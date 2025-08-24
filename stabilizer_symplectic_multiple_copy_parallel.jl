using QuantumClifford
using LinearAlgebra
using JSON
using Combinatorics
using Base.Threads
using Dates

include("src/Config.jl")
using .Config  

struct Comb
    stab::Stabilizer
    enc::Vector{Tuple{String, Vararg{Int}}}
    dec::Vector{Tuple{String, Vararg{Int}}}

    function Comb(stab, enc, dec)
        new(deepcopy(stab), deepcopy(enc), deepcopy(dec))
    end
end

include("post_processing_results.jl")

# Dynamic creation of initial stabilizers
function build_stabilizer(d::Int, k::Int)
    lines = String[]
    total = 4 * k * d + 1

    for i in 1:2:(4 * k * d)
        push!(lines, "+" * repeat("I", i - 1) * repeat("X", d) * repeat("I", total - (i + d)))
        push!(lines, "+" * repeat("I", i - 1) * repeat("Z", d) * repeat("I", total - (i + d)))
    end

    stab_str = join(lines, "\n")
    return eval(Meta.parse("S\"\"\"$stab_str\"\"\""))
end

function read_sequences_from_file(filename::String)
    open(filename, "r") do io
        raw = JSON.parse(io)
        return [Tuple.(s) for s in raw]
    end
end

const CLIFFORD_CACHE = Dict{Int, Vector{Vector{Tuple{String, Vararg{Int}}}}}()

function get_cliffords(k::Int)
    if haskey(CLIFFORD_CACHE, k)
        return CLIFFORD_CACHE[k]
    else
        filename = "cliffords_generation_outputs/cliffords_depth=4_#qubits=$(k).json"
        sequences = read_sequences_from_file(filename)
        CLIFFORD_CACHE[k] = sequences
        return sequences
    end
end

function get_1q_cliffords()
    filename = "cliffords_generation_outputs/cliffords_depth=6_#qubits=1.json"
    sequences = read_sequences_from_file(filename)
    return sequences
end

function build_meas_stabilizer(total::Int, indexes::Vector{Int64}, meas::Char)
    stab_char = fill('I', total)
    @inbounds for idx in indexes
        stab_char[idx] = meas
    end
    stab_str = String(stab_char)
    return stab_str
end

function get_choi(A::AbstractMatrix, d::Int)
    # |phi+> = 1/sqrt(d)sum(|ii>)i=0..d-1
    phi_plus = zeros(ComplexF64, d^2)
    for i in 1:(d+1):d^2
        phi_plus[i] = 1 / sqrt(d)
    end

    # |phi+><phi+|
    density_matrix_phi_plus = phi_plus * phi_plus'

    return d * kron(Matrix{ComplexF64}(LinearAlgebra.I, d, d), A) * density_matrix_phi_plus * kron(Matrix{ComplexF64}(LinearAlgebra.I, d, d), A')
end

function build_Bell_meas_2q(first_qubit, second_qubit, total)
    indexes = Vector{Int64}()
        
    push!(indexes, first_qubit)
    push!(indexes, second_qubit)

    P_XX = build_meas_stabilizer(total, indexes, 'X')
    P_ZZ = build_meas_stabilizer(total, indexes, 'Z')

    return P_XX, P_ZZ
end

function extract_partial_stabilizers(stab, qubit_indices::Vector{Int})
    partials = String[]
    for g in stab
        pstr = string(g)
        selected = join(pstr[i + 2] for i in qubit_indices)
        push!(partials, pstr[1] * selected)
    end
    return unique(partials)
end

const PAULI_CACHE = PauliOperator[]

function build_pauli_cache!(strings::Vector{String})
    empty!(PAULI_CACHE)
    for s in strings
        push!(PAULI_CACHE, eval(Meta.parse("P\"$s\"")))
    end
    return nothing
end

function build_stabilizer_cache(d, k)
    all_strings = Vector{String}()

    first_order = 2
    second_order = 3
    total = 4 * k * d
    indexes = Vector{Int64}()
    for qubit in 1:k
        first_qubit = 1 + (first_order - 1) * d * k + d * (qubit - 1)
        second_qubit = 1 + (second_order - 1) * d * k + d * (qubit - 1)
        push!(indexes, first_qubit)
        push!(indexes, second_qubit)
    end

    P_XX = build_meas_stabilizer(total, indexes, 'X')
    P_ZZ = build_meas_stabilizer(total, indexes, 'Z')

    push!(all_strings, P_XX)
    push!(all_strings, P_ZZ)

    first_order = 2
    second_order = 4
    for qubit in 1:k
        first_qubit = 1 + (first_order - 1) * d * k + d * (qubit - 1) + 1
        second_qubit = 1 + (second_order - 1) * d * k + d * (qubit - 1) + 1
        
        P_XX, P_ZZ = build_Bell_meas_2q(first_qubit, second_qubit, total)
        
        push!(all_strings, P_XX)
        push!(all_strings, P_ZZ)
    end

    # Bell measure
    first_order = 1
    second_order = 4
    for qubit in 1:k
        first_qubit = 1 + (first_order - 1) * d * k + d * (qubit - 1) + 1
        second_qubit = 1 + (second_order - 1) * d * k + d * (qubit - 1)
        
        P_XX, P_ZZ = build_Bell_meas_2q(first_qubit, second_qubit, total)
    
        push!(all_strings, P_XX)
        push!(all_strings, P_ZZ)
    end

    build_pauli_cache!(all_strings)
end

@inline pauli_from_cache(i::Int) = @inbounds PAULI_CACHE[i]

function apply_comb(U, choi_func, d, k)

    #amount of threads
    nt = nthreads()
    # per-thread buffers
    buffers = [[] for _ in 1:nt]

    # Get generated gate sequences
    cliffords = get_cliffords(k)

    # Initialize stabilizer state
    global_stab = build_stabilizer(d, k)

    # Selecting qubits that contain output
    selected_qubits = Vector{Int64}()
    first_order = 1
    second_order = 3
    for qubit in 1:k
        first_qubit = 1 + (first_order - 1) * d * k + d * (qubit - 1)
        second_qubit = 1 + (second_order - 1) * d * k + d * (qubit - 1) + 1
        push!(selected_qubits, first_qubit)
        push!(selected_qubits, second_qubit)
    end

    @threads for enc_seq in cliffords
        t_id = threadid()

        stab = deepcopy(global_stab)

        # Encode
        order = 3
        for (g, qubits...) in enc_seq
            if length(qubits) == 1
                qubit = qubits[1]
                apply!(stab, CLIFFORD_GATES[g](1 + 2 * d * k + d * (qubit - 1)))
            else
                ctrl = qubits[1]
                trgt = qubits[2]
                apply!(stab, CLIFFORD_GATES[g](1 + 2 * d * k + d * (ctrl - 1), 1 + 2 * d * k + d * (trgt - 1)))
            end
        end

        # Bell measure 3rd and 5th qubits / 5th, 7th, 9th, 11th
        stab = project!(stab, pauli_from_cache(1))[1]
        stab = project!(stab, pauli_from_cache(2))[1]

        # Decode
        for dec_seq in cliffords

            dec_stab = deepcopy(stab)

            order = 1
            for (g, qubits...) in dec_seq
                if length(qubits) == 1
                    qubit = qubits[1]
                    apply!(dec_stab, CLIFFORD_GATES[g](1 + d * k * (order - 1) + d * (qubit - 1)))
                else
                    ctrl = qubits[1]
                    trgt = qubits[2]
                    apply!(dec_stab, CLIFFORD_GATES[g](1 + d * k * (order - 1) + d * (ctrl - 1), 1 + d * k * (order - 1) + d * (trgt - 1)))
                end
            end

            # Inject U - Clifford
            order = 4
            for qubit in 1:k
                qubit = 1 + (order - 1) * d * k + d * (qubit - 1)
                for (g, index) in U
                    apply!(dec_stab, CLIFFORD_GATES[g](qubit))
                end
            end

            # Bell measure 
            for i in 3:length(PAULI_CACHE)
                dec_stab = project!(dec_stab, pauli_from_cache(i))[1]
            end
            # dec_stab = project!(dec_stab, pauli_from_cache(3))[1]
            # dec_stab = project!(dec_stab, pauli_from_cache(4))[1]
            # dec_stab = project!(dec_stab, pauli_from_cache(5))[1]
            # dec_stab = project!(dec_stab, pauli_from_cache(6))[1]

            # println("maes ", dec_stab)
        
            # Bell measure
            # comb.stab = project!(comb.stab, pauli_from_cache(7))[1]
            # comb.stab = project!(comb.stab, pauli_from_cache(8))[1]
            # comb.stab = project!(comb.stab, pauli_from_cache(9))[1]
            # comb.stab = project!(comb.stab, pauli_from_cache(10))[1]

            # Reduce stabilizers
            canonicalize!(dec_stab)

            #Checking stabilizers
            partial_strs = extract_partial_stabilizers(dec_stab, selected_qubits)
            nontrivial = filter(x -> x[2:length(selected_qubits)+1] != repeat("_", length(selected_qubits)), partial_strs)

            flag = 0
            for out_stab in nontrivial
                out_stab = replace(out_stab, '_' => 'I')
                sign = out_stab[1] == '-' ? -1 : 1
                body = out_stab[2:end]
                ops = [PAULI_GATES[c] for c in body]
                operator = foldl(kron, ops)
                transformed = operator * choi_func
                if (isapprox(transformed, sign * choi_func))
                    flag += 1 
                end
            end
            if (flag == length(nontrivial))
                comb = Comb(dec_stab, enc_seq, dec_seq)
                push!(buffers[t_id], comb)
            end
            
        end
        
    end
    return vcat(buffers...)
end

# BLAS.set_num_threads(1)
# versioninfo()

start_time = time_ns()

d = 2
k = 2
func = conj

println(Threads.nthreads())

file = open("outputs/k=$(k)/output_func=$(func)_d=$(d)_k=$(k)_depth=4.txt", "w")

println(file, "Start:")
flush(file)

println("Started at: ", now())

build_stabilizer_cache(d, k)

arrays = []
U_cliffords = get_1q_cliffords()
global i = 1
for U in U_cliffords
    combs = []

    U_arr = [CLIFFORD_MATRICES[g] for (g, qubit) in U]
    U_combined = reduce(*, U_arr)
    choi_func_U = get_choi(func(kron(U_combined, U_combined)), d^2)
    # choi_func_U = get_choi(func(U_combined), d)

    println("Processing ", i, "/", length(U_cliffords))

    start_time_comb = time_ns()
    combs = apply_comb(U, choi_func_U, d, k)
    end_time_comb = time_ns()

    println(file, "Done U in seconds: ", (end_time_comb - start_time_comb) / 1e9)
    flush(file)

    global i += 1

    println(file, length(combs))
    flush(file)
    println(file, U)
    flush(file)
    seen = Set{Tuple{Vector{Tuple{String, Vararg{Int}}}, Vector{Tuple{String, Vararg{Int}}}}}()
    for comb in combs
        key = (comb.enc, comb.dec)
        if key ∉ seen
            push!(seen, key)
        end
    end
    for comb in seen
        println(file, comb)
        flush(file)
    end

    push!(arrays, combs)
end

println(file)
println(file, "---------")
println(file)

println(file, "arr length: ", length(arrays))
flush(file)

# FIND ALL COMBS WHICH GIVES TRANSFORMSTIONS
all = count_enc_dec_presence(arrays)

println(file, "all: ", length(all))
flush(file)

global ten = 0

ten_arr = []

for ((enc, dec), count, indices) in all
    if count == 10
        # println("enc: ", enc, ", dec: ", dec, " → appears in $count arrays: ", indices)
        push!(ten_arr, indices)
        global ten += 1
    end 
    println(file, "enc: ", enc, ", dec: ", dec, " → appears in $count arrays: ", indices)
    flush(file)
end
println(file, "Appeared in 10 arr amount:", ten)
flush(file)

# -------

elements = collect(1:24)
println(file, ten_arr)
flush(file)

element_map = compute_element_array_map(elements, ten_arr)
intersection_matrix = compute_array_intersections(ten_arr)

println(file, "Each element belongs to:")
flush(file)
for (e, arrs) in element_map
    println(file, "$e → Arrays: ", sort(arrs))
    flush(file)
end

println(file, "\nPairwise intersection matrix between arrays:")
flush(file)
for row in eachrow(intersection_matrix)
    println(file, row)
    flush(file)
end

for i in 1:10
    println(file, intersection_lengths(ten_arr, i))
    flush(file)
end

end_time = time_ns()

println(file, "Done in seconds: ", (end_time - start_time) / 1e9)
flush(file)


# FIND COMBS UNTIL FIND TRANSFORMATION FOR ALL CLIFFORD
# top_combs = find_top_shared_until_unrepresented(arrays)

# for (comb, count, indices) in find_top_shared_until_unrepresented(arrays)
#     println("enc: ", comb.enc, ", dec: ", comb.dec,
#             " → shared in $count arrays: ", sort(indices))

#     println("enc: ")
#     enc = reduce(*, reverse([CLIFFORD_MATRICES[c] for c in comb.enc]))
#     show(stdout, "text/plain", enc)
#     println()

#     println("dec: ")
#     dec = reduce(*, reverse([CLIFFORD_MATRICES[c] for c in comb.dec]))
#     show(stdout, "text/plain", dec)
#     println()
# end
# -------

# FIND COMBS WHICH GIVES MOST CLIFFORDS
# top_shared = find_top_shared_elements(arrays, top_n = 5)

# for (comb, count) in top_shared
#     println("Appears in $count arrays: ", comb.enc, " → ", comb.dec)
# end
# -------