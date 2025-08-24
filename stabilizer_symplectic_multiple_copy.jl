using QuantumClifford
using LinearAlgebra
using JSON
using Combinatorics
using Base.Threads

# include("clifford_generation.jl")  # gate sequence generator

include("src/Config.jl")
using .Config  
# include("src/CliffordData.jl")
# using .CliffordData  

mutable struct Comb
    stab::Stabilizer
    enc::Vector{Tuple{String, Vararg{Int}}}
    dec::Vector{Tuple{String, Vararg{Int}}}
end

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

const CLIFFORD_CACHE = Dict{Int, Any}()

function get_cliffords(k::Int)
    if haskey(CLIFFORD_CACHE, k)
        return CLIFFORD_CACHE[k]
    else
        filename = "cliffords_generation_outputs/cliffords_depth=6_#qubits=$(k).json"
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
    return eval(Meta.parse("P\"$stab_str\""))
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

function apply_comb(U, func::Function, d, k)
    # Initialize stabilizer state
    stab = build_stabilizer(d, k)

    comb = Comb(stab, [], [])

    # Get generated gate sequences
    # cliffords = generate_cliffords(6)
    cliffords = get_cliffords(k)

    # enc_target_qubit = 5 #TODO HARD CODED HERE
    # enc_target_qubit_1 = 9
    # enc_target_qubit_2 = 11

    order = 3
    enc_comb_array = []
    for seq in cliffords 
        copy_stab = deepcopy(stab)
        for (g, qubits...) in seq
            if length(qubits) == 1
                qubit = qubits[1]
                apply!(copy_stab, CLIFFORD_GATES[g](1 + 2 * d * k + d * (qubit - 1)))
            else
                ctrl = qubits[1]
                trgt = qubits[2]
                apply!(copy_stab, CLIFFORD_GATES[g](1 + 2 * d * k + d * (ctrl - 1), 1 + 2 * d * k + d * (trgt - 1)))
            end
        end
        push!(enc_comb_array, Comb(copy_stab, seq, []))
    end
    enc_comb_array = enc_comb_array[1:2]

    # Bell measure 3rd and 5th qubits / 5th, 7th, 9th, 11th
    # 3 = 1 + 1 * d * k + d * (qubit - 1)
    # 5 = 1 + 2 * d * k + d * (qubit - 1)

    # 5 = 1 + 1 * d * k + d * (qubit - 1)
    # 7 = 1 + 1 * d * k + d * (qubit - 1)
    # 9 = 1 + 2 * d * k + d * (qubit - 1)
    # 11 = 1 + 2 * d * k + d * (qubit - 1)    

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
    
    for comb in enc_comb_array
        comb.stab = project!(comb.stab, P_XX)[1]
        comb.stab = project!(comb.stab, P_ZZ)[1]
    end

    #Decode
    # dec_target_qubit = 1
    # dec_target_qubit_1 = 1
    # dec_target_qubit_2 = 3
    order = 1
    dec_comb_array = []
    for comb in enc_comb_array
        for seq in cliffords
            copy_stab = deepcopy(comb.stab)
            for (g, qubits...) in seq
                if length(qubits) == 1
                    qubit = qubits[1]
                    apply!(copy_stab, CLIFFORD_GATES[g](1 + d * k * (order - 1) + d * (qubit - 1)))
                else
                    ctrl = qubits[1]
                    trgt = qubits[2]
                    apply!(copy_stab, CLIFFORD_GATES[g](1 + d * k * (order - 1) + d * (ctrl - 1), 1 + d * k * (order - 1) + d * (trgt - 1)))
                end
            end
            push!(dec_comb_array, Comb(copy_stab, comb.enc, seq))
        end
    end

    #choi function state
    # U_combned = Config.I

    U_arr = [CLIFFORD_MATRICES[g] for (g, qubit) in U]
    U_combined = reduce(*, U_arr)

    # for i in length(U)
        
    #     # println(U_arr)
    #     # println(U_arr)
    #     # tensor = reduce(kron, U_arr)
    #     # println(U)
        
    # end
    # FUNCTION DEFINITION HERE
    # choi_func = get_choi(func(kron(U_combined, U_combined)), d^2)
    choi_func = get_choi(func(U_combined), d)

    # selected_qubits = [1, 6]
    # selected_qubits = [1, 3, 10, 12]
    selected_qubits = Vector{Int64}()
    first_order = 1
    second_order = 3
    for qubit in 1:k
        first_qubit = 1 + (first_order - 1) * d * k + d * (qubit - 1)
        second_qubit = 1 + (second_order - 1) * d * k + d * (qubit - 1) + 1
        push!(selected_qubits, first_qubit)
        push!(selected_qubits, second_qubit)
    end

    ans_combs = []

    for comb in dec_comb_array

        # Reduce stabilizers
        canonicalize!(comb.stab)

        # Inject U - Clifford
        # 7
        # 13 & 15
        order = 4
        for qubit in 1:k
            qubit = 1 + (order - 1) * d * k + d * (qubit - 1)
            for (g, ind) in U
                apply!(comb.stab, CLIFFORD_GATES[g](qubit))
            end
        end

        # Bell measure 4th and 8th qubits
        # 6, 14 / Bell measure 8, 16

        first_order = 2
        second_order = 4
        for qubit in 1:k
            first_qubit = 1 + (first_order - 1) * d * k + d * (qubit - 1) + 1
            second_qubit = 1 + (second_order - 1) * d * k + d * (qubit - 1) + 1
            # println(first_qubit)
            # println(second_qubit)
            
            P_XX, P_ZZ = build_Bell_meas_2q(first_qubit, second_qubit, total)
            # println(P_XX)
            # println(P_ZZ)
        
            comb.stab = project!(comb.stab, P_XX)[1]
            comb.stab = project!(comb.stab, P_ZZ)[1]
        end
        
        # Bell measure 2nd and 7th qubits
        # 2, 13 / Bell measure 4, 15
        first_order = 1
        second_order = 4
        qubit = 1
        for qubit in 1:k
            first_qubit = 1 + (first_order - 1) * d * k + d * (qubit - 1) + 1
            second_qubit = 1 + (second_order - 1) * d * k + d * (qubit - 1)
            # println(first_qubit)
            # println(second_qubit)
            
            P_XX, P_ZZ = build_Bell_meas_2q(first_qubit, second_qubit, total)
            # println(P_XX)
            # println(P_ZZ)
        
            comb.stab = project!(comb.stab, P_XX)[1]
            comb.stab = project!(comb.stab, P_ZZ)[1]
        end

        canonicalize!(comb.stab)

        #Checking stabilizers
        partial_strs = extract_partial_stabilizers(comb.stab, selected_qubits)
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
            push!(ans_combs, comb)
        end
    end

    return ans_combs
end

import Base: isequal, hash

function isequal(a::Comb, b::Comb)
    return a.enc == b.enc && b.dec == a.dec
end

function hash(c::Comb, h::UInt)
    return hash((c.enc, c.dec), h)
end

function find_top_shared_elements(arrays; top_n::Int = 3)
    # Create a dictionary to count in how many arrays each Comb appears
    counts = Dict{Comb, Int}()

    for arr in arrays
        seen = Set{Comb}() 
        for comb in arr
            if comb ∉ seen
                counts[comb] = get(counts, comb, 0) + 1
                push!(seen, comb)
            end
        end
    end

    # Sort by how many arrays the comb appears in
    sorted = sort(collect(counts), by = x -> -x[2])

    # Take top N
    top_shared = first(sorted, top_n)

    return top_shared
end

function find_top_shared_until_unrepresented(arrays; top_n::Int = 3)
    counts = Dict{Comb, Int}()
    presence = Dict{Comb, Set{Int}}()

    # Track which arrays each Comb appears in
    for (i, arr) in enumerate(arrays)
        seen = Set{Comb}()
        for comb in arr
            if comb ∉ seen
                counts[comb] = get(counts, comb, 0) + 1
                if !haskey(presence, comb)
                    presence[comb] = Set{Int}()
                end
                push!(presence[comb], i)
                push!(seen, comb)
            end
        end
    end

    # Sort combinations by how many arrays they appear in
    sorted = sort(collect(counts), by = x -> -x[2])

    top_comb_set = Set{Comb}()
    covered_arrays = Set{Int}()
    result = []

    for (comb, _) in sorted
        push!(top_comb_set, comb)
        union!(covered_arrays, presence[comb])
        push!(result, comb)

        if length(covered_arrays) < length(arrays)
            continue
        end

        if length(result) >= top_n
            break
        end
    end

    result_with_counts = [(comb, length(presence[comb]), collect(presence[comb])) for comb in result]
    return result_with_counts
end

function count_enc_dec_presence(arrays)
    presence = Dict{Tuple{Vector{Tuple{String, Int64}}, Vector{Tuple{String, Int64}}}, Set{Int}}()

    for (i, arr) in enumerate(arrays)
        seen = Set{Tuple{Vector{Tuple{String, Int64}}, Vector{Tuple{String, Int64}}}}()
        for comb in arr
            key = (comb.enc, comb.dec)
            if key ∉ seen
                if !haskey(presence, key)
                    presence[key] = Set{Int}()
                end
                push!(presence[key], i)
                push!(seen, key)
            end
        end
    end

    return [(key, length(presence[key]), sort(collect(presence[key]))) for key in keys(presence)]
end

function equal_up_to_global_phase(A::AbstractMatrix, B::AbstractMatrix; atol=1e-8)
    # Flatten and find first nonzero element in B to estimate phase
    idx = findfirst(!≈(0), B)
    λ = A[idx] / B[idx]  # estimated global phase
    return isapprox(A, λ * B; atol=atol)
end

start_time = time_ns()

d = 2
k = 1
func = conj

println(Threads.nthreads())

file = open("output_func=$(func)_d=$(d)_k=$(k)_depth=6_old_version.txt", "w")

println(file, "Start:")
flush(file)

arrays = []
# cliffords = generate_cliffords(6)
# cliffords = read_sequences_from_file("cliffords_depth=8_#qubits=2.json")
U_cliffords = get_1q_cliffords()

global i = 1
for U in U_cliffords
    if (i==1)
        start_time_comb = time_ns()
        combs = apply_comb(U, func, d, k)
        end_time_comb = time_ns()

        println(file, "Done U in seconds: ", (end_time_comb - start_time_comb) / 1e9)
        flush(file)

        println(file, length(combs))
        flush(file)
        println(file, U)
        flush(file)
        seen = Set{Tuple{Vector{Tuple{String, Int64}}, Vector{Tuple{String, Int64}}}}()
        for comb in combs
            key = (comb.enc, comb.dec)
            if key ∉ seen
                push!(seen, key)
            end
        end
        for comb in seen
            println(comb)
            flush(file)
        end

        global i += 1

        push!(arrays, combs)
    end
end

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

function compute_element_array_map(elements, arrays)
    # Map from element => Set of array indices
    element_to_arrays = Dict(e => Int[] for e in elements)

    for (i, arr) in enumerate(arrays)
        for el in arr
            push!(element_to_arrays[el], i)
        end
    end

    return element_to_arrays
end

function compute_array_intersections(arrays)
    n = length(arrays)
    intersection_matrix = zeros(Int, n, n)

    for i in 1:n
        for j in i+1:n
            common = length(intersect(arrays[i], arrays[j]))
            intersection_matrix[i, j] = common
            intersection_matrix[j, i] = common
        end
    end

    return intersection_matrix
end


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

function intersection_lengths(arrays, group_size::Int)
    sets = [Set(a) for a in arrays]
    lens = Set{Int}()
    
    for idx_group in combinations(1:length(sets), group_size)
        common = reduce(intersect, sets[idx_group])
        push!(lens, length(common))
    end

    return lens
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