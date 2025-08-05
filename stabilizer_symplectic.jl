using QuantumClifford
using LinearAlgebra
include("test.jl")  # gate sequence generator

# Gate map
const CLIFFORD_GATES = Dict(
    "H" => sHadamard,
    "S" => sPhase,
    "I" => sId1
)

X = [0 1; 1 0]
Y = [0 -im; im 0]
Z = [1 0; 0 -1]
I = [1 0; 0 1]
S = [1 0; 0 im]
H = (1 / sqrt(2)) * [1  1; 1 -1]

const PAULI_GATES = Dict(
    'I' => I,
    'X' => X,
    'Y' => Y,
    'Z' => Z
)

const CLIFFORD_MATRICES = Dict(
    "H" => H,
    "S" => S,
    "I" => I
)

mutable struct Comb
    stab::Stabilizer
    enc::Vector{String}
    dec::Vector{String}
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

function extract_partial_stabilizers(stab, qubit_indices::Vector{Int})
    partials = String[]
    for g in stab
        pstr = string(g)
        selected = join(pstr[i + 2] for i in qubit_indices)
        push!(partials, pstr[1] * selected)
    end
    return unique(partials)
end

function apply_comb(U, func::Function, d)
    # Initialize stabilizer state
    stab = S"XXIIIIII
        ZZIIIIII
        IIXXIIII
        IIZZIIII
        IIIIXXII
        IIIIZZII
        IIIIIIXX
        IIIIIIZZ"
    comb = Comb(stab, [], [])

    # Get generated gate sequences
    results = generate_cliffords_full(6)

    enc_target_qubit = 5
    enc_comb_array = []
    for ((px, pz), seq) in results
        copy_stab = deepcopy(stab)
        for g in seq
            apply!(copy_stab, CLIFFORD_GATES[g](enc_target_qubit))
        end
        push!(enc_comb_array, Comb(copy_stab, seq, []))
    end

    # Bell measure 3rd and 5th qubits
    for comb in enc_comb_array
        comb.stab = project!(comb.stab, P"IIXIXIII")[1]
        comb.stab = project!(comb.stab, P"IIZIZIII")[1]
    end

    #Decode
    dec_target_qubit = 1
    dec_comb_array = []
    for comb in enc_comb_array
        for ((px, pz), seq) in results
            copy_stab = deepcopy(comb.stab)
            for g in seq
                apply!(copy_stab, CLIFFORD_GATES[g](dec_target_qubit))
            end
            push!(dec_comb_array, Comb(copy_stab, comb.enc, seq))
        end
    end

    #choi function state
    U_arr = [CLIFFORD_MATRICES[c] for c in U]
    U_combined = reduce(*, U_arr)
    # FUNCTION DEFINITION HERE
    choi_func = get_choi(func(U_combined), d)

    selected_qubits = [1, 6]

    ans_combs = []

    for comb in dec_comb_array

        # Reduce stabilizers
        canonicalize!(comb.stab)

        # Inject U - Clifford
        for g in U
            apply!(comb.stab, CLIFFORD_GATES[g](7))
        end

        # Bell measure 4th and 8th qubits
        comb.stab = project!(comb.stab, P"IIIXIIIX")[1]
        comb.stab = project!(comb.stab, P"IIIZIIIZ")[1]

        # Bell measure 2nd and 7th qubits
        comb.stab = project!(comb.stab, P"IXIIIIXI")[1]
        comb.stab = project!(comb.stab, P"IZIIIIZI")[1]

        # println("Full circuit:\n", z2z7)

        canonicalize!(comb.stab)

        #Checking stabilizers
        partial_strs = extract_partial_stabilizers(comb.stab, selected_qubits)
        nontrivial = filter(x -> x[2:3] != "__", partial_strs)

        flag = 0
        for out_stab in nontrivial
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
    presence = Dict{Tuple{Vector{String}, Vector{String}}, Set{Int}}()

    for (i, arr) in enumerate(arrays)
        seen = Set{Tuple{Vector{String}, Vector{String}}}()
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



d = 2

arrays = []
results = generate_cliffords_full(6)
for ((px, pz), seq) in results
    push!(arrays, apply_comb(seq, inv, d))
end
println(length(arrays))

# FIND ALL COMBS WHICH GIVES TRANSFORMSTIONS
all = count_enc_dec_presence(arrays)

println(length(all))

global ten = 0

ten_arr = []

for ((enc, dec), count, indices) in all
    if count == 10
        println("enc: ", enc, ", dec: ", dec, " → appears in $count arrays: ", indices)
        push!(ten_arr, indices)
        global ten += 1
    end 
    # println("enc: ", enc, ", dec: ", dec, " → appears in $count arrays: ", indices)
end
println("Appeared in 10 arr amount:", ten)
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
println(ten_arr)

element_map = compute_element_array_map(elements, ten_arr)
intersection_matrix = compute_array_intersections(ten_arr)

println("Each element belongs to:")
for (e, arrs) in element_map
    println("$e → Arrays: ", sort(arrs))
end

println("\nPairwise intersection matrix between arrays:")
for row in eachrow(intersection_matrix)
    println(row)
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
    println(intersection_lengths(ten_arr, i))
end


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