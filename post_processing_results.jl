using LinearAlgebra
using Combinatorics

import Base: isequal, hash

# struct Comb
#     stab::Stabilizer
#     enc::Vector{Tuple{String, Vararg{Int}}}
#     dec::Vector{Tuple{String, Vararg{Int}}}
# end

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
    presence = Dict{Tuple{Vector{Tuple{String, Vararg{Int}}}, Vector{Tuple{String, Vararg{Int}}}}, Set{Int}}()

    for (i, arr) in enumerate(arrays)
        seen = Set{Tuple{Vector{Tuple{String, Vararg{Int}}}, Vector{Tuple{String, Vararg{Int}}}}}()
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

function intersection_lengths(arrays, group_size::Int)
    sets = [Set(a) for a in arrays]
    lens = Set{Int}()
    
    for idx_group in combinations(1:length(sets), group_size)
        common = reduce(intersect, sets[idx_group])
        push!(lens, length(common))
    end

    return lens
end