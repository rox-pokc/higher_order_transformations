using LinearAlgebra
using Combinatorics

# Define basic Clifford generators
const H = ComplexF64.(1/sqrt(2) * [1 1; 1 -1])
const S = ComplexF64.([1 0; 0 im])
const I = ComplexF64.([1 0; 0 1])

const GENERATORS = Dict(
    "I" => I,
    "H" => H,
    "S" => S
)

# Normalize global phase
function normalize_global_phase(U::Matrix{ComplexF64})
    phase = det(U)^(1/2)
    return U / phase
end

# Main function
function generate_cliffords(depth::Int)
    seen = Dict{Vector{ComplexF64}, Vector{String}}()
    gate_labels = collect(keys(GENERATORS))
    
    for seq in Iterators.product(ntuple(_ -> gate_labels, depth)...)
        U = I
        for g in seq
            U = GENERATORS[g] * U
        end
        U_norm = normalize_global_phase(U)
        key = vec(round.(U_norm; digits=6))
        if !haskey(seen, key)
            seen[key] = collect(seq)
        end
    end
    return seen
end


cliffs = generate_cliffords(6)

for (k, v) in cliffs
    println(k)
    println(join(v, " → "))
end

println("Total unique Clifford unitaries found: ", length(cliffs))
