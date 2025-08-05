using LinearAlgebra
using Combinatorics

# Define actual 1-qubit Clifford gates
const H = ComplexF64.(1/sqrt(2) * [1 1; 1 -1])
const S = ComplexF64.([1 0; 0 im])
const I = ComplexF64.([1 0; 0 1])
const X = ComplexF64.([0 1; 1 0])
const Y = ComplexF64.([0 -im; im 0])
const Z = ComplexF64.([1 0; 0 -1])

# Gate dictionary
const unitary_gates = Dict(
    "I" => I,
    "H" => H,
    "S" => S,
    "X" => X,
    "Y" => Y,
    "Z" => Z,
)

# Compose unitary matrix from a sequence
function compose_unitary(sequence::Vector{String})
    result = I
    for g in reverse(sequence)  # apply right-to-left
        result = unitary_gates[g] * result
    end
    return result
end

function normalize_global_phase(U::Matrix{ComplexF64})
    for i in 1:size(U,1), j in 1:size(U,2)
        if abs(U[i,j]) > 1e-8
            phase = angle(U[i,j])
            return U * exp(-im * phase)
        end
    end
    return U  # fallback for zero matrix (should not happen with valid unitaries)
end

function equivalent_up_to_global_phase(U1, U2; atol=1e-8)
    U1n = normalize_global_phase(U1)
    U2n = normalize_global_phase(U2)
    return isapprox(U1n, U2n; atol=atol)
end

# Generate unique Clifford unitaries
function generate_unique_unitaries(depth::Int)
    gate_names = collect(keys(unitary_gates))
    seen = Vector{Matrix{ComplexF64}}()
    circuits = Dict{Matrix{ComplexF64}, Vector{String}}()

    for seq in Iterators.product(ntuple(_ -> gate_names, depth)...)
        gatelist = collect(seq)
        U = compose_unitary(gatelist)
        
        if !any(existing -> equivalent_up_to_global_phase(existing, U), seen)
            push!(seen, U)
            circuits[U] = gatelist
        end
    end

    return circuits
end

# Run
depth = 4
cliffords = generate_unique_unitaries(depth)

# Print results
for (U, seq) in cliffords
    println("Sequence: ", join(seq, " → "))
    println("Matrix:")
    display(U)
    println()
end

println("Total unique Clifford unitaries found: ", length(cliffords))
