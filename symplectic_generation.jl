using LinearAlgebra
using Combinatorics

# Pauli basis order: X = (1,0), Z = (0,1)
# Pauli vector: [x; z] ∈ F₂²
MOD2 = x -> mod.(x, 2)

# Symplectic action: maps [x; z] → A * [x; z]
# For 1-qubit Clifford, A is 2x2 over F₂
const H_sym = [0 1; 1 0]  # Hadamard swaps X ↔ Z
const S_sym = [1 0; 1 1]  # Phase: X → Y, Z unchanged
const I_sym = [1 0; 0 1]  # Identity
# const X_sym = [1 0; 0 1]        
# const Z_sym = [1 0; 0 1]
# const Y_sym = [1 0; 0 1]

# Actual 2x2 unitaries
const H = 1/sqrt(2) * [1 1; 1 -1]
const S = [1 0; 0 im]
const I = [1 0; 0 1]
# const X = [0 1; 1 0]
# const Z = [1 0; 0 -1]
# const Y = [0 -im; im 0]

# Named dictionary of generators
const generators = Dict(
    "I" => I_sym,
    "H" => H_sym,
    "S" => S_sym,
    # "X" => X_sym,
    # "Y" => Y_sym,
    # "Z" => Z_sym,
)

const unitary_gates = Dict(
    "I" => I,
    "H" => H,
    "S" => S,
    # "X" => X,
    # "Y" => Y,
    # "Z" => Z,
)

# Function to compose symplectic matrices mod 2
function compose_symplectic(sequence)
    result = I_sym
    for g in sequence
        result = MOD2(g * result)
    end
    return result
end

# Compose actual unitary matrix from gate sequence
function compose_unitary(sequence::Vector{String})
    result = I
    for g in reverse(sequence)  # apply right-to-left
        result = unitary_gates[g] * result
    end
    return result
end

# Generate all sequences of length `depth` over `generators`
function generate_unique_symplectic(depth)
    gate_names = collect(keys(generators))
    seen = Set{Matrix{Int}}()
    unique_circuits = Dict{Matrix{Int}, Vector{String}}()

    for seq in Iterators.product(ntuple(_ -> gate_names, depth)...)
        gate_matrices = [generators[g] for g in seq]
        composed = compose_symplectic(reverse(gate_matrices))  # reverse to apply in correct order
        if !any(x -> composed == x, seen)
            push!(seen, composed)
            unique_circuits[composed] = collect(seq)
        end
    end

    return unique_circuits
end

# Example usage: all unique 1-qubit Cliffords of depth 3
depth = 5
unique_circs = generate_unique_symplectic(depth)

for (sym, seq) in unique_circs
    println("Sequence: ", join(seq, " → "))
    println("Matrix:")
    display(compose_unitary(seq))
    println()
end

println("Total unique 1-qubit Cliffords found: ", length(unique_circs))
