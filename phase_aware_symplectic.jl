using LinearAlgebra
using Combinatorics

# --- Pauli Binary Vectors and Phase-Aware Tableau Representation ---

# Pauli operators are represented as (x, z) ∈ F₂²
# Each row is a stabilizer generator: [x z r] where r ∈ Z₄ encodes phase:
# 0 → +1, 1 → +i, 2 → -1, 3 → -i

struct Tableau
    x::Vector{Int}
    z::Vector{Int}
    r::Int  # phase in Z₄
end

# Define basic Clifford gates and their symplectic + phase update rules

const H = 1/sqrt(2) * [1 1; 1 -1]
const S = [1 0; 0 im]
const I = [1 0; 0 1]

# Initialize stabilizers for |0⟩ state: stabilized by Z
function initial_tableau()
    Tableau([0], [1], 0)  # corresponds to Z with +1 eigenvalue
end

# Function to apply H gate to 1-qubit tableau
function apply_H(t::Tableau)
    x, z, r = t.x[1], t.z[1], t.r
    new_x, new_z = z, x
    new_r = (r + 2 * x * z) % 4  # phase update rule for H
    return Tableau([new_x], [new_z], new_r)
end

# Function to apply S gate to 1-qubit tableau
function apply_S(t::Tableau)
    x, z, r = t.x[1], t.z[1], t.r
    new_x, new_z = x, (z + x) % 2
    new_r = (r + x * z) % 4  # phase update rule for S
    return Tableau([new_x], [new_z], new_r)
end

# Compose gate sequence and return final tableau
function apply_sequence_to_tableau(seq::Vector{String})
    t = initial_tableau()
    for g in reverse(seq)
        if g == "H"
            t = apply_H(t)
        elseif g == "S"
            t = apply_S(t)
        elseif g == "I"
            continue
        else
            error("Unsupported gate: $g")
        end
    end
    return t
end

# Generate unique tableaus of fixed depth
function generate_unique_tableaus(depth)
    gate_names = ["H", "S", "I"]
    seen = Vector{Tableau}()
    results = Dict{Tableau, Vector{String}}()

    function tableau_equal(t1::Tableau, t2::Tableau)
        return t1.x == t2.x && t1.z == t2.z && t1.r % 4 == t2.r % 4
    end

    for seq in Iterators.product(ntuple(_ -> gate_names, depth)...)
        gatelist = collect(seq)
        t = apply_sequence_to_tableau(gatelist)
        if !any(existing -> tableau_equal(existing, t), seen)
            push!(seen, t)
            results[t] = gatelist
        end
    end

    return results
end

# Run the generator
depth = 4
unique_tableaus = generate_unique_tableaus(depth)

# Print the results
for (stab, seq) in unique_tableaus
    println("Sequence: ", join(seq, " → "))
    println("x = ", stab.x, ", z = ", stab.z, ", phase r = ", stab.r)
    println()
end

println("Total unique stabilizer transformations (with signs): ", length(unique_tableaus))