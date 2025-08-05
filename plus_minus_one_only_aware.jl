using LinearAlgebra
using Combinatorics

# Stabilizer generator: (x, z) for Pauli, s ∈ Z₂ for ± sign
struct Tableau
    x::Int
    z::Int
    s::Int  # 0 = +1, 1 = -1
end

# Pretty-print Pauli with ±
function Base.show(io::IO, t::Tableau)
    pauli = if t.x == 0 && t.z == 0
        "I"
    elseif t.x == 1 && t.z == 0
        "X"
    elseif t.x == 0 && t.z == 1
        "Z"
    elseif t.x == 1 && t.z == 1
        "Y"
    end
    sign = t.s == 0 ? "+" : "-"
    print(io, "$sign$pauli")
end

# Apply H gate
function apply_H(t::Tableau)
    x, z, s = t.x, t.z, t.s
    new_x, new_z = z, x
    new_s = (s + 2 * x * z) % 2
    return Tableau(new_x, new_z, new_s)
end

# Apply S gate
function apply_S(t::Tableau)
    x, z, s = t.x, t.z, t.s
    new_x, new_z = x, (z + x) % 2
    new_s = (s + x * z) % 2
    return Tableau(new_x, new_z, new_s)
end

# Apply sequence of gates to initial Z
function apply_sequence(seq::Vector{String})
    t = Tableau(0, 1, 0)  # Start with +Z
    for g in reverse(seq)
        if g == "H"
            t = apply_H(t)
        elseif g == "S"
            t = apply_S(t)
        elseif g == "I"
            continue
        else
            error("Unknown gate: $g")
        end
    end
    return t
end

# Check tableau equality
equal_t(t1::Tableau, t2::Tableau) = t1.x == t2.x && t1.z == t2.z && t1.s == t2.s

# Generate all distinct ±1 stabilizers under Clifford gates
function generate_unique_pm_stabilizers(depth::Int)
    gate_names = ["H", "S", "I"]
    seen = Tableau[]
    results = Dict{Tableau, Vector{String}}()

    for seq in Iterators.product(ntuple(_ -> gate_names, depth)...)
        gatelist = collect(seq)
        t = apply_sequence(gatelist)
        if !any(existing -> equal_t(existing, t), seen)
            push!(seen, t)
            results[t] = gatelist
        end
    end

    return results
end

# Run it
depth = 4
results = generate_unique_pm_stabilizers(depth)

for (t, seq) in results
    println("Sequence: ", join(seq, " → "), "  →  ", t)
end
