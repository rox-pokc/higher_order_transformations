using Combinatorics

# include("src/Config.jl")
# using .Config  

# Represent a signed Pauli operator
struct Pauli1Q
    x::Int  # 0 or 1
    z::Int  # 0 or 1
    r::Int  # 0: +1, 1: -1
end

Base.isequal(a::Pauli1Q, b::Pauli1Q) = a.x == b.x && a.z == b.z && a.r == b.r
Base.hash(p::Pauli1Q, h::UInt) = hash((p.x, p.z, p.r), h)

# Pretty print
function string_of(p::Pauli1Q)
    s = p.r == 0 ? "+" : "-"
    op = p.x == 0 && p.z == 0 ? "I" :
         p.x == 1 && p.z == 0 ? "X" :
         p.x == 0 && p.z == 1 ? "Z" :
         p.x == 1 && p.z == 1 ? "Y" : "?"
    return s * op
end

# Conjugate Pauli by gate
function conjugate(p::Pauli1Q, gate::String)
    x, z, r = p.x, p.z, p.r
    if gate == "I"
        return Pauli1Q(x, z, r)
    elseif gate == "H"
        return Pauli1Q(z, x, (r + x*z) % 2)
    elseif gate == "S"
        return Pauli1Q(x, (z + x) % 2, (r + x*z) % 2)
    else
        error("Unknown gate $gate")
    end
end

# Track conjugation of X and Z under a sequence
function act_on_XZ(seq::Vector{String})
    px = Pauli1Q(1, 0, 0)  # X
    pz = Pauli1Q(0, 1, 0)  # Z
    for g in seq
        px = conjugate(px, g)
        pz = conjugate(pz, g)
    end
    return (px, pz)
end

# Generate unique Clifford actions up to global phase
function generate_cliffords(depth::Int)
    seen = Dict{Tuple{Pauli1Q, Pauli1Q}, Vector{String}}()
    gates = ["I", "H", "S"]

    for seq in Iterators.product(ntuple(_ -> gates, depth)...)
        px, pz = act_on_XZ(collect(seq))
        key = (px, pz)
        if !haskey(seen, key)
            seen[key] = collect(seq)
        end
    end
    return seen
end

# GENERATION EXAMPLE

# cliffs = generate_cliffords(6)

# global yes = 0
# global no = 0
# for ((px, pz), seq) in cliffs
#     println("X ↦ $(string_of(px)), Z ↦ $(string_of(pz))  ←  ", seq)
#     # println("X ↦ $(string_of(px)), Z ↦ $(string_of(pz))  ←  ", join(seq, " → "))
#     U_arr = [CLIFFORD_MATRICES[c] for c in seq]
#     U_combined = reduce(*, U_arr)
#     show(stdout, "text/plain", U_combined)
#     println()
#     if conj(U_combined) == U_combined
#         global yes += 1
#     else
#         global no += 1
#     end
# end
# println("YES: ", yes)
# println("NO: ", no)

# println("\nTotal unique Clifford operations: ", length(cliffs))