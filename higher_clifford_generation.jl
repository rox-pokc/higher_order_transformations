using Combinatorics
using JSON

# Multi-qubit Pauli represented as (x_vec, z_vec, r)
struct PauliKQ
    x::Vector{Int}
    z::Vector{Int}
    r::Int  # 0: +1, 1: -1
end

Base.isequal(a::PauliKQ, b::PauliKQ) = a.x == b.x && a.z == b.z && a.r == b.r
Base.hash(p::PauliKQ, h::UInt) = hash((p.x, p.z, p.r), h)

function string_of(p::PauliKQ)
    rstr = p.r == 0 ? "+" : "-"
    opstr = ""
    for (xi, zi) in zip(p.x, p.z)
        op = xi == 0 && zi == 0 ? "I" :
             xi == 1 && zi == 0 ? "X" :
             xi == 0 && zi == 1 ? "Z" :
             xi == 1 && zi == 1 ? "Y" : "?"
        opstr *= op
    end
    return rstr * opstr
end

# Conjugate by 1-qubit gate on qubit q
function conjugate_1q(p::PauliKQ, gate::String, q::Int)
    x, z, r = copy(p.x), copy(p.z), p.r
    if gate == "H"
        r = (r + x[q]*z[q]) % 2
        x[q], z[q] = z[q], x[q]
    elseif gate == "S"
        r = (r + x[q]*z[q]) % 2
        z[q] = (z[q] + x[q]) % 2
    elseif gate == "I"
        return PauliKQ(x, z, r)
    else
        error("Unknown 1q gate $gate")
    end
    return PauliKQ(x, z, r)
end

# Conjugate by CNOT q1→q2
function conjugate_CNOT(p::PauliKQ, q1::Int, q2::Int)
    x, z, r = copy(p.x), copy(p.z), p.r
    r = (r + x[q1]*z[q2] * (1 + x[q2] + z[q1])) % 2
    x[q2] = (x[q2] + x[q1]) % 2
    z[q1] = (z[q1] + z[q2]) % 2
    return PauliKQ(x, z, r)
end

# Apply gate sequence to all basis Paulis: X1,...,Xk and Z1,...,Zk
function act_on_generators(seq, k::Int)
    paulis = [PauliKQ([i == j ? 1 : 0 for j in 1:k], zeros(Int, k), 0) for i in 1:k]  # Xi
    paulis_append = [PauliKQ(zeros(Int, k), [i == j ? 1 : 0 for j in 1:k], 0) for i in 1:k]  # Zi
    generators = vcat(paulis, paulis_append)

    for gate in seq
        if gate[1] == "CNOT"
            for i in eachindex(generators)
                generators[i] = conjugate_CNOT(generators[i], gate[2], gate[3])
            end
        else
            for i in eachindex(generators)
                generators[i] = conjugate_1q(generators[i], gate[1], gate[2])
            end
        end
    end
    return generators
end

# Generate all Clifford gates for k qubits using depth of gates
function generate_cliffords_kq(depth::Int, k::Int)
    gates_1q = ["I", "H", "S"]

    # Generate all possible 1-qubit gate applications to k qubits
    function one_qubit_steps()
        gate_sets = [gates_1q for _ in 1:k]
        return Iterators.product(gate_sets...)
    end

    # Generate all CNOT gates between different pairs of qubits
    function cnot_steps()
        cnots = []
        for ctrl in 1:k
            for tgt in 1:k
                if ctrl != tgt
                    push!(cnots, [("CNOT", ctrl, tgt)])
                end
            end
        end
        return cnots
    end

    # Combine all instructions (1-qubit + 1 CNOT) for each depth layer
    function all_layers()
        layers = []
        for step in one_qubit_steps()
            layer = [(g, i) for (i, g) in enumerate(step)]
            append!(layers, [layer])
        end
        append!(layers, cnot_steps())
        return layers
    end

    # Build all possible sequences of instructions up to the given depth
    layer_choices = all_layers()
    gate_sequences = Iterators.product(ntuple(_ -> layer_choices, depth)...) 

    seen = Dict{Vector{PauliKQ}, Vector{Any}}()
    for seq_tuple in gate_sequences
        seq = reduce(vcat, seq_tuple)
        generators = act_on_generators(seq, k)
        if !haskey(seen, generators)
            seen[generators] = seq
        end
    end

    return seen
end

function save_cliffords_to_file(cliffords, filename)
    json_array = [
        val 
        for (key, val) in cliffords
    ]

    open(filename, "w") do io
        JSON.print(io, json_array)
    end
end


# Set number of qubits and gate depth
k = 1          # number of qubits
depth = 6     # number of gate layers

# Generate the Clifford dictionary
# cliffords = generate_cliffords_kq(depth, k)

elapsed = @elapsed cliffords = generate_cliffords_kq(depth, k)
save_cliffords_to_file(cliffords, "cliffords_generation_outputs/cliffords_depth=$(depth)_#qubits=$(k).json")
println("Elapsed time: $(round(elapsed, digits=3)) seconds")

# for (generators, sequence) in cliffords
#     println("Gate sequence: ", sequence)
#     println("Pauli generator action:")
#     for (i, gen) in enumerate(generators)
#         label = i <= k ? "X[$i]" : "Z[$(i - k)]"
#         println("  $label → ", string_of(gen))
#     end
#     println("-"^40)
# end

println("Number of distinct Clifford operations generated: ", length(cliffords))