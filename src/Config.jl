module Config

using QuantumClifford

export CLIFFORD_GATES, PAULI_GATES, CLIFFORD_MATRICES, X, Y, Z, I, S, H

X = [0 1; 1 0]
Y = [0 -im; im 0]
Z = [1 0; 0 -1]
I = [1 0; 0 1]
S = [1 0; 0 im]
H = (1 / sqrt(2)) * [1  1; 1 -1]
CNOT = [1 0 0 0;
        0 1 0 0;
        0 0 0 1;
        0 0 1 0]

const CLIFFORD_GATES = Dict(
    "H" => sHadamard,
    "S" => sPhase,
    "I" => sId1,
    "CNOT" => (ctrl::Int, trgt::Int) -> sCNOT(ctrl, trgt)
)

const PAULI_GATES = Dict(
    'I' => I,
    'X' => X,
    'Y' => Y,
    'Z' => Z
)

const CLIFFORD_MATRICES = Dict(
    "H" => H,
    "S" => S,
    "I" => I,
    "CNOT" => CNOT
)

end