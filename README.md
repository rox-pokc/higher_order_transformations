# Higher-Order Quantum Transformations (Julia)
Clifford-structured and SDP-based search for quantum combs/superinstruments that implement functions of unknown unitaries (e.g., inversion, conjugation). This repo contains small, focused scripts rather than a monolithic package.

---

## What’s inside
- **SDP approach (single copy):** exact probabilistic superinstruments with no ancillas.
- **Stabilizer-symplectic approaches:** Clifford-restricted searches (single and two copies).
- **Clifford generators:** enumerate 1-qubit and multi-qubit Clifford groups to feed the searches.

---

## Scripts 

| Script | Method | Qubits (dim) | Copies | Aux | Target / Notes |
|---|---|---:|---:|---:|---|
| `sdp_single_copy.jl` | SDP (exact, probabilistic) | 2 | 1 | 0 | Solve for a heralded superinstrument; maximizes success under exactness constraints. |
| `stabilizer_formalism_single_copy.jl` | Stabilizer formalism | 2 | 1 | 0 | Clifford-induced combs for **conjugation**. |
| `clifford_generation.jl` | Enumeration | 1-qubit | — | — | Generate all single-qubit Clifford gates (sequences + matrices). |
| `higher_clifford_generation.jl` | Enumeration | many-qubit | — | — | Generate multi-qubit Cliffords (e.g., via symplectic tableaux). |
| `stabilizer_symplectic_search_single_copy.jl` |  Stabilizer formalism | 2 | 1 | 0 | Search for Clifford-induced combs. |
| `stabilizer_symplectic_search_multiple_copy_parallel.jl` |  Stabilizer formalism | 2 | 2 | 0 | Search for Clifford-induced combs for parallel two copies of unitary. |

---

## Requirements
- **Julia** ≥ 1.10 (tested on 1.11.x)
- Common packages:
  - `LinearAlgebra`, `SparseArrays`, `Random`, `Statistics`
  - `JuMP`, `MosekTools` 
  - `QuantumClifford`, `Combinatorics`
  - `JSON`

Install (example):
```julia
import Pkg
Pkg.activate(".")
Pkg.add([
  "LinearAlgebra","Statistics",
  "JuMP","MosekTools",
  "QuantumClifford","Combinatorics","JSON"
])
