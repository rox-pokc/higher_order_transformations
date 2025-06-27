using JuMP
using MosekTools
using LinearAlgebra
using Random

sigma_x = [0 1; 1 0]
sigma_y = [0 -im; im 0]
sigma_z = [1 0; 0 -1]

function get_choi(A::AbstractMatrix, d::Int)
    # |phi+> = 1/sqrt(d)sum(|ii>)i=0..d-1
    phi_plus = zeros(ComplexF64, d^2)
    for i in 1:(d+1):d^2
        phi_plus[i] = 1 / sqrt(d)
    end

    # |phi+><phi+|
    density_matrix_phi_plus = phi_plus * phi_plus'

    return d * kron(I(d), A) * density_matrix_phi_plus * kron(I(d), A')
end

function partial_trace(A::AbstractMatrix, trace_out_dim::Vector{Int}; total_dim::Int = round(Int, log2(size(A, 1))))
    @assert size(A, 1) == size(A, 2) "Matrix must be square"
    @assert 2^total_dim == size(A, 1) "Matrix size doesn't match number of qubits"

    keep_dim = setdiff(1:total_dim, trace_out_dim)
    n_keep = length(keep_dim)

    dim_keep = 2^n_keep
    reduced = zeros(eltype(A), dim_keep, dim_keep)

    for i in 0:(dim_keep - 1), j in 0:(dim_keep - 1)
        sum_expr = zero(eltype(A))

        for trace_vals in Iterators.product(fill(0:1, length(trace_out_dim))...)
            i_bits = zeros(Int, total_dim)
            j_bits = zeros(Int, total_dim)

            # Set bits for kept dimensions
            i_keep_bits = reverse(digits(i, base=2, pad=n_keep))
            j_keep_bits = reverse(digits(j, base=2, pad=n_keep))

            for (bitpos, val) in zip(keep_dim, i_keep_bits)
                i_bits[bitpos] = val
            end
            for (bitpos, val) in zip(keep_dim, j_keep_bits)
                j_bits[bitpos] = val
            end

            # Set bits for traced-out dimensions
            for (bitpos, val) in zip(trace_out_dim, trace_vals)
                i_bits[bitpos] = val
                j_bits[bitpos] = val
            end

            i_idx = foldl((acc, b) -> 2*acc + b, i_bits, init=0) + 1
            j_idx = foldl((acc, b) -> 2*acc + b, j_bits, init=0) + 1

            sum_expr += A[i_idx, j_idx]
        end

        reduced[i+1, j+1] = sum_expr
    end

    return reduced
end

function partial_transpose_second(A::AbstractMatrix, dims::Tuple)
    reshaped = reshape(A, dims)
    permuted = permutedims(reshaped, (1, 4, 3, 2))
    return reshape(permuted, 4, 4)
end

function threshold_zero(A::AbstractArray, ϵ=1e-5)
    return map(x -> abs(x) < ϵ ? 0.0 : x, A)
end

function diagonalize(A::AbstractMatrix)
    # eigen(Hermitian(Matrix(A))) if A is sparse?...
    E = eigen(A)    
    D = Diagonal(E.values)    
    return D
end

function random_unitary(dim::Int; real::Bool=false) 
    gin = real ? randn(dim, dim) : randn(dim, dim) + im * randn(dim, dim)
    Q, R = qr(gin)
    d = diag(R)
    phases = d ./ abs.(d)
    phases[isnan.(phases)] .= 1
    Q = Q * Diagonal(phases)

    detQ = det(Q)
    correction = (detQ + 0im)^(-1/dim)
    Q = Q * correction

    return Q
end

function tensor_power(U::AbstractMatrix, k::Int)
    result = U
    for _ in 2:k
        result = kron(result, U)
    end
    return result
end

function unitary_linear_independent(list_U::Vector{Matrix{ComplexF64}}, d::Int)
    N = length(list_U)
    A = Matrix{ComplexF64}(undef, d^4, N)
    for i in 1:N
        choi = get_choi(list_U[i], d)
        A[:, i] = vec(choi)
    end
    return rank(A, atol=1e-8)
end

function find_dimension(d::Int, k::Int)
    n = 1000
    nLI = 0
    while true
        list_U_k = [tensor_power(random_unitary(d), k) for _ in 1:n]
        nLI_new = unitary_linear_independent(list_U_k, d^k)
        if nLI_new < n
            return nLI_new
        else
            n += 100
        end
    end
end

function minimal_stinespring_from_choi(C::AbstractMatrix)
    # Eigendecomposition of Choi matrix
    vals, vecs = eigen(Hermitian(C))
    # Filter nonzero eigenvalues
    keep = findall(>(1e-10), vals)
    eigen_vals = vals[keep]
    eigen_vecs = vecs[:, keep]

    r = length(eigen_vals)  # rank of the Choi matrix
    d = Int(sqrt(size(C, 1)))  # assume C ∈ L(Ho ⊗ Hi), dim(Ho) * dim(Hi) = size(C, 1)
    d_i = d_o = Int(d)   # assume square dims

    # Build the purification |ψ⟩ = ∑ √λ_k |e_k⟩ ⊗ |k⟩
    pure_vecs = zeros(ComplexF64, d_i, d_o, r)
    for (k, λ) in enumerate(eigen_vals)
        e_k = eigen_vecs[:, k]
        e_k_tensor = reshape(e_k, d_o, d_i)
        pure_vecs[:, :, k] .= sqrt(λ) * permutedims(e_k_tensor, (2, 1))  # reorder to H_i ⊗ H_o
    end

    # Reshape into a matrix: V: H_i → H_o ⊗ H_a
    V = reshape(permutedims(pure_vecs, (2, 3, 1)), d_o * r, d_i)


    C_reconstructed = zeros(ComplexF64, d_o * d_i, d_o * d_i)
    for k in 1:r
        v = vec(pure_vecs[:, :, k])
        C_reconstructed += v * v'  # outer product
    end

    println(isapprox(C, C_reconstructed; atol=1e-10))


    Kraus_ops = [pure_vecs[:, :, k] / sqrt(eigen_vals[k]) for k in 1:r]

    return pure_vecs, V, Kraus_ops
end

function svd_decomposition(V::AbstractMatrix)
    svd_result = svd(V)
    return (U = svd_result.U, S = svd_result.S, Vt = svd_result.Vt)
end



d = 2  # 1 qubit
k = 1 # 1 copy

dim = find_dimension(d, k)
println("dim = ", dim)

list_U = Vector{Matrix{ComplexF64}}(undef, dim)
list_U_k = Vector{Matrix{ComplexF64}}(undef, dim)
Jin = Vector{Matrix{ComplexF64}}(undef, dim)
Jout = Vector{Matrix{ComplexF64}}(undef, dim)

for i in 1:dim
    U = random_unitary(d, real=false)
    list_U[i] = U
    list_U_k[i] = tensor_power(U, k)

    # Choi matrix for Φ_U
    Jin[i] = get_choi(U, d)

    # Define target U
    func_U = conj(U)
    Jout[i] = get_choi(func_U, d)
end

function solve_sample()
    model = Model(Mosek.Optimizer)
    # set_optimizer_attribute(model, MOI.Silent(), true)

    @variable(model, S[1:d^4, 1:d^4], PSD)
    @variable(model, F[1:d^4, 1:d^4], PSD)
    @variable(model, p >= 0)

    # C = S + F
    @expression(model, C[i=1:d^4, j=1:d^4], S[i, j] + F[i, j])
    
    for i in 1:dim
        @constraint(model, partial_trace(S * kron(I(d), transpose(Jin[i]), I(d)), [2, 3]) == p * Jout[i])
    end

    # First constraint: Tr_4(C) == kron(Tr_{3,4}(C), I/d)
    PT_4 = partial_trace(C, [4])          # shape: d^3 × d^3
    PT_34 = partial_trace(C, [3, 4])      # shape: d^2 × d^2
    kron_target = kron(PT_34, I(d) / d)    # shape: d^3 × d^3

    @constraint(model, PT_4 == kron_target)

    # Second constraint: Tr_{2,3,4}(C) == d * I_d
    PT_234 = partial_trace(C, [2, 3, 4])  # shape: d × d
    @constraint(model, PT_234 == I(d) * d)
    
    @objective(model, Max, p)

    optimize!(model)
    println(termination_status(model))
    println(primal_status(model))
    println(dual_status(model))

    return value.(S), value.(C), value(p)
end

S_opt, C_opt, p_star = solve_sample()
println("Optimal p = ", p_star)
println("Optimal S = ")
show(stdout, "text/plain", threshold_zero(S_opt))
println()

println("Optimal C = ")
show(stdout, "text/plain", threshold_zero(C_opt))
println()

println("Trace C_opt: ", tr(C_opt))

evals_opt = eigen(S_opt).values
rank_S_opt = count(abs.(evals_opt) .> 1e-6)
println("Optimal rank = $rank_S_opt")

rank_C_opt = count(abs.(eigen(C_opt).values) .> 1e-6)
println("C_opt rank = $rank_C_opt")

println("Optimal F = ")
show(stdout, "text/plain", threshold_zero(C_opt - S_opt))
println()

## Checking on random unitary 

U_rand = random_unitary(d, real=false)
map_rand = C_opt * kron(I(d), transpose(get_choi(U_rand, d)), I(d))
show(stdout, "text/plain", threshold_zero(partial_trace(map_rand, [2, 3])))
println()
show(stdout, "text/plain", threshold_zero(get_choi(conj(U_rand), d)))
println()

sigma_y_choi = get_choi(sigma_y, d)
circuit = sigma_y * U_rand * sigma_y
choi_circuit = get_choi(circuit, d)
println(choi_circuit ≈ get_choi(conj(U_rand), d))

instrument = kron(sigma_y_choi, I(d^2)) * kron(I(d^2), transpose(sigma_y_choi))

show(stdout, "text/plain", threshold_zero(instrument))
println("Trace instrument: ", tr(instrument))
println(size(instrument))
println()
show(stdout, "text/plain", threshold_zero(partial_trace(instrument * kron(I(d), transpose(get_choi(U_rand, d)), I(d)), [2, 3])))
println()
show(stdout, "text/plain", get_choi(conj(U_rand), d))

println("1st: ", isapprox(partial_trace(instrument, [4]), kron(partial_trace(instrument, [3, 4]),  I(d) / d)))
println("2nd: ", isapprox(partial_trace(instrument, [2, 3, 4]), I(d) * d))

pure_vecs, V, Kraus_ops = minimal_stinespring_from_choi(S_opt)

println("Purification vector shape: ", size(pure_vecs))   # (d_in, d_out, d_a)
show(stdout, "text/plain", threshold_zero(pure_vecs))
println()
println("Isometry shape: ", size(V))                # (d_out * d_a, d_in)
show(stdout, "text/plain", threshold_zero(V))
println()

println("Kraus: ", size(Kraus_ops))                
show(stdout, "text/plain", threshold_zero(Kraus_ops[1]))
println()

println("Choi of sigma_y: ", size(sigma_y_choi))
show(stdout, "text/plain", threshold_zero(sigma_y_choi))
println()

Q, R = qr(V)  
println("Q of V: ")
show(stdout, "text/plain", threshold_zero(Matrix(Q)))
println()

println("R of V: ")
show(stdout, "text/plain", threshold_zero(R))
println()

H = (sigma_x + sigma_z) / sqrt(2)

components = svd_decomposition(V)

# SVD decomposition
U = components.U
S = components.S
Vt = components.Vt

println("U: ")
show(stdout, "text/plain", threshold_zero(U))
println()

println("S: ")
show(stdout, "text/plain", threshold_zero(S))
println()

println("Vt: ")
show(stdout, "text/plain", threshold_zero(Vt))
println()


# instrument

S1 = kron(I(d), I(d), sigma_x, sigma_x)
println(isapprox(S1 * C_opt * transpose(S1), C_opt))

S2 = kron(I(d), I(d), sigma_z, sigma_z)
println(isapprox(S2 * C_opt * transpose(S2), C_opt))

S3 = kron(sigma_x, sigma_x, I(d), I(d))
println(isapprox(S3 * C_opt * transpose(S3), C_opt))

S4 = kron(sigma_z, sigma_z, I(d), I(d))
println(isapprox(S4 * C_opt * transpose(S4), C_opt))


H = (1 / sqrt(2)) * [1  1; 1 -1]
println(det(H))
choi_H_func = get_choi(conj(H), d)
show(stdout, "text/plain", threshold_zero(choi_H_func))
println()

stab1 = kron(sigma_x, sigma_z)
println(isapprox(stab1 * choi_H_func * transpose(stab1), choi_H_func))

stab2 = kron(sigma_z, sigma_x)
println(isapprox(stab2 * choi_H_func * transpose(stab2), choi_H_func))

S = [1 0; 0 im]
println(det(S))

choi_S_func = get_choi(conj(S), d)
show(stdout, "text/plain", threshold_zero(choi_S_func))
println()

stab3 = kron(sigma_y, sigma_x)
println(isapprox(stab3 * choi_S_func * transpose(stab3), (-1) * choi_S_func))
show(stdout, "text/plain", threshold_zero(stab3 * choi_S_func * transpose(stab3)))
println()

stab4 = kron(sigma_z, sigma_z)
println(isapprox(stab4 * choi_S_func * transpose(stab4), choi_S_func))

show(stdout, "text/plain", threshold_zero(sigma_y * (im * H) * sigma_y))
println()

show(stdout, "text/plain", conj(im * H))
println()

println(det(im * H))