using LinearAlgebra
using JuMP, MosekTools
using TensorOperations

function partial_trace(A::AbstractMatrix, dims::Tuple, trace_out::Tuple)
    reshaped = reshape(A, dims)
    
    show(stdout, "text/plain", reshaped)

    traced = zeros(2, 2)

    for j=1:2
        for k=1:2
            for i=1:2
                traced[j,k] += reshaped[i,j,i,k]
            end
        end
    end

    show(stdout, "text/plain", traced)
    
    return traced
end 

function trace_two_out_of_four(A::AbstractMatrix)
    tensor = A
    reduced = zeros(GenericAffExpr{ComplexF64,VariableRef}, 4, 4)
    for i1 in 1:2, i4 in 1:2, j1 in 1:2, j4 in 1:2
        idx_reduced = (i1 - 1) * 2 + i4
        jdx_reduced = (j1 - 1) * 2 + j4

        sum_expr = 0
        for a in 1:2, b in 1:2  # qubits 2 and 3
            # compute flat indices in original 4-qubit system (16×16)
            i_full = (i1 - 1) * 8 + (a - 1) * 4 + (b - 1) * 2 + i4
            j_full = (j1 - 1) * 8 + (a - 1) * 4 + (b - 1) * 2 + j4

            sum_expr += tensor[i_full, j_full]
        end

        reduced[idx_reduced, jdx_reduced] = sum_expr
    end
    # @tensor reduced[i1, i4, j1, j4] := tensor[i1, a, b, i4, j1, a, b, j4]
   
    return reduced
end

function trace_second_out_of_three(A::AbstractMatrix)
    tensor = reshape(A, 2, 2, 2, 2, 2, 2)
    @tensor reduced[i1, i3, j1, j3] := tensor[i1, a, i3, j1, a, j3]
    matrix = reshape(reduced, 4, 4)
    return matrix
end



function partial_transpose_second(A::AbstractMatrix, dims::Tuple)
    reshaped = reshape(A, dims)
    permuted = permutedims(reshaped, (1, 4, 3, 2))
    return reshape(permuted, 4, 4)
end

# Pauli matrices
sigma_x = [0 1; 1 0]
sigma_y = [0 -im; im 0]
sigma_z = [1 0; 0 -1]
paulis = [sigma_x, sigma_y, sigma_z]

# Function to generate a random SU(2) matrix
function random_su2()
    # Generate random real coefficients for Pauli matrices
    θ = 2π * rand()  # Random angle
    n = randn(3)     # Random direction
    n /= norm(n)     # Normalize to get unit vector

    # Construct the unitary matrix U = exp(-i θ n ⋅ σ / 2)
    H = -im * θ / 2 * (n[1] * sigma_x + n[2] * sigma_y + n[3] * sigma_z)
    U = exp(H)

    return U
end

# Generate a random SU(2) matrix
U = random_su2()
# U = [0 1; -1 0]

# Verify properties
print("U = \n")
show(stdout, "text/plain", U)
println()
# println("Determinant: ", det(U))  # Should be 1
# println("U * U† ≈ I: ", U * U' ≈ I)  # Should be identity
# println(conj(U))
# test = sigma_y * U * sigma_y
# println(test ≈ conj(U))

d = 2

phi_plus = zeros(ComplexF64, d^2)
phi_plus[1] = 1 / sqrt(d)
phi_plus[d^2] = 1 / sqrt(d)
println(phi_plus)

density_matrix_phi_plus = phi_plus * phi_plus'
# println(density_matrix_phi_plus)
# show(stdout, "text/plain", density_matrix_phi_plus)
# println()

choi_U = d * kron(I(d), U) * density_matrix_phi_plus * kron(I(d), U')
print("choi_U = \n")
show(stdout, "text/plain", choi_U)
println()

func_U=conj(U);
choi_func_U = d * kron(I(d), func_U) * density_matrix_phi_plus * kron(I(d), func_U')
print("choi_func_U = \n")
show(stdout, "text/plain", choi_func_U)
println()

model = Model(Mosek.Optimizer)
@variable(model, S[1:d^4, 1:d^4], PSD)
@variable(model, p >= 0)  

transformed_matrix = S * kron(I(d), transpose(choi_U), I(d))
println(length(S))
println(length(transformed_matrix))

channel = trace_two_out_of_four(transformed_matrix)
print(length(channel))

@constraint(model, channel == p * choi_func_U)
@objective(model, Max, p)
optimize!(model)

println("Optimal success probability p: ", value(p))
println("Optimized Choi matrix of S:\n", value.(S))


# A = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]

# traced_out_A = partial_trace(A, (2, 2, 2, 2), tuple(1, 4))


circuit = sigma_y * U * sigma_y
choi_circuit = d * kron(I(d), circuit) * density_matrix_phi_plus * kron(I(d), circuit')
show(stdout, "text/plain", choi_circuit)
println()

println(choi_func_U ≈ choi_circuit)

choi_y = d * kron(I(d), sigma_y) * density_matrix_phi_plus * kron(I(d), sigma_y')
partial_T_y = partial_transpose_second(choi_y, (2, 2, 2, 2))

first = kron(I(d), partial_T_y) * kron(choi_U, I(d))
first = trace_second_out_of_three(first)
yu = U * sigma_y
# choi_yu = d * kron(I(d), yu) * density_matrix_phi_plus * kron(I(d), yu')
partial_T_first = partial_transpose_second(first, (2, 2, 2, 2))
second = kron(I(d), partial_T_first) * kron(choi_y, I(d))
second = trace_second_out_of_three(second)
show(stdout, "text/plain", second)
println()
println(choi_circuit ≈ second)