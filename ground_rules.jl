using LinearAlgebra

d = 2

sigma_x = [0 1; 1 0]
sigma_y = [0 -im; im 0]
sigma_z = [1 0; 0 -1]
H = (1 / sqrt(2)) * [1  1; 1 -1]
S = [1 0; 0 im]

println(H * sigma_x * adjoint(H))
println(H * sigma_y * adjoint(H))
println(H * sigma_z * adjoint(H))

println(S * sigma_x * adjoint(S))
println(S * sigma_y * adjoint(S))
println(S * sigma_z * adjoint(S))

println(sigma_y * sigma_x * adjoint(sigma_y))
println(sigma_y * sigma_y * adjoint(sigma_y))
println(sigma_y * sigma_z * adjoint(sigma_y))

# ["S", "S", "H", "S", "I", "I"], ["S", "S", "H", "S", "I", "I"]) 
# println(I(d) * I(d) * S * H * S * S)
# res = S * S * H * S * I(d) * I(d)
# res = I(d) * I(d) * S * H * S * S
# println(S * S * H * S * I(d) * I(d))

# ["S", "H", "S", "H", "I", "I"], ["S", "H", "S", "H", "I", "I"]
# "S", "S", "H", "S", "S", "H"
enc = H * S * S * H * S * S
dec = H * S * S * H * S * S

show(stdout, "text/plain", enc)
println()
show(stdout, "text/plain", dec)
println()

println("H and S shitshow")

show(stdout, "text/plain", conj(H))
println()
show(stdout, "text/plain", dec * (H) * enc)
println()

show(stdout, "text/plain", conj(S))
println()
show(stdout, "text/plain", dec * (S) * enc)
println()

show(stdout, "text/plain", conj(H * S))
println()
show(stdout, "text/plain", dec * (H * S) * enc)
println()

println(isapprox(conj(H), dec * H * enc))
println(isapprox((-im) * conj(S), dec * S * enc))
println(isapprox((im) * conj(H * S), dec * (H * S) * enc))

println(isapprox(sigma_y * (im) * H * sigma_y, conj((im) * H)))
println(isapprox(sigma_y * (-im) * S * sigma_y, conj((-im) * S)))
println((-im) * sigma_y * H * (-im) * sigma_y)
println(((-im) * sigma_y * S * (-im) * sigma_y)/(-im))

println("transpose & inverse")
U2 = H * S * S * H * S
# U2 = S * H * S * S * H * I
println(U2)
println(inv(U2))
println(S * H * U2 * S * H)
println(isapprox(inv(U2), S * H * U2 * S * H))
# println((-im) * sigma_y * H * S * H * S * S * S * U2 * S * H * S * S * (-im) * sigma_y)
println(transpose(U2))
# enc: ["S", "S", "S", "H", "S", "H"], dec: ["S", "S", "H", "S", "I", "I"] - transposition
# enc: ["H", "S", "I", "I", "I", "I"], dec: ["S", "H", "S", "H", "I", "I"] 
# dec = S * H * S * S
dec = H * S * H * S
# enc = H * S * H * S * S * S
enc = S * H
println(dec * U2 * enc)

function equal_up_to_global_phase(A::AbstractMatrix, B::AbstractMatrix; atol=1e-8)
    # Flatten and find first nonzero element in B to estimate phase
    idx = findfirst(!≈(0), B)
    λ = A[idx] / B[idx]  # estimated global phase
    println(λ)
    return isapprox(A, λ * B; atol=atol)
end

println(equal_up_to_global_phase(inv(U2), (-im) * sigma_y * dec * U2 * enc * (-im) * sigma_y))