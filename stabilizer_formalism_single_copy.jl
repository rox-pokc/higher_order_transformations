using QuantumClifford

# Initialize сircuit
stab = S"XXIIIIII
        ZZIIIIII
        IIXXIIII
        IIZZIIII
        IIIIXXII
        IIIIZZII
        IIIIIIXX
        IIIIIIZZ"

# Encode
enc = sY
apply!(stab, enc(5))
println(stab)

# Bell measure 3rd and 5th qubits
x3x5 = project!(stab, P"IIXIXIII")[1]
z3z5 = project!(x3x5, P"IIZIZIII")[1]
println(z3z5)

# Decode
dec = sY
apply!(z3z5, dec(1))

# Reduce stabilizers
canonicalize!(z3z5)
println("Reduced quantum comb:\n", z3z5)

# Inject U - Clifford
U = sPhase
apply!(z3z5, U(7))

# Bell measure 4th and 8th qubits
x4x8 = project!(z3z5, P"IIIXIIIX")[1]
z4z8 = project!(x4x8, P"IIIZIIIZ")[1]

# Bell measure 2nd and 7th qubits
x2x7 = project!(z4z8, P"IXIIIIXI")[1]
z2z7 = project!(x2x7, P"IZIIIIZI")[1]

println("Full circuit:\n", z2z7)

canonicalize!(z2z7)
println("Reduced full circuit:\n", z2z7)