using Printf

# EA in units of the Boltzmann constant k_B
EA   = 300
A(T) = 5 + 1 / (1 + 1e-2(T - EA)^2)

k(T) = A(T) * exp(-EA / T)

println("# Temperature(K)  Rate(1/s)")
for T in 250:10:350
    @printf "  %5.1f           %.5f\n" T k(T)
end

println()
println()
println("# Fitting curve")
println("arrhenius(T) = $(A(250)) * exp(-$EA / T)")
