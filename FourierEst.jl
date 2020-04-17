
#--------------------------------------------------
# Script for Fourier estimation method            #
# UMN, Pharmacometrics © 2020                     #
#--------------------------------------------------
# Author: Mutaz Jaber                             #
# GitHub: Mutaz94                                 #
#--------------------------------------------------
# different harmonics will be used and it depends #
# on the L2-approximation algorithm for the       #
# Adjustment                                      #
#--------------------------------------------------
# Equation:                                       #
# R(t)= a₀ + ∑ a_n *cos(2nπt/24)+b_n*sin(2nπt/24) #
# n: harmonic number, a,b: Fourier coeff          #
#-------------------------------------------------$

using LinearAlgebra, Optim, Plots
using CSV

# Importing some data to test
df = CSV.read("../data.csv")

time = convert(Array, df.TIME)
Conc = convert(Array, df.AVG)

Conc_norm = normalize(Conc) #to preserve the shape 

using LsqFit


 function fourEst!(t,p)
            a₀, a_1, b_1, a_2, b_2 = p
            a₀ .+ a_1*cos.((2*π.*t)/24) .+ b_1*sin.((2*π.*t)/24) +
                 a_2*cos.((2*π.*t)/12) .+ b_2*sin.((2*π.*t)/12)
end

init = [0.1, 0.2, 0.1, 1, 1]

fit_norm = curve_fit(fourEst!, time, Conc_norm, init)
fit2 = curve_fit(fourEst!, time, Conc, init) #this is not normalized

est = fit_norm.param
res = fit_norm.residual
pred = fourEst!(time, est)

using Plots; plotly()
plot(time, pred)
plot!(time, Conc_norm)
scatter(time, res)

# L2 approximation function to check the percentage
# of each harmonic (Contribution).

function L2_approx(est)
    a₀, a₁, b₁, a₂, b₂ = est
    R₀ = a₀
    R₁ = sqrt((a₁)^2 + (b₁)^2)
    R₂ = sqrt((a₂)^2 + (b₂)^2)
    Rn = R₀ + R₁ + R₂

    H₀ = (2 *(R₀)^2/(2 *(R₀)^2 + Rn^2))*100
    H₁ = ((R₁)^2/(2 *(R₀)^2 + Rn^2))*100
    H₂ = ((R₂)^2/(2 *(R₀)^2 + Rn^2))*100
    println("p₀ = $H₀", "\n p₁ = $H₁", "\n p₂ = $H₂")
    return  nothing
end

percentage = L2_approx(est)

