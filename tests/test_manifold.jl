"""
testing manifold
"""


using LinearAlgebra
using DifferentialEquations
using BenchmarkTools
using Plots

include("../R3BP/src/R3BP.jl")

println("Running test for manifold!")

## define parameters
paramEM = R3BP.get_cr3bp_param(399, 301)
mu = paramEM.mu
println("mu: $mu")
lp = R3BP.lagrangePoints(mu)

# initial condition of halo
guess_halo = R3BP.halo_analytical_construct(mu, 2, 4500, paramEM.lstar, 1)
res_iter = R3BP.ssdc_periodic_xzplane([mu,], guess_halo.x0, guess_halo.period;
        fix="z", method=Tsit5(), reltol=1e-12, abstol=1e-12);

X0 = res_iter.x0 # [1.176924090973164, 0.0, -0.060210863312217, 0.0, -0.173836346247689, 0.0]
T = res_iter.period #3.385326412831325

## Setup for manifold call
function condition(u,t,integrator)
  u[2] - 0.3   # when y-value hits xx
end
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)

# parameters for manifolds
stability = true;
tf = -15.0
N = 50
direction = "positive"

## generate manifolds
sim = R3BP.get_manifold(mu, X0, T, tf, stability, N, direction, cb)

## create plot
plot_manif = plot(;linealpha=0.4, vars=(1,2), flip=false, aspect_ratio=:equal, size=(800,650),
     frame_style=:box, gridalpha=0.4, xlabel="x", ylabel="y")
for sol in sim
    plot!(plot_manif, Array(sol)[1,:], Array(sol)[2,:], label=false,
      linealpha=0.4, vars=(1,2), flip=false, c=:dodgerblue)
end
display(plot_manif)
println("Done!")
