"""
testing manifold
"""


using LinearAlgebra
using DifferentialEquations
using BenchmarkTools
using Plots

pyplot()

include("../R3BP/src/R3BP.jl")

println("Running test for manifold!")

## define parameters
params = R3BP.get_cr3bp_param(399, 301)
mu = params.mu
println("mu: $mu")
lp = R3BP.lagrangePoints(mu)

# initial condition of halo
X0 = [1.176924090973164, 0.0, -0.060210863312217, 0.0, -0.173836346247689, 0.0]
T = 3.385326412831325

## Setup for manifold call
function condition(u,t,integrator)
  u[1] - 1.25   # when y-value hits xx
end

affect!(integrator) = terminate!(integrator)

# assign callback
cb = ContinuousCallback(condition,affect!)

# parameters for manifolds
stability = true;
tf = -10.0

## generate manifolds
outsim, WarmStart= R3BP.get_manifold(mu, X0, T, tf, stability; lstar=params.lstar, callback=nothing, xdir="negative");

## create plot
display(plot(outsim, linealpha=0.4, vars=(1,2), flip=false, aspect_ratio=:equal, size=(800,650), c=:orangered,
     frame_style=:box, gridalpha=0.4, xlabel="x", ylabel="y"))
println("Done!")
