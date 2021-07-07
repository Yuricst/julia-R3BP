"""
Test julia script to be used on PACE
"""

using Distributed
#using DifferentialEquations

#include("../R3BP/src/R3BP.jl")

println("Testing on PACE...")

np = nprocs()
println("Number of procs ..... $np")


println("Done!")
