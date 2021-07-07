"""
Test for multiple shooting
"""


using LinearAlgebra
using DifferentialEquations
using BenchmarkTools
using Plots

pyplot()

include("../R3BP/src/R3BP.jl")

println("Running test for multiple shooting!")

## Initialize ODE settings
# ODE settings throughout this notebook
reltol = 1.e-13
abstol = 1.e-13
method = Tsit5()

## construct halo initial guess
CR3BP_param = R3BP.get_cr3bp_param(10, 299)
lp = 2
Az_km = 260*1e3
println("Halo guess Az_km: $Az_km")
northsouth = 3   # 1 or 3
guess0 = R3BP.halo_analytical_construct(CR3BP_param.mu, lp, Az_km, CR3BP_param.lstar, northsouth)


## prepare initial guess for multiple shooting
# number of data in initial guess array
n_data, _ = size(guess0.fullstates)
# number of nodes to use
n = 6
# get index to extract from 3rd order initial guess state-history
skip_idx = Int(floor(500/(n-1)))
extract_idx = []
for i = 1:skip_idx:n_data
    push!(extract_idx, i)
end
# append final node for periodic orbit
push!(extract_idx, 1)

# get x0s, states along the trajectory
x0s = []
for i = 1:n
    extract_idx[i]
    push!(x0s, guess0.fullstates[extract_idx[i], 1:6])
end

# get tofs, propagation of the first ~ (n-1)th state
tofs = zeros(n-1)
for i = 1:n-1
    if i != n-1
        tof = guess0.fullstates[extract_idx[i+1], 7] - guess0.fullstates[extract_idx[i], 7]
    else
        tof = guess0.period - guess0.fullstates[extract_idx[i], 7]
    end
    tofs[i] = tof
end

## create ODE Problem
# construct initial ODE problem
μ = CR3BP_param.mu
prob_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!,
    vcat(guess0.x0, reshape(I(6), (36,)))[:], guess0.period, (μ));
# propagate initial guess for plotting later
sol = solve(prob_stm, method, reltol=reltol, abstol=abstol)
solig = R3BP.sol_to_arrays(sol);   # obtain state-history as array for plotting


## Run multiple shooting
tolDC = 1.e-13   # tolerance to be achieved by multiple-shooting algorithm
x0_conv, tof_conv = R3BP.multiple_shooting(prob_stm, x0s, tofs, tolDC;
    periodic=true, reltol=reltol, abstol=abstol, use_ensemble=true);

# propagate converged result of multiple-shooting
prob = ODEProblem(R3BP.rhs_cr3bp_sv!, x0_conv[1:6], sum(tof_conv), (μ))
sol = solve(prob, Tsit5(), reltol=reltol, abstol=abstol);
solmsdc = R3BP.sol_to_arrays(sol.u);   # obtain state-history as array for plotting

## Benchmarking without / with ensemble
println("\n*** without ensemble ***")
@benchmark x0_conv, tof_conv = R3BP.multiple_shooting(prob_stm, x0s, tofs, tolDC;
    periodic=true, reltol=reltol, abstol=abstol, use_ensemble=false);

println("\n*** with ensemble ***")
@benchmark x0_conv, tof_conv = R3BP.multiple_shooting(prob_stm, x0s, tofs, tolDC;
    periodic=true, reltol=reltol, abstol=abstol, use_ensemble=true);


## create plot of trajectory
# idx1, idx2 = 1, 2
# ptraj = plot(flip=false, aspect_ratio=:equal, size=(800,650),
#              xlabel="state $idx1", ylabel="state $idx2", gridalpha=0.4, frame_style=:box)
#
# plot!(ptraj, solig[idx1,:], solig[idx2,:], label="Initial guess propagated")
# plot!(ptraj, solmsdc[idx1,:], solmsdc[idx2,:], label="Multiple-shooting")
#
# # plot nodes used by multiple-shooting algorithm
# for i = 1:n
#     scatter!(ptraj, [x0_conv[(i-1)*6+1:6i][idx1]], [x0_conv[(i-1)*6+1:6i][idx2]],
#         label="node $i", marker=(:cross, 3.5, 3.5))
# end
# display(ptraj)

println("Done!")
