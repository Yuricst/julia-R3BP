"""
Test for multiple shooting
"""


using LinearAlgebra
using DifferentialEquations
using BenchmarkTools
using Plots

include("../R3BP/src/R3BP.jl")

println("Running test for multiple shooting!")

## Initialize ODE settings
# ODE settings throughout this notebook
reltol = 1.e-13
abstol = 1.e-14
method = Vern7()

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
p = (CR3BP_param.mu)
prob_stm = ODEProblem(R3BP.rhs_cr3bp_svstm!,
    vcat(guess0.x0, reshape(I(6), (36,)))[:], guess0.period, p,
    method=method, reltol=reltol, abstol=abstol
);
# propagate initial guess for plotting later
sol = solve(prob_stm)
solig = R3BP.sol_to_arrays(sol);   # obtain state-history as array for plotting


## Run all setting combinations of multiple shooting
tolDC = 1.e-12  # tolerance to be achieved by multiple-shooting algorithm
# x0_conv, tof_conv, convflag = R3BP.multiple_shooting(prob_stm, x0s, tofs, tolDC;
#     periodic=true, fix_time=false, fix_x0=false, fix_xf=false,
#     use_ensemble=false, rhs=R3BP.rhs_cr3bp_svstm!, p=p,
#     maxiter=10, reltol=reltol, abstol=abstol, method=method)

x0_conv, tof_conv, convflag = R3BP.multiple_shooting(prob_stm, x0s, tofs, tolDC;
    periodic=true, fix_time=false, fix_x0=false, fix_xf=false,
    use_ensemble=false, rhs=R3BP.rhs_cr3bp_svstm!, p=p,
    maxiter=8, reltol=reltol, abstol=abstol, method=method)

# propagate converged result of multiple-shooting
prob = ODEProblem(R3BP.rhs_cr3bp_sv!, x0_conv[1:6], sum(tof_conv), p, method=method, reltol=reltol, abstol=abstol)
sol = solve(prob)
solmsdc = R3BP.sol_to_arrays(sol.u);   # obtain state-history as array for plotting


## Benchmarking without / with ensemble
# println("\n*** Benchmarking without ensemble ***")
# b1 = @benchmark x0_conv, tof_conv = R3BP.multiple_shooting(prob_stm, x0s, tofs, tolDC;
#         periodic=true, reltol=reltol, abstol=abstol, use_ensemble=false, verbose=false)
# println(mean(b1))
#
# println("\n*** Benchmarking with ensemble ***")
# b2 = @benchmark x0_conv, tof_conv = R3BP.multiple_shooting(prob_stm, x0s, tofs, tolDC;
#     periodic=true, reltol=reltol, abstol=abstol, use_ensemble=true, verbose=false)
# println(mean(b2))

## create plot of trajectory
function get_plot(idx1, idx2, aspect_ratio=:equal)
    ptraj = plot(flip=false, aspect_ratio=aspect_ratio, size=(800,650),
                 xlabel="state $idx1", ylabel="state $idx2", gridalpha=0.4, frame_style=:box)

    plot!(ptraj, solig[idx1,:], solig[idx2,:], label="Initial guess propagated")
    plot!(ptraj, solmsdc[idx1,:], solmsdc[idx2,:], label="Multiple-shooting")

    # plot nodes used by multiple-shooting algorithm
    for i = 1:n
        scatter!(ptraj, [x0s[i][idx1]], [x0s[i][idx2]],
            label="Initial node $i", marker=(:circle, 3.5, 3.5))
        scatter!(ptraj, [x0_conv[(i-1)*6+1:6i][idx1]], [x0_conv[(i-1)*6+1:6i][idx2]],
            label="Final node $i", marker=(:cross, 3.5, 3.5))
    end
    return ptraj
end
ptraj1 = get_plot(1, 2, :none)
ptraj2 = get_plot(1, 3, :none)
ptraj3 = get_plot(2, 3, :none)
display(plot(ptraj1, ptraj2, ptraj3, layout = (1,3)))

println("Done!")
