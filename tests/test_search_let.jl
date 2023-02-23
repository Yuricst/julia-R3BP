"""
test for perilune targeting in CR3BP
"""

using LinearAlgebra
using Roots
using Plots
using DifferentialEquations
using Dierckx
using Suppressor: @suppress_err
using ProgressMeter
using Printf

using AstrodynamicsBase

plotly()

include("../R3BP/src/R3BP.jl")

params = R3BP.get_bcr4bp_param(399, 301)

# periodic orbit
x0 = [  0.9828163970888667
 -0.0002378763303112657
  6.2209559202675e-5
  0.10193123254290165
 -2.1469314506800314
  0.037255023498971654];
tof = -150*86400/params.tstar
radius_target = 6578/params.lstar

# create LET
lets, pt0, cbs = R3BP.grid_search_let(
    params, x0, tof, radius_target, 40, 1e-6, true, true
);

# propagate trajectory
sims_let = []
@showprogress for traj in lets
    rp1, t0, _tof = traj
    @printf("r_perigee = %4.4f, t0 = %3.4f\n", rp1*params.lstar, t0*180/π)
    # p[1] = μ, p[2] = μ_3, p[3] = t0, p[4] = a, p[5] = ω_s
    p_ode = (params.mu, params.μ_3, t0, params.a, params.ω_s)
    prob = ODEProblem(R3BP.rhs_bcr4bp_sv!, x0, (0,_tof), p_ode, callback=cbs,
        method=Tsit5(), reltol=1e-12, abstol=1e-12)
    sol = @suppress_err solve(prob);
    push!(sims_let, sol)
end

# trajectory plot
ptraj = plot(aspect_ratio=:equal, flip=false, legend=false, frame_style=:box)
for (idx,sol) in enumerate(sims_let)
    plot!(ptraj, Array(sol)[1,:], Array(sol)[2,:], flip=false,
        line_z=[el[2] for el in lets][idx], cmap=:winter)
end
display(ptraj)