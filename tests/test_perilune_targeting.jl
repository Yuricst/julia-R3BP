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

params = R3BP.get_cr3bp_param(399, 301)
mu = params.mu

# periodic orbit
X0 = [1.176924090973164, 0.0, -0.060210863312217, 0.0, -0.173836346247689, 0.0];
T = 3.385326412831325;

# compact generator function
i_strip    = 2
n_interval = 30
n_strip    = 50
tf_fwd     = 15.0
states_llo = R3BP.lpo2llo_target(
    mu, X0, T, i_strip, n_interval, n_strip, tf_fwd);

# propagate found trajectories backward
tf_bck = -3.0
prob_base = ODEProblem(R3BP.rhs_cr3bp_sv!, [1,0,0,0,1,0], (0, tf_bck), [mu,],
    method=Tsit5(), reltol=1e-12, abstol=1e-12,
)

traj_llo_targeting = []
@showprogress for state in states_llo  # length == n_strip
    _prob = remake(prob_base; u0=state)
    _sol = @suppress_err solve(_prob)
    push!(traj_llo_targeting, _sol)
end

# trajectory plot
ptraj = plot(frame_style=:box, aspect_ratio=:equal)
for sol in traj_llo_targeting
    plot!(ptraj, Array(sol)[1,:], Array(sol)[2,:], label=false)
end
display(ptraj)
