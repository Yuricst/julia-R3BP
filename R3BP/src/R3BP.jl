"""
Module related to restricted three body problem
"""
module R3BP

using Roots
using DifferentialEquations
using DifferentialEquations.EnsembleAnalysis
using ForwardDiff
using Printf
using DataFrames
using ProgressMeter
using JSON
using Dierckx
using Suppressor: @suppress_err
using Plots

using AstrodynamicsBase
#using Plots

# defining system parameters
#include("get_gm.jl")
#include("get_semiMajorAxis.jl")
include("get_cr3bp_param.jl")

# basic functions for CR3BP
include("lagrangePoints.jl")
include("jacobiConstant.jl")
include("ode/equationsOfMotion.jl")
include("librationpointorbit/lpo_stability.jl")

# LPO & manifold
include("librationpointorbit/lpo_family.jl")
include("librationpointorbit/manifold.jl")
include("librationpointorbit/analyticalCollinearHalo.jl")
include("librationpointorbit/stretching.jl")

# differential correction
include("differential_correction/singleshooting.jl")
include("differential_correction/multipleshooting.jl")

# equations of motion for trajectory design
include("ode/deltaV_transcription.jl")
include("ode/equationsOfMotionWithThrust.jl")

# LET design
include("let/perilune_targeting.jl")
include("let/search_bcr4bp.jl")

# miscellaneous
include("misc_tools.jl")
include("plot_support.jl")


# utility in R3BP
export lagrangePoint, get_cr3bp_param

# RHS
export rhs_cr3bp_sv, rhs_cr3bp_svstm, rhs_pcr3bp_sv, rhs_pcr3bp_svstm, rhs_pcr3bp_svstm

# manifold function
export get_stm, scale_ϵ, get_manifold

# differential correction
export halo_analytical_construct, ssdc_periodic_xzplane
export multiple_shooting

# miscellaneous
export sol_to_arrays

end
