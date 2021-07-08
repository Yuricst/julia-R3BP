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

# defining system parameters
include("get_gm.jl")
include("get_semiMajorAxis.jl")
include("get_cr3bp_param.jl")

# basic functions for CR3BP
include("lagrangePoints.jl")
include("jacobiConstant.jl")
include("equationsOfMotion.jl")
include("lpo_stability.jl")

# LPO & manifold
include("lpo_family.jl")
include("manifold.jl")
include("analyticalCollinearHalo.jl")

# differential correction
include("differential_correction/differentialCorrection_singleshoot.jl")
include("differential_correction/multipleshooting.jl")

# equations of motion for trajectory design
include("deltaV_transcription.jl")
include("equationsOfMotionWithThrust.jl")

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
