"""
Module related to restricted three body problem
"""
module R3BP

using Roots
using DifferentialEquations
using ForwardDiff
using Printf
using DataFrames

# basic functions for CR3BP
include("lagrangePoints.jl")
include("jacobiConstant.jl")
include("equationsOfMotion.jl")
include("lpo_stability.jl")

# manifold, halo initial guess, differential correction
include("manifold.jl")
include("analyticalCollinearHalo.jl")
include("differentialCorrection_singleshoot.jl")

# LPO family handling
include("lpo_family.jl")

# defining system parameters
include("get_gm.jl")
include("get_semiMajorAxis.jl")
include("get_cr3bp_param.jl")

# for trajectory design
include("deltaV_transcription.jl")
include("equationsOfMotionWithThrust.jl")

# Misc
include("unpack_kwargs.jl")

# utility in R3BP
export lagrangePoint, get_cr3bp_param

# RHS
export rhs_cr3bp_sv, rhs_cr3bp_svstm, rhs_pcr3bp_sv, rhs_pcr3bp_svstm, rhs_pcr3bp_svstm

# manifold function
export get_stm, scale_ϵ, get_manifold

# differential correction
export halo_analytical_construct, ssdc_periodic_xzplane

end
