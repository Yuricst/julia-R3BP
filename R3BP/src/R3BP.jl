"""
Module related to restricted three body problem
"""
module R3BP

using Roots
using DifferentialEquations

# basic functions for CR3BP
include("lagrangePoints.jl")
include("jacobiConstant.jl")
include("equationsOfMotion.jl")

# manifold, halo initial guess, differential correction
include("manifold.jl")
include("analyticalCollinearHalo.jl")
include("differentialCorrection_singleshoot.jl")

# defining system parameters
include("get_gm.jl")
include("get_semiMajorAxis.jl")
include("get_cr3bp_param.jl")


#include("get_poincareSection_from_manifold.jl")


#export lagrangePoint, rhs_cr3bp_sv, rhs_cr3bp_svstm, rhs_pcr3bp_sv, rhs_pcr3bp_svstm, rhs_pcr3bp_svstm, get_stm, scale_ϵ, get_manifold

end
