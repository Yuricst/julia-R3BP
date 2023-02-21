"""
    dv_lvlh2inertial(state0, vinf_params)

Convert velocity vector in LVLH frame to inertial frame

Args:
    `state0::Vector{Float}`: state before delta v is applied, should have at least length 6 or 7 (rx,ry,rz,vx,vy,vz,m)
    `vinf_params::Vector{Float}`: [τ, θ, β]

Returns:
   ::Vector{Float}`: delta-v vector direction scaled by tau, in heliocentric frame
"""
function dv_lvlh2inertial(state0, vinf_params)
    # get delta-V in LVLH frame
    τ, θ, β = vinf_params[1], vinf_params[2], vinf_params[3]
    dv_lvlh = τ*[cos(θ)*cos(β), sin(θ)*cos(β), sin(β)]
    # conversion matrix following Markley and Crassidis 2014 pg.36
    r, v = state0[1:3], state0[4:6]
    o3I = -r / norm(r)
    o2I = -cross(r,v)/norm(cross(r,v))
    o1I = cross(o2I, o3I)
    A_IO = reshape(vcat(o1I, o2I, o3I), 3,3)
    dv_inertial = A_IO*dv_lvlh
    return dv_inertial
end
