"""
equations of motion with thrust
"""


"""
    rhs_bcr4bp_thrust!(du,u,p,t)

BCR4BP equation of motion

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector, including mass
    - `p`: parameters, where p[1] = μ, p[2] = μ_3, p[3] = t0, p[4] = a, p[5] = ω_s,
                             p[6] = τ, p[7] = θ, p[8] = β, p[9] = mdot, p[10] = tmax
    - `t`: time
"""
function rhs_bcr4bp_thrust!(du,u,p,t)
    # unpack arguments
    μ, μ_3, t0, a_s, ω_s, τ, θ, β, mdot, tmax = p
    # decompose state
    x, y, z, vx, vy, vz, mass = u

    # calculate radii
    r1 = sqrt( (x+μ)^2 + y^2 + z^2 )
    r2 = sqrt( (x-1+μ)^2 + y^2 + z^2 )

    # sun position
    xs = a_s*cos(ω_s*t + t0)
    ys = a_s*sin(ω_s*t + t0)
    zs = 0.0
    r3 = sqrt( (x-xs)^2 + (y-ys)^2 + (z-zs)^2 )

    # compute delta-V vector
    dir_v = dv_lvlh2inertial(u[1:6], [τ, θ, β])
    ts = tmax*dir_v

    # position-state derivative
    du[1] = vx
    du[2] = vy
    du[3] = vz

    # velocity derivatives
    du[4] =  2*vy + x - ((1-μ)/r1^3)*(μ+x)  + (μ/r2^3)*(1-μ-x)  + ( -(μ_3/r3^3)*(x-xs) - (μ_3/a_s^3)*xs ) + ts[1]/mass
    du[5] = -2*vx + y - ((1-μ)/r1^3)*y      - (μ/r2^3)*y        + ( -(μ_3/r3^3)*(y-ys) - (μ_3/a_s^3)*ys ) + ts[2]/mass
    du[6] =           - ((1-μ)/r1^3)*z      - (μ/r2^3)*z        + ( -(μ_3/r3^3)*(z)    - (μ_3/a_s^3)*zs ) + ts[3]/mass

    # mass derivative
    du[7] = -mdot*τ
end