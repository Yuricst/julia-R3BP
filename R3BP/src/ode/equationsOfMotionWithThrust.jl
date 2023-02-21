"""
Equations of motion with thrust
"""


"""
    rhs_pcr3bp_thrust_m1dir!(du,u,p,t)

Right-hand side expression for state-vector in PCR3BP with thrust in m1 direction

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector
    - `p`: parameters, where p[1] = μ, p[2] = τ, p[3] = mdot, p[4] = tmax
    - `t`: time
"""
function rhs_pcr3bp_thrust_m1dir!(du,u,p,t)
    # unpack parameters
    μ, τ, mdot, tmax = p
    # unpack state
    x, y = u[1], u[2]
    vx, vy = u[3], u[4]
    # compute distances
    r1 = sqrt( (x+μ)^2 + y^2 );
    r2 = sqrt( (x-1+μ)^2 + y^2 );

    # compute thrust in m1 direction
    thrust_vec = τ*tmax * [x+μ, y] / r1

    # derivatives of positions
    du[1] = u[3]
    du[2] = u[4]
    # derivatives of velocities
    du[3] =  2*vy + x - ((1-μ)/r1^3)*(μ+x) + (μ/r2^3)*(1-μ-x) + thrust_vec[1]/u[5]
    du[4] = -2*vx + y - ((1-μ)/r1^3)*y      - (μ/r2^3)*y        + thrust_vec[2]/u[5]
    # mass derivative
    du[5] = -mdot*abs(τ)
end



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



"""
    rhs_bcr4bp_sunenhanced!(du,u,p,t)

BCR4BP equation of motion

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector, including mass
    - `p`: parameters, where p[1] = μ, p[2] = μ_3, p[3] = t0, p[4] = a, p[5] = ω_s,
                             p[6] = τ, p[7] = mdot, p[8] = tmax
    - `t`: time
"""
function rhs_bcr4bp_sunenhanced!(du,u,p,t)
    # unpack arguments
    μ, μ_3, t0, a_s, ω_s, τ, mdot, tmax = p
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
    sun_v = [
        -(μ_3/r3^3)*(x-xs) - (μ_3/a_s^3)*xs,
        -(μ_3/r3^3)*(y-ys) - (μ_3/a_s^3)*ys,
        -(μ_3/r3^3)*(z)    - (μ_3/a_s^3)*zs
    ]
    thrust_vec = τ*tmax*sun_v/norm(sun_v)

    # position-state derivative
    du[1] = vx
    du[2] = vy
    du[3] = vz

    # velocity derivatives
    du[4] =  2*vy + x - ((1-μ)/r1^3)*(μ+x)  + (μ/r2^3)*(1-μ-x)  + sun_v[1] + thrust_vec[1]/mass
    du[5] = -2*vx + y - ((1-μ)/r1^3)*y      - (μ/r2^3)*y        + sun_v[2] + thrust_vec[2]/mass
    du[6] =           - ((1-μ)/r1^3)*z      - (μ/r2^3)*z        + sun_v[3] + thrust_vec[3]/mass

    # mass derivative
    du[7] = -mdot*τ
end


# --------------------------------------------------------------------------- #
# EoM with control-law given as callable

"""
    rhs_pcr3bp_control!(du,u,p,t)

Right-hand side expression for state-vector in PCR3BP with control.
Expect control-law function in the form:

    ```julia
        tx, ty = control_law(du,u,p,t)
    ```

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector
    - `p`: parameters, where p[1] = μ
    - `t`: time
"""
function rhs_pcr3bp_control!(du,u,p,t)
    # unpack parameters
    μ, mdot, τ = p[1:3]
    control_law = p[end]
    # get thrust direction
    tx, ty, mdot_scaled = control_law(du,u,p[1:end-1],t)
    # unpack state
    x, y = u[1], u[2]
    vx, vy = u[3], u[4]
    m = u[5]
    # compute distances
    r1 = sqrt( (x+μ)^2 + y^2 );
    r2 = sqrt( (x-1+μ)^2 + y^2 );
    # derivatives of positions
    du[1] = u[3]
    du[2] = u[4]
    # derivatives of velocities
    du[3] =  2*vy + x - ((1-μ)/r1^3)*(μ+x) + (μ/r2^3)*(1-μ-x) + tx/m
    du[4] = -2*vx + y - ((1-μ)/r1^3)*y     - (μ/r2^3)*y       + ty/m
    # mass derivative
    du[5] = -mdot_scaled
end



"""
    rhs_cr3bp_control!(du,u,p,t)

Right-hand side expression for state-vector in CR3BP with control.
Expect control-law function in the form:

    ```julia
        tx, ty, tz = control_law(du,u,p,t)
    ```

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector
    - `p`: parameters, where p[1] = μ, p[2] = mdot, p[3] = control-law, p[4:6] = τ, θ, β
    - `t`: time
"""
function rhs_cr3bp_control!(du,u,p,t)
    # unpack parameters
    μ, mdot, τ = p[1:3]
    control_law = p[end]
    # get thrust direction
    tx, ty, tz = control_law(du,u,p[1:end-1],t)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    m = u[7]
    # compute distances
    r1 = sqrt( (x+μ)^2 + y^2 + z^2 );
    r2 = sqrt( (x-1+μ)^2 + y^2 + z^2 );
    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]
    # derivatives of velocities
    du[4] = 2*vy + x - ((1-μ)/r1^3)*(μ+x) + (μ/r2^3)*(1-μ-x) + tx/m
    du[5] = -2*vx + y - ((1-μ)/r1^3)*y - (μ/r2^3)*y + ty/m
    du[6] = -((1-μ)/r1^3)*z - (μ/r2^3)*z + tz/m
    # mass derivative
    du[7] = -mdot*τ
end



"""
    rhs_bcr4bp_control!(du,u,p,t)

BCR4BP equation of motion with control.
Expect control-law function in the form:

    ```julia
        tx, ty, tz = control_law(du,u,p,t)
    ```

# Arguments
    - `du`: cache array of duative of state-vector
    - `u`: state-vector
    - `p`: parameters, where p[1] = μ, p[2] = μ_3, p[3] = t0, p[4] = a, p[5] = ω_s
    - `t`: time
"""
function rhs_bcr4bp_control!(du,u,p,t)
    # unpack parameters
    μ, μ_3, t0, a_s, ω_s, mdot, τ = p[1:7]
    control_law = p[end]
    # get thrust direction
    tx, ty, tz = control_law(du,u,p[1:end-1],t)

    # unpack state
    x, y, z, vx, vy, vz, m = u

    # calculate radii
    r1 = sqrt( (x+μ)^2 + y^2 + z^2 )
    r2 = sqrt( (x-1+μ)^2 + y^2 + z^2 )

    # sun position
    xs = a_s*cos(ω_s*t + t0)
    ys = a_s*sin(ω_s*t + t0)
    zs = 0.0
    r3 = sqrt( (x-xs)^2 + (y-ys)^2 + (z-zs)^2 )

    # position-state derivative
    du[1] = vx
    du[2] = vy
    du[3] = vz

    # velocity derivatives
    du[4] =  2*vy + x - ((1-μ)/r1^3)*(μ+x) + (μ/r2^3)*(1-μ-x) + ( -(μ_3/r3^3)*(x-xs) - (μ_3/a_s^3)*xs ) + tx/m
    du[5] = -2*vx + y - ((1-μ)/r1^3)*y      - (μ/r2^3)*y        + ( -(μ_3/r3^3)*(y-ys) - (μ_3/a_s^3)*ys ) + ty/m
    du[6] =           - ((1-μ)/r1^3)*z      - (μ/r2^3)*z        + ( -(μ_3/r3^3)*(z)    - (μ_3/a_s^3)*zs ) + tz/m
    # mass derivative
    du[7] = -mdot*τ
end
