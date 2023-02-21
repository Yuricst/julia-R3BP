"""
Function for computing Jacobi constant
"""


"""
    jacobiConstant(μ::Float64, sv)

Compute Jacobi constant

# Arguments
    - `μ::Float64`: CR3BP parameter
    - `sv::Vector`: state-vector, planar or spatial
"""
function jacobiConstant(μ::Float64, sv)
    if length(sv)==4
        x, y, vx, vy = sv
        z, vz = 0.0, 0.0;
    elseif length(sv)==6
        x, y, z, vx, vy, vz = sv
    end
    r1 = sqrt( (x+μ)^2 + y^2 + z^2 );
    r2 = sqrt( (x-1+μ)^2 + y^2 + z^2 );
    μ1 = 1.0-μ;
    μ2 = μ;
    # compute augmented potential
    ubar = 0.5*(x^2 + y^2) + μ1/r1 + μ2/r2;
    jc   = 2*ubar - (vx^2 + vy^2 + vz^2);
    return jc
end


"""
Compute state from Jacobi-constant
"""
function velocity_from_jacobiConstant(μ::Float64, jc::Float64; kwargs...)
    kwargs_dict = Dict(kwargs)
    # main manifold options
    x = assign_from_kwargs(kwargs_dict, :x, 1.0)
    y = assign_from_kwargs(kwargs_dict, :y, 0.0)
    z = assign_from_kwargs(kwargs_dict, :z, 0.0)
    vx = assign_from_kwargs(kwargs_dict, :vx, 0.0)
    vy = assign_from_kwargs(kwargs_dict, :vy, nothing)
    vz = assign_from_kwargs(kwargs_dict, :vz, 0.0)
    # compute velocity from jacobi constant
    r1 = sqrt( (x+μ)^2 + y^2 + z^2 );
    r2 = sqrt( (x-1+μ)^2 + y^2 + z^2 );
    if isnothing(vy) == true
        kinetic = vx^2 + vz^2
    end

    calc = sqrt((x^2 + y^2) + 2((1-μ)/r1 + μ/r2) - kinetic - jc)

    return calc
end
