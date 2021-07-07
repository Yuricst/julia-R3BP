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
