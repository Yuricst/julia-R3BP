"""
Function for computing Jacobi constant
"""


function jacobiConstant(mu, x0)
    """Compute Jacobi constant"""
    if length(x0)==4
        x, y, vx, vy = x0[1], x0[2], x0[3], x0[4];
        z, vz = 0.0, 0.0;
    elseif length(x0)==6
        x, y, z    = x0[1], x0[2], x0[3];
        vx, vy, vz = x0[4], x0[5], x0[6];
    end
    r1 = sqrt( (x+mu)^2 + y^2 + z^2 );
    r2 = sqrt( (x-1+mu)^2 + y^2 + z^2 );
    mu1 = 1.0-mu;
    mu2 = mu;
    # compute augmented potential
    ubar = 0.5*(x^2 + y^2) + mu1/r1 + mu2/r2;
    jc   = 2*ubar - (vx^2 + vy^2 + vz^2);
    return jc
end
