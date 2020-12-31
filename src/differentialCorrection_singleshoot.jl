"""
Single shooting differential correction
"""

using DifferentialEquations
using LinearAlgebra
using Printf


# -------------------------------------------------------------------------------------------- #
# pseudo-event function for easy handling
function pseudo_rhs_cr3bp_sv(u,p,t)
    """Right-hand side expression for state-vector in CR3BP
    Input:
        du : cache array of derivative of state-vector
        u : state-vector
        p : parameters, where p[1] = mu
        t : time
    """
    du = zeros(6)
    # unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]
    # compute distances
    r1 = sqrt( (x+p[1])^2 + y^2 + z^2 );
    r2 = sqrt( (x-1+p[1])^2 + y^2 + z^2 );
    # derivatives of positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]
    # derivatives of velocities
    du[4] = 2*vy + x - ((1-p[1])/r1^3)*(p[1]+x) + (p[1]/r2^3)*(1-p[1]-x);
    du[5] = -2*vx + y - ((1-p[1])/r1^3)*y - (p[1]/r2^3)*y;
    du[6] = -((1-p[1])/r1^3)*z - (p[1]/r2^3)*z;
    return du
end


function pseudo_rhs_pcr3bp_sv(u,p,t)
    """Right-hand side expression for state-vector in CR3BP
    Input:
        du : cache array of derivative of state-vector
        u : state-vector
        p : parameters, where p[1] = mu
        t : time
    """
    du = zeros(4)
    # unpack state
    x, y = u[1], u[2]
    vx, vy = u[3], u[4]
    # compute distances
    r1 = sqrt( (x+p[1])^2 + y^2 );
    r2 = sqrt( (x-1+p[1])^2 + y^2 );
    # derivatives of positions
    du[1] = u[3]
    du[2] = u[4]
    # derivatives of velocities
    du[3] = 2*vy + x - ((1-p[1])/r1^3)*(p[1]+x) + (p[1]/r2^3)*(1-p[1]-x);
    du[4] = -2*vx + y - ((1-p[1])/r1^3)*y - (p[1]/r2^3)*y;
    return du
end


# -------------------------------------------------------------------------------------------- #
# single=shooting differnetial correction
struct Struct_out_ssdc
    x0
    period
    sol
    flag
end


function ssdc_periodic_xzplane(mu, x0, period0; kwargs...)
    """Single-shooting differential correction for periodic trajectory with symmetry across xz-plane"""

    # ------- unpack kwargs ----- #
    if :maxiter in keys(kwargs)
        maxiter = kwargs[:n];
    else
        maxiter = 15
    end

    if :reltol in keys(kwargs)
        reltol = kwargs[:reltol];
    else
        reltol = 1e-12
    end

    if :abstol in keys(kwargs)
        abstol = kwargs[:abstol];
    else
        abstol = 1e-12
    end

    if :method in keys(kwargs)
        method = kwargs[:method];
    else
        method = Tsit5()
    end

    if :fix in keys(kwargs)
        fix = kwargs[:fix]
    else
        fix = "period"
    end

    if :tolDC in keys(kwargs)
        tolDC = kwargs[:tolDC]
    else
        tolDC = 1e-11
    end

    x0iter = deepcopy(x0)
    period = deepcopy(period0)

    # ----- iterate until convergence ----- #
    idx = 1
    flag = 0
    while idx < maxiter+1

        # define ODE problem
        if length(x0iter)==4
            x0_stm = hcat(x0iter, [1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1]);
            prob = ODEProblem(R3BP.rhs_pcr3bp_svstm!, x0_stm, period/2, (mu));
        elseif length(x0iter)==6
            x0_stm = hcat(x0iter, [1 0 0 0 0 0  0 1 0 0 0 0  0 0 1 0 0 0  0 0 0 1 0 0  0 0 0 0 1 0  0 0 0 0 0 1]);
            prob = ODEProblem(R3BP.rhs_cr3bp_svstm!, x0_stm, period/2, (mu));
        else
            error("x0 should be length 4 or 6")
        end
        # solve ODE problem
        sol = solve(prob, method, reltol=reltol, abstol=abstol);
        # final state and STM
        statef = sol.u[end][1:length(x0iter)]

        if length(x0iter)==4
            dstatef = pseudo_rhs_pcr3bp_sv(statef, (mu), 1.0)
        elseif length(x0iter)==6
            dstatef = pseudo_rhs_cr3bp_sv(statef, (mu), 1.0)
        end
        stm = reshape(sol.u[end][length(x0iter)+1:end], (length(x0iter),length(x0iter)))';

        # residual vector
        if length(x0iter)==4
            #        y         vx
            ferr = [ statef[2] statef[3] ]'
        elseif length(x0iter)==6
            #        y         vx        vz
            ferr = [ statef[2] statef[4] statef[6] ]'
        end

        # create free parameters vector and DF
        if fix=="vy" && length(x0iter)==4
            #      dy / dx    dvx / d(P/2)
            xi = [ x0iter[1]  period/2 ]'
            df = [ [ stm[2,1] dstatef[2] ];
                   [ stm[3,1] dstatef[3] ] ];
        elseif fix=="z" && length(x0iter)==6
            #      / dx      / dvy     / d(P/2)
            xi = [ x0iter[1] x0iter[5] period/2 ]'
            df = [ [ stm[2,1] stm[2,5] dstatef[2] ];
                   [ stm[4,1] stm[4,5] dstatef[4] ];
                   [ stm[6,1] stm[6,5] dstatef[6] ] ];
        elseif fix=="period" && length(x0iter)==4
            #      dy / dx    dvx / dvy
            xi = [ x0iter[1]  x0iter[4] ]'
            df = [ [ stm[2,1] stm[2,4] ];
                   [ stm[3,1] stm[3,4] ] ];
        elseif fix=="period" && length(x0iter)==6
            #      / dx      / dz      / dvy
            xi = [ x0iter[1] x0iter[3] x0iter[5] ]'
            df = [ [ stm[2,1] stm[2,3] stm[2,5] ];
                   [ stm[4,1] stm[4,3] stm[4,5] ];
                   [ stm[6,1] stm[6,3] stm[6,5] ] ];
        end
        # correct state
        xii = xi - inv(df)*ferr

        # update state
        if fix=="vy" && length(x0iter)==4
            #      x         P/2
            x0iter[1] = xii[1];
            period    = xii[2];
        elseif fix=="z" && length(x0iter)==6
            #      x         vy        P/2
            x0iter[1] = xii[1];
            x0iter[5] = xii[2];
            period    = xii[3];
        elseif fix=="period" && length(x0iter)==4
            #      x         vy
            x0iter[1] = xii[1];
            x0iter[4] = xii[2];
        elseif fix=="period" && length(x0iter)==6
            #      x         z         vy
            x0iter[1] = xii[1];
            x0iter[3] = xii[2];
            x0iter[5] = xii[3];
        end

        # check convergence
        if norm(ferr) < tolDC
            println("Converged!")
            flag = 1
            break
        else
            @printf("Current residual: %s\n", norm(ferr))
            idx += 1
        end
    end

    # return result
    if length(x0iter)==4
        prob = ODEProblem(R3BP.rhs_pcr3bp_sv!, x0iter, period, (mu));
    else
        prob = ODEProblem(R3BP.rhs_cr3bp_sv!, x0iter, period, (mu));
    end
    return Struct_out_ssdc(x0iter, period, solve(prob, method, reltol=reltol, abstol=abstol), flag);
end
