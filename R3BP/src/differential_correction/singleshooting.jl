"""
Single shooting differential correction algorithm for periodic orbits in R3BP systems
"""


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


"""
Right-hand side expression for state-vector in CR3BP

Input:
    du : cache array of derivative of state-vector
    u : state-vector
    p : parameters, where p[1] = mu
    t : time
"""
function pseudo_rhs_pcr3bp_sv(u,p,t)
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
# single-shooting differnetial correction
struct Struct_out_ssdc
    x0
    period::Float64
    sol::ODESolution
    flag::Int
    fiters::Vector
end

struct SingleShootingSolution
    x0::Vector
    period::Float64
    flag::Int
    fiters::Vector
end

"""
    ssdc_periodic_xzplane(p, x0, period0; kwargs...)

Single-shooting differential correction for periodic trajectory with symmetry across xz-plane

# Arguments
    p (tuple): parameters for DifferentialEquations
    x0 (Array):
    period0 (float): period of LPO; shooting is done based on perpendicular plane crossing at period0/2
    kwargs:
        maxiter
        reltol
        abstol
        method
        fix
        tolDC
        system (str): "cr3bp" or "er3bp"
        stm_option (str): "analytical" or "ad"

# Returns
    (struct): struct with fields: x0, period, sol, flag, fiters
"""
function ssdc_periodic_xzplane(p, x0::Vector, period0::Real; kwargs...)
    # unpack kwargs
    kwargs_dict = Dict(kwargs)
    maxiter = assign_from_kwargs(kwargs_dict, :maxiter, 15, Int)
    method  = assign_from_kwargs(kwargs_dict, :method, Tsit5())
    reltol  = assign_from_kwargs(kwargs_dict, :reltol, 1.e-13, Float64)
    abstol  = assign_from_kwargs(kwargs_dict, :abstol, 1.e-14, Float64)
    fix     = assign_from_kwargs(kwargs_dict, :fix, "period", String)

    tolDC      = assign_from_kwargs(kwargs_dict, :tolDC, 1.0e-12, Float64)
    system     = assign_from_kwargs(kwargs_dict, :system, "cr3bp", String)
    verbose    = assign_from_kwargs(kwargs_dict, :verbose, false, Bool)
    stm_option = assign_from_kwargs(kwargs_dict, :stm_option, "analytical", String)
    maxiter    = assign_from_kwargs(kwargs_dict, :maxiter, 15, Int)

    # unpack mu
    mu = p[1]

    # initialize with array and period
    x0iter = deepcopy(x0)
    period = deepcopy(period0)

    # if stm_option is AD, create a base problem for later
    #if cmp(stm_option, "ad")==0
    if cmp(system, "cr3bp")==0
        if cmp(stm_option, "ad")==0
            prob_base = ODEProblem(
                R3BP.rhs_cr3bp_sv!, x0iter, period/2, p,
                reltol=reltol, abstol=abstol
            );
        else
            prob_base = ODEProblem(
                R3BP.rhs_cr3bp_svstm!, x0iter, period/2, p,
                reltol=reltol, abstol=abstol
            );
        end
    elseif cmp(system, "er3bp")==0
        if cmp(stm_option, "ad")==0
            prob_base = ODEProblem(
                R3BP.rhs_er3bp_sv!, x0iter, period/2, p,
                reltol=reltol, abstol=abstol
            );
        else
            prob_base = ODEProblem(
                R3BP.rhs_er3bp_svstm!, x0iter, period/2, p,
                reltol=reltol, abstol=abstol
            );
        end
    end
    #end

    # ----- iterate until convergence ----- #
    idx = 1
    flag = 0
    fiters = []
    while idx < maxiter+1
        # remake ODE problem
        if cmp(stm_option, "ad")==0
            prob = remake(prob_base; tspan=(0.0, period/2), u0=x0iter, p=p)
        else
            if length(x0iter)==4
                x0_stm = vcat(x0iter, reshape(I(4), (16,)))[:]
            elseif length(x0iter)==6
                x0_stm = vcat(x0iter, reshape(I(6), (36,)))[:]
            end
            prob = remake(prob_base; tspan=(0.0, period/2), u0=x0_stm, p=p)
        end

        # solve ODE problem
        sol = solve(prob, method);
        # final state and STM
        statef = sol.u[end][1:length(x0iter)]

        # get derivative of state
        if length(x0iter)==4
            dstatef = zeros(4)
            if cmp(system, "cr3bp")==0
                _ = R3BP.rhs_pcr3bp_sv!(dstatef, statef, p, sol.t[end])
            else
                error("Planar ER3BP not implemented!")
            end
        elseif length(x0iter)==6
            dstatef = zeros(6)
            if cmp(system, "cr3bp")==0
                _ = R3BP.rhs_cr3bp_sv!(dstatef, statef, p, sol.t[end])
            else
                _ = R3BP.rhs_er3bp_sv!(dstatef, statef, p, sol.t[end])
            end
        end

        # get state-transition matrix
        if cmp(stm_option, "analytical")==0
            # STM using analytical expression
            stm = reshape(sol.u[end][length(x0iter)+1:end], (length(x0iter),length(x0iter)))';
        else
            # STM using AD
            stm = ForwardDiff.jacobian(x0iter -> get_statef(prob_base, x0iter, period/2, p, method), x0iter)
        end

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
        push!(fiters, norm(ferr))

        # update state
        if fix=="vy" && length(x0iter)==4
            #      x         P/2
            x0iter[1] = xii[1];
            period    = 2xii[2];
        elseif fix=="z" && length(x0iter)==6
            #      x         vy        P/2
            x0iter[1] = xii[1];
            x0iter[5] = xii[2];
            period    = 2xii[3];
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
            if verbose==true
                @printf("Converged at iteration %i\n", idx)
            end
            flag = 1
            break
        else
            if verbose==true
                @printf("Iteration %i: residual: %s\n", idx, norm(ferr))
            end
            idx += 1
        end
    end

    # construct problem
    return SingleShootingSolution(x0iter, period, flag, fiters)
    #sol_final = solve(prob, method, reltol=reltol, abstol=abstol)
    #return Struct_out_ssdc(x0iter, period, sol_final, flag, fiters)
end



"""
Function to remake ODE problem and get final state
"""
function get_statef(prob, x0, tf, p, method)
    _prob = remake(prob, u0=x0, p=p)
    sol = solve(_prob, method, saveat=tf)
    sol.u[end]
end


function print_stm(stm)
    nrow, ncol = size(stm)
    for i = 1:nrow
        for j = 1:ncol
            @printf("% 1.6e   ", stm[i,j])
            if j == ncol
                print("\n")
            end
        end
    end
end
