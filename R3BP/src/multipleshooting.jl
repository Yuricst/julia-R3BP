"""
Multiple-shooting algorithm
"""


"""
    multiple_shooting(prob_stm, x0s, tofs, tolDC; kwargs...)

Multiple shooting to correct trajectory

# Arguments
    - `prob_stm::ODEProblem`: ODEProblem for propagating state+STM
    - `x0s::Array`: list of n initial guess states
    - `tofs::Array`: propagation times of the first (n-1) initial guesses
    - `tolDC::Float64`: tolerance on multiple-shooting
"""
function multiple_shooting(prob_stm::ODEProblem, x0s::Array, tofs::Array, tolDC::Float64; kwargs...)
    # unpack kwargs
    kwargs_dict = Dict(kwargs)
    fix_time = assign_from_kwargs(kwargs_dict, :fix_time, true)
    periodic = assign_from_kwargs(kwargs_dict, :periodic, false)
    fix_idx = assign_from_kwargs(kwargs_dict, :fix_idx, nothing)
    maxiter = assign_from_kwargs(kwargs_dict, :maxiter, 15)
    method  = assign_from_kwargs(kwargs_dict, :method, Tsit5())
    reltol  = assign_from_kwargs(kwargs_dict, :reltol, 1.e-12)
    abstol  = assign_from_kwargs(kwargs_dict, :abstol, 1.e-12)
    verbose = assign_from_kwargs(kwargs_dict, :verbose, true)
    p       = assign_from_kwargs(kwargs_dict, :p, [])
    rhs!    = assign_from_kwargs(kwargs_dict, :rhs!, nothing)

    # exception if rhs! is not provided but fix_time == false
    if fix_time == false && isnothing(rhs!) == true
        error("Provide callable rhs! if fix_time==false")
    end

    # get state-vector length
    n_sv = length(x0s[1])

    # get number of nodes from x0s
    n = length(x0s)
    __print_verbose("Multiple-shooting using $n nodes, periodic: $periodic", verbose)
    if length(tofs) != n-1
        error("Must follow length(x0s)-1 == length(tofs)!")
    elseif n == 1
        error("Provide at least 2 nodes for multiple-shooting!")
    end

    # prepare initial guess vector
    if fix_time==true
        x0_vec = zeros(n_sv*n)
    else
        x0_vec = zeros(n_sv*n + n-1)
    end
    # append state-vector into initial guess
    for i = 1:n
        # append state-vector
        x0_vec[(i-1)*n_sv + 1 : i*(n_sv)] = x0s[i]
    end
    # append integration times into initial guess
    if fix_time == false
        if i != n
            x0_vec[n_sv*n+i] = tofs[i]
        end
    end

    # prepare DF matrix
    if periodic==false
        df = zeros(n_sv*(n-1), length(x0_vec))   # DF matrix
    else
        df = zeros(n_sv*(n-1)+n_sv, length(x0_vec))   # DF matrix
    end
    # fill-in negative identity into DF matrix
    for j = 1:n-1
        df[(j-1)*n_sv+1:j*n_sv, j*n_sv+1:(j+1)*n_sv] = -I(n_sv)
    end
    # connect initial and final states if periodic solution is sought
    if periodic==true
        df[end-n_sv+1:end, 1:n_sv]              = -I(n_sv)
        df[end-n_sv+1:end, n_sv*(n-1)+1:n*n_sv] =  I(n_sv)
    end

    # initialize storages
    convflag = 0                  # convergence flag
    xfs  = zeros(n_sv*(n-1))      # derivatives
    if periodic == false
        ferr = zeros(n_sv*(n-1))  # error vector
    else
        ferr = zeros(n_sv*n)      # error vector
    end

    # iterate multiple-shooting method
    for iter = 1:maxiter

        # re-make ODE problem and integrate until tofs[ix]
        for i = 1:n-1
            x0_stm = vcat(x0_vec[(i-1)*n_sv + 1 : i*(n_sv)], reshape(I(n_sv), (n_sv*n_sv,)))[:]
            if fix_time == true
                tprop = tofs[i]
            else
                tprop = x0_vec[n_sv*n+i]
            end
            _prob = remake(prob_stm; tspan=(0.0, tprop), u0=x0_stm)
            sol = DifferentialEquations.solve(_prob, method, reltol=reltol, abstol=abstol)
            # construct vector of propagation results
            xfs[(i-1)*n_sv + 1 : i*(n_sv)] = sol.u[end][1:n_sv]
            # fill in STM into DF matrix
            df[(i-1)*n_sv+1:i*n_sv, (i-1)*n_sv+1:i*n_sv] = transpose( reshape(sol.u[end][n_sv+1:end], (n_sv,n_sv)) )
            if fix_time == false
                duf = rhs!(zeros(length(x0_stm)), sol.u[end], p, 0.0)
            end
        end

        # compute error vector
        if periodic == false
            ferr = xfs - x0_vec[n_sv+1:n_sv*n]
        else
            ferr[1:n_sv*(n-1)] = xfs - x0_vec[n_sv+1:n_sv*n]
            ferr[n_sv*(n-1)+1:end] = x0_vec[end-n_sv+1:end] - x0_vec[1:n_sv]     # glueing final - initial states
            #ferr[n_sv*(n-1)+1:end] = x0_vec[end-n_sv+1:end] - x0s[1]   # glueing final position
            #ferr[n_sv*(n-1)+1:end] = x0_vec[1:n_sv] - x0s[1]   # glueing initial position
        end

        # check if tolerance is achieved
        if norm(ferr) < tolDC
            err = norm(ferr)
            convflag = 1
            __print_verbose("Tolerance achieved at iteration $iter with error = $err", verbose)
            break
        else
            err = norm(ferr)
            __print_verbose("Multiple shooting iteration $iter, error = $err", verbose)
        end

        # apply update on initial guess
        x0_vec = multiple_shooting_update!(x0_vec, ferr, df)
    end

    return x0_vec, tofs, convflag
end



"""
    multiple_shooting_update!(x, ferr, df)

Apply multiple-shooting update
"""
function multiple_shooting_update!(x, ferr, df)
    nf = length(ferr)
    nx = length(x)
    if nf == nx
        x = x - pinv(df)*ferr
    elseif nx > nf      # under-determined system: minimum-norm solution
        dfT = transpose(df)
        x = x - dfT * pinv(df * dfT) * ferr
    else                # over-determined system: least-square solution
        dfT = transpose(df)
        x = x - pinv(dfT * df) * dfT * ferr
    end
    return x
end



"""
Print message if verbose==true
"""
function __print_verbose(message, verbose::Bool)
    if verbose==true
        println(message)
    end
end
