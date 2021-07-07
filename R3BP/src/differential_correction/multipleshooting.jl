"""
Multiple-shooting differential correction algorithm
"""


"""
    multiple_shooting(prob_stm, x0s::Vector, tofs, tolDC; kwargs...)

Multiple shooting to correct trajectory

# Arguments
    - `prob_stm::ODEProblem`: ODEProblem for propagating state+STM
    - `x0s::Array`: list of n initial guess states
    - `tofs::Array`: propagation times of the first (n-1) initial guesses
    - `tolDC::Float64`: tolerance on multiple-shooting
"""
function multiple_shooting(prob_stm::ODEProblem, x0s::Array, tofs::Array, tolDC::Float64; kwargs...)
    # ------------- unpack kwargs ------------- #
    kwargs_dict = Dict(kwargs)

    # options for multiple shooting
    maxiter  = assign_from_kwargs(kwargs_dict, :maxiter, 15, Int)
    fix_time = assign_from_kwargs(kwargs_dict, :fix_time, true, Bool)
    periodic = assign_from_kwargs(kwargs_dict, :periodic, false, Bool)
    fix_x0   = assign_from_kwargs(kwargs_dict, :fix_x0, false, Bool)
    fix_xf   = assign_from_kwargs(kwargs_dict, :fix_xf, false, Bool)

    # ODE settings
    method   = assign_from_kwargs(kwargs_dict, :method, Tsit5())
    reltol   = assign_from_kwargs(kwargs_dict, :reltol, 1.e-12, Float64)
    abstol   = assign_from_kwargs(kwargs_dict, :abstol, 1.e-12, Float64)
    p        = assign_from_kwargs(kwargs_dict, :p, [])
    rhs!     = assign_from_kwargs(kwargs_dict, :rhs!, nothing)
    use_ensemble = assign_from_kwargs(kwargs_dict, :use_ensemble, false, Bool)

    # misc function settings
    verbose  = assign_from_kwargs(kwargs_dict, :verbose, true, Bool)

    # exception if rhs! is not provided but fix_time == false
    if fix_time == false && isnothing(rhs!) == true
        error("Provide callable rhs! if fix_time==false")
    end

    # ------------- storage setups ------------- #
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

    # ------------- iterate multiple-shooting ------------- #
    for iter = 1:maxiter

        # re-make ODE problem and integrate until tofs[ix]
        if use_ensemble==false
            multiple_shooting_propagate_sequential!(df, xfs, prob_stm, n, n_sv, x0_vec, tofs, fix_time, method, reltol, abstol)
        else
            println("FIXME!")
            multiple_shooting_propagate_ensemble!(df, xfs, prob_stm, n, n_sv, x0_vec, tofs, fix_time, method, reltol, abstol)
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
        err = norm(ferr)
        if err < tolDC
            convflag = 1
            __print_verbose("Tolerance achieved at iteration $iter with error = $err", verbose)
            break
        else
            __print_verbose("Multiple shooting iteration $iter, error = $err", verbose)
        end

        # apply update on initial guess
        x0_vec = multiple_shooting_update!(x0_vec, ferr, df)
    end
    return x0_vec, tofs, convflag
end


"""
    multiple_shooting_propaagte_sequential!(df::Matrix{Float64}, xfs::Vector{Float64}, prob_stm::ODEProblem, n::Int, n_sv::Int, x0_vec::Vector{Float64}, tofs::Array, fix_time::Bool, method, reltol::Float64, abstol::Float64)

Propagates nodes and mutates DF and xfs; solved sequentially using `remake()`
"""
function multiple_shooting_propagate_sequential!(df::Matrix{Float64}, xfs::Vector{Float64}, prob_stm::ODEProblem, n::Int, n_sv::Int, x0_vec::Vector{Float64}, tofs::Array, fix_time::Bool, method, reltol::Float64, abstol::Float64)
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
        # if fix_time == false, fill in df/dx into DF matrix
        duf = zeros(length(x0_stm))
        if fix_time == false
            _ = rhs!(duf, sol.u[end], p, sol.t[end])
        end
    end
    return
end


"""
    multiple_shooting_propaagte_ensemble!(df::Matrix{Float64}, xfs::Vector{Float64}, prob_stm::ODEProblem, n::Int, n_sv::Int, x0_vec::Vector{Float64}, tofs::Array, fix_time::Bool, method, reltol::Float64, abstol::Float64)

Propagates nodes and mutates DF and xfs; solved sequentially using `ensembleSimulation()`
"""
function multiple_shooting_propagate_ensemble!(df::Matrix{Float64}, xfs::Vector{Float64}, prob_stm::ODEProblem, n::Int, n_sv::Int, x0_vec::Vector{Float64}, tofs::Array, fix_time::Bool, method, reltol::Float64, abstol::Float64)
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
        # if fix_time == false, fill in df/dx into DF matrix
        duf = zeros(length(x0_stm))
        if fix_time == false
            _ = rhs!(duf, sol.u[end], p, sol.t[end])
        end
    end
    return
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
