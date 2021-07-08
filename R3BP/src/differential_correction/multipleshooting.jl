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
	- `maxiter::Int`: max iteration used by multiple-shooting algorithm
	- `fix_time::Bool`: whether to use tofs as part of the free-parameters
	- `fix_x0::Bool`: whether to freeze initial state-vector to provided value
	- `fix_xf::Bool`: whether to freeze final state-vector to provided value
	- `method::ODESolver`: `DifferentialEquations` ODE algorithm
	- `reltol::Float64`: `DifferentialEquations` relative tolerance
	- `abstol::Float64`: `DifferentialEquations` absolute tolerance
	- `p::Array`: `DifferentialEquations` parameters passed to rhs
	- `rhs::callable`: `DifferentialEquations`'s `rhs!(du,u,p,t)` function
	- `use_ensemble::Bool`: whether to use `DifferentialEquations`'s `EnsembleProblem`
	- `verbose::Bool`: whether to run algorithm in verbose mode

# Returns
	(`Array`): x0_vec, tofs, convflag
"""
function multiple_shooting(prob_stm::ODEProblem, x0s::Array, tofs::Array, tolDC::Float64; kwargs...)
    # ------------- unpack kwargs ------------- #
    kwargs_dict = Dict(kwargs)

    # options for multiple shooting
    maxiter  = assign_from_kwargs(kwargs_dict, :maxiter,  15, Int)
    fix_time = assign_from_kwargs(kwargs_dict, :fix_time, false, Bool)
    periodic = assign_from_kwargs(kwargs_dict, :periodic, false, Bool)
    fix_x0   = assign_from_kwargs(kwargs_dict, :fix_x0,   false, Bool)
    fix_xf   = assign_from_kwargs(kwargs_dict, :fix_xf,   false, Bool)

    # ODE settings
    method          = assign_from_kwargs(kwargs_dict, :method, Tsit5())
    reltol          = assign_from_kwargs(kwargs_dict, :reltol, 1.e-12, Float64)
    abstol          = assign_from_kwargs(kwargs_dict, :abstol, 1.e-12, Float64)
    p               = assign_from_kwargs(kwargs_dict, :p, [])
    rhs!            = assign_from_kwargs(kwargs_dict, :rhs, nothing)
    use_ensemble    = assign_from_kwargs(kwargs_dict, :use_ensemble, false, Bool)
	ensemble_method = assign_from_kwargs(kwargs_dict, :ensemble_method, EnsembleThreads())

    # misc function settings
    verbose  = assign_from_kwargs(kwargs_dict, :verbose, true, Bool)

	# ------------- multiple shooting size parameters ------------- #
    n_sv = length(x0s[1])	# get state-vector length
    n    = length(x0s)      # get number of nodes from x0s

	# ------------- check for exceptions ------------- #
	multiple_shooting_exceptions(n, tofs, fix_time, rhs!, fix_x0, fix_xf)

    # ------------- storage setups ------------- #
	__print_verbose(" ---------------- Multiple shooting algorithm ---------------- ", verbose)
	__print_verbose("    Multiple shooting solver settings", verbose)
	__print_verbose("        nodes         ..... $n", verbose)
	__print_verbose("        periodic      ..... $periodic", verbose)
	__print_verbose("        fix_time      ..... $fix_time", verbose)
	__print_verbose("        fix_x0        ..... $fix_x0", verbose)
	__print_verbose("        fix_xf        ..... $fix_xf", verbose)
	__print_verbose("        tolDC         ..... $tolDC", verbose)
	__print_verbose("        maxiter       ..... $maxiter", verbose)
	__print_verbose("    DifferentialEquations settings", verbose)
	__print_verbose("        ODE method    ..... $method", verbose)
	__print_verbose("        ODE reltol    ..... $reltol", verbose)
	__print_verbose("        ODE abstol    ..... $abstol", verbose)
	__print_verbose("        ODE Ensembles ..... $use_ensemble", verbose)
	__print_verbose("    Starting iterations", verbose)

    # prepare initial guess vector
	x0_vec, n_first_x, n_last_x = initialize_multiple_shooting_x(x0s, tofs, n, n_sv, fix_time, fix_x0, fix_xf)

    # prepare DF matrix
    if periodic == false
        df = zeros(n_sv*(n-1), length(x0_vec))   # DF matrix
    else
        df = zeros(n_sv*(n-1)+n_sv, length(x0_vec))   # DF matrix
    end
    # fill-in negative identity into DF matrix
    for j = 1:n-1
        df[(j-1)*n_sv+1:j*n_sv, j*n_sv+1:(j+1)*n_sv] = -I(n_sv)
    end
    # connect initial and final states if periodic solution is sought
    if periodic == true
        df[end-n_sv+1:end, 1:n_sv]              = -I(n_sv)
        df[end-n_sv+1:end, n_sv*(n-1)+1:n*n_sv] =  I(n_sv)
    end

    # initialize storages
    convflag = 0                  # convergence flag
    svf  = zeros(n_sv*(n-1))      # derivatives
    if periodic == false
        ferr = zeros(n_sv*(n-1))  # error vector
    else
        ferr = zeros(n_sv*n)      # error vector
    end

	# create problem-remake function for ensemble simulation
	if use_ensemble == true
		function prob_func(prob, i, repeat)
			if fix_time == true
				tprop = tofs[i]
			else
				tprop = x0_vec[n_sv*n+i]
			end
			remake(prob, u0=vcat(x0_vec[(i-1)*n_sv + 1 : i*(n_sv)], reshape(I(n_sv), (n_sv*n_sv,)))[:], tspan=(0.0, tprop))
		end
	end

    # ------------- iterate multiple-shooting ------------- #
    for iter = 1:maxiter

        # propagate ODE nodes and update DF
        if use_ensemble == false
            multiple_shooting_propagate_sequential!(df, svf, prob_stm, rhs!, p, n, n_sv, x0_vec, tofs, fix_time, fix_x0, fix_xf, method, reltol, abstol)
        else
            multiple_shooting_propagate_ensemble!(df, svf, prob_stm, rhs!, p, prob_func, ensemble_method, n, n_sv, x0_vec, tofs, fix_time, fix_x0, fix_xf, method, reltol, abstol)
        end

        # update error vector F
		multiple_shooting_update_f!(n, n_sv, ferr, svf, x0_vec, periodic, fix_x0, fix_xf)

        # check if tolerance is achieved
        err = norm(ferr)
        if err < tolDC
            convflag = 1
			conv_err_rounded = round(err, sigdigits=7)
			__print_verbose("        Iteration $iter, error = $conv_err_rounded < $tolDC :)", verbose)
            break
        else
            __print_verbose("        Iteration $iter, error = $err", verbose)
        end

        # update X
		#x0_vec = multiple_shooting_update_x!(x0_vec, ferr, df)
        multiple_shooting_update_x!(x0_vec, ferr, df)
    end
	__print_verbose(" ------------------------------------------------------------ ", verbose)
    return x0_vec, tofs, convflag
end


"""
    multiple_shooting_propaagte_sequential!(df::Matrix{Float64}, svf::Vector{Float64}, prob_stm::ODEProblem, n::Int, n_sv::Int, x0_vec::Vector{Float64}, tofs::Array, fix_time::Bool, fix_x0::Bool, fix_xf::Bool, method, reltol::Float64, abstol::Float64)

Propagates nodes and mutates DF and xfs; solved sequentially using `remake()`
"""
function multiple_shooting_propagate_sequential!(df::Matrix{Float64}, svf::Vector{Float64}, prob_stm::ODEProblem, rhs!, p, n::Int, n_sv::Int, x0_vec::Vector{Float64}, tofs::Array, fix_time::Bool, fix_x0::Bool, fix_xf::Bool, method, reltol::Float64, abstol::Float64)
    # re-make ODE problem and integrate until tofs[ix]
    for i = 1:n-1

		# integrate state-vector and STM
        x0_stm = vcat(x0_vec[(i-1)*n_sv + 1 : i*(n_sv)], reshape(I(n_sv), (n_sv*n_sv,)))[:]
        if fix_time == true
            tprop = tofs[i]
        else
            tprop = x0_vec[n_sv*n+i]
        end
        _prob = remake(prob_stm; tspan=(0.0, tprop), u0=x0_stm)
        sol = DifferentialEquations.solve(_prob, method, reltol=reltol, abstol=abstol)

        # construct vector of propagation results
        svf[(i-1)*n_sv + 1 : i*(n_sv)] = sol.u[end][1:n_sv]

        # fill in STM into DF matrix
        df[(i-1)*n_sv+1:i*n_sv, (i-1)*n_sv+1:i*n_sv] = transpose( reshape(sol.u[end][n_sv+1:end], (n_sv,n_sv)) )

		# if fix_time == false, fill in df/dx into DF matrix
		if fix_time == false   # FIXME
			duf = zeros(n_sv + n_sv*n_sv)
			_ = rhs!(duf, sol.u[end], p, sol.t[end])
			df[(i-1)*n_sv+1:i*n_sv, n_sv*n + i] = duf[1:n_sv]
		end
    end
    return
end



"""
    multiple_shooting_propaagte_ensemble!(df::Matrix{Float64}, svf::Vector{Float64}, prob_stm::ODEProblem, n::Int, n_sv::Int, x0_vec::Vector{Float64}, tofs::Array, fix_time::Bool, fix_x0::Bool, fix_xf::Bool, method, reltol::Float64, abstol::Float64)

Propagates nodes and mutates DF and svf; solved sequentially using `ensembleSimulation()`
"""
function multiple_shooting_propagate_ensemble!(df::Matrix{Float64}, svf::Vector{Float64}, prob_stm::ODEProblem, rhs!, p, prob_func, ensemble_method, n::Int, n_sv::Int, x0_vec::Vector{Float64}, tofs::Array, fix_time::Bool, fix_x0::Bool, fix_xf::Bool, method, reltol::Float64, abstol::Float64)
	# solve ensemble problem
	ensemble_prob = EnsembleProblem(prob_stm, prob_func=prob_func)
	sim = solve(ensemble_prob, method, ensemble_method, trajectories=n-1, reltol=reltol, abstol=abstol);

	# get iterable of final states
	#xf_iter = get_timestep(sim, END)
	# for xf in xf_iter
	# 	xfs[(i-1)*n_sv + 1 : i*(n_sv)] = xf_iter[1:n_sv]
	# end

	# process result
	for (i,sol) in enumerate(sim)

		# construct vector of propagation results
		svf[(i-1)*n_sv + 1 : i*(n_sv)] = sol.u[end][1:n_sv]

		# fill in STM into DF matrix
		df[(i-1)*n_sv+1:i*n_sv, (i-1)*n_sv+1:i*n_sv] = transpose( reshape(sol.u[end][n_sv+1:end], (n_sv,n_sv)) )

		# if fix_time == false, fill in df/dx into DF matrix
		if fix_time == false
			duf = zeros(n_sv + n_sv*n_sv)
			_ = rhs!(duf, sol.u[end], p, sol.t[end])
			df[(i-1)*n_sv+1:i*n_sv , n_sv*n + i] = duf[1:n_sv]
		end
	end
    return
end


"""
	multiple_shooting_update_f!(n::Int, n_sv::Int, ferr::Vector{Float64}, svf::Vector{Float64}, x0_vec::Vector{Float64}, periodic::Bool, fix_x0::Bool, fix_xf::Bool)
	
Mutate error-vector f
"""
function multiple_shooting_update_f!(n::Int, n_sv::Int, ferr::Vector{Float64}, svf::Vector{Float64}, x0_vec::Vector{Float64}, periodic::Bool, fix_x0::Bool, fix_xf::Bool)
	# update error vector F
	ferr[1:n_sv*(n-1)] = svf - x0_vec[n_sv+1:n_sv*n]
	if periodic == true
		ferr[n_sv*(n-1)+1:end] = x0_vec[n_sv*n-n_sv+1:n_sv*n] - x0_vec[1:n_sv]     # glueing final - initial states
		#ferr[n_sv*(n-1)+1:end] = x0_vec[end-n_sv+1:end] - x0s[1]   # glueing final position
		#ferr[n_sv*(n-1)+1:end] = x0_vec[1:n_sv] - x0s[1]   # glueing initial position
	end
	return
end


"""
    multiple_shooting_update_x!(x::Vector{Float64}, ferr::Vector{Float64}, df::Matrix{Float64})

Apply multiple-shooting update to mutate x
"""
function multiple_shooting_update_x!(x::Vector{Float64}, ferr::Vector{Float64}, df::Matrix{Float64})
    nf = length(ferr)
    nx = length(x)
    if nf == nx
        x[:] = x - pinv(df)*ferr
    elseif nx > nf      # under-determined system: minimum-norm solution
        dfT = transpose(df)
        x[:] = x - dfT * pinv(df * dfT) * ferr
    else                # over-determined system: least-square solution
        dfT = transpose(df)
        x[:] = x - pinv(dfT * df) * dfT * ferr
    end
    return
end


"""
Check multiple shooting exceptions
"""
function multiple_shooting_exceptions(n::Int, tofs::Array, fix_time::Bool, rhs!, fix_x0::Bool, fix_xf::Bool)
	# exception if rhs! is not provided but fix_time == false
	if fix_time == false && isnothing(rhs!) == true
        error("Provide callable rhs! if fix_time==false")
    end

	# exception if n < 3 but x0 or xf is fixed
	if (fix_x0 == true && fix_xf == false) || (fix_x0 == false && fix_xf == true)
		if n < 3
			error("Provide at least 3 nodes if fixing x0 or xf")
		end

	# exception if n < 4 but x0 and xf is fixed
	elseif fix_x0 == true && fix_xf == true
		if n < 4
			error("Provide at least 4 nodes if fixing x0 or xf")
		end
	end

	# exception if length(tof) != length(x0) - 1
	if length(tofs) != n-1
        error("Must follow length(x0s)-1 == length(tofs)")

	# exception if less than 2 nodes supplied
    elseif n == 1
        error("Provide at least 2 nodes for multiple-shooting")
    end
end


"""
	initialize_multiple_shooting_x(x0s::Array, n::Int, n_sv::Int, fix_time::Bool, fix_x0::Bool, fix_xf::Bool)

Initialize initial-guess vector x
"""
function initialize_multiple_shooting_x(x0s::Array, tofs::Array, n::Int, n_sv::Int, fix_time::Bool, fix_x0::Bool, fix_xf::Bool)
	# ------------- get number of free x's ------------- #
	# check if final sv is used in initial-guess vector
	if fix_x0 == true
		n_first_x = 2
	else
		n_first_x = 1
	end
	# check if last sv is used in initial-guess vector
	if fix_xf == true
		n_last_x = n-1
	else
		n_last_x = n
	end

	# ------------- initialize space ------------- #
	if fix_time == true
        x0_vec = zeros(n_sv*n)
    else
        x0_vec = zeros(n_sv*n + n-1)
    end

	# ------------- append state-vector ------------- #
    # append state-vector into initial guess
    for i = n_first_x:n_last_x
        x0_vec[(i-1)*n_sv + 1 : i*(n_sv)] = x0s[i]
    end

	# ------------- append time of flight ------------- #
    # append integration times into initial guess
	#println("n_first_x: $n_first_x, n_last_x: $n_last_x")
    if fix_time == false
		for i = n_first_x:n_last_x
	        if i != n_last_x
	            x0_vec[n_sv*n+i] = tofs[i]
	        end
		end
    end
	#println("x0_vec: $x0_vec")
	return x0_vec, n_first_x, n_last_x
end
