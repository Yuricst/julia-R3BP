"""
Multiple-shooting differential correction algorithm
"""


struct MultipleShootingOptions
	verbose::Bool
	n_sv::Int
	n::Int
	fix_time::Bool
	periodic::Bool
	fix_x0::Bool
	fix_xf::Bool
end


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

	# construct struct of options
	MSOptions = MultipleShootingOptions(verbose, n_sv, n, fix_time, periodic, fix_x0, fix_xf)

	# ------------- check for exceptions ------------- #
	multiple_shooting_exceptions(tofs, fix_time, rhs!, MSOptions)

    # ------------- storage setups ------------- #
	__print_verbose(" -------------- Multiple shooting algorithm -------------- ", verbose)
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

    # initialize free-parameter guess vector
	x0_vec, svf, n_propagated = initialize_multiple_shooting_x(x0s, tofs, MSOptions)

	# initialize error vector
	ferr = initialize_multiple_shooting_f(MSOptions)

    # initialize DF matrix
	df = initialize_multiple_shooting_df(length(ferr), length(x0_vec), MSOptions)#n_sv, n_propagated, periodic)

    # initialize storages
    convflag = 0    # convergence flag

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
            multiple_shooting_propagate_sequential!(
				df, svf, prob_stm, rhs!, p, n, n_sv, n_propagated, x0_vec, tofs, fix_time
			)
        else
            multiple_shooting_propagate_ensemble!(
				df, svf, prob_stm, rhs!, p, prob_func, ensemble_method, n, n_sv,
				n_propagated, x0_vec, tofs, fix_time
			)
        end

        # update error vector F
		multiple_shooting_update_f!(n, n_sv, ferr, svf, x0_vec, periodic, n_propagated)

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
        #multiple_shooting_update_x!(x0_vec, ferr, df, n, n_sv, fix_time, fix_x0, fix_xf)
		apply_update_law!(x0_vec, ferr, df)
    end
	__print_verbose(" -------------------------------------------------------- ", verbose)
    return x0_vec, tofs, convflag
end


"""
    multiple_shooting_propaagte_sequential!(df::Matrix{Float64}, svf::Vector{Float64}, prob_stm::ODEProblem, n::Int, n_sv::Int, n_propagated::Int, x0_vec::Vector{Float64}, tofs::Array, fix_time::Bool, method, reltol::Float64, abstol::Float64)

Propagates nodes and mutates DF and xfs; solved sequentially using `remake()`
"""
function multiple_shooting_propagate_sequential!(
	df::Matrix{Float64}, svf::Vector{Float64}, prob_stm::ODEProblem,
	rhs!, p, n::Int, n_sv::Int, n_propagated::Int, x0_vec::Vector{Float64},
	tofs::Array, fix_time::Bool,
)
    # re-make ODE problem and integrate until tofs[ix]
    for i = 1:n_propagated

		# integrate state-vector and STM
        x0_stm = vcat(x0_vec[(i-1)*n_sv + 1 : i*(n_sv)], reshape(I(n_sv), (n_sv*n_sv,)))[:]
        if fix_time == true
            tprop = tofs[i]
        else
            tprop = x0_vec[n_sv*n+i]
        end
        _prob = remake(prob_stm; tspan=(0.0, tprop), u0=x0_stm)
        sol = DifferentialEquations.solve(_prob)

        # append to vector of propagation results
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
    multiple_shooting_propaagte_ensemble!(df::Matrix{Float64}, svf::Vector{Float64}, prob_stm::ODEProblem, n::Int, n_sv::Int, n_propagated::Int, x0_vec::Vector{Float64}, tofs::Array, fix_time::Bool, method, reltol::Float64, abstol::Float64)

Propagates nodes and mutates DF and svf; solved sequentially using `ensembleSimulation()`
"""
function multiple_shooting_propagate_ensemble!(
	df::Matrix{Float64}, svf::Vector{Float64}, prob_stm::ODEProblem,
	rhs!, p, prob_func, ensemble_method, n::Int, n_sv::Int, n_propagated::Int,
	x0_vec::Vector{Float64}, tofs::Array, fix_time::Bool)
	# solve ensemble problem
	ensemble_prob = EnsembleProblem(prob_stm, prob_func=prob_func)
	sim = solve(ensemble_prob, method, ensemble_method; trajectories=n_propagated);

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
function multiple_shooting_update_f!(n::Int, n_sv::Int, ferr::Vector{Float64}, svf::Vector{Float64}, x0_vec::Vector{Float64}, periodic::Bool, n_propagated::Int)
	# update error vector F
	ferr[1:n_sv*n_propagated] = svf - x0_vec[n_sv+1:n_sv*(n_propagated+1)]
	if periodic == true
		ferr[end-n_sv+1:end] = x0_vec[n_sv*n-n_sv+1:n_sv*n] - x0_vec[1:n_sv]
	end
	return
end


# """
#     multiple_shooting_update_x!(x::Vector{Float64}, ferr::Vector{Float64}, df::Matrix{Float64})
#
# Apply multiple-shooting update to mutate x
# """
# function multiple_shooting_update_x!(x::Vector{Float64}, ferr::Vector{Float64}, df::Matrix{Float64}, n::Int, n_sv::Int, fix_time, fix_x0, fix_xf)
# 	# reduce x acoordingly
# 	if fix_x0 == false && fix_xf == false
# 		x_red = x
# 		df_red = df
#
# 	elseif fix_x0 == false && fix_xf == true
# 		if fix_time == true
# 			x_red  = x[1:end-n_sv]
# 			df_red = df[:, 1:(n-1)*n_sv]
# 		else  # fix_time == false
# 			x_red = vcat(x[1:(n-1)*n_sv], x[n*n_sv+1:end])
# 			df_red = hcat(df[:, 1:(n-1)*n_sv], df[:, n*n_sv+1:end])
# 		end
#
# 	elseif fix_x0 == true  && fix_xf == false
# 		x_red = x[n_sv+1:end]
# 		df_red= df[:, n_sv+1:end]
#
# 	elseif fix_x0 == true && fix_xf == true
# 		if fix_time == true
# 			x_red = x[n_sv+1:end-n_sv]
# 			df_red = df[:, n_sv+1:(n-1)*n_sv]
# 		else  # fix_time == false
# 			x_red = vcat(x[n_sv+1:(n-1)*n_sv], x[n*n_sv+1:end])
# 			df_red = hcat(df[:, n_sv+1:(n-1)*n_sv], df[:, n*n_sv+1:end])
# 		end
# 	end
#
# 	print("ferr: "); println(length(ferr))
# 	print("x_red: "); println(length(x_red))
# 	print("df_red: "); println(size(df_red))
#
# 	# apply mutliple-shooting update law
# 	apply_update_law!(x_red, ferr, df_red)
#
# 	# re-construct x
# 	if fix_x0 == false && fix_xf == false
# 		x[:] = x_red
# 	elseif fix_x0 == false && fix_xf == true
# 		x[:] = vcat(x_red, x[end-n_sv+1:end])
# 		#xred =  x[1:end-n_sv]
# 	elseif fix_x0 == true  && fix_xf == false
# 		#xred = x[n_sv+1:end]
# 		x[:] = vcat(x[1:n_sv], x_red)
#
# 	elseif fix_x0 == true && fix_xf == true
# 		x[:] = vcat(x[1:n_sv], x_red, x[end-n_sv+1:end])
# 	end
#     return
# end


"""
	apply_update_laws!(x::Vector{Float64}, ferr::Vector{Float64}, df::Matrix{Float64})

Apply update-laws to mutate x (or reduced x)
"""
function apply_update_law!(x::Vector{Float64}, ferr::Vector{Float64}, df::Matrix{Float64})
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
	multiple_shooting_exceptions(tofs::Array, fix_time::Bool, rhs!, MSOptions::MultipleShootingOptions)

Check multiple shooting exceptions
"""
function multiple_shooting_exceptions(tofs::Array, fix_time::Bool, rhs!, MSOptions::MultipleShootingOptions)
	# exception if rhs! is not provided but fix_time == false
	if MSOptions.fix_time == false && isnothing(rhs!) == true
        error("Provide callable rhs! if fix_time==false")
    end

	# exception if n < 3 but x0 or xf is fixed
	if (MSOptions.fix_x0 == true && MSOptions.fix_xf == false) || (MSOptions.fix_x0 == false && MSOptions.fix_xf == true)
		if MSOptions.n < 3
			error("Provide at least 3 nodes if fixing x0 or xf")
		end

	# exception if n < 4 but x0 and xf is fixed
	elseif MSOptions.fix_x0 == true && MSOptions.fix_xf == true
		if MSOptions.n < 4
			error("Provide at least 4 nodes if fixing x0 or xf")
		end
	end

	# exception if length(tof) != length(x0) - 1
	if length(tofs) != MSOptions.n-1
        error("Must follow length(x0s)-1 == length(tofs)")

	# exception if less than 2 nodes supplied
	elseif MSOptions.n == 1
        error("Provide at least 2 nodes for multiple-shooting")
    end
end


"""
	initialize_multiple_shooting_x(x0s::Array, n::Int, n_sv::Int, fix_time::Bool, fix_x0::Bool, fix_xf::Bool)

Initialize initial-guess vector and propagation output vector
"""
function initialize_multiple_shooting_x(x0s::Array, tofs::Array, MSOptions::MultipleShootingOptions)
	# ------------- get number of free x's ------------- #
	n_propagated = MSOptions.n-1  # number of nodes to be propagated

	# ------------- initialize space ------------- #
	# nodes and decision vector
	if MSOptions.fix_time == true
        x0_vec = zeros(MSOptions.n_sv*MSOptions.n)   # length = (sv-elements)*(number of nodes)
    else
        x0_vec = zeros(MSOptions.n_sv*MSOptions.n + MSOptions.n-1)  # length = (sv-elements + 1)*(number of nodes) - 1
    end
	# propaated outcome of nodes
	svf = zeros(MSOptions.n_sv*(MSOptions.n-1))

	# ------------- append state-vector ------------- #
    # append state-vector into initial guess
    for i = 1:MSOptions.n
        x0_vec[(i-1)*MSOptions.n_sv + 1 : i*(MSOptions.n_sv)] = x0s[i]
    end

	# ------------- append time of flight ------------- #
    # append integration times into initial guess
    if MSOptions.fix_time == false
		for i = 1:MSOptions.n-1
			x0_vec[MSOptions.n_sv*MSOptions.n+i] = tofs[i]
		end
    end
	return x0_vec, svf, n_propagated
end


"""
	initialize_multiple_shooting_f(MSOptions::MultipleShootingOptions)

Initialize array for error vector F
"""
function initialize_multiple_shooting_f(MSOptions::MultipleShootingOptions)
	# initialize length of F vector
	n_f = MSOptions.n_sv*(MSOptions.n-1)
	if MSOptions.periodic == true
		n_f += MSOptions.n_sv
	end
	return zeros(n_f)
end



"""
	initialize_multiple_shooting_df(nf::Int, nx::Int, n_sv::Int, n_propagated::Int)

Initialize array for DF = dF/dX
"""
function initialize_multiple_shooting_df(nf::Int, nx::Int, MSOptions::MultipleShootingOptions)# n_sv::Int, n_propagated::Int, periodic::Bool)
	# initialize DF matrix
	df = zeros(nf, nx)

	# fill-in negative identity into DF matrix
    for j = 1:MSOptions.n-1  #n_propagated
        df[(j-1)*MSOptions.n_sv+1:j*MSOptions.n_sv, j*MSOptions.n_sv+1:(j+1)*MSOptions.n_sv] = -I(MSOptions.n_sv)
    end

    # fill-in identity to inital and final states
    if MSOptions.periodic == true
		df[(MSOptions.n-1)*MSOptions.n_sv+1:MSOptions.n*MSOptions.n_sv, 1:MSOptions.n_sv] = -I(MSOptions.n_sv)
		df[(MSOptions.n-1)*MSOptions.n_sv+1:MSOptions.n*MSOptions.n_sv, MSOptions.n_sv*(MSOptions.n-1)+1:MSOptions.n_sv*MSOptions.n] = I(MSOptions.n_sv)
    end

	# # fill-in identity to initial state
	# if fix_x0 == true
	# 	df[n_propagated*n_sv+1+n_f_additional:(n_propagated+1)*n_sv+n_f_additional, 1:n_sv] = I(n_sv)
	# 	n_f_additional += n_sv
	# end
	#
	# # fill-in identity to final state
	# if fix_xf == true
	# 	df[n_propagated*n_sv+1+n_f_additional:(n_propagated+1)*n_sv+n_f_additional, n_sv*n_propagated+1:n_sv*(n_propagated+1)] = I(n_sv)
	# end

	return df
end
