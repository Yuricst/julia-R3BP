"""
Search through grid in BCR4BP
"""

function find_perigee(sol, mu)
    r1min = 100.0
    for u in sol.u
        r1 = sqrt((u[1]+ mu)^2 + u[2]^2 + u[3]^2)
        if r1 < r1min
            r1min = r1
        end
    end
    return r1min
end


function find_perigee_tof(sol, mu)
    r1min, tof = 100.0, 0.0
    for (idx,u) in enumerate(sol.u)
        r1 = sqrt((u[1]+ mu)^2 + u[2]^2 + u[3]^2)
        if r1 < r1min
            r1min = r1
            tof = sol.t[idx]
        end
    end
    return r1min, tof
end



function get_callbacks(params)
    # affects
    affect_terminate!(integrator) = terminate!(integrator)
    affect_nothing! = function (integrator)
        return
    end

	# callbacks
	rmax = 2.0e6 / params.lstar   # apogee max
	function condition_rmax(u,t,integrator)
	    return norm(u[1:3]) - rmax
	end
	cb_rmax = ContinuousCallback(condition_rmax, affect_terminate!, rootfind=false);

	r1min = 300/params.lstar
	function condition_r1min(u,t,integrator)
	    return sqrt((u[1]+ params.mu)^2 + u[2]^2 + u[3]^2) - r1min
	end
	cb_r1min = ContinuousCallback(condition_r1min, affect_terminate!, rootfind=false);

    function condition_r1apse(u,t,integrator)
        return (u[1]+ params.mu)*u[4] + u[2]*u[5] + u[3]*u[6]
    end
    cb_r1apse = ContinuousCallback(condition_r1apse, affect_nothing!, rootfind=true)

	# collect callbacks
	cbs = CallbackSet(cb_rmax, cb_r1min, cb_r1apse);
	return cbs
end


"""
Zoom in on grid in terms of t0
"""
function zoom_grid(params, x0, tof, cbs, n_init::Int=30, rp_threshold::Real=0.75, Δt0::Real=π/12)
    t0s = LinRange(0, 2π, n_init+1)[1:n_init]
    
    # create base ODE
    p_ode = (params.mu, params.μ_3, 0.0, params.a, params.ω_s)
    prob_base = ODEProblem(R3BP.rhs_bcr4bp_sv!, x0, (0.0, tof), p_ode, callback=cbs,
        method=Tsit5(), reltol=1e-12, abstol=1e-12, save_everystep=false)

    perigees = Float64[]
    sims = []
    for t0 in t0s
        # p[1] = μ, p[2] = μ_3, p[3] = t0, p[4] = a, p[5] = ω_s
        p_ode = (params.mu, params.μ_3, t0, params.a, params.ω_s)
        _prob = remake(prob_base, p=p_ode)
        sol = @suppress_err solve(_prob)
        # find perigee
        push!(perigees, find_perigee(sol, params.mu))
    end
    
    # iterate through perigee and find local minima
    t0_interest = Real[]
    Δdir = 0.0
    for idx = 1:length(perigees)-1
        if (perigees[idx+1] - perigees[idx])*Δdir < 0.0 && perigees[idx] < rp_threshold
            push!(t0_interest, t0s[idx])
        end
        Δdir = perigees[idx+1] - perigees[idx]
    end
    
    # research through regions of interest
    _perigees = Float64[]
    _t0s = zeros(2n_init)
    for (idx,t0) in enumerate(t0_interest)
        _t0s_iter = LinRange(t0-Δt0, t0+Δt0, n_init+1)[1:n_init]
        _t0s[(idx-1)*n_init+1:idx*n_init] = _t0s_iter
        for t0 in _t0s_iter
            # p[1] = μ, p[2] = μ_3, p[3] = t0, p[4] = a, p[5] = ω_s
            p_ode = (params.mu, params.μ_3, t0, params.a, params.ω_s)
            _prob = remake(prob_base, p=p_ode)
            sol = @suppress_err solve(_prob)
            # find perigee
            push!(_perigees, find_perigee(sol, params.mu))
        end
    end
    return _t0s, _perigees, prob_base
end



"""
Conduct root-solving against t0 to match perigee radius
"""
function rootsolve_t0(params, prob_base, t0s, perigees, x0, tof, radius_target, cbs, tol::Real=1e-6, verbose::Bool=false)
    # base ODE
    p_ode = (params.mu, params.μ_3, 0.0, params.a, params.ω_s)
    
    # residual function
    res_func = function (t0, full_return::Bool=false)
        _prob = remake(prob_base; p=(params.mu, params.μ_3, t0, params.a, params.ω_s))
        sol = @suppress_err solve(_prob);
        if full_return == false
            return find_perigee(sol, params.mu) - radius_target
        else
            return find_perigee(sol, params.mu), sol
        end
    end
    
    # find places to iterate
    lets = Vector[]
    for idx = 1:length(perigees)-1
        if (perigees[idx+1] - radius_target)*(perigees[idx] - radius_target) < 0.0
            # solve bisection problem
            t0_solved = find_zero(res_func, (t0s[idx], t0s[idx+1]), Bisection())
            rp, sol = res_func(t0_solved, true)
            _, tof = find_perigee_tof(sol, params.mu)
            if abs(rp - radius_target) < tol
                push!(lets, [rp, t0_solved, tof])
                if verbose
	            	@printf("LET t0 = %3.4f, rp = %4.4f\n", t0_solved*180/π, rp*params.lstar)
	            end
            end
        end
    end
    
    return lets
end



"""
Grid-search LET from perilune state via grid-search and bisection method in BCR4BP

# Arguments
- `params::BCR4BP_param`: BCR4BP parameters from `R3BP.get_bcr4bp_param(399, 301)`
- `x0::Vector`: perilune state
- `tof::Real`: time of flight; set to negative for Earth -> perilune transfer
- `radius_target::Real`: target radius at Earth
- `n::Int`: number of points on grid of t0
- `verbose::Bool`: flag on printing info
- `full_return::Bool`: whether to return all information

# Returns
- List of LET info, consisting of [perigee, t0, tof] for input x0
"""
function grid_search_let(
	params,
	x0::Vector,
	tof::Real,
	radius_target::Real,
	n::Int=40,
	tol::Real=1e-6,
	verbose::Bool=true,
	full_return::Bool=true,
	rp_threshold::Real=0.75,
	Δt0::Float64=deg2rad(15),
)

	# callbacks
	cbs = get_callbacks(params)

	# grid-search
	if verbose
		println("Constructing grid search with n = $n")
	end
	t0_interest, _perigees, prob_base = zoom_grid(params, x0, tof, cbs, n, rp_threshold, Δt0);

	# plot
	if full_return
		pt0 = plot(frame_style=:box, size=(600,400), ylabel="Perigee, km",
	    yscale=:log10, legend=:bottom)
		plot!(pt0, t0_interest, _perigees*params.lstar, marker=:circle)
		hline!(pt0, [radius_target*params.lstar,])
	end

	# root-solve
	lets = rootsolve_t0(params, prob_base, t0_interest, _perigees, x0, tof, radius_target, cbs, tol, verbose);
	if full_return
		return lets, pt0, cbs
	else
		return lets
	end
end