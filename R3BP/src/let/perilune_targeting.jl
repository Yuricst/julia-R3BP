"""
Functions for perilune targeting via interpolation
"""


"""
Process poincare section to separate top and bottom surfaces
"""
function separate_top_bottom(xs, ys)
    n = size(xs,1)
    direction = 0
    flipped = Int[]
    for i = 1:n-1
        if direction * sign(xs[i+1] - xs[i]) < 0
            push!(flipped, i)
        end
        direction = sign(xs[i+1] - xs[i])  # update
    end
    if length(flipped)==1
        side1 = hcat(xs[1:flipped[1]], ys[1:flipped[1]])
        side2 = hcat(xs[flipped[1]+1:end], ys[flipped[1]+1:end])

    elseif length(flipped)==2
        side1 = hcat(
            vcat(xs[flipped[2]+1:end], xs[1:flipped[1]])[:],
            vcat(ys[flipped[2]+1:end], ys[1:flipped[1]])[:]
        )
        side2 = hcat(xs[flipped[1]+1:flipped[2]], ys[flipped[1]+1:flipped[2]])
    end
    # order them
    side1 = side1[sortperm(side1[:,1]), :]
    side2 = side2[sortperm(side2[:,1]), :]
    return flipped, side1, side2
end



"""
Interpolate Poincare-Section

Example: `ps = hcat([sol.u[end] for sol in sim]...)` where `sim` is the manifold.
"""
function interpolate_ps(ps, i_strip::Int=2, n_interval::Int=20, n_strip::Int=8)
	strips_per_state = Dict()
	anchors = Float64[]
	for i = 1:6
	    flipped, side1, side2 = separate_top_bottom(ps[i_strip, :], ps[i, :]);

	    # construct spline & evaluate
	    anchors = LinRange(side1[1,1], side1[end,1], n_interval)
	    spl1 = Spline1D(side1[:,1], side1[:,2])
	    y1s = evaluate(spl1, anchors)
	    
	    spl2 = Spline1D(side2[:,1], side2[:,2])
	    y2s = evaluate(spl2, anchors)
	    
	    # get strips
	    strips = Dict()
	    for j = 1:n_strip
	        strips[j] = (y2s-y1s) * j/(n_strip+1) + y1s
	    end
	    strips_per_state[i] = strips
	end

	strips_list = []
	for j = 1:n_strip
	    push!(strips_list,
	        [
	            strips_per_state[1][j],
	            strips_per_state[2][j],
	            strips_per_state[3][j],
	            strips_per_state[4][j],
	            strips_per_state[5][j],
	            strips_per_state[6][j],
	        ]
	    )
	end
    return anchors, strips_list
end


"""
Function to find perilune
"""
function find_perilune(sol, mu::Real, get_state::Bool=false)
    r2min = 100.0
    u_min = [1.0,0.0,0.0,0.0,0.0,0.0]
    for u in sol.u
        r2 = sqrt((u[1]-(1-mu))^2 + u[2]^2 + u[3]^2)
        if r2 < r2min
            r2min = r2
            u_min = u
        end
    end
    if get_state == false
        return r2min
    else
        return r2min, u_min
    end
end


"""
Root-solve strip of Poincare-section for intersection
"""
function get_strips(
    strips_list, mu::Real, r2_threshold_min::Real, tf_fwd::Real=5.0,
    method=Tsit5(), reltol=1e-12, abstol=1e-12,
) 
	# create events
	condition_r2hit = function (u,t,integrator)
	    r2 = sqrt((u[1]-(1-mu))^2 + u[2]^2 + u[3]^2)
	    return r2 - r2_threshold_min
	end
	affect1!(integrator) = terminate!(integrator)
	cb1 = ContinuousCallback(condition_r2hit, affect1!);

	condition_apse = function (u,t,integrator)
	    return (u[1]-(1-mu))*u[4] + u[2]*u[5] + u[3]*u[6]
	end
	affect2! = function (integrator) end
	cb2 = ContinuousCallback(condition_apse, affect2!);

	cbs = CallbackSet(cb1,cb2);

    # create base ODE problem
    prob_base = ODEProblem(
        R3BP.rhs_cr3bp_sv!, [1.1,0.0,0.0,0.0,1.0,0.0], tf_fwd, [mu,],
        method=method, reltol=reltol, abstol=abstol, callback=cbs,
    );

    # prepare strips
    strip_traj_list = []
	perilunes_per_strip = []
	@showprogress for _strip in strips_list  # length == n_strip
	    traj_per_strip = []
	    perilunes = Real[]
	    for j = 1:length(strips_list[1][1])  # n_interval
	        _x0 = [
	            _strip[1][j], _strip[2][j], _strip[3][j],
	            _strip[4][j], _strip[5][j], _strip[6][j],
	        ]
	        _prob = remake(prob_base; u0=_x0)
	        sol = @suppress_err solve(_prob)
	        push!(traj_per_strip, sol)
	        push!(perilunes, find_perilune(sol, mu))
	    end
	    
	    push!(strip_traj_list, traj_per_strip)
	    push!(perilunes_per_strip, perilunes)
	end

	return strips_list, perilunes_per_strip
end



"""
Interpolate strip to find targeting trajectory
"""
function interpolate_strip(
    mu::Real, anchors::Union{Vector,LinRange}, strip_states::Vector, perilunes::Vector, 
    target_radius::Real, i_strip::Int, r2_threshold_min, tf_fwd::Real=5.0,
    verbose::Bool=false, prob_base=nothing,
    method=Tsit5(), reltol=1e-12, abstol=1e-12,
)
    # interpolate along strip
    spl_state = [Spline1D(anchors, strip_states[i_state]) for i_state=1:6]  # state at manifold

	# create events
	condition_r2hit = function (u,t,integrator)
	    r2 = sqrt((u[1]-(1-mu))^2 + u[2]^2 + u[3]^2)
	    return r2 - r2_threshold_min
	end
	affect1!(integrator) = terminate!(integrator)
	cb1 = ContinuousCallback(condition_r2hit, affect1!);

	condition_apse = function (u,t,integrator)
	    return (u[1]-(1-mu))*u[4] + u[2]*u[5] + u[3]*u[6]
	end
	affect2! = function (integrator) end
	cb2 = ContinuousCallback(condition_apse, affect2!);

	cbs = CallbackSet(cb1,cb2);

    if isnothing(prob_base)
	    # create base ODE problem
	    prob_base = ODEProblem(
	        R3BP.rhs_cr3bp_sv!, [1.1,0.0,0.0,0.0,1.0,0.0], tf_fwd, [mu,],
	        method=method, reltol=reltol, abstol=abstol, callback=cbs,
	    );
	end
    
    residual = function (anchor_var, get_state::Bool=false)
        #println("anchor_var: $anchor_var")
        _x0_interp = [evaluate(spl_state[j], anchor_var) for j = 1:6]
        _prob = remake(prob_base; u0=_x0_interp)
        sol_anchor = @suppress_err solve(_prob)
        if get_state == true
            return find_perilune(sol_anchor, mu, get_state)
        else
            return find_perilune(sol_anchor, mu, get_state) - target_radius
        end
    end
    
    # iterate through sampled perilune-going locations
    state_target = Vector[]
    for i = 1:length(perilunes)-1
        if (perilunes[i]-target_radius) * (perilunes[i+1]-target_radius) < 0
        	if verbose == true
	            println("Detected change in sign at index $i")
	            #println("anchors[i]   => ", strip_states[i_strip][i])
	            #println("anchors[i+1] => ", strip_states[i_strip][i+1])
	            #println("Stored: ", perilunes[i]-target_radius)
	            println("Computed: ", residual(strip_states[i_strip][i]))
	            println("Computed: ", residual(strip_states[i_strip][i+1]))
	        end
            res_anchor = find_zero(
                residual, 
                (strip_states[i_strip][i], strip_states[i_strip][i+1]),
                Bisection(), atol=1e-8, rtol=1e-8
            )
            # store
            r2min, state_llo = residual(res_anchor, true)
            if abs(r2min - target_radius) < 1e-8
            	push!(state_target, state_llo)
            end
            if verbose == true
            	println(" ..... residual: ", r2min - target_radius)
            end
        end
    end
    return state_target
end



function lpo2llo_target(
	mu::Real,
	X0_lpo,
	period::Real,
	i_strip    = 2,
	n_interval = 30,
	n_strip    = 50,
	tf_fwd     = 5.0,
	r2_threshold_min = 0.004518,
	target_radius    = 0.005819,
	N_manif::Int=50,
	tf_manif::Real=-10.0,
	manif_cb=nothing,
	verbose::Bool=true,
)
	if isnothing(manif_cb)
		r2_threshold = 1.0
		function condition(u,t,integrator)
		    r2 = sqrt((u[1]-(1-mu))^2 + u[2]^2 + u[3]^2)
		    return r2 - r2_threshold
		end
		affect!(integrator) = terminate!(integrator)
		manif_cb = ContinuousCallback(condition,affect!);
	end
	# generate manifolds
	if verbose
		println("Generating manifolds...")
	end
	stability = true
	direction = "positive"
	sim = R3BP.get_manifold(mu, X0_lpo, period, tf_manif, stability, N_manif, direction, manif_cb);

	# get poincare section
	if verbose
		println("Preparing Poincare Section...")
	end
	ps = hcat([sol.u[end] for sol in sim]...)

	# process PS
	anchors, strips_list = R3BP.interpolate_ps(ps, i_strip, n_interval, n_strip);
	strips_list, perilunes_per_strip = R3BP.get_strips(strips_list, mu, r2_threshold_min, tf_fwd);

	# get states at LLO
	states_llo = []
	@showprogress for istrip = 1:length(perilunes_per_strip)
	    _states_llo = R3BP.interpolate_strip(
	        mu, anchors, strips_list[istrip], perilunes_per_strip[istrip], target_radius,
	        i_strip, r2_threshold_min, tf_fwd, verbose
	    );
	    if length(states_llo) == 0
	        states_llo = _states_llo
	    else
	        states_llo = vcat(states_llo, _states_llo)
	    end
	end
	return states_llo
end