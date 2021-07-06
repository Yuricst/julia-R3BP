"""
"""

using LinearAlgebra
using DifferentialEquations

using Dates
using JSON

include("../R3BP/src/R3BP.jl")


"""
	get_initial_condition(mu::Float64, rp::Float64, ra::Float64, theta::Float64, m::Int=2, velocity_dir::String="positive")

Get initial condition, from Belbruno 2004
    
Args:
    mu (float): R3BP system parameter
    rp (float): periapsis radius
    ra (float): apoapsis radius
    theta (float): angle w.r.t. x-axis, radian
    m (int): center mass, 1 for m1, 2 for m2
    velocity_dir (str): "positive" or "negative"

Returns:
    (Array): initial condition
"""
function get_initial_condition(mu::Float64, rp::Float64, ra::Float64, theta::Float64, m::Int=2, velocity_dir::String="positive")
    # compute e
    e = (ra-rp)/(ra+rp)
    # compute v
    if m ==2
        v = sqrt( mu*(1+e)/rp )
    elseif m ==1
        v = sqrt( (1-mu)*(1+e)/rp ) 
    end
    # construct state
    if m==2
        x = (1 - mu) + rp*cos(theta)
    elseif m==1
        x = (-mu) + rp*cos(theta)   # FIXME?
    end
    y  =  rp*sin(theta)
    if velocity_dir=="positive"
        vx =  rp*sin(theta) - v*sin(theta)
        vy = -rp*cos(theta) + v*cos(theta)
    elseif velocity_dir=="negative"
        vx =  rp*sin(theta) + v*sin(theta)
        vy = -rp*cos(theta) - v*cos(theta)
    end
    return [ x, y, vx, vy ]
end


"""
	check_apsis(sv, mu::Float64, m::Int=2)

Check if apsis is hit for planar problem
"""
function check_apsis(sv, mu::Float64, m::Int=2)
    if m == 2
        x = sv[1] - (1-mu)
    end
    return dot([x, sv[2]], sv[3:4])
end


"""
Main function
"""
function main()
	# set-up Sun-Earth system
	params = R3BP.get_cr3bp_param(10, 399)
	mu = params.mu
	println("mu: $mu")

	# callback event based on having perigee around 0.9~1.1 * sma of the moon
	moon_sma = 384748.0 / params.lstar
	function condition(u,t,integrator)
	    r_local = sqrt((u[1] - (1-mu))^2 + u[2]^2)
	    if 0.9moon_sma < r_local < 1.1moon_sma
	        return check_apsis(u, mu)   # when hitting apsis
	    else
	        return NaN
	    end
	end
	affect!(integrator) = terminate!(integrator)
	cb = ContinuousCallback(condition,affect!)

	# tolerance of ODE problem
	reltol = 1.0e-13
	abstol = 1.0e-13

	# perigee km
	rp = (6378.0 + 185.0)/params.lstar

	# setup storage
	x0s = []
	tfs = []
	sim_info = []
	n = 360
	n_ra = 20

	thetas = LinRange(0.0, 2π, n+1)[1:end-1]
	ras = LinRange(0.8e6/params.lstar, 1.6e6/params.lstar, n_ra)

	for ra in ras
	    sma = (rp+ra)/2
	    period = period = 2π*sqrt((sma*params.lstar)^3/R3BP.get_gm("399")[1]) / params.tstar
	    for theta_iter in thetas
	        push!(x0s, get_initial_condition(mu, rp, ra, theta_iter))
	        push!(tfs, period)
	        push!(sim_info, Dict("theta"=>theta_iter, "rp"=>rp, "ra"=>ra))
	    end
	end

	# check number of initial guesses
	nic = length(x0s)
	ntf = length(tfs)
	println("Using $nic - $ntf initial conditions")

	# set-up initial problem
	prob = ODEProblem(R3BP.rhs_pcr3bp_sv!, x0s[1], (0.0, 2.0*tfs[1]), (mu))

	# ---------- ensemble simulation ---------- #
	function prob_func(prob, i, repeat)
	    remake(prob, u0=x0s[i], tspan=(0.0, 2.0tfs[i]))
	end
	ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
	sim = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories=length(x0s), callback=cb, reltol=reltol, abstol=abstol);

	# create storage for terminated solutions
	out_term = []
	for (idx, sol) in enumerate(sim)
	    if sol.retcode == :Terminated
	        push!(out_term, Dict(
	            "x0" => sol.u[1],
	            "xf" => sol.u[end],
	            "tf" => sol.t[end],
	            "theta" => sim_info[idx]["theta"],
	            "rp" => sim_info[idx]["rp"],
	            "ra" => sim_info[idx]["ra"],
	        ))
	    end
	end
	n_found = length(out_term)
	println("Found $n_found solutions!")

	# export data
	timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMM")
	flename = "se_ballistic_" * timestamp
	flepath = "./let-SunEarthSystem-data/" * flename * ".json"
	open(flepath,"w") do f
	    JSON.print(f, out_term)
	end
	print("Save data at: ")
	println(flepath)
end


# run main analysis
main()
println("Done!")