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
function check_apsis(sv, mu::Float64)
    return dot([sv[1] - (1-mu), sv[2]], sv[3:4])
end


"""
Main function
"""
function main(params, thetas, ras, verbose=false)
	# set-up Sun-Earth system
	mu = params.mu
	println("mu: $mu")

	# ---------- callbacks ---------- #
	affect!(integrator) = terminate!(integrator)

	# callback event based on having perigee around 0.9~1.1 * sma of the moon
	moon_sma = 384748.0 / params.lstar
	function condition_perigee_radius(u,t,integrator)
	    r_local = sqrt((u[1] - (1-mu))^2 + u[2]^2)
	    if 0.9moon_sma < r_local < 1.1moon_sma
	        return check_apsis(u[1:4], mu)  # when hitting apsis
	    else
	        return NaN
	    end
	end
	cb1 = ContinuousCallback(condition_perigee_radius, affect!, abstol=1.0e-12,reltol=1.0e-12)

	# callback event based on intersecting with n*earth radius with vr<0
	nr_earth = moon_sma/2  # 10*6378.0 / params.lstar
	function condition_nearth_Intersect(u,t,integrator)
		vr_sign = check_apsis(u, mu)
		if vr_sign < 0.0
			r_local = sqrt((u[1] - (1-mu))^2 + u[2]^2)
			return r_local - nr_earth
		else
			return NaN
		end
	end
	cb2 = ContinuousCallback(condition_nearth_Intersect, affect!, abstol=1.0e-12,reltol=1.0e-12)

	# callback event when leaving n*earth_SOI
	nearth_soi = 3*0.929e6 / params.lstar
	function condiction_leave_nSOI(u,t,integrator)
		r_local = sqrt((u[1] - (1-mu))^2 + u[2]^2)
		return nearth_soi - r_local
	end
	cb3 = ContinuousCallback(condiction_leave_nSOI, affect!, abstol=1.0e-12,reltol=1.0e-12)

	# create call back set
	cbs = CallbackSet(cb1, cb2, cb3)

	# ---------- setup initial guess ---------- #
	# perigee km
	rp = (6378.0 + 185.0)/params.lstar

	# setup storage
	x0s = []
	tfs = []
	sim_info = []

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

	# ---------- thruster parameter ---------- #
	mstar = 4100.0       # 4100 in kg for bepi-colombo
	tmax_newtons = 0.0   # N, used to be 0.4
	tmax = tmax_newtons *(1/mstar)*(params.tstar^2/(1e3*params.lstar))
	isp = 3500.0   # seconds
	mdot = 0.0  #tmax_newtons/(isp*9.80665) *(params.tstar/mstar)
	println("tmax: $tmax")
	println("mdot: $mdot")

	# ---------- setup on ODE problem ---------- #
	# tolerance of ODE problem
	reltol = 1.0e-13
	abstol = 1.0e-13

	# set-up initial problem
	out_term = []
	for (idx_ensemble, τ_iter) in enumerate([-1.0, 1.0])
		println("\nEnsemble sim # $idx_ensemble ..... using τ = $τ_iter")
		#τ = -1.0
		p = (mu, τ_iter, mdot, tmax)
		m0 = 1.0
		prob = ODEProblem(R3BP.rhs_pcr3bp_thrust_m1dir!, vcat(x0s[1], m0), (0.0, 2.0*tfs[1]), p)

		# ---------- ensemble simulation ---------- #
		function prob_func(prob, i, repeat)
			if verbose==true
				print("\rproblem # $i / $nic")
			end
		    remake(prob, u0=vcat(x0s[i], m0), tspan=(0.0, 2.0tfs[i]))
		end

		ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
		sim = solve(ensemble_prob, Tsit5(), EnsembleThreads(), 
			trajectories=length(x0s), callback=cbs, reltol=reltol, abstol=abstol);

		# ---------- post-process result ---------- #
		# create storage for terminated solutions
		for (idx, sol) in enumerate(sim)
		    if sol.retcode == :Terminated
		    	# check which event terminated the propagation
		    	if isnan(condition_perigee_radius(sol.u[end], 0.0, 0.0)) == false
			        push!(out_term, Dict(
			            "x0" => sol.u[1],
			            "xf" => sol.u[end],
			            "tf" => sol.t[end],
			            "theta" => sim_info[idx]["theta"],
			            "rp" => sim_info[idx]["rp"],
			            "ra" => sim_info[idx]["ra"],
			            "τ" => τ_iter,
			            "p" => p,
			            "m0" => m0,
			            "mf" => sol.u[end][end],
			            "mstar" => mstar,
			            "isp" => isp, 
			            "mdot" => mdot,
			            "tmax" => tmax,
			            "mstar" => mstar,
			        ))
			    end
		    end
		end
	end
	n_found = length(out_term)
	println("\nFound $n_found solutions!")

	# export data
	timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMM")
	flename = "se_enhanced_" * timestamp
	flepath = "./let-SunEarthSystem-data/" * flename * ".json"
	open(flepath,"w") do f
	    JSON.print(f, out_term)
	end
	print("Save data at: ")
	println(flepath)
	return
end


# ----------------------------------------- #
# run main analysis
params = R3BP.get_cr3bp_param(10, 399)
n_theta = 360
n_ra = 50
thetas = LinRange(0.0, 2π, n_theta+1)[1:end-1]
ras = LinRange(1.0e6/params.lstar, 2.0e6/params.lstar, n_ra)

main(params, thetas, ras)
println("Done!")
