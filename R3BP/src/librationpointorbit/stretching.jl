"""
Computation of stretching directions for stable orbits
"""


"""
Get SVD decomposition of monodromy along LPO

# Arguments
    - `p::List`: arguments to be passed to ODE
    - `x0_lpo::Vector`: initial condition of LPO
    - `period::Float64`: period
    - `n::Int`: period
    - `direction::String`: period

# Returns
    - (tuple): x0_departures, Fs
"""
function get_lpo_svd(p, x0_lpo::Vector, period::Float64, n::Int, direction::String; kwargs...)
    # ---------- extract arguments ---------- #
    kwargs_dict = Dict(kwargs)

    # ODE settings
    reltol = assign_from_kwargs(kwargs_dict, :reltol, 1.e-12)
    abstol = assign_from_kwargs(kwargs_dict, :abstol, 1.e-12)
    method = assign_from_kwargs(kwargs_dict, :method, Vern7())

    # get state length
    nsv = length(x0_lpo)
    if nsv==4
        rhs!     = rhs_pcr3bp_sv!
        rhs_stm! = rhs_pcr3bp_svstm!
    elseif nsv==6
        rhs!     = rhs_cr3bp_sv!
        rhs_stm! = rhs_cr3bp_svstm!
    else
        error("x0 should be length 4 or 6!")
    end

    # propagate LPO by full period and get initial states
    prob_lpo = ODEProblem(rhs!, x0_lpo, period, p)
    sol_x0s = solve(prob_lpo, method, reltol=reltol, abstol=abstol, saveat=LinRange(0, period, n+1))

    x0_stm = vcat(x0_lpo, reshape(I(nsv), (nsv^2,)))[:]
    prob_lpo_stm = ODEProblem(rhs_stm!, x0_stm, period, p)
    Fs = []
    x0_departures = []

    @showprogress for i = 1:n
        x0_departure = sol_x0s.u[i][1:nsv]
        if i != 1
            x0_stm = vcat(x0_departure, reshape(I(nsv), (nsv^2,)))[:]
            prob_lpo_stm = remake(prob_lpo_stm; u0=x0_stm)
        end
        # solve problem
        sol = solve(prob_lpo_stm, method, reltol=reltol, abstol=abstol)
        monodromy = R3BP.get_stm(sol, nsv)
        # perform SVD & store
        if direction == "forward"
            push!(Fs, svd(monodromy))
        else
            push!(Fs, svd(inv(monodromy)))
        end
        push!(x0_departures, x0_departure)
    end
    return x0_departures, Fs
end
