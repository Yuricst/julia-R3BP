"""
Function associated with manifold
"""


# -------------------------------------------------------------------------------- #
# methods associated with STM
"""
    get_stm(sol::ODESolution, nsv::Int)

Extract STM from solution of ODEProblem at final index
"""
function get_stm(sol, nsv::Int)
    return transpose( reshape(sol.u[end][nsv+1:end], (nsv, nsv)) )
end


"""
    get_stm(sol::ODESolution, nsv::Int, idx::Int)

Extract STM from solution of ODEProblem at index=idx
"""
function get_stm(sol, nsv::Int, idx::Int)
    return transpose( reshape(sol.u[idx][nsv+1:end], (nsv, nsv)) )
end


# -------------------------------------------------------------------------------- #
# methods associated with manifold computation
"""
    get_stable_unstable_eigvecs(λs, vs)

Identidy stable and unstable eigenvalue-eigenvector pair
"""
function get_stable_unstable_eigvecs(λs, vs)
    λs_real = []
    for (idx,λ) in enumerate(λs)
        if imag(λ) == 0
            push!(λs_real, [idx, real(λ)])
        end
    end
    λs_real = hcat(λs_real...)'
    idx_λmin = Int(λs_real[argmin(λs_real[:,2]), 1])
    idx_λmax = Int(λs_real[argmax(λs_real[:,2]), 1])
    # get real parts
    emin = real(λs[idx_λmin])
    emax = real(λs[idx_λmax])
    vmin = real(vs[:,idx_λmin])
    vmax = real(vs[:,idx_λmax])

    if (emin < 1.0) && (emax > 1.0)
        return emax, emin, vmax, vmin
    else
        return 1.0, 1.0, 0.0, 0.0
    end
    # idx_stb_unstb = [];
    # for idx in 1:length( λs )
    #     if abs(real(λs[idx]) - 1) > 1e-4 && imag(λs[idx])==0
    #         push!(idx_stb_unstb, idx);
    #     end
    # end
    #
    # if length(idx_stb_unstb) > 0
    #     if abs(λs[ idx_stb_unstb[1] ]) > abs(λs[ idx_stb_unstb[2] ])
    #         eig_unstb = real( λs[ idx_stb_unstb[1] ] );
    #         v_unstb   = real( vs[:, idx_stb_unstb[1] ] );
    #         eig_stb   = real( λs[ idx_stb_unstb[2] ] );
    #         v_stb     = real( vs[:, idx_stb_unstb[2] ] );
    #     else
    #         eig_unstb = real( λs[ idx_stb_unstb[2] ] );
    #         v_unstb   = real( vs[:, idx_stb_unstb[2] ] );
    #         eig_stb   = real( λs[ idx_stb_unstb[1] ] );
    #         v_stb     = real( vs[:, idx_stb_unstb[1] ] );
    #     end
    #     return eig_unstb, eig_stb, v_unstb, v_stb
    # else
    #     return 1.0, 1.0, 0.0, 0.0
    # end
end


"""
    get_eigenvector(monodromy, stable::Bool=true)

Get eigenvectors from monodromy matrix
"""
function get_eigenvector(monodromy, stable::Bool, verbosity::Int=0)
    λs = eigvals(monodromy)
    vs = eigvecs(monodromy)
    eig_unstb, eig_stb, v_unstb, v_stb = get_stable_unstable_eigvecs(λs, vs)
    nu = stability(eig_unstb, eig_stb)
    __print_verbosity("Linear stability ν = $nu", verbosity, 0)
    #@printf("Eigenvalue unstable: %1.4f, stable: %1.4f, product: %1.4f\n", eig_unstb, eig_stb, eig_unstb*eig_stb)
    if stable==true
        return v_stb
    else
        return v_unstb
    end
end


"""
    scale_ϵ(μ::Float64, x0, period::Float64, stable::Bool, monodromy, y0, lstar::Float64, relative_tol_manifold::Float64=0.1, absolute_tol_manifold_km::Float64=100.0)

Obtain linear perturbation ϵ magnitude for manifolds

# Arguments
    - `μ::Float64`: CR3BP parameter
"""
function scale_ϵ(μ::Float64, x0, period::Float64, stable::Bool, monodromy, y0, lstar::Float64,
    relative_tol_manifold::Float64=0.1, absolute_tol_manifold_km::Float64=100.0)
    if length(x0)==4
        idx_pos_last = 2
    elseif length(x0)==6
        idx_pos_last = 3
    else
        error("Expected x0 to be length 4 or 6!")
    end
    # define absolute tolerance in canonical units
    absolute_tol_manifold = absolute_tol_manifold_km/lstar

    # list of perturbation sizes to try
    perturbation_km_lst = [0.01 0.05 0.1 0.5 1.0 5.0 10.0 50.0 100.0 500.0 1000.0]
    perturbation_lst = perturbation_km_lst/lstar

    # compute eigenvectors error after 1 rev using STM
    if stable==true   # stable case
        ef = inv(monodromy) * y0
    else              # unstbale case
        ef = monodromy * y0
    end
    # compute position error
    efPosition = norm( ef[1:idx_pos_last] )

    # initialise epsilon with smallest epsilon
    ϵ = deepcopy( perturbation_lst[1] )

    for ϵ_try in perturbation_lst
        # predicted error based on STM
        predictedError = ϵ_try * efPosition

        # compute error by propagating one period
        x0_manifold = x0 + ϵ_try*y0#reshape(y0, (1,length(x0)))

        if stable==true   # stable case
            #prop_manif_out = propagate_cr3bp(μ, x0_manifold, -period)
            prob_estErr = ODEProblem(rhs_pcr3bp_sv!, x0_manifold, (0, -period), (μ),
                method=Tsit5(), reltol=1e-11, abstol=1e-11
            );
        else               # unstable case
            #prop_manif_out = propagate_cr3bp(μ, x0_manifold, period)
            prob_estErr = ODEProblem(rhs_pcr3bp_sv!, x0_manifold, (0, period), (μ),
                method=Tsit5(), reltol=1e-11, abstol=1e-11
            );
        end
        sol = solve(prob_estErr);

        # actual position error
        propPositionError = norm(sol.u[end][1:idx_pos_last] - x0[1:idx_pos_last])
        #la.norm( prop_manif_out["statef"][0:3] - stateP[0:3] )

        if abs( (propPositionError - predictedError)/propPositionError ) > relative_tol_manifold && abs( propPositionError - predictedError ) > absolute_tol_manifold
            break
        else
            # else store solution largest epsilon allowed
            ϵ = deepcopy( ϵ_try )
        end
    end
    return ϵ
end


## manifold functions

# struct WarmManifold
#     x0_ptb_vec::Array{Float64,1}
#     n::Int
#     nsv::Int
#     μ::Float64
#     rhs!
# end

"""
    get_manifold(
            μ::Float64,
            x0::Array,
            period::Float64,
            tf::Float64,
            stable::Bool;
            kwargs...)

Function to obtain manifold of LPO

# Arguments
    - `μ::Float64`:: R3BP parameter
    - `x0::Array`: LPO state
    - `period::Float64`: LPO period
    - `tf::Float64`: time of flight of manifold
    - `stable::Bool`: true for stable manifold, false for unstable manifold
    - `xdir::String`: direction of manifold, "positive" or "negative" along x-axis

# Returns
    `EnsembleSolution`: ODE solution of `n` discrete manifold branches
    `array`: x0_ptb_vec
"""
function get_manifold(
            μ::Float64,
            x0::Array,
            period::Float64,
            tf::Float64,
            stable::Bool;
            kwargs...)
    # ---------- extract arguments ---------- #
    kwargs_dict = Dict(kwargs)
    # main manifold options
    xdir     = assign_from_kwargs(kwargs_dict, :xdir, "positive")
    n        = assign_from_kwargs(kwargs_dict, :n, 50)
    callback = assign_from_kwargs(kwargs_dict, :callback, nothing)
    ϵ      = assign_from_kwargs(kwargs_dict, :ϵ, nothing)
    lstar    = assign_from_kwargs(kwargs_dict, :lstar, 1.0)
    relative_tol_manifold    = assign_from_kwargs(kwargs_dict, :relative_tol_manifold, 0.1)
    absolute_tol_manifold_km = assign_from_kwargs(kwargs_dict, :absolute_tol_manifold_km, 100.0)
    verbosity = assign_from_kwargs(kwargs_dict, :verbosity, 0)
    detailed_output = assign_from_kwargs(kwargs_dict, :detailed_output, false)
    ensemblealg = assign_from_kwargs(kwargs_dict, :ensemblealg, EnsembleThreads())

    # ODE settings
    reltol = assign_from_kwargs(kwargs_dict, :reltol, 1.e-12)
    abstol = assign_from_kwargs(kwargs_dict, :abstol, 1.e-12)
    method = assign_from_kwargs(kwargs_dict, :method, Tsit5())

    # warm-start options
    x0_ptb_vec = assign_from_kwargs(kwargs_dict, :wm, nothing)

    __print_verbosity("========== Manifold Setting ==========", verbosity, 0)
    __print_verbosity("   xdir: $xdir\n   Stable: $stable", verbosity, 0)

    # ---------- propagate c0 by one full period with STM ---------- #
    nsv = length(x0)
    if nsv==4
        rhs!     = rhs_pcr3bp_sv!
        rhs_stm! = rhs_pcr3bp_svstm!
    elseif nsv==6
        rhs!     = rhs_cr3bp_sv!
        rhs_stm! = rhs_cr3bp_svstm!
    else
        error("x0 should be length 4 or 6!")
    end

    # if no WarmManifold is provided, compute parameters
    if isnothing(x0_ptb_vec) == true

        # obtain STM history
        x0_stm = vcat(x0, reshape(I(nsv), (nsv^2,)))[:]
        prob_lpo = ODEProblem(rhs_stm!, x0_stm, period, (μ), method=method, reltol=reltol, abstol=abstol)
        sol = solve(prob_lpo, saveat=LinRange(0, period, n+1))

        # get monodromy matrix
        monodromy = get_stm(sol, nsv)

        # get eigenvectors at initial state
        y0 = get_eigenvector(monodromy, stable, verbosity)

        # define ϵ (linear perturbation)
        if isnothing(ϵ)
            ϵ = scale_ϵ(μ, x0, period, stable, monodromy, y0, lstar, relative_tol_manifold, absolute_tol_manifold_km)
        end
        # decide sign of ϵ based on xdir
        if cmp(xdir, "positive")==0
            if y0[1] < 0.0
                ϵ_corr = -abs(ϵ)
            else
                ϵ_corr =  abs(ϵ)
            end
        elseif cmp(xdir, "negative")==0
            if y0[1] > 0.0
                ϵ_corr = -abs(ϵ)
            else
                ϵ_corr =  abs(ϵ)
            end
        else
            error("xdir should be \"positive\" or \"negative\"")
        end
        __print_verbosity("Using linear perturbation ϵ = $ϵ", verbosity, 0)

        # ---------- construct and append initial condition ---------- #
        x0_ptb_vec = zeros(n*nsv)
        for i in 1:n
            # map eigenvector (y = stm*y0)
            y_transcribed = get_stm(sol, nsv, i) * y0
            # construct linearly perturbed state
            x0_ptb_vec[1+(i-1)*nsv : i*nsv] = sol.u[i][1:nsv] + ϵ_corr*y_transcribed/norm(y_transcribed)
        end
    else
        __print_verbosity("Warm-starting EnsembleProblem", verbosity, 0)
    end

    # define base ODE problem for manifold branch
    prob_branch = ODEProblem(rhs!, x0_ptb_vec[1:nsv], tf, (μ))

    # define remake function
    function prob_func(prob, i, repeat)
        remake(prob, u0=x0_ptb_vec[1+(i-1)*nsv : i*nsv])
    end

    # construct Ensemble Problem
    ensemble_prob = EnsembleProblem(prob_branch, prob_func=prob_func)

    # return output of EnsembleProblem and perturbed ic vector
    sim = solve(
        ensemble_prob, method, ensemblealg, trajectories=n,
        #callback=callback, method=method, reltol=reltol, abstol=abstol
    )
    if detailed_output == false
        return sim, x0_ptb_vec
    else
        return sim, monodromy
    end
end



"""
    get_manifold(
        x0_ptb_vec::Array{Float64,1},
        μ::Float64,
        nsv::Int,
        tf::Float64,
        kwargs...)

Function to obtain manifold of LPO.
This dispatch utilizes pre-computed perturbation x's.

# Arguments
    - `x0_ptb_vec::Array{Float64,1}`: array of branch-points x0's, length (num. of branches)*(state-vector length)
    - `μ::Float64`:: R3BP parameter
    - `nsv::Int`: number of elements in state-vector
    - `tf::Float64`: time of flight of manifold

# Returns
    `EnsembleSolution`: ODE solution of `n` discrete manifold branches
"""
function get_manifold(
    x0_ptb_vec::Array{Float64,1},
    μ::Float64,
    nsv::Int,
    tf::Float64,
    kwargs...
)
    # ---------- extract arguments ---------- #
    kwargs_dict = Dict(kwargs)
    # main manifold options
    callback = assign_from_kwargs(kwargs_dict, :callback, nothing)
    # ODE settings
    reltol = assign_from_kwargs(kwargs_dict, :reltol, 1.e-12)
    abstol = assign_from_kwargs(kwargs_dict, :abstol, 1.e-12)
    method = assign_from_kwargs(kwargs_dict, :method, Tsit5())
    verbosity = assign_from_kwargs(kwargs_dict, :verbosity, 0)

    __print_verbosity("Warm-starting EnsembleProblem", verbosity, 0)

    # get problem rhs function and dimension
    if nsv==4
        rhs!     = rhs_pcr3bp_sv!
        rhs_stm! = rhs_pcr3bp_svstm!
    elseif nsv==6
        rhs!     = rhs_cr3bp_sv!
        rhs_stm! = rhs_cr3bp_svstm!
    else
        error("x0 should be length 4 or 6!")
    end
    n = length(x0_ptb_vec) ÷ nsv

    # define base ODE problem for manifold branch
    prob_branch = ODEProblem(rhs!, x0_ptb_vec[1:nsv], tf, (μ))

    # define remake function
    function prob_func(prob, i, repeat)
        remake(prob, u0=x0_ptb_vec[1+(i-1)*nsv : i*nsv])
    end

    # construct Ensemble Problem
    ensemble_prob = EnsembleProblem(prob_branch, prob_func=prob_func)

    # return output of EnsembleProblem
    return solve(ensemble_prob, method, EnsembleThreads(), trajectories=n, callback=callback, method=method, reltol=reltol, abstol=abstol), x0_ptb_vec
end


# ----------------------------------------------------------------------------------- #
# function for extracting poincare section
struct Struct_out_PoincareSection
    u::Matrix{Float64}
    t::Vector{Float64}
end


"""
    get_manifold_ps(sim_manifold)

Get manifold poincare section
"""
function get_manifold_ps(sim_manifold::EnsembleSolution)
    # initialize array of poincare section
    t_ps = zeros(length(sim_manifold))
    x_ps, y_ps   = zeros(length(sim_manifold)), zeros(length(sim_manifold))
    vx_ps, vy_ps = zeros(length(sim_manifold)), zeros(length(sim_manifold))
    if length(sim_manifold[1].u[1])==6
        z_ps, vz_ps = zeros(length(sim_manifold)), zeros(length(sim_manifold))
    end
    # populate poincare section
    for idx in 1:length(sim_manifold)
        t_ps[idx] = sim_manifold[idx].t[end]
        if length(sim_manifold[idx].u[1])==4
            x_ps[idx]  = sim_manifold[idx].u[end][1]
            y_ps[idx]  = sim_manifold[idx].u[end][2]
            vx_ps[idx] = sim_manifold[idx].u[end][3]
            vy_ps[idx] = sim_manifold[idx].u[end][4]

        elseif length(sim_manifold[idx].u[1])==6
            x_ps[idx]  = sim_manifold[idx].u[end][1]
            y_ps[idx]  = sim_manifold[idx].u[end][2]
            z_ps[idx]  = sim_manifold[idx].u[end][3]
            vx_ps[idx] = sim_manifold[idx].u[end][4]
            vy_ps[idx] = sim_manifold[idx].u[end][5]
            vz_ps[idx] = sim_manifold[idx].u[end][6]
        end
    end
    # concatenate into single array
    if length(sim_manifold[1].u[1])==4
        u = cat(x_ps, y_ps, vx_ps, vy_ps, dims=(2,2))
    elseif length(sim_manifold[1].u[1])==6
        u = cat(x_ps, y_ps, z_ps, vx_ps, vy_ps, vz_ps, dims=(2,2))
    end
    return Struct_out_PoincareSection(u, t_ps)
end
