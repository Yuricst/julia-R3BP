"""
Function associated with manifold
"""


# -------------------------------------------------------------------------------- #
# methods associated with STM
"""
    get_stm(sol, idx::Int)

Extract STM from solution of ODEProblem
"""
function get_stm(sol, idx::Int)
    if length(sol.u[1])==20
        stm = reshape(sol.u[idx][4+1:end], (4, 4))'
    elseif length(sol.u[1])==42
        stm = reshape(sol.u[idx][6+1:end], (6, 6))'
    end
    return stm
end


# -------------------------------------------------------------------------------- #
# methods associated with manifold computation
"""
    get_stable_unstable_eigvecs(λs, vs)

Identidy stable and unstable eigenvalue-eigenvector pair
"""
function get_stable_unstable_eigvecs(λs, vs)
    idx_stb_unstb = [];
    for idx in 1:length( λs )
        if abs(real(λs[idx]) - 1) > 1e-4 && imag(λs[idx])==0
            push!(idx_stb_unstb, idx);
        end
    end

    if abs(λs[ idx_stb_unstb[1] ]) > abs(λs[ idx_stb_unstb[2] ])
        eig_unstb = real( λs[ idx_stb_unstb[1] ] );
        v_unstb   = real( vs[:, idx_stb_unstb[1] ] );
        eig_stb   = real( λs[ idx_stb_unstb[2] ] );
        v_stb     = real( vs[:, idx_stb_unstb[2] ] );
    else
        eig_unstb = real( λs[ idx_stb_unstb[2] ] );
        v_unstb   = real( vs[:, idx_stb_unstb[2] ] );
        eig_stb   = real( λs[ idx_stb_unstb[1] ] );
        v_stb     = real( vs[:, idx_stb_unstb[1] ] );
    end
    return eig_unstb, eig_stb, v_unstb, v_stb
end


"""
    get_eigenvector(monodromy, stable::Bool=true)

Get eigenvectors from monodromy matrix
"""
function get_eigenvector(monodromy, stable::Bool=true)
    λs = eigvals(monodromy)
    vs = eigvecs(monodromy)
    eig_unstb, eig_stb, v_unstb, v_stb = get_stable_unstable_eigvecs(λs, vs)
    @printf("Linear stability ν = %1.4f \n", 0.5*(abs(eig_unstb) + abs(eig_stb)))
    @printf("Eigenvalue unstable: %1.4f, stable: %1.4f, product: %1.4f\n", eig_unstb, eig_stb, eig_unstb*eig_stb)
    if stable==true
        return v_stb
    else
        return v_unstb
    end
end


"""
    scale_ϵ(μ, x0, period, stable, monodromy, y0, lstar::Float64, relative_tol_manifold::Float64=0.1, absolute_tol_manifold_km::Float64=100.0)

Obtain linear perturbation ϵ magnitude for manifolds

# Arguments
    - `μ::Float64`: CR3BP parameter
"""
function scale_ϵ(μ, x0, period, stable, monodromy, y0, lstar::Float64, relative_tol_manifold::Float64=0.1, absolute_tol_manifold_km::Float64=100.0)
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
            prob_estErr = ODEProblem(rhs_pcr3bp_sv!, x0_manifold, (0, -period), (μ));
        else               # unstable case
            #prop_manif_out = propagate_cr3bp(μ, x0_manifold, period)
            prob_estErr = ODEProblem(rhs_pcr3bp_sv!, x0_manifold, (0, period), (μ));
        end
        sol = solve(prob_estErr, Tsit5(), reltol=1e-11, abstol=1e-11);

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


## manifold function
"""
    get_manifold(
           μ::Float64,
           x0::Array,  # ::Array{Float64,2},
           period::Float64,
           tf::Float64,
           stable::Boolean;
           kwargs...)

Function to obtain manifold of LPO

# Optional arugments:
    - ϵ
    - xdir (default: "positive")
    - n
    - callback
    - lstar
    - relative_tol_manifold
    - absolute_tol_manifold_km
    - reltol
    - abstol
    - method
    """
function get_manifold(
            μ::Float64,
            x0::Array,
            period::Float64,
            tf::Float64,
            stable::Bool;
            kwargs...)
    # ---------- extract arguments ---------- #
    if :xdir in keys(kwargs)
        xdir = kwargs[:xdir];
    else
        xdir = "positive";
    end

    if :n in keys(kwargs)
        n = kwargs[:n];
    else
        n = 50
    end

    if :callback in keys(kwargs)
        callback = kwargs[:callback];
    else
        callback = nothing;
    end

    if :lstar in keys(kwargs)
        lstar = kwargs[:lstar];
    else
        lstar=384400.0;
    end

    if :relative_tol_manifold in keys(kwargs)
        relative_tol_manifold = kwargs[:relative_tol_manifold];
    else
        relative_tol_manifold=0.1;
    end

    if :absolute_tol_manifold_km in keys(kwargs)
        absolute_tol_manifold_km = kwargs[:absolute_tol_manifold_km];
    else
        absolute_tol_manifold_km=100.0;
    end

    if :reltol in keys(kwargs)
        reltol = kwargs[:reltol];
    else
        reltol=1e-11;
    end

    if :abstol in keys(kwargs)
        abstol = kwargs[:abstol];
    else
        abstol=1e-11;
    end

    if :method in keys(kwargs)
        method = kwargs[:method];
    else
        method=Tsit5();
    end

    if :verbosity in keys(kwargs)
        verbosity = kwargs[:verbosity]
    else
        verbosity = 0
    end

    if verbosity > 0
        println("========== Manifold Setting ==========")
        print("xdir: ")
        println(xdir)
        print("Stable: ")
        println(stable)
    end

    # ---------- propagate c0 by one full period with STM ---------- #
    if length(x0)==4
        #x0_stm = vcat(x0, [1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1]);
        x0_stm = vcat(x0iter, reshape(I(4), (16,)))[:]
        prob_lpo = ODEProblem(rhs_pcr3bp_svstm!, x0_stm, period, (μ));
    elseif length(x0)==6
        #x0_stm = vcat(x0, [1 0 0 0 0 0  0 1 0 0 0 0  0 0 1 0 0 0  0 0 0 1 0 0  0 0 0 0 1 0  0 0 0 0 0 1]);
        x0_stm = vcat(x0, reshape(I(6), (36,)))[:]
        prob_lpo = ODEProblem(rhs_cr3bp_svstm!, x0_stm, period, (μ));
    else
        error("x0 should be length 4 or 6")
    end
    ts_lpo = LinRange(0, period, n+1)
    sol = solve(prob_lpo, method, reltol=reltol, abstol=abstol, saveat=ts_lpo);

    # get monodromy matrix (careful of order from reshape function!)
    monodromy = reshape(sol.u[end][length(x0)+1:end], (length(x0),length(x0)))';

    # get eigenvectors at initial state
    y0 = get_eigenvector(monodromy, stable);

    # iterate over each branch to be generated
    #ts_branch = LinRange(0, tf, steps);

    # define ϵ (linear perturbation)
    if :ϵ in keys(kwargs)
        ϵ = kwargs[:ϵ];
    else
        ϵ = scale_ϵ(μ, x0, period, stable, monodromy, y0, lstar, relative_tol_manifold, absolute_tol_manifold_km);
    end

    # decide sign of ϵ based on xdir
    if cmp(xdir,"positive")==0
        if y0[1] < 0.0
            ϵ_corr = -abs(ϵ)
        else
            ϵ_corr =  abs(ϵ)
        end
    elseif cmp(xdir,"negative")==0
        if y0[1] > 0.0
            ϵ_corr = -abs(ϵ)
        else
            ϵ_corr =  abs(ϵ)
        end
    else
        error("xdir should be \"positive\" or \"negative\"")
    end
    if verbosity > 0
        @printf("Using linear perturbation ϵ = %s \n", ϵ_corr)
    end

    # ---------- construct and append initial condition ---------- #
    x0_ptrbs = []
    for idx_x0 in 1:n
        # map eigenvector (y = stm*y0)
        y_transcribed = reshape(sol.u[idx_x0][length(x0)+1:end], (length(x0), length(x0)))' * reshape(y0, (length(x0), 1));

        # construct perturbed state
        x0_ptrb   = sol.u[idx_x0][1:length(x0)] + ϵ_corr*y_transcribed/norm(y_transcribed);
        push!(x0_ptrbs, x0_ptrb)
    end

    # define base ODE problem for manifold branch
    if length(x0)==4
        prob_branch = ODEProblem(rhs_pcr3bp_sv!, reshape(x0_ptrbs[1],(1,4)), tf, (μ));
    elseif length(x0)==6
        prob_branch = ODEProblem(rhs_cr3bp_sv!, reshape(x0_ptrbs[1],(1,6)), tf, (μ));
    end

    # ---------- ensemble siμlation ---------- #
    function prob_func(prob, i, repeat)
        remake(prob, u0=reshape(x0_ptrbs[i],(1,length(x0))))
    end

    ensemble_prob = EnsembleProblem(prob_branch, prob_func=prob_func)
    if isnothing(callback)
        #println("No callback function")
        sim = solve(ensemble_prob, method, EnsembleThreads(), trajectories=n, reltol=reltol, abstol=abstol)
    else
        #println("Using callback function")
        sim = solve(ensemble_prob, method, EnsembleThreads(), trajectories=n, callback=callback, reltol=reltol, abstol=abstol)
    end
    return sim
end


# ----------------------------------------------------------------------------------- #
# function for extracting poincare section
# solution output
struct Struct_out_PoincareSection
    u
    t
end


"""
    get_manifold_ps(sim_manifold)

Get manifold poincare section
"""
function get_manifold_ps(sim_manifold)
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
