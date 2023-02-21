"""
LPO family handling
"""


"""
    lpo2df!(mu::Float64, x0::Vector, period::Float64, m1::Int, m2::Int; kwargs...)

Create dataframe entry from LPO characteristics

# Optional keyword arguments
    - `store_stability::Bool`: whether to compute & store stability
"""
function lpo2df!(mu::Float64, x0::Vector, period::Float64, m1::Int, m2::Int, lp::Int, familyname::String="default"; kwargs...)
    # unpack keyword arguments
    kwargs_dict = Dict(kwargs)
    store_stability = assign_from_kwargs(kwargs_dict, :store_stability, true)
    get_dict = assign_from_kwargs(kwargs_dict, :store_stability, false)

    # construct dictionary
    if length(x0)==6
        dict_entry = Dict(
            "mu" => mu,
            "sv_rx" => x0[1],
            "sv_ry" => x0[2],
            "sv_rz" => x0[3],
            "sv_vx" => x0[4],
            "sv_vy" => x0[5],
            "sv_vz" => x0[6],
            "period" => period,
            "family" => familyname,
            "jacobi" => jacobiConstant(mu, x0),
        )
    elseif length(x0) == 4
        dict_entry = Dict(
            "mu" => mu,
            "sv_rx" => x0[1],
            "sv_ry" => x0[2],
            "sv_vx" => x0[3],
            "sv_vy" => x0[4],
            "period" => period,
            "family" => familyname,
            "jacobi" => jacobiConstant(mu, x0),
        )
    else
        throw(ArgumentError(x0, "State-vector should be length 4 or 6!"))
    end
    if store_stability == true
        # compute stability
        stability = R3BP.stability(mu, x0, period)
        dict_entry["stability"] = stability
    end
    if get_dict == false
        return DataFrame(dict_entry)
    else
        return dict_entry
    end
end


"""
    lpo2df!(df::DataFrame, mu::Float64, x0::Vector, period::Float64, m1::Int, m2::Int, lp::Int, familyname::String="default"; kwargs...)

Mutate DataFrame by appending new member
"""
function lpo2df!(df::DataFrame, mu::Float64, x0::Vector, period::Float64, m1::Int, m2::Int, lp::Int, familyname::String="default"; kwargs...)
    # unpack keyword arguments
    kwargs_dict = Dict(kwargs)
    store_stability = assign_from_kwargs(kwargs_dict, :store_stability, true)
    dict_new = lpo2df!(mu, x0, period, m1, m2, lp, familyname, store_stability=store_stability, get_dict=true)
    # combine to existing dataframe
    push!(df, dict_new)
    return df
end



"""
    find_closest_member(df::DataFrame, column::String, value)

Get index of closest member
"""
function find_closest_member(df::DataFrame, column::String, value)
    diffsqr = (df[!, column] .- value).^2
    idx = argmin(diffsqr)
    return df[idx,:]
end


"""
    get_state_from_df(df::DataFrame, idx::Int)

Get state-vectort from DataFrame
"""
function get_state_from_df(df::DataFrame, idx::Int)
    return [df[idx,:]["sv_rx"], df[idx,:]["sv_ry"], df[idx,:]["sv_rz"],
            df[idx,:]["sv_vx"], df[idx,:]["sv_vy"], df[idx,:]["sv_vz"]]
end


"""
    get_state_from_df(df::DataFrame, idx::Int)

Get state-vectort from DataFrame
"""
function get_state_from_df(df::DataFrameRow)
    return [df["sv_rx"], df["sv_ry"], df["sv_rz"],
            df["sv_vx"], df["sv_vy"], df["sv_vz"]]
end


"""
Construct LPO with closest period
"""
function df_lpo_specify_period(df::Union{DataFrame, DataFrameRow}, period_target::Float64)
    closest = R3BP.find_closest_member(df, "period", period_target)
    return R3BP.ssdc_periodic_xzplane(closest["mu"], R3BP.get_state_from_df(closest),
        period_target, fix="period")
end


"""
    get_family!(p, x0_guess::Vector, period_guess::Float64, dperiod::Float64, maxiter::Int; kwargs...)

Construct family of LPO based on natural parameter continuation in period and single-shooting with xz-symmetry.
"""
function get_family!(
    p, x0_guess::Vector, period_guess::Float64, dperiod::Float64, maxiter::Int; kwargs...
)

    # ---------- extract arguments ---------- #
    kwargs_dict = Dict(kwargs)
    # main manifold options
    max_not_improved = R3BP.assign_from_kwargs(kwargs_dict, :max_not_improved, 5)
    m1 = R3BP.assign_from_kwargs(kwargs_dict, :m1, 0)
    m2 = R3BP.assign_from_kwargs(kwargs_dict, :m2, 0)
    lp = R3BP.assign_from_kwargs(kwargs_dict, :lp, 0)
    fix = R3BP.assign_from_kwargs(kwargs_dict, :fix, "period")
    familyname = R3BP.assign_from_kwargs(kwargs_dict, :familyname, "default")
    period_min = R3BP.assign_from_kwargs(kwargs_dict, :period_min, 0.0)
    period_max = R3BP.assign_from_kwargs(kwargs_dict, :period_max, 100.0)
    method  = assign_from_kwargs(kwargs_dict, :method, Vern7())
    reltol  = assign_from_kwargs(kwargs_dict, :reltol, 1.e-12, Float64)
    abstol  = assign_from_kwargs(kwargs_dict, :abstol, 1.e-12, Float64)
    planar_check = R3BP.assign_from_kwargs(kwargs_dict, :planar_check, true)
    x_mag_max = R3BP.assign_from_kwargs(kwargs_dict, :x_mag_max, 3.0)

    mu = p[1]

    # initialize counters and storage
    i_not_improved = 0
    i_family = 0
    df = DataFrame()

    @showprogress for i = 1:maxiter
        # run new shooting
        res_iter = R3BP.ssdc_periodic_xzplane(p, x0_guess, period_guess;
            fix=fix, method=method, reltol=reltol, abstol=abstol);

        # store & update guess if converged
        if res_iter.flag == 1
            # store into dataframe
            if i_family == 0
                df = R3BP.lpo2df!(mu, res_iter.x0, res_iter.period, m1, m2, lp, familyname)
            else
                df = R3BP.lpo2df!(df, mu, res_iter.x0, res_iter.period, m1, m2, lp, familyname)
            end

            # update initial guess
            x0_guess  = res_iter.x0
            period_guess = res_iter.period + dperiod
            i_family += 1

        # else decrease natural parameter step
        else
            println("Not converged at iteration $i !")
            dperiod = dperiod*0.8
            i_not_improved += 1
        end

        # break if not imporoved too many times
        if (i_not_improved == max_not_improved)
            println("Stopping condition met at $i, period_guess = $period_guess ! (Max number of iteration for single solution)")
            break
        end

        if (period_min > period_guess)
            println("Stopping condition met at $i, period_guess = $period_guess ! (< Min period)")
            break
        end

        if (period_max < period_guess)
            println("Stopping condition met at $i, period_guess = $period_guess ! (> Max period)")
            break
        end

        if (x_mag_max < abs(res_iter.x0[1]))
            println("Stopping condition met at $i, |x| > $x_mag_max ! (> x max magnitude)")
            break
        end

        if planar_check == true
            if (abs(res_iter.x0[3]) < 1.e-8) && (abs(res_iter.x0[6]) < 1.e-8)
                println("Stopping condition met at $i, period_guess = $period_guess ! (Planar hit)")
                break
            end
        end
    end
    return df
end