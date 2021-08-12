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
            "m1" => m1,
            "m2" => m2,
            "family" => familyname
        )
    elseif length(x0) == 4
        dict_entry = Dict(
            "mu" => mu,
            "sv_rx" => x0[1],
            "sv_ry" => x0[2],
            "sv_vx" => x0[3],
            "sv_vy" => x0[4],
            "period" => period,
            "m1" => m1,
            "m2" => m2,
            "family" => familyname
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


# function lpo2df(df::DataFrame, mu::Float64, x0::Vector, period::Float64, m1::Int, m2::Int, lp::Int, familyname::String="default"; kwargs...)
#     # unpack keyword arguments
#     kwargs_dict = Dict(kwargs)
#     store_stability = assign_from_kwargs(kwargs_dict, :store_stability, true)
#     dict_new = lpo2df(mu, x0, period, m1, m2, lp, familyname, store_stability=store_stability, get_dict=true)
#     # combine to existing dataframe
#     push!(df, dict_new)
#     return df
# end


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
