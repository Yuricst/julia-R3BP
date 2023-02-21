"""
    assign_from_kwargs(keyargs::Dict, keyname, default_value, type=nothing)

Utility to unpack kwargs
"""
function assign_from_kwargs(keyargs::Dict, keyname, default_value, type=nothing)
    if haskey(keyargs, keyname)==true
        var = keyargs[keyname]
    else
        var = default_value
    end
    # check type if provided
    if isnothing(type) == false
        given_type = typeof(var)
        if given_type != type
            error("kwarg $keyname got $given_type, should be $type")
        end
    end
    return var
end


"""
    __print_verbose(message::String, verbose::Bool)

Print message if verbose==true
"""
function __print_verbose(message::String, verbose::Bool)
    if verbose==true
        println(message)
    end
end


"""
    __print_verbosity(message::String, verbosity::Int, threshold::Int)

Print message if verbosity > threshold
"""
function __print_verbosity(message::String, verbosity::Int, threshold::Int)
    if verbosity > threshold
        println(message)
    end
end



"""
ODE settings class for use with DifferentialEquations
"""
struct ODESettings
    reltol::Float64
    abstol::Float64
    method
end


"""
Compute M-N resonance period
"""
function resonance_period(M,N)
    return 2*N*π/M
end


"""
Load JSON file to dataframe
"""
function json_to_df(path::String)
    dict = Dict()
    open(path, "r") do f
        dict = JSON.parse(f)  # parse and transform data
    end
    return DataFrame(dict)
end


"""
Save DataFrame as JSON
"""
function df_to_json(path::String, df::Union{DataFrame, DataFrameRow})
    dict = Dict(pairs(eachcol(dfcomb)))
    open(path,"w") do f
        JSON.print(f, dict)
    end
end
