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
ODE settings class for use with DifferentialEquations
"""
struct ODESettings
    reltol::Float64
    abstol::Float64
    method
end
