"""
    assign_from_kwargs(keyargs::Dict, keyname, default)

Utility to unpack kwargs
"""
function assign_from_kwargs(keyargs::Dict, keyname, default)
    if haskey(keyargs, keyname)==true
        var = keyargs[keyname]
    else
        var = default
    end
    return var
end
