"""
Typical callback functions
"""


"""
    check_cr3bp_apsis(sv, mu::Float64, m::Int=2)

Check if in w.r.t. primary m=1 or m=2 by returning dot(r,v)
"""
function check_cr3bp_apsis(sv, mu::Float64, m::Int=2)
    # relocate x-axis value
    if m == 2
        x = sv[1] - (1-mu)
    else
        x = sv[1] + mu
    end
    # re-construct r and v vectors
    if length(sv) == 4
        r = [x, sv[2]]
        v = sv[3:4]
    elseif length(sv) == 6
        r = [x, sv[2], sv[3]]
        v = sv[4:6]
    end
    return dot(r,v)
end
