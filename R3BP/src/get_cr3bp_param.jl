"""
Function define CR3BP system parameters
"""

struct CR3BP_param
    mu::Float64
    lstar::Float64
    tstar::Float64
    mstar::Float64
    m2_soi::Float64
end


"""
    get_cr3bp_param(m1_naifID::Int, m2_naifID::Int)

Obtain CR3BP parameters mu, Lstar, Tstar, soi of m2

Args:
    m1_naifID (Int): mass of first primary
    m2_naifID (Int): mass of second primary

Returns:
    (obj): object with fields:
        mu (float): mass-parameter
        lstar (float): non-dimensional distance
        tstar (float): non-dimensional time
        m2_soi (float): sphere of influence of second mass, in km
"""
function get_cr3bp_param(m1_naifID::Int, m2_naifID::Int)
    return get_cr3bp_param(string(m1_naifID), string(m2_naifID))
end



"""
    get_cr3bp_param(m1_naifID::String, m2_naifID::String)

Obtain CR3BP parameters mu, Lstar, Tstar, soi of m2

Args:
    m1_naifID (str): mass of first primary
    m2_naifID (str): mass of second primary

Returns:
    (obj): object with fields:
        mu (float): mass-parameter
        lstar (float): non-dimensional distance
        tstar (float): non-dimensional time
        m2_soi (float): sphere of influence of second mass, in km
"""
function get_cr3bp_param(m1_naifID::String, m2_naifID::String)
    # list of semi-major axis
    if m1_naifID =="10"
        a2 = get_semiMajorAxes(m2_naifID)[1]
    elseif m2_naifID =="10"
        a2 = get_semiMajorAxes(m1_naifID)[1]
    else
        semiMajorAxes = get_semiMajorAxes(m1_naifID, m2_naifID)
        a1, a2 = semiMajorAxes[1], semiMajorAxes[2]
    end
    # list of gm
    gmlst = get_gm(m1_naifID, m2_naifID)
    m1_gm, m2_gm = gmlst[1], gmlst[2]
    if m1_gm < m2_gm
        error("Excepting m1 > m2!")
    end
    # create list of parameters
    mu     = m2_gm / (m1_gm + m2_gm)
    lstar  = a2
    tstar  = sqrt( ( a2 )^3 / ( m1_gm + m2_gm ) )
    mstar  = m1_gm + m2_gm
    m2_soi = a2 * (m2_gm/m1_gm)^(2/5)
    return CR3BP_param(mu, lstar, tstar, mstar, m2_soi)
end
