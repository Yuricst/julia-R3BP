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

Obtain CR3BP parameters mu, lstar, tstar, soi of m2

# Arguments
    m1_naifID (str): mass of first primary
    m2_naifID (str): mass of second primary

# Returns
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



struct BCR4BP_param
    mu::Float64
    lstar::Float64
    tstar::Float64
    mstar::Float64
    m2_soi::Float64
    μ_3::Float64
    a::Float64
    ω_s::Float64
    tsyn::Float64
end



"""
    get_bcr4bp_param(m1_naifID::Int, m2_naifID::Int)

Obtain BCR4BP parameters mu, Lstar, Tstar, soi of m2, from m1-m2 planet-moon-sun system

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
function get_bcr4bp_param(m1_naifID::Int, m2_naifID::Int)
    return get_bcr4bp_param(string(m1_naifID), string(m2_naifID))
end


"""
    get_bcr4BP_param(m1_naifID::String, m2_naifID::String)

Obtain BCR4BP parameters mu, lstar, tstar, soi of m2, from m1-m2 planet-moon-sun system

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
function get_bcr4bp_param(m1_naifID::String, m2_naifID::String)
    CR3BP_param = get_cr3bp_param(m1_naifID, m2_naifID)

    # get information of Sun
    gm_sun     = get_gm("10")[1]
    a_m1_sun   = get_semiMajorAxes(m1_naifID)[1]   # semi-major axis of first (larger) body
    println("a_m1_sun: $a_m1_sun, gm_sun: $gm_sun")
    period_sun = 2π*sqrt(a_m1_sun^3/gm_sun)   # [sec]

    # compute rotation rate of sun about m1-m2 barycenter
    tsyn = get_synodic_period(period_sun, 2π*CR3BP_param.tstar)  # [sec]
    ω_s  = -2π / (tsyn / CR3BP_param.tstar)                    # [rad/canonical time]
    println("ω_s: $ω_s")
    # compute scaled gm and length of sun
    a_sun = a_m1_sun / CR3BP_param.lstar
    μ_3   = gm_sun   / CR3BP_param.mstar

    # return structure
    return BCR4BP_param(CR3BP_param.mu, CR3BP_param.lstar, CR3BP_param.tstar,
    CR3BP_param.mstar, CR3BP_param.m2_soi, μ_3, a_sun, ω_s, tsyn)
end



"""

Compute synodic period between two systems with periods p1 and p2

# Arguments
    p1 (float): period of first system
    p2 (float): period of second system

# Returns
    (float): synodic period
"""
function get_synodic_period(p1, p2)
    return 1/abs( 1/p1 - 1/p2 )
end
