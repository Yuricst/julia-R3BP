"""
Stability functions
"""


"""
    stability(λ_unstb::Float64, λ_stb::Float64)

Compute linear stability given stable and unstable eigenvalues
"""
function stability(λ_unstb::Float64, λ_stb::Float64)
    return 0.5*(abs(λ_unstb) + abs(λ_stb))
end


"""
    stability(λ_unstb::Float64, λ_stb::Float64)

Compute linear stability given stable and unstable eigenvalues
"""
function stability(monodromy)
    λs = eigvals(monodromy)
    vs = eigvecs(monodromy)
    λ_unstb, λ_stb, _, _ = get_stable_unstable_eigvecs(λs, vs)
    return stability(λ_unstb, λ_stb)
end


"""
    stability(x::Vector, period::Float64)

Compute linear stability given stable and unstable eigenvalues
"""
function stability(μ::Float64, x0::Vector, period::Float64; kwargs...)
    # unpack arguments
    kwargs_dict = Dict(kwargs)
    reltol = assign_from_kwargs(kwargs_dict, :reltol, 1.e-12)
    abstol = assign_from_kwargs(kwargs_dict, :abstol, 1.e-12)
    method = assign_from_kwargs(kwargs_dict, :method, Tsit5())

    # integrate x0 by full period
    if length(x0)==4
        #x0_stm = vcat(x0, [1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1]);
        x0_stm = vcat(x0iter, reshape(I(4), (16,)))[:]
        prob_lpo = ODEProblem(rhs_pcr3bp_svstm!, x0_stm, period, (μ),
            method=method, reltol=reltol, abstol=abstol
        );
    elseif length(x0)==6
        #x0_stm = vcat(x0, [1 0 0 0 0 0  0 1 0 0 0 0  0 0 1 0 0 0  0 0 0 1 0 0  0 0 0 0 1 0  0 0 0 0 0 1]);
        x0_stm = vcat(x0, reshape(I(6), (36,)))[:]
        prob_lpo = ODEProblem(rhs_cr3bp_svstm!, x0_stm, period, (μ),)
            #method=method, reltol=reltol, abstol=abstol);
    else
        error("x0 should be length 4 or 6")
    end
    sol = solve(prob_lpo, method, reltol=reltol, abstol=abstol)
    # get monodromy matrix (careful of order from reshape function!)
    monodromy = reshape(sol.u[end][length(x0)+1:end], (length(x0),length(x0)))';
    # call stability function with dispatch for monodromy
    return stability(monodromy)
end
