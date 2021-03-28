# manifold_sampling
Sampling from manifold in CR3BP

### Dependencies
Dependencies used in this repo are: `DifferentialEquations`, `LinearAlgebra`, `Distributed`, `Printf`, `Plots`, `DataFrames`, `CSV`, `CSVFiles`, `Statistics`. 


### Integration of PCR3BP / CR3BP states
For integrating a trajectory, use the DifferentialEquations module together with R3BP: 

```julia
using DifferentialEquations
using Plots
include("R3BP.jl")

# define parameters
mu = 0.01215058426994;
X0 = [1.176924090973164 0.0 -0.060210863312217 0.0 -0.173836346247689 0.0];
T = 3.385326412831325;

tspan = (0.0, T)
p = (mu)
prob = ODEProblem(R3BP.rhs_cr3bp_sv!, X0, tspan, p)
sol = DifferentialEquations.solve(prob, Tsit5(), reltol=1e-11, abstol=1e-11)

# plot using the sol object
plot(sol, vars=(1,2,3))
```

### Manifold
A manifold of a LPO is usually propagated until a Poincare section, which is to be defined by a callback function:
```julia
function condition(u,t,integrator)
  u[1] - 1.25   # when y-value hits xx
end

affect!(integrator) = terminate!(integrator)

# assign callback
cb = ContinuousCallback(condition,affect!)
```
then define parameters for generating manifold
```julia
# parameters for manifolds
num_branch = 100;
stability = true;
epsilon = 1e-5
lstar = 384400.
relative_tol_manifold = 0.1
absolute_tol_manifold_km = 100.0
tf = -7.0
```
and finally the function wraps DifferentialEquations' `EnsembleProblem()`
```julia
# generate manifolds
outsim = R3BP.get_manifold(mu, X0, T, tf, num_branch, stability, epsilon, cb, 
    "positive", lstar, relative_tol_manifold, absolute_tol_manifold_km, 1e-11, 1e-11, Tsit5())
```


