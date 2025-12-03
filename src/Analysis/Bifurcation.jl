module Bifurcation

using ..System
using BifurcationKit
using Accessors
using StaticArrays
using Parameters



export get_bifurcation_problem, continuation_opts

"""
    get_bifurcation_problem(p::DuffingParams)

Creates a BifurcationProblem for the autonomous system.
Parameter to vary: ϵ (linear damping).
"""
function get_bifurcation_problem(p::DuffingParams)
    # Define the vector field for BifurcationKit
    # F(u, p) where p is the parameter vector (or struct)
    function F_duffing(u, par)
        x, y = u
        @unpack ϵ, μ, A, ω = par
        dx = y
        dy = x - x^3 - ϵ*y - μ*x^2*y # Autonomous part (A=0 assumed for bifurcation of equilibria)
        return [dx, dy]
    end

    # Initial guess (fixed point at x=1)
    u0 = [1.0, 0.0]

    
    # Parameter axis: Lens for ϵ
    param_axis = (@optic _.ϵ)
    
    return BifurcationProblem(F_duffing, u0, p, param_axis)

end

"""
    continuation_opts()

Returns standard continuation options for PALC.
"""
function continuation_opts()
    return ContinuationPar(
        p_min = -0.5,
        p_max = 0.5,
        ds = 0.001,        # Initial step size
        dsmax = 0.01,      # Max step size
        dsmin = 1e-4,      # Min step size
        max_steps = 10000,
        detect_bifurcation = 3, # Detect all bifurcations
        n_inversion = 6,
        max_bisection_steps = 25,
        nev = 2            # Number of eigenvalues to compute
    )
end

end # module Bifurcation
