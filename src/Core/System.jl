module System

using Parameters
using StaticArrays

export DuffingParams, duffing!

"""
    DuffingParams

Parameters for the Modified Duffing Oscillator.
- `ϵ`: Linear damping coefficient.
- `μ`: Nonlinear damping coefficient (Fixed at -0.225).
- `A`: Forcing amplitude.
- `ω`: Forcing frequency.
"""
@with_kw struct DuffingParams{T<:Real}
    ϵ::T = 0.1
    μ::T = -0.225
    A::T = 0.0
    ω::T = 1.0
end

"""
    duffing!(du, u, p, t)

In-place definition of the Modified Duffing Oscillator ODEs.
Equations:
    dx/dt = y
    dy/dt = x - x^3 - ϵ*y - μ*x^2*y + A*cos(ω*t)
"""
function duffing!(du, u, p::DuffingParams, t)
    x, y = u
    @unpack ϵ, μ, A, ω = p
    
    du[1] = y
    du[2] = x - x^3 - ϵ*y - μ*x^2*y + A*cos(ω*t)
    return nothing
end

end # module System
