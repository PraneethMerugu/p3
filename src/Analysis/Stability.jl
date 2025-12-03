module Stability

using ..System
using LinearAlgebra
using StaticArrays
using Parameters


export fixed_points, jacobian, eigenvalues, stability_type

"""
    fixed_points(p::DuffingParams)

Returns the analytical fixed points of the unforced system (A=0).
Points: (0,0), (1,0), (-1,0).
"""
function fixed_points(p::DuffingParams)
    return [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(-1.0, 0.0)]
end

"""
    jacobian(u, p::DuffingParams)

Calculates the Jacobian matrix at state `u`.
J = [0 1; 1-3x^2-2μxy -ϵ-μx^2]
"""
function jacobian(u, p::DuffingParams)
    x, y = u
    @unpack ϵ, μ = p
    
    j11 = 0.0
    j12 = 1.0
    j21 = 1.0 - 3*x^2 - 2*μ*x*y
    j22 = -ϵ - μ*x^2
    
    return SMatrix{2,2}(j11, j21, j12, j22)
end

"""
    eigenvalues(u, p::DuffingParams)

Calculates the eigenvalues of the Jacobian at state `u`.
"""
function eigenvalues(u, p::DuffingParams)
    J = jacobian(u, p)
    return eigen(J).values
end

"""
    stability_type(u, p::DuffingParams)

Classifies the stability of a fixed point.
Returns: :Stable, :Unstable, :Saddle, :Center
"""
function stability_type(u, p::DuffingParams)
    evals = eigenvalues(u, p)
    r1, r2 = real.(evals)
    
    if r1 < 0 && r2 < 0
        return :Stable
    elseif r1 > 0 && r2 > 0
        return :Unstable
    elseif r1 * r2 < 0
        return :Saddle
    else
        return :Center # Pure imaginary or zero real part
    end
end

end # module Stability
