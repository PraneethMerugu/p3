module Manifolds

using ..System
using ..Solvers
using ..Stability
using DifferentialEquations
using LinearAlgebra
using CairoMakie


using Optim

export plot_homoclinic_loop

"""
    plot_homoclinic_loop(output_path::String="homoclinic_loop.png")

Finds and plots the homoclinic orbit connecting the saddle point at the origin to itself.
Uses the shooting method to minimize the return distance.
"""
function plot_homoclinic_loop(output_path::String="homoclinic_loop.png")
    # Saddle point is at (0,0)
    # We need to find ϵ such that the unstable manifold W^u(0) returns to W^s(0).
    # For the unforced Duffing with negative stiffness, the origin is a saddle.
    # The unstable eigenvector direction depends on ϵ.
    
    # Objective function to minimize: Distance to origin after time T
    function objective(ϵ_val)
        p = DuffingParams(ϵ=ϵ_val)
        
        # Calculate unstable eigenvector at (0,0)
        # J = [0 1; 1 -ϵ]
        # λ = (-ϵ + sqrt(ϵ^2 + 4))/2  (Positive eigenvalue)
        λ_u = (-ϵ_val + sqrt(ϵ_val^2 + 4)) / 2
        v_u = [1.0, λ_u]
        v_u = v_u / norm(v_u)
        
        # Initialize small perturbation along unstable manifold
        δ = 1e-4
        u0 = δ * v_u
        
        # Integrate
        # We need to detect when it crosses x=0 again or gets close to origin?
        # Actually, for the loop, it goes out to x ~ 1 and comes back.
        # We integrate for a fixed time or until x crosses 0 from positive to negative?
        
        # Let's integrate for a sufficient time and find the minimum distance to origin
        tspan = (0.0, 50.0)
        prob = ODEProblem(duffing!, u0, tspan, p)
        sol = solve(prob, get_solver(), reltol=1e-6, abstol=1e-6)
        
        # Find min distance after leaving the neighborhood
        # Skip first few points
        dists = [norm(u) for u in sol.u[10:end]]
        return minimum(dists)
    end
    
    # Optimize ϵ
    # We expect ϵ_c < 0 for the homoclinic orbit in this modified Duffing?
    # Or maybe it exists for specific parameters.
    # The PDF mentions "Homoclinic Bifurcation point (Hc)".
    # We scan or optimize.
    
    res = optimize(objective, -0.5, 0.5)
    ϵ_c = Optim.minimizer(res)
    
    # Plotting
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], title = "Homoclinic Loop (ϵ ≈ $(round(ϵ_c, digits=4)))")
    
    # Re-run simulation at optimal ϵ
    p = DuffingParams(ϵ=ϵ_c)
    λ_u = (-ϵ_c + sqrt(ϵ_c^2 + 4)) / 2
    v_u = [1.0, λ_u] / norm([1.0, λ_u])
    u0 = 1e-4 * v_u
    prob = ODEProblem(duffing!, u0, (0.0, 100.0), p)
    sol = solve(prob, get_solver(), reltol=1e-8, abstol=1e-8)
    
    lines!(ax, sol, color=:red, linewidth=2, label="Homoclinic Orbit")
    scatter!(ax, [0.0], [0.0], color=:black, markersize=15, label="Saddle (0,0)")
    
    axislegend(ax)
    save(output_path, fig)
end

end # module Manifolds
