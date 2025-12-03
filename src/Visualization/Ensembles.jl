module Ensembles

using ..System
using ..Solvers
using DifferentialEquations
using CairoMakie
using LinearAlgebra


using Distributions
using Colors

export smale_blob, fractal_basins

"""
    smale_blob(output_path::String="smale_blob.mp4")

Visualizes the stretching and folding of a blob of initial conditions.
Demonstrates positive Lyapunov exponent and mixing.
"""
function smale_blob(output_path::String="smale_blob.mp4")
    # Parameters for Chaos
    p = DuffingParams(A=0.02, ω=1.0, ϵ=0.05) # As per PDF Part II
    
    # Find a point on the attractor (transient removal)
    # Start near origin to avoid escaping to infinity due to negative μ
    prob_transient = ODEProblem(duffing!, [0.1, 0.0], (0.0, 100.0), p)
    sol_transient = solve(prob_transient, Tsit5(), reltol=1e-6, abstol=1e-6)
    attractor_point = sol_transient[end]
    println("Attractor point: ", attractor_point)

    # Initialize Blob around the attractor point
    μ_pos = attractor_point
    σ = 1e-3 # Smaller blob
    d = MvNormal(μ_pos, σ^2 * I)
    
    N_particles = 2000 
    u0s = [rand(d) for _ in 1:N_particles]
    
    # Callback to terminate escaping trajectories
    is_out_of_domain(u, t, integrator) = abs(u[1]) > 10.0 || abs(u[2]) > 10.0
    terminate_cb = DiscreteCallback(is_out_of_domain, terminate!)

    # Ensemble Problem
    prob = ODEProblem(duffing!, [0.0, 0.0], (0.0, 20.0), p) # Shorter duration for demo
    ensemble_prob = EnsembleProblem(prob, prob_func = (prob, i, repeat) -> remake(prob, u0=u0s[i]))
    
    # Solve
    solver = Tsit5()
    sim = solve(ensemble_prob, solver, EnsembleThreads(), trajectories=N_particles, saveat=0.1, reltol=1e-6, abstol=1e-6, callback=terminate_cb)

    
    # Filter successful trajectories or find common length
    valid_sims = [s for s in sim if s.retcode == ReturnCode.Success || s.retcode == ReturnCode.Terminated]
    if isempty(valid_sims)
        println("All trajectories failed!")
        return
    end
    
    # Use the minimum length to avoid bounds errors
    min_len = minimum(length(s) for s in valid_sims)
    times = valid_sims[1].t[1:min_len]
    
    fig = Figure(size = (800, 800))
    ax = Axis(fig[1, 1], title = "Smale Blob Evolution", limits = (-2, 2, -2, 2))
    
    points_obs = Observable(Point2f[])
    scatter!(ax, points_obs, color=:green, markersize=2, transparency=true)
    
    record(fig, output_path, 1:min_len; framerate=30) do i
        current_points = Point2f[]
        for sol in valid_sims
            if i <= length(sol)
                push!(current_points, Point2f(sol[i][1], sol[i][2]))
            end
        end
        points_obs[] = current_points
    end
end

"""
    fractal_basins(output_path::String="fractal_basins.png")

Generates a high-resolution map of the basins of attraction.
Grid: 1000x1000 (Reduced to 200x200 for demo speed).
"""
function fractal_basins(output_path::String="fractal_basins.png")
    p = DuffingParams(ϵ=0.1, A=0.0) # Autonomous multistable regime
    
    # Grid
    res = 200
    x_range = range(-2, 2, length=res)
    y_range = range(-2, 2, length=res)
    
    basins = zeros(Int, res, res)
    
    # We can use DynamicalSystems.basins_of_attraction if available, 
    # but let's do a simple endpoint check.
    
    # Attractors:
    # 1: Right Well (x > 0)
    # 2: Left Well (x < 0)
    # 3: Origin (Unlikely for stable)
    
    # Parallel loop?
    # Simple serial for now to avoid complexity with shared arrays in this script
    
    tspan = (0.0, 100.0)
    prob = ODEProblem(duffing!, [0.0, 0.0], tspan, p)
    
    for (i, x) in enumerate(x_range)
        for (j, y) in enumerate(y_range)
            _prob = remake(prob, u0=[x, y])
            sol = solve(_prob, get_solver(), reltol=1e-4, abstol=1e-4, save_everystep=false)
            final_x = sol[end][1]
            
            if final_x > 0.5
                basins[i, j] = 1 # Right
            elseif final_x < -0.5
                basins[i, j] = 2 # Left
            else
                basins[i, j] = 0 # Boundary/Origin
            end
        end
    end
    
    # Plot
    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1], title = "Fractal Basins of Attraction")
    heatmap!(ax, x_range, y_range, basins, colormap = [:black, :red, :blue])
    
    save(output_path, fig)
end

end # module Ensembles
