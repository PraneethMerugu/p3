module PhaseSpace

using ..System
using ..Solvers
using ..Stability
using CairoMakie
using ColorSchemes



using DifferentialEquations

export hopf_scanner

"""
    hopf_scanner(output_path::String="hopf_scanner.mp4")

Generates an animation of the phase portrait as ϵ varies from 0.5 to -0.1.
Visualizes the Supercritical Hopf Bifurcation at ϵ = 0.225.
"""
function hopf_scanner(output_path::String="hopf_scanner.mp4")
    # Parameter range
    ϵ_range = range(0.5, -0.1, length=300)
    
    # Observable for current ϵ
    ϵ_obs = Observable(ϵ_range[1])
    
    # Setup Figure
    fig = Figure(resolution = (800, 800))
    ax = Axis(fig[1, 1], 
        title = @lift("Hopf Scanner: ϵ = $(round($ϵ_obs, digits=3))"),
        xlabel = "x", ylabel = "y",
        limits = (-2, 2, -2, 2)
    )
    
    # Vector Field Grid
    xs = range(-2, 2, length=20)
    ys = range(-2, 2, length=20)
    points_flat = [Point2f(x, y) for x in xs, y in ys] |> vec
    
    # Observables for vector field (Vector of Point2f)
    vectors_obs = Observable(Vector{Point2f}(undef, length(points_flat)))
    
    # Initial calculation
    function get_field(ϵ)
        p = DuffingParams(ϵ=ϵ)
        vecs = Vector{Point2f}(undef, length(points_flat))
        for (k, pt) in enumerate(points_flat)
            x, y = pt[1], pt[2]
            dx = y
            dy = x - x^3 - ϵ*y - p.μ*x^2*y
            # Normalize for visualization
            mag = sqrt(dx^2 + dy^2)
            if mag > 0
                vecs[k] = Point2f(dx / mag, dy / mag)
            else
                vecs[k] = Point2f(0.0, 0.0)
            end
        end
        return vecs
    end
    
    initial_vectors = get_field(ϵ_range[1])
    
    # Plot Arrows (Static initial call)
    arrow_plt = arrows!(ax, points_flat, initial_vectors, arrowsize=10, lengthscale=0.1, color=:gray)
    
    # Streamlines - seed points in a grid
    streamline_seeds = [
        Point2f(-1.5, -1.5), Point2f(-1.5, 0.0), Point2f(-1.5, 1.5),
        Point2f(0.0, -1.5), Point2f(0.0, 1.5),
        Point2f(1.5, -1.5), Point2f(1.5, 0.0), Point2f(1.5, 1.5)
    ]
    
    # Function to compute streamlines with velocity
    function compute_streamlines(ϵ, seeds)
        p = DuffingParams(ϵ=ϵ)
        streamlines = Vector{Vector{Point2f}}()
        velocities = Vector{Vector{Float64}}()
        
        for seed in seeds
            prob = ODEProblem(duffing!, [seed[1], seed[2]], (0.0, 10.0), p)
            sol = solve(prob, Tsit5(), reltol=1e-4, abstol=1e-4, saveat=0.1)
            if !isempty(sol)
                points = Point2f[]
                vels = Float64[]
                for u in sol.u
                    if abs(u[1]) < 2.5 && abs(u[2]) < 2.5
                        push!(points, Point2f(u[1], u[2]))
                        # Compute velocity magnitude
                        dx = u[2]
                        dy = u[1] - u[1]^3 - ϵ*u[2] - p.μ*u[1]^2*u[2]
                        vel = sqrt(dx^2 + dy^2)
                        push!(vels, vel)
                    end
                end
                if length(points) > 1
                    push!(streamlines, points)
                    push!(velocities, vels)
                end
            end
        end
        return streamlines, velocities
    end
    
    # Initial streamlines
    initial_streamlines, initial_velocities = compute_streamlines(ϵ_range[1], streamline_seeds)
    
    # Create colored streamlines
    streamline_plots = []
    for (sl, vels) in zip(initial_streamlines, initial_velocities)
        # Normalize velocities for color mapping
        vmin, vmax = 0.0, 2.0  # Adjust range as needed
        colors = [get(ColorSchemes.viridis, clamp((v - vmin) / (vmax - vmin), 0, 1)) for v in vels]
        plt = lines!(ax, sl, color=colors, linewidth=2.5)
        push!(streamline_plots, plt)
    end
    
    # Fixed Points
    # We want to color them based on stability
    # This is tricky in Makie with Observables for color, so we'll just plot markers
    # and update their color in the loop if possible, or just static for now.
    # Better: Re-plot markers on every frame? Efficient enough for 3 points.
    
    fp_scatter = scatter!(ax, [0.0, 1.0, -1.0], [0.0, 0.0, 0.0], 
        color = [:red, :green, :green], markersize = 20, strokewidth=2)

    # Trajectory Tracer
    # We simulate a short trajectory for the current ϵ
    traj_points = Observable(Point2f[])
    lines!(ax, traj_points, color=:cyan, linewidth=2)

    # Animation Loop
    record(fig, output_path, ϵ_range; framerate = 30) do ϵ_val
        ϵ_obs[] = ϵ_val
        
        # Update Vector Field manually
        new_vecs = get_field(ϵ_val)
        arrow_plt[2] = new_vecs # Update the second argument (vectors)
        
        # Update Streamlines
        new_streamlines, new_velocities = compute_streamlines(ϵ_val, streamline_seeds)
        for (i, (sl, vels)) in enumerate(zip(new_streamlines, new_velocities))
            if i <= length(streamline_plots)
                streamline_plots[i][1] = sl
                # Update colors based on velocity
                vmin, vmax = 0.0, 2.0
                colors = [get(ColorSchemes.viridis, clamp((v - vmin) / (vmax - vmin), 0, 1)) for v in vels]
                streamline_plots[i].color = colors
            end
        end
        
        # Update Fixed Point Colors
        p = DuffingParams(ϵ=ϵ_val)
        stabilities = [stability_type(fp, p) for fp in [[0.0,0.0], [1.0,0.0], [-1.0,0.0]]]
        colors = [s == :Stable ? :green : :red for s in stabilities]
        fp_scatter.color = colors
        
        # Update Trajectory
        # Start near x=1
        u0 = [1.01, 0.01]
        tspan = (0.0, 50.0)
        prob = ODEProblem(duffing!, u0, tspan, p)
        sol = solve(prob, get_solver(), reltol=1e-6, abstol=1e-6, saveat=0.1)
        traj_points[] = Point2f.(sol[1,:], sol[2,:])
    end
end

end # module PhaseSpace
