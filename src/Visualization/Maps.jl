module Maps

using ..System
using ..Solvers
using ..Chaos
using DynamicalSystemsBase
using CairoMakie
using DifferentialEquations




using Statistics
using LinearAlgebra

export plot_poincare, plot_cobweb, plot_complexity_spectrum, plot_markov_heatmap

"""
    plot_poincare(output_path::String="poincare_map.png")

Generates Poincaré sections for Periodic (ω=0.9) and Chaotic (ω=1.0) regimes.
"""
function plot_poincare(output_path::String="poincare_map.png")
    # Periodic Regime
    p_periodic = DuffingParams(A=0.02, ω=0.9, ϵ=0.05)
    ds_periodic = CoupledODEs(duffing!, [0.1, 0.1], p_periodic)
    
    # Chaotic Regime
    p_chaotic = DuffingParams(A=0.02, ω=1.0, ϵ=0.05)
    ds_chaotic = CoupledODEs(duffing!, [0.1, 0.1], p_chaotic)
    
    # Compute Poincaré Maps
    # Plane: t % (2π/ω) = 0
    # DynamicalSystems.poincaremap uses a plane definition.
    # For stroboscopic map of non-autonomous system, we can use trajectory with sampling.
    
    # Callback to terminate escaping trajectories
    is_out_of_domain(u, t, integrator) = abs(u[1]) > 10.0 || abs(u[2]) > 10.0
    terminate_cb = DiscreteCallback(is_out_of_domain, terminate!)

    # Ensemble parameters
    N_traj = 100
    T_max_periods = 50
    
    function get_poincare_points(p, N_traj, T_max_periods)
        T = 2π / p.ω
        points = Point2f[]
        
        # Grid of ICs
        x_range = range(-0.5, 0.5, length=isqrt(N_traj))
        y_range = range(-0.5, 0.5, length=isqrt(N_traj))
        
        for x0 in x_range, y0 in y_range
            prob = ODEProblem(duffing!, [x0, y0], (0.0, T_max_periods * T), p)
            sol = solve(prob, Tsit5(), callback=terminate_cb, saveat=T)
            
            for u in sol.u
                if abs(u[1]) < 2.0 && abs(u[2]) < 2.0 # Only plot if within bounds
                    push!(points, Point2f(u[1], u[2]))
                end
            end
        end
        return points
    end

    # Periodic
    println("Generating Periodic Ensemble...")
    points_p = get_poincare_points(p_periodic, 200, 50)
    println("Periodic points: ", length(points_p))
    
    # Chaotic
    println("Generating Chaotic Ensemble...")
    points_c = get_poincare_points(p_chaotic, 200, 50)
    println("Chaotic points: ", length(points_c))
    
    # Plot
    fig = Figure(size = (1000, 500))
    
    ax1 = Axis(fig[1, 1], title = "Poincaré Map: Periodic (ω=0.9)", xlabel="x", ylabel="y", limits=(-2, 2, -2, 2))
    if !isempty(points_p)
        scatter!(ax1, points_p, markersize=2, color=:blue)
    end
    
    ax2 = Axis(fig[1, 2], title = "Poincaré Map: Chaotic (ω=1.0)", xlabel="x", ylabel="y", limits=(-2, 2, -2, 2))
    if !isempty(points_c)
        scatter!(ax2, points_c, markersize=1, color=:red, transparency=true)
    end
    
    save(output_path, fig)
end

"""
    plot_cobweb(output_path::String="cobweb_map.png")

Generates a Cobweb plot from the x-coordinates of the Poincaré map.
"""
function plot_cobweb(output_path::String="cobweb_map.png")
    p = DuffingParams(A=0.02, ω=1.0, ϵ=0.05)
    
    # Callback
    is_out_of_domain(u, t, integrator) = abs(u[1]) > 10.0 || abs(u[2]) > 10.0
    terminate_cb = DiscreteCallback(is_out_of_domain, terminate!)
    
    # Callback
    is_out_of_domain(u, t, integrator) = abs(u[1]) > 10.0 || abs(u[2]) > 10.0
    terminate_cb = DiscreteCallback(is_out_of_domain, terminate!)
    
    # Ensemble parameters
    N_traj = 50
    T_max_periods = 50
    T = 2π / p.ω
    
    all_x_n = Vector{Float64}[]
    
    # Grid of ICs
    x_range = range(-0.5, 0.5, length=isqrt(N_traj))
    y_range = range(-0.5, 0.5, length=isqrt(N_traj))
    
    for x0 in x_range, y0 in y_range
        prob = ODEProblem(duffing!, [x0, y0], (0.0, T_max_periods * T), p)
        sol = solve(prob, Tsit5(), callback=terminate_cb, saveat=T)
        
        # Extract x_n
        x_n = [u[1] for u in sol.u if abs(u[1]) < 2.0]
        if length(x_n) > 1
            push!(all_x_n, x_n)
        end
    end
    
    if isempty(all_x_n)
        println("Cobweb simulation failed/empty.")
        return
    end
    
    fig = Figure(size = (600, 600))
    ax = Axis(fig[1, 1], title = "Cobweb Plot (Chaotic Ensemble)", xlabel="x_n", ylabel="x_{n+1}", limits=(-2, 2, -2, 2))
    
    # Plot y=x line
    lines!(ax, [-2, 2], [-2, 2], color=:gray, linestyle=:dash)
    
    # Plot Cobweb lines for a few trajectories
    for i in 1:min(5, length(all_x_n))
        x_n = all_x_n[i]
        points = Point2f[]
        for j in 1:length(x_n)-1
            push!(points, Point2f(x_n[j], x_n[j]))     # Start at diagonal
            push!(points, Point2f(x_n[j], x_n[j+1]))   # Go to function
            push!(points, Point2f(x_n[j+1], x_n[j+1])) # Go to diagonal
        end
        if !isempty(points)
            lines!(ax, points, color=:purple, linewidth=0.5, transparency=true)
        end
    end
    
    # Scatter points for ALL trajectories
    for x_n in all_x_n
        if length(x_n) > 1
            scatter!(ax, x_n[1:end-1], x_n[2:end], color=:black, markersize=2, transparency=true)
        end
    end
    
    save(output_path, fig)
end

"""
    plot_complexity_spectrum(output_path::String="complexity_spectrum.png")

Plots Lyapunov Exponent vs Lempel-Ziv Complexity as ω varies.
"""
function plot_complexity_spectrum(output_path::String="complexity_spectrum.png")
    ω_range = range(0.5, 1.5, length=50) # Reduced for speed
    
    lyaps = Float64[]
    lzs = Float64[]
    
    for ω_val in ω_range
        p = DuffingParams(A=0.02, ω=ω_val, ϵ=0.05)
        
        # Lyapunov
        λ = compute_lyapunov(p; T=500.0)
        push!(lyaps, λ)
        
        # Lempel-Ziv
        ds = CoupledODEs(duffing!, [0.1, 0.1], p)
        tr, t = trajectory(ds, 1000.0; Δt=0.1)
        lz = compute_lempel_ziv(tr[:, 1])
        push!(lzs, lz)
    end
    
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], title = "Complexity Spectrum", xlabel="Frequency ω")
    
    lines!(ax, ω_range, lyaps, color=:blue, label="Max Lyapunov Exponent")
    lines!(ax, ω_range, lzs, color=:orange, label="Lempel-Ziv Complexity")
    
    axislegend(ax)
    save(output_path, fig)
end

"""
    plot_markov_heatmap(output_path::String="markov_heatmap.png")

Generates Markov transition matrices for Periodic and Chaotic regimes.
"""
function plot_markov_heatmap(output_path::String="markov_heatmap.png")
    # Helper to compute matrix
    function get_markov_matrix(p)
        ds = CoupledODEs(duffing!, [0.1, 0.1], p)
        tr, t = trajectory(ds, 2000.0; Δt=0.1)
        x_data = tr[:, 1]
        symbols = [x < 0 ? 0 : 1 for x in x_data]
        
        M = zeros(2, 2)
        for i in 1:length(symbols)-1
            s_curr = symbols[i] + 1 # 1-based index
            s_next = symbols[i+1] + 1
            M[s_curr, s_next] += 1
        end
        
        # Normalize rows
        for i in 1:2
            row_sum = sum(M[i, :])
            if row_sum > 0
                M[i, :] ./= row_sum
            end
        end
        return M
    end
    
    # Periodic
    p_p = DuffingParams(A=0.02, ω=0.9, ϵ=0.05)
    M_p = get_markov_matrix(p_p)
    
    # Chaotic
    p_c = DuffingParams(A=0.02, ω=1.0, ϵ=0.05)
    M_c = get_markov_matrix(p_c)
    
    fig = Figure(resolution = (800, 400))
    
    ax1 = Axis(fig[1, 1], title = "Markov Matrix (Periodic)", xticks=(1:2, ["0", "1"]), yticks=(1:2, ["0", "1"]))
    heatmap!(ax1, M_p, colormap=:blues, colorrange=(0, 1))
    text!(ax1, "0->0: $(round(M_p[1,1], digits=2))", position=(1, 1), color=:black)
    
    ax2 = Axis(fig[1, 2], title = "Markov Matrix (Chaotic)", xticks=(1:2, ["0", "1"]), yticks=(1:2, ["0", "1"]))
    heatmap!(ax2, M_c, colormap=:reds, colorrange=(0, 1))
    
    save(output_path, fig)
end

end # module Maps
