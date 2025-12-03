println("Starting Full Generation Pipeline...")
push!(LOAD_PATH, joinpath(@__DIR__, "../src"))


include("01_autonomous.jl")
include("02_nonautonomous.jl")
include("03_global_dynamics.jl")

println("All visualizations generated successfully.")
