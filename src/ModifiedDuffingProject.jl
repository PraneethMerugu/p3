module ModifiedDuffingProject

# Include Core Modules
include("Core/System.jl")
include("Core/Solvers.jl")

# Export Core
using .System
using .Solvers
export DuffingParams, duffing!
export get_solver, get_tolerances

# Analysis Modules
include("Analysis/Stability.jl")
include("Analysis/Bifurcation.jl")
include("Analysis/Chaos.jl")

# Visualization Modules
include("Visualization/PhaseSpace.jl")
include("Visualization/Diagrams.jl")
include("Visualization/Manifolds.jl")
include("Visualization/Ensembles.jl")
include("Visualization/Maps.jl")

# Create Namespaces for Scripts
module Analysis
    using ..Stability
    using ..Bifurcation
    using ..Chaos
    export Stability, Bifurcation, Chaos
end

module Visualization
    using ..PhaseSpace
    using ..Diagrams
    using ..Manifolds
    using ..Ensembles
    using ..Maps
    export PhaseSpace, Diagrams, Manifolds, Ensembles, Maps
end

export Analysis, Visualization


end # module ModifiedDuffingProject
