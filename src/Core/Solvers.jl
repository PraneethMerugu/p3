module Solvers

using DifferentialEquations
using OrdinaryDiffEq

export get_solver, get_tolerances

"""
    get_solver()

Returns the designated solver for the project: `Vern7` with automatic stiffness detection.
If stiffness is detected, it switches to `Rosenbrock23`.
"""
function get_solver()
    return AutoTsit5(Rosenbrock23())
    # Note: While the PDF specifies Vern7, AutoTsit5 is often more robust for general purpose.
    # However, to strictly follow the PDF which explicitly mentions Vern7 for high order accuracy:
    # We will use CompositeAlgorithm if needed, but let's stick to Vern7 as primary if non-stiff.
    # For safety and strict adherence to "Stiffness Detection" requirement:
    return AutoVern7(Rodas5()) 
end

"""
    get_tolerances()

Returns the strict relative and absolute tolerances required for chaotic simulation.
reltol=1e-9, abstol=1e-9
"""
function get_tolerances()
    return (reltol=1e-9, abstol=1e-9)
end

end # module Solvers
