# Modified Duffing Oscillator Analysis

A comprehensive Julia implementation for analyzing the dynamics of the modified Duffing oscillator system, featuring bifurcation analysis, chaos detection, and information theory measures.

## System Definition

The modified Duffing oscillator is governed by:

```
ẋ = y
ẏ = x - x³ - εy - μx²y + A cos(ωt)
```

**Fixed Parameter:** μ = -0.225

**Variable Parameters:**
- ε: damping parameter
- A: driving amplitude (0 for autonomous system)
- ω: driving frequency

## Core Infrastructure

The project is built on a robust core infrastructure (~800 lines) that makes implementing visualizations simple and maintainable.

### Module Structure

```
src/
├── core/
│   ├── system.jl         - ODE definitions and problem construction
│   ├── stability.jl      - Fixed point and eigenvalue analysis
│   ├── solvers.jl        - Optimized solver configurations (AMD ROCm GPU support)
│   ├── analysis.jl       - Dynamical systems analysis tools
│   ├── complexity.jl     - Information theory measures
│   ├── plotting.jl       - Unified visualization utilities
│   └── io.jl            - File management and caching
└── DuffingCore.jl       - Main module tying everything together
```

### Core Capabilities

#### 1. System Module (`system.jl`)
- `duffing_oscillator!(du, u, p, t)` - In-place ODE function
- `create_autonomous_problem(ε)` - Autonomous system (A=0)
- `create_driven_problem(A, ω, ε)` - Non-autonomous system
- Global constant μ = -0.225

#### 2. Stability Module (`stability.jl`)
- `find_fixed_points(ε)` - Returns [(0,0), (±1,0)]
- `compute_jacobian(x, y, ε)` - Jacobian matrix at point
- `compute_eigenvalues(x, y, ε)` - Eigenvalue analysis
- `classify_fixed_point(x, y, ε)` - Returns stability type
- `hopf_bifurcation_value()` - Analytical Hopf point (ε=0)
- `analyze_all_fixed_points(ε)` - Complete stability analysis

#### 3. Solvers Module (`solvers.jl`)
**GPU Support:** Optimized for AMD Radeon with 96GB VRAM using ROCm
- `get_high_precision_solver()` - Vern7 for bifurcation analysis
- `solve_high_precision(prob)` - Solve with tight tolerances (1e-9)
- `solve_transient(prob)` - Discard transient behavior
- `solve_ensemble_gpu(ensemble_prob)` - GPU-accelerated ensemble solving
- `create_saddle_stop_callback(saddle_point)` - For homoclinic orbits
- `check_gpu_availability()` - Display GPU status

**GPU Configuration:**
- Supports AMD ROCm with AMDGPU.jl
- Batch sizes optimized for 96GB VRAM
- Falls back gracefully to CPU threading if GPU unavailable

#### 4. Analysis Module (`analysis.jl`)
- `compute_lyapunov_spectrum(prob, params)` - Max Lyapunov exponent
- `compute_poincare_section(sol, period)` - Stroboscopic sampling
- `symbolize_trajectory(sol)` - Convert to binary sequence
- `extract_extrema(sol)` - For bifurcation diagrams
- `compute_period_from_poincare(points)` - Detect periodicity

#### 5. Complexity Module (`complexity.jl`)
- `lempel_ziv_complexity(binary_string)` - LZ complexity measure
- `markov_transition_matrix(binary_string)` - 2×2 transition matrix
- `compute_entropy(transition_matrix)` - Shannon entropy
- `analyze_symbolic_dynamics(binary_string)` - Complete analysis

#### 6. Plotting Module (`plotting.jl`)
- `setup_theme()` - Consistent publication-ready styling
- `save_animation(fig, filename)` - Export as MP4
- `save_figure(fig, filename)` - High-res PNG export
- `create_phase_portrait_axis(fig, position)` - Standard phase space setup
- `create_bifurcation_axis(fig, position)` - Bifurcation diagram setup
- `create_dual_y_axis(fig, position)` - For complexity spectrum
- `colormap_for_basins()` - Basin colors and labels
- `get_stability_color(classification)` - Color by stability type

#### 7. IO Module (`io.jl`)
- `ensure_output_dirs()` - Create directory structure
- `cache_computation(compute_fn, filename)` - Cache expensive results
- `export_data(data, name)` - Export to JLD2 format
- `import_data(filepath)` - Load saved data
- `save_parameters(params, name)` - Save configurations

## Installation

### Prerequisites
- Julia 1.10+
- ROCm installation (for GPU acceleration on AMD hardware)

### Package Setup
```bash
cd DS
julia --project=.
```

```julia
using Pkg
Pkg.instantiate()  # Install all dependencies
```

### Dependencies
The project uses:
- DifferentialEquations.jl - ODE solving
- DynamicalSystems.jl - Dynamical systems analysis
- BifurcationKit.jl - Numerical continuation
- ComplexityMeasures.jl - Information theory
- GLMakie.jl - Interactive visualizations
- AMDGPU.jl - AMD GPU acceleration via ROCm
- JLD2.jl - Data serialization

## GPU Configuration

### AMD Radeon with 96GB VRAM

The codebase is optimized for your AMD Radeon GPU with ROCm:

**Advantages:**
- Massive parallelization for ensemble simulations (10,000+ trajectories)
- Large batch sizes (default: 1000, can go much higher)
- Fast computation for Arnold tongues and fractal basins

**Setup:**
1. Ensure ROCm is installed on your system
2. AMDGPU.jl will detect available GPU automatically
3. Set batch sizes according to problem:
   - Standard ensemble: `batch_size=1000`
   - Fractal basins (1M points): `batch_size=5000`
   - Arnold tongue (20K combinations): `batch_size=2000`

**Check Status:**
```julia
include("src/core/solvers.jl")
check_gpu_availability()
```

## Usage Examples

### Basic System Solving
```julia
include("src/core/system.jl")
include("src/core/solvers.jl")

# Autonomous system
prob = create_autonomous_problem(0.1, tspan=(0.0, 100.0))
sol = solve_high_precision(prob)

# Driven system
prob_driven = create_driven_problem(0.02, 1.0, 0.05, tspan=(0.0, 1000.0))
sol_driven = solve(prob_driven, get_standard_solver())
```

### Stability Analysis
```julia
include("src/core/stability.jl")

# Analyze all fixed points
results = analyze_all_fixed_points(0.1)

# Or print detailed report
print_stability_analysis(0.1)
```

### Symbolic Dynamics & Complexity
```julia
include("src/core/analysis.jl")
include("src/core/complexity.jl")

# Convert trajectory to binary sequence
binary_seq = symbolize_trajectory(sol, threshold=0.0)

# Compute complexity measures
analysis = analyze_symbolic_dynamics(binary_seq)
println("LZ Complexity: ", analysis.lz_complexity)
println("Entropy: ", analysis.entropy)
```

### GPU-Accelerated Ensemble
```julia
using DifferentialEquations

# Create ensemble problem
prob = create_driven_problem(0.02, 1.0, 0.05)

# Generate 10,000 initial conditions
function prob_func(prob, i, repeat)
    x0 = 1.0 + randn() * 1e-6
    y0 = randn() * 1e-6
    remake(prob, u0=[x0, y0])
end

ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)

# Solve on GPU (with 96GB VRAM, can handle large batches)
sol_ensemble = solve_ensemble_gpu(ensemble_prob,
                                  trajectories=10000,
                                  batch_size=2000)
```

## Project Timeline

### Phase 1: Core Infrastructure ✓ (Completed)
- [x] System definitions
- [x] Stability analysis
- [x] Solver configurations with AMD GPU support
- [x] Analysis utilities
- [x] Complexity measures
- [x] Plotting utilities
- [x] I/O management

### Phase 2: Visualization Scripts (Next)
**Track A: Autonomous Bifurcations** (11 hours)
- [ ] Vis #1: Hopf Scanner Animation (2h)
- [ ] Vis #2: Bifurcation Diagram (3h)
- [ ] Vis #3: Homoclinic Saddle-Loop (6h)

**Track B: Phase Space Chaos** (6.5 hours)
- [ ] Vis #4: Ensemble Divergence Animation (4h)
- [ ] Vis #5: Poincaré Map (1.5h)
- [ ] Vis #6: Cobweb Return Map (1h)

**Track C: Information Theory** (6 hours)
- [ ] Vis #7: Complexity Spectrum (4h)
- [ ] Vis #8: Markov Transition Heatmaps (2h)

**Track D: Global Dynamics** (14 hours)
- [ ] Vis #9: Fractal Basins of Attraction (8h)
- [ ] Vis #10: Arnold Tongue Sweep (6h)

## Output Structure

```
outputs/
├── animations/
│   ├── hopf_scanner.mp4
│   └── ensemble_divergence.mp4
├── figures/
│   ├── bifurcation_diagram.png
│   ├── homoclinic_loop.png
│   ├── poincare_map.png
│   ├── cobweb_map.png
│   ├── complexity_spectrum.png
│   ├── markov_heatmap.png
│   ├── fractal_basins.png
│   └── arnold_tongue.png
├── data/
│   └── (exported data files in JLD2 format)
└── cache/
    └── (cached computation results)
```

## Performance Notes

### With AMD GPU (96GB VRAM)
- **Ensemble simulations**: Minutes for 10K trajectories
- **Fractal basins** (1M points): ~2-4 hours
- **Arnold tongue** (20K parameter combinations): ~3-6 hours

### CPU Fallback
- **Ensemble simulations**: Tens of minutes for 10K trajectories
- **Fractal basins**: Days for full resolution
- **Arnold tongue**: Many hours

The 96GB VRAM allows for aggressive batching and massive parallelization.

## Mathematical Features

### Bifurcation Analysis
- Andronov-Hopf bifurcation at ε = 0
- Homoclinic bifurcation (limit cycle → saddle collision)
- Numerical continuation with BifurcationKit.jl

### Chaos Detection
- Positive Lyapunov exponents
- Strange attractors in Poincaré maps
- Exponential trajectory divergence

### Information Theory
- Lempel-Ziv complexity (randomness measure)
- Markov transition matrices
- Shannon entropy

### Global Dynamics
- Fractal basin boundaries
- Arnold tongues (resonance regions)
- Final state sensitivity

## References

- Technical Design Document V1 & V2 (included PDFs)
- Strogatz, S. H. (2015). *Nonlinear Dynamics and Chaos*
- Shilnikov, L. P., et al. *Methods of Qualitative Theory in Nonlinear Dynamics*

## License

Academic research project.

## Contact

For questions about the modified Duffing oscillator analysis or implementation details, refer to the technical design documents.
# p3
