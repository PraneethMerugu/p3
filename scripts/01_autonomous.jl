using ModifiedDuffingProject
using ModifiedDuffingProject.Visualization.PhaseSpace
using ModifiedDuffingProject.Visualization.Diagrams
using ModifiedDuffingProject.Visualization.Manifolds

println("Generating Figure 1: Hopf Scanner Animation...")
hopf_scanner("hopf_scanner.mp4")
println("Done.")

println("Generating Figure 2: Global Bifurcation Diagram...")
plot_bifurcation_diagram("bifurcation_diagram.png")
println("Done.")

println("Generating Figure 3: Homoclinic Saddle-Loop...")
plot_homoclinic_loop("homoclinic_loop.png")
println("Done.")
