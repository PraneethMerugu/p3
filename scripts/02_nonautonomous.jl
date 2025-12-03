using ModifiedDuffingProject
using ModifiedDuffingProject.Visualization.Ensembles
using ModifiedDuffingProject.Visualization.Maps

println("Generating Figure 4: Smale Blob Animation...")
smale_blob("smale_blob.mp4")
println("Done.")

println("Generating Figure 5: Poincaré Maps...")
plot_poincare("poincare_map.png")
println("Done.")

println("Generating Figure 6: Cobweb Plot...")
plot_cobweb("cobweb_map.png")
println("Done.")

println("Generating Figure 7: Complexity Spectrum...")
plot_complexity_spectrum("complexity_spectrum.png")
println("Done.")

println("Generating Figure 8: Markov Heatmaps...")
plot_markov_heatmap("markov_heatmap.png")
println("Done.")
