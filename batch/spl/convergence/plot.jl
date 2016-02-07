#!/usr/bin/env julia

using CoordinatedPrecoding, OptimalBaseStationClustering
using JLD, LaTeXStrings

sim_name = "convergence"
data = load("$(sim_name).jld")

plot_params_bounds = Dict(
    "objective" => :none,

    "methods" => Dict(
        "BranchAndBoundClustering" => [
            ("upper_bound_evolution",),
            ("lower_bound_evolution",),
        ],
    )
)
_, processed_data_bounds_mean, _ = postprocess_assignment_convergence(data["raw_results"], data["simulation_params"], plot_params_bounds)

plot_params_fathom = Dict(
    "objective" => :none,

    "methods" => Dict(
        "BranchAndBoundClustering" => [
            ("num_iters",),
            ("num_bounded_nodes",),
            ("fathoming_evolution",),
        ],
    )
)
_, processed_data_fathom_mean, _ = postprocess_assignment_convergence(data["raw_results"], data["simulation_params"], plot_params_fathom)

# 8-class Set1
colours = Dict(
    :red => "#e41a1c",
    :blue => "#377eb8",
    :green => "#4daf4a",
    :purple => "#984ea3",
    :orange => "#ff7f00",
    :yellow => "#ffff33",
    :brown => "#a65628",
    :pink => "#f781bf",
)

# Plot defaults
PyPlot.rc("lines", linewidth=1., markersize=3, markeredgewidth=0.5)
PyPlot.rc("font", size=6, family="serif", serif="Computer Modern Sans Serif")
PyPlot.rc("text", usetex=true)
PyPlot.rc("text.latex", preamble="\\usepackage{amsmath}")
PyPlot.rc("axes", linewidth=0.5, labelsize=6)
PyPlot.rc("xtick", labelsize=6)
PyPlot.rc("ytick", labelsize=6)
PyPlot.rc("legend", fancybox=true, fontsize=6)
PyPlot.rc("figure", figsize=(3.50,1.2), dpi=125)

len = length(processed_data_bounds_mean["BranchAndBoundClustering"]["upper_bound_evolution"])

##########################################################################
# bounds
plot_name = "bounds"
fig = PyPlot.figure()
ax = fig[:add_axes]((0.12,0.24,0.88-0.12,0.95-0.24))

ax[:plot](1:len, processed_data_bounds_mean["BranchAndBoundClustering"]["upper_bound_evolution"],
    color=colours[:blue], linestyle=":", marker="v", markeredgecolor=colours[:blue], markevery=50,
    label="Best upper bound")
ax[:plot](1:len, processed_data_bounds_mean["BranchAndBoundClustering"]["lower_bound_evolution"],
    color=colours[:blue], linestyle="-", marker="^", markeredgecolor=colours[:blue], markevery=50,
    label="Best lower bound (incumbent)")
ax[:set_xlabel]("Iterations")
ax[:set_ylabel]("Sum t'put [nats/s/Hz]")
ax[:set_xlim](0,201)
ax[:set_ylim](30,40)
legend = ax[:legend](loc="upper right")
legend_frame = legend[:get_frame]()
PyPlot.setp(legend_frame, linewidth=0.5)
fig[:savefig]("$(sim_name)_$(plot_name).eps")

##########################################################################
# fathom
plot_name = "fathom"
fig = PyPlot.figure()
ax1 = fig[:add_axes]((0.12,0.24,0.88-0.12,0.95-0.24))

tree_size = OptimalBaseStationClustering.subtree_size(1, 1, data["simulation_params"]["I"])
fathomed_subtree_sizes_cum = cumsum(processed_data_fathom_mean["BranchAndBoundClustering"]["fathoming_evolution"])
len = length(fathomed_subtree_sizes_cum)

line1 = ax1[:plot](1:len, fathomed_subtree_sizes_cum, color=colours[:blue], linestyle="-", marker="o", markeredgecolor=colours[:blue], markevery=50)
ax1[:set_xlabel]("Iterations")
ax1[:set_ylabel]("\\# of nodes pruned")
ax1[:set_yscale]("log")

ax2 = ax1[:twinx]()
line2 = ax2[:plot](1:len, 100*fathomed_subtree_sizes_cum/tree_size, color=colours[:blue], linestyle="--", marker="d", markeredgecolor=colours[:blue], markevery=50)
ax2[:plot](len, 100*fathomed_subtree_sizes_cum[end]/tree_size, color=colours[:blue], linestyle="--", marker="d", markeredgecolor=colours[:blue])
ax2[:set_ylabel]("\\% of tree pruned")
ax2[:set_ylim](-5,105)
ax2[:set_xlim](0,201)

legend = ax1[:legend]([ line1[1], line2[1] ], [ "Absolute number", "Percentage" ], bbox_to_anchor=(0.8, 0.6))
legend_frame = legend[:get_frame]()
PyPlot.setp(legend_frame, linewidth=0.5)

fig[:savefig]("$(sim_name)_$(plot_name).eps")

max_cluster_size = (data["simulation_params"]["M"] + data["simulation_params"]["N"] - data["simulation_params"]["d"])/(data["simulation_params"]["Kc"]*data["simulation_params"]["d"])
stirlings = collect(CoordinatedPrecoding.Stirling2NumberIterator(data["simulation_params"]["I"]))
naive_search_complexity = sum(stirlings[ceil(Int, data["simulation_params"]["I"]/max_cluster_size):end])
num_bounded_nodes = processed_data_fathom_mean["BranchAndBoundClustering"]["num_bounded_nodes"][1]

println("Number of iterations: ", len)
println("Number of nodes bounded: ", num_bounded_nodes)
println("Fraction of search tree bounded: ", num_bounded_nodes/tree_size)
println("Fraction of search tree pruned: ", 1 - num_bounded_nodes/tree_size)
println("Complexity of naive search: ", naive_search_complexity)
println("Complexity reduction: ", naive_search_complexity/num_bounded_nodes)
