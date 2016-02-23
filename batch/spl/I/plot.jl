#!/usr/bin/env julia

using CoordinatedPrecoding, DistributedBaseStationClustering
using Compat, JLD, LaTeXStrings

sim_name = "I"
data = load("$(sim_name).jld")

# Naive exhaustive search
max_cluster_size = (data["simulation_params"]["M"] + data["simulation_params"]["N"] - data["simulation_params"]["d"])/(data["simulation_params"]["Kc"]*data["simulation_params"]["d"])
exhaustive_search_complexity = zeros(Int, length(data["Is"]))
worst_case_bnb_complexity = zeros(Int, length(data["Is"]))
for i in data["Is"]
    stirlings = collect(CoordinatedPrecoding.Stirling2NumberIterator(i))
    exhaustive_search_complexity[i] = sum(stirlings[cld(i, max_cluster_size):end])
    worst_case_bnb_complexity[i] = sum(exhaustive_search_complexity[1:i])
end

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
PyPlot.rc("legend", fancybox=true, fontsize=5)
PyPlot.rc("figure", figsize=(3.50,1.2), dpi=125)

# Plot it
fig = PyPlot.figure()
ax = fig[:add_axes]((0.12,0.24,0.88-0.12,0.95-0.24))

ax[:plot](data["Is"], worst_case_bnb_complexity,
    color=colours[:blue], linestyle=":", marker="x", markeredgecolor=colours[:blue], markevery=4,
    label="Branch and bound (worst case)")
ax[:plot](data["Is"], exhaustive_search_complexity,
    color=colours[:orange], linestyle="-", marker="v", markeredgecolor=colours[:orange], markevery=4,
    label="Exhaustive search (exact)")
ax[:plot](data["Is"], mean(data["results"][:,:,1], 2),
    color=colours[:blue], linestyle="-", marker="o", markeredgecolor=colours[:blue], markevery=4,
    label="Branch and bound (average)")
ax[:plot](data["Is"], mean(data["results"][:,:,2], 2),
    color=colours[:green], linestyle="-", marker="^", markeredgecolor=colours[:green], markevery=4,
    label="Heuristic (average)")
ax[:set_xlabel](L"Number of BSs $I$")
ax[:set_ylabel]("\\# iterations")
ax[:set_yscale]("log")
legend = ax[:legend](loc="upper left")
legend_frame = legend[:get_frame]()
PyPlot.setp(legend_frame, linewidth=0.5)
fig[:savefig]("I.eps")
