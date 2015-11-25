#!/usr/bin/env julia

using CoordinatedPrecoding, DistributedBaseStationClustering
using Compat, JLD

plot_params = Dict(
    "objective" => :sum,

    "methods" => Dict(
        "BranchAndBoundClustering" => [
            ("throughputs",),
        ],

        "GreedyClustering" => [
            ("throughputs",),
        ],

        "GrandCoalitionClustering" => [
            ("throughputs",),
        ],

        "NoClustering" => [
            ("throughputs",),
        ],
    )
)

##########################################################################
# Load data
sim_name = "SNR"
seeds = split(readchomp("../seeds.txt"), '\n'); Nseeds = length(seeds)
file_names = [ "$(sim_name)-seed$n.jld" for n in seeds ]

# Load first
println("Initial loading from $(file_names[1])")
data = load(file_names[1])
simulation_params = data["simulation_params"]; Ndrops = simulation_params["Ndrops"]
assn_siz = size(data["raw_assignment_results"].simulation_results)

raw_assignment_results = CoordinatedPrecoding.MultipleSimulationResults(Nseeds*assn_siz[1], assn_siz[2:end]...)

for (file_idx, file_name) in enumerate(file_names)
    println("Loading from $(file_name)")
    data = load(file_name)
    raw_assignment_results.simulation_results[(file_idx - 1)*Ndrops + 1:file_idx*Ndrops, :, :, :] = data["raw_assignment_results"].simulation_results
    data = []
end

_, processed_data_mean, _ = postprocess(raw_assignment_results, simulation_params, plot_params)

##########################################################################
# Save merged and post processed data
println("-- Saving merged results")
save("$(sim_name).jld",
     "simulation_params", simulation_params,
     "processed_data_mean", processed_data_mean)
