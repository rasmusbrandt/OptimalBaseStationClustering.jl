#!/usr/bin/env julia

using CoordinatedPrecoding, OptimalBaseStationClustering, DistributedBaseStationClustering
using JLD

# Parameters
include("../simulation_params.jl")
const Is = 1:16

# General
simulation_params = Dict(
    "Ndrops" => Ndrops, "Nsim" => Nsim,
    "I" => 1, "Kc" => Kc,
    "M" => M, "N" => N, "d" => d,
    "geography_size" => (geography_width, geography_width),
    "MS_serving_BS_distance" => MS_serving_BS_distance,
    "precoding_methods" => [ NoPrecoding ],
    "assignment_methods" => [ BranchAndBoundClustering, GreedyClustering ],
    "aux_assignment_params" => Dict(
        "max_num_MSs_per_BS" => Kc,

        "GeneralBranchAndBoundClustering:branching_rule" => :dfs,
        "GeneralBranchAndBoundClustering:improve_initial_incumbent" => true,
        "GeneralBranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:store_evolution" => false,
    ),
    "independent_variable" => (set_average_SNRs_dB!, [SNR_dB]),
)

# Simulate and save
srand(725242)
results = zeros(Float64, length(Is), simulation_params["Ndrops"], 2)
for (idx, It) in enumerate(Is)
    network =
        setup_random_large_scale_network(It,
            simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
            num_streams=simulation_params["d"],
            geography_size=(sqrt(It*BS_density), sqrt(It*BS_density)),
            MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

    _, raw_assignment_results =
        simulate(network, simulation_params, loop_over=:assignment_methods)
    results[idx, :, 1] = [ sum(raw_assignment_results[ii_Ndrop]["BranchAndBoundClustering"]["num_bounded_nodes"]) for ii_Ndrop = 1:simulation_params["Ndrops"] ]
    results[idx, :, 2] = [ sum(raw_assignment_results[ii_Ndrop]["GreedyClustering"]["num_constraint_checks"]) for ii_Ndrop = 1:simulation_params["Ndrops"] ]
end

println("-- Saving results")
save("I.jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "Is", Is,
     "results", results)
