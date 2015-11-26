#!/usr/bin/env julia

using CoordinatedPrecoding, OptimalBaseStationClustering, DistributedBaseStationClustering
using JLD

# Parameters
include("../simulation_params.jl")

# General
simulation_params = Dict(
    "Ndrops" => 1, "Nsim" => Nsim,
    "I" => I, "Kc" => Kc,
    "M" => M, "N" => N, "d" => d,
    "geography_size" => (geography_width, geography_width),
    "MS_serving_BS_distance" => MS_serving_BS_distance,
    "assignment_methods" => [ BranchAndBoundClustering ],
    "aux_assignment_params" => Dict(
        "max_num_MSs_per_BS" => Kc,

        "GeneralBranchAndBoundClustering:branching_rule" => :bfs,
        "GeneralBranchAndBoundClustering:improve_initial_incumbent" => true,
        "GeneralBranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:store_evolution" => true
    ),
)

# Simulate and save
srand(338822)
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
set_average_SNRs_dB!(network, SNR_dB)

raw_results =
    simulate_assignment_convergence(network, simulation_params)

println("-- Saving results")
save("convergence.jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
