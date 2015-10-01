#!/usr/bin/env julia

using CoordinatedPrecoding, OptimalBaseStationClustering, DistributedBaseStationClustering
using JLD

# Parameters
include("../simulation_params.jl")

# General
simulation_params = Dict(
    "Ndrops" => Ndrops, "Nsim" => Nsim,
    "I" => I, "Kc" => Kc,
    "M" => M, "N" => N, "d" => d,
    "geography_size" => (geography_width, geography_width),
    "MS_serving_BS_distance" => MS_serving_BS_distance,
    "assignment_methods" => [
        BranchAndBoundClustering,
        GreedyClustering,

        CoalitionFormationClustering_AttachOrSupplant,

        GrandCoalitionClustering,
        NoClustering,
    ],
    "precoding_methods" => [ NoPrecoding ],
    "aux_assignment_params" => Dict(
        "max_num_MSs_per_BS" => Kc,

        "GeneralBranchAndBoundClustering:branching_rule" => :dfs,
        "GeneralBranchAndBoundClustering:improve_initial_incumbent" => true,
        "GeneralBranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:store_evolution" => false,
    ),

    # Needed for CoalitionFormationClustering_AttachOrSupplant
    "aux_network_params" => Dict(
        "num_coherence_symbols" => num_coherence_symbols,
        "beta_network_sdma" => beta_network_sdma,
    ),
    "independent_variable" => (set_average_SNRs_dB!, -10:5:40),
)
