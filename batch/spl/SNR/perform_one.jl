#!/usr/bin/env julia

using CoordinatedPrecoding, OptimalBaseStationClustering, DistributedBaseStationClustering
using Compat, JLD, ArgParse
Lumberjack.remove_truck("default")

include(joinpath(dirname(@__FILE__), "../simulation_params.jl"))

function perform_one()
    # Get seed
    s = ArgParseSettings()
    @add_arg_table s begin
        "--seed", "-s"
            help = "RNG seed"
            arg_type = Int
            required = true
    end
    args = parse_args(s)
    seed = args["seed"]

    # Simulation
    srand(seed)
    simulation_params = Dict(
        "Ndrops" => fld(Ndrops, 25), "Nsim" => Nsim,
        "I" => I, "Kc" => Kc,
        "M" => M, "N" => N, "d" => d,
        "geography_size" => (geography_width, geography_width),
        "MS_serving_BS_distance" => MS_serving_BS_distance,
        "assignment_methods" => [
            BranchAndBoundClustering,
            GreedyClustering,

            NoClustering,
        ],
        "precoding_methods" => [ NoPrecoding ],
        "aux_network_params" => Dict(
            "beta_network_sdma" => beta_network_sdma,
            "num_coherence_symbols" => num_coherence_symbols,
        ),
        "aux_assignment_params" => Dict(
            "max_num_MSs_per_BS" => Kc,

            "GeneralBranchAndBoundClustering:branching_rule" => :dfs,
            "GeneralBranchAndBoundClustering:improve_initial_incumbent" => true,
            "GeneralBranchAndBoundClustering:max_abs_optimality_gap" => 0.,
            "GeneralBranchAndBoundClustering:max_rel_optimality_gap" => 0.,
            "GeneralBranchAndBoundClustering:store_evolution" => false,
        ),

        "independent_variable" => (set_average_SNRs_dB!, -10:5:40),
    )
    simulation_params["simulation_name"] = "SNR-seed$seed"
    network =
        setup_random_large_scale_network(simulation_params["I"],
            simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
            num_streams=simulation_params["d"],
            geography_size=simulation_params["geography_size"],
            MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])
    raw_precoding_results, raw_assignment_results =
        simulate(network, simulation_params, loop_over=:assignment_methods)

    # Save for publication plots
    save("$(simulation_params["simulation_name"]).jld",
         "simulation_params", clean_simulation_params_for_jld(simulation_params),
         "raw_assignment_results", raw_assignment_results)
end

perform_one()
