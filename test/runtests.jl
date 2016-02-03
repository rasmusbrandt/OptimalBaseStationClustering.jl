using Base.Test
using CoordinatedPrecoding, OptimalBaseStationClustering, DistributedBaseStationClustering

const fc = 2e9 # GHz
const Wc = 300e3 # kHz
const c = 300e6 # m/s
const λ = c/fc # m
const v_kmh = 30 # km/h
const v = v_kmh*(1e3/3600) # m/s
const fd = v/(λ*Wc)
const num_coherence_symbols = 1/(2*fd)
const beta_network_sdma = 0.5

const design_ISD = 500
const BS_density = sqrt(3)/2*design_ISD^2 # BSs per m^2, for hexagonal cells
const I = 8; const Kc = 2
const M = 8; const N = 2; const d = 1
const geography_width = sqrt(I*BS_density) # sqrt(I*BS_density)
const MS_serving_BS_distance = Nullable(150.) # geography_width/10. = 161.1854897735313

const SNR_dB = 30

const Ndrops = 100

# Utility model and specialized branch and bound
f1 = 1/I - (M + Kc*(N + d))/num_coherence_symbols
f2 = -Kc*M/num_coherence_symbols
alpha(C) = (f1 + f2*C)*C
r1(SNR) = (1 - beta_network_sdma)*exp_times_E1(SNR)
t(C, SNR, SINR) = d*(alpha(C)*r1(SNR) + beta_network_sdma*exp_times_E1(SINR))
D = floor(Int, (M + N - d)/(Kc*d)) # from Liu2013
_, B = findmax([ t(b, 1., 1.) for b in 1:I ]) # Optimality of B is independent of rho, for this utility model
f(ts) = sum(ts) # sum throughput

BranchAndBoundClustering(channel, network) = GeneralBranchAndBoundClustering(channel, network, f, t, D, B)
ExhaustiveSearchClustering(channel, network) = GeneralExhaustiveSearchClustering(channel, network, f, t, D, B)
GreedyClustering(channel, network) = GeneralGreedyClustering(channel, network, f, t, D, B)

# General
simulation_params = Dict(
    "Ndrops" => Ndrops, "Nsim" => 1,
    "I" => I, "Kc" => Kc,
    "M" => M, "N" => N, "d" => d,
    "geography_size" => (geography_width, geography_width),
    "MS_serving_BS_distance" => MS_serving_BS_distance,
    "assignment_methods" => [ BranchAndBoundClustering, ExhaustiveSearchClustering, GreedyClustering ],
    "precoding_methods" => [ NoPrecoding ],
    "independent_variable" => (set_average_SNRs_dB!, [SNR_dB]),
)

settings = [
    Dict("max_num_MSs_per_BS" => Kc,
        "GeneralBranchAndBoundClustering:branching_rule" => :bfs,
        "GeneralBranchAndBoundClustering:improve_initial_incumbent" => true,
        "GeneralBranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:store_evolution" => false
    ),
    Dict("max_num_MSs_per_BS" => Kc,
        "GeneralBranchAndBoundClustering:branching_rule" => :dfs,
        "GeneralBranchAndBoundClustering:improve_initial_incumbent" => true,
        "GeneralBranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:store_evolution" => false
    ),
    Dict("max_num_MSs_per_BS" => Kc,
        "GeneralBranchAndBoundClustering:branching_rule" => :bfs,
        "GeneralBranchAndBoundClustering:improve_initial_incumbent" => false,
        "GeneralBranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:store_evolution" => false
    ),
    Dict("max_num_MSs_per_BS" => Kc,
        "GeneralBranchAndBoundClustering:branching_rule" => :dfs,
        "GeneralBranchAndBoundClustering:improve_initial_incumbent" => false,
        "GeneralBranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:store_evolution" => false
    ),
]

for s in settings
    simulation_params["aux_assignment_params"] = s

    network =
        setup_random_large_scale_network(simulation_params["I"],
            simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
            num_streams=simulation_params["d"],
            geography_size=simulation_params["geography_size"],
            MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

    _, raw_assignment_results =
        simulate(network, simulation_params, loop_over=:assignment_methods)

    for i in eachindex(raw_assignment_results.simulation_results)
        b1 = raw_assignment_results.simulation_results[i]["BranchAndBoundClustering"]["throughputs"]
        b2 = raw_assignment_results.simulation_results[i]["BranchAndBoundClustering"]["a"]
        e1 = raw_assignment_results.simulation_results[i]["ExhaustiveSearchClustering"]["throughputs"]
        e2 = raw_assignment_results.simulation_results[i]["ExhaustiveSearchClustering"]["a"]
        g = raw_assignment_results.simulation_results[i]["GreedyClustering"]["throughputs"]
        @test_approx_eq b1 e1
        @test all(b2 .== e2)
        @test f(g) <= f(e1)
    end
end
