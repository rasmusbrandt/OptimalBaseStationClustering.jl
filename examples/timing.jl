#!/usr/bin/env julia

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
const I = 10; const Kc = 2
const M = 8; const N = 2; const d = 1
const geography_width = 2000. # sqrt(I*BS_density)
const MS_serving_BS_distance = Nullable(150.) # geography_width/10. = 161.1854897735313

const SNR_dB = 20.

const Ndrops = 1
const Nsim = 1

# Utility model and specialized branch and bound
f1 = 1/I - (M + Kc*(N + d))/num_coherence_symbols
f2 = -Kc*M/num_coherence_symbols
alpha(C) = (f1 + f2*C)*C
r1 = (1 - beta_network_sdma)*exp_times_E1(1 + 10^(SNR_dB/10))
t(C, SINR) = alpha(C)*r1 + beta_network_sdma*exp_times_E1(1 + SINR)
D = floor(Int, (M + N - d)/(Kc*d)) # from Liu2013
_, B = findmax([ t(b, 1.) for b in 1:I ]) # Optimality of B is independent of rho, for this utility model
f(ts) = sum(ts) # sum throughput

BranchAndBoundClustering(channel, network) = GeneralBranchAndBoundClustering(channel, network, f, t, D, B)
ExhaustiveSearchClustering(channel, network) = GeneralExhaustiveSearchClustering(channel, network, f, t, D, B)
GreedyClustering(channel, network) = GeneralGreedyClustering(channel, network, f, t, D, B)

# General
simulation_params = Dict(
    "Ndrops" => 1, "Nsim" => 1,
    "Ntest" => 10,
    "I" => I, "Kc" => Kc,
    "M" => M, "N" => N, "d" => d,
    "geography_size" => (geography_width, geography_width),
    "MS_serving_BS_distance" => MS_serving_BS_distance,
    "assignment_methods" => [ BranchAndBoundClustering, ExhaustiveSearchClustering, GreedyClustering ],
    "aux_assignment_params" => Dict(
        "max_num_MSs_per_BS" => Kc,

        "GeneralBranchAndBoundClustering:branching_rule" => :bfs,
        "GeneralBranchAndBoundClustering:improve_initial_incumbent" => true,
        "GeneralBranchAndBoundClustering:max_abs_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:max_rel_optimality_gap" => 0.,
        "GeneralBranchAndBoundClustering:store_evolution" => true
    ),
)

network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

timing(network, simulation_params, loop_over=:assignment_methods)
