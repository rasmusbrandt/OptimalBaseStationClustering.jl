#!/usr/bin/env julia

##########################################################################
# run_SNR.jl
#
# Performance as a function of transmit power, comparing different
# cluster assignment methods.
##########################################################################

using CoordinatedPrecoding, OptimalBaseStationClustering, DistributedBaseStationClustering
using Compat, JLD

##########################################################################
# General settings
srand(725242)
start_time = @compat Libc.strftime("%Y%m%dT%H%M%S", time())

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
const geography_width = 2000. # sqrt(I*BS_density)
const MS_serving_BS_distance = Nullable(150.) # geography_width/10. = 161.1854897735313

const SNRs_dB = -10:5:60

const Ndrops = 1
const Nsim = 1

# Utility model and specialized branch and bound
f1 = 1/I - (M + Kc*(N + d))/num_coherence_symbols
f2 = -Kc*M/num_coherence_symbols
alpha(C) = (f1 + f2*C)*C
r1(SNR) = (1 - beta_network_sdma)*exp_times_E1(SNR)
t(C, SNR, SINR) = d*(alpha(C)*r1(SNR) + beta_network_sdma*exp_times_E1(SINR))
D = floor(Int, (M + N - d)/(Kc*d)) # from Liu2013
_, B = findmax([ t(b, 1., 1.) for b in 1:I ]) # Optimality of B is independent of SNR and SINR, for this utility model
f(ts) = sum(ts) # sum throughput

BranchAndBoundClustering(channel, network) = GeneralBranchAndBoundClustering(channel, network, f, t, D, B)
GreedyClustering(channel, network) = GeneralGreedyClustering(channel, network, f, t, D, B)

##########################################################################
# Simulation settings
simulation_params = Dict(
    "Ndrops" => Ndrops, "Nsim" => Nsim,
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

    "independent_variable" => (set_average_SNRs_dB!, SNRs_dB),
)
network =
    setup_random_large_scale_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"],
        geography_size=simulation_params["geography_size"],
        MS_serving_BS_distance=simulation_params["MS_serving_BS_distance"])

##########################################################################
# Run and process results
raw_precoding_results, raw_assignment_results =
    simulate(network, simulation_params, loop_over=:assignment_methods)
plot_params = Dict(
    "objective" => :sum,

    "methods" => Dict(
        "BranchAndBoundClustering" => [ ("throughputs",), ],
        "GreedyClustering" => [ ("throughputs",), ],
        "NoClustering" => [ ("throughputs",), ],
    )
)
_, processed_data_mean, _ = postprocess(raw_assignment_results, simulation_params, plot_params)

##########################################################################
# Plot
idp_vals = simulation_params["independent_variable"][2]

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

# Plot it
fig = PyPlot.figure()
ax = fig[:add_axes]((0.12,0.24,0.88-0.12,0.95-0.24))

ax[:plot](idp_vals, processed_data_mean["BranchAndBoundClustering"]["throughputs"],
    color=colours[:blue], linestyle="-", marker="o", markeredgecolor=colours[:blue], markevery=3,
    label="Branch and bound")
ax[:plot](idp_vals, processed_data_mean["GreedyClustering"]["throughputs"],
    color=colours[:green], linestyle="-", marker="^", markeredgecolor=colours[:green], markevery=3,
    label="Heuristic")
ax[:plot](idp_vals, processed_data_mean["NoClustering"]["throughputs"],
    color=colours[:pink], linestyle="-", marker="x", markeredgecolor=colours[:pink], markevery=3,
    label="No clustering")
ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel]("Sum t'put [nats/s/Hz]")
legend = ax[:legend](loc="upper left")
legend_frame = legend[:get_frame]()
PyPlot.setp(legend_frame, linewidth=0.5)
fig[:savefig]("SNR.eps")
