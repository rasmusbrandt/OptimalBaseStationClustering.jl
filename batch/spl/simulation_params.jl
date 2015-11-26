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
const I = 16; const Kc = 2
const M = 8; const N = 2; const d = 1
const geography_width = 2000. # sqrt(I*BS_density)
const MS_serving_BS_distance = Nullable(150.) # geography_width/10. = 161.1854897735313

const SNR_dB = 20.

const Ndrops = 250
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
