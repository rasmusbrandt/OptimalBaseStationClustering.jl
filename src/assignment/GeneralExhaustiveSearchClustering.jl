##########################################################################
# Optimal base station clustering based on exhaustive search.
#
# All possible restricted growth strings (and thus set partitions) are
# enumerated, and the best is picked.

function GeneralExhaustiveSearchClustering(channel, network, f, t, D, B = 0)
    # Static params
    LargeScaleFadingCellAssignment!(channel, network)
    assignment = get_assignment(network)
    require_equal_num_MSs_per_cell(assignment)
    I = get_num_BSs(network); K = get_num_MSs(network); Kc = convert(Int, K/I)
    require_equal_num_BS_antennas(network); M = get_num_BS_antennas(network)[1]
    require_equal_num_MS_antennas(network); N = get_num_MS_antennas(network)[1]
    require_equal_num_streams(network); d = get_num_streams(network)[1]
    Ps = get_transmit_powers(network); sigma2s = get_receiver_noise_powers(network)
    static_params = I::Int, K::Int, Kc::Int, M::Int, N::Int, d::Int, D::Int, B::Int, Ps::Vector{Float64}, sigma2s::Vector{Float64}, assignment

    # Utility params
    utility_params = f::Function, t::Function

    # Warn if this will be slow...
    I >= 12 && Lumberjack.warn("ExhaustiveSearchClustering will be slow since I = $I.")

    # Exhaustive search over all partitions
    num_iters = 0
    best_objective = 0.; best_throughputs = Array(Float64, K, d)
    best_partition = Partition(collect(0:(I-1)))
    for partition in PartitionIterator(I)
        num_iters += 1

        # Calculate throughputs
        throughputs_ = throughputs(partition, channel, static_params, utility_params)

        objective = f(throughputs_)
        if objective >= best_objective
            best_objective = objective
            best_throughputs = throughputs_
            best_partition = partition
        end
    end
    a = restricted_growth_string(best_partition)
    Lumberjack.info("GeneralExhaustiveSearchClustering finished.",
        Dict(:best_objective => best_objective,
             :num_iters => num_iters,
             :a => a))

    # Return results
    results = AssignmentResults()
    results["throughputs"] = best_throughputs
    results["a"] = a
    results["num_clusters"] = reshape([1 + maximum(a)], 1, 1)
    results["avg_cluster_size"] = reshape([average_cluster_size(a)], 1, 1)
    results["num_iters"] = reshape([num_iters], 1, 1)
    return results
end
