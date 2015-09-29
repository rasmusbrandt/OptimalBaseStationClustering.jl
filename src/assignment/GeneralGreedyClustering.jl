##########################################################################
# Greedy base station clustering.
#
# A simple greedy clustering method based on path loss and transmit powers.
# In each step, the method finds the strongest interfering link, and matches
# the corresponding cells into a cluster, if that is IA feasible. This is
# done either by putting the offending BS into the cluster of the victim
# BS (GeneralGreedyClustering_Single), or by trying to merge the respective
# clusters, if possible (GeneralGreedyClustering_Multiple). The bound on the
# cluster sizes will be given by IA feasibility, and not by the overhead
# pre-log factor.

function GeneralGreedyClustering(channel, network, f, t, D)
    LargeScaleFadingCellAssignment!(channel, network)
    assignment = get_assignment(network)
    require_equal_num_MSs_per_cell(assignment)
    I = get_num_BSs(network); K = get_num_MSs(network); Kc = convert(Int, K/I)
    require_equal_num_BS_antennas(network); M = get_num_BS_antennas(network)[1]
    require_equal_num_MS_antennas(network); N = get_num_MS_antennas(network)[1]
    require_equal_num_streams(network); d = get_num_streams(network)[1]
    Ps = get_transmit_powers(network); sigma2s = get_receiver_noise_powers(network)
    static_params = I::Int, K::Int, Kc::Int, M::Int, N::Int, d::Int, Ps::Vector{Float64}, sigma2s::Vector{Float64}

    # Cluster assignment matrix
    partition_matrix = eye(Int, I, I)

    # Clustering metric
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    F = zeros(Float64, I, I)
    for i = 1:I; for j = 1:I
        if i == j
            F[i,j] = -Inf
            continue
        end

        # Sum interference from BS j to MSs served by BS i
        F[i,j] = sum([ log2(1 + (channel.large_scale_fading_factor[k,j]^2)*Ps[j]/sigma2s[k]) for k in served_MS_ids(i, assignment) ])
    end; end

    # Greedily build clusters based on strongest sum interference between cells
    num_objective_calculations = 0
    while !all(F .== -Inf)
        # Find strongest interfering link that is still active
        _, idx = findmax(F)
        i, j = ind2sub((I, I), idx)

        # Join clusters of BS i and BS j
        i_cluster = find(partition_matrix[i,:] .== 1)
        j_cluster = find(partition_matrix[j,:] .== 1)

        if length(i_cluster) + length(j_cluster) <= D
            partition_matrix[i_cluster,j_cluster] = 1
            partition_matrix[j_cluster,i_cluster] = 1
        end

        # Never consider this link again
        F[i,j] = -Inf
    end
    partition = Partition(partition_matrix)
    a = restricted_growth_string(partition_matrix)
    desired_powers, interfering_powers =
        precompute_powers(channel, assignment, static_params)
    throughputs = zeros(Float64, K, d)
    for block in partition.blocks
        outside_all_BSs = setdiff(IntSet(1:I), block.elements)
        cluster_size = length(block.elements)
        for i in block.elements; for k in served_MS_ids(i, assignment)
            irreducible_interference_power = Float64(0)
            for j in outside_all_BSs
                irreducible_interference_power += interfering_powers[j,k]
            end
            rho = desired_powers[k]/(sigma2s[k] + irreducible_interference_power)
            throughputs[k,:] = t(cluster_size, rho)
        end; end
    end
    objective = f(throughputs)
    Lumberjack.info("GeneralGreedyClustering finished.",
        Dict(:objective => objective, :a => a))

    # Return results
    results = AssignmentResults()
    results["throughputs"] = throughputs
    results["a"] = a
    results["num_clusters"] = 1 + maximum(a)
    results["avg_cluster_size"] = avg_cluster_size(a)
    results["num_objective_calculations"] = num_objective_calculations
    return results
end
