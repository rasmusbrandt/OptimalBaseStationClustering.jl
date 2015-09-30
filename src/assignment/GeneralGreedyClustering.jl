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
    # Static params
    LargeScaleFadingCellAssignment!(channel, network)
    assignment = get_assignment(network)
    require_equal_num_MSs_per_cell(assignment)
    I = get_num_BSs(network); K = get_num_MSs(network); Kc = convert(Int, K/I)
    require_equal_num_BS_antennas(network); M = get_num_BS_antennas(network)[1]
    require_equal_num_MS_antennas(network); N = get_num_MS_antennas(network)[1]
    require_equal_num_streams(network); d = get_num_streams(network)[1]
    Ps = get_transmit_powers(network); sigma2s = get_receiver_noise_powers(network)
    static_params = I::Int, K::Int, Kc::Int, M::Int, N::Int, d::Int, Ps::Vector{Float64}, sigma2s::Vector{Float64}, assignment

    # Scenario params
    _, B = findmax([ t(b, 1., 1.) for b in 1:I ])
    scenario_params = f::Function, t::Function, D::Int, B::Int

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
    num_constraint_checks = 0
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
        num_constraint_checks += 1

        # Never consider this link again
        F[i,j] = -Inf
    end
    a = restricted_growth_string(partition_matrix)
    partition = Partition(partition_matrix)
    throughputs_ = throughputs(partition, channel, static_params, scenario_params)
    Lumberjack.info("GeneralGreedyClustering finished.", Dict(:objective => f(throughputs_), :a => a))

    # Return results
    results = AssignmentResults()
    results["throughputs"] = throughputs_
    results["a"] = a
    results["num_clusters"] = reshape([1 + maximum(a)], 1, 1)
    results["average_cluster_size"] = reshape([average_cluster_size(a)], 1, 1)
    results["num_constraint_checks"] = reshape([num_constraint_checks], 1, 1)
    return results
end
