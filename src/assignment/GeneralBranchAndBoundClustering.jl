##########################################################################
# Optimal base station clustering based on branch and bound.
#
# A branch and bound tree is developed which describes all possible
# restricted growth strings that describe all possible base station
# clusters.

function GeneralBranchAndBoundClustering(channel, network, f, t, D, B = 0)
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

    # Algorithm parameters
    aux_params = get_aux_assignment_params(network)
    @defaultize_param! aux_params "GeneralBranchAndBoundClustering:branching_rule" :dfs
    branching_rule = aux_params["GeneralBranchAndBoundClustering:branching_rule"]
    @defaultize_param! aux_params "GeneralBranchAndBoundClustering:improve_initial_incumbent" true
    improve_initial_incumbent = aux_params["GeneralBranchAndBoundClustering:improve_initial_incumbent"]
    @defaultize_param! aux_params "GeneralBranchAndBoundClustering:max_abs_optimality_gap" 0.
    max_abs_optimality_gap = aux_params["GeneralBranchAndBoundClustering:max_abs_optimality_gap"]
    @defaultize_param! aux_params "GeneralBranchAndBoundClustering:max_rel_optimality_gap" 0.
    max_rel_optimality_gap = aux_params["GeneralBranchAndBoundClustering:max_rel_optimality_gap"]
    @defaultize_param! aux_params "GeneralBranchAndBoundClustering:store_evolution" false
    store_evolution = aux_params["GeneralBranchAndBoundClustering:store_evolution"]

    # Lumberjack.debug("GeneralBranchAndBoundClustering started.")

    # Precompute desired and interfering powers
    desired_powers, interfering_powers = compute_powers(channel, static_params)
    scratch = Array(Float64, I)

    # Initial bounds
    best_upper_bound = Inf
    incumbent_throughputs, incumbent_a, incumbent_objective =
        initial_incumbent(channel, network, static_params, utility_params, improve_initial_incumbent)

    # Eager branch and bound
    num_iters = 0; num_bounded_nodes = 0
    abs_conv_crit = 0.; premature_ending = false
    lower_bound_evolution = Float64[]; upper_bound_evolution = Float64[]; fathoming_evolution = Int[]
    live = initialize_live(channel, network, static_params, utility_params,
                           desired_powers, interfering_powers)
    while length(live) > 0
        num_iters += 1

        # Select next node to be processed.
        if branching_rule == :bfs
            # Best first, i.e. the highest (best) upper bound
            parent = Base.Collections.heappop!(live, Base.Order.Reverse)
            best_upper_bound = parent.upper_bound
        elseif branching_rule == :dfs
            # Depth first.
            best_upper_bound = maximum([ node.upper_bound for node in live ])
            parent = pop!(live)
        end

        if store_evolution
            # Store bound evolution per iteration
            push!(lower_bound_evolution, incumbent_objective)
            push!(upper_bound_evolution, best_upper_bound)
        end

        # Check convergence
        abs_conv_crit = best_upper_bound - incumbent_objective
        rel_conv_crit = abs_conv_crit/incumbent_objective
        if abs_conv_crit <= max_abs_optimality_gap || rel_conv_crit <= max_rel_optimality_gap
            # Lumberjack.debug("Converged.", Dict(:abs_conv_crit => abs_conv_crit, :max_abs_optimality_gap => max_abs_optimality_gap))
            premature_ending = true
            break
        end

        fathomed_subtree_size = 0
        for child in branch(parent)
            throughput_bounds =
                bound!(child, channel, network, static_params, utility_params,
                       desired_powers, interfering_powers, scratch)
            num_bounded_nodes += 1

            # Worthwhile investigating this subtree/leaf more?
            if child.upper_bound > incumbent_objective
                if is_leaf(child, I)
                    # For leaves, the upper bound is tight. Thus, we
                    # have found a new incumbent!
                    incumbent_objective = child.upper_bound
                    incumbent_throughputs = throughput_bounds
                    incumbent_a = child.a

                    # Lumberjack.debug("Found new incumbent solution.", Dict(:node => child, :incumbent_objective => incumbent_objective))
                else
                    # Lumberjack.debug("Keeping node since upper bound is above incumbent value.", Dict(:node => child, :incumbent_objective => incumbent_objective))
                    if branching_rule == :bfs
                        Base.Collections.heappush!(live, child, Base.Order.Reverse)
                    else
                        push!(live, child)
                    end
                end
            else
                # Lumberjack.debug("Discarding node since upper bound is below incumbent value.", Dict(:node => child, :incumbent_objective => incumbent_objective))

                if store_evolution
                    fathomed_subtree_size += subtree_size(child, I) - 1 # minus one since we already explored child
                end
            end
        end

        if store_evolution
            push!(fathoming_evolution, fathomed_subtree_size)
        end
    end

    # Add remaining subtrees that were implicitly fathomed
    if store_evolution
        push!(fathoming_evolution, sum([ subtree_size(node, I) for node in live ]))
    end

    # Did we find the global optimum?
    if (abs_conv_crit < 0.) || !premature_ending
        abs_conv_crit = 0.
    end

    Lumberjack.info("GeneralBranchAndBoundClustering finished.",
        Dict(:objective => f(incumbent_throughputs),
             :num_bounded_nodes => num_bounded_nodes,
             :abs_conv_crit => abs_conv_crit,
             :max_abs_optimality_gap => max_abs_optimality_gap,
             :a => incumbent_a))

    # Store a for next branch and bound run
    set_aux_network_param!(network, incumbent_a, "GeneralBranchAndBoundClustering:cache:optimal_a")

    # Return results
    results = AssignmentResults()
    results["throughputs"] = incumbent_throughputs
    results["a"] = incumbent_a
    results["num_clusters"] = reshape([1 + maximum(incumbent_a)], 1, 1)
    results["average_cluster_size"] = reshape([average_cluster_size(incumbent_a)], 1, 1)
    results["num_iters"] = reshape([num_iters], 1, 1)
    results["num_bounded_nodes"] = reshape([num_bounded_nodes], 1, 1)
    results["lower_bound_evolution"] = reshape(lower_bound_evolution, (1, 1, length(lower_bound_evolution)))
    results["upper_bound_evolution"] = reshape(upper_bound_evolution, (1, 1, length(upper_bound_evolution)))
    results["fathoming_evolution"] = reshape(fathoming_evolution, (1, 1, length(fathoming_evolution)))
    return results
end

function initial_incumbent(channel, network,
    static_params, utility_params, improve_initial_incumbent)

    I, K, Kc, M, N, d, D, B, Ps, sigma2s, assignment = static_params
    f, t = utility_params
    incumbent_throughputs = zeros(Float64, K, d)
    incumbent_a = zeros(Int, I)
    incumbent_objective = 0.

    if improve_initial_incumbent
        # Solution from heuristic
        heuristic_results = GeneralGreedyClustering(channel, network, f, t, D)
        heuristic_objective = f(heuristic_results["throughputs"])
        if heuristic_objective > incumbent_objective
            incumbent_throughputs = heuristic_results["throughputs"]
            incumbent_a = heuristic_results["a"]
            incumbent_objective = heuristic_objective
        end
        # Lumberjack.debug("Potential incumbent throughputs calculated.", Dict(:incumbent_objective => incumbent_objective))

        # Optimal solution from previous branch and bound round
        if has_aux_network_param(network, "GeneralBranchAndBoundClustering:cache:optimal_a")
            previous_a = get_aux_network_param(network, "GeneralBranchAndBoundClustering:cache:optimal_a")
            previous_partition = Partition(previous_a)
            previous_throughputs = throughputs(previous_partition, channel, static_params, utility_params)
            previous_objective = f(previous_throughputs)
            if previous_objective > incumbent_objective
                incumbent_throughputs = previous_throughputs
                incumbent_a = previous_a
                incumbent_objective = previous_objective
            end
        end
    end

    incumbent_throughputs, incumbent_a, incumbent_objective
end

type BranchAndBoundNode
    a::Vector{Int} # restricted growth string describing the (partial) partition
    upper_bound::Float64
end

# For sorting of the list of live nodes
Base.isless(N1::BranchAndBoundNode, N2::BranchAndBoundNode) = (N1.upper_bound < N2.upper_bound)

# Helper functions
m(node) = 1 + maximum(node.a)
depth(node) = size(node.a, 1)
is_leaf(node, I) = (depth(node) == I)
num_children(node) = 1 + m(node)

subtree_size(node, I) = subtree_size(depth(node), m(node), I)
function subtree_size(depth, m, I)
    if depth == I
        return 1
    else
        return 1 + m*subtree_size(depth+1, m, I) + subtree_size(depth+1, m+1 ,I)
    end
end

# Initialize the live structure by creating the root node.
function initialize_live(channel, network, static_params, utility_params,
                         desired_powers, interfering_powers)
    root = BranchAndBoundNode([0], Inf)
    bound!(root, channel, network, static_params, utility_params,
           desired_powers, interfering_powers,
           zeros(Float64, get_num_BSs(network)))
    # Lumberjack.debug("Root created.", Dict(:node => root))
    return [ root ]
end

# Branch a node by creating a number of descendants, corresponding to putting
# the BS at this depth in different clusters. Branched nodes inherit the
# parent bound, until their bounds are updated.
function branch(parent)
    num_children_ = num_children(parent)
    children = Array(BranchAndBoundNode, num_children_)

    for p = 1:num_children_
        child_a = push!(copy(parent.a), p - 1)
        child = BranchAndBoundNode(child_a, parent.upper_bound)
        children[p] = child

        # Lumberjack.debug("Branching.", Dict(:node => child))
    end

    return children
end

# Bound works by optimistically removing interference for unclustered BSs.
function bound!(node, channel, network, static_params, utility_params,
    desired_powers::Vector{Float64}, interfering_powers::Matrix{Float64},
    scratch::Vector{Float64})

    I, K, Kc, M, N, d, D, B, Ps, sigma2s, assignment = static_params
    f, t = utility_params

    # Preallocate aggregated interfering powers
    aggregated = fill(0., 2)

    # The BSs are grouped based on their position in the graph. If they are put
    # in a cluster, they are 'clustered', otherwise they are 'unclustered'.
    all_BSs = IntSet(1:I)
    N_clustered = length(node.a); N_unclustered = I - N_clustered
    clustered_BSs = IntSet(1:N_clustered); unclustered_BSs = IntSet(N_clustered+1:I)

    # For looping over clusters, we create a pseudo partition, where the
    # unclustered BSs belong to singleton blocks.
    pseudo_partition_a = Array(Int, I)
    @simd for i = 1:N_clustered
        @inbounds pseudo_partition_a[i] = node.a[i]
    end
    m = 1 + maximum(node.a)
    @simd for i = (N_clustered+1):I
        @inbounds pseudo_partition_a[i] = m + (i - N_clustered - 1)
    end
    pseudo_partition = Partition(pseudo_partition_a, skip_check=true)

    # The number of slots available will be important in the bounds.
    # We also store the BSs that are in the 'full' clusters,
    # i.e. clusters that cannot accept any more BSs.
    BSs_in_full_clusters = IntSet()
    for block in pseudo_partition.blocks
        N_available_ = D - length(block.elements)
        if N_available_ <= 0
            union!(BSs_in_full_clusters, block.elements)
        end
    end
    BSs_in_nonfull_clusters = setdiff(all_BSs, BSs_in_full_clusters)
    N_BSs_in_nonfull_clusters = length(BSs_in_nonfull_clusters)

    # Prelog, rate, and throughput bounds
    node_is_leaf = is_leaf(node, I)
    throughput_bounds = zeros(Float64, K, d)
    for block in pseudo_partition.blocks
        cluster_size = length(block.elements)
        N_available_slots = D - cluster_size

        # Enforce hard cluster size constraint
        if cluster_size > D
            for i in block.elements; for k in served_MS_ids(i, assignment)
                @simd for n = 1:d
                    @inbounds throughput_bounds[k,n] = -Inf
                end
            end; end
            break
        end

        # Stuff needed for rate bounds
        outside_all_BSs = setdiff(all_BSs, block.elements)
        outside_clustered_BSs = setdiff(clustered_BSs, block.elements)
        outside_BSs_in_nonfull_clusters = setdiff(BSs_in_nonfull_clusters, block.elements)
        N_outside_BSs_in_nonfull_clusters = length(outside_BSs_in_nonfull_clusters)

        # Get appropriate bounds for all MSs in this cluster.
        for i in block.elements; for k in served_MS_ids(i, assignment)
            @inbounds aggregated[1] = 0.; @inbounds aggregated[2] = 0.

            # SINR bound (solution to supermodular optimization problem)
            if node_is_leaf
                # All BSs outside my cluster contribute irreducible interference.
                sum_irreducible_interference_power!(aggregated, k, outside_all_BSs, interfering_powers)
            else
                # The bounds depends on if this BS is clustered or not.
                if i <= N_clustered
                    # Clustered BSs outside my cluster contribute irreducible interference.
                    sum_irreducible_interference_power!(aggregated, k, outside_clustered_BSs, interfering_powers)

                    # The bound is now due to picking the N_available_slots strongest interferers (that are not clustered)
                    # and assuming that this interference is reducible. This is a bound since we are not ensuring the
                    # disjointness of the clusters here.
                    sum_reducible_interference_power!(aggregated, k, unclustered_BSs, N_unclustered - N_available_slots,
                                                      interfering_powers, sub(scratch, 1:N_unclustered))
                else
                    # I cannot join any BS that belong to a full cluster, so
                    # those BSs contribute irreducible interference.
                    sum_irreducible_interference_power!(aggregated, k, BSs_in_full_clusters, interfering_powers)

                    # We now pick the N_available_slots strongest interferers (which do not belong to full clusters),
                    # and assume that this interference is reducible.
                    sum_reducible_interference_power!(aggregated, k, outside_BSs_in_nonfull_clusters, N_outside_BSs_in_nonfull_clusters - N_available_slots,
                                                      interfering_powers, sub(scratch, 1:N_outside_BSs_in_nonfull_clusters))
                end
            end
            @inbounds SNR = desired_powers[k]/sigma2s[k]
            @inbounds SINR_bound = desired_powers[k]/(sigma2s[k] + sum(aggregated))

            # Cluster size bound
            if node_is_leaf
                cluster_size_bound = cluster_size
            else
                # Use the globally optimal B, if there is one. Otherwise, find
                # it, given the current SNR and SINR bound.
                if B != 0
                    B_local = B
                else
                    @simd for b in 1:I
                        @inbounds scratch[b] = t(b, SNR, SINR_bound)::Float64
                    end
                    _, B_local = findmax(scratch)
                end
                if cluster_size < B_local
                    # There is still space, so add more BSs to this cluster. We never want more than B
                    # BSs in our clusters, because that is when the prelog starts going down again.
                    if i <= N_clustered
                        cluster_size_bound = min(cluster_size + N_unclustered, B_local)
                    else
                        cluster_size_bound = min(cluster_size + N_BSs_in_nonfull_clusters, B_local)
                    end
                else
                    # We (potentially) already have too many BSs in this cluster,
                    # thus our prelog can only go down by adding more. Bound the
                    # prelog by our current value.
                    cluster_size_bound = cluster_size
                end
            end

            # Throughput bound (though unimodality and increasingness)
            throughput_bound = t(cluster_size_bound, SNR, SINR_bound)::Float64
            @simd for n = 1:d
                @inbounds throughput_bounds[k,n] = throughput_bound
            end
        end; end
    end

    # Final system bound
    node.upper_bound = f(throughput_bounds)

    # Lumberjack.debug("Bounding.", Dict(:node => node))

    return throughput_bounds
end

function sum_irreducible_interference_power!(aggregated, k, from_BSs, interfering_powers)
    for j in from_BSs
        @inbounds aggregated[1] += interfering_powers[j,k]
    end
end

function sum_reducible_interference_power!(aggregated, k, from_BSs, N_available, interfering_powers, scratch)
    idx = 1
    for j in from_BSs
        @inbounds scratch[idx] = interfering_powers[j,k]
        idx += 1
    end
    sort!(scratch)
    @simd for idx in 1:N_available
        @inbounds aggregated[2] += scratch[idx]
    end
end
