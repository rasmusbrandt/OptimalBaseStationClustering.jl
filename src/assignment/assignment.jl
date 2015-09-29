function throughputs(partition, channel, static_params, scenario_params)
    I, K, Kc, M, N, d, Ps, sigma2s, assignment = static_params
    f, t, D, B = scenario_params

    desired_powers, interfering_powers = compute_powers(channel, static_params)

    throughputs = zeros(Float64, K, d)
    for block in partition.blocks
        outside_all_BSs = setdiff(IntSet(1:I), block.elements)
        cluster_size = length(block.elements)
        for i in block.elements; for k in served_MS_ids(i, assignment)
            irreducible_interference_power = Float64(0)
            for j in outside_all_BSs
                irreducible_interference_power += interfering_powers[j,k]
            end
            SNR = desired_powers[k]/sigma2s[k]
            rho = desired_powers[k]/(sigma2s[k] + irreducible_interference_power)
            throughputs[k,:] = t(cluster_size, SNR, rho)
        end; end
    end
    throughputs
end

function compute_powers(channel, static_params)
    I, K, Kc, M, N, d, Ps, sigma2s, assignment = static_params
    desired_powers = Array(Float64, K)
    interfering_powers = Array(Float64, I, K)
    for i = 1:I; for k in served_MS_ids(i, assignment)
        desired_powers[k] = channel.large_scale_fading_factor[k,i]*channel.large_scale_fading_factor[k,i]*(Ps[i]/(Kc*d))
        for j = 1:I
            interfering_powers[j,k] = channel.large_scale_fading_factor[k,j]*channel.large_scale_fading_factor[k,j]*Ps[j]
        end
    end; end
    desired_powers, interfering_powers
end

function average_cluster_size(a)
    num_clusters = 1 + maximum(a)

    num_members = zeros(Int, num_clusters)
    for i = 1:num_clusters
        num_members[i] = length(find(a .== (i-1)))
    end

    mean(num_members)
end
