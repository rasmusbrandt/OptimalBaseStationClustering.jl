function precompute_powers(channel, assignment, static_params)
    I, K, Kc, M, N, d, Ps, sigma2s = static_params
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

function avg_cluster_size(a)
    num_clusters = 1 + maximum(a)

    num_members = zeros(Int, num_clusters)
    for i = 1:num_clusters
        num_members[i] = length(find(a .== (i-1)))
    end

    return mean(num_members)
end
