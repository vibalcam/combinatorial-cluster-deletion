"""
Same methods as in cd_det_mfp.jl but without the initialization

This is useful for a more modular approach, where the lower bound
is computed for all MFP methods, then the graph is created and 
passed to each method.
"""

function det_pivot_degree(Anew)

    # Get the number of rows in matrix A
    n = size(Anew,1)

    # RUN THE DETERMINISTIC PIVOT ROUNDING ALGORITHM
    ################################################
    counter = 0
    V = trues(n) # Unclustered nodes
    clus = zeros(Int, n)  # Cluster assignment vector
    nc = 1  # Cluster counter

    # compute degrees of nodes
    INVALID_NODE = 0
    degrees = sum(Anew, dims=1)[:]
    degrees = convert(Vector{Int}, degrees)
    curr_bin = Int(maximum(degrees))
    # initalize bins for each degree
    bins = [Int[] for _ in 1:curr_bin]
    node_to_bin = zeros(Int, n)

    # # function to set index as singleton
    # function set_singleton(idx)
    #     # node_to_bin[idx] = INVALID_NODE
    #     # cluster
    #     nc += 1
    #     clus[idx] = nc
    #     V[idx] = false
    # end

    # assign nodes to bins
    for (idx,d) in enumerate(degrees)
        # if degree 0, cluster singleton and continue
        if d == 0
            clus[idx] = nc
            V[idx] = false
            nc += 1

            # node_to_bin[idx] = INVALID_NODE
            # nc += 1
            # clus[idx] = nc
            # V[idx] = false
            continue
        end

        push!(bins[d], idx)
        node_to_bin[idx] = length(bins[d])
        # if isdefined(bins, d)
        #     push!(bins[d], idx)
        #     node_to_bin[idx] = length(bins)
        # else
        #     bins[d] = [idx]
        #     node_to_bin[idx] = 1
        # end
    end

    # transform to adjacency list
    Anew = ConstructAdj(Anew, n)
    
    idx_bin = 1
    p = INVALID_NODE
    while any(V)
        while true
            # check if we need to move to next bin
            while idx_bin > length(bins[curr_bin])
                curr_bin -= 1
                idx_bin = 1
            end
            # choose highest degree node
            p = bins[curr_bin][idx_bin]
            idx_bin += 1
            # check if node is valid
            if p != INVALID_NODE && V[p]
                break
            end
        end

        # form cluster with p and its neighbors
        p_neigh = Anew[p]
        idx = V[p_neigh]
        cluster = p_neigh[idx]
        # cluster = [k for k in Anew[p] if V[k]]
        # cluster = collect(rows[nzrange(Anew, p)])
        push!(cluster, p)

        # @assert(all(V[cluster]))
        clus[cluster] .= nc
        nc += 1

        # counter += length(cluster)
        counter = sum(V)
        print("\r$(100 - counter*100/n)")

        # remove cluster from V and G
        V[cluster] .= false

        # update degrees
        for k in cluster
            bins[degrees[k]][node_to_bin[k]] = INVALID_NODE

            # ignore pivot
            if k == p
                continue
            end
            # check neighbors of nodes
            for j in Anew[k]
                # ignore clustered nodes
                if !V[j]
                    continue
                end
                # update degree
                deg_j = degrees[j]
                degrees[j] = deg_j - 1
                # update old bin
                bin_j = node_to_bin[j]
                bins[deg_j][bin_j] = INVALID_NODE
                # check if new degree is 0
                if deg_j == 1
                    # set_singleton(j)

                    # cluster singleton
                    clus[j] = nc
                    V[j] = false
                    nc += 1
                    continue
                end
                # move to new bin
                push!(bins[deg_j-1], j)
                node_to_bin[j] = length(bins[deg_j-1])
            end
        end
    end
    
    return clus
end

function det_pivot_ratio(Anew)
    
    n = size(Anew,1)

    # transform to adjacency list
    Anew = ConstructAdj(Anew, n)

    # RUN THE DETERMINISTIC PIVOT ROUNDING ALGORITHM
    ################################################
    counter = 0
    # V = ones(Bool, n) # Unclustered nodes
    V = trues(n) # Unclustered nodes
    clus = zeros(Int, n)  # Cluster assignment vector
    nc = 0  # Cluster counter

    numerator = zeros(Int, n)
    denominator = zeros(Int, n)
    scores = zeros(Float64, n)

    # dictionary of wedges where the node is the center
    wedges_tmp = [Dict{Int, Set{Int}}() for _ in 1:n]
    # dictionary of centers of the wedges a node i belongs to (not including itself)
    centers_tmp = [Set{Int}() for _ in 1:n]

    # get initial scores
    # rows = rowvals(Anew)
    n_wedges = 0
    # numerator[V] .= 0
    # denominator[V] .= 0

    # Threads.@threads for k in 1:n
    for k in 1:n
        # find neighbors of k
        neighbors_k = Anew[k]
        for idx_neigh_i in eachindex(neighbors_k)
            i = neighbors_k[idx_neigh_i]
            # no need to check values since we are removing
            # zeros from Anew, thus it must be a nonzero edge

            for idx_neigh_j in idx_neigh_i+1:length(neighbors_k)
                j = neighbors_k[idx_neigh_j]
                if i != j
                    # open wedge centered at k found
                    n_wedges += 1
                    denominator[k] += 1
                    numerator[i] += 1
                    numerator[j] += 1

                    # add wedges to each node (i, center, j)
                    # s1, s2 = i < j ? (i,j) : (j,i)
                    push!(centers_tmp[i], k)
                    push!(centers_tmp[j], k)
                    # not including itself as center
                    # push!(centers[k], k)

                    # push!(wedges[k], (s1, k, s2))
                    # push!(wedges[k], (s1, s2))
                    push!(get!(wedges_tmp[k], i, Set{Int}()), j)
                    push!(get!(wedges_tmp[k], j, Set{Int}()), i)
                end
            end
        end
    end

    # change sets to arrays
    wedges = Vector{Dict{Int, Vector{Int}}}(undef, length(wedges_tmp))
    centers = Vector{Vector{Int}}(undef, length(centers_tmp))
    for (i, dict_k) in enumerate(wedges_tmp)
        d = Dict{Int, Vector{Int}}()
        for (j, set_k) in dict_k
            d[j] = collect(set_k)
        end
        wedges[i] = d
    end
    for (i, set_k) in enumerate(centers_tmp)
        centers[i] = collect(set_k)
    end
    wedges_tmp = nothing
    centers_tmp = nothing
    
    while any(V)
        nc += 1

        if n_wedges == 0
            # get the first nonzero value
            p = findfirst(V)
        else
            # scores[V] = numerator[V] ./ denominator[V]
            # scores[denominator .== 0] .= Inf
            # p = argmin(scores)
            p = argmin(ifelse.(denominator .== 0, Inf, numerator ./ denominator))
        end

        # form cluster with p and its neighbors
        p_neigh = Anew[p]
        idx = V[p_neigh]
        cluster = p_neigh[idx]
        # cluster = [k for k in Anew[p] if V[k]]
        # cluster = collect(rows[nzrange(Anew, p)])
        push!(cluster, p)

        # @assert(all(V[cluster]))
        clus[cluster] .= nc

        # counter += length(cluster)
        counter = sum(V)
        print("\r$(100 - counter*100/n)")

        # compute changes to scores for next iteration
        if n_wedges > 0
            # remove wedges from each node in cluster
            # Threads.@threads for i in cluster
            for i in cluster
                # remove all wedges it is the center of
                for (s1, set_s2) in wedges[i]
                    for s2 in set_s2
                        # to ensure we do not double count wedges
                        if s1 > s2
                            continue
                        end

                        # update scores
                        n_wedges -= 1
                        denominator[i] -= 1
                        numerator[s1] -= 1
                        numerator[s2] -= 1
                    end
                end
                # cleanup vector
                empty!(wedges[i])
                # check scores for correctness
                @assert(denominator[i] == 0)

                # remove wedges it is not the center of
                for k in centers[i]
                    # if k is in cluster, it will be removed later
                    if k in cluster
                        continue
                    end
                    wedges_k = wedges[k]
                    # check if that specific wedge has been removed
                    if !haskey(wedges_k, i)
                        continue
                    end
                    # otherwise, remove wedge
                    for j in wedges_k[i]
                        # update scores
                        n_wedges -= 1
                        denominator[k] -= 1
                        numerator[i] -= 1
                        numerator[j] -= 1
                        
                        # remove wedge from other node
                        deleteat!(wedges_k[j], findfirst(wedges_k[j] .== i))
                    end
                    empty!(wedges_k[i])
                end
            end
        end

        # remove cluster from V and G
        V[cluster] .= false
        # Anew[cluster, :] .= 0
        # Anew[:, cluster] .= 0
        # dropzeros!(Anew)
    end
    
    println()
    return clus
end