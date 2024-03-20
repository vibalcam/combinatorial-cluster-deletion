## Functions for merging clustering to improve the results of MFP

# Please use this function to merge clusters.
# Input:    A: just A in context
#           clus: the original solution clusters
# Output:   A feasible cluster deletion solution which cannot be furter merged.
# In each loop the function tries to find two clusters that can be merged. If found, 
# merge them and start a new loop. Otherwise, return the merged clusters.
# Time complexity:  O(|clus|(|E| + |clus|^2)). There are at most |clus| loops. In each loop, once
# a non-edge is found between two clusters, it will stop checking the edges between 
# them and starts to check the edges between other clusters.
function mfpmerge(A, clus)
    new_clus = copy(clus)
    # println(clus)
    while(true)
        new_clus, cls_update = try_merge_cd_clus(A, new_clus)
        if cls_update == 0
            break
        end
    #    println(new_clus)
    end
    return new_clus
end

# Try to merge clusters if possible
function try_merge_cd_clus(A, clus)
    C = Vector{Vector}()
    for i = 1:(findmax(clus)[1])
        push!(C, Vector{Int64}())
    end
    for i in eachindex(clus)
        push!(C[clus[i]], i)
    end
    for i = 1:length(C)
        for j = i+1:length(C)
            flag = 1
            for u in C[i]       # for each vertex u in cluster C[i]
                for v in C[j]   # for each vertex v in cluster C[j]
                    if A[u, v] <= 0
                        flag = 0
                        break
                    end
                end
                if flag == 0
                    break
                end
            end
            if flag == 1
                return merge_clus(clus, C, i, j)
            end
        end
    end
    return clus, 0
end

# Combine clusters
function merge_clus(clus, C, x, y)
    i = min(x, y)
    j = max(x, y)
    for v in C[j]
        clus[v] = i
    end
    for ii in eachindex(clus)
        if clus[ii] > j
            clus[ii] -= 1
        end
    end
    return clus, 1
end