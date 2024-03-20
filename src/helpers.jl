using Random

"""
Just a linearization of the node pairs. 
"""
function cc_elist(n)
    Elist = zeros(Int64,round(Int64,n*(n-1)/2),2)
    nvars = 0
    for j=1:n
        for i=1:j-1
            nvars += 1
            Elist[nvars,1] = i
            Elist[nvars,2] = j
        end
    end
    return Elist
end


"""
For graph with adjacency matrix A, extract triplet indices
for open wedges and triangles.

For each entry in W, the first index is the center of the wedge
each triangle is listed 3 times, once for each "center"
"""
function get_wedges_triangles(A)
    n = size(A,1)
    Neighbs = ConstructAdj(A,n)
    T = Vector{Tuple{Int,Int,Int}}()
    W = Vector{Tuple{Int,Int,Int}}()
    for i = 1:n
        N = Neighbs[i]
        for jj = 1:length(N)
            j = N[jj]
            for kk = jj+1:length(N)
                k = N[kk]
                if A[j,k] == 1
                    push!(T,(i,j,k))
                else
                    push!(W,(i,j,k))
                end
            end
        end
    end
    return T, W
end

# From the adjacency matrix, build an adjacency list for the graph
function ConstructAdj(C::SparseMatrixCSC,n::Int64)
    rp = C.rowval
    ci = C.colptr
    Neighbs = Vector{Vector{Int64}}()
    for i = 1:n
        # chop up the rp vector and put it in Neighbs
        push!(Neighbs,rp[ci[i]:ci[i+1]-1])
    end
    return Neighbs

end

function count_wedges_triangles(A)
    n = size(A,1)
    Neighbs = ConstructAdj(A,n)
    T = 0
    W = 0
    for i = 1:n
        N = Neighbs[i]
        # println("$i")
        for jj = 1:length(N)
            j = N[jj]
            for kk = jj+1:length(N)
                k = N[kk]
                if A[j,k] == 1
                    T += 1
                else
                    W += 1
                end
            end
        end
    end
    return T, W
end




"""
Confirm that a clustering of a graph is indeed a disjoint
union of cliques
"""
function isCDfeasible(A,c)
    for j = 1:maximum(c)
        S = findall(x->x==j,c)
        if length(S) == 1
            continue
        end
        for ii = 1:length(S)
            Ind = S[ii]
            for jj = ii+1:length(S)
                Jnd = S[jj]
                if A[Ind,Jnd] == 0
                    return false
                end
            end
        end
    end
    return true
end


"""
Given a clustering, compute the
unweighted correlation clustering objective.
This may be a cluster deletion or a
cluster editing solution.

Naive, slow code.
"""
function cc_cd_obj(A,c)
    n = size(A,1)
    mistakes = 0
    for i = 1:n
        for j = i+1:n
            if A[i,j] == 0 && c[i] == c[j]
                mistakes += 1
                
                # for cluster deletion this should not happen
                @assert(false) 
            end
            if A[i,j] == 1 && c[i] != c[j]
                mistakes += 1
            end
        end
    end
    return mistakes
end

"""
Faster checking of objective
"""
function check_cc_obj_fastish(A::SparseArrays.SparseMatrixCSC,Elist,c,m)
    n = length(c)
    pos_mis = 0
    for k = 1:size(Elist,1)
        if c[Elist[k,1]] != c[Elist[k,2]]
            pos_mis += 1
        end
    end

    clus_num = maximum(c)
    clus_sizes = zeros(Int64,clus_num)
    for i = 1:n
        clus_sizes[c[i]] += 1
    end

    neg_mis = pos_mis - m
    for k = 1:maximum(clus_num)
        S = clus_sizes[k]
        neg_mis += S*(S-1)/2
    end
    return round(Int64,neg_mis) + pos_mis
end

"""
Fast version of the pivot algorithm: a random permutation decides the ordering of pivots.
"""
function permutation_pivot(A)
    n = size(A,1)
    p = randperm(n)         # specify pivots in advance--equivalent to uniform random selection at each step
    c = zeros(Int64,n)      # current cluster
    Neighbs = ConstructAdj(A,n)
    clusnum = 1
    for i = 1:n
        K = p[i]        # current pivot candidate
        if c[K] > 0
            # already clustered, move to next pivot
            continue
        end
        nbs = Neighbs[K]
        c[K] = clusnum
        for j in nbs
            # Rather than deleting nodes everywhere at each step, just check them and
            # ignore if they are already clustered
            if c[j] == 0
                c[j] = clusnum
            end
        end
        clusnum += 1
    end
    return c
end


"""
Permutation pivot, in the case where the graph is given as an adjacency list
"""
function permutation_pivot(Neighbs::Vector{Vector{Int64}},n)
    # specify pivots in advance--equivalent to uniform random selection at each step
    p = randperm(n)
    c = zeros(Int64,n)      # current cluster
    clusnum = 1
    for i = 1:n
        K = p[i]        # current pivot candidate
        if c[K] > 0
            # already clustered, move to next pivot
            continue
        end
        nbs = Neighbs[K]
        c[K] = clusnum
        for j in nbs
            # Rather than deleting nodes everywhere at each step, just check them and
            # ignore if they are already clustered
            if c[j] == 0
                c[j] = clusnum
            end
        end
        clusnum += 1
    end
    return c
end


"""
find_violations
Given a candidate distance matrix D, iterate through all 3-tuples of nodes
and return the tuples where triangle inequality constraints have been violated.

Output is stored in vector 'violations'

Note that we want D to be lower triangular here, though if it is symmetric
this will work fine as well. We just need to make sure the distance information
in D is not just stored in the upper triangular portion
"""
function find_violations!(D::Matrix{Float64}, violations::Vector{Tuple{Int,Int,Int}})
  n = size(D,1)

  # We only need this satisfied to within a given tolerance, since the
  # optimization software will only solve it to within a certain tolerance
  # anyways. This can be tweaked if necessary.
  epsi = 1e-8
  @inbounds for i = 1:n-2
       for j = i+1:n-1
          a = D[j,i]
           for k = j+1:n
              b = D[k,i]
              c = D[k,j]
        if a - b > epsi && a - c > epsi && a-b-c > epsi
            push!(violations, (i,j,k))
                # @constraint(m, x[i,j] - x[i,k] - x[j,k] <= 0)
        end
        if b - a > epsi && b - c > epsi && b-a-c > epsi
            push!(violations, (i,k,j))
            # @constraint(m, x[i,k] - x[i,j] - x[j,k] <= 0)
        end

        if c - a > epsi && c-b>epsi && c-a-b > epsi
            push!(violations, (j,k,i))
            # @constraint(m, x[j,k] - x[i,k] - x[i,j] <= 0)
        end
      end
    end
  end
end

