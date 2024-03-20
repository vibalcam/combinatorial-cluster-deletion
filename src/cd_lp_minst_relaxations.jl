"""
Combinatorial solver for STC-LP relying in PushRelabelMaxFlow
"""

using SparseArrays
using LinearAlgebra

include("PushRelabelMaxflow.jl")


function CDSTC_MINST(A, Elist, use_min_source=false)
    n = size(A,1)
    # adjacency list
    Neighbs = ConstructAdj(A,n)

    # create matrix B for min-st, two nodes (++, --) for each edge
    n_B = round(Int64, sum(A) + 2)
    # B = sparse(Float64[], Float64[], Float64[], n_B, n_B)
    # B = zeros(n_B, n_B)
    B_i = Vector{Int64}()
    B_j = Vector{Int64}()
    B_v = Vector{Float64}()
    function add_edge_B(i,j,v)
        push!(B_i, i)
        push!(B_j, j)
        push!(B_v, v)
    end

    # given edge (i,j) with j > i, ij-- given by idx, ij++ idx+1
    index = Dict()
    idx = 0

    # loop over all vertices to build s-t graph
    for i = 1:n
        # get neighbors of i
        N = @views Neighbs[i]

        # loop over neighbors of i (edge ij)
        for jj = 1: length(N)
            j = @views N[jj]

            # unseen edge
            if j > i && A[i, j] > 0
                idx += 2
                index[(i, j)] = idx

                # connect s to node corresponding to edge ij-- with 1
                # B[1,idx] = 1
                add_edge_B(1, idx, 1)
                # connects node corresponding to ij++ to t with 1
                # B[idx+1, n_B] = 1
                add_edge_B(idx+1, n_B, 1)
            end

            # loop over unseen neigbors of i (edge ik)
            for kk = 1:jj-1
                # k < j
                k = @views N[kk]

                # triangle (edge jk exists)
                if @views A[k,j] > 0
                    continue
                end

                # open wedge centered at i with inf
                a1,a2 = minmax(i, j)
                b1,b2 = minmax(i, k)
                idx1 = @views index[(a1,a2)]
                idx2 = @views index[(b1,b2)]
                # connect ij-- to ik++
                # B[idx1, idx2+1] = Inf
                add_edge_B(idx1, idx2+1, Inf)
                # connect ik-- to ij++ with inf
                # B[idx2, idx1+1] = Inf
                add_edge_B(idx2, idx1+1, Inf)
            end
        end
    end

    B = sparse(B_i, B_j, B_v, n_B, n_B)
    s = 1
    t = n_B
    tic = time()
    F = maxflow(B,s,t)
    if use_min_source
        S = source_nodes_min(F)
    else
        S = source_nodes(F)
    end
    # get array with 1 if in S, 0 otherwise
    x_h = @views [i in S ? 1 : 0 for i in 1:n_B]
    time_st = time()-tic
    # get array of -- and of ++ without s and t
    x_minus = @views x_h[2:2:end-1]
    x_plus = @views x_h[3:2:end-1]

    # calculate back the variables x = 1/2 + (x++ - x--)/2
    x = @views 1/2 .+ (x_plus - x_minus) ./ 2

    # ensure that x is returned in the same order as Elist
    idx = @views [index[(Elist[i,1], Elist[i,2])] for i in 1:size(Elist, 1)] .÷ 2
    sol = @views x[idx]

    return sum(sol), sol, time_st
    # sum(sol): objective value
    # sol: solution
    # time_st: runtime for min-st
end
