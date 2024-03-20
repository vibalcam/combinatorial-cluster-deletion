## This file originally written and shared by Nate Veldt under the following license:
# MIT License

# Copyright (c) 2022 Nate Veldt

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
##################################################################################


# This is even faster code for many of the functions that rely on running pivot multiple times
# to either the original graph or some type of derived graph constructed based on
# a certain type of lower bound.

"""
Faster version of the pivot algorithm
"""
function permutation_pivot_faster(A)
    n = size(A,1)
    p = randperm(n)         # specify pivots in advance--equivalent to uniform random selection at each step
    c = zeros(Int64,n)      # current cluster
    Neighbs = ConstructAdj(A,n)
    clusnum = 1

    # Sbin is used to more quickly compute the correlation cluster objective
    # sum of binomials: sum_{clusters S} binomial(S,2)
    Sbin = 0

    for i = 1:n
        K = p[i]        # current pivot candidate
        if c[K] > 0
            # already clustered, move to next pivot
            continue
        end
        nbs = Neighbs[K]
        c[K] = clusnum
        S = 1           # cluster size

        for j in nbs
            # Rather than deleting nodes everywhere at each step, just check them and
            # ignore if they are already clustered
            if c[j] == 0
                c[j] = clusnum
                S += 1
            end
        end
        Sbin += S*(S-1)/2
        clusnum += 1
    end
    return c, round(Int64,Sbin)
end


"""
Faster checking of objective.

Elist here is the actual edge list of a graph, don't double count.

It's not usually all pairs (unless this is a complete graph)
"""
function check_cc_obj_faster(Elist,c,m,Sbin)
    n = length(c)
    pos_mis = 0
    for k = 1:size(Elist,1)
        if c[Elist[k,1]] != c[Elist[k,2]]
            pos_mis += 1
        end
    end

    # clus_num = maximum(c)
    # clus_sizes = zeros(Int64,clus_num)
    # for i = 1:n
    #     clus_sizes[c[i]] += 1
    # end

    neg_mis = pos_mis - m + Sbin
    return round(Int64,neg_mis) + pos_mis
end

function many_pivot_faster(A,Elist,pivtimes = 10)
    m = round(Int64,sum(A)/2)
    clus, Sbin = permutation_pivot_faster(A)
    cc_obj = check_cc_obj_faster(Elist,clus,m,Sbin)
    # tst = check_cc_obj(A,clus)
    # @assert(cc_obj == tst)

    # Pivot is fast, so we can run it multiple times
    for jj = 1:pivtimes
        clusnew, Sbin = permutation_pivot_faster(A)
        objnew = check_cc_obj_faster(Elist,clusnew,m,Sbin)
        # tst = check_cc_obj(A,clusnew)
        # @assert(objnew == tst)
        if objnew < cc_obj
            cc_obj = objnew
            clus = clusnew
        end
    end
    return cc_obj, clus
end


"""
Rounds a feasible STC+ solution into a solution
    for a cluster editing problem.

Input:
    A = adjacency matrix

    Welist = node pairs to be flipped

Output:
    Approximation for CC.

    Step 1: Construct derived graph which is A with some flips
    Step 2: Apply pivot to the new graph
    Step 3: Compute CC objective score
    Step 4: Check how far this is from the lower bound

    pivtimes is the number of times to run the pivot
    method on the derived graph.

    Also returns the average time for a pivot step and the average
    quality solution
"""
function STCplus_to_CC_round_rand_faster(A,Eadd,Edel,pivtimes = 5)

    # Elist is list of edges for graph A
    n = size(A,1)
    tic = time()
    Anew = flip_graph(A,Eadd,Edel)
    setuptime = time()-tic
    # @show setuptime
    Is, Js = findnz(triu(A))
    ElistA = [Is Js]
    m = sum(A)/2

    NeighbsNew = ConstructAdj(Anew,n)

    # Apply pivot on this graph multiple times,
    # returning the best output

    # clus = permutation_pivot(Anew) # Apply pivot on one graph
    # cc_obj = check_cc_obj(A,clus)  # Check output on the other graph

    clus = zeros(n)
    Sbin = permutation_pivot_fastest!(Anew,NeighbsNew,clus)
    cc_obj = check_cc_obj_faster(ElistA,clus,m,Sbin)
    # cc_obj2 = check_cc_obj_fastish(A,Elist,clus,m)
    # objnew = check_cc_obj(A,clus)
    # @show objnew, cc_obj, cc_obj2
    # @assert(cc_obj == objnew)
    clusnew = zeros(n)
    objs = cc_obj
    tic = time()
    for jj = 1:pivtimes
        # Pivot is fast, so we can run it multiple times
        # clusnew, Sbin = permutation_pivot_faster(Anew)
        # clusnew, Sbin = permutation_pivot_fastest(Anew,NeighbsNew)
        Sbin = permutation_pivot_fastest!(Anew,NeighbsNew,clusnew)
        objnew = check_cc_obj_faster(ElistA,clusnew,m,Sbin)
        # objnew2 = check_cc_obj(A,clusnew) # Check output on original graph
        # @show objnew, objnew2
        # @assert(objnew2 == objnew)

        objs += objnew
        if objnew < cc_obj
            cc_obj = objnew
            clus = clusnew
        end
    end
    totaltime = time()-tic
    avg_time = totaltime/pivtimes + setuptime
    avg_obj = objs/pivtimes
    return round(Int64,cc_obj), clus, avg_obj, avg_time
end


"""
Fastest version of the pivot algorithm
"""
function permutation_pivot_fastest(A,Neighbs)
    n = size(A,1)
    p = randperm(n)         # specify pivots in advance--equivalent to uniform random selection at each step
    c = zeros(Int64,n)      # current cluster
    clusnum = 1

    # Sbin is used to more quickly compute the correlation cluster objective
    # sum of binomials: sum_{clusters S} binomial(S,2)
    Sbin = 0

    for i = 1:n
        K = p[i]        # current pivot candidate
        if c[K] > 0
            # already clustered, move to next pivot
            continue
        end
        nbs = Neighbs[K]
        c[K] = clusnum
        S = 1           # cluster size

        for j in nbs
            # Rather than deleting nodes everywhere at each step, just check them and
            # ignore if they are already clustered
            if c[j] == 0
                c[j] = clusnum
                S += 1
            end
        end
        Sbin += S*(S-1)/2
        clusnum += 1
    end
    return c, round(Int64,Sbin)
end


"""
Fastest version of the pivot algorithm
"""
function permutation_pivot_fastest!(A,Neighbs,c)
    n = size(A,1)
    p = randperm(n)         # specify pivots in advance--equivalent to uniform random selection at each step
    # c = zeros(Int64,n)      # current cluster
    clusnum = 1

    for i = 1:n
        c[i] = 0
    end
    # Sbin is used to more quickly compute the correlation cluster objective
    # sum of binomials: sum_{clusters S} binomial(S,2)
    Sbin = 0

    for i = 1:n
        K = p[i]        # current pivot candidate
        if c[K] > 0
            # already clustered, move to next pivot
            continue
        end
        nbs = Neighbs[K]
        c[K] = clusnum
        S = 1           # cluster size

        for j in nbs
            # Rather than deleting nodes everywhere at each step, just check them and
            # ignore if they are already clustered
            if c[j] == 0
                c[j] = clusnum
                S += 1
            end
        end
        Sbin += S*(S-1)/2
        clusnum += 1
    end
    return round(Int64,Sbin)
end


"""
Rounds a lower bound for STC+ into a feasible solution
    for a cluster editing problem.

This works for:
    * rounding the CC LP relaxation (2-approx)
    * rounding the STC+ LP relaxation (4-approx)
    * rounding a feasible solution to the STC+ ILP (6-approx)

Input:
    A = adjacency matrix

    Elist = linear ordering of node-pairs

    soln = feasible solution to the CC LP relaxation
            or the STC+ LP relaxation
            or the STC+ ILP

Output:
    Approximation for CC.

    Step 1: Construct derived graph where edges are pairs with x[e] < 1/2
    Step 2: Apply pivot to the new graph
    Step 3: Compute CC objective score
    Step 4: Check how far this is from the lower bound

    pivtimes is the number of times to run the pivot 
    method on the derived graph.
""" 
function STCplus_to_CC_round_LP(A,Elist,soln,pivtimes = 1)
    n = size(A,1)

    # Elist is a linearization of all node pairs in A,
    # not just the edges in A.

    # Find all the node pairs Elist[k] that have soln[k] < 1/2
    tic = time()
    tol = 1e-8
    keep = findall(x->x<0.5-tol,soln)

    # These all become edges in a new derived graph Anew
    Inew = Elist[keep,1]
    Jnew = Elist[keep,2]
    Anew = sparse(Inew,Jnew,ones(length(Jnew)),n,n)
    Anew = Anew + Anew'

    Is, Js = findnz(triu(A))
    ElistA = [Is Js]
    m = sum(A)/2

    NeighbsNew = ConstructAdj(Anew,n)
    setuptime = time()-tic
    # Apply pivot on this graph multiple times,
    # returning the best output

    # clus = permutation_pivot(Anew) # Apply pivot on one graph
    # cc_obj = check_cc_obj(A,clus)  # Check output on the other graph
    clus, Sbin = permutation_pivot_fastest(Anew,NeighbsNew)
    cc_obj = check_cc_obj_faster(ElistA,clus,m,Sbin)
    # cc_obj2 = check_cc_obj_fastish(A,Elist,clus,m)
    # objnew = check_cc_obj(A,clus)
    # @show objnew, cc_obj, cc_obj2
    # @assert(cc_obj == objnew)
    objs = cc_obj
    tic = time()
    for jj = 1:pivtimes
        # Pivot is fast, so we can run it multiple times
        # clusnew, Sbin = permutation_pivot_faster(Anew)
        # clusnew, Sbin = permutation_pivot_fastest(Anew,NeighbsNew)
        clusnew,Sbin = permutation_pivot_fastest(Anew,NeighbsNew)
        objnew = check_cc_obj_faster(ElistA,clusnew,m,Sbin)
        # objnew2 = check_cc_obj(A,clusnew) # Check output on original graph
        # @show objnew, objnew2
        # @assert(objnew2 == objnew)

        objs += objnew
        if objnew < cc_obj
            cc_obj = objnew
            clus = clusnew
        end
    end
    totaltime = time()-tic
    avg_time = totaltime/pivtimes + setuptime
    avg_obj = objs/pivtimes
    return round(Int64,cc_obj), clus, avg_obj, avg_time
end