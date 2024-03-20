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
#####################################################
#
# Several updates made for algorithms accompanying paper: 
# "Combinatorial Approximations for Cluster Deletion: Simpler, Faster, and Better"


# Linear programming relaxations for cluster deletions
using Gurobi
using SparseArrays
using LinearAlgebra
gurobi_env = Gurobi.Env()

include("helpers.jl")


"""
CDSTC_Gurobi (CDSTC = cluster deletion & strong triadic closure)

Use Gurobi WITHOUT JuMP to solve cluster deletion or 
minimum weakness strong triadic closure labeling,
either optimally, or just the canonical LP relaxation.

This is quite a bit faster than using JuMP.

Input:
    A = adjacency matrix (unsigned, unweighted)

Optional parameters:
    clusterdeletion = true:
        solve cluster deletion objective and not just stc

    LP = true:
        just solve the LP relaxation, not the exact ILP

    Timelimit: upper limit on how long to let Gurobi run

    Outputflag: whether or not to show Gurobi solver status during solve

Output:
    obj = output objective
    Elist = ordered list of edges in graph
    Evals = variable value for each edge
"""
function CDSTC_Gurobi(A; clusterdeletion = true, outputflag = false, LP = false, timelimit = Inf)
    n = size(A,1)
    Is, Js, Vs = findnz(triu(A))
    # build varmap and obj
    vmap = -ones(Int,n,n)
    obj = Float64[]
    nvars = 0
    Elist = zeros(Int64,length(Vs),2)
    for k = 1:length(Vs)
        # variable for each edge
        nvars += 1
        i = Is[k]
        j = Js[k]
        Elist[k,1] = i
        Elist[k,2] = j
        vmap[i,j] = nvars-1
        vmap[j,i] = nvars-1
        push!(obj, 1.0)
    end
    aptr = Ref{Ptr{Cvoid}}()
    # GRBnewmodel (	GRBenv	*env,
 	#  	GRBmodel	**modelP,
 	#  	const char	*Pname,
 	#  	int	numvars,
 	#  	double	*obj,
 	#  	double	*lb,
 	#  	double	*ub,
 	#  	char	*vtype,
 	#  	const char	**varnames )
    if LP
        vtypes = repeat(GRB_CONTINUOUS, nvars)
        ub = ones(nvars)
        err = GRBnewmodel(gurobi_env, aptr, "CDLP", nvars, obj, C_NULL, ub, vtypes, C_NULL)

    else
        vtypes = repeat(GRB_BINARY, nvars)
        err = GRBnewmodel(gurobi_env, aptr, "ExactCD", nvars, obj, C_NULL, C_NULL, vtypes, C_NULL)
    end

    m = aptr[]
    GRBsetdblparam(GRBgetenv(m), GRB_DBL_PAR_TIMELIMIT, timelimit)
    GRBsetintparam(GRBgetenv(m), "OutputFlag",outputflag)

    T,W = get_wedges_triangles(A)

    try
        cind2 = Int32[0,0]
        cval2 = Float64[0,0]
        for t in W
            i = t[1]    # This is always the center
            j = t[2]
            k = t[3]
            cind2[1] = vmap[min(i,j),max(i,j)] # xij
            cind2[2] = vmap[min(i,k),max(i,k)] # xik
            # x[i,j] + x[i,k] >= 1, so -xij - xjk <= -1
            cval2[1] = -1  
            cval2[2] = -1
            error = GRBaddconstr(m, 2, cind2, cval2, GRB_LESS_EQUAL, -1.0, C_NULL)
        end
    
        if clusterdeletion
            # If we are solving cluster deletion, need more triangle inequality constraints
            cind = Int32[0,0,0]
            cval = Float64[0,0,0]
            # println("Adding extra triangle constraints.")
            for t in T
                i = t[1]
                j = t[2]
                k = t[3]
                cind[1] = vmap[min(i,j),max(i,j)] # xij
                cind[2] = vmap[min(i,k),max(i,k)] # xik
                cind[3] = vmap[min(j,k),max(j,k)] # xjk
                # @constraint(m, x[i,j] + x[i,k] >= x[j,k])
                # @constraint(m, x[i,j] + x[j,k] >= x[i,k])
                # @constraint(m, x[j,k] + x[i,k] >= x[i,j])

                # xij - xik - xjk <= 0.0
                cval[1] = 1  
                cval[2] = -1
                cval[3] = -1
                error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)

                # -xij + xik - xjk <= 0.0
                cval[1] = -1  
                cval[2] = 1
                cval[3] = -1
                error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)

                # -xij - xik + xjk <= 0.0
                cval[1] = -1  
                cval[2] = -1
                cval[3] = 1
                error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
            end
        end

    GRBoptimize(m)
    stat = Ref{Int32}(0)
    GRBgetintattr(m, GRB_INT_ATTR_STATUS, stat)
    #Status codes: https://www.gurobi.com/documentation/9.5/refman/optimization_status_codes.html#sec:StatusCodes
    status = stat[]
    # println("Status = $status")
    if status == 2
        optimal = true
    else 
        optimal = false
    end
    robj = Ref{Float64}(0.0)
    GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, robj)

    obj = robj[]
    soln = zeros(nvars)
    GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, nvars, soln)

    # return D, obj
    return optimal, obj, Elist, soln
 finally
     #@show "freemodel here!"
     GRBfreemodel(m)
 end
end


"""
Given a graph A, return a maximal set of edge-disjoint open wedges.
"""
function maximal_disjoint_openwedge(A)
    n = size(A,1)
    Neighbs = ConstructAdj(A,n)
    Is = Vector{Int64}()
    Js = Vector{Int64}()
    B = copy(A)
    for i = 1:n
        N = Neighbs[i]
        for jj = 1:length(N)
            j = N[jj]
            ij1 = max(i,j)
            ij2 = min(i,j)
            if B[ij1,ij2] == 2
                continue
            end
            for kk = jj+1:length(N)
                k = N[kk]
                ik1 = max(i,k)
                ik2 = min(i,k)
                if B[ik1,ik2] == 2 || A[k,j] > 0
                    continue
                else
                    # if we reach here: (i,j) and (i,k) are edges, but (j,k) is not
                    # So (i,j,k) is an open wedge centered at i.
                    # Also, neither (i,j) nor (i,k) already appears in our open wedge set.
                    push!(Is,ij2)
                    push!(Js,ij1)
                    push!(Is,ik2)
                    push!(Js,ik1)
                    B[ij1,ij2] = 2
                    B[ik1,ik2] = 2

                    # Now exit this and start with a new j, as we have already
                    # used edge (i,j) and we don't need to check for more
                    # edges involving (i,j)
                    break
                end
            end
        end
    end
    Welist = [Is Js] # list of the edges in the wedgelist
    return Welist
end


"""
Same function, with output in a more convenient format for rounding
"""
function maximal_disjoint_openwedge_rounding(A)
    n = size(A,1)
    Neighbs = ConstructAdj(A,n)
    Is = Vector{Int64}()
    Js = Vector{Int64}()
    B = copy(A)
    for i = 1:n
        N = Neighbs[i]
        for jj = 1:length(N)
            j = N[jj]
            ij1 = max(i,j)
            ij2 = min(i,j)
            if B[ij1,ij2] == 2
                continue
            end
            for kk = jj+1:length(N)
                k = N[kk]
                ik1 = max(i,k)
                ik2 = min(i,k)
                if B[ik1,ik2] == 2 || A[k,j] > 0
                    continue
                else
                    # if we reach here: (i,j) and (i,k) are edges, but (j,k) is not
                    # So (i,j,k) is an open wedge centered at i.
                    # Also, neither (i,j) nor (i,k) already appears in our open wedge set.
                    push!(Is,ij2)
                    push!(Js,ij1)
                    push!(Is,ik2)
                    push!(Js,ik1)
                    B[ij1,ij2] = 2
                    B[ik1,ik2] = 2

                    # Now exit this and start with a new j, as we have already
                    # used edge (i,j) and we don't need to check for more
                    # edges involving (i,j)
                    break
                end
            end
        end
    end
    Ib, Jb, Vb = findnz(triu(B'))
    Elist = [Ib Jb]

    # Entry will be 1 if the corresponding edge in Elist 
    # is deleted. It will be zero otherwise.
    soln = Vb .- 1         
    bd = sum(soln)/2
    return round(Int64,bd), Elist, soln
end

"""
Checks whether a solution to the STC LP
satisfies the constraints of the CD LP

Elist is an m x 2 array of edges

soln[k] is the variable for edge Elist[k,:].

"""
function check_cd_constraints(A,Elist,soln)
    
    # map from edge to edge index
    edgeDict = Dict()
    for k = 1:size(Elist,1)
        @assert(Elist[k,1]< Elist[k,2])
        edgeDict[(Elist[k,1],Elist[k,2])] = k
    end
    T,W = get_wedges_triangles(A)

    for t in T
        i = t[1]
        j = t[2]
        k = t[3]

        ij = edgeDict[(min(i,j),max(i,j))]
        ik = edgeDict[(min(i,k),max(i,k))]
        jk = edgeDict[(min(j,k),max(j,k))]

        tol = 1e-8
        if soln[ij] + soln[ik] - soln[jk] < -tol
    
            return false
        end
        if soln[ij] + soln[jk] - soln[ik] < -tol
            return false
        end
        if soln[jk] + soln[ik] - soln[ij] < -tol
            return false
        end
    end
    return true
end