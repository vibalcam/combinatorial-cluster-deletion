using MAT
using SparseArrays
using LinearAlgebra
include("../../src/cd_lp_minst_relaxations.jl")
include("../../src/cd_det_mfp.jl")
include("../../src/cd_lp_relaxations.jl")
include("../../src/cd_rounding.jl")
include("../../src/faster_many_pivot.jl")
include("../run-cd-functions.jl")
include("mfp-merge.jl")
include("../../src/modular_deterministic_pivots.jl")

# choose path to data and saving results
save_path = "fb-results"

# In order to run these experiments, access to Facebook100 datasets is needed
fbs = readlines("Facebook_Sets.txt")
graphs = split(fbs[1])
data_path = "../../data/Facebook100" # assumes access to Facebook100 datasets
mvals = zeros(100)
nvals = zeros(100)
for i = 1:length(graphs)
    graph = graphs[i]
    F = matread("$data_path/$graph.mat")
    A = F["A"]
    n = size(A,1)
    m = sum(A)/2
    mvals[i] = m
    nvals[i] = n
end
p = sortperm(mvals)
## Run experiments in order of size
for i = 1:length(graphs)
    graph = graphs[p[i]]
    F = matread("$data_path/$graph.mat")
    A = F["A"]
    n = size(A,1)
    m = sum(A)/2
    mvals[i] = m
    nvals[i] = n
    
    Is, Js, Vs = findnz(triu(A))
    Elist = [Is Js]

    seed = 37

    # Get the lower bound: (MATCH STEP)
    Random.seed!(seed)
    lbtime = @elapsed bd, Elist_match, soln_match = maximal_disjoint_openwedge_rounding(A)

    # Construct the new graph: (FLIP STEP)
    fliptime = @elapsed Anew = flip_graph(n,Elist,soln_match)

    # Run Pivot to get clustering: (PIVOT STEP)


    # degree-based rounding
    degtime = @elapsed degclus = det_pivot_degree(Anew)
    degobj = check_cd_obj(Elist,degclus)
    degrat = degobj/bd
    
    mergetime = @elapsed mergeclus = mfpmerge(A, degclus)
    mergeobj = check_cd_obj(Elist, mergeclus)
    mergerat = mergeobj/bd
 
    println("")
    println("$(graph): \n \t\t $mergetime $degtime $lbtime $fliptime \t\t $mergerat $degrat")
    matwrite("$(save_path)/$(graph)-merge-results.mat", Dict("m"=>m,"n"=>n,"mergetime"=>mergetime,"degtime"=>degtime,"mergeclus"=>mergeclus,"degclus"=>degclus,"degrat"=>degrat,"mergerat"=>mergerat))
end
