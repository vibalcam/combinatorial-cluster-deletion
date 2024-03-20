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
save_path = "othergraph-results"

graphs = [
    "polbooksA";
    "adjnounA";
    "footballA";
    "Harvard500A";
    "Erdos991A"; 
    "celegansneuralA";
    "Netscience";
    "celegansmetabolicA";
    "RogetA";
    "SmaGriA";
    "emailA";
    "polblogsA"
]
data_path = "../../data/smallgraphs"

## Run experiments in order of size
for i = 1:length(graphs)
    graph = graphs[i]
    graph = graphs[i]
    if isfile("$data_path/$graph.mat")
        F = matread("$data_path/$graph.mat")
    else
        println("File $graph not exist. Continuing...")
    end
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
    matwrite("$(save_path)/$graph-merge-results.mat", Dict("m"=>m,"n"=>n,"mergetime"=>mergetime,"degtime"=>degtime,"mergeclus"=>mergeclus,"degclus"=>degclus,"degrat"=>degrat,"mergerat"=>mergerat))
end
