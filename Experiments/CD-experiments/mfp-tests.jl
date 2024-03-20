"""
Run the MFP methods in a modular and faster manner
Reuses the lower bound to not have to recompute it
"""

using MAT
using SparseArrays
using LinearAlgebra
include("../../src/cd_lp_minst_relaxations.jl")
include("../../src/cd_det_mfp.jl")
include("../../src/cd_lp_relaxations.jl")
include("../../src/cd_rounding.jl")
include("../../src/faster_many_pivot.jl")
include("../run-cd-functions.jl")
# include("mfp-merge.jl")
include("../../src/modular_deterministic_pivots.jl")

tinytestgraphs = [
    "KarateA";
    "dolphinsA";
    "lesmisA";
    "polbooksA";
    "adjnounA";
    "footballA";
]

medgraphs = [
        "Harvard500A";
        "Erdos991A"; 
        "celegansneuralA";
        "Netscience";
        "celegansmetabolicA";
        "RogetA";
        "SmaGriA";
        "emailA"
]

biggraphs = [

    # larger graphs
    "simple-amazon0302";
    "simple-amazon0312";
    "simple-amazon0505";
    "simple-amazon0601";
    "simple-cit-Patents";
    "simple-com-Amazon";
    "simple-com-DBLP";
    "simple-com-LiveJournal";
    "simple-com-Youtube";
    "simple-email-EuAll";
    "simple-loc-Brightkite";
    "simple-loc-Gowalla";
    "simple-roadNet-CA";
    "simple-roadNet-PA";
    "simple-roadNet-TX";
    "simple-soc-Epinions1";
    "simple-soc-Slashdot0811";
    "simple-soc-Slashdot0902";
    "simple-web-BerkStan";
    "simple-web-Google";
    "simple-web-NotreDame";
    "simple-web-Stanford";
    "simple-wiki-Talk";
    "simple-wiki-topcats";
    # ########################################
        # "simple-soc-LiveJournal1";
]

# choose path to data and saving results
# data_path = "../../data/smallgraphs"
data_path = "../../data/simple-snap"
save_path = "cd-experiments-local"

## Test just one

# try to merge clusters
# graphs = tinytestgraphs
# graphs = medgraphs
graphs = biggraphs
time_limit = 172800
pivtimes = 100

Objectives = zeros(length(graphs),5)
Runtimes = zeros(length(graphs),6)
# Load graph
for i = 1:length(graphs)
    graph = graphs[i]
    if isfile("$data_path/$graph.mat")
        F = matread("$data_path/$graph.mat")
    else
        println("File $graph not exist. Continuing...")
    end
    A = F["A"]
    n = size(A,1)
    A = sparse(A)
    Is, Js, Vs = findnz(triu(A))
    Elist = [Is Js]

    seed = 37

    # Get the lower bound: (MATCH STEP)
    Random.seed!(seed)
    lbtime = @elapsed bd, Elist_match, soln_match = maximal_disjoint_openwedge_rounding(A)

    # Construct the new graph: (FLIP STEP)
    fliptime = @elapsed Anew = flip_graph(n,Elist,soln_match)

    # Run Pivot to get clustering: (PIVOT STEP)

    # Random pivot, two ways
    Random.seed!(seed)
    pivobjs = zeros(pivtimes)
    pivruns = zeros(pivtimes)
    checktimes = zeros(pivtimes)
    for jj = 1:pivtimes
        ptime = @elapsed clus = permutation_pivot(Anew)
        checktime = @elapsed objnew = check_cd_obj(Elist,clus)
        pivruns[jj] = ptime
        pivobjs[jj] = objnew
        checktimes[jj] = checktime
    end
    pivmanytime = sum(checktimes) + sum(pivruns) + fliptime + lbtime
    pivavgtime = sum(pivruns)/pivtimes + fliptime + lbtime

    pivmanyobj = minimum(pivobjs)
    pivavgobj = sum(pivobjs)/pivtimes

    pivmanyrat = pivmanyobj/bd
    pivavgrat = pivavgobj/bd

    # degree-based rounding
    degtime = @elapsed clus_deg = det_pivot_degree(Anew)
    degobj = check_cd_obj(Elist,clus_deg)
    rattime = @elapsed clus_rat = det_pivot_ratio(Anew)
    ratobj = check_cd_obj(Elist,clus_rat)
    degrat = degobj/bd
    ratrat = ratobj/bd

    println("")
    println("$(graph): $pivmanytime $degtime $rattime $pivmanyrat $degrat $ratrat ; $lbtime $fliptime")
    Runtimes[i,:] = [lbtime, fliptime, pivavgtime, pivmanytime, degtime,rattime]
    Objectives[i,:] = [bd, pivavgobj, pivmanyobj, degobj, ratobj]
    matwrite("$(save_path)/sofar_mfp.mat", Dict("Runtimes"=>Runtimes,"Objectives"=>Objectives))
end

## Get graph stats quick
mvals = zeros(length(graphs))
nvals = zeros(length(graphs))
for i = 1:length(graphs)
    graph = graphs[i]
    if isfile("$data_path/$graph.mat")
        F = matread("$data_path/$graph.mat")
    else
        println("File $graph not exist. Continuing...")
    end
    A = F["A"]
    n = size(A,1)
    m = nnz(A)/2
    mvals[i] = m
    nvals[i] = n
end

## Quick plot output
using Plots
f = scatter(mvals,Objectives[:,2]./Objectives[:,1],yscale = :normal, xscale = :log10,grid = false, label = "RanMFP-Avg")
scatter!(mvals,Objectives[:,3]./Objectives[:,1],label = "RanMFP-100")
scatter!(mvals,Objectives[:,4]./Objectives[:,1],label = "DegMFP")
scatter!(mvals,Objectives[:,5]./Objectives[:,1],label = "RatMFP")

