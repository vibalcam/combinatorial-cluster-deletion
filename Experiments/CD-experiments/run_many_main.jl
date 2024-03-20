"""
Run CD methods one by one and save results for each method and graph
Uses run-cd-functions.jl
"""

using MAT
using SparseArrays

include("../../src/cd_lp_minst_relaxations.jl")
include("../../src/cd_det_mfp.jl")
include("../../src/cd_lp_relaxations.jl")
include("../../src/cd_rounding.jl")
include("../../src/faster_many_pivot.jl")
include("../run-cd-functions.jl")

graphs = [
    "simple-amazon0302";
    "simple-amazon0312";
    "simple-amazon0505";
    "simple-amazon0601";
    "simple-ca-AstroPh";
    "simple-ca-CondMat";
    "simple-ca-GrQc";
    "simple-ca-HepPh";
    "simple-ca-HepTh";
    "simple-cit-HepPh";
    "simple-cit-HepTh";
    "simple-cit-Patents";
    "simple-com-Amazon";
    "simple-com-DBLP";
    "simple-com-LiveJournal";
    "simple-com-Youtube";
    "simple-email-Enron";
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
    ########################################
    "simple-soc-LiveJournal1";
    "simple-com-Orkut";
]

# choose path to data and saving results
# data_path = "../../data/graphs"
data_path = "../../data/simple-snap"
# save_path = "cd-experiments-server"
save_path = "cd-experiments-local"

# time_limit = 1800
# pivot_times = 100
time_limit = 172800
pivot_times = 100

###################################
# MFP 
###################################
for graph in graphs
    tl = time_limit
    pivtimes = pivot_times

    # Load graph
    if isfile("$data_path/$graph.mat")
        F = matread("$data_path/$graph.mat")
    else
        println("File $graph not exist. Continuing...")
        continue
    end
    A = F["A"]
    n = size(A,1)
    A = sparse(A)
    Is, Js, Vs = findnz(triu(A))
    Elist = [Is Js]
    println("$graph mfp")
    try
        clus, runtimes, results = run_matchflippivot_cd(A,Elist,pivtimes)
        matwrite("$(save_path)/$(graph)_mfp.mat", Dict("clus"=>clus,
            "runtimes"=>runtimes,"results"=>results, "finished"=>true))
    catch e
        if isa(e,OutOfMemoryError)
            println("$(graph) mfp failed")
            matwrite("$(save_path)/$(graph)_mfp.mat", Dict("finished"=>false))
        else
            rethrow(e)
        end
    end
end

###################################
# Deterministic MFP Faster RatMFP
###################################
for graph in graphs
    tl = time_limit
    pivtimes = 1

    # Load graph
    if isfile("$data_path/$graph.mat")
        F = matread("$data_path/$graph.mat")
    else
        println("File $graph not exist. Continuing...")
        continue
    end
    A = F["A"]
    n = size(A,1)
    A = sparse(A)
    Is, Js, Vs = findnz(triu(A))
    Elist = [Is Js]
    println("$graph deterministic RatMFP")

    try
        clus, runtimes, results = run_det_matchflippivot_cd(A,Elist,pivtimes, false)
        matwrite("$(save_path)/$(graph)_detmfpfaster.mat", Dict("clus"=>clus,
            "runtimes"=>runtimes,"results"=>results, "finished"=>true))
    catch e
        if isa(e,OutOfMemoryError)
            println("$(graph) deterministic RatMFP failed")
            matwrite("$(save_path)/$(graph)_detmfpfaster.mat", Dict("finished"=>false))
        else
            rethrow(e)
        end
    end
end

###################################
# Deterministic MFP Faster DegMFP
###################################
for graph in graphs
    tl = time_limit
    pivtimes = 1

    # Load graph
    if isfile("$data_path/$graph.mat")
        F = matread("$data_path/$graph.mat")
    else
        println("File $graph not exist. Continuing...")
        continue
    end
    A = F["A"]
    n = size(A,1)
    A = sparse(A)
    Is, Js, Vs = findnz(triu(A))
    Elist = [Is Js]
    println("$graph deterministic DegMFP")

    try
        clus, runtimes, results = run_det_matchflippivot_cd(A,Elist,pivtimes, true)
        matwrite("$(save_path)/$(graph)_detdegmfp.mat", Dict("clus"=>clus,
            "runtimes"=>runtimes,"results"=>results, "finished"=>true))
    catch e
        if isa(e,OutOfMemoryError)
            println("$(graph) deterministic DegMFP failed")
            matwrite("$(save_path)/$(graph)_detdegmfp.mat", Dict("finished"=>false))
        else
            rethrow(e)
        end
    end
end

###################################
# STC LP Min s-t
###################################
for graph in graphs
    tl = time_limit; pivtimes = pivot_times; 

    if isfile("$data_path/$graph.mat")
        F = matread("$data_path/$graph.mat")
    else
        println("File $graph not exist. Continuing...")
        continue
    end
    
    # read matrix and convert to sparse
    A = F["A"]
    A = sparse(A)
    # find the nonzero elements of the upper triangular part and retrieve row, column, and value arrays
    Is, Js, Vs = findnz(triu(A))
    # create edge list (edges, 2) concatenating rows and column arrays
    Elist = [Is Js]
    println("$graph stc min s-t")

    try
        clus, runtimes, results = run_stcminst(A,Elist,pivtimes, tl, false)
        matwrite("$(save_path)/$(graph)_stcminst.mat", Dict("clus"=>clus,
            "runtimes"=>runtimes,"results"=>results, "finished"=>true))
    catch e
        if isa(e,OutOfMemoryError)
            println("$(graph) stc min s-t failed")
            matwrite("$(save_path)/$(graph)_stcminst.mat", Dict("finished"=>false))
        else
            rethrow(e)
        end
    end
end

###################################
# STC LP
###################################
for graph in graphs
    tl = time_limit; pivtimes = pivot_times; 

    if isfile("$data_path/$graph.mat")
        F = matread("$data_path/$graph.mat")
    else
        println("File $graph not exist. Continuing...")
        continue
    end
    A = F["A"]; A = sparse(A); Is, Js, Vs = findnz(triu(A)); Elist = [Is Js]
    println("$graph stc lp")
    try
        clus, runtimes, results = run_stclp(A,Elist,pivtimes,tl, false)
        matwrite("$(save_path)/$(graph)_stclp.mat", Dict("clus"=>clus,
            "runtimes"=>runtimes,"results"=>results, "finished"=>true))
    catch e
        if isa(e,OutOfMemoryError)
            println("$(graph) stc lp failed")
            matwrite("$(save_path)/$(graph)_stclp.mat", Dict("finished"=>false))
        else
            rethrow(e)
        end
    end
end

###################################
# STC LP Min s-t
###################################
for graph in graphs
    tl = time_limit; pivtimes = pivot_times; 

    if isfile("$data_path/$graph.mat")
        F = matread("$data_path/$graph.mat")
    else
        println("File $graph not exist. Continuing...")
        continue
    end
    
    # read matrix and convert to sparse
    A = F["A"]
    A = sparse(A)
    # find the nonzero elements of the upper triangular part and retrieve row, column, and value arrays
    Is, Js, Vs = findnz(triu(A))
    # create edge list (edges, 2) concatenating rows and column arrays
    Elist = [Is Js]
    println("$graph stc min s-t")

    try
        clus, runtimes, results = run_stcminst(A,Elist,pivtimes, tl, true)
        matwrite("$(save_path)/$(graph)_stcminst_det.mat", Dict("clus"=>clus,
            "runtimes"=>runtimes,"results"=>results, "finished"=>true))
    catch e
        if isa(e,OutOfMemoryError)
            println("$(graph) stc min s-t failed")
            matwrite("$(save_path)/$(graph)_stcminst_det.mat", Dict("finished"=>false))
        else
            rethrow(e)
        end
    end
end

###################################
# STC LP DETERMINISTIC ROUNDING
###################################
for graph in graphs
    tl = time_limit; pivtimes = pivot_times; 

    if isfile("$data_path/$graph.mat")
        F = matread("$data_path/$graph.mat")
    else
        println("File $graph not exist. Continuing...")
        continue
    end
    A = F["A"]; A = sparse(A); Is, Js, Vs = findnz(triu(A)); Elist = [Is Js]
    println("$graph stc lp")
    try
        clus, runtimes, results = run_stclp(A,Elist,pivtimes,tl, true)
        matwrite("$(save_path)/$(graph)_stclp_det.mat", Dict("clus"=>clus,
            "runtimes"=>runtimes,"results"=>results, "finished"=>true))
    catch e
        if isa(e,OutOfMemoryError)
            println("$(graph) stc lp failed")
            matwrite("$(save_path)/$(graph)_stclp_det.mat", Dict("finished"=>false))
        else
            rethrow(e)
        end
    end
end

# ###################################
# # STC LP Min s-t (Min source version)
# ###################################
# for graph in graphs
#     tl = time_limit; pivtimes = pivot_times; 

#     if isfile("$data_path/$graph.mat")
#         F = matread("$data_path/$graph.mat")
#     else
#         println("File $graph not exist. Continuing...")
#         continue
#     end
    
#     # read matrix and convert to sparse
#     A = F["A"]
#     A = sparse(A)
#     # find the nonzero elements of the upper triangular part and retrieve row, column, and value arrays
#     Is, Js, Vs = findnz(triu(A))
#     # create edge list (edges, 2) concatenating rows and column arrays
#     Elist = [Is Js]
#     println("$graph stc min s-t")

#     try
#         clus, runtimes, results = run_stcminst(A,Elist,pivtimes, tl, true)
#         matwrite("$(save_path)/$(graph)_stcminst_minsource.mat", Dict("clus"=>clus,
#             "runtimes"=>runtimes,"results"=>results, "finished"=>true))
#     catch e
#         if isa(e,OutOfMemoryError)
#             println("$(graph) stc min s-t failed")
#             matwrite("$(save_path)/$(graph)_stcminst_minsource.mat", Dict("finished"=>false))
#         else
#             rethrow(e)
#         end
#     end
# end
