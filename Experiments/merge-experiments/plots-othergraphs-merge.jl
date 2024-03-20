using Plots
using MAT
using SparseArrays
using LinearAlgebra

##
include("../../src/cd_lp_relaxations.jl")

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


mergerun = Vector{Float64}()
degrun = Vector{Float64}()
mergerat = Vector{Float64}()
degrat = Vector{Float64}()
mvals = Vector{Float64}()
nvals = Vector{Float64}()
g = Vector{String}()
L = Vector{Float64}()
for i = 1:length(graphs)
    graph = graphs[i]
        
    try
        
        M = matread("$(save_path)/$graph-merge-results.mat")
        F = matread("$data_path/$graph.mat")
        A = F["A"]
        n = size(A,1)
        m = sum(A)/2
                
        Is, Js, Vs = findnz(triu(A))
        Elist = [Is Js]
    
        seed = 37
    
        # Get the lower bound: (MATCH STEP)
        Random.seed!(seed)
        lbtime = @elapsed bd, Elist_match, soln_match = maximal_disjoint_openwedge_rounding(A)
        push!(L, bd)
        push!(mvals,M["m"])
        push!(nvals,M["n"])
        push!(mergerun, M["mergetime"])
        push!(degrun, M["degtime"])
        push!(mergerat,M["mergerat"])
        push!(degrat,M["degrat"])
        if graph[end] == 'A'
            push!(g,graph[1:end-1])
        else
            push!(g,graph)
        end
    catch
        # println("\t no results for $graph")
    end
end

## Plot

p = sortperm(degrat)
stepy = 1
ln = g[p]
x = mvals
s1 = 500
s2 = 400
title = ""
xl = "Problem instance"
yl = "\nApprox Ratio"
x = 1:(stepy):length(g)
f = plot(title = title,grid = false,legend = :bottomright,legendbox = false, xlabel = xl, ylabel = yl,foreground_color_legend = false)
scatter!(f,x,degrat[p],label = "RanMFP-Deg",size = (s1,s2),markerstrokecolor = :lightblue,background_color_legend = nothing,
markerstrokewidth = 6;color = :green,markershape = :cross, markersize = 8)
scatter!(f,x,mergerat[p],label = "   + Merge",markerstrokecolor = :red,
markerstrokewidth = 0;color = :red,markershape = :diamond, markersize = 6,
xticks = (1:stepy:length(ln), "        ".*ln[1:stepy:length(ln)]),xrotation = 50,
xtickfont=font(10),
ylim = [1.2,2.3],
ytickfont=font(11),
guidefont=font(13),
titlefont=font(12),
legendfont=font(10)
)
savefig("Figures/merge-other-graphs-ratios.pdf")


## Table of graph stats

p = sortperm(L./mvals)
for i = 1:length(g)
    m = round(Int64,mvals[p[i]])
    n = round(Int64,nvals[p[i]])
    w = round(200*L[p[i]]/m,digits = 2)
    dr = round(degrat[p[i]],digits = 2)
    mr = round(mergerat[p[i]],digits = 2)
    println("$(g[p[i]]) & $n & $m & $w & $dr& $mr \\\\")
end

# Rm = (2*L)./mvals
