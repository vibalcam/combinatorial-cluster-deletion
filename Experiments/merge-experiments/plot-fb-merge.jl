using Plots
using MAT
using SparseArrays
using LinearAlgebra
include("../../src/cd_lp_relaxations.jl")

##

fbs = readlines("Facebook_Sets.txt")
graphs = split(fbs[1])

save_path = "fb-results"
num = 46
mergerun = Vector{Float64}()
degrun = Vector{Float64}()
mergerat = Vector{Float64}()
degrat = Vector{Float64}()
mvals = Vector{Float64}()
nvals = Vector{Float64}()
g = Vector{String}()
L = Vector{Float64}()
data_path = "../../data/Facebook100" # assumes access to Facebook100 datasets
for i = 1:length(graphs)
    graph = graphs[i]
        
    try
        println("have results for $graph")
        
        M = matread("$(save_path)/($graph)-merge-results.mat")
        F = matread("$data_path/$graph.mat")
        println("yep")
        A = F["A"]
        n = size(A,1)
        m = sum(A)/2
                
        Is, Js, Vs = findnz(triu(A))
        Elist = [Is Js]
    
        seed = 37
    
        # Get the lower bound: (MATCH STEP)
        Random.seed!(seed)
        lbtime = @elapsed bd, Elist_match, soln_match = maximal_disjoint_openwedge_rounding(A)
        @show bd
        push!(L, bd)
        push!(mvals,M["m"])
        push!(nvals,M["n"])
        push!(mergerun, M["mergetime"])
        push!(degrun, M["degtime"])
        push!(mergerat,M["mergerat"])
        push!(degrat,M["degrat"])
        push!(g,graph)
    catch
        # println("\t no results for $graph")
    end
end

## Plot

p = sortperm(degrat)
stepy = 1
ln = g[p]
x = mvals
s1 = 600
s2 = 250
title = ""
xl = "Problem instance"
yl = "\nApprox Ratio"
x = 1:(stepy):length(g)
f = plot(title = title,grid = false,legend = :bottomright,legendbox = false, xlabel = xl, ylabel = yl,foreground_color_legend = false)
scatter!(f,x,degrat[p],label = "RanMFP-Deg",size = (s1,s2),markerstrokecolor = :lightblue,background_color_legend = nothing,
markerstrokewidth = 6;color = :green,markershape = :cross, markersize = 8)
scatter!(f,x,mergerat[p],label = "   + Merge",markerstrokecolor = :red,
markerstrokewidth = 0;color = :red,markershape = :diamond, markersize = 6,
xticks = (1:stepy:length(ln), "                           ".*ln[1:stepy:length(ln)]),xrotation = 50,
xtickfont=font(7),
ylim = [1.83,2.01],
ytickfont=font(11),
guidefont=font(13),
titlefont=font(12),
legendfont=font(10)
)
savefig("Figures/merge-ratios.pdf")


## Show lower bounds versus edges

for i = 1:length(g)
    println("$(mvals[i]) $(L[i]) $(mvals[i]/L[i])")
end

Rm = (2*L)./mvals
