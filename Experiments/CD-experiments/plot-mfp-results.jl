"""
Prints the results from mfp-tests.jl
"""

using Plots
using MAT

M = matread("cd-experiments-local/snap-mfp-results.mat")
Run = M["Runtimes"]
Obj = M["Objectives"]
mvals = M["mvals"]
nvals = M["nvals"]
graphs = M["graphs"]
##

g = Vector{String}()
for i = 1:length(graphs)
    push!(g,graphs[i][8:end])
end
g

##
# Runtimes[i,:] = [lbtime, fliptime, pivavgtime, pivmanytime, degtime,rattime]
# Objectives[i,:] = [bd, pivavgobj, pivmanyobj, degobj, ratobj]

Lb = Objectives[:,1]
Obs = Objectives[:,2:5]./Lb
##

p = sortperm(Obs[:,1])
Os = Obs[p,:]
stepy = 1
ln = g
x = mvals
s1 = 500
s2 = 250
ms = 6
title = ""
xl = "Problem instance"
yl = "Approx Ratio"
f = plot(title = title,grid = false,legend = :bottomleft,legendbox = false, xlabel = xl, ylabel = yl,foreground_color_legend = nothing)
scatter!(f,1:length(graphs),Os[:,1],label = "RanMFP-Avg",size = (s1,s2),markerstrokecolor = :lightblue,
markerstrokewidth = 0;color = :lightblue,markershape = :star, markersize = ms)
scatter!(f,1:length(graphs),Os[:,2],label = "RanMFP-100",markerstrokecolor = :blue,
markerstrokewidth = 0;color = :blue,markershape = :star, markersize = ms)
scatter!(f,1:length(graphs),Os[:,3],label = "DegMFP",markerstrokecolor = :green,
markerstrokewidth = 8;color = :green,markershape = :cross, markersize = 8)
scatter!(f,1:length(graphs),Os[:,4],label = "RatMFP",markerstrokecolor = :orange,
markerstrokewidth = 0;color = :orange,markershape = :circle, markersize = 4,
xticks = (1:stepy:length(ln), "            ".*ln[1:stepy:length(ln)]),xrotation = 40,
xtickfont=font(8),
ytickfont=font(9),
guidefont=font(10),
titlefont=font(10),
legendfont=font(9)
)

savefig("Figures/mfp-ratios.pdf")
# 

## Runtime Plot


x = mvals
s1 = 420
s2 = 300
ms = 6
title = ""
xl = "Number of edges |E|"
yl = "Runtime (s)"
f = plot(title = title,grid = false,legend = :topleft,legendbox = false, xlabel = xl, ylabel = yl,foreground_color_legend = nothing)
scatter!(f,x,Run[:,3],label = "RanMFP-Avg",size = (s1,s2),yscale = :log10, xscale = :log10,background_color_legend = nothing,
markerstrokewidth = 0;color = :lightblue,markershape = :star, markersize = ms,markerstrokecolor = :lightblue)
scatter!(f,x,Run[:,4],label = "RanMFP-100",markerstrokecolor = :blue,
markerstrokewidth = 0;color = :blue,markershape = :star, markersize = ms)
scatter!(f,x,Run[:,5],label = "DegMFP",markerstrokecolor = :green,
markerstrokewidth = 3;color = :green,markershape = :cross, markersize = 6)
scatter!(f,x,Run[:,6],label = "RatMFP",markerstrokecolor = :orange,
markerstrokewidth = 0;color = :orange,markershape = :circle, markersize = 4,
# xticks = ((1e5, 1e6, 1e7),[1e5, 1e6, 1e7]),
# xlim = [1e5,1e8],
ylim = [0.01,10^4],
# yticks = [0 1 10 100 1000],
# yformatter = :scientific,
# yticks = ((10^-1, 1, 10, 100,1000),[0.1, 10^0, 10^1, 100,1000]),
xtickfont=font(15),
ytickfont=font(15),
guidefont=font(16),
# titlefont=font(10),
legendfont=font(10)
)

savefig("Figures/mfp-runtimes.pdf")
