"""
Prints the results from run_many_main.jl

Code to print the results of all the cluster deletion experiments in a latex table.
It also creates the figures for the paper.
"""
#####################################
########### INCLUDES ################
#####################################

using MAT
using Plots
using SparseArrays
using LinearAlgebra
using Statistics
using Serialization
using Printf

include("../../src/helpers.jl")

#####################################
########### HELPER FUNCTIONS ########
#####################################

# Initialize dictionaries
function initialize_dict()
    return Dict(
        :lb => [],
        # lower bound
        :ub => [],
        # upper bound
        :rat => [],
        # ratio (upper/lower)
        :run_total => [],
        # total runtime (lower + rounding)
        :run => [],
        # runtime of lower bound
        :run_round => [],
        # runtime of rounding
        :status => [],  
        # if false, we didn't converge
        :iscd => [],    
        # if true, this is the cd solution too!
        :time_st => [],
        # time to run min st cut
        :avg_run_round => [],
        # average runtime of rounding
        :avg_ub => [],
        # average objective value
        :avg_ratio => [],
        # average ratio

        :perc_deleted => [],
        # percentage of delted edges
        :perc_singletons => [],
        # percentage of singleton clusters
        :avg_cluster_size => [],
        # average cluster size
        :m => [],
        # number of edges
        :n => [],
        # number of nodes
        :name => [],
        # name of the graph
    )
end

# Load graph data
function load_graph_data(graph, data_path)
    if isfile("$data_path/$graph.mat")
        F = matread("$data_path/$graph.mat")

        A = F["A"]
        n = size(A,1)
        # m = round(Int64,sum(A)/2)
        m = Int(sum(A)/2)
        
        A = sparse(A)
        Is, Js, Vs = findnz(triu(A))
        Elist = [Is Js]
        
        return A, Elist, n, m
    else
        return nothing
    end
end

function get_cluster_info(d, G, clus)
    A, Elist, n, m = G

    # clusters info
    # -- count size of clusters
    clusters = zeros(Int, maximum(clus))
    for i in clus
        clusters[i] += 1
    end

    obj = check_cc_obj_fastish(A,Elist,clus,m)
    # @assert(obj == cc_cd_obj(A,clus))

    # -- percentage of edges deleted
    perc_deleted = obj / m
    # -- number of singletons
    perc_singletons = count(x -> x == 1, clusters) / length(clusters)
    # -- avg cluster size
    avg_cluster_size = mean(clusters)

    # push values into dictionary
    push!(d[:perc_deleted], perc_deleted)
    push!(d[:perc_singletons], perc_singletons)
    push!(d[:avg_cluster_size], avg_cluster_size)
    push!(d[:m], m)
    push!(d[:n], n)

    return perc_deleted, perc_singletons, avg_cluster_size
end

# Process data and add it to dictionary
function process_data(type, M, gname, d, G, finished_key="finished", results_key="results", runtimes_key="runtimes", clusters_key="clus") #, dig = 3, dig2 = 4)
    # G = (A, Elist, n, m)

    push!(d[:name], gname)
    is_valid = false

    l_runtimes = missing
    l_obj = missing

    if !isnothing(M) && M[finished_key]
        res = M[results_key]
        time = M[runtimes_key]

        if length(time) > 2
            @assert(type == "mfp")
            # l_runtimes = mean(time[3:end])
            l_runtimes = mean(time[3])
            l_obj = mean(res[4])

            push!(d[:status], true)
            push!(d[:iscd], missing)
        else
            if length(res) >= 4
                push!(d[:status], Bool(res[4]))
                push!(d[:iscd], res[5])
            else
                push!(d[:status], true)
                push!(d[:iscd], missing)
            end
        end

        if length(res) >= 6
            push!(d[:time_st], res[6])
        else
            push!(d[:time_st], missing)
        end

        is_valid = d[:status][end]
    else
        push!(d[:status], true)
        push!(d[:iscd], missing)
        push!(d[:time_st], missing)

        is_valid = false
    end

    push!(d[:avg_run_round], l_runtimes)
    push!(d[:avg_ub], l_obj)

    if is_valid
        # approximation ratio info
        push!(d[:lb], res[1])
        push!(d[:ub], res[2])
        push!(d[:rat], d[:ub][end]/d[:lb][end])

        if !ismissing(l_obj)
            push!(d[:avg_ratio], d[:avg_ub][end]/d[:lb][end])
        end
        
        # runtime info
        push!(d[:run], time[1])
        push!(d[:run_round], time[2])
        push!(d[:run_total], time[1] + time[2])
        # runtime is the sum of lower bound and rounding

        # clusters info
        get_cluster_info(d, G, M[clusters_key])
    else
        push!(d[:lb], missing)
        push!(d[:ub], missing)
        push!(d[:rat], missing)
        push!(d[:run_total], missing)
        push!(d[:run], missing)
        push!(d[:run_round], missing)
        push!(d[:avg_ratio], missing)
    end

    return is_valid
end

#####################################
########### PARAMETERS ##############
#####################################

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
    "simple-com-Orkut";
    "simple-com-Youtube";
    "simple-email-Enron";
    "simple-email-EuAll";
    "simple-loc-Brightkite";
    "simple-loc-Gowalla";
    "simple-roadNet-CA";
    "simple-roadNet-PA";
    "simple-roadNet-TX";
    "simple-soc-Epinions1";
    "simple-soc-LiveJournal1";
    "simple-soc-Slashdot0811";
    "simple-soc-Slashdot0902";
    "simple-web-BerkStan";
    "simple-web-Google";
    "simple-web-NotreDame";
    "simple-web-Stanford";
    "simple-wiki-Talk";
    "simple-wiki-topcats";
]

# save_path = "cd-experiments-server"
save_path = "cd-experiments-local"

models_names = [
    ("mfp", "RanMFP"),
    ("detdegmfp", "DegMFP"),
    ("detmfpfaster", "RatMFP"),
    # ("stcminst", "Comb-LP"),
    ("stcminst_det", "Comb-LP"),
    # ("stcminst_minsource", "Comb-LP-MinSource"),
    # ("stclp", "Gurobi-LP"),
    ("stclp_det", "Gurobi-LP"),
]

# symbols for not finished or not converged
symbol_not_finished = "***"
symbol_not_converged = "--"
# data_path = "../../data/graphs"
data_path = "../../data/simple-snap"

dig = 3
dig2 = 4

sample_dict = initialize_dict()
models_data = Dict()


#####################################
########### PRINT TABLE ##############
#####################################

function commas(num)
    str = string(num)
    parts = split(str, ".")
    
    if length(parts) == 2
        # Apply commas to the left part only
        left_part = replace(parts[1], r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
        result_str = join([left_part, parts[2]], ".")
    else
        # No decimal point, apply commas to the whole string
        result_str = replace(str, r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
    end
    
    return result_str
end

function commas_print(num)
    print(commas(num))
end

# print table
print("\\begin{tabular}{ll")
for (i, j) in models_names
    print("r")
end
println("}")

# header for LaTeX table
println("\\toprule")
print("\\textbf{Graph} & ")
for (p, name) in models_names
    print("& \\textbf{$(name)}")

    # initialize dictionary for this model
    models_data[p] = initialize_dict()
end
println("\\\\*")

# check if data exists
file_name = "cd-experiments-data.jls"
if !isfile(file_name)
    for i = 1:length(graphs)
        graph = graphs[i]

        # Load graph data
        graph_data = load_graph_data(graph, data_path)
        if graph_data === nothing
            continue
        end

        A, Elist, n, m = graph_data
        # graph information
        gname = graph[8:end]

        # get data per model
        anyvalid = false
        for (p, name) in models_names
            if isfile("$(save_path)/$(graph)_$(p).mat")
                M = matread("$(save_path)/$(graph)_$(p).mat")
                # process_data(M, models_data[p])
            else
                M = nothing
            end
            anyvalid = process_data(p, M, gname, models_data[p], (A, Elist, n, m)) || anyvalid
        end

        # if no data available, skip
        if !anyvalid
            continue
        end

        # node and edges information
        N = commas(n)
        # if n > 10000
        #     N = round(Int64,n/1000)
        #     N = "$(N)k"
        # end
        M = commas(m)
        # if m > 10000
        #     M = round(Int64,m/1000)
        #     M = "$(M)k"
        # end

        println("\\midrule")
        
        # print name and lb
        print("\\textsc{$gname} & LB") 
        for (p, name) in models_names
            print(" & ")
            if !models_data[p][:status][end]
                print("$symbol_not_converged")
            elseif ismissing(models_data[p][:lb][end])
                print("$symbol_not_finished")
            else
                commas_print(models_data[p][:lb][end])
            end
        end
        println("\\\\*")

        # print ub
        print("& UB")
        for (p, name) in models_names
            print(" & ")
            if !models_data[p][:status][end]
                print("$symbol_not_converged")
            elseif ismissing(models_data[p][:ub][end])
                print("$symbol_not_finished")
            else
                commas_print(models_data[p][:ub][end])
            end
        end
        println("\\\\*")

        # print ratio
        print("\$n =\$ $N & Ratio ")
        for (p, name) in models_names
            print(" & ")
            if !models_data[p][:status][end]
                print("$symbol_not_converged")
            elseif ismissing(models_data[p][:rat][end])
                print("$symbol_not_finished")
            else
                print(round(models_data[p][:rat][end], sigdigits = dig2))
            end
        end
        println("\\\\* [2mm]")

        # print run lower bound
        print("\$m =\$ $M & Run LB ")
        for (p, name) in models_names
            print(" & ")
            if !models_data[p][:status][end]
                print("$symbol_not_converged")
            elseif ismissing(models_data[p][:run][end])
                print("$symbol_not_finished")
            else
                commas_print(@sprintf("%.3f", models_data[p][:run][end]))
            end
        end
        println("\\\\*")

        # print run rounding
        print("& Run rounding ")
        for (p, name) in models_names
            print(" & ")
            if !models_data[p][:status][end]
                print("$symbol_not_converged")
            elseif ismissing(models_data[p][:run_round][end])
                print("$symbol_not_finished")
            else
                commas_print(@sprintf("%.3f", models_data[p][:run_round][end]))
            end
        end
        println("\\\\*")

        # print run total
        print("& Run total ")
        for (p, name) in models_names
            print(" & ")
            if !models_data[p][:status][end]
                print("$symbol_not_converged")
            elseif ismissing(models_data[p][:run_total][end])
                print("$symbol_not_finished")
            else
                commas_print(@sprintf("%.3f", models_data[p][:run_total][end]))
            end
        end
        println("\\\\*")
    end

    println("\\bottomrule")
    println("\\end{tabular}")

    # # save julia data
    # open(file_name, "w") do file
    #     serialize(file, models_data)
    # end
else
    # load julia data
    models_data =  deserialize(open(file_name, "r"))
end

#####################################
########### PLOT FIGURES ############
#####################################

# helper function to get style for each model
function get_style(model_key)
    if model_key == "mfpavg"
        d =  Dict(
            :markerstrokecolor => :skyblue,
            :markerstrokewidth => 0,
            :color => :skyblue,
            :markershape => :star,
            :label => "RanMFP-Avg",
            :markersize => 6,
        )
    elseif model_key == "mfp"
        d =  Dict(
            :markerstrokecolor => :blue,
            :markerstrokewidth => 0,
            :color => :blue,
            :markershape => :star,
            :label => "RanMFP-100",
            :markersize => 6,
        )
    elseif model_key == "detmfpfaster"
        d =  Dict(
            :markerstrokecolor => :orange,
            :markerstrokewidth => 0,
            :color => :orange,
            :markershape => :circle,
            :markersize => 4,
            :label => "RatMFP",
        )
    elseif model_key in ["stcminst", "stcminst_det"]
        d =  Dict(
            :markerstrokecolor => :purple,
            :markerstrokewidth => 0,
            :color => :purple,
            :markershape => :utriangle,
            :label => "Comb-LP",
        )
    elseif model_key in ["stclp", "stclp_det"]
        d =  Dict(
            :markerstrokecolor => :red,
            :markerstrokewidth => 0,
            :color => :red,
            :markershape => :hexagon,
            :label => "Gurobi-LP",
        )
    elseif model_key == "detdegmfp"
        d =  Dict(
            :markerstrokecolor => :green,
            :markerstrokewidth => 8,
            :color => :green,
            :markershape => :cross,
            :label => "DegMFP",
            :markersize => 6,
        )
    else
        error("Unknown model key: $model_key")
    end

    return d
end

# graph information
m = models_data["detdegmfp"][:m]
n = models_data["detdegmfp"][:n]
gnames = models_data["detdegmfp"][:name]
println("Number of graphs: $(length(gnames))")

# log n and m values
n_values = log10.(n)
m_values = log10.(m)


############ Sorted ratio vs graphs
y_key = :rat
sort_idx = sortperm(models_data["mfp"][y_key])
title = ""
x = gnames[sort_idx]
xl = "Problem instance"
yl = "Approx Ratio"
ln = x
s1 = 500
s2 = 250
stepy = 1
f = plot(
    title = title, 
    grid = false, 
    legend = :best, 
    legendbox = false,
    xlabel = xl, 
    ylabel = yl,
    foreground_color_legend = nothing,
)
for (model_key, model_name) in models_names[1:3]
    scatter!(
        f,1:length(x),models_data[model_key][y_key][sort_idx],
        size = (s1,s2);
        get_style(model_key)...
    )
end
model_key = "mfpavg"
scatter!(
    f,1:length(x),models_data["mfp"][:avg_ratio][sort_idx],
    size = (s1,s2);
    xticks = (1:stepy:length(ln), "            ".*ln[1:stepy:length(ln)]),xrotation = 40,
    xtickfont=font(8),
    ytickfont=font(9),
    guidefont=font(10),
    titlefont=font(10),
    legendfont=font(9),
    get_style(model_key)...
)
savefig("Figures/mfp-ratios.pdf")

############ Runtime plots MFP

title = ""
x = m
y_key = :run_total
xl = "Number of edges |E|"
yl = "Runtime (s)"
s1 = 420
s2 = 300
stepy = 1
f = plot(
    title = title, 
    grid = false, 
    legend = :topleft, 
    legendbox = false,
    xlabel = xl, 
    ylabel = yl,
    foreground_color_legend = nothing,
)
for (model_key, model_name) in models_names[1:3]
    scatter!(
        f,x,models_data[model_key][y_key],
        size = (s1,s2);
        yscale = :log10, xscale = :log10,
        background_color_legend = nothing,
        get_style(model_key)...
    )
end
model_key = "mfpavg"
scatter!(
    f,x, models_data["mfp"][:avg_run_round] .+ models_data["mfp"][:run],
    label = "RanMFP-Avg",
    size = (s1,s2);
    yscale = :log10, xscale = :log10,
    background_color_legend = nothing,
    ylim = [0.01,10^4],
    xtickfont=font(15),
    ytickfont=font(15),
    guidefont=font(16),
    legendfont=font(12),
    get_style(model_key)...
)
savefig("Figures/mfp-runtimes.pdf") 

############ Runtime plots LP

title = ""
x = m
y_key = :run
xl = "Number of edges |E|"
yl = "Runtime (s)"
s1 = 420
s2 = 300
stepy = 1
global missing_fill = 4
f = plot(
    title = title, 
    grid = false, 
    legend = :topleft, 
    legendbox = false,
    xlabel = xl, 
    ylabel = yl,
    foreground_color_legend = nothing,
)
for (model_key, model_name) in models_names[4:5]
    global missing_fill
    missing_fill += 0.5

    run_values = models_data[model_key][y_key]
    missing_indices = findall(ismissing, run_values)
    non_missing_indices = findall(!ismissing, run_values)

    # plot vertical line for maximum non-missing value
    if length(missing_indices) > 0
        d = get_style(model_key)
        # remove label
        d[:label] = ""

        vline!(
            [maximum(x[non_missing_indices])], 
            linestyle = :dot, linewidth = 4, 
            color = d[:color], 
            label = "", 
            alpha=0.5, 
            xscale=:log10, yscale=:log10
        )

        # Identify missing values and plot them as crosses
        scatter!(
            f,x[missing_indices], fill(10^missing_fill, length(missing_indices)),
            size = (s1,s2);
            yscale = :log10, xscale = :log10,
            background_color_legend = nothing,
            d...
        )
    end

    scatter!(
        f,x[non_missing_indices],run_values[non_missing_indices],
        size = (s1,s2);
        yscale = :log10, xscale = :log10,
        background_color_legend = nothing,
        ylim = [0.01,10^5],
        xtickfont=font(15),
        ytickfont=font(15),
        guidefont=font(16),
        legendfont=font(12),
        get_style(model_key)...
    )
end
savefig("Figures/lp-runtimes.pdf") 
