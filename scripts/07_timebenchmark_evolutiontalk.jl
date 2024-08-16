using Pesto
using Distributions
using ProgressMeter
using CSV
using DataFrames
using JLD2
using StatsPlots
using Statistics

#### load simulated trees
## hard code tree heights
tree_heights = [25.0, 50.0, 75.0, 100.0, 125.0]
n_heights = length(tree_heights)
n_starting = 1
starting = collect(1:n_starting)
n_models = 4
n_trees = 100


fpaths = Array{String, 4}(undef, n_trees, n_heights, n_starting, n_models)
datasets = Array{SSEdata, 4}(undef, n_trees, n_heights, n_starting, n_models)
taxa = zeros(n_trees, n_heights, n_starting, n_models)
treelengths = zeros(n_trees, n_heights, n_starting, n_models)

for i in 1:n_trees
    for j in 1:n_heights
        for k in 1:n_starting
            for m in 1:n_models
                s = starting[k]
                h = tree_heights[j]
                fpath = string("simulated_trees/s",s,"_h",h,"_m",m,"_",i,".tre")
                fpaths[i,j,k,m] = fpath
                phy = readtree(fpath)
                ρ = 1.0
                data = SSEdata(phy, ρ)
                datasets[i,j,k,m] = data
                taxa[i,j,k,m] = length(data.tiplab)
                treelengths[i,j,k,m] = sum(data.branch_lengths)
            end
        end
    end
end

datasets[1,1,1,1]
ntips_goal = Pesto.lrange(100.0, 20000.0, 10)

coords = [
    argmin((ntips .- taxa) .^ 2 ) for ntips in ntips_goal
]

benchmark_datasets = [
    datasets[c] for c in coords
]

n_benchmarks = length(benchmark_datasets)
n_replicates = 5

times1 = load("output/benchmark_times.jld2")["times1"]

ntaxa = [
    length(data.tiplab) for data in benchmark_datasets
]

p = StatsPlots.plot(grid = false, 
        #xlabel = "number of taxa", ylabel = "time (s)", 
            #yscale = :log10,
            size = (400, 300),
            xticks = (collect(1:n_benchmarks), ntaxa))
labels = ["likelihood", "Pesto", "shifts", "shift rate"]
cls = palette(:tab10, 4)
offsets = [-0.15, -0.05, 0.05, 0.15] .* 1.0


ys = Statistics.mean(times1, dims = 3)[:,:,1]
ys = Statistics.median(times1, dims = 3)[:,:,1]
xs = hcat([collect(1:n_benchmarks) .+ 0 for offset in offsets]...)
vars = zeros(n_benchmarks, 4)
quants = zeros(n_benchmarks, 4, 2)
for i in 1:10
    for j in 1:4
        vars[i,j] = sum((times1[i,j,:] .- ys[i,j]).^2) / (n_replicates-1)
        quants[i,j,:] = quantile(times1[i,j,:], [0.25, 0.75])
    end
end


#for q in 1:4
for q in [2]
    StatsPlots.plot!(p, xs[:,q], ys[:,q], yscale = :log10, 
            label = labels[q], color = cls[q], linewidth = 2)
    StatsPlots.scatter!(p, xs[:,q], ys[:,q], yscale = :log10, 
            label = "", color = cls[q])
end
yt = Pesto.lrange(0.01, 10_000_000.0, 10)
ytl = ["0.01 s", "0.1 s", "1 s", "10 s", "", "", "3 hour", "1 day", "2 weeks", "4 months"]
StatsPlots.plot!(p, legend = :right,
    yticks = (yt, ytl),
    ylabel = "run time", xlabel = "number of taxa", 
     xrotation = 90)


revbayes_times = [
    30*60*60, ## about 30 hrs
    55*60*60, ## about 55 hrs
    72.5*60*60, ## about 72.5 hrs
    125.0*60*60, ## about 125 hr
    198*60*60 ## about 200 hrs
]
revbayes_ntaxa = [1,2,3,4,5]
StatsPlots.plot!(p, revbayes_ntaxa, revbayes_times, color = cls[1], linewidth = 2, label = "RevBayes")
StatsPlots.scatter!(p, revbayes_ntaxa, revbayes_times, color = cls[1], linewidth = 2, label = "")
StatsPlots.plot!(p, [5,10], [198*60*60, 144*24*60*60], color = cls[1], linewidth = 2, label = "", linestyle = :dash)

using Measures

StatsPlots.plot!(p,
    ylims = (0.01, 150*60^3),
    bottom_margin = 3mm
)
p
    
StatsPlots.savefig(p, "/home/bkopper/evolutiontalk_pesto_revbayes_runtime.pdf")














