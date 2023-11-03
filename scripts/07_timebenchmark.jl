using Pesto
using Distributions
using ProgressMeter
using CSV
using DataFrames
using JLD2
using StatsPlots

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

parms = zeros(
    n_benchmarks,
    3, ## lambda, mu, eta
)

## fit parameters once
bmodels = Array{SSEconstant, 1}(undef, n_benchmarks)
@showprogress for i in 1:n_benchmarks
    data = benchmark_datasets[i]
    ρ = 1.0
    λml, μml, = estimate_constant_bdp(data)

    H = 0.587405
    d1 = LogNormal(log(λml), H)
    d2 = LogNormal(log(μml), H)
    n = 6
    speciation = make_quantiles(d1, n)
    extinction = make_quantiles(d2, n)
    
    λ, μ = allpairwise(speciation, extinction)
    η = optimize_eta(λ, μ, data)
    parms[i,:] .= [λml, μml, η]
    model = SSEconstant(λ, μ, η)
    bmodels[i] = model
end


n_replicates = 5
times1 = zeros(
    n_benchmarks,
    4, ## D(t), F(t), N(t) and η estimation
    n_replicates
)


#i = 1
@showprogress for i in 1:n_benchmarks
    data = benchmark_datasets[i]
    model = bmodels[i]
    Ds, Fs = backwards_forwards_pass(model, data)
    Ss = ancestral_state_probabilities(data, Ds, Fs)
    λ = model.λ
    μ = model.μ

    for j in 1:n_replicates
        ## logL
        GC.gc()
        t = @elapsed logL_root(model, data)
        times1[i,1,j] = t 

        ## D+F
        GC.gc()
        t = @elapsed backwards_forwards_pass(model, data)
        times1[i,2,j] = t

        ## N
        GC.gc()
        t = @elapsed compute_nshifts(model, data, Ds, Ss)
        times1[i,3,j] = t

        ## eta
        GC.gc()
        t = @elapsed optimize_eta(λ, μ, data)
        times1[i,4,j] = t
    end
end

#save("output/benchmark_times.jld2", "times1", times1)
times1 = load("output/benchmark_times.jld2")["times1"]

ntaxa = [
    length(data.tiplab) for data in benchmark_datasets
]

p = StatsPlots.plot(grid = false, 
        #xlabel = "number of taxa", ylabel = "time (s)", 
            #yscale = :log10,
            xticks = (collect(1:n_benchmarks), ntaxa))
labels = ["logL", "D(t)+F(t)", "N(t)", "η"]
cls = palette(:tab10, 4)
offsets = [-0.15, -0.05, 0.05, 0.15] .* 1.0

using Statistics

ys = Statistics.mean(times1, dims = 3)[:,:,1]
ys = Statistics.median(times1, dims = 3)[:,:,1]
xs = hcat([collect(1:n_benchmarks) .+ 0 for offset in offsets]...)
vars = zeros(n_benchmarks, 4)
quants = zeros(n_benchmarks, 4, 2)
interquart = zeros(n_benchmarks, 4)
for i in 1:10
    for j in 1:4
        vars[i,j] = sum((times1[i,j,:] .- ys[i,j]).^2) / (n_replicates-1)
        quants[i,j,:] = quantile(times1[i,j,:], [0.25, 0.75])
        interquart[i,j] = abs(diff(quants[i,j,:])[1])
    end
end
sds = sqrt.(vars)


for q in 1:4
    StatsPlots.plot!(p, xs[:,q], ys[:,q], yscale = :log10, 
            yerror = interquart[:,q], 
            label = labels[q], color = cls[q], linewidth = 2)
end
yt = Pesto.lrange(0.01, 1000.0, 6)
StatsPlots.plot!(p, legend = :topleft,
    yticks = (yt, round.(yt; digits = 2)),
    ylabel = "time (s)", xlabel = "number of taxa", 
    size = (300, 300), xrotation = 90)
    
StatsPlots.savefig(p, "figures/timebenchmark_algorithms.pdf")












## Set up the model
model = SSEconstant(λ, μ, η)
## Compute the results
results = Dict()
times = Dict()
@showprogress for (ntaxa, d) in data
    t1 = time()
    Ds, Fs = backwards_forwards_pass(model, d; verbose = false)
    res = calculate_tree_rates(d, model, Ds, Fs; verbose = false);
    t2 = time()
    times[ntaxa] = t2 - t1

    results[ntaxa] = res
end

rdf = CSV.read("output/runtimes_r.csv", DataFrame)

jdf = DataFrame(time = collect(values(times)), ntaxa = collect(keys(times)))
sort!(jdf, :ntaxa)

ps = []
for i in 1:2
    p = plot(xlab = "ntaxa", ylab = "time (seconds)", legend = :bottomright, title = "Arithmetic")

    if i == 1
        plot!(p, xscale = :log10, yscale = :log10, title = "Log-scale")

        for i in 1:length(rdf[:,:ntaxa])
            annotate!(p, rdf[i, :ntaxa], rdf[i, :tnativeR], text(rdf[:,:ntaxa][i], 6, :left, :top))
        end
        for i in 1:length(jdf[:,:ntaxa])
            annotate!(p, jdf[i, :ntaxa], jdf[i, :time], text(jdf[:,:ntaxa][i], 6, :left, :top))
        end
    end
    plot!(p, jdf[:, :ntaxa], jdf[:, :time], label = "", color = "black")
    scatter!(p, jdf[:, :ntaxa], jdf[:, :time], label = "Julia", color = "black")

    ## R native
    plot!(p, rdf[:, :ntaxa], rdf[:, :tnativeR], label = "", color = "orange")
    scatter!(p, rdf[:, :ntaxa], rdf[:, :tnativeR], label = "R", color = "orange")

    ## Rcpp
    plot!(p, rdf[:, :ntaxa], rdf[:, :tRcpp], label = "", color = "green")
    scatter!(p, rdf[:, :ntaxa], rdf[:, :tRcpp], label = "Rcpp (postorder)", color = "green")


    append!(ps, [p])
end

runtimeplot = plot(ps...)

savefig(runtimeplot, "figures/runtime_36states_100knotsperbranch.pdf")








