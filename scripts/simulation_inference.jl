using Pesto
using JLD2
using Glob
using ProgressMeter
using Distributions

dat = load("simulated_trees/N.jld2")
N = dat["N"]

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
Nt = sum(N, dims = (5,6))[1:n_trees,:,:,:,1,1] ./ treelengths
Nsum = sum(N, dims = (5,6))[:,:,:,:,1,1]


N_estimates = zeros(n_trees, n_heights, n_starting, n_models, 4)
η_estimates = zeros(n_trees, n_heights, n_starting, n_models, 4)


r0 = 0.04
ηtrue = r0 / 50
βs = [0.01, 0.02, 0.03, 0.04]
models = SSEconstant[]
for β in βs
    r = [
        r0,
        r0 + β,
        r0 + 2*β
    ]
    ϵ = 2/3
    λ = r ./ (1 - ϵ)
    μ = λ .- r
    model = SSEconstant(λ, μ, ηtrue)
    append!(models, [model])
end
models


params = collect(
    Iterators.product(
        1:n_trees,
        1:n_heights,
        1:n_starting,
        1:n_models
    )
)

#######################################
##
## Test true λ,μ and η
##
#######################################
io = open("output/prog.txt","w")
prog = ProgressMeter.Progress(n_trees * n_heights * n_starting * n_models; desc = "Inference 0: ", output = io)
#Threads.@threads for (i,j,k,m) in params
for (i,j,k,m) in params
    data = datasets[i,j,k,m]
    λtrue = models[m].λ
    μtrue = models[m].μ
    try
        rates = birth_death_shift(models[m], data; shift_bayes_factor = true)
        fpath = string("simulated_trees/results/true_all/s1_h", tree_heights[j], "_m", m, "_", i, ".tre")
        Pesto.writenewick(fpath, data, rates)
    catch

    end
    ProgressMeter.next!(prog)
end
close(io)

#######################################
##
## Test true λ,μ, and estimating η
##
#######################################
io = open("output/prog.txt","w")
prog = ProgressMeter.Progress(n_trees * n_heights * n_starting * n_models; desc = "Inference 1: ", output = io)
#Threads.@threads for (i,j,k,m) in params
for (i,j,k,m) in params
    data = datasets[i,j,k,m]
    λtrue = models[m].λ
    μtrue = models[m].μ
    try
	    η_estimates[i,j,k,m,1] = optimize_eta(λtrue, μtrue, data)
        model = SSEconstant(λtrue, μtrue, η_estimates[i,j,k,m,1])
        rates = birth_death_shift(model, data; shift_bayes_factor = true)
        fpath = string("simulated_trees/results/true_div/s1_h", tree_heights[j], "_m", m, "_", i, ".tre")
        Pesto.writenewick(fpath, data, rates)
    catch
        η_estimates[i,j,k,m,1] = -1
    end
    ProgressMeter.next!(prog)
end
close(io)


#######################################
##
## Test unknown λ,μ, and shift rate
##
#######################################
cbdp = zeros(n_trees, n_heights, n_starting, n_models, 2)
io = open("output/prog.txt","w")
prog = ProgressMeter.Progress(n_trees * n_heights * n_starting * n_models, desc = "Inference 2...", output = io)
#Threads.@threads for (i,j,k,m) in params
for (i,j,k,m) in params
    println(i,"\t", j, "\t", k, "\t", m)
    data = datasets[i,j,k,m]
    λml, μml = estimate_constant_bdp(data)
    cbdp[i,j,k,m,:] .= [λml, μml]

    H = 0.587
    n = 6
    dλ = LogNormal(log(λml), H)
    dμ = LogNormal(log(µml), H)
    λquantiles = make_quantiles(dλ, n)
    µquantiles = make_quantiles(dμ, n)
    λ, μ = allpairwise(λquantiles, µquantiles)
   
    try
        η_estimates[i,j,k,m,2] = optimize_eta(λ, μ, data)
        model = SSEconstant(λ, μ, η_estimates[i,j,k,m,2])
        rates = birth_death_shift(model, data; shift_bayes_factor = true)
        fpath = string("simulated_trees/results/unknown_rates/s1_h", tree_heights[j], "_m", m, "_", i, ".tre")
        Pesto.writenewick(fpath, data, rates)
    catch
        η_estimates[i,j,k,m,2] = -1
    end
    ProgressMeter.next!(prog)
end
close(io)

save("simulated_trees/eta_estimates.jld2", "eta", η_estimates)
save("simulated_trees/cbdp_estimates.jld2", "cbdp", cbdp)
