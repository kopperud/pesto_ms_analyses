using Revise
using BirthDeathSimulation
using ProgressMeter
using JLD2
using Distributions
using Random

Random.seed!(1234)

βs = [0.01, 0.02, 0.03, 0.04]
r0 = 0.04
models = bdsmodel[]
for β in βs
    r = [
        r0,
        r0 + β,
        r0 + 2*β
    ]
    ϵ = 2/3
    λ = r ./ (1 - ϵ)
    μ = λ .- r
    η = r0 / 50 ## fifty branching events per rate shift events
    model = bdsmodel(λ, μ, η)
    append!(models, [model])
end
models


maxtaxa = 50000
starting = [1]
n_trees = 100
n_starting = 1#length(λ)

tree_heights = [25.0, 50.0, 75.0, 100.0, 125.0]
n_heights = length(tree_heights)
n_models = length(models)

trees = Array{Tree, 4}(undef, n_trees, n_heights, n_starting, n_models)
N = Array{Float64, 6}(undef, n_trees, n_heights, n_starting, n_models, length(models[1].λ), length(models[1].λ))
for j in 1:n_heights
    h = tree_heights[j]
    for m in 1:n_models
        model = models[m]

        for k in 1:n_starting
            state = 1
            prog = ProgressUnknown("Complete trees sim (tree height = $h, state = $state, model = $m):")
            i = 1
            while i <= n_trees
                maxtime = tree_heights[j]
                tree = sim_bdshift(model, maxtime, maxtaxa, state)
                if length(tree.Leaves) < maxtaxa ## reject complete trees that termined when too many taxa
                    if length(tree.Leaves) > 1 ## Reject completete trees where all taxa went extinct
                        prune_extinct!(tree)
                        if abs(treeheight(tree) - maxtime) < 0.001 ## reject complete trees where both root children did not survive
                            N0 = +([branch.N for (idx, branch) in tree.Branches]...)
                            N[i,j,k,m,:,:] .= N0
                            trees[i,j,k,m] = tree
                            i += 1
                        end
                    end
                end            
                ProgressMeter.next!(prog)
            end
            sleep(0.1)
            ProgressMeter.finish!(prog)
        end
    end
end

taxa = zeros(Int64, n_trees, n_heights, n_starting, n_models)
treelengths = zeros(n_trees, n_heights, n_starting, n_models)
for i in 1:n_trees
    for j in 1:n_heights
        for k in 1:n_starting
            for m in 1:n_models
                tree = trees[i,j,k,m]
                taxa[i,j,k,m] = PhylogeneticTrees.ntaxa(tree)
                treelengths[i,j,k,m] = PhylogeneticTrees.treelength(tree)
            end
        end
    end
end
sum(taxa, dims = 1)[1,:,:,:] ./ n_trees


prog = ProgressMeter.Progress(n_trees * n_heights * n_starting * n_models, "Writing trees...")
for i in 1:n_trees
    for j in 1:n_heights
        for k in 1:n_starting
            for m in 1:n_models
                tree = trees[i,j,k,m]
                h = tree_heights[j]
                state = collect(1:n_starting)[k]
                fpath = string("simulated_trees/s", state, "_h", h, "_m", m, "_", i, ".tre")
                writenewick(fpath, tree, models[m])
                ProgressMeter.next!(prog)
            end
        end
    end
end

save("simulated_trees/results/N.jld2", "N", N)