using Random
using Distributions
using BirthDeathSimulation
using StatsPlots

Random.seed!(1234)

#λ = [0.3, 0.5] ## speciation rates
#µ = [0.05, 0.15] ## extinction rates
λ = [0.2, 0.4]
μ = [0.05, 0.15]
η = 0.05 ## shift rate
model = bdsmodel(λ, µ, η)



max_time = 25.0
max_taxa = 10_000

trees = []
for i in 1:2500
    starting_state = 1
    tree = sim_bdshift(model, max_time, max_taxa, starting_state)

    if ntaxa(tree) < max_taxa
        append!(trees, [tree])
    end    
end

event_rates = Float64[]
event_probs = Float64[]
for tree in trees
    sumN = sum([sum(branch.N) for (i, branch) in tree.Branches])
    tree_length = sum(branch.bl for (i, branch) in tree.Branches)
    event_rate = sumN / tree_length
    append!(event_rates, event_rate)

    d = Poisson(η * tree_length)
    event_prob = exp(logpdf(d, sumN))
    append!(event_probs, event_prob)
end

p1 = histogram(event_rates, bins = 55, grid = false, 
            label = "N/treelength", title = "complete trees",
            xlabel = "events per time", ylabel = "frequency")
plot!(p1, [η, η], [0.0, 200], linestyle = :dash, color = :black, linewidth = 3, label = "η")

#### now reconstructed trees 

trees2 = []
event_rates2 = Float64[]
#for i in 1:500
while length(trees2) < 2500
    starting_state = 1
    tree = sim_bdshift(model, max_time, max_taxa, starting_state)

    if ntaxa(tree) < max_taxa
        if ntaxa(tree) > 1 ## reject small trees
            prune_extinct!(tree)
            if abs(treeheight(tree) - max_time) < 0.001 ## reject complete trees where both root children did not survive
                append!(trees2, [tree])

                sumN = sum([sum(branch.N) for (i, branch) in tree.Branches])
                tree_length = sum(branch.bl for (i, branch) in tree.Branches)
                event_rate = sumN / tree_length
                append!(event_rates2, event_rate)
            end
        end
    end    
end

length(trees2)

p2 = histogram(event_rates2, bins = 12, grid = false, 
    label = "N/treelength", title = "reconstructed trees",
    xlabel = "events per time", ylabel = "frequency")
plot!(p2, [η, η], [0.0, 100], linestyle = :dash, color = :black, linewidth = 3, label = "η")

plot(p1, p2, xlim = (0.0, 0.2))

p3 = plot(grid = false, xlabel = "events per time", ylabel = "density")
StatsPlots.density!(p3, event_rates, label = "complete trees")
StatsPlots.density!(p3, event_rates2, label = "reconstructed trees")
d = Poisson(η)
logpdf.(d, collect(range(0, 1.0; length=100)))
StatsPlots.plot!(p3, )

plot(p3, xlim = (0.0, 0.12))
savefig(p3, "figures/N_prior_distribution.pdf")


mean(event_rates2)
mean(event_rates)

using Statistics
using HypothesisTests

ApproximateTwoSampleKSTest(event_rates, event_rates2)








