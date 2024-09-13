using Revise
using Pesto
using CairoMakie


#tree = readtree("simulated_trees/results/unknown_rates/s1_h125.0_m4_6.tre")
tree = readtree("simulated_trees/s1_h100.0_m4_22.tre")
ρ = 1.0

data = SSEdata(tree, ρ)


## the true rate categories
ϵ = 2/3
r_true = [0.04, 0.08, 0.12]

μ_true = r_true * ϵ / (1 - ϵ)
λ_true = μ_true .+ r_true


## mis-specified models
r_low = r_true .- 0.035
μ_low = r_low * ϵ / (1 - ϵ)
λ_low = μ_low .+ r_low

r_high = r_true .+ 0.05
μ_high = r_high * ϵ / (1 - ϵ)
λ_high = μ_high .+ r_high

## fit the models

η1 = optimize_eta(λ_true, μ_true, data)
η2 = optimize_eta(λ_low, μ_low, data)
η3 = optimize_eta(λ_high, μ_high, data)

model1 = SSEconstant(λ_true, μ_true, η1)
model2 = SSEconstant(λ_low, μ_low, η1)
model3 = SSEconstant(λ_high, μ_high, η1)
models = [model1, model2, model3]

df_rates = []
for model in models
    rates = birth_death_shift(model, data)
    push!(df_rates, rates)
end

[sum(rates[!,:nshift]) for rates in df_rates]

treeplot(data, df_rates[1])
treeplot(data, df_rates[2])
treeplot(data, df_rates[3])

7 / sum(data.branch_lengths)

writenewick("/tmp/true.tre", data, df_rates[1])
writenewick("/tmp/low.tre", data, df_rates[2])
writenewick("/tmp/high.tre", data, df_rates[3])


using RCall


R"""
library(ggtree)
library(treeio)
library(ggplot2)
library(patchwork)

tr1 <- read.beast.newick("/tmp/true.tre")
tr2 <- read.beast.newick("/tmp/low.tre")
tr3 <- read.beast.newick("/tmp/high.tre")
trees <- list(tr1, tr2, tr3)
#tr1@data$nshift <- sapply(tr1@data$N, sum)
rate_limits <- c(
    min(sapply(trees, function(tree) min(tree@data$mean_netdiv))),
    max(sapply(trees, function(tree) max(tree@data$mean_netdiv)))
    )
shift_limits <- c(
    min(sapply(trees, function(tree) min(tree@data$nshift))),
    max(sapply(trees, function(tree) max(tree@data$nshift)))
    )

p1a <- ggtree(tr1, aes(color = nshift)) +
    scale_colour_gradient("", low = "black", high = "red", limits = shift_limits) +
    ggtitle("a) true rate categories, r = [0.04,0.08,0.12]") +
    theme(legend.position = "none")
p2a <- ggtree(tr2, aes(color = nshift)) +
    scale_colour_gradient("", low = "black", high = "red", limits = shift_limits) +
    ggtitle("b) too low netdiv categories, r = [0.005,0.045,0.085]") +
    theme(legend.position = "none")
p3a <- ggtree(tr3, aes(color = nshift)) +
    scale_colour_gradient("no. rate shifts (N)", low = "black", high = "red", limits = shift_limits) +
    ggtitle("c) too high netdiv categories, r = [0.09,0.13,0.17]")

    
p1b <- ggtree(tr1, aes(color = mean_netdiv)) +
    scale_colour_gradient("", low = "gray", high = "#619CFF", limits = rate_limits) +
    theme(legend.position = "none")
p2b <- ggtree(tr2, aes(color = mean_netdiv)) +
    scale_colour_gradient("", low = "gray", high = "#619CFF", limits = rate_limits) +
    theme(legend.position = "none")
p3b <- ggtree(tr3, aes(color = mean_netdiv)) +
    scale_colour_gradient("mean netdiv", low = "gray", high = "#619CFF", limits = rate_limits)

p_shift <- p1a | p2a | p3a
p_rates <- p1b | p2b | p3b
p <- p_shift / p_rates
ggsave("figures/misspecied_rate_categories.pdf", p, width = 400, height = 315, units = "mm")
p
"""