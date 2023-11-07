using DataFrames
using CSV
using Statistics
using Pesto
using CairoMakie
using LaTeXStrings

df = CSV.read("output/primates_LSBDS_rates.log", DataFrame)

nrows = size(df)[1]

n_edges = 465
rb_rates = zeros(nrows, n_edges)
rb_shifts = zeros(nrows, n_edges)

for branch_index in 1:n_edges
    rb_rates[:,branch_index] = df[!, "avg_lambda[$branch_index]"]
    rb_shifts[:,branch_index] = df[!, "num_shifts[$branch_index]"]
end


mean_lambda = Statistics.mean(rb_rates, dims = 1)[1,:]
mean_shifts = Statistics.mean(rb_shifts, dims = 1)[1,:]
mean_lambda[end] = NaN
mean_shifts[end] = NaN

## Pesto Inference
treefile = "data/primates.tre"
phy = Pesto.readtree(treefile)
num_total_species = 367
ρ = length(phy.tip_label) / num_total_species
primates = SSEdata(phy, ρ)

using Pesto
λ = [0.1, 0.2, 0.4, 0.1]
μ = [0.05, 0.15, 0.10, 0.2]
η = optimize_eta(λ, μ, primates)

model = SSEconstant(λ, μ, η)
pesto_rates = birth_death_shift(model, primates)
sort!(pesto_rates, :node) ## sort by node index

R"""
source("scripts/matchnodes.R")
mn = matchNodes(phy)
""";
@rget mn
mn[!, :Rev] = [Int64(x) for x in mn[!,:Rev]]

rev_to_ape = Dict()
ape_to_rev = Dict()
for row in eachrow(mn)
    rev_to_ape[row[:Rev]] = row[:R]
    ape_to_rev[row[:R]] = row[:Rev]
end

rev_to_ape

new_df = DataFrame(
    "pesto_lambda" => pesto_rates[!,:mean_lambda],
    "pesto_nshifts" => pesto_rates[!,:nshift],
    "revbayes_lambda" => [
        mean_lambda[ape_to_rev[i]] for i in 1:n_edges
    ],
    "revbayes_nshifts" => [
        mean_shifts[ape_to_rev[i]] for i in 1:n_edges
    ]
)


fig = Figure(resolution=(550, 250), fontsize = 14,
            figure_padding = (0,2,5,0))

ax1 = Axis(fig[1,1],
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false, 
        ygridvisible = false,
        #xscale = Makie.log10,
        #yscale = Makie.log10,
        title = "a) mean speciation rates",
        titlealign = :left,
        ylabel = L"\text{Pesto}",
        xlabel = L"\text{RevBayes}")

ax2 = Axis(fig[1,2],
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false, 
        ygridvisible = false,
        #xscale = Makie.log10,
        #yscale = Makie.log10,
        title = "b) number of shifts",
        titlealign = :left,
        #ylabel = L"\text{Pesto (deterministic alg.)}",
        xlabel = L"\text{RevBayes}")

CairoMakie.plot!(
    ax1,
    new_df[!,:pesto_lambda],
    new_df[!,:revbayes_lambda],
    label = "branch rate",
    strokewidth = 1,
    color = :white)
CairoMakie.lines!(
    ax1,
    [0.18, 0.4],
    [0.18, 0.4],
    label = "one-to-one",
    linestyle = :dash,
    color = :gray
)

CairoMakie.plot!(
    ax2, 
    new_df[!,:pesto_nshifts],
    new_df[!,:revbayes_nshifts],
    label = "branch rate",
    strokewidth = 1,
    color = :white)
CairoMakie.lines!(
        ax2,
        [0.0, 1.0],
        [0.0, 1.0],
        label = "one-to-one",
        linestyle = :dash,
        color = :gray
    )
fig
CairoMakie.save("figures/rb_vs_pesto.pdf", fig)


## plot tree for fun
using RCall

@rput primates
@rput pesto_rates
R"""
library(ape)
library(ggtree)
library(tidytree)
library(ggplot2)
library(dplyr)

x <- as_tibble(primates)
phydf <- merge(x, pesto_rates, by = "node")
th <- max(node.depth.edgelength(primates))

td_phy <- as.treedata(phydf)

td_phy@phylo$tip.label <- gsub("_", " ", td_phy@phylo$tip.label)

scalexformat <- function(x) sprintf("%.0f Ma", th - abs(round(x, 0)))

bgcolor="white"
fgcolor="black"

p1 <- ggtree(td_phy, aes(color = nshift), size = 1) +
    theme(legend.position = c(0.15, 0.8)) +
    scale_colour_gradient(low = "black", high = "red", name = "a) Rate shifts") +
    scale_x_continuous(labels = scalexformat, breaks=seq(0,th, length.out = 5), 
                        limits = c(0, th+3))
                       

p2 <- ggtree(td_phy, aes(color = `mean_lambda`), size = 1) +
    theme(legend.position = c(0.15, 0.8),
    legend.background = element_blank()) +
    scale_colour_gradient(low = "black", high = "#56B1F7", name = "b) Speciation rate") +
    ggtree::geom_cladelab(node = 386, label = "Cercopithecidae", 
                          offset = 0.5, offset.text = 1, angle = 90, hjust = 0.5)  +
    scale_x_continuous(labels = scalexformat, breaks=seq(0,th, length.out = 5), 
                        limits = c(0, th+4)) 

pc <- (p1 | p2) &
        theme(
    axis.line.x=element_line(color=fgcolor),
    axis.text.x=element_text(size=12),
    legend.text = element_text(size=15),
    legend.title = element_text(size=15),
    legend.background = element_blank()
        )

#ggsave("figures/primates_shift_lambda.pdf", pc, width = 200, height = 120, units = "mm")
"""