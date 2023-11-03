using Distributions
using RCall
using StatsPlots
using Pesto
using Optim
using LaTeXStrings

##############################
##
##   Load the data files
##
###############################
treefile = "data/primates.tre"
phy = Pesto.readtree(treefile)
num_total_species = 367
ρ = length(phy.tip_label) / num_total_species
primates = SSEdata(phy, ρ)

##############################
##
##   Set up the model
##
###############################

H = 0.587405

λml, μml = Pesto.estimate_constant_bdp(primates)
d1 = LogNormal(log(λml), H)
d2 = LogNormal(log(μml), H)

n = 10
speciation = make_quantiles(d1, n)
extinction = make_quantiles(d2, n)

k = n ^2
λ, μ = allpairwise(speciation, extinction)

η = optimize_eta(λ, μ, primates)
println(λml, "\t", μml, "\t", η)

model = SSEconstant(λ, μ, η)

rates = birth_death_shift(model, primates)

##############################
##
## Plot the tree using some R code
##
###############################

@rput primates
@rput rates
R"""
library(ape)
library(ggtree)
library(tidytree)
library(ggplot2)
library(dplyr)

x <- as_tibble(primates)
phydf <- merge(x, rates, by = "node")
th <- max(node.depth.edgelength(primates))

td_phy <- as.treedata(phydf)

td_phy@phylo$tip.label <- gsub("_", " ", td_phy@phylo$tip.label)

scalexformat <- function(x) sprintf("%.0f Ma", th - abs(round(x, 0)))

bgcolor="white"
fgcolor="black"

p1 <- ggtree(td_phy, aes(color = nshift), size = 1) +
    #geom_tiplab(size = 2) +
    theme(legend.position = c(0.15, 0.8)) +
    scale_colour_gradient(low = "black", high = "red", name = "a) Rate shifts") +
    #ggtree::geom_cladelab(node = 386, label = "Cercopithecidae", 
    #                      offset = 0.5, offset.text = 1, angle = 90, hjust = 0.5)  +
    scale_x_continuous(labels = scalexformat, breaks=seq(0,th, length.out = 5), 
                        limits = c(0, th+3))
                       

p2 <- ggtree(td_phy, aes(color = `mean_lambda`), size = 1) +
    #geom_tiplab(size = 3, fontface = 3) +
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

ggsave("figures/primates_shift_lambda.pdf", pc, width = 200, height = 120, units = "mm")
""";
