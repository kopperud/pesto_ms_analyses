## this script calculates the probability densities
## for the schematic in the ms (Fig. 1), and prints them

# libraries
using Plots
using CairoMakie
using RCall
using Pesto

R"""       
library(ape)
treestring <- "((A:0.1,B:0.1):0.1,C:0.2):0.0;"
phy <- read.tree(text = treestring)

nde <- node.depth.edgelength(phy)
node_depths <- max(nde) - nde
phy$node_depths <- node_depths
phy$branching_times <- branching.times(phy)

po <- postorder(phy)
phy$po <- po
"""
@rget phy
tree = Pesto.phylo(phy[:edge],
    phy[:edge_length],
    phy[:Nnode],
    phy[:tip_label],
    phy[:node_depths],
    phy[:branching_times],
    phy[:po])
ρ = 1.0
faketree = SSEdata(tree, ρ)

## set up a two-state model
lambda_vec = [1.2, 0.5]
mu_vec = [1.0, 0.3]
η = 0.01
model = SSEconstant(
    lambda_vec,
    mu_vec,
    0.01
)

## Need to hack the source code in Pesto to remove re-scaling before 
## running this, otherwise it will give the re-scaled densities.
Ds, Fs = backwards_forwards_pass(model, faketree);
Ss = ancestral_state_probabilities(faketree, Ds, Fs);

for (i, D) in Ds
    anc, dec = tree.edge[i,:]
    if dec > 3
        label = "internal"
    else
        label = dec
    end
    a, b = extrema(D.t)
    r = D(a)
    println("D[$label](a) = $r")
    r = D(b)
    println("D[$label](b) = $r")
    #println(b)
end

for (i, F) in Fs
    anc, dec = tree.edge[i,:]
    if dec > 3
        label = "internal"
    else
        label = dec
    end
    a, b = extrema(F.t)
    r = F(a)
    println("F[$label](a) = $r")
    r = F(b)
    println("F[$label](b) = $r")
    #println(b)
end


Fs[4](0.2) .* Ds[4](0.2)





















