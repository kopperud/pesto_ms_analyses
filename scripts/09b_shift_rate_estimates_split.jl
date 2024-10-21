using Revise
using Pesto
using JLD2
using Glob
using ProgressMeter
using Distributions
#using StatsPlots
using CairoMakie
using LaTeXStrings


dat = load("simulated_trees/results/N.jld2")
N = dat["N"]

dat2 = load("simulated_trees/results/eta_estimates.jld2")
η_estimates = dat2["eta"]

cbdp = load("simulated_trees/results/cbdp_estimates.jld2")["cbdp"]

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


η_estimates[:,4,:,4,1]
size(η_estimates)
βs = [0.01, 0.02, 0.03, 0.04]
ηtrue = 0.0008


## Some trees where it crashed
sum(η_estimates .== -1)
argmax(η_estimates .== -1)


fig = Figure(size=(650, 500), fontsize = 14,
            figure_padding = (0,0,0,0))
axs = []

titles = [
    L"\text{tiny var.}",
    L"\text{small var.}",
    L"\text{moderate var.}",
    L"\text{large var.}"
]

#xt = [10.0, 100.0, 1_000.0, 10_000.0]

xlabel = Label(fig[1, 2:5], L"\text{a) Simulated trees without rate shifts}", justification = :left)
xlabel = Label(fig[4, 2:5], L"\text{b) Simulated trees with }\geq \text{1 rate shifts}", justification = :left)

for q in 1:2
    for m in 1:n_models
        if q == 1
            title = titles[m]
        else
            title = ""
        end
        ax = Axis(fig[q+1,m+1], 
                  yscale = Makie.log10,
                  xscale = Makie.log10,
                  topspinevisible = false,
                  xgridvisible = false, 
                  ygridvisible = false,
                  rightspinevisible = false,
                  #title = title,
                  #xticks = xt,
                  #xlabel = L"\text{number of taxa}", 
                  )
        if m > 1
            hideydecorations!(ax, ticks = false)
        end
        if q == 1
            hidexdecorations!(ax, ticks = false)
        end

        for j in 1:5
            y = η_estimates[:,j,1,m,q]
            is_pos = y .> 0
            no_shifts = Nsum[:,j,1,m] .== 0
            x = taxa[:,j,1,m]
           
            # filter out those that i) crashed and ii) those that did not experience a shift
            y = y[is_pos .& no_shifts]
            x = x[is_pos .& no_shifts]

            CairoMakie.scatter!(ax, x, y, 
                                label = LaTeXString(string(Int64(tree_heights[j])) * raw" Ma"),
                                markersize = 7,
                                strokewidth = 0.5)
        end
        CairoMakie.lines!(ax, [extrema(taxa[:,:,1,m])...], 
                        [ηtrue, ηtrue], label = L"\text{true}~\eta",
                         color = :red, linewidth = 2, linestyle = :dash)
        append!(axs, [ax])
    end
end

#############################
##
## for those that experienced >= 1 rate shifts
##
#######################

for q in 1:2
    for m in 1:n_models
        if q == 1
            title = titles[m]
        else
            title = ""
        end
        ax = Axis(fig[q+4,m+1], 
                  yscale = Makie.log10,
                  xscale = Makie.log10,
                  topspinevisible = false,
                  xgridvisible = false, 
                  ygridvisible = false,
                  rightspinevisible = false,
                  title = title,
                  #xticks = xt,
                  #xlabel = L"\text{number of taxa}", 
                )
        if m > 1
            hideydecorations!(ax, ticks = false)
        end
        if q == 1
            hidexdecorations!(ax, ticks = false)
        end

        for j in 1:5
            y = η_estimates[:,j,1,m,q]
            is_pos = y .> 0
            minimum_one_shift = Nsum[:,j,1,m] .> 0
            x = taxa[:,j,1,m]
           
            # filter out those that i) crashed and ii) those that did not experience a shift
            y = y[is_pos .& minimum_one_shift]
            x = x[is_pos .& minimum_one_shift]
            
            CairoMakie.scatter!(ax, x, y, 
                                label = LaTeXString(string(Int64(tree_heights[j])) * raw" Ma"),
                                markersize = 7,
                                strokewidth = 0.5)
        end
        CairoMakie.lines!(ax, [extrema(taxa[:,:,1,m])...], 
                        [ηtrue, ηtrue], label = L"\text{true}~\eta",
                         color = :red, linewidth = 2, linestyle = :dash)
        append!(axs, [ax])
    end
end
#CairoMakie.legend(fig[1:2,5])
#ylabel = Label(fig[1, 6], L"\text{true}~𝛌,𝛍,\eta", rotation = π/2)
## fake LABELS
xlabel1 = Label(fig[7, 2:5], L"\text{number of taxa}", justification = :center)

ylabel = Label(fig[2:3, 1], L"\text{shift rate~}(\eta)", rotation = π/2)
ylabel = Label(fig[2, 6], L"\text{true}~𝛌,𝛍", rotation = π/2)
ylabel = Label(fig[3, 6], L"\text{unknown}~𝛌,𝛍", rotation = π/2)

ylabel = Label(fig[5:6, 1], L"\text{shift rate~}(\eta)", rotation = π/2)
ylabel = Label(fig[5, 6], L"\text{true}~𝛌,𝛍", rotation = π/2)
ylabel = Label(fig[6, 6], L"\text{unknown}~𝛌,𝛍", rotation = π/2)

fig[3:5, 7] = Legend(fig, axs[1], "", framevisible = false, 
                patchsize = (30, 30))


for row in [2,3,5,6]
    rowsize!(fig.layout, row, Relative(0.19))
end
for row in [1,4,7]
    rowsize!(fig.layout, row, Relative(0.08))
end


rowgap!(fig.layout, 2.0)
colgap!(fig.layout, 2.0)
linkaxes!(axs...)
fig

CairoMakie.save("figures/eta-estimates-sim-split.pdf", fig)

