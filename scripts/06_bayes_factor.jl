using Pesto
using Distributions
using ProgressMeter
using Printf
using LaTeXStrings
using Measures
using CairoMakie

##############################
##
##   Load the data files
##
###############################
treefile = "data/primates.tre"
phy = readtree(treefile)
ρ = 0.635
primates = SSEdata(phy, ρ)

##############################
##
##   Set up the model
##
###############################
λml, μml = estimate_constant_bdp(primates)

H = 0.587405
d1 = LogNormal(log(λml), H)
d2 = LogNormal(log(μml), H)

n = 10
speciation = make_quantiles(d1, n)
extinction = make_quantiles(d2, n)

k = n ^2
λ, μ = allpairwise(speciation, extinction)
η = optimize_eta(λ, μ, primates) 
model = SSEconstant(λ, μ, η)

###############################
##
## Calculate probability of there being a shift
##
################################
posterior_prior_odds = posterior_prior_shift_odds(model, primates)

n_etas = 50
ηs = Pesto.lrange(0.0001, 0.7, n_etas)
models = [SSEconstant(λ, μ, rate) for rate in ηs]
n_edges = length(primates.branch_lengths)

ps = zeros(n_etas, n_edges)
prog = ProgressMeter.Progress(n_etas, "calculating...")
for i in 1:n_etas
    ps[i,:] = posterior_prior_shift_odds(models[i], primates; n_knots=50)
    next!(prog)
end

significance_level = 10.0
supported_shifts = sum(ps .> significance_level, dims = 2)[:,1]



fig = Figure(resolution=(850, 300), fontsize = 16,
figure_padding = (5,5,5,5))


ax1 = Axis(
        fig[1,1],
        xlabel = L"\text{branch index}", 
        ylabel = L"\text{Bayes factor (>0 shifts)}",
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false,
        ygridvisible = false,
        yscale = CairoMakie.log10,
        title = "a)",
        titlealign = :left,
        #xticks = [0.0, 0.25, 0.5, 0.75, 1.0],
        )

n_edges = length(primates.branch_lengths)
CairoMakie.plot!(ax1, 1:n_edges, posterior_prior_odds, label = "", color = :white, strokewidth = 1)

levels = [3.2, 10.0]
labels = [L"\text{substantial support}", L"\text{strong support}"]
linestyles = [:dash, :dashdot]
colors = [:gray, :black]
ls = []
for (level, label, c) in zip(levels, labels, colors)
    l = CairoMakie.lines!(ax1, [1,n_edges], [level, level], label = label, linestyle = :dash, 
            linewidth = 2, color = c)
    append!(ls, [l])
end
        

ax2 = Axis(
        fig[1,2],
        xlabel = L"\text{shift rate}~(\eta)", 
        ylabel = L"\text{no. supported branches}",
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xscale = CairoMakie.log10,
        title = "b)",
        titlealign = :left,
        yticks = [0,2,4,6,8,10],
        )

l3 = CairoMakie.lines!(ax2, ηs, supported_shifts,
    label = L"\text{alternative }\eta", 
    color = "black")
l4 = CairoMakie.plot!(ax2, [η], [sum(posterior_prior_odds .> significance_level)], 
        label = L"\text{MLE }\eta", color = :white, strokewidth = 1)

labels = [
    L"\text{strong support}",
    L"\text{substantial support}",
    L"\text{alternative }~\eta",
    L"\text{MLE}~\eta"]
Legend(fig[1, 3], [reverse(ls)..., l3, l4], labels, nbanks = 1, framevisible = false)
colsize!(fig.layout, 1, Relative(0.37))
colsize!(fig.layout, 2, Relative(0.37))
colgap!(fig.layout, 2)
fig
save("figures/shift_bayes_factor.pdf", fig)
