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

#################################
##
##  helper function
##
############################

function foobaz(data, η)
    n = 10
    λml, μml = estimate_constant_bdp(data)
    H = 0.587
    dλ = LogNormal(log(λml), H)
    dμ = LogNormal(log(μml), H)

    λquantiles = make_quantiles(dλ, n)
    µquantiles = make_quantiles(dμ, n)

    λ, μ = allpairwise(λquantiles, µquantiles)
    model = SSEconstant(λ, μ, η)
    N = compute_nshifts(model, data)
    Nsum = sum(N)
    #r = rates[!, :mean_lambda]
    #r = r[.! isnan.(r)]
    return(Nsum)
end

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

ηtl = Pesto.lrange(0.05, 100.0, 9) |> collect
tl = sum(primates.branch_lengths)
Nhats = [foobaz(primates, ηt/tl) for ηt in ηtl]

n_etas = 50
ηs = Pesto.lrange(0.05, 100.0, n_etas) ./ tl
models = [SSEconstant(λ, μ, rate) for rate in ηs]
n_edges = length(primates.branch_lengths)

ps = zeros(n_etas, n_edges)
prog = ProgressMeter.Progress(n_etas, "calculating...")
for i in 1:n_etas
    ps[i,:] = posterior_prior_shift_odds(models[i], primates)
    next!(prog)
end

significance_level = 10.0
supported_shifts = sum(ps .> significance_level, dims = 2)[:,1]

## logL surface
logLs = zeros(50)
@showprogress for (i, eta) in enumerate(ηs)
    m = SSEconstant(λ, μ, eta)
    logLs[i] = logL_root(m, primates)
end



fig = Figure(resolution=(850, 350), fontsize = 16,
figure_padding = (5,5,5,5))


ax1 = Axis(
        fig[1,1],
        #xlabel = L"\text{branch index}", 
        ylabel = L"\text{Bayes factor (}\geq 1~\text{shifts)}",
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false,
        ygridvisible = false,
        yscale = CairoMakie.log10,
        title = "a)",
        titlealign = :left,
        #xticks = [0.0, 0.25, 0.5, 0.75, 1.0],
        yticks = [0.3, 1.0, 3.2, 10.0]
        )
ylims!(ax1, minimum(posterior_prior_odds).*0.5, maximum(posterior_prior_odds)*10)

n_edges = length(primates.branch_lengths)
CairoMakie.scatter!(ax1, 1:n_edges, posterior_prior_odds, 
                label = "", color = :black, 
                markersize = 6,
                strokewidth = 1)

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
CairoMakie.axislegend(ax1, reverse(ls), reverse(labels), position = :lt, framevisible = false)#; merge = true, unique = true, position = :lt)
#Legend(fig[1,1], [reverse(ls)...], labels, nbanks = 1, framevisible = false)
       

ax2 = Axis(
        fig[1,2],
        #xlabel = L"\text{shift rate}~(\eta)", 
        ylabel = L"\text{no. supported branches}",
        #xlabel = L"\text{expected number of shifts}~(E[N] = \eta t)",
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xscale = CairoMakie.log10,
        xticklabelrotation = π/2,
        title = "b)",
        titlealign = :left,
        yticks = [0,2,4,6,8,10],
        xticks = round.(ηtl; digits = 2),
        )

l3 = CairoMakie.lines!(ax2, ηs .* tl, supported_shifts,
    label = L"\text{alternative }\eta", 
    color = "black")
l4 = CairoMakie.plot!(ax2, [η * tl], [sum(posterior_prior_odds .> significance_level)], 
        label = L"\text{MLE }\eta", color = :white, strokewidth = 1)
CairoMakie.axislegend(ax2, merge = true, unique = true, position = :rt, framevisible = false)


ax3 = Axis(fig[1,3],
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false, 
        ygridvisible = false,
        xscale = Makie.log10,
        yscale = Makie.log10,
        title = "c)",
        titlealign = :left,
        xticklabelrotation = π/2,
        xticks = round.(ηtl; digits = 2),
        yticks = round.(ηtl; digits = 2),
        ylabel = L"\text{posterior number of shifts}~(\hat{N})",
        #xlabel = L"\text{expected number of shifts}~(E[N] = \eta t)"
        )
#CairoMakie.scatter!(ax3, ηtl, Nhats, label = "shifts", color = :black, markerstrokewidth = 1)
CairoMakie.lines!(ax3, [1.0, 100.0], [1.0, 100.0], label = L"\text{identity line}", color = :black, linestyle = :dash)
CairoMakie.lines!(ax3, ηtl, Nhats, label = L"\text{alternative}~\eta", color = :black)
CairoMakie.plot!(ax3, [η*tl], [foobaz(primates, η)], label = L"\text{MLE}~\eta", 
                color = :white, strokewidth = 1)
CairoMakie.axislegend(ax3, merge = true, unique = true, position = :lt, framevisible = false)

#= 
ax4 = Axis(fig[1,4],
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false, 
        ygridvisible = false,
        xscale = Makie.log10,
        #yscale = Makie.log10,
        title = "c)",
        titlealign = :left,
        xticklabelrotation = π/2,
        #yticks = round.(ηtl; digits = 2),
        xticks = round.(ηtl; digits = 2),
        #yticks = round.(ηtl; digits = 2),
        ylabel = L"\text{estimated number of shifts}~(\hat{N})",
        )
CairoMakie.lines!(ax4, ηs .* tl, logLs, label = L"\text{alternative}~\eta", color = :black, linestyle = :dash)
 =#
xlabel1 = Label(fig[2, 1], L"\text{branch index}")
xlabel2 = Label(fig[2, 2:3], L"\text{prior number of shifts}~(E[N] = \eta t)")

for i in 1:3
    colsize!(fig.layout, i, Relative(0.33))
end
colgap!(fig.layout, 2)
rowgap!(fig.layout, 0)
fig

fig
save("figures/shift_bayes_factor.pdf", fig)

fig