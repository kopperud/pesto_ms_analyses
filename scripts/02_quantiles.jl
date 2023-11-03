## plot the quantiles
using Distributions
using LaTeXStrings
using CairoMakie
import Pesto

H = 0.587
dλ = LogNormal(log(0.30), H)
dμ = LogNormal(log(0.22), H)
n = 6

fig = Figure(resolution = (320, 220), 
    fontsize = 14,
    figure_padding = (17.0, 17.0, 2.0, 0.0))

## TOP PANEL
xs = collect(range(0.0, 0.8; length = 100))
λquantiles = Pesto.make_quantiles(dλ, n)
ax1 = Axis(fig[1,1],
        xticks = (λquantiles, [L"\lambda_1", L"\lambda_2", L"\lambda_3", L"\lambda_4", L"\lambda_5", L"\lambda_6"]),
        xticklabelsize = 12)
for q in λquantiles
    lines!(ax1, [q, q], [0.0, exp(loglikelihood(dλ, q))], label = "", color = "gray", alpha = 1)
end
lines!(ax1, xs, exp.(loglikelihood.(dλ, xs)), color = "black")

### RIGHT PANEL
ys = collect(range(0.0, 0.6; length = 100))
μquantiles = Pesto.make_quantiles(dμ, n)
ax2 = Axis(fig[2,2],
        xgridvisible = false, 
        ygridvisible = false,
        yticks = (μquantiles, [L"\mu_1", L"\mu_2", L"\mu_3", L"\mu_4", L"\mu_5", L"\mu_6"]),
        xticklabelrotation = 0,
        yticklabelsize = 12)

for q in μquantiles
    lines!(ax2, [0.0, exp(loglikelihood(dμ, q))], [q, q], label = "", color = "gray", alpha = 1)
end
lines!(ax2, exp.(logpdf.(dμ,ys)), ys, color = "black")

## BOTTOM PANEL
ax3 = Axis(fig[2,1],
        xgridvisible = false, 
        ygridvisible = false,
        xlabel = L"\text{speciation rate }(\lambda)",
        ylabel = L"\text{extinction rate }(\mu)",
        topspinevisible = false,
        rightspinevisible = false,
        xticklabelrotation = π/2,
        xticklabelsize = 9,
        yticklabelsize = 9)

λs, μs = Pesto.allpairwise(λquantiles, μquantiles)
CairoMakie.plot!(ax3, λs, μs, color = "black", markersize = 8)
CairoMakie.plot!(ax3, [0.30], [0.22], color = "red", 
                label = "MLE simple model",
                markersize = 8)

elements = [MarkerElement(marker = :circle, color = :black),
            MarkerElement(marker = :circle, color = :red)]
labels = [L"\text{States}", L"\text{Simple model}"]
leg = Legend(fig[1,2], elements, labels, labelsize=9, 
                #patchsize = (10.0f0, 10.0f0),
                framevisible = false)

#leg.tellheight = true
#leg.tellwidth = true

linkxaxes!(ax1, ax3)
linkyaxes!(ax2, ax3)

hideydecorations!(ax1)
hidexdecorations!(ax2)
hideydecorations!(ax3, label = false, ticks = false, ticklabels = false)
hidexdecorations!(ax3, label = false, ticks = false, ticklabels = false)
hidespines!(ax1)
hidespines!(ax2)
hidexdecorations!(ax1, ticklabels = false)
hideydecorations!(ax2, ticklabels = false)

colsize!(fig.layout, 1, Relative(0.7))
rowsize!(fig.layout, 2, Relative(0.7))
rowgap!(fig.layout, -8.0)
colgap!(fig.layout, -8.0)

fig
CairoMakie.save("figures/model_makie.pdf", fig)