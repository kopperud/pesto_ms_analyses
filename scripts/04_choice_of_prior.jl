## choice of prior distribution on the rates
using Distributions
using CairoMakie
using ProgressMeter
using Pesto
using Plots.PlotMeasures

primates_tree = readtree(Pesto.path("primates.tre"))
primates = SSEdata(primates_tree, 0.635)

H = 0.587405

λml, μml = estimate_constant_bdp(primates)
v = (exp(H^2) - 1)* exp(2*log(λml) + H^2) ## variance of the lognormal
θλ = v / λml
αλ = λml / θλ
v = (exp(H^2) - 1)* exp(2*log(μml) + H^2) ## variance of the lognormal
θμ = v / μml
αμ = μml / θμ


order = 1
aλ = λml * log(10^order)/(10^order - 1)
bλ = aλ * 10^order
aμ = μml * log(10^order)/(10^order-1)
bμ = aμ * 10^order

λdistributions = [
    LogNormal(log(λml), H),
    LogUniform(aλ, bλ),
    Exponential(λml),
    Gamma(αλ, θλ)
]
μdistributions = [
    LogNormal(log(μml), H),
    LogUniform(aμ, bμ),
    Exponential(μml),
    Gamma(αµ, θµ)
]

names_simple = [
    "LogNormal",
    "LogUniform",
    "Exponential",
    "Gamma"
]

names_latex = [
    L"\text{log normal}",
    L"\text{log uniform}",
    L"\text{exponential}",
    L"\text{gamma}"
]


n = 10

λml, μml = estimate_constant_bdp(primates)
v = (exp(H^2) - 1)* exp(2*log(λml) + H^2)
θλ = v / λml
αλ = λml / θλ
v = (exp(H^2) - 1)* exp(2*log(μml) + H^2)
θμ = v / μml
αμ = μml / θμ


order = 10
aλ = λml * log(order)/(order-1)
bλ = aλ * order
aμ = μml * log(order)/(order-1)
bμ = aμ * order

λdistributions = [
    LogNormal(log(λml), H),
    LogUniform(aλ, bλ),
    Exponential(λml),
    Gamma(αλ, θλ)
]

μdistributions = [
    LogNormal(log(μml), H),
    LogUniform(aμ, bμ),
    Exponential(μml),
    Gamma(αµ, θµ)
]

models = []
for (dλ, dμ) in zip(λdistributions, μdistributions)
    λquantiles = make_quantiles(dλ, n)
    µquantiles = make_quantiles(dμ, n)

    λ, μ = allpairwise(λquantiles, µquantiles)

    η = optimize_eta(λ, µ, primates)
    model = SSEconstant(λ, μ, η)
    append!(models, [model])
end


## Number of shifts
rates = Dict()
for (i, model) in enumerate(models)
    rates[i] = birth_death_shift(model, primates)
end
nshifts = [
    sum(rates[i][!,:nshift]) for i in 1:4
]
bfactors = [
    rates[i][!,:shift_bf] for i in 1:4
]
[maximum(x[1:end-1]) for x in bfactors]
strong_supported_bf = [
    sum(bfactors[i] .> 10) for i in 1:4
]

rates[4][rates[4][!,:shift_bf] .> 10, :]

fig = Figure(resolution=(650, 300), fontsize = 14,
            figure_padding = (5,15,5,5))

ax = Axis(
        fig[1,1],
        xlabel = L"\text{speciation rate}~(\lambda)", 
        ylabel = L"\text{density,}~\text{f}(\lambda)",
        #title = "\text{a) Priors",
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false,
        ygridvisible = false,
        title = "a)",
        titlealign = :left,
        xticks = [0.0, 0.25, 0.5, 0.75, 1.0],
        )
CairoMakie.ylims!(ax, 0.0, 8.0)
CairoMakie.xlims!(ax, 0.0, 1.0)


x = range(0.0, 2.0; length = 400)
colors = [:black, Makie.wong_colors()[1:3]...]
for (name, rdist, c) in zip(names_simple, λdistributions, colors)
    CairoMakie.lines!(ax, x, exp.(logpdf.(rdist, x)), label = name, linewidth = 3, color = c)
end
fig


##Loglikelihoods
logLs = zeros(length(λdistributions))
for (j, model) in enumerate(models)
    logLs[j] = logL_root(model, primates)
end
ΔlogLs = logLs .- (maximum(logLs))

ax1 = Axis(
    fig,
    bbox = CairoMakie.BBox(250, 400, 160, 280),
    rightspinevisible = false,
    bottomspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    title = "b)",
    titlealign = :left,
    ylabel = L"\Delta \log L",
    xticks = ([1,2,3,4], names_latex),
    xticklabelrotation=-π/4,
)
for (i, ΔlogL) in enumerate(ΔlogLs)
    CairoMakie.barplot!([i], [ΔlogLs[i]], color = colors[i])
end
fig

ax2 = Axis(
    fig,
    bbox = CairoMakie.BBox(450, 600, 160, 280),
    rightspinevisible = false,
    topspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
    title = "c)",
    titlealign = :left,
    ylabel = L"\text{number of shifts}~(\hat{N})",
    xticks = ([1,2,3,4], names_latex),
    xticklabelrotation=-π/4,
)


for (i, ΔlogL) in enumerate(ΔlogLs)
    CairoMakie.barplot!([i], [nshifts[i]], color = colors[i])
end

fig
#save("figures/choice-of-prior-makie-order100.pdf", fig)
save("figures/choice-of-prior.pdf", fig)