## How much does the number of shifts depend on the η parameter?
## How much does the number of shifts depend on the size of the state space?

using Distributions
using Pesto
using ProgressMeter
using CairoMakie

phy = readtree(Pesto.path("primates.tre"))
ρ = 0.635
primates = SSEdata(phy, ρ)

function foobar(data, n)
    λml, μml = estimate_constant_bdp(data)
    H = 0.587
    dλ = LogNormal(log(λml), H)
    dμ = LogNormal(log(μml), H)

    λquantiles = make_quantiles(dλ, n)
    µquantiles = make_quantiles(dμ, n)

    λ, μ = allpairwise(λquantiles, µquantiles)
    η = optimize_eta(λ, μ, data)
    model = SSEconstant(λ, μ, η)

    Ds, Fs = backwards_forwards_pass(model, data);

    N = compute_nshifts(model, data, Ds, Fs)
    Nsum = sum(N)
    return(Nsum)
end


ns = [2, 3, 4, 5, 6, 7, 8, 9, 10]
Ns = zeros(length(ns))

@showprogress for (i, n) in enumerate(ns)
    Ns[i] = foobar(primates, n)
end


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
    Ds, Fs = backwards_forwards_pass(model, data);

    N = compute_nshifts(model, data, Ds, Fs)
    Nsum = sum(N)
    #r = rates[!, :mean_lambda]
    #r = r[.! isnan.(r)]
    return(Nsum)
end

function branchrates(data, n)
    λml, μml = estimate_constant_bdp(data)
    H = 0.587
    dλ = LogNormal(log(λml), H)
    dμ = LogNormal(log(μml), H)

    λquantiles = make_quantiles(dλ, n)
    µquantiles = make_quantiles(dμ, n)

    λ, μ = allpairwise(λquantiles, µquantiles)
    ηml = optimize_eta(λ, μ, data)
    model = SSEconstant(λ, μ, ηml)
    rates = birth_death_shift(model, data)
    r = rates[!, :mean_netdiv]
    return(r[1:end-1])
end

ηtl = Pesto.lrange(0.05, 100.0, 9) |> collect
tl = sum(primates.branch_lengths)
Nhats = [foobaz(primates, ηt/tl) for ηt in ηtl]

branch_rates = zeros(length(ns), 464)
@showprogress for (i, n) in enumerate(ns)
    branch_rates[i,:] .= branchrates(primates, n)
end

fig = Figure(size=(550, 300), fontsize = 14,
            figure_padding = (0,0,0,0))

ax1 = Axis(fig[1,1],
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false,
        ygridvisible = false,
        title = "a)",
        titlealign = :left,
        xscale = Makie.sqrt,
        xticklabelrotation = π/2,
        xticks = ns.^2,
        ylabel = L"\text{posterior mean no. shifts}~(\hat{N})",
        #xlabel = L"\text{no. rate categories}~(n^2 = K)"
        )
CairoMakie.ylims!(ax1, 0.0, 6.2)
CairoMakie.scatter!(ax1, ns.^2, Ns[:,1], label = "primates", color = :black, markerstrokewidth = 1)
CairoMakie.lines!(ax1, ns.^2, Ns[:,1], label = "", color = :black)

ax2 = Axis(fig[1,2],
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false, 
        ygridvisible = false,
        title = "b)",
        titlealign = :left,
        #xscale = Makie.sqrt,
        #yscale = Makie.log10,
        xticklabelrotation = π/2,
        xticks = (ns, string.(ns.^2)),
#        xlabel = L"\text{no. rate categories}~(n^2 = K)",
        ylabel = L"\text{net diversification rate}~(\bar{r})")

xlabel = Label(fig[2,1:2], L"\text{no. rate categories}~(n^2 = K)")
rowgap!(fig.layout, 0)

colors = ["black", "gray", "black", "gray","black", "gray","black", "gray", "black"]
#colors = repeat([:steelblue], 9)
xs = vcat([repeat([i], 464) for i in ns]...)
cs = vcat([repeat([i], 464) for i in colors]...)

CairoMakie.rainclouds!(ax2, xs, vcat(branch_rates'...), 
            color = cs, label = "Posterior",
                clouds=nothing,
                plot_boxplots = false,
                jitter_width=0.3,
                side_nudge = 0.0,
                markersize=4)
fig

for i in 1:2
    colsize!(fig.layout, i, Relative(1/2))
end
fig

#CairoMakie.save("figures/state_and_eta_sensitivity.pdf", fig)
