using Distributions
using Pesto
using LaTeXStrings
using CairoMakie
using CSV
using ProgressMeter

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

λ_twostep, μ_twostep = Pesto.estimate_constant_bdp(primates)
d1 = LogNormal(log(λ_twostep), H)
d2 = LogNormal(log(μ_twostep), H)

n = 10
speciation = make_quantiles(d1, n)
extinction = make_quantiles(d2, n)

k = n ^2
λ, μ = allpairwise(speciation, extinction)

η_twostep = optimize_eta(λ, μ, primates)
println(λ_twostep, "\t", μ_twostep, "\t", η_twostep)

model = SSEconstant(λ, μ, η_twostep)

#############################
#
#  calculate likelihood grids
#
#############################
function foo_make_model2(netdiv, relext, η)
    #d1 = LogNormal(log(netdiv), H)
    d1 = Normal(netdiv, H)
    #d2 = LogNormal(log(relext), H)
    β = 4
    α = relext * β / (1 - relext)
    d2 = Beta(α, β)

    n = 20
    q1 = make_quantiles(d1, n)
    q2 = make_quantiles(d2, n)

    r, ϵ = allpairwise(q1, q2)

    λ = r ./ (1 .- ϵ)
    μ = λ .- r

    μ = μ[λ .> 0]
    λ = λ[λ .> 0]

    model = SSEconstant(λ, μ, η)
    return(model)
end


function lrange(from::Float64, to::Float64, length::Int64 = 6)
    exp.(collect(range(log(from), log(to); length = length)))
end






#############################
#
#  keep eta constant
#
#############################
n = 65
m = 70

logls1_unscaled = zeros(n,m)
rate_for_one_shift = 1 / sum(primates.branch_lengths)

rs = range(-0.1, 0.3; length = n)
ϵs = range(0.01,0.9; length = m)

models = Array{SSEconstant}(undef, n, m)
prog = ProgressMeter.Progress(n*m, "calculating likelihoods...");
for (i, rmean) in enumerate(rs)
    for (j, ϵmean) in enumerate(ϵs)
        model = foo_make_model2(rmean, ϵmean, rate_for_one_shift*10)
        logls1_unscaled[i,j] = logL_root(model, primates)
        models[i,j] = model
        ProgressMeter.next!(prog)
    end
end

## make it top at 0
ml1 = maximum(logls1_unscaled)
logls1 = logls1_unscaled .- ml1

## make it bottom out at Δ
Δ = -6
logls1[logls1 .< Δ] .= Δ

###################
##
##   assemble figure
##
#################

fig = Figure(figure_padding = 0, size = (650,500))

#(left, right, bottom, top)
protrusions = (65, 30, 60, 25)

## panel a)
ax1 = Axis3(fig[1,1],
           xlabel = L"\text{net div rate}",
           ylabel = L"\text{rel ext rate}",
           zlabel = L"\text{Δlogl}",
           azimuth = 1.35*π,
          xgridvisible = false, 
          ygridvisible = false,
          zgridvisible = false,
          viewmode = :stretch,
          title = L"\text{a) shift rate fixed}",
          protrusions = protrusions,
          )

CairoMakie.surface!(
    ax1,
    rs,
    ϵs, logls1,
    colormap = :viridis)

mle_index = argmax(logls1)
r_mle = rs[mle_index[1]]; ϵ_mle = ϵs[mle_index[2]];
z_mle = logls1[mle_index] 
#z_twostep = logL_root(SSEconstant(λ, μ, rate_for_one_shift), primates) - ml1
kwargs = (; markersize = 10, strokecolor = :black, strokewidth = 1)
sc1 = CairoMakie.scatter!(ax1, r_mle, ϵ_mle, z_mle; kwargs..., color = :gray, label = L"\text{joint maximum likelihood}")
#sc2 = CairoMakie.scatter!(ax1, λ_twostep, μ_twostep, z_twostep; kwargs..., color = :white, label = L"\text{joint maximum likelihood}")

## create the legends
empty_ax = Axis(fig[2,1],
    topspinevisible = false,
    leftspinevisible = false,
    bottomspinevisible = false,
    rightspinevisible = false,
) 
hidedecorations!(empty_ax)
labels = [
    L"\text{joint maximum likelihood}", 
    L"\text{two step approach}", 
]
axislegend(empty_ax, [sc1, sc2], labels, position = :cc)
#=
Colorbar(fig[1, 2], limits = (-2, 0),
    vertical = false,
    colormap = :viridis, size = 25,
    label = L"\text{\Delta logl}",
)
=#

rowsize!(fig.layout, 1, Relative(0.8))
rowgap!(fig.layout, 0.0)
colgap!(fig.layout, 0.0)

fig

#CairoMakie.save("figures/likelihood-surface.pdf", fig) ## for some reason the pdf is 60 megabytes
#CairoMakie.save("figures/likelihood-surface.png", fig)
