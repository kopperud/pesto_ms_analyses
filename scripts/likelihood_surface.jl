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
function foo_make_model(λmean, μmean, η)
    d1 = LogNormal(log(λmean), H)
    d2 = LogNormal(log(μmean), H)

    n = 10
    speciation = make_quantiles(d1, n)
    extinction = make_quantiles(d2, n)

    λ, μ = allpairwise(speciation, extinction)

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
n = 150
m = 160

logls1_unscaled = zeros(n,m)
rate_for_one_shift = 1 / sum(primates.branch_lengths)

λs1 = range(0.1,0.5; length = n)
μs1 = range(0.01,0.4; length = m)

prog = ProgressMeter.Progress(n*m, "calculating likelihoods...");
for (i, λmean) in enumerate(λs1)
    for (j, μmean) in enumerate(μs1)
        model = foo_make_model(λmean, μmean, rate_for_one_shift)
        logls1_unscaled[i,j] = logL_root(model, primates)
        ProgressMeter.next!(prog)
    end
end
#############################
#
#  keep mu constant
#
############################
logls2_unscaled = zeros(n,m)

λs2 = range(0.1,0.45; length = n)
#ηs2 = range(0.00001,0.1; length = m)
ηs2 = lrange(0.00001,0.1, m)

prog = ProgressMeter.Progress(n*m, "calculating likelihoods...");
for (i, λmean) in enumerate(λs2)
    for (j, η) in enumerate(ηs2)
        model = foo_make_model(λmean, 0.22, η)
        logls2_unscaled[i,j] = logL_root(model, primates)
        ProgressMeter.next!(prog)
    end
end

ml2 = maximum(logls2_unscaled)
logls2 = logls2_unscaled .- ml2
logls2[logls2 .< -6] .= -6
fig = Figure(); ax = Axis3(fig[1,1], azimuth = 1.6*π); CairoMakie.surface!(ax, λs2, ηs2, logls2); fig
#############################
#
#  keep lambda constant
#
############################
logls3_unscaled = zeros(n,m)

μs3 = range(0.01,0.3; length = n)
ηs3 = lrange(0.00001,0.05, m)

prog = ProgressMeter.Progress(n*m, "calculating likelihoods...");
for (i, μmean) in enumerate(μs3)
    for (j, η) in enumerate(ηs3)
        model = foo_make_model(0.30, μmean, η)
        logls3_unscaled[i,j] = logL_root(model, primates)
        ProgressMeter.next!(prog)
    end
end


## make it top at 0
ml1 = maximum(logls1_unscaled)
ml2 = maximum(logls2_unscaled)
ml3 = maximum(logls3_unscaled)
logls1 = logls1_unscaled .- ml1
logls2 = logls2_unscaled .- ml2
logls3 = logls3_unscaled .- ml3

## make it bottom out at Δ
Δ = -2
logls1[logls1 .< Δ] .= Δ
Δ = -6
logls2[logls2 .< Δ] .= Δ
logls3[logls3 .< Δ] .= Δ


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
           xlabel = L"\text{speciation rate }(\hat{\lambda})",
           ylabel = L"\text{extinction rate }(\hat{\mu})",
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
    λs1,
    μs1, logls1,
    colormap = :viridis)
mle_index = argmax(logls1)
λ_mle = λs1[mle_index[1]]; μ_mle = μs1[mle_index[2]];
z_mle = logls1[mle_index] 
z_twostep = logL_root(SSEconstant(λ, μ, rate_for_one_shift), primates) - ml1
kwargs = (; markersize = 10, strokecolor = :black, strokewidth = 1)
sc1 = CairoMakie.scatter!(ax1, λ_mle, μ_mle, z_mle; kwargs..., color = :gray, label = L"\text{joint maximum likelihood}")
sc2 = CairoMakie.scatter!(ax1, λ_twostep, μ_twostep, z_twostep; kwargs..., color = :white, label = L"\text{joint maximum likelihood}")

## create the legends
empty_ax = Axis(fig[1,2],
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

## keep mu fixed
## panel b)
ax2 = Axis3(fig[2,1],
           xlabel = L"\text{speciation rate }(\hat{\lambda})",
           ylabel = L"\text{shift rate }(\eta)",
           zlabel = L"\text{Δlogl}",
           azimuth = 1.35*π,
          xgridvisible = false, 
          ygridvisible = false,
          zgridvisible = false,
          viewmode = :stretch,
          title = L"\text{b) extinction rate fixed}",
          protrusions = protrusions,
          )

CairoMakie.surface!(
    ax2,
    λs2,
    ηs2, 
    logls2,
    colormap = :viridis)
mle_index = argmax(logls2)
λ_mle = λs2[mle_index[1]]; η_mle = ηs2[mle_index[2]];
z_mle = logls2[mle_index] 
z_twostep = logL_root(SSEconstant(λ, μ, η_twostep), primates) - ml2
kwargs = (; markersize = 10, strokecolor = :black, strokewidth = 1)
sc1 = CairoMakie.scatter!(ax2, λ_mle, η_mle, z_mle; kwargs..., color = :gray, label = L"\text{joint maximum likelihood}")
sc2 = CairoMakie.scatter!(ax2, λ_twostep, η_twostep, z_twostep; kwargs..., color = :white, label = L"\text{joint maximum likelihood}")


## keep lambda fixed
## panel c)
ax3 = Axis3(fig[2,2],
           xlabel = L"\text{extinction rate }(\hat{\mu})",
           ylabel = L"\text{shift rate }(\eta)",
           zlabel = L"\text{Δlogl}",
           azimuth = 1.35*π,
          xgridvisible = false, 
          ygridvisible = false,
          zgridvisible = false,
          viewmode = :stretch,
          title = L"\text{c) speciation rate fixed}",
          protrusions = protrusions,
          #yreversed = true
          )

CairoMakie.surface!(
    ax3,
    μs3,
    ηs3, 
    logls3,
    colormap = :viridis)
mle_index = argmax(logls3)
μ_mle = μs3[mle_index[1]]; η_mle = ηs3[mle_index[2]];
z_mle = logls3[mle_index] 
z_twostep = logL_root(SSEconstant(λ, μ, η_twostep), primates) - ml3
kwargs = (; markersize = 10, strokecolor = :black, strokewidth = 1)
sc1 = CairoMakie.scatter!(ax3, μ_mle, η_mle, z_mle; kwargs..., color = :gray, label = L"\text{joint maximum likelihood}")
sc2 = CairoMakie.scatter!(ax3, μ_twostep, η_twostep, z_twostep; kwargs..., color = :white, label = L"\text{joint maximum likelihood}")



for i in 1:2
    rowsize!(fig.layout, i, Relative(0.5))
    colsize!(fig.layout, i, Relative(0.5))
end

rowgap!(fig.layout, 0.0)
colgap!(fig.layout, 0.0)

fig

CairoMakie.save("figures/likelihood-surface.pdf", fig) ## for some reason the pdf is 60 megabytes
CairoMakie.save("figures/likelihood-surface.png", fig)
