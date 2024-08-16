using CSV
using DataFrames
using CairoMakie
using LaTeXStrings

function lrange(from::Float64, to::Float64, length::Int64 = 6)
    exp.(collect(range(log(from), log(to); length = length)))
end


df = CSV.read("output/branch_specific_estimation_error.csv", DataFrame,
                missingstring=["NA", "NAN", "NULL"])


tree_heights = [25, 50, 75, 100, 125]
rate_variation = ["tiny", "small", "moderate", "large"]
inferences = ["true_all", "true_div", "unknown_rates"]

n_models = length(rate_variation)
n_heights = length(tree_heights)
n_inferences = length(inferences)

titles = [
    L"\text{tiny}",
    L"\text{small}",
    L"\text{moderate}",
    L"\text{large}"
]

yt = lrange(0.25, 4.0, 5)
fig = Figure(resolution = (650, 370), fontsize = 14, 
            figure_padding = (1,1,1,1))
axs = []
for (i, inference) in enumerate(inferences)
    for j in 1:n_models
        if i == 1
            title = titles[j]
        else
            title = ""
        end
        ax = Axis(fig[i,j+1], 
                    xscale = Makie.log10,
                    yscale = Makie.log10,
                    xgridvisible = false,
                    ygridvisible = false,
                    topspinevisible = false,
                    title = title,
                    yticks = yt,
                    rightspinevisible = false)
        CairoMakie.ylims!(ax, 0.2, 7.0)
        if j > 1
            hideydecorations!(ax, ticks = false)
        end
        if i < 3
            hidexdecorations!(ax, ticks = false)
        end

        for (q, h) in enumerate(tree_heights)
            #println("$i $j $h")
            df4 = subset(
                df,
                :criterion => x -> x .== "N_over_half",
                :inference => x -> x .== inference,
                :height =>    x -> x .== h,
                :model  =>    x -> x .== j,
                :N_true_sum => x -> x .> 0,
            )
            x = collect(df4[!,:ntaxa])
            y = collect(df4[!,:prop_error_geomean_lambda])
            scatter!(ax, x, y, markersize = 7,
                    strokewidth = 0.5,
                    label = LaTeXString(string(Int64(tree_heights[q])) * raw" Ma"))
        end
        CairoMakie.lines!(ax, [extrema(df[!,:ntaxa])...], 
                        [1.0, 1.0], label = L"\text{no error}",
                         color = :red, linewidth = 2, linestyle = :dash)

        append!(axs, [ax])
    end
end
linkaxes!(axs...)

## fake LABELS
ylabel = Label(fig[1:3, 1], L"\text{proportional error (speciation rate)}", rotation = π/2)
xlabel = Label(fig[4, 2:5], L"\text{number of taxa}")


xlabel = Label(fig[0, 2:5], L"\text{allowed rate variation}")
ylabel = Label(fig[1, 6], L"\text{true}~𝛌,𝛍,\eta", rotation = π/2)
ylabel = Label(fig[2, 6], L"\text{true}~𝛌,𝛍", rotation = π/2)
ylabel = Label(fig[3, 6], L"\text{unknown}~𝛌,𝛍,\eta", rotation = π/2)
fig

labels = [
    L"\text{25 Ma}",
    L"\text{50 Ma}",
    L"\text{75 Ma}",
    L"\text{100 Ma}",
    L"\text{125 Ma}"
]
title = "Shift in"
fig[1:3, 7] = Legend(fig, axs[1], "", framevisible = false, patchsize = (30, 30))
colsize!(fig.layout, 7, Relative(0.2))
colgap!(fig.layout, 7)
rowgap!(fig.layout, 7)

rowsize!(fig.layout, 0, Relative(0.05))
#rowsize!(fig.layout, 1, Relative(0.05))
for i in 1:3
    rowsize!(fig.layout, i, Relative(0.85/3))
end
rowsize!(fig.layout, 4, Relative(0.05))
fig

CairoMakie.save("figures/proportional-error-lambda.pdf", fig)

##################################
##
##   Do the equivalent but for extinction rates
##
####################################

#yt = lrange(0.125, 8.0, 5)
#yt = round.(lrange(0.1, 10.0, 5); digits = 2)
yt = round.(lrange(0.111111111111, 9.0, 5); digits = 2)
fig = Figure(size = (650, 370), fontsize = 14, 
            figure_padding = (1,1,1,1))
axs = []
for (i, inference) in enumerate(inferences)
    for j in 1:n_models
        if i == 1
            title = titles[j]
        else
            title = ""
        end
        ax = Axis(fig[i,j+1], 
                    xscale = Makie.log10,
                    yscale = Makie.log10,
                    xgridvisible = false,
                    ygridvisible = false,
                    topspinevisible = false,
                    title = title,
                    yticks = yt,
                    rightspinevisible = false)
        CairoMakie.ylims!(ax, 0.05, 16.0)
        if j > 1
            hideydecorations!(ax, ticks = false)
        end
        if i < 3
            hidexdecorations!(ax, ticks = false)
        end

        for (q, h) in enumerate(tree_heights)
            #println("$i $j $h")
            df4 = subset(
                df,
                :criterion => x -> x .== "N_over_half",
                :inference => x -> x .== inference,
                :height =>    x -> x .== h,
                :model  =>    x -> x .== j
            )
            x = collect(df4[!,:ntaxa])
            y = collect(df4[!,:prop_error_geomean_mu])
            scatter!(ax, x, y, markersize = 7,
                    strokewidth = 0.5,
                    label = LaTeXString(string(Int64(tree_heights[q])) * raw" Ma"))
        end
        CairoMakie.lines!(ax, [extrema(df[!,:ntaxa])...], 
                        [1.0, 1.0], label = L"\text{no error}",
                         color = :red, linewidth = 2, linestyle = :dash)

        append!(axs, [ax])
    end
end
linkaxes!(axs...)

## fake LABELS
ylabel = Label(fig[1:3, 1], L"\text{proportional error (extinction rate)}", rotation = π/2)
xlabel = Label(fig[4, 2:5], L"\text{number of taxa}")


xlabel = Label(fig[0, 2:5], L"\text{allowed rate variation}")
ylabel = Label(fig[1, 6], L"\text{true}~𝛌,𝛍,\eta", rotation = π/2)
ylabel = Label(fig[2, 6], L"\text{true}~𝛌,𝛍", rotation = π/2)
ylabel = Label(fig[3, 6], L"\text{unknown}~𝛌,𝛍,\eta", rotation = π/2)

##
fig[1:3, 7] = Legend(fig, axs[1], "", framevisible = false, patchsize = (30, 30))
colsize!(fig.layout, 7, Relative(0.2))
colgap!(fig.layout, 7)
rowgap!(fig.layout, 7)

rowsize!(fig.layout, 0, Relative(0.05))
#rowsize!(fig.layout, 1, Relative(0.05))
for i in 1:3
    rowsize!(fig.layout, i, Relative(0.85/3))
end
rowsize!(fig.layout, 4, Relative(0.05))
fig

CairoMakie.save("figures/proportional-error-mu.pdf", fig)


###########################################
##
##    Figure 11 -- false error ratios/accuracy
## 
###########################################


fig = Figure(size = (450, 150), fontsize = 14, 
            figure_padding = (1,1,1,1))

axs = []
xlabels = [
    "accuracy",
    "false positive ratio",
    #"false omission ratio",
    "false negative ratio"
]
yt = [0, 400, 800, 1200, 1600]
for i in 1:3
    ax = Axis(fig[1,i], 
                #yscale = Makie.log10,
                xgridvisible = false,
                ygridvisible = false,
                topspinevisible = false,
                xlabel = xlabels[i],
                ylabel = "frequency (trees)",
                #yticks = yt,
                rightspinevisible = false)
    xlims!(ax, -0.1, 1.1)
    append!(axs, [ax])
end
ax1, ax2, ax3 = axs

df_empirical_bayes = subset(
    df,
    :inference => x -> x .== "unknown_rates",
    #:criterion => x -> x .== "N_over_half",
    #:criterion => x -> x .== "bayes_factor_10",
    :criterion => x -> x .== "N_half_and_bayes_factor",
)
fp = df_empirical_bayes[!,"false positive"]
fn = df_empirical_bayes[!,"false negative"]
tp = df_empirical_bayes[!,"true positive"]
tn = df_empirical_bayes[!,"true negative"]
fpr = fp ./ (fp .+ tn)
fnr = fn ./ (tp .+ fn)

acc = (tp .+ tn) ./ (tp .+ tn .+ fp .+ fn)
false_ommission_ratio = fn ./ (fn .+ tn)

import StatsPlots
p1 = StatsPlots.violin(fpr1[fpr1 .< 0.01], title = "bayes factor > 10")
p2 = StatsPlots.violin(fpr2[fpr2 .< 0.01], title = "N > 10")
StatsPlots.plot(p1, p2, ylim = (0.0, 0.01), xlabel = "density", ylabel = "false positive ratio", grid = false)


bins = 15
CairoMakie.hist!(ax1, acc, color = :black)#, bins = bins)
CairoMakie.hist!(ax2, fpr, color = :black, bins = 1)
#CairoMakie.hist!(ax3, false_ommission_ratio, color = :black, bins = bins)
CairoMakie.hist!(ax3, fnr[.!isnan.(fnr)], color = :black)#, bins = bins)

for i in 2:3
    hideydecorations!(axs[i], ticks = false)
end
colgap!(fig.layout, 7)
linkaxes!(axs...)
fig
CairoMakie.save("figures/number-of-shifts-confusion_bayesfactor_and_N_half.pdf", fig)


function mean(x::Vector{Float64})
    s = sum(x)
    res = s / length(x)
    return(res)
end

metrics = [
    "mean false positive ratio            " => mean(fpr),
    "mean false negative ratio            " => mean(fnr[.!isnan.(fnr)]),
    "mean false ommission ratio           " => mean(false_ommission_ratio),
    "mean accuracy                        " => mean(acc),
    "frequency at least one false positive" => sum(fpr .> 0) / length(fpr)
]

using Printf

println("Shift criterion: Bayes factor > 10 and N > 0.5")
#println("Shift criterion: N > 0.5")
for (name, metric) in metrics
    println(name, "\t", @sprintf "%f" metric)
end

fnr

fig = Figure()
ax = Axis(fig[1,1],
        xscale = Makie.log10,
        xlabel = "number of taxa",
        ylabel = "false positive ratio")
        #yscale = Makie.log10)
scatter!(ax, df_empirical_bayes[!,:ntaxa], fpr)
fig


## plot metrics by tree size


fig5 = Figure(size = (600, 400));

yt = [0.0, 0.25, 0.50, 0.75, 1.0]

xs = [1,2,3,4]
xtl = ["(0,50]", "(50,250]", "(250,1000]", "(1000,∞]"]

ax1 = Axis(fig5[1,1], 
        ylabel = L"\text{mean accuracy}",
        xgridvisible = false,
        ygridvisible = false,
        topspinevisible = false,
        yticks = yt,
        xticks = (xs, xtl),
        xticklabelrotation = pi/2,
        rightspinevisible = false)
ax2 = Axis(fig5[1,2], 
        ylabel = L"\text{mean false positive ratio}",
        xgridvisible = false,
        ygridvisible = false,
        topspinevisible = false,
        xticks = (xs, xtl),
        #yticks = yt,
        xticklabelrotation = pi/2,
        rightspinevisible = false)
ax3 = Axis(fig5[1,3], 
        ylabel = L"\text{mean false negative ratio}",
        xgridvisible = false,
        ygridvisible = false,
        topspinevisible = false,
        yticks = yt,
        xticks = (xs, xtl),
        xticklabelrotation = pi/2,
        rightspinevisible = false)
        
ylims!(ax1, (0.0, 1.0))
#ylims!(ax2, (0.0, 0.001))
ylims!(ax3, (0.0, 1.0))


function metrics(df, ntaxa_lower, ntaxa_upper)
    df = deepcopy(df)
    df = subset(df,
        :ntaxa =>  x -> x .> ntaxa_lower,
    )
    df = subset(df,
        :ntaxa =>  x -> x .<= ntaxa_upper,
    )

    fp = df[!,"false positive"]
    fn = df[!,"false negative"]
    tp = df[!,"true positive"]
    tn = df[!,"true negative"]
    fpr = fp ./ (fp .+ tn)
    fnr = fn ./ (tp .+ fn)

    fnr = fnr[.!isnan.(fnr)]

    acc = (tp .+ tn) ./ (tp .+ tn .+ fp .+ fn)

    mean_acc = mean(acc)
    mean_fpr = mean(fpr)
    mean_fnr = mean(fnr)

    return(mean_acc, mean_fpr, mean_fnr)
end

ymetrics = zeros(4,3)

ymetrics[1,:] .= metrics(df_empirical_bayes, 0.0, 50.0)
ymetrics[2,:] .= metrics(df_empirical_bayes, 50.0, 250.0)
ymetrics[3,:] .= metrics(df_empirical_bayes, 250.0, 1000.0)
ymetrics[4,:] .= metrics(df_empirical_bayes, 1000.0, 1e20)

barplot!(ax1, xs, ymetrics[:,1], color = "gray")
barplot!(ax2, xs, ymetrics[:,2], color = "gray")
barplot!(ax3, xs, ymetrics[:,3], color = "gray")

ylabel = Label(fig5[2, 1:3], L"\text{tree size (tips)}")
fig5
CairoMakie.save("figures/metrics-with-tree-size.pdf", fig5)
