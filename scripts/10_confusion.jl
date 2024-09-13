###########################################
##
##    Figure 10 -- false error ratios/accuracy
## 
###########################################

using CSV
using DataFrames
using CairoMakie
using LaTeXStrings

function lrange(from::Float64, to::Float64, length::Int64 = 6)
    exp.(collect(range(log(from), log(to); length = length)))
end


df = CSV.read("output/branch_specific_estimation_error.csv", DataFrame,
                missingstring=["NA", "NAN", "NULL"])

################################
##
## prepare metrics
##
################################
fp = df[!,"false positive"]
fn = df[!,"false negative"]
tp = df[!,"true positive"]
tn = df[!,"true negative"]
fpr = fp ./ (fp .+ tn)
fnr = fn ./ (tp .+ fn)
#fnr = fnr[.!isnan.(fnr)]

acc = (tp .+ tn) ./ (tp .+ tn .+ fp .+ fn)

df[!,"accuracy"] = acc
df[!,"false_negative_ratio"] = fnr
df[!,"false_positive_ratio"] = fpr

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

##################
##
## prepare the figure
##
##################
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


function mymean(x::Vector{Float64})
    s = sum(x)
    res = s / length(x)
    return(res)
end

function standard_error(x::Vector{Float64})
    m = mymean(x)

    v = sum((x .- m) .^2) / (length(x)-1)
    var_mean = v / length(x)
    se_mean = sqrt(var_mean)

    return(se_mean)
end

mymetrics = [
    "mean false positive ratio            " => mymean(fpr),
    "mean false negative ratio            " => mymean(fnr[.!isnan.(fnr)]),
    "mean false ommission ratio           " => mymean(false_ommission_ratio),
    "mean accuracy                        " => mymean(acc),
    "frequency at least one false positive" => sum(fpr .> 0) / length(fpr)
]

using Printf

println("Shift criterion: Bayes factor > 10 and N > 0.5")
#println("Shift criterion: N > 0.5")
for (name, metric) in mymetrics
    println(name, "\t", @sprintf "%f" metric)
end


## make bins for trees as a function of number of taxa
tree_size_bins = Int64[]
for nt in df[!,:ntaxa]
    if (nt <= 50)
        push!(tree_size_bins, 1)
    elseif (nt > 50) & (nt <= 250)
        push!(tree_size_bins, 2)
    elseif (nt > 250) & (nt <= 1000)
        push!(tree_size_bins, 3)
    elseif  nt > 1000
        push!(tree_size_bins, 4)
    end
end
df[!,:tree_size_bin] = tree_size_bins


dodge = Int64[]
for N in df[!,:N_true_sum]
    if N > 0
        push!(dodge, 1)
    else
        push!(dodge, 2)
    end
end
df[!,:dodge] = dodge



fig = Figure(size = (600, 400));

xt = [1,2,3,4]
titles = [L"\text{(0,50]}", L"\text{(50,250]}", L"\text{(250,1000]}",L"\text{(1000,∞]}"]
xtl = [
    L"\text{tiny}",
    L"\text{small}",
    L"\text{moderate}",
    L"\text{large}"
]

shift_or_no_shift = [
    x -> x .== 0,
    x -> x .> 0,
]

summaries1 = []

for metric in [:accuracy, :false_positive_ratio, :false_negative_ratio]
    for (i, crit) in enumerate(shift_or_no_shift)
        for taxa_category in 1:4
            df1 = subset(
                    df,
                    :inference => x -> x .== "unknown_rates",
                   # :criterion => x -> x .== "bayes_factor_10",
                    :criterion => x -> x .== "N_half_and_bayes_factor",
                    :tree_size_bin => x -> x .== taxa_category,
                    :N_true_sum => crit,
                    )
            y = df1[!,metric]
            y = y[.!isnan.(y)]
            value = mymean(y)
            se = standard_error(y)


            summary_df = DataFrame(
                                :value => value,
                                :se => se,
                                :metric => metric,
                                :n => size(df1)[1],
                                :tree_size_bin => taxa_category,
                                :shifts => i-1,
                                  )
            push!(summaries1, summary_df)
        end
    end
end

summaries2 = []

for metric in [:accuracy, :false_positive_ratio, :false_negative_ratio]
    for (i, crit) in enumerate(shift_or_no_shift)
        for model_index in 1:4
            df1 = subset(
                    df,
                    :inference => x -> x .== "unknown_rates",
                    :criterion => x -> x .== "N_half_and_bayes_factor",
                    :model => x -> x .== model_index,
                    :N_true_sum => crit,
                    )
            y = df1[!,metric]
            y = y[.!isnan.(y)]
            value = mymean(y)
            se = standard_error(y)


            summary_df = DataFrame(
                                :value => value,
                                :se => se,
                                :metric => metric,
                                :n => size(df1)[1],
                                :model => model_index,
                                :shifts => i-1,
                                  )
            push!(summaries2, summary_df)
        end
    end
end


df5a = reduce(vcat, summaries1)
df5b = reduce(vcat, summaries2)


fig4 = Figure(size = (700, 400))
ms = [:accuracy, :false_positive_ratio, :false_negative_ratio]


axs = []
for (row, metric) in enumerate(ms) 
    df6 = subset(df5a,
             :metric => x -> x .== metric,
            )
    xt = 1:4
    xtl = [L"(0,50]", L"(50,250]",L"(250,1000]",L"(1000,∞)"]
    ax = Axis(
            fig4[row,1],
            topspinevisible = false,
            rightspinevisible = false,
            xgridvisible = false,
            ygridvisible = false,
            xticks = (xt, xtl),
            xlabel = L"\text{number of taxa}",
        )
    if row < 3
        hidexdecorations!(ax, ticks = false)
    end
    push!(axs, ax)

    df_shift = subset(df6, :shifts => x -> x .> 0)
    df_no_shift = subset(df6, :shifts => x -> x .== 0)

    barplot!(ax, 
             df_shift[!,:tree_size_bin] .- 0.25, 
             df_shift[!,:value],
             color = :white, 
             strokecolor = :black,
             strokewidth = 1,
             width = 0.5,
            )
    errorbars!(ax,
               df_shift[!,:tree_size_bin] .- 0.25,
               df_shift[!,:value],
               df_shift[!,:se],
               color = :black
              )
    barplot!(ax, 
             df_no_shift[!,:tree_size_bin] .+ 0.25, 
             df_no_shift[!,:value],
             color = Pattern("/", background_color = :white, linecolor = :black),
             strokecolor = :black,
             strokewidth = 1,
             width = 0.5,
            )
    errorbars!(ax,
               df_no_shift[!,:tree_size_bin] .+ 0.25,
               df_no_shift[!,:value],
               df_no_shift[!,:se],
               color = :black
              )
end



for (row, metric) in enumerate(ms) 
    df6 = subset(df5b,
             :metric => x -> x .== metric,
            )
    xt = 1:4
    xtl = [L"\text{tiny}", L"\text{small}",L"\text{moderate}",L"\text{large}"]
    ax = Axis(
            fig4[row,2],
            topspinevisible = false,
            rightspinevisible = false,
            xgridvisible = false,
            ygridvisible = false,
            xticks = (xt, xtl),
            xlabel = L"\text{allowed rate variation}",
        )
    if row < 3
        hidexdecorations!(ax, ticks = false)
    end
    push!(axs, ax)

    df_shift = subset(df6, :shifts => x -> x .> 0)
    df_no_shift = subset(df6, :shifts => x -> x .== 0)

    barplot!(ax, 
             df_shift[!,:model] .- 0.25, 
             df_shift[!,:value],
             color = :white, 
             strokecolor = :black,
             strokewidth = 1,
             width = 0.5,
            )
    errorbars!(ax,
               df_shift[!,:model] .- 0.25,
               df_shift[!,:value],
               df_shift[!,:se],
               color = :black
              )
    barplot!(ax, 
             df_no_shift[!,:model] .+ 0.25, 
             df_no_shift[!,:value],
             color = Pattern("/", background_color = :white, linecolor = :black),
             strokecolor = :black,
             strokewidth = 1,
             width = 0.5,
            )
    errorbars!(ax,
               df_no_shift[!,:model] .+ 0.25,
               df_no_shift[!,:value],
               df_no_shift[!,:se],
               color = :black
              )
              
end

if false
    ylims!(axs[1], (0.96, 1.01))
    ylims!(axs[4], (0.96, 1.01))

    ylims!(axs[2], (0.00, 0.00015))
    ylims!(axs[5], (0.00, 0.00015))

    ylims!(axs[3], (0.85, 1.01))
    ylims!(axs[6], (0.85, 1.01))
else
    for i in 1:6
        ylims!(axs[i], (-0.03, 1.03))
    end
end
fig4

linkxaxes!(axs[1:3]...)
linkxaxes!(axs[4:6]...)

linkyaxes!(axs[1], axs[4])
linkyaxes!(axs[2], axs[5])
linkyaxes!(axs[3], axs[6])


Label(fig4[0,1], L"\text{a) split by tree size}")
Label(fig4[0,2], L"\text{b) split by model}")
Label(fig4[1,0], L"\text{accuracy}", rotation = π/2)
Label(fig4[2,0], L"\text{FPR}", rotation = π/2)
Label(fig4[3,0], L"\text{FNR}", rotation = π/2)
groups = [
      PolyElement(color = :white, strokecolor = :black, strokewidth = 1),
    PolyElement(color = Pattern("/", background_color = :white, linecolor = :black), strokecolor = :black, strokewidth = 1),
   ]
Legend(fig4[4,1:2], 
       groups, 
       [L"\geq 1 \text{ rate shifts}", L"\text{no rate shifts}"],
      )

for row in 1:3
    rowsize!(fig4.layout, row, Relative(0.28))
end

for col in 1:2
    colsize!(fig4.layout, col, Relative(0.45))
end

fig4

#CairoMakie.save("figures/confusion-split-simple.pdf", fig4)
#CairoMakie.save("figures/confusion-split-simple_ylimit_01.pdf", fig4)
         
