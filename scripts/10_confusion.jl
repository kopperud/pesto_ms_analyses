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
ylabel = Label(fig5[2, 1:3], L"\text{tree size (tips)}")
fig5
CairoMakie.save("figures/metrics-with-tree-size.pdf", fig5)

## make bin size
tree_size_bins = Int64[]
for nt in df[!,:ntaxa]
    if (nt < 50)
        push!(tree_size_bins, 1)
    elseif (nt >= 50) & (nt < 250)
        push!(tree_size_bins, 2)
    elseif (nt >= 250) & (nt < 1000)
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


axs = []
for (j, metric) in enumerate([:accuracy, :false_positive_ratio, :false_negative_ratio])
    for i in 1:4
        if i == 1
            ylabel = L"\text{accuracy}"
        else
            ylabel = ""
        end

        if j == 3
            xtl1 = xtl
        else
            xtl1 = ["" for _ in 1:4]
        end

        if j == 1 
            title = titles[j]
        else
            title = ""
        end

        ax = Axis(
            fig[j+1,i+1],
            xgridvisible = false,
            ygridvisible = false,
            topspinevisible = false,
            #ylabel = ylabel,
            xticklabelrotation = pi/2,
            #xlabel = L"\text{rate variation}",
            title = title,
            xticks = (xt, xtl1),
            rightspinevisible = false,
        )
        
        df1 = subset(
            df,
            :tree_size_bin => x -> x .== i,
            :inference => x -> x .== "unknown_rates",
        )
        
        values = df1[!,metric]
        categories = df1[!,:model]
        d = df1[!,:dodge]

        ## remove NaN, happens for fnr
        categories = categories[.!isnan.(values)] 
        d = d[.!isnan.(values)] 
        values = values[.!isnan.(values)] 
        
        CairoMakie.boxplot!(ax, categories, values, 
            dodge = d,
            color = map(x->x==1 ? :gray : :orange, d),
        )

        #ylims!(ax, [0.9, 1.01])
        push!(axs, ax)
    end
end

## extra labels
overall_title = Label(fig[1, 2:5], L"\text{tree size (number of tips)}")
xlabel = Label(fig[5, 2:5], L"\text{allowed rate variation}")
ylabel = Label(fig[2, 1], L"\text{accuracy}", rotation = π/2)
ylabel = Label(fig[3, 1], L"\text{FPR}", rotation = π/2)
ylabel = Label(fig[4, 1], L"\text{FNR}", rotation = π/2)


for ax in axs[1:4]
    ylims!(ax, [0.9, 1.01])
end

for ax in axs[5:8]
    ylims!(ax, [-0.01,0.01])
end

for ax in axs[1:4]
    ylims!(ax, [0.49, 1.01])
end

linkaxes!(axs[1:4]...)
linkaxes!(axs[5:8]...)
linkaxes!(axs[9:12]...)


for i in [2,3,4]
    rowsize!(fig.layout, i, Relative(0.3))
end

fig
CairoMakie.save("figures/number-of-shifts-confusion-split.pdf", fig)


fig2 = Figure(size = (500, 800));
xtl = [2, 20, 200, 2000, 20_000]
xt = log10.(xtl)
xtl = string.(xtl)

cl_maps = [:blues, :reds, :jet, :greens]
axs = []
for i in 1:4
        
    df1 = subset(
        df,
        #:tree_size_bin => x -> x .== i,
        :model => x -> x .== i,
        :inference => x -> x .== "unknown_rates",
        :criterion => x -> x .== "N_half_and_bayes_factor",
        #:N_true_sum => x -> x .== 0,
    )


    ax = Axis(fig2[i,1], 
                    #yscale = Makie.log10,
                    #xscale = Makie.log10,
                    xgridvisible = false,
                    ygridvisible = false,
                    topspinevisible = false,
                    xlabel = L"\text{number of taxa}",
                    ylabel = L"\text{accuracy}",
                    xticks = (xt, xtl),
                    #yticks = yt,
                    rightspinevisible = false)
    xlims!(ax, (log10(1), log10(100_000)))
    #ylims!(ax, (0.0, 1.0))
    push!(axs, ax)
    hb = CairoMakie.hexbin!(
        ax, 
        log10.(df1[!,:ntaxa]), 
        df1[!,:accuracy],
        colorscale = log10,
        colormap = :greens,
        )
    Colorbar(fig2[i, 2], hb,
        label = L"\text{number of trees}",
        height = Relative(0.5)
    )    
end
linkaxes!(axs...)
fig2

log10.(df1[!,:ntaxa])

xt = [1, 10, 100, 1000, 10_000]
fig3 = Figure(size = (800, 600)); 
xtl = string.(yt)

textpos = [0.25, -0.1, -0.1]
for (row, metric) in zip(1:3, [:accuracy, :false_positive_ratio, :false_negative_ratio])
    axs = []
    ax_counter = 1
    for category in 1:2
        for i in 1:4
            if category == 1
                df1 = subset(
                    df,
                    :tree_size_bin => x -> x .== i,
                    :inference => x -> x .== "unknown_rates",
                    #:criterion => x -> x .== "N_half_and_bayes_factor",
                    :criterion => x -> x .== "bayes_factor_10",
                )
            else
                df1 = subset(
                    df,
                    :model => x -> x .== i,
                    :inference => x -> x .== "unknown_rates",
                    :criterion => x -> x .== "bayes_factor_10",
                )
            end

            for j in 1:2
                if ax_counter == 1
                    #ylabel = L"\text{accuracy}"
                    ylabel = ""
                else
                    ylabel = ""
                end

                ax = Axis(
                        fig3[row,ax_counter], 
                        xscale = Makie.pseudolog10,
                        ylabel = ylabel,
                        xgridvisible = false,
                        ygridvisible = false,
                        topspinevisible = false,
                        bottomspinevisible = true,
                        leftspinevisible = ax_counter == 1,
                        rightspinevisible = false,
                        xticks = (xt, xtl),
                )

                if ax_counter > 1
                    hideydecorations!(ax)
                end

                hidexdecorations!(ax, ticks = false)
                push!(axs, ax)

                if j == 1
                    criterion = x -> x .> 0
                    color = "gray"
                else
                    criterion = x -> x .== 0
                    color = "orange"
                end

                df2 = subset(
                    df1,
                    :N_true_sum => criterion 
                )


                y = df2[!,metric]
                y = y[.!isnan.(y)]

                println(length(y))

                if length(y) > 20
                    hist!(ax, y, direction = :x, color = color)
                    #mean_acc = string.(mean(y) * 100)[1:4] * "%"
                    #text!(ax, 2.0, 1.0, text = "mean = $mean_acc")
                    #text!(ax, 1.0, textpos[row], text = "$mean_acc")
                end

                ax_counter += 1
                if ax_counter == 9
                    ax_counter += 1
                end
            end
        end

    end
    linkaxes!(axs...)

end

#label = LaTeXString(string(Int64(tree_heights[q])) * raw" Ma"))
Label(fig3[0, 1:8], L"\text{a) split by tree size}")
Label(fig3[4, 1:2], L"(2,50]")
Label(fig3[4, 3:4], L"(50,250]")
Label(fig3[4, 5:6], L"(250,1000]")
Label(fig3[4, 7:8], L"(1000,∞)")
Label(fig3[5, 1:8], L"\text{number of taxa}")

Label(fig3[0, 9:17], L"\text{b) split by model}")
Label(fig3[4, 10:11], L"\text{tiny}")
Label(fig3[4, 12:13], L"\text{small}")
Label(fig3[4, 14:15], L"\text{moderate}")
Label(fig3[4, 16:17], L"\text{large}")
Label(fig3[5, 9:17], L"\text{allowed rate variation}")

Label(fig3[1, 0], L"\text{accuracy}", rotation = π/2)
Label(fig3[2, 0], L"\text{false positive ratio}", rotation = π/2)
Label(fig3[3, 0], L"\text{false negative ratio}", rotation = π/2)

rowgap!(fig3.layout, 15)
colgap!(fig3.layout, 3)

for row in 1:3
    rowsize!(fig3.layout, row, Relative(0.3))
end

groups = [PolyElement(color = color, strokecolor = :transparent) for color in [:gray, :orange]]
Legend(fig3[1,5:8], groups, 
       [L"\geq 1 \text{ rate shifts}", L"\text{no rate shifts}"],
      )

fig3
CairoMakie.save("figures/confusion-split.pdf", fig3)


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
            value = mean(y)


            summary_df = DataFrame(
                                :value => value,
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
            value = mean(y)


            summary_df = DataFrame(
                                :value => value,
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
    barplot!(ax, 
             df_no_shift[!,:tree_size_bin] .+ 0.25, 
             df_no_shift[!,:value],
             color = Pattern("/", background_color = :white, linecolor = :black),
             strokecolor = :black,
             strokewidth = 1,
             width = 0.5,
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
    barplot!(ax, 
             df_no_shift[!,:model] .+ 0.25, 
             df_no_shift[!,:value],
             color = Pattern("/", background_color = :white, linecolor = :black),
             strokecolor = :black,
             strokewidth = 1,
             width = 0.5,
            )
end

ylims!(axs[1], (0.96, 1.01))
ylims!(axs[4], (0.96, 1.01))

ylims!(axs[3], (0.85, 1.01))
ylims!(axs[6], (0.85, 1.01))

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
CairoMakie.save("figures/confusion-split-simple.pdf", fig4)
         
