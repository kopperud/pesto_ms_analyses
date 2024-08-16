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
    L"\text{tiny var.}",
    L"\text{small var.}",
    L"\text{moderate var.}",
    L"\text{large var.}"
]

yt = lrange(0.25, 4.0, 5)
fig = Figure(size = (650, 620), fontsize = 14, 
            figure_padding = (1,1,1,1))
axs = []

####################
##
## trees without any rate shifts
##
###################
for (i, inference) in enumerate(inferences)
    for j in 1:n_models
        if i == 1
            title = titles[j]
        else
            title = ""
        end
        ax = Axis(fig[i+1,j+1], 
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
            df4 = subset(
                df,
                :criterion => x -> x .== "N_over_half",
                :inference => x -> x .== inference,
                :height =>    x -> x .== h,
                :model  =>    x -> x .== j,
                :N_true_sum => x -> x .== 0,
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

####################
##
## trees with at least one rate shifts
##
###################
for (i, inference) in enumerate(inferences)
    for j in 1:n_models
    
        ax = Axis(fig[i+5,j+1], 
                    xscale = Makie.log10,
                    yscale = Makie.log10,
                    xgridvisible = false,
                    ygridvisible = false,
                    topspinevisible = false,
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
ylabel = Label(fig[2:8, 1], L"\text{proportional error (speciation rate)}", rotation = π/2)
title1 = Label(fig[1, 2:5], L"\text{a) Simulated trees without rate shifts}")
title2 = Label(fig[5, 2:5], L"\text{b) Simulated trees with }\geq \text{1 rate shifts}", justification = :bottom)

xlabel = Label(fig[9, 2:5], L"\text{number of taxa}")


#xlabel = Label(fig[1, 2:5], L"\text{allowed rate variation}")
ylabel = Label(fig[2, 6], L"\text{true}~𝛌,𝛍,\eta", rotation = π/2)
ylabel = Label(fig[3, 6], L"\text{true}~𝛌,𝛍", rotation = π/2)
ylabel = Label(fig[4, 6], L"\text{unknown}~𝛌,𝛍,\eta", rotation = π/2)

ylabel = Label(fig[6, 6], L"\text{true}~𝛌,𝛍,\eta", rotation = π/2)
ylabel = Label(fig[7, 6], L"\text{true}~𝛌,𝛍", rotation = π/2)
ylabel = Label(fig[8, 6], L"\text{unknown}~𝛌,𝛍,\eta", rotation = π/2)




labels = [
    L"\text{25 Ma}",
    L"\text{50 Ma}",
    L"\text{75 Ma}",
    L"\text{100 Ma}",
    L"\text{125 Ma}"
]
title = "Shift in"
fig[4:6, 7] = Legend(fig, axs[1], "", framevisible = false, patchsize = (30, 30))
colsize!(fig.layout, 7, Relative(0.2))
colsize!(fig.layout, 1, Relative(0.04))
colgap!(fig.layout, 7)
rowgap!(fig.layout, 7)

#rowsize!(fig.layout, 0, Relative(0.05))
for row in [1,5,9]
    rowsize!(fig.layout, row, Relative(0.04))
end

for i in [2,3,4,6,7,8]
    rowsize!(fig.layout, i, Relative(0.14667))
end
#rowsize!(fig.layout, 4, Relative(0.05))
fig

CairoMakie.save("figures/proportional-error-lambda-split.pdf", fig)

