using CSV
using DataFrames
using CairoMakie
using LaTeXStrings

function lrange(from::Float64, to::Float64, length::Int64 = 6)
    exp.(collect(range(log(from), log(to); length = length)))
end


df = CSV.read("output/branch_specific_estimation_error.csv", DataFrame,
                missingstring=["NA", "NAN", "NULL"])

df = subset(
    df, 
    :inference => x -> x .== "unknown_rates",
    :criterion => x -> x .== "N_half_and_bayes_factor"
    )



fig = CairoMakie.Figure(size = (400, 400))

heights = [25, 50, 75, 100, 125]

dfs = [
    subset(df, :height => x -> x .== height) for height in heights
]


function proportion_shifts(df, ntaxa_lower, ntaxa_upper)
    df = deepcopy(df)
    df = subset(df,
        :ntaxa =>  x -> x .> ntaxa_lower,
    )
    df = subset(df,
        :ntaxa =>  x -> x .<= ntaxa_upper,
    )

    p = sum(df[!,:N_true_sum] .> 0) ./ size(df)[1]
    return(p)
end


props = zeros(4)
props[1] = proportion_shifts(df, 0.0, 50.0)
props[2] = proportion_shifts(df, 50.0, 250.0)
props[3] = proportion_shifts(df, 250.0, 1000.0)
props[4] = proportion_shifts(df, 1000.0, 1e20)

## of the treed that experienced at least one shift, how many ?
N_true_distribution = subset(df, :N_true_sum => x -> x .> 0)[!,:N_true_sum]
mean(N_true_distribution)
median(N_true_distribution)
hist(N_true_distribution)





props