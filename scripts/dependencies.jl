import Pkg

Pkg.add("Distributions")
#Pkg.add("StatsPlots")
#Pesto")
#Pkg.add("Optim")
#Pkg.add("BenchmarkTools")
Pkg.add("ForwardDiff")
Pkg.add("ProgressMeter")
Pkg.add("RCall")
Pkg.build("RCall")
#Pkg.add("Measures")
Pkg.add("JLD2")
Pkg.add("Glob")
Pkg.add("CSV")
Pkg.add("Revise")
Pkg.add("DataFrames")
Pkg.add(url="https://github.com/kopperud/Pesto.jl")

exit()
