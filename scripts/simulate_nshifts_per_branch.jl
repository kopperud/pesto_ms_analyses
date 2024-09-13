using Revise
using Pesto
using DataFrames
using CSV
using Statistics
using RCall
using CairoMakie
using LaTeXStrings
using Distributions

## Pesto Inference
treefile = "data/primates.tre"
phy = Pesto.readtree(treefile)
num_total_species = 367
ρ = length(phy.tip_label) / num_total_species
primates = SSEdata(phy, ρ)

λ = [0.1, 0.2, 0.4, 0.1]
μ = [0.05, 0.15, 0.10, 0.2]
η = optimize_eta(λ, μ, primates)

model = SSEconstant(λ, μ, η)
pesto_rates = birth_death_shift(model, primates)


edge_index = 304


Ds, Fs = backwards_forwards_pass(model, primates);
E = extinction_probability(model, primates);
Ss = ancestral_state_probabilities(Ds, Fs);


analytical_shift_prob = posterior_shift_prob(model, primates)[304]



function simulate_nshifts_branch(D, St, E, n_times)
    t0, t1 = extrema(D.t)
    
    n_steps = n_times-1
    Δt = -(t1 - t0)/n_steps
   
    d = Categorical(St)

    states = Int64[] 
    push!(states, rand(d))

    ## Initialize first so there is less memory allocation
    N = Int64[0]
    F = zeros(4)
    dF = zeros(4)
    P = zeros(4)

    for step in 1:n_steps
        current_time = t1 + (step-1) * Δt
        current_state = states[end]
        
        F[:] .= 0
        F[current_state] = 1
        r = η/3
       
        ## a few Euler's steps
        tiny_steps = 2
        for j in 1:tiny_steps
            ΔΔt = Δt / tiny_steps
            Et = E(current_time + ΔΔt * (j-1))
            dF[:] = (-1) .* (-(λ .+ μ .+ η) .* F .+ 2 .* λ .* F .* Et .+ r .* (sum(F) .- F))
            F += ΔΔt .* dF
        end

        P[:] = F .* D(current_time + Δt)
        P[:] = P ./ sum(P)

        new_state_d = Categorical(P)
        new_state = rand(new_state_d)

        if new_state != current_state
            N[1] += 1
        end

        push!(states, new_state)
    end

    return(N[1])
end

## wrapper that does replicates
function simulate_nshifts_branch(edge_index, Ds, Fs, E, n, n_replicates)
    N = zeros(Int64, n_replicates)

    D = Ds[edge_index]
    F = Fs[edge_index]

    t0, t1 = extrema(D.t)

    St = Pesto.ancestral_state_probability(D(t1), F(t1), t1)

    for i in 1:n_replicates
    #Threads.@threads for i in 1:n_replicates
        N[i] = simulate_nshifts_branch(D, St, E, n)
    end

    return(N)
end



n_steps = [10, 25, 50, 100, 150, 300, 500, 1000, 2500, 5000, 10_000]

n_replicates = 100_000

N = zeros(Int64, length(n_steps), n_replicates)
for (i, n) in enumerate(n_steps)
    print(".")
    #N[i,:] = [simulate_nshifts_branch(304, Ds, Ss, E, n) for _ in 1:n_replicates]
    res = simulate_nshifts_branch(304, Ds, Fs, E, n, n_replicates)
    N[i,:] .= res
end

fig = Figure();
xtl = ["$n" for n in n_steps]
ax = Axis(fig[1,1], 
    xlabel = "number of steps",
    ylabel = "posterior prob ≥1 rate shifts",
    title = "ODE vs stochastic map (1 branch, 100,000 runs)",
    xgridvisible = false, 
    ygridvisible = false,
    xscale = log10,
    xticks = (n_steps, xtl)
)
p_rate_shift = sum(N .> 0, dims = 2)[:,1] ./ n_replicates
sc = scatter!(ax, 
    n_steps,
    p_rate_shift,
    label = "stochastic mapping"
)
li = CairoMakie.lines!(
    ax, [extrema(n_steps)...], [analytical_shift_prob, analytical_shift_prob],
    linestyle = :dash,
    color = :gray,
    label = "differential equation"
)
axislegend(ax, position = :rt)
ylims!(ax, (0.96, 0.98))
fig

save("figures/stochastic_mapping_vs_differential_equation.pdf", fig)
fig



p = Pesto.posterior_shift_prob(model, primates)[304]

Nhat = pesto_rates[304,:nshift]
#(1.0 - p)^2 + (0.0 - (1-p))^2
v = (1.0 - Nhat)^2 * p  + (0.0 - Nhat)^2 * (1-p)
sqrt(v)

p

hist(df[!,"num_shifts[461]"])

for i in 1:5
    y = df[!, "num_shifts[461]"]
    r = (y .== i) / size(df)[1]
    println(y .== i)
end

px = [sum(df[!, "num_shifts[461]"] .== i) / size(df)[1] for i in 0:5]

px
sqrt(sum([(i -1 - Nhat)^2*px[i] for i in 1:6]))
sqrt(sum([(i -1 - Nhat)^2*px[i] for i in 1:2]))


## variance of N

#N grouped by departure category j
function ode(du, u, p, t)
    η, K, D, F = p

    Dt = D(t)
    Ft = F(t)
    St = Pesto.ancestral_state_probability(Dt, Ft, t)

    Dsum = sum(Dt)
    r = -η / (K- 1.0)

    du[:] = r .* St .* (Dsum .- Dt) ./ Dt
    #for j in 1:K
        #du[j] = r * (Dsum - Dt[j]) / Dt[j]
    #end
end

using OrdinaryDiffEq
D = Ds[edge_index]
F = Fs[edge_index]

K = length(model.λ)
t0, t1 = extrema(D.t);
tspan = (t1, t0)
params = (model.η, K, D, F);

u0 = [0.0, 0.0, 0.0, 0.0]
prob = OrdinaryDiffEq.ODEProblem(ode, u0, tspan, params)

St = Pesto.ancestral_state_probability(D(t1), F(t1), t1)
sol = solve(prob, Tsit5())
N = sol.u[end]

X = Pesto.posterior_shift_prob_categories(model, D, K, Pesto.no_shifts_problem(model), Tsit5())

p_no_shift = X
p_one_shift = 1 .- X

v = (1 .- N).^2 .* p_one_shift .+ (0 .- N).^2 .* p_no_shift

sqrt(sum(v .* St))

px = [sum(rb_shifts[:,461] .== i-1) / size(rb_shifts)[1] for i in 1:6]

mu = mean_shifts[461]
v_approx = [(0 - mu)^2 * px[1], (1 - mu)^2 * (1 - px[1])]
std_approx = sqrt(v_approx)

v_numerical = sum([(i-1 - mean_shifts[461])^2 * px[i] for i in 1:6])
std_numerical = sqrt(v_numerical)

std_shifts[461]
v


sol.u[end]

sol.u[end] .* St



