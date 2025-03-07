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

    #states = Int64[] 
    states = zeros(Int64, n_times)
    #push!(states, rand(d))
    states[1] = rand(d)

    ## Initialize first so there is less memory allocation
    N = Int64[0]
    F = zeros(4)
    dF = zeros(4)
    P = zeros(4)

    for step in 1:n_steps
        current_time = t1 + (step-1) * Δt
        current_state = states[step]
        
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

        P[:] = F .* D(current_time + Δt)[:,2]
        P[:] = P ./ sum(P)

        new_state_d = Categorical(P)
        new_state = rand(new_state_d)

        if new_state != current_state
            N[1] += 1
        end

        #push!(states, new_state)
        states[step+1] = new_state
    end

    return(states, N[1])
end

## wrapper that does replicates
function simulate_nshifts_branch(edge_index, Ds, Fs, E, n, n_replicates)
    N = zeros(Int64, n_replicates)
    states = zeros(Int64, n, n_replicates)

    D = Ds[edge_index]
    F = Fs[edge_index]

    t0, t1 = extrema(D.t)

    St = Pesto.ancestral_state_probability(D(t1)[:,2], F(t1)[:,2], t1)
    println(St)

    for i in 1:n_replicates
    #Threads.@threads for i in 1:n_replicates
        states1, N1 = simulate_nshifts_branch(D, St, E, n)

        N[i] = N1
        states[:,i] = states1

    end

    return(states, N)
end

n_time_steps = 100 ## not so important
n_replicates = 30_000 ## important 
states, N = simulate_nshifts_branch(304, Ds, Fs, E, n_time_steps, n_replicates)

## variance per time slice/bin
sim_vars = zeros(n_time_steps)
for i in 1:n_time_steps
    sim_vars[i] = var([model.λ[states[i,j]] for j in 1:n_replicates])
end
sort(sim_vars)

D = Ds[304]
F = Fs[304]

theoretical_vars = zeros(n_time_steps)
times = collect(range(D.t[end], D.t[1]; length = n_time_steps))
for i in 1:n_time_steps
    St = Pesto.ancestral_state_probability(D(times[i])[:,2], F(times[i])[:,2], 0.0)
    first_moment = sum(model.λ .* St)
    second_moment = sum(model.λ .* model.λ .* St)
    theoretical_vars[i] = second_moment - first_moment ^2
end

fig = Figure(size = (800, 350));
ax1 = Axis(fig[1,1], xlabel = "simulated Var[X(t)]", ylabel = "theoretical Var[X(t)]")
ax2 = Axis(fig[1,2], xlabel = "time along branch (Ma)", ylabel = "abs(Var (theoretical) - Var (simulated))")
ax3 = Axis(fig[1,3], xlabel = "time along branch (Ma)", ylabel = "Var[X(t)]")
lines!(ax1, sim_vars, theoretical_vars, label = "Var[X(t)]")
lines!(ax2, times, abs.(sim_vars .- theoretical_vars))
lines!(ax3, times, sim_vars, label = "simulated")
lines!(ax3, times, theoretical_vars, label = "theoretical")
ax2.xreversed = true
ax3.xreversed = true
axislegend(ax1, position = :lt)
axislegend(ax3, position = :lt)

fig

means, vars, prob = brvar(model, primates, 304, 8)
sum(vars .* prob)

## variance of the branch-specific mean
sim_vars_branch = zeros(n_replicates)
for i in 1:n_replicates
    branch_mean = mean([model.λ[states[j,i]] for j in 1:n_time_steps])
    sim_vars_branch[i] = branch_mean
end
sim_vars_branch

theoretical_vars_branch = 0.0
for j in 1:n_time_steps
    t = times[j]
    St = Pesto.ancestral_state_probability(D(t)[:,2], F(t)[:,2], 0.0)
    ## 1st moment
    m1 = sum(model.λ .* St)
    ## 2nd moment
    m2 = sum(model.λ .* model.λ .* St)
    
    ## var (assuming independence)
    v = m2 - m1^2
    theoretical_vars_branch += v / (n_time_steps^2)
end
theoretical_vars_branch
var(sim_vars_branch)

var(sim_vars_branch) / theoretical_vars_branch
theoretical_vars_branch / var(sim_vars_branch)

## calculate covariance
t1 = 30.0
t2 = 25.0
S1 = Pesto.ancestral_state_probability(D(t1)[:,2], F(t1)[:,2], 0.0)
S2 = Pesto.ancestral_state_probability(D(t2)[:,2], F(t2)[:,2], 0.0)

E1 = sum(S1 .* model.λ)
E2 = sum(S2 .* model.λ)

sum([(model.λ[i] - E1)*(model.λ[j] - E2)*S1[i]*S2[j] for i in 1:4, j in 1:4])



n_steps = [10, 25, 50, 100, 150, 300, 500, 1000, 2500, 5000, 10_000]

n_replicates = 10_000

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

#save("figures/stochastic_mapping_vs_differential_equation.pdf", fig)
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



