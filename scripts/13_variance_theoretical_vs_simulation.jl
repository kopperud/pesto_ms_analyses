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

model = BDSconstant(λ, μ, η)
pesto_rates = birth_death_shift(model, primates)
edge_index = 304

Ds, Fs = backwards_forwards_pass(model, primates);
E = extinction_probability(model, primates);
Ss = ancestral_state_probabilities(Ds, Fs);

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

n_time_steps = 300 ## not so important
n_replicates = 50_000 ## important 
states, N = simulate_nshifts_branch(304, Ds, Fs, E, n_time_steps, n_replicates)

## variance per time slice/bin
sim_means = zeros(n_time_steps)
sim_vars = zeros(n_time_steps)
for i in 1:n_time_steps
    sim_vars[i] = var([model.λ[states[i,j]] for j in 1:n_replicates])
    sim_means[i] = mean([model.λ[states[i,j]] for j in 1:n_replicates])
end
sort(sim_vars)

D = Ds[304]
F = Fs[304]

theoretical_vars = zeros(n_time_steps)
theoretical_means = zeros(n_time_steps)

times = collect(range(D.t[end], D.t[1]; length = n_time_steps))
for i in 1:n_time_steps
    St = Pesto.ancestral_state_probability(D(times[i])[:,2], F(times[i])[:,2], 0.0)
    first_moment = sum(model.λ .* St)
    second_moment = sum(model.λ .* model.λ .* St)
    theoretical_means[i] = first_moment 
    theoretical_vars[i] = second_moment - first_moment ^2
end

fig = Figure(size = (650, 280));

ax1 = Axis(fig[1,1], 
    xlabel = L"\text{time before the present (Ma)}",
    ylabel = L"\mathbb{E}[\lambda(t)]",
    xgridvisible = false, 
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
    ) 
ylims!(ax1, (0.0, maximum(theoretical_means)+0.1))
    
ax2 = Axis(fig[1,2], 
    xlabel = L"\text{time before the present (Ma)}",
    ylabel = L"\text{Var}[\lambda(t)]",
    xgridvisible = false, 
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
    ) 

lines!(ax1, times, sim_means, color = :black, label = L"\text{simulated}")
lines!(ax1, times, theoretical_means, color = :black, linestyle = :dash, label = L"\text{theoretical}")
ax1.xreversed = true

lines!(ax2, times, sim_vars, color = :black, label = L"\text{simulated}")
lines!(ax2, times, theoretical_vars, color = :black, linestyle = :dash, label = L"\text{theoretical}")
ax2.xreversed = true
axislegend(ax2, position = :lt)

fig

#save("figures/speciation_rate_variance.pdf", fig)

function covariance(t1, t2, model, D, F)
    D1 = D(t1)[:,2]
    D2 = D(t2)[:,2]

    F1 = F(t1)[:,2]
    F2 = F(t2)[:,2]

    S1 = D1 .* F1
    S1 = S1 / sum(S1)

    S2 = D2 .* F2
    S2 = S2 / sum(S2)

    P = Pesto.Pmatrix_precise(model, D, E, t1, t2)

    E1 = sum(model.λ .* S1)
    E2 = sum(model.λ .* S2)

    ## covariance
    co = 0.0
    for i in 1:4
        for j in 1:4
            co += S1[j] * P[i,j] * (model.λ[j] - E1) * (model.λ[i] - E2)
        end
    end
    return(co)
end

l = model.λ[states]
covariance(times[1], times[2], model, D, F)
cov(l[1,:], l[2,:])

Pesto.Pmatrix_precise(model, D, E, 32.0, 31.0)

Δt = times[2] - times[1]
P = Pmatrix(model, D, E, times[1], times[60] - times[1])
P = Pesto.Pmatrix_precise(model, D, E, times[1], times[60])

theoretical_covars = Float64[]
times2 = collect(range(D.t[end], D.t[1]; length = 300))
for i in 2:length(times2) 
    push!(theoretical_covars, covariance(times2[1], times2[i], model, D, F))
end
simulation_covars = Float64[]
for i in 2:n_time_steps
    push!(simulation_covars, cov(l[1,:], l[i,:]))
end
#scatter(1.0, 1.0)

fig = Figure()
ax = Axis(fig[1,1], xlabel = "time along branch", ylabel = "cov(x(0), x(t1))")
lines!(ax, times2[2:end], theoretical_covars, label = "theoretical")
lines!(ax, times[2:end], simulation_covars, label = "simulation")
axislegend(ax, position = :lb)
ax.xreversed = true
fig


function foobaz(c::Int64)
    
    var_barlambda = 0.0
    times2 = collect(range(D.t[end], D.t[1]; length = c))

    for a in 1:c
        t = times2[a]
        
        Dt = D(t)[:,2]
        Ft = F(t)[:,2]
        St = Dt .* Ft
        St = St / sum(St)

        first_moment = sum(model.λ .* St)
        second_moment = sum(model.λ .* model.λ .* St)

        va = second_moment - first_moment^2
        var_barlambda += (1/c^2) * va

        for b in 2:c
            if b > a
                var_barlambda += 2*(1/c^2) * covariance(t, times2[b], model, D, F)
            end
        end
    end
    return(var_barlambda)
end

using ProgressMeter

function foobaz2(times, states, model)
    var_barlambda = 0.0
    l = model.λ[states]

    n_steps = length(times)


    @showprogress for a in 1:n_steps
        Et = mean(l[a,:]) / n_steps
        Vt = var(l[a,:]) / (n_steps^2)

        var_barlambda += Vt
    end
    return(var_barlambda)
end

foobaz2(times, states, model)
foobaz(300)

ns = collect(range(2, 50))
vars = [foobaz(n) for n in ns]

fig2 = Figure(size = (350, 250))
ax = Axis(fig2[1,1], 
    xlabel = L"\text{number of time bins}", 
    ylabel = L"\text{Var}[\bar{\lambda}]",
    xgridvisible = false, 
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
    )
lines!(ax, ns, vars, label = "theoretical", color = :black)
#scatter!(ax, ns, vars, color = :black)
vs = var(mean(l, dims = 1)[1,:])
lines!(ax, [extrema(ns)...], [vs, vs], label = "simulation", linestyle = :dash, color = :black)
axislegend(ax, position = :rb)
fig2

save("figures/speciation_rate_branchvariance.pdf", fig2)




