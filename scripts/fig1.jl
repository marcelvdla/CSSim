# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.6
#   kernelspec:
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# +
# Produce figure 1, phase portrait of Antonovsky model from Cantin2020
using DrWatson
using Distributions
using Random
using Plots
using ChaosTools
@quickactivate :CSSim
include(joinpath(srcdir(), "one_forest_system.jl"))

# +
# Define parameters
ρ = 4.2
f = 1
h = 2

# Determine number of timesteps and initial states
T = 100
n_initial_points = 50
n_states = 2

# Create initial states
rng = MersenneTwister(42)
u0s = rand(rng, Uniform(0, 5), (n_initial_points, n_states))
# -

# symbolic jacobian for antonvosky model
antonovsky_sym_jacob()

# +
# Numerically calculate fixed points
low = interval(0, 5)
high = interval(0, 5)
box = low × high
ds = one_forest_system(u0s[1, :]; ρ = ρ, f = f, h = h)

fps, eigs, stable = fixedpoints(ds, box, antonovsky_jacob)
# -

# Iterate through initial states and create phase portrait
nrows = size(u0s)[1]
P = plot(title = "Phase Portrait of Antonovsky Model", legend = false)
for i in 1:nrows
    ds = one_forest_system(u0s[i,:], ρ = ρ, f = f, h = h)
    trajectory_ds = trajectory(ds, T)
    X, _ = trajectory_ds
    young_tree_density, old_tree_density = X[:, 1], X[:, 2]
    plot!(P, young_tree_density, old_tree_density, color = :grey)
end
scatter!([fps[1, 1]], [fps[1, 2]], markersize = 7, markercolor=:black)
scatter!([fps[2, 1]], [fps[2, 2]], markersize = 7, markercolor=:black)
scatter!([fps[3, 1]], [fps[3, 2]], markersize = 7, markercolor=:black)
display(P)


