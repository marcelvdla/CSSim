# Produce figure 1, phase portrait of Antonovsky model from Cantin2020

using DrWatson
using Distributions
using Random
@quickactivate :CSSim

# Define parameters
ρ = 4.2
f = 1
h = 2

# Determine number of timesteps and initial states
T = 100
n_initial_points = 100
n_states = 2

# Create initial states
rng = MersenneTwister(42)
u0s = rand(rng, Uniform(0, 5), (n_initial_points, n_states))

# Iterate through initial states and create phase portrait
# TODO: finish
for i in 1:size(arr)[1]
    ds = one_forest_system(u0s[i,:], ρ = ρ, f = f, h = h)
    trajectory_ds = trajectory(ds, T, )
end 


