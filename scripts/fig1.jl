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
nrows = size(u0s)[1]
for i in 1:nrows
    ds = one_forest_system(u0s[i,:], ρ = ρ, f = f, h = h)
    trajectory_ds = trajectory(ds, T)
    young_tree_density_x, old_tree_density_y = trajectory_ds
    
    break
end

