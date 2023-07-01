module CSSim

using DynamicalSystems
using DrWatson

# Functions common to models 
export α, w, γ, εₖ, B, θ, ecosystems_times_to_deforest
export iskilled
export EcosystemDeforestTime

# Models and model jacobians
export one_forest_system
export two_forest_system 
export n_forest_rule!
export n_forest_system
export n_forest_system_oop
export n_forest_tds
export antonovsky_jacob, two_forest_jacob, n_forest_jacob!

# Plotting 
export phase_portrait

# Simulation
export deforested_position_sim

include("common.jl") 
include("n_forest.jl")
include(joinpath(scriptsdir(), "plots.jl"))
include(joinpath(scriptsdir(), "deforested_position_sim.jl"))

end # module CSSim
