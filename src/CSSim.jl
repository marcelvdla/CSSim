module CSSim

using DynamicalSystems
using DrWatson

# Functions common to models 
export α, w, γ, εₖ, B, θ, ecosystems_times_to_deforest, EcosystemDeforestTime

# Models and model jacobians
export one_forest_system, antonovsky_jacob
export two_forest_system, two_forest_jacob
export n_forest_system, n_forest_system_oop, n_forest_tds, n_forest_jacob!

# Plotting 
export phase_portrait

include("common.jl") 
include("n_forest.jl")
include("chaos.jl")
include(joinpath(scriptsdir(), "plots.jl"))

end # module CSSim
