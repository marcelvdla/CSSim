module CSSim

using DynamicalSystems
using DrWatson

export α, w, γ, εₖ, B, θ, ecosystems_times_to_deforest, EcosystemDeforestTime
export one_forest_system, antonovsky_jacob
export two_forest_system, two_forest_jacob
export n_forest_system, n_forest_system_oop

hello_world() = println("Hello world!")

include("common.jl") 
include("n_forest.jl")
include("fixedpoints.jl")
include(joinpath(scriptsdir(), "plots.jl"))

end # module CSSim
