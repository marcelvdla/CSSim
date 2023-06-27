module CSSim

using DynamicalSystems
using DrWatson

hello_world() = println("Hello world!")

# include("common.jl")  # 
include("n_forest.jl")
include(joinpath(scriptsdir(), "plots.jl"))

end # module CSSim
