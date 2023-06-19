module CSSim

using DynamicalSystems

hello_world() = println("Hello world!")

include("common.jl")
include("one_forest_system.jl")
include("two_forest_system.jl")
include("n_forest_system.jl")

end # module CSSim
