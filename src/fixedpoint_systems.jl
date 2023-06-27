# Rules for n-forests but using state vectors instead of state matrices
using DynamicalSystems
using UnPack

include("common.jl")

"""
    n_forest_rule(u::SVector, params::Dict{Symbol, Any}, t)

OOP `n`-forest rule for biotic pump system with `u` as a vector of length 
`n * 2` such that `[x1, x2, ..., xn, y1, y2, ..., yn]`
"""
function n_forest_rule(u::SVector, params::Dict{Symbol, Any}, t) 
    d::Union{Vector{Any}, Vector{Int}} = [] # distances vector
    @unpack ρ, f, a₁, h, a₂, d, l, 
            α₀, w₀, P₀, β₁, β₂, 
            n  = params
    @assert n >= 1 "At least 1 forest ecosystem"
    @assert length(d) == (n-1) "n-1 distances in distance vector `d`"

    n_states = 2
    n = length(u) ÷ n_states

    # Get all x and y state variables
    x = u[1:n]
    y = u[n+1:end]

    # For storing the computed time derivatives
    du = zeros(n*n_states)

    for i in 1:n 
        # Get state variables 
        xᵢ = x[i]
        yᵢ = y[i]

        # Biotic pump
        wᵢ = n > 1 ? w(i, x, y, d, l, P₀, β₁, β₂) : 0
        αᵢ = n > 1 ? α(wᵢ, α₀, w₀) : 0

        # Compute time derivatives
        ẋᵢ = ρ*yᵢ - γ(yᵢ)*xᵢ - f*xᵢ + a₁*αᵢ*xᵢ
        ẏᵢ = f*xᵢ - h*yᵢ + a₂*αᵢ*yᵢ

        # update time derivative vector
        du[i] =  ẋᵢ
        du[i+n] = ẏᵢ
    end 

    return SVector{n*n_states}(du)
end 

"""
    n_forest_jacob 

1D jacobian. 
"""
function n_forest_jacob(J, u, params::Dict{Symbol, Any}, t)

end 

"""
    n_forest_system_oop(u0, params::Dict{Symbol, Any})

# Examples 
```julia-repl
julia> n = 2
julia> params = Dict{Symbol,Any}(
        :ρ  => 4.2, 
        :f  => 1.0,
        :α₀ => -1.0, 
        :w₀ =>  1.0,
        :a₁ => 1.0, 
        :h  => 2.0, 
        :a₂ => 0.0, 
        :d  => [42], 
        :l  => 600.0, 
        :P₀ => 1.0, 
        :β₁ => 0.0, 
        :β₂ => 1.0,
        :n  => n,
        :ecosystems_to_deforest => [],
    )
julia> u0 = 4*rand(n, 2)
julia> ds = n_forest_system_oop(u0, params)
```
"""
function n_forest_system_oop(u0, params::Dict{Symbol, Any})
    return CoupledODEs(n_forest_rule, u0, params)
end 

""" 
    n_forest_rule(u::SVector{<:Real}, params::Dict{Symbol, Any}, t)

OOP rule for `n`-forest system.

The vector `u` is a vector of length `n * 2` such that 
`[x1, x2, y1, y2, ..., xi, xi+1, yi, yi+1]` for i in 1 to n-1.
function n_forest_rule(u, params::Dict{Symbol, Any}, t)
    # distance between `i` and `i+1` forest vector
    d::Union{Vector{Any}, Vector{Int}} = []

    ecosystems_to_deforest::Union{
        Vector{Any}, Vector{EcosystemDeforestTime}} = []

    @unpack ρ, f, a₁, h, a₂, d, l, 
            α₀, w₀, P₀, β₁, β₂, 
            ecosystems_to_deforest,
            n  = params
    @assert n >= 1 "At least 1 forest ecosystem"
    @assert length(d) == (n-1) "n-1 distances in distance vector `d`"

    # two state variable indices
    x_ix = 1
    y_ix = 2

    # reshape the input vector
    n_states = 2
    u = reshape(u, n, n_states)

    # Make a time derivative matrix
    du = zeros(n, n_states)

    # x and y for all ecosystems
    x = u[:, x_ix]
    y = u[:, y_ix]

    # For all ecosystems, define the linear complex network of dynamical systems
    ecosystem_ids_to_deforest = getproperty.(ecosystems_to_deforest, :i)
    tstars = getproperty.(ecosystems_to_deforest, :tstar)
    for i in 1:n 
        xᵢ = x[i]
        yᵢ = y[i]

        # Biotic pump 
        wᵢ = n > 1 ? w(i, x, y, d, l, P₀, β₁, β₂) : 0
        αᵢ = n > 1 ? α(wᵢ, α₀, w₀) : 0

        # Possible deforestation of ecosystem i by finding the index of the first 
        # deforested ecosystem id that matches the current ecosystem id 
        # and then indexing the parallel `tstars` list
        εᵢ = εₖ(i, ecosystem_ids_to_deforest)
        tstar_ix = findfirst(id -> id == i, ecosystem_ids_to_deforest)
        tstar = εᵢ ? tstars[tstar_ix] : 0
        θᵢ = εᵢ ? θ(t, tstar) : 0

        # forest dynamical system for `i^th` ecosystem
        du[i, x_ix] = ρ*yᵢ - γ(yᵢ)*xᵢ - f*xᵢ + a₁*αᵢ*xᵢ - εᵢ*θᵢ*xᵢ
        du[i, y_ix] = f*xᵢ - h*yᵢ + a₂*αᵢ*yᵢ - εᵢ*θᵢ*yᵢ
    end 

    return SVector{n*n_states}(du)
end 
"""


"""
    n_forest_rule(u, params::Dict{Symbol, Any}, t)

DEPRECATED

OOP `n_forest_rule` using [`n_forest_rule!`](@ref).

function n_forest_rule(u, params::Dict{Symbol, Any}, t)
    n_states = 2
    n = params[:n]
    du = zeros(n, n_states)

    # The in-place rule requires an abstract matrix
    u_reshaped = reshape(u, n, n_states)
    n_forest_rule!(du, u_reshaped, params, t)

    # return a static vector and permute the dims of du since
    # du is of the shape [n, 2] such that [x1 y1; ... ; xi yi]
    # so e.g., [x1 x2; y1 y2] --> [x1 y1; x2 y2]
    return SVector{n*n_states}(permutedims(du))
end 
"""
