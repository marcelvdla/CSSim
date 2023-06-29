# Functions for the n-forest model with/without biotic/deforestation components
using DynamicalSystems
using Symbolics
using UnPack
using CSSim

"""
    n_forest_system(u0, params::Dict{Symbol, Any})
"""
function n_forest_system(u0, params::Dict{Symbol, Any})
    return CoupledODEs(n_forest_rule!, u0, params)
end

"""
    n_forest_system(
        u0::AbstractMatrix, params::Dict{Symbol, Any}, diffeq::NamedTuple)

`n_forest_system` with parameter `diffeq` for solver arguments.
"""
function n_forest_system(
    u0::AbstractMatrix, params::Dict{Symbol, Any}, diffeq::NamedTuple)
    return CoupledODEs(n_forest_rule!, u0, params; diffeq=diffeq)
end 

"""
    n_forest_rule!(du, u::AbstractMatrix, params::Dict{Symbol, Any}, t)

[1] : Equation (10) from Cantin2020
"""
function n_forest_rule!(du, u::AbstractMatrix, params::Dict{Symbol, Any}, t)
    # distance between `i` and `i+1` forest vector
    d::Union{Vector{Any}, Vector{<:Real}} = []

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

    return nothing
end 

"""
    n_forest_sym_rule(n::Int)

Return symbolic repr. for `n`-forest complex network without perturbations.

[1] : Inspiration from the implementation of the Jacobian for the Lorenz96
    system based on "Ordinal pattern-based complexity analysis of 
    high-dimensional chaotic time series" (https://doi.org/10.1063/5.0147219).
"""
function n_forest_sym_rule(n::Int)
    # state variables  
    @variables x[1:n] y[1:n]

    # penalty variables 
    @variables α[1:n]

    # time derivatives 
    dx = Matrix{Num}(undef, n, 2)

    # parameters 
    @variables ρ f h a₁ a₂ a b c 

    x_ix = 1
    y_ix = 2
    for i in 1:n 
        γ = a*(y[i] - b)^2 + c
        dx[i, x_ix] = ρ*y[i] - γ*x[i] - f*x[i] + a₁*α[i]*x[i]
        dx[i, y_ix] = f*x[i] - h*y[i] + a₂*α[i]*y[i]
    end 

    # reshape to an [x1 y1, ..., xi yi] dx and a [x1, y1, ..., xi, yi] x
    return (reshape(permutedims(dx), :), 
            reshape(permutedims(reshape([x; y], n, 2)), :))
end 

"""
    n_forest_sym_jacob(n::Int)

Return symbolic jacobian for `n`-forest complex network without perturbations.

Symbolic jacobian will be of shape `2*n × 2*n`. See also 
[`n_forest_sym_rule`](@ref).
"""
function n_forest_sym_jacob(n::Int)
    dx, x = n_forest_sym_rule(n)
    return Symbolics.jacobian(dx, x)
end 

"""
    one_forest_system(u0; ρ, f, h, a = 1, b = 1, c = 1)

Coupled ODE system for modeling the dynamics of a single forest ecosystem's
young and old tree density.
"""
function one_forest_system(u0; ρ, f, h, a = 1, b = 1, c = 1)
    params = [ρ, f, h, a, b, c]
    return CoupledODEs(antonovsky_rule, u0, params)
end 

"""
    antonovsky_rule(u, params, t)

Return the density of young trees (`x`) and old trees (`y`) for a one-species,
two sub-population forest ecosystem. 

[1] : Equation (1) from Cantin2020

[2] : Antonovsky1990
"""
function antonovsky_rule(u, params, t)
    x, y = u
    ρ, f, h, a, b, c = params
    xdot = ρ*y - γ(y, a, b, c)*x - f*x
    ydot = f*x - h*y
    du = SVector(xdot, ydot)
    #@show typeof(du) typeof(t)
    return du
end 

"""
    antonovsky_jacob(u, params, t)

Jacobian for the Antonovsky forest ecosystem model.
"""
function antonovsky_jacob(u, params, t) 
    x, y = u
    ρ, f, h, a, b, c = params
    return SMatrix{2,2}(-c - f - a*((y - b)^2), f, ρ - 2a*x*(y - b), -h)
end

"""
    antonovsky_sym_rule()

Return the symbolic representation of the Antonovsky dynamical system
and the state variables of that system.
"""
function antonovsky_sym_rule()
    # state variables
    @variables x y

    # time derivatives
    @variables xdot ydot

    # parameters
    @variables a b c ρ f h γ
    
    # antonovsky dynamical system
    γ = a*(y-b)^2 + c
    xdot = ρ*y - γ*x - f*x
    ydot = f*x - h*y

    return([xdot, ydot], [x, y])
end 

"""
    antonovsky_sym_jacob()
    
Return the symbolic Jacobian of the Antonovsky model.
"""
function antonovsky_sym_jacob() 
    dx, x = antonovsky_sym_rule()
    return Symbolics.jacobian(dx, x)
end 
