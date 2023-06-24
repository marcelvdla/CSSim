using DynamicalSystems
using Symbolics
using UnPack
include("common.jl")

"""
    n_forest_system(u0::Matrix, params::Dict{Symbol, Any})
"""
function n_forest_system(u0::Matrix, params::Dict{Symbol, Any})
    return CoupledODEs(n_forest_rule, u0, params)
end 

"""
    n_forest_rule(u::Matrix, params::Dict{Symbol, Any}, t)

[1] : Equation (10) from Cantin2020
"""
function n_forest_rule(u::Matrix, params, t)
    @unpack ρ, f, a₁, h, a₂, d, l, α₀, P₀, β₁, β₂, n = params
    @assert n >= 1 "At least 1 forest ecosystem"
    @assert length(d) == (n-1), "n-1 distances in distance vector `d`"

    n_state_vars = 2
    x_ix = 1
    y_ix = 2

    x = u[:, x_ix]
    y = u[:, y_ix]

    du = zeros(n, n_state_vars)
    for i in 1:n 
        xᵢ = u[i, x_ix] 
        yᵢ = u[i, y_ix]
        wᵢ = n > 1 ? w(i, x, y, d, l, P₀, β₁, β₂) : 0
        αᵢ = n > 1 ? α(wᵢ, α₀, w₀) : 0
        du[i, x_ix] = ρ*yᵢ - γ(yᵢ)*xᵢ - f*xᵢ + a₁*αᵢ*xᵢ
        du[i, y_ix] = f*xᵢ - h*yᵢ + a₂*αᵢ*yᵢ
    end 

    return du
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
    xdot = ρ*y - gamma(y, a, b, c)*x - f*x
    ydot = f*x - h*y
    return SVector(xdot, ydot)
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
