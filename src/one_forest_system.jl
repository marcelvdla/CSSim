using DynamicalSystems
using CSSim

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
