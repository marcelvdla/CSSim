using DynamicalSystems

"""
    one_forest_system(u0; p, f, h)

Coupled ODE system for modeling the dynamics of a single forest ecosystem's
young and old tree density.
"""
function one_forest_system(u0; p, f, h, a, b, c)
    params = [p, f, h, a, b, c]
    return CoupledODEs(antonovsky_rule, u0, params)
end 

"""
    antonovsky_rule(u, params, t)

Return the density of young trees (`x`) and old trees (`y`) for a one-species,
two sub-population forest ecosystem. 

[1] : Equation (1) from [Cantin2020](https://www.sciencedirect.com/science/article/pii/S1476945X20300386)

[2] : [Antonosky1990](https://www.sciencedirect.com/science/article/abs/pii/004058099090043U?via%3Dihub)
"""
function antonovsky_rule(u, params, t)
    x, y = u
    ρ, f, h, a, b, c = params
    xdot = ρ*y - gamma(y, a, b, c)*x - f*x
    ydot = f*x - h*y
    return SVector(xdot, ydot)
end 
