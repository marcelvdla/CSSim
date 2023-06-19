using DynamicalSystems

"""
    one_forest_system(u0; p, f, h)

Coupled ODE system for modeling the dynamics of a single forest ecosystem's
young and old tree density.
"""
function one_forest_system(u0; p, f, h)
    params = [p, f, h]
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
    ρ, f, h = params
    xdot = ρ*y - gamma(y)*x - f*x
    ydot = f*x - h*y
    return SVector(xdot, ydot)    
end 

"""
    gamma(y, a, b, c)

Return mortality rate of young trees

TODO: Define default params

[1] : Equation (2) from [Cantin2020](https://www.sciencedirect.com/science/article/pii/S1476945X20300386)
"""
gamma(y, a, b, c) = a*(y - b)^2 + c
