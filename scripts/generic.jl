# TODO: Search for generic lorenz96 or something with n-states 
# and which also implements this OOP for reference
using DrWatson
@quickactivate :CSSim

using ChaosTools
using LinearAlgebra: I
using IntervalArithmetic: Interval


function generic_system(u0, n_states) 
    @show "generic system $(n_states)"
    CoupledODEs(generic_rule, u0, n_states)
end 

function two_state_system(u0) 
    @show "two state system"
    CoupledODEs(two_state_rule, u0, 2)
end 

function three_state_system(u0) 
    @show "three state system"
    CoupledODEs(three_state_rule, u0, 3)
end

function generic_rule(u, n_states, t)
    #du = Vector{Union{Float64, Interval}}(undef, n_states) # PROBLEM HERE!!!!
    #du = similar(u)
    du = undef
    if typeof(u) == SVector{Float64} || typeof(u) == Vector{Float64}
        @show "make float vector"
        du = Vector{Float64}(undef, n_states)
    else 
        @show "make interval vector"
        du = Vector{Interval{Float64}}(undef, n_states)
    end 
    
    for i in 1:n_states
        du[i] = u[i]
    end 
    @show typeof(u) typeof(du)
    return SVector{n_states}(du)
end 

function generic_jacob(u, n_states, t)
    return SMatrix{n_states, n_states}(I)
end 

function two_state_rule(u, n_states, t)
    @show typeof(u) typeof(t)
    SVector(u[1], u[2])
end 

function three_state_rule(u, n_states, t)
    @show typeof(u) typeof(t)
    SVector(u[1], u[2], u[3])
end 

"""
How to construct a variable length system here using a for loop ....
"""
function lorenz_rule(u, p, t)
    σ = p[1]; ρ = p[2]; β = p[3]
    du1 = σ*(u[2]-u[1])
    du2 = u[1]*(ρ-u[3]) - u[2]
    du3 = u[1]*u[2] - β*u[3]
    du =  SVector{3}(du1, du2, du3)
    @show typeof(u) typeof(du)
    return du
end

function lorenz_jacob(u, p, t)
    σ, ρ, β = p
    return SMatrix{3,3}(-σ, ρ - u[3], u[2], σ, -1, u[1], 0, -u[1], -β)
end

# calling the two state function (this will work)
u0 = rand(2)
ds = two_state_system(u0)
box = interval(0,1) × interval(0,1)

fps, eigs, stable = fixedpoints(ds, box, generic_jacob)

# calling the 3 state function 
u0 = rand(3)
ds = three_state_system(u0)
box = interval(0,1) × interval(0,1) × interval(0,1)

fps, eigs, stable = fixedpoints(ds, box, generic_jacob)

# calling the lorenz function
ρ, β = 30.0, 10/3
@show "lorenz system"
lorenz = CoupledODEs(lorenz_rule, 10ones(3), [10.0, ρ, β])
x = y = interval(-20, 20)
z = interval(0, 40)
box = x × y × z

fp, eigs, stable = fixedpoints(lorenz, box, lorenz_jacob)


# calling the generic function with 3 states (this will fail)
n_states = 3
u0 = rand(n_states)
state_i_interval = interval(0, 1)
box = IntervalBox(repeat([state_i_interval], n_states)...)
ds = generic_system(u0, n_states)

fps, eigs, stable = fixedpoints(ds, box, generic_jacob)