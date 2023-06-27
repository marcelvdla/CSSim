using DrWatson
@quickactivate :CSSim

using ChaosTools
using LinearAlgebra: I
using IntervalArithmetic: Interval


generic_system(u0, n_states) = CoupledODEs(generic_rule, u0, n_states)

two_state_system(u0) = CoupledODEs(two_state_rule, u0, 2)

function generic_rule(u, n_states, t)
    du = Vector{Union{Float64, Interval}}(undef, n_states) # PROBLEM HERE!!!!
    for i in 1:n_states
        du[i] = u[i]
    end 
    return SVector{n_states}(du)
end 

function generic_jacob(u, n_states, t)
    return SMatrix{n_states, n_states}(I)
end 

function two_state_rule(u, n_states, t)
    @show typeof(u) typeof(t)
    SVector(u[1], u[2])
end 

# calling the two state function (this will work)
u0 = rand(2)
ds = two_state_system(u0)
box = interval(0,1) × interval(0,1)

fps, eigs, stable = fixedpoints(ds, box, generic_jacob)

# calling the generic function with 3 states (this will fail)
n_states = 3
u0 = rand(n_states)
state_i_interval = interval(0, 1)
box = IntervalBox(repeat([state_i_interval], n_states)...)
ds = generic_system(u0, n_states)

fps, eigs, stable = fixedpoints(ds, box, generic_jacob)