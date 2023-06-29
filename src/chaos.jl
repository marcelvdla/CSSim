# Functions needed for chaotic analysis of n-forest systems
using DynamicalSystems
using ChaosTools
using UnPack
using IntervalArithmetic: Interval, IntervalBox
using LinearAlgebra: I
using CSSim

"""
    n_forest_tds(u0::AbstractVector, params::Dict{Symbol, Any})

Tangent dynamical system for `n`-forest system without perturbation.

# Examples 
```jldoctest 
julia> n = 2
julia> params = Dict{Symbol,Any}(
    :ρ  => 4.2, :f  => 1.0, :α₀ => -1.0, :w₀ =>  1.0, :a₁ => 1.0, :h  => 2.0, 
    :a₂ => 0.0, :d  => [42],  :l  => 600.0,  :P₀ => 1.0, :β₁ => 0.0, :β₂ => 1.0,
    :n  => n, :ecosystems_to_deforest => [])
julia> u0 = 4*rand(n, 2)
julia> tds = n_forest_tds(u0, params)
```
"""
function n_forest_tds(u0::AbstractVector, params::Dict{Symbol, Any})
    ds = n_forest_system(u0, params)
    tds = TangentDynamicalSystem(ds; J = n_forest_jacob!)
    return tds 
end

"""
    n_forest_rule!(du, u::SVector, params::Dict{Symbol, Any}, t)

IIP `n`-forest rule for use with [`n_forest_tds`](@ref) since this requires
a 1D `u` IIP rule and 1D `u` IIP Jacobian.
"""
function n_forest_rule!(du, u, params::Dict{Symbol, Any}, t) 
    # Distance between `i` and `i+1` forest vector
    d::Union{Vector{Any}, Vector{<:Real}} = [] 

    # not used
    ecosystems_to_deforest::Union{
        Vector{Any}, Vector{EcosystemDeforestTime}} = []

    @unpack ρ, f, a₁, h, a₂, d, l, 
            α₀, w₀, P₀, β₁, β₂, 
            ecosystems_to_deforest,
            n  = params
    @assert n >= 1 "At least 1 forest ecosystem"
    @assert length(d) == (n-1) "n-1 distances in distance vector `d`"

    n_states = 2
    n = length(u) ÷ n_states

    # Get all x and y state variables
    x = u[1:n]
    y = u[n+1:length(u)]

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
        du[i] = ẋᵢ
        du[i+n] = ẏᵢ
    end 

    return nothing
end 

""" 
    n_forest_jacob!(J, u::AbstractVector, params::Dict{Symbol, Any}, t)

Jacobian for `n`-forest system without perturbation and `u` as a vector of 
length `n * 2` such that `[x1, x2, ..., xn, y1, y2, ..., yn]`. The rule for 
the Jacobian was determined by visual inspection of the outputs of 
[`n_forest_sym_jacob`](@ref). Note that this will be used with 
[`ChaosTools.fixedpoints`](@ref).

# Examples 
```jldoctest
julia> # Test jacobian for one forest system
julia> n = 1;
julia> params = Dict(
    :ρ => 1, :f => 1, :a₁ => 1, :h => 1, :a₂ => 1, 
    :d => [], :l => 1, :α₀ => 1, :w₀ => 1, :P₀ => 1, 
    :β₁ => 1, :β₂ => 1, :ecosystems_to_deforest => [], 
    :n => n);
julia> J = zeros(n*2, n*2);
julia> u = SVector{n*2}(ones(n, 2));
julia> n_forest_jacob!(J, u, params, 1);
julia> J
2×2 Matrix{Float64}:
 -2.0   1.0
  1.0  -1.0
julia> # Test jacobian for two forest system
julia> n = 2
julia> params[:n] = n
julia> params[:d] = [1]
julia> J = zeros(n*2, n*2);
julia> u = SVector{n*2}(ones(n, 2));
julia> n_forest_jacob!(J, u, params, 1);
julia> J
4×4 Matrix{Float64}:
 -2.0   1.0   0.0       0.0
  1.0  -1.0   0.0       0.0
  0.0   0.0  -2.10364   1.0
  0.0   0.0   1.0      -1.10364
```
"""
function n_forest_jacob!(J, u::AbstractVector, params::Dict{Symbol, Any}, t)
    # Distance between `i` and `i+1` forest vector
    d::Union{Vector{Any}, Vector{<:Real}} = [] 

    # Not used
    ecosystems_to_deforest::Union{
        Vector{Any}, Vector{EcosystemDeforestTime}} = []
    
    @unpack ρ, f, a₁, h, a₂, d, l, 
            α₀, w₀, P₀, β₁, β₂, 
            ecosystems_to_deforest,
            n  = params
    @assert n >= 1 "At least 1 forest ecosystem"
    @assert length(d) == (n-1) "n-1 distances in distance vector `d`"
        
    a = b = c = 1
    J[:, :] .= 0
 
    n_states = 2

    # Get all x and y state variables
    x = u[1:n]
    y = u[n+1:n*n_states]
    
    # define the jacobian inplace using the rule determined 
    ecosystem_id = 1
    for i in 1:2:2*n
        xᵢ = x[ecosystem_id]
        yᵢ = y[ecosystem_id]
        wᵢ = w(ecosystem_id, x, y, d, l, P₀, β₁, β₂)
        αᵢ = α(wᵢ, α₀, w₀)
        
        J[i, i] = -c - f + a₁*αᵢ - a*(-b + yᵢ)^2
        J[i, i+1] = ρ - 2*a*(-b + yᵢ)*xᵢ
        J[i+1, i] = f
        J[i+1, i+1] = -h + a₂*αᵢ
        ecosystem_id += 1
    end 
    
    return nothing
end


"""
    n_forest_system_oop(u0, params::Dict{Symbol, Any})

OOP `n`-forest CoupledODEs system without perturbation.
"""
function n_forest_system_oop(u0, params::Dict{Symbol, Any})
    return CoupledODEs(n_forest_rule, u0, params)
end 

"""
    n_forest_rule(u, params::Dict{Symbol, Any}, t)

OOP `n`-forest rule without perturbation with `u` as a vector of length 
`n * 2` such that `[x1, x2, ..., xn, y1, y2, ..., yn]`. This rule is an
out-of-place creation of the time derivative vector `du` since this was 
originally required for use with [`ChaosTools.fixedpoints`](@ref)
"""
function n_forest_rule(u::SVector, params::Dict{Symbol, Any}, t) 
    # Distance between `i` and `i+1` forest vector
    d::Union{Vector{Any}, Vector{<:Real}} = [] 

    # not used
    ecosystems_to_deforest::Union{
        Vector{Any}, Vector{EcosystemDeforestTime}} = []

    @unpack ρ, f, a₁, h, a₂, d, l, 
            α₀, w₀, P₀, β₁, β₂, 
            ecosystems_to_deforest,
            n  = params
    @assert n >= 1 "At least 1 forest ecosystem"
    @assert length(d) == (n-1) "n-1 distances in distance vector `d`"

    n_states = 2
    n = length(u) ÷ n_states

    @show typeof(u)

    # Get all x and y state variables
    x = u[1:n]
    y = u[n+1:length(u)]

    # For storing the computed time derivatives
    du = similar(u) # problem here

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
    two_forest_system
"""
function two_forest_system(u0, params::Dict{Symbol, Any})
    return CoupledODEs(two_forest_rule, u0, params)
end 

function two_forest_rule(u, params::Dict{Symbol, Any}, t)
    @unpack ρ, f, a₁, h, a₂, d, l, 
            α₀, w₀, P₀, β₁, β₂, 
            ecosystems_to_deforest, # not used
            n  = params
    @assert n >= 1 "At least 1 forest ecosystem"
    @assert length(d) == (n-1) "n-1 distances in distance vector `d`"

    # state variables 
    x₁ = u[1]
    x₂ = u[2]
    x = SVector(x₁, x₂)

    y₁ = u[3]
    y₂ = u[4]
    y = SVector(y₁, y₂)

    # Biotic pump
    w₁ = n > 1 ? w(1, x, y, d, l, P₀, β₁, β₂) : 0
    α₁ = n > 1 ? α(w₁, α₀, w₀) : 0

    w₂ = n > 1 ? w(2, x, y, d, l, P₀, β₁, β₂) : 0
    α₂ = n > 1 ? α(w₂, α₀, w₀) : 0

    # forest systems 
    ẋ₁ = ρ*y₁ - γ(y₁)*x₁ - f*x₁ + a₁*α₁*x₁
    ẏ₁ = f*x₁ - h*y₁ + a₂*α₁*y₁

    ẋ₂ = ρ*y₂ - γ(y₂)*x₂ - f*x₂ + a₁*α₂*x₂
    ẏ₂ = f*x₂ - h*y₂ + a₂*α₂*y₂

    return SVector(ẋ₁, ẋ₂, ẏ₁, ẏ₂)
end 

function two_forest_jacob(u, params::Dict{Symbol, Any}, t)
    @unpack ρ, f, a₁, h, a₂, d, l, 
            α₀, w₀, P₀, β₁, β₂, 
            ecosystems_to_deforest, # not used
            n  = params
    a = b = c = 1
    # state variables 
    x₁ = u[1]
    x₂ = u[2]
    x = SVector(x₁, x₂)

    y₁ = u[3]
    y₂ = u[4]
    y = SVector(y₁, y₂)

    # Biotic pump
    w₁ = n > 1 ? w(1, x, y, d, l, P₀, β₁, β₂) : 0
    α₁ = n > 1 ? α(w₁, α₀, w₀) : 0

    w₂ = n > 1 ? w(2, x, y, d, l, P₀, β₁, β₂) : 0
    α₂ = n > 1 ? α(w₂, α₀, w₀) : 0
    
    J = SMatrix{4, 4}(
        [
            (a₁*α₁ - c - f - a*((y[1] - b)^2)) (ρ - 2a*(y[1] - b)*x[1]) 0 0
            (f) (a₂*α₁ - h) 0 0
            0 0 (a₁*α₂ - c - f - a*((y[2] - b)^2)) (ρ - 2a*(y[2] - b)*x[2])
            0 0 (f) (a₂*α₂ - h)
        ]
    )
    return J
end 

""" 
    boxes(n::Int, interval_ui)

Return `[interval]^2n` interval for an `n`
"""
function boxes(n::Int, interval_ui::Interval)::IntervalBox
    IntervalBox(repeat([interval_ui, interval_ui], n)...)
end 