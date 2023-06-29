# Small functions used in the modeling of the n-forest system
using StatsBase
using Random

"""
    w(i, x, y, d, l, P₀, β₁ = 0, β₂ = 1)

Return water quantity received by ecosystem `i` from all other ecosystems.

[1] : Equation (3), (5), and (6) from Cantin2020
"""
function w(i, x, y, d, l, P₀, β₁, β₂)
    n = length(x)
    if i == 1
        return P₀
    else
        w_sum = (P₀ + B(x[1], y[1], β₁, β₂))*exp(-sum(d[1:n-1]) / l)
        for i in 2:n-1
            w_sum += B(x[i], y[i], β₁, β₂)*exp(-sum(d[i:n-1]) / l)
        end
        return w_sum
    end 
end 

"""
    B(xᵢ, yᵢ, β₁ = 0, β₂ = 1)

Return quantity of water evaporated of `i^th` forest assuming constant surface.

[1] : Equation (7) from Cantin2020
"""
B(xᵢ, yᵢ, β₁, β₂) = β₁*xᵢ + β₂*yᵢ

"""
    B(xᵢ, yᵢ, Sᵢ)

Return quantity of water evaporated of `i^th` forest assuming variable surface.

[1] : Equation (7) from Cantin2020
"""
B(xᵢ, yᵢ, Sᵢ) = β₁(Sᵢ)*xᵢ + β₂(Sᵢ)*yᵢ

β₁(Sᵢ) = throw("unimplemented")

β₂(Sᵢ) = throw("unimplemented")

"""
    α(wᵢ, α₀, w₀)

Return penalty value for ecosystems receiving optimal/suboptimal water.

`α` can be typed by `\\alpha<tab>`

[1] : Equation (8) from Cantin2020
"""
α(wᵢ, α₀, w₀) = α₀*(1 - wᵢ / w₀)

"""
    γ(y, a = 1, b = 1, c = 1)

Return mortality rate of young trees.

`γ` can be typed by `\\gamma<tab>`

[1] : Equation (2) from Cantin2020
"""
γ(y, a = 1, b = 1, c = 1) = a*(y - b)^2 + c

""" 
    θ(t, tstar, A = 2)

Return deforestation coefficient.

`θ` can be typed by `\\theta<tab>`

[1] : Equation (16) from Cantin2020
"""
function θ(t, tstar, A = 2)
    if t <= tstar
        return 0
    elseif tstar < t < tstar + 1
        return (A/2) - (A/2)*cos(π*(t - tstar))
    else
        return A 
    end
end 

"""
    εₖ(k, ecosystems_ids_to_deforest::AbstractVector)

Return deforestation bool from ecosystem id `k` in `ecosystems_ids_to_deforest`.

`εₖ` can be typed by `\\varepsilon<tab>\\_k<tab>`

[1] : Equation (17) from Cantin2020
"""
function εₖ(k, ecosystems_ids_to_deforest::AbstractVector)
    return k in collect(ecosystems_ids_to_deforest)
end 

"""
    ecosystems_times_to_deforest(
        T::Int, N::Int, n::Int, seed = 42)::Vector{EcosystemDeforestTime}

Return vector of elements `(i, tstar)` to randomly deforest ecosystem `i` at 
time `tstar`.

# Arguments 
- `T::Int`: Timesteps from which tstar will be sampled.
- `N::Int`: Number of ecosystems to deforest. 
- `n::Int`: Total number of ecosystems in the complex network.
- `seed=42`: Random seed for sampling.

[1] : Expression (15) from Cantin2020
"""
function ecosystems_times_to_deforest(
    T::Int, N::Int, n::Int, seed = 42)::Vector{EcosystemDeforestTime}
    @assert N <= n "Number of ecosystems to deforest is leq total ecosystems"
    rng = MersenneTwister(seed)
    ecosystem_ids_k = sample(rng, 1:n, N, replace=false)
    tstars = sample(rng, 1:T, N, replace=true)
    tups = [
            EcosystemDeforestTime(
                ecosystem_ids_k[i], 
                tstars[i]) 
            for i in 1:N
        ]

    return tups
end

"""
    ecosystems_times_to_deforest(
        T::Int, 
        deforest_ecosystems::Union{Vector{Int}, Vector{Any}}, 
        deforest_times::Union{Vector{Int}, Vector{Any}}, 
        n::Int, 
        seed = 42)::Vector{EcosystemDeforestTime}

Return vector of `(i, tstar)` with selected ecosystems/times to deforest.

# Arguments 
- `T::Int`: Timesteps from which tstar will be sampled.
- `deforest_ecosystems::Union{Vector{Int}, Vector{Any}}`: Vector of ecosystem 
    ids to deforest. Pass an empty vector `[]` if you want random ecosystem ids 
    to be used.
- `deforest_times::Union{Vector{Int}, Vector{Any}}`: Vector of times to 
        deforest. Pass an empty vector `[]` if you want random times to be used.
- `n`: Total number of ecosystems in the complex network.
- `seed=42`: Random seed for sampling.
"""
function ecosystems_times_to_deforest(
    T::Int, 
    deforest_ecosystems::Union{Vector{Int}, Vector{Any}}, 
    deforest_times::Union{Vector{Int}, Vector{Any}}, 
    n::Int, 
    seed = 42)
    # Compute length of args
    len_deforest_ecosystems = length(deforest_ecosystems)
    len_deforest_times = length(deforest_times)
    
    # Validate args
    at_least_one_vector = len_deforest_ecosystems > 0 || len_deforest_times > 0
    @assert at_least_one_vector "Length of one of `deforest_ecosystems`" 
        " or `deforest_times` is at least 1"

    @assert len_deforest_ecosystems <= n && len_deforest_times <= n "Number of" 
        " ecosystems to deforest is leq total ecosystems:"
        " $(len_deforest_ecosystems), $(len_deforest_times)"

    # Generate random samples as needed
    N = max(len_deforest_ecosystems, len_deforest_times)
    rng = MersenneTwister(seed)
    ecosystem_ids_k = len_deforest_ecosystems > 0 ? deforest_ecosystems :
        sample(rng, 1:n, N, replace=false)
    tstars = len_deforest_times > 0 ? deforest_times : 
        sample(rng, 1:T, N, replace=true)

    @assert length(ecosystem_ids_k) == length(tstars) "Valid deforestation pairs"

    # get deforestation vector
    tups = [
        EcosystemDeforestTime(
            ecosystem_ids_k[i], 
            tstars[i]) 
        for i in 1:N
    ]

    return tups
end 

function iskilled(
    y_last, 
    kiled_threshold = 0.0, 
    atol::AbstractFloat = 1e-3)
    return isapprox(y_last, kiled_threshold; atol = atol)
end 

struct EcosystemDeforestTime 
    "ecosystem id"
    i::Int
    "deforestation time during integration"
    tstar::Int 
end 