using StatsBase
using Random

"""
    w(i, x, y, d, l, P₀, β₁ = 0, β₂ = 1)

Return water quantity received by ecosystem `i` from all other ecosystems.

[1] : Equation (3), (5), and (6) from Cantin2020
"""
function w(i, x, y, d, l, P₀, β₁ = 0, β₂ = 1)
    n = length(x)
    if i == 1
        return P₀
    else
        w_sum = (P₀ + B(x[1], y[1]))*exp(-sum(d[1:n-1]) / l)
        for i in 2:n-1
            w_sum += B(x[i], y[i], β₁, β₂)*exp(-sum(d[i:n-1]) / l)
        end
        return w_sum
    end 
end 

"""
    B(xᵢ, yᵢ, β₁ = 0, β₂ = 1)

Return quantity of water evaporated of `i^th` forest.

[1] : Equation (7) from Cantin2020
"""
B(xᵢ, yᵢ, β₁ = 0, β₂ = 1) = β₁*xᵢ + β₂*yᵢ

"""
    α(wᵢ, α₀, w₀)

Return penalty value for ecosystems receiving optimal/suboptimal water.

[1] : Equation (8) from Cantin2020
"""
α(wᵢ, α₀ = -1, w₀ = 1) = α₀*(1 - wᵢ / w₀)

"""
    γ(y, a = 1, b = 1, c = 1)

Return mortality rate of young trees.

[1] : Equation (2) from Cantin2020
"""
γ(y, a = 1, b = 1, c = 1) = a*(y - b)^2 + c

""" 
    θ(t, tstar)

Return deforestation coefficient.

[1] : Equation (16) from Cantin2020
"""
function θ(t, tstar)
    A = 2
    if t <= tstar
        return 0
    elseif tstar < t < tstar + 1
        return (A/2) - (A/2)cos(π*(t-tstar))
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
    random_ecosystems_times_to_deforest(
        T, N, n, seed = 42)::Vector{EcosystemDeforestTime}

Return vector of named tuples `(i, tstar)` to deforest ecosystem `i` at time `tstar`.

# Arguments 
- `T`: Total time system system is evolved (e.g., T = 50 for 50 years)
- `N`: Number of ecosystems to deforest. 
- `n`: Total number of ecosystems in the complex network.
- `seed=42`: Random seed for sampling.

[1] : Expression (15) from Cantin2020
"""
function random_ecosystems_times_to_deforest(
    T, N, n, seed = 42)::Vector{EcosystemDeforestTime}
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

struct EcosystemDeforestTime 
    "ecosystem id"
    i::Int
    "deforestation time during integration"
    tstar::Int 
end 