"""
    w(i, x, y, d, l, P₀, β₁ = 0, β₂ = 1)

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

[1] : Equation (7) from Cantin2020
"""
B(xᵢ, yᵢ, β₁ = 0, β₂ = 1) = β₁*xᵢ + β₂*yᵢ

"""
    α(wᵢ, α₀, w₀)

Return penalty value α_i(w), using α_0, w_0

TODO: α_0, w_0 should be global for the system, could get rid of them somehow in the function call?

[1] : Equation (8) from Cantin2020
"""
α(wᵢ, α₀ = -1, w₀ = 1) = α₀*(1 - wᵢ / w₀)

"""
    γ(y, a = 1, b = 1, c = 1)

Return mortality rate of young trees

[1] : Equation (2) from Cantin2020
"""
γ(y, a = 1, b = 1, c = 1) = a*(y - b)^2 + c