
"""
    gamma(y, a = 1, b = 1, c = 1)

Return mortality rate of young trees

[1] : Equation (2) from Cantin2020
"""
gamma(y, a = 1, b = 1, c = 1) = a*(y - b)^2 + c
