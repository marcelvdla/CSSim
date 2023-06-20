
"""
    gamma(y, a = 1, b = 1, c = 1)

Return mortality rate of young trees

[1] : Equation (2) from [Cantin2020](https://www.sciencedirect.com/science/article/pii/S1476945X20300386)
"""
gamma(y, a = 1, b = 1, c = 1) = a*(y - b)^2 + c
