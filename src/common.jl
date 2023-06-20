
"""
    
α(w_i, α_0, w_0)

Return penalty value α_i(w), using α_0, w_0

TODO: α_0, w_0 should be global for the system, could get rid of them somehow in the function call?

[1] : Equation (8) from [Cantin2020]
"""
α(w_i, α_0, w_0) = α_0 * (1 - w_i / w_0)