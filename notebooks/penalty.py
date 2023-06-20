import numpy as np

def alpha(x, y, w_0=1, alpha_0=-1, P_0=1.05):
    """ Computes the penalty values and returns list of penalty per forest
    """
    d = len(x) * [42]       # distance hardcoded to 42km

    def B(xi, yi, beta_1=0, beta_2=1):
        return beta_1 * xi + beta_2 * yi

    def w_i(i, x, y, d):
        """ i::forest index
            x::array of densities of young trees
            y::array of densities of old trees
            d::distances between forests
        """
        
        if i == 0:
            return P_0  # the first forest receives only the base level of precipitation
        else:
            w = P_0 * np.exp(-np.sum(d[:i]) / sum(d))
            for j in range(i, i+1):
                w += B(x[j], y[j]) * np.exp(-np.sum(d[j-1:i]) / sum(d))
            return w
    
    return [alpha_0 * (1 - w_i(i, x, y, d) / w_0) for i in range(len(x))]