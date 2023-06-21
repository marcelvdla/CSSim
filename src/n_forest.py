from typing import List, Set, Callable
from typing import NamedTuple

from collections import namedtuple 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def B(xi, yi, beta_1=0, beta_2=1):
        return beta_1 * xi + beta_2 * yi

def w_i(i, x, y, d, l=600, P_0=1.05):
    """ i::forest index
        x::array of densities of young trees
        y::array of densities of old trees
        d::distances between forests
    """
    if i == 0:
        return P_0  # the first forest receives only the base level of precipitation
    else:
        w = P_0 * np.exp(-np.sum(d[:i]) / l)
        for j in range(i, i+1):
            w += B(x[j], y[j]) * np.exp(-np.sum(d[j-1:i]) / l)
        return w

def alpha(x, y, dist, w_0=1, alpha_0=-1):
    """ Computes the penalty values and returns list of penalty per forest
    """
    d = len(x) * [dist]
    
    return [alpha_0 * (1 - w_i(i, x, y, d) / w_0) for i in range(len(x))]

def deriv_forest(x, y, penalty_rate, args):
    """
    Calculate derivatives of x and y over time as per the Antonovsky & Korzukhin rule.
    
    INPUT:
    x            = value for x; density of young trees in ecosystem, typically between 0 and 4
    y            = value for y; density of old trees in ecosystem, typically between 0 and 4
    penalty_rate = penalty rate calculated by penalty function
    args         = tuple of 6 arguments
    
    Returns derivatives of x and y.
    """
    # Unpack arguments rho, gamma (not used), f, h, a1 and a2
    fertility, mortality_young, aging_rate, biotic_pump_young, mortality_old, biotic_pump_old = args
    
    # Calculate derivatives
    dx = fertility * y - ((y - 1)**2 + 1) * x - aging_rate * x + biotic_pump_young * penalty_rate * x
    dy = aging_rate * x - mortality_old * y + biotic_pump_old * penalty_rate * y
    
    return dx, dy

def system_n_forests(x0s, y0s, args, timesteps = 100, dt = 0.01, dist=42):
    """
    Solves a system of ODEs
    INPUT:
    x0s       = array of values for x, the density of young trees in the ecosystem; typically between 0 and 4
    y0s       = array of values for y, the density of old trees in the ecosystem; typically between 0 and 4
    args      = tuple of arguments needed for the derivation function
    timesteps = int, timesteps to iterate over; default 100
    dt        = delta time, default 0.01
    
    Returns two arrays of x and y values per time
    """
       
    n = len(x0s)
    x_vals = np.empty((n, int(timesteps)))
    y_vals = np.empty((n, int(timesteps)))
    penalty = 0.5
    for t in range(timesteps):
        penalties = alpha(x0s, y0s, dist)
        for i in range(n):
            x_vals[i, t] = x0s[i]
            y_vals[i, t] = y0s[i]
            # implement solve ivp
            dx, dy = deriv_forest(x0s[i], y0s[i], penalties[i], args)
            x0s[i] += dx * dt
            y0s[i] += dy * dt
            
    return x_vals, y_vals


def theta(t: int, t_star: int):
    """Models beginning of deforestation process at time t*

    Args: 
        t: Timestep in integration of dynamical system.
        t_star: Timestep at which an ecosystem becomes deforested.

    Returns:
        The deforestation coefficient theta to be used in Equation (18)
        of Cantin2020.

    References:
        Equation (16) from Cantin2020
    """

    A = 2 # constant defined by Cantin2020
    deforestation_coefficient: float

    if t <= t_star:
        deforestation_coefficient = 0
    elif t_star < t and t < t_star + 1:
        deforestation_coefficient = A/2 - (A/2)*np.cos(np.pi*(t - t_star))
    else:
        deforestation_coefficient = A
    
    return deforestation_coefficient 


def epsilon_k(k: int, random_distinct_integers_i: Set[int]):
    """Boolean integer for deforestation effect on forest dynamical system.

    Args:
        k: The integer of the `k^th` ecosystem. 
        random_distinct_integers_i: Integers corresponding to the i^th
            ecosystem that will undergo deforestation. This is 
            `{i_1, ..., i_N}` in Equation (17) Cantin2020.

    Returns:
        1 if k in the set, 0 otherwise. 

    References:
        Equation (17) from Cantin2020
    """

    return k in random_distinct_integers_i 


def randomly_perturb_ecosystems(
        t_timesteps: int, 
        n_deforested_ecosystems: int, 
        n_total_ecosystems: int,
        seed: int = 42) -> List[NamedTuple]:
    """Randomly perturb ecosystems at random timesteps

    Example:

    ```python
    >>> perturbed_ecosystems = randomly_perturb_ecosystems(50, 2, 10)
    >>> perturbed_ecosystems
    [EcosystemDeforestTime(ecosystem_id=8, t_star=18),
     EcosystemDeforestTime(ecosystem_id=1, t_star=22)] 
    >>> perturbed_ecosystems[0].ecosystem_id
    8
    >>> perturbed_ecosystems[0].t_star
    18
    ```
    
    Args:
        t_timesteps: Number of timesteps over which integration will occur.
            For example, one might compute trajectories of a dynamical system
            over 50 timesteps (e.g., years), so then `t_timesteps = 50`.
        n_deforested_ecosystems: 
        n_total_ecosystems: Total number of ecosystems.    

    Returns:
        A list of named tuples. The named tuple has a `ecosystem_id` and `t_star`
        property. The `ecosystem_id` is deforested at time `t_star` according to 
        the `epsilon_k` and `theta` functions defined in this file.

    References:
        Equation (15) in Cantin2020
    """
    np.random.seed(seed)
    assert n_total_ecosystems >= 2, "there must be 2 or more ecosystems"
    assert n_deforested_ecosystems >= 0, \
        "number of ecosystems to perturb must be >= 0"
    
    # randomly generate ids for ecosystems that will be deforested 
    ecosystem_ids_to_deforest = np.random.choice(
        n_total_ecosystems, size=(n_deforested_ecosystems,), replace=False)

    # randomly generate time steps at which ecosystem ids will be deforested
    timesteps_to_deforest = np.random.choice(
        t_timesteps, size=(n_deforested_ecosystems))
    
    # combine the ecosystems and timesteps into a list of tuples
    ecosystem_time_tuples = []
    EcosystemDeforestTime = namedtuple(
        "EcosystemDeforestTime", ["ecosystem_id", "t_star"])

    for i in range(len(ecosystem_ids_to_deforest)):
        ecosystem_time_tuple = EcosystemDeforestTime(
            ecosystem_ids_to_deforest[i],
            timesteps_to_deforest[i])

        ecosystem_time_tuples.append(ecosystem_time_tuple)

    return ecosystem_time_tuples
