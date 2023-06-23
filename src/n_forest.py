from typing import List, Set, Callable, Union
from typing import NamedTuple

from collections import namedtuple

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def B(xi: float, yi: float, beta_1=0, beta_2=1):
    """Return quantity of water evaporated over `i^th` forest.

    Default values from Table 3 in Cantin2020

    Args:
        xi: Density of young tree species in `i^th` forest.
        yi: Density of old tree spcies in `i^th` forest.
        beta_1: Water evaporation coefficient of young trees (`mm * year**-1`).
        beta_2: Water evaporation coefficient of old trees (`mm * year**-1``).

    References:
        Equation (7) in Cantin2020.
    """
    return beta_1 * xi + beta_2 * yi


def w_i(
    i: int,
    x: List[float],
    y: List[float],
    d: List[float],
    beta_2: float=1,
    l: float=600,
    P_0: float=1.0,):
    """Water received by ecosystems in network.

    This is the biotic pump mechanism. Default values from Table 3 in
    Cantin2020.

    Args:
        i: `i^th` forest.
        x: List of densities of young trees for N ecosystems.
        y: List of densities of old trees for N ecosystems.
        d: List of distances where `d_i` is distance between forest `i` and `i+1`.
        beta_2: Water evaporation coefficient of old trees.
        l: Positive normalization coefficient from size of forested area.
        P_0: Average water quantity evaporated over nearby maritime zone.
            For figure 8, this is 1.05, otherwise it is 1.00.

    References:
        Equations (3) and (6) in Cantin2020
    """
    if i == 0:
        # the first forest receives only the base level of precipitation
        # equation (3)
        return P_0
    elif i == 1:
        # equation (5)
        w = (P_0 + B(x[0], y[0], beta_2=beta_2)) * np.exp(-d[0] / l)
        return w
    else:
        # equation (6)
        w = (P_0 + B(x[0], y[0], beta_2=beta_2)) * np.exp(-np.sum(d[0:i]) / l)
        for j in range(1, i+1):
            w += B(x[j], y[j], beta_2=beta_2) * np.exp(-np.sum(d[j:i]) / l)
        return w


def alpha(
    x: List[float],
    y: List[float],
    dist: Union[float, List[float]],
    w_0=1,
    alpha_0=-1,
    beta_2=1,
    P_0=1.00):
    """Penalty function for quantity of water received by forest ecosystems.

    Computes the penalty values and returns list of penalty per forest.

    TODO: This function should likely be simplified because this might introduce
    complications to the dynamical system (solve_ivp).

    Args:
        x: List of young tree species densities for each ecosystem.
        y: List of old tree species densities for each ecosystem
        dist: The distance the `i^th` ecosystem is from the `i+1` ecosystem.
        w_0: Threshold above which densities of trees are positively affected.
        alpha_0: Negative penalty coefficient.

    References:
        Equation (8) and (9) in Cantin2020
    """
    assert alpha_0 < 0, "alpha_0 must be negative."

    # If not a arraylike, then forests all same distance from each other
    if not isinstance(dist, (np.ndarray, list)):
        d = len(x) * [dist]
    else:
        assert len(dist) == len(x), \
            "number of distances should match number of forests"
        d = dist

    return [alpha_0 * (1 - w_i(i, x, y, d, beta_2=beta_2, P_0=P_0) / w_0)
            for i in range(len(x))]


def deriv_forest(x, y, penalty_rate, args):
    """
    Calculate derivatives of x and y over time as per the Antonovsky &
    Korzukhin rule.

    INPUT:
    x            = value for x; density of young trees in ecosystem, typically
                   between 0 and 4
    y            = value for y; density of old trees in ecosystem, typically
                   between 0 and 4
    penalty_rate = penalty rate calculated by penalty function
    args         = tuple of 9 arguments

    Returns derivatives of x and y.
    """
    # Unpack arguments rho, gamma (not used), f, h, a1 and a2
    fertility, mortality_young, aging_rate, biotic_pump_young, \
    mortality_old, biotic_pump_old, dist, beta_2, P_0 = args

    # Calculate derivatives
    dx = fertility * y - ((y - 1)**2 + 1) * x - aging_rate * x \
        + biotic_pump_young * penalty_rate * x
    dy = aging_rate * x - mortality_old * y + biotic_pump_old * penalty_rate * y

    return dx, dy


def system_n_forests(x0s, y0s, args, timesteps = 100, dt = 0.01, dist=42, beta_2=1):
    """Solve the system of ODEs of the Antonovosky rule using forward euler.

    Args:
        x0s: array of values for x, the density of young trees in the ecosystem;
            typically between 0 and 4.
        y0s: array of values for y, the density of old trees in the ecosystem;
            typically between 0 and 4.
        args: tuple of arguments needed for the derivation function
        timesteps: int, timesteps to iterate over; default 100
        dt: Delta time, default 0.01

    Returns:
        Two arrays of x and y values per time
    """

    assert (isinstance(args, (tuple, list, np.ndarray))),(
            "Your args variable should be a list-like.")
    assert (len(args)== 9), ("Make sure you have all arguments needed included:"
                             " fertility, mortality_young, aging_rate,"
                             "biotic_pump_young,"
                             "mortality_old, biotic_pump_old, dist, beta_2, P_0.")
    assert len(x0s) == len(y0s), ("The input vector for the young and old trees"
                                   "should be the same length.")

    n = len(x0s)
    x_vals = np.empty((n, int(timesteps/dt)))
    y_vals = np.empty((n, int(timesteps/dt)))

    for t in range(int(timesteps / dt)):
        penalties = alpha(x0s, y0s, dist=args[6]/(n-1), beta_2 = args[7],
                          P_0 = args[8])
        for i in range(n):

            x_vals[i, t] = x0s[i]
            y_vals[i, t] = y0s[i]

            dx, dy = deriv_forest(x0s[i], y0s[i], penalties[i], args)
            x0s[i] += dx * dt
            y0s[i] += dy * dt


    return x_vals, y_vals
