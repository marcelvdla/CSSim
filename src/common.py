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

