from typing import List, Set, Callable, Union
from typing import NamedTuple

from collections import namedtuple

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from src.perturbation import (
    get_ecosystems_to_perturb_at_times_tstar,
    theta,
    epsilon_k,
    EcosystemDeforestTime
)

from src.common import *


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

def deriv_forest_perturbed(x, y, penalty_rate, epsilon, theta, args):
    """
    Calculate derivatives of x and y over time as per the Antonovsky & Korzukhin rule.

    INPUT:
    x            = value for x; density of young trees in ecosystem, typically between 0 and 4
    y            = value for y; density of old trees in ecosystem, typically between 0 and 4
    penalty_rate = penalty rate calculated by penalty function
    args         = tuple of 9 arguments

    Returns derivatives of x and y.
    """
    # Unpack arguments rho, gamma (not used), f, h, a1 and a2
    fertility, aging_rate, biotic_pump_young, mortality_old, biotic_pump_old, dist, beta_2, P_0 = args

    # Calculate derivatives
    dx = fertility * y - ((y - 1)**2 + 1) * x - aging_rate * x + biotic_pump_young * penalty_rate * x - epsilon * theta * x
    dy = aging_rate * x - mortality_old * y + biotic_pump_old * penalty_rate * y - epsilon * theta * y

    return dx, dy


def system_n_forests(x0s, y0s, args, perturbed = False, timesteps = 100, dt = 0.01, ecosystems_to_deforest = 2):
    """Solve the system of ODEs of the Antonovosky rule using forward Euler.

    Args:
        x0s: array of values for x, the density of young trees in the ecosystem;
            typically between 0 and 4.
        y0s: array of values for y, the density of old trees in the ecosystem;
            typically between 0 and 4.
        args: tuple of arguments needed for the derivation function
        timesteps: int, timesteps to iterate over; default 100
        dt: Delta time, default 0.01
        ecosystems_to_deforest: int, how many ecosystems should be deforested initially

    Returns:
        Two arrays of x and y values per time
    """

    assert (isinstance(args, (tuple, list, np.ndarray))),(
            "Your args variable should be a list-like.")
    assert (len(args)== 11), ("Make sure you have all arguments needed included: "
                             "fertility, aging_rate, biotic_pump_young,"
                             "mortality_old, biotic_pump_old, dist,"
                              "beta_2, P_0, w_0, alpha_0, beta_1.")
    assert len(x0s) == len(y0s), ("The input vector for the young and old trees"
                                   "should be the same length.")

    n = len(x0s)
    x_vals = np.empty((n, int(timesteps/dt)))
    y_vals = np.empty((n, int(timesteps/dt)))



    if perturbed:

        # Choose networks to deforest; gets a list of instances the deforested ecosystems
        deforest_ecosystem_at_t : List[EcosystemDeforestTime] = \
            get_ecosystems_to_perturb_at_times_tstar(
                t_timesteps=int(timesteps/dt),
                n_deforested_ecosystems=ecosystems_to_deforest,
                n_total_ecosystems=n)

        # Get a list with the ids of the ecosystems that will be deforested
        deforest_ecosystem_ids = [ecosystem.ecosystem_id for ecosystem in deforest_ecosystem_at_t]

    for t in range(int(timesteps / dt)):

        penalties = alpha(x0s, y0s, dist=args[5]/(n-1), beta_2 = args[6],
                          P_0 = args[7], w_0=args[8], alpha_0=args[9], beta_1= args[10])
        for i in range(n):

            x_vals[i, t] = x0s[i]
            y_vals[i, t] = y0s[i]

            if perturbed:

                # Checks if this ecosystem i will be deforested
                epsilon_i = epsilon_k(k=i, ecosytem_ids=deforest_ecosystem_ids)


                if epsilon_i:
                    # Get index of this ecosystem i in the list of deforested ecosystems to get the t_star
                    t_star_ix = deforest_ecosystem_ids.index(i)
                    used_t_star = deforest_ecosystem_at_t[t_star_ix].t_star

                    # Calculate theta using this t_star
                    theta_i = theta(t=t, t_star=used_t_star)

                else:
                    theta_i = 0

                # Use epsilon and theta in the other derivative function when perturbing the system
                dx, dy = deriv_forest_perturbed(x0s[i], y0s[i], penalties[i], epsilon_i, theta_i, args)

            else:
                dx, dy = deriv_forest(x0s[i], y0s[i], penalties[i], args)

            x0s[i] += dx * dt
            y0s[i] += dy * dt

    if perturbed:
        return x_vals, y_vals, deforest_ecosystem_ids

    else:
        return x_vals, y_vals
