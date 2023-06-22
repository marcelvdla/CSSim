"""Perturb the complex network with deforestation."""

from typing import List, Optional, Union, Set, NamedTuple, Dict 
from collections import namedtuple

import numpy as np
from scipy.integrate import solve_ivp

from src.n_forest import w_i


def perturbation_solve():
    """"""
    raise NotImplementedError


def perturbation_rule(t, state, params: List):
    """General form of perturbed `i^th` dynamical ecosystem.

    NOTE: In principle, setting the epsilon parameter always to 0 would
    reproduce the general system without perturbation.

    Args:
        t: Timestep argument required by scipy.
        state: A state vector whose first element is density of young
            trees and second element is density of old trees.
        params: Parameters for ecosystem dynamics. Parameters should be passed
            in the same order of appearance in equation (18) read from left
            to right and then top to bottom.
            - rho: Fertility of tree species.
            - f: Aging rate of young trees.
            - a_1: Biotic pump weight of young trees.
            - alpha_i: Penalty weight for high/low water received by forest `i`.
            - epsilon_i: Perturbation of system bool at timestep t.
            - theta_i: Models beginning of deforestation of forest `i` at `t`.
            - h: Mortality rate of old trees.
            - a_2: Biotic pump weight of old trees.

    Returns:
        Density of young and old species of trees, respectively.

    References:
        Equation (18) in Cantin2020
    """
    rho, f, a_1, alpha_i, epsilson_i, theta_i, h, a_2 = params 
    x_i, y_i = state 
    xdot_i = rho*y_i - gamma(y_i)*x_i - f*x_i + a_1*alpha_i*x_i \
        - epsilson_i*theta_i*x_i
    ydot_i = f*x_i - h*y_i + a_2*alpha_i*y_i - epsilson_i*theta_i*y_i
    return xdot_i, ydot_i


def gamma(y, a=1, b=1, c=1):
    """"""
    return a*(y-b)**2 - c


def alpha(xs_to_i, ys_to_i, i):
    """"""
    raise NotImplementedError


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
