"""Perturb the complex network with deforestation."""

from typing import List, Optional, Union, Set, NamedTuple, Dict 
from collections import namedtuple

import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import solve_ivp
from scipy.integrate._ivp.ivp import OdeResult

from src.n_forest import w_i


def perturbation_solve(T: int, y0: np.ndarray, n_ecosystems: int):
    """

    Args:
        y0: Array of shape (n, 2)  n-ecosystems and 2 state variables.

    The first update of states 
    """
    assert T > 0, "evoluation of system occurs for more than 0 timesteps"
    assert isinstance(y0, np.ndarray), "initial values are stored in ndarray"
    assert len(y0) == n_ecosystems, "initial values for each ecosystem"

    t_eval = range(T+1)
    solutions = np.empty() # TODO
    for initial_ecosystem_state_i in y0:
        pass 

    solution_i: OdeResult = solve_ivp(
        perturbation_rule, [t, T], y0, t_eval=t_eval)
    
    for t_ix, t in range(t_eval):
        current_solution = None
        for forest_i in range(1, n_ecosystems):
            pass 

    return 


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


def gamma(y, a = 1, b = 1, c = 1):
    return a*(y-b)**2 - c


def alpha(
    xs_to_i: List[float], 
    ys_to_i: List[float], 
    d_to_i: List[float],
    i: int, 
    w_0: float = 1, 
    alpha_0: float = -1,
    beta_2: float = 1,
    l: float = 600,
    P_0: float = 1):
    """

    Args:
        xs_to_i: Densities of young trees up to ecosystem `i`.
        ys_to_i: Densities of old trees up to ecosystem `i`.
        d_to_i: Distances up to ecosystem `i`.
        i: The integer `i` denoting the `i^th` ecosystem.
        w_0: Threshold for water needed for postive effect on ecosystem.
        alpha_0: Initial penalty coefficient.
        beta_2: Water evaporation coefficient of old trees.
        l: Positive normalization coefficient from size of forested area.
        P_0: Average water quantity evaporated over nearby maritime zone.
            For figure 8, this is 1.05, otherwise it is 1.00.

    Returns:
        TODO
    
    References:
        Equation (8) and (9) in Cantin2020
    """
    return alpha_0 *\
        (1 - w_i(i=i, x=xs_to_i, y=ys_to_i, d=d_to_i, beta_2=beta_2, l=l, P_0=P_0)/w_0)


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


def epsilon_k(
        k: int, 
        ecosytem_ids: Union[Set[int], List[NamedTuple]]):
    """Boolean integer for deforestation effect on forest dynamical system.

    Args:
        k: The integer of the `k^th` ecosystem. 
        ecosystem_ids: Integers corresponding to the i^th
            ecosystem that will undergo deforestation. This is 
            `{i_1, ..., i_N}` in Equation (17) Cantin2020. This could also
            just be a list of EcosystemDeforestTime named tuples.

    Returns:
        1 if k in the set, 0 otherwise. 

    References:
        Equation (17) from Cantin2020
    """
    if isinstance(ecosytem_ids, list) \
        and isinstance(ecosytem_ids[0], namedtuple):
        ecosytem_ids = {ecosystem.ecosystem_id for ecosystem in ecosytem_ids}
    return k in ecosytem_ids


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
