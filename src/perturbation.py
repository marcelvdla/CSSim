"""Perturb the complex network with deforestation."""

from typing import List, Optional, Union, Set, NamedTuple, Dict 
from collections import namedtuple

import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import solve_ivp
from scipy.integrate._ivp.ivp import OdeResult

from src.n_forest import w_i


def perturbation_rule(t, u, params: List):
    """General form of perturbed `i^th` dynamical ecosystem.

    NOTE: In principle, setting the epsilon parameter to none
    reproduce the general system without perturbation.

    TODO: If state \in R^{n x 2}, then this drastically simplifies
    the interdependency problem ( i think .... )

    Args:
        t: Timestep argument required by scipy.
        u: A state vector of shape [n * 2]. Element `i` is density of 
            young trees and element `i+1` is the density of old trees.
        params: Parameters for ecosystem dynamics. Parameters should be
            passed exactly in this order.
            - rho: Fertility of tree species.
            - f: Aging rate of young trees.
            - a_1: Biotic pump weight of young trees.
            - h: Mortality rate of old trees.
            - a_2: Biotic pump weight of old trees.
            - dists: Distance between `i` and `i+1` ecosystem.
            - w_0: Threshold for water needed for postive effect on ecosystem.
            - alpha_0: Initial penalty coefficient.
            - beta_2: Water evaporation coefficient of old trees.
            - l: Positive normalization coefficient from size of forested area.
            - P_0: Average water quantity evaporated over nearby maritime zone.
                For figure 8, this is 1.05, otherwise it is 1.00.
            - ecosystem_id_t_star: List of namedtuple's with `ecosystem_id` and
                `t_star` fields indicating deforestation time of ecosystem.
            - n_ecosystems: Number of ecosystems.

    Returns:
        A vector where each even index is `x_i` and each odd index is `y_i`
        for `i in range(n_ecosystems)`.

    References:
        Equation (18) in Cantin2020
    """
    # Extract parameters
    rho,  \
    f,  \
    a_1, \
    h, \
    a_2, \
    dists, \
    w_0, \
    alpha_0, \
    beta_2, \
    l, \
    P_0, \
    ecosystem_id_t_star, \
    n_ecosystems = params 
    
    # iterate through ecosystems and compute new values
    n_state_vars = 2
    du = np.empty(shape=(n_ecosystems * n_state_vars,))
    du_counter = 0
    for i in range(n_ecosystems):
        # TODO: is this getting the appropriate ecosystem vectors...
        raise NotImplementedError
        x_i, y_i = u[du_counter], u[du_counter+1]
        xs_to_i = u[0:(2*i)+1:n_state_vars]
        ys_to_i = u[1:(2*i)+2:n_state_vars]

        try: 
            # never prints...
            print(i, len(x_i), len(y_i))
        except:
            pass 

        alpha_i = alpha(
            xs_to_i=xs_to_i, 
            ys_to_i=ys_to_i,
            d_to_i=dists[:i],
            i=i,
            w_0=w_0,
            alpha_0=alpha_0,
            beta_2=beta_2,
            l=l,
            P_0=P_0)
        
        # Determine deforestation phenomena
        if ecosystem_id_t_star is None:
            epsilon_i = 0
            theta_i = 0
        else:
            epsilon_i = epsilon_k(
                k=i, ecosytem_ids=ecosystem_id_t_star.ecoystem_id)
            theta_i = theta(t=t, t_star=ecosystem_id_t_star.t_star)

        # Compute change in densities of forest ecosystem i
        xdot_i = rho*y_i - gamma(y_i)*x_i - f*x_i \
            + a_1*alpha_i*x_i \
            - epsilon_i*theta_i*x_i
        ydot_i = f*x_i - h*y_i + a_2*alpha_i*y_i - epsilon_i*theta_i*y_i

        # Update derivative of state vector
        du[du_counter] = xdot_i
        du[du_counter+1] = ydot_i
        du_counter += n_state_vars

    return du


def gamma(y, a = 1, b = 1, c = 1):
    """Mortality rate of young trees."""
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
    """Penalty function for ecosystems receiving optimal/suboptimal water.

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
        Penalty coefficient.
    
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
