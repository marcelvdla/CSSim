#!/usr/bin/env python3
"""
Note: this file can only be run when the seed in the function 
get_ecosystems_to_perturb_at_times_tstar is commented out,
otherwise it will only do the same combination, when calling
this file also provide the distance (e.g. python3 killedheatmap.py 50)
"""

# Imports
from src.n_forest import *
from itertools import combinations
import sys

if __name__=="__main__":
    rho = 4.2     # fertility
    f = 1         # aging_rate
    a_1 = 1       # biotic_pump_young
    a_2 = 0       # biotic_pump_old
    h = 2         # mortality_old
    dist = int(sys.argv[1])

    beta_2 = 0.15
    P_0 = 1.05
    w_0 = 1
    alpha_0 = -1.0
    beta_1 = 0
    
    args = (rho,f,a_1,h,a_2,dist,beta_2,P_0,w_0,alpha_0,beta_1)
    killgrid = np.zeros((10,10))
    combs = list(combinations(range(10),2))

    # Find number of killed trees for all combinations
    # This is like rolling a dice, since the implementation
    # didn't support choosing which ecosystems to deforest
    while len(combs) > 0:
        # Random initial densities
        x0s, y0s = np.random.uniform(0,4,10), np.random.uniform(0,4,10) 
        x, y, ids = system_n_forests(x0s, y0s, args, perturbed = True)

        if tuple(ids) in combs:
            combs.remove(tuple(ids))
            print('removed', ids, 'remaining: ', len(combs))

        # Add number of killed forests to the grid
        num_kills = sum([1 for forest in y if np.isclose(forest[-1], 0.0)])
        killgrid[ids[0], ids[1]] = num_kills

    print(killgrid)

    fig, ax = plt.subplots()
    im = ax.imshow(killgrid, origin='lower', cmap='Reds')
    ax.set_xlabel('ID of deforest 1')
    ax.set_xticks(range(10))
    ax.set_xticklabels(range(1,11))
    ax.set_ylabel('ID of deforest 2')
    ax.set_yticks(range(10))
    ax.set_yticklabels(range(1,11))
    fig.colorbar(im)
    plt.show()