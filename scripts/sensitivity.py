#!/usr/bin/env python3
"""
This script computes total order Sobol indices for the N-forest network with 
biotic pump. This script adapts code that was provided during the ABM course
for sensitivity analysis to work with the biotic pump model. 
"""

# Imports
from SALib.sample import saltelli
from SALib.analyze import sobol

import sys
import pdb

import pandas as pd
import numpy as np
from itertools import combinations

from src.n_forest import *


def get_data_sobol(problem, replicates, distinct_samples):
    # Sample the parameter space
    param_values = saltelli.sample(problem, 
                                   distinct_samples,
                                   calc_second_order=False)
    print(param_values)

    # Create dataframe
    count = 0
    data = pd.DataFrame(index=range(replicates*len(param_values)), 
                                    columns=['alpha_0',
                                             'P_0',
                                             'w_0',
                                             'distance',
                                             'beta_1',
                                             'beta_2'])
    data['x1'], data['y1'], data['x2'], data['y2'] = None, None, None, None

    # Fixed variable values
    rho = 4.2     # fertility
    f = 1         # aging_rate
    a_1 = 0.1     # biotic_pump_young
    a_2 = 0.05    # biotic_pump_old
    h = 2         # mortality_old

    for i in range(replicates):
        for vals in param_values: 
            # Transform to dict with parameter names and their values
            variable_parameters = {}
            for name, val in zip(problem['names'], vals):
                variable_parameters[name] = val

            # Sample parameters
            dist = variable_parameters['distance']
            beta_1 = variable_parameters['beta_1']
            beta_2 = variable_parameters['beta_2']
            P_0 = variable_parameters['P_0']
            w_0 = variable_parameters['w_0']
            alpha_0 = variable_parameters['alpha_0']

            arguments = (rho,f,a_1,h,a_2,dist,beta_2,P_0,w_0,alpha_0,beta_1)
            
            # Run model
            x, y = system_n_forests(np.random.uniform(0,5,2),
                                    np.random.uniform(0,5,2),
                                    arguments, timesteps=800, 
                                    dt = 0.01)
            densities = np.array([x[0][-1],y[0][-1],x[1][-1],y[1][-1]])

            # Add data to the dataframe
            data.iloc[count, 0:6] = vals
            data.iloc[count, 6] = densities[0]
            data.iloc[count, 7] = densities[1]
            data.iloc[count, 8] = densities[2]
            data.iloc[count, 9] = densities[3]
            count += 1

            # print progress
            print(f'{count/(len(param_values)*(replicates))*100:.2f}% done',
                  end='\r')
    
    return data 


def plot_index(s, params, i, title=''):
    """
    Creates a plot for Sobol sensitivity analysis that shows the 
    contributions of each parameter to the global sensitivity.

    Args:
        s (dict): dictionary {'S#': dict, 'S#_conf': dict} of dicts that hold
            the values for a set of parameters
        params (list): the parameters taken from s
        i (str): string that indicates what order the sensitivity is.
        title (str): title for the plot
    """

    if i == '2':
        p = len(params)
        params = list(combinations(params, 2))
        indices = s['S' + i].reshape((p ** 2))
        indices = indices[~np.isnan(indices)]
        errors = s['S' + i + '_conf'].reshape((p ** 2))
        errors = errors[~np.isnan(errors)]
    else:
        indices = s['S' + i]
        errors = s['S' + i + '_conf']
        plt.figure()

    l = len(indices)

    plt.title(title)
    plt.ylim([-0.2, len(indices) - 1 + 0.2])
    plt.yticks(range(l), params)
    plt.errorbar(indices, range(l), xerr=errors, linestyle='None', marker='o')
    plt.axvline(0, c='k')


if __name__ == "__main__":
    replicates = int(sys.argv[1])
    distinct_samples = int(sys.argv[2])

    problem = {
        'num_vars': 6,
        'names': ['alpha_0', 'P_0', 'w_0', 'distance', 'beta_1', 'beta_2'],
        'bounds': [[-2.0, -1.0], [0.8, 1.2], [0.75, 1.25], [10, 900], [0.0, 0.4], [0.0, 0.4]]
    }

    data = get_data_sobol(problem, replicates, distinct_samples)
    data.to_csv(f'sobol_{replicates}_{distinct_samples}.csv')

    # pdb.set_trace()
    
    Si_density = sobol.analyze(problem, data['y2'].values, print_to_console=True, calc_second_order=False)
    
    # First order, shouldnt be needed
    # plot_index(Si_density, problem['names'], '1', 'First order sensitivity')
    # plt.show()

    # Total order
    plot_index(Si_density, problem['names'], 'T', 'Total order sensitivity')
    plt.show()