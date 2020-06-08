''' This module contains functions necessary to fit a negative binomial
using the maximum likelihood estimator and some numerical analysis
@author: Peter Xenopoulos
@website: http://www.peterxeno.com
@downloaded 17 oct 2019 from https://github.com/pnxenopoulos/negative_binomial/blob/master/negative_binomial/core.py
'''

import math
import numpy as np

from scipy.optimize import newton
from scipy.special import digamma

def r_derv(r_var, vec):
    ''' Function that represents the derivative of the neg bin likelihood wrt r
    @param r: The value of r in the derivative of the likelihood wrt r
    @param vec: The data vector used in the likelihood
    '''
    if not r_var or not vec:
        raise ValueError("r parameter and data must be specified")

    if r_var <= 0:
        raise ValueError("r must be strictly greater than 0")

    total_sum = 0
    obs_mean = np.mean(vec)  # Save the mean of the data
    n_pop = float(len(vec))  # Save the length of the vector, n_pop

    for obs in vec:
        total_sum += digamma(obs + r_var)

    total_sum -= n_pop*digamma(r_var)
    total_sum += n_pop*math.log(r_var / (r_var + obs_mean))

    return total_sum

def p_equa(r_var, vec):
    ''' Function that represents the equation for p in the neg bin likelihood wrt p
    @param r: The value of r in the derivative of the likelihood wrt p
    @param vec: Te data vector used in the likelihood
    '''
    if not r_var or not vec:
        raise ValueError("r parameter and data must be specified")

    if r_var <= 0:
        raise ValueError("r must be strictly greater than 0")

    data_sum = np.sum(vec)
    n_pop = float(len(vec))
    p_var = 1 - (data_sum / (n_pop * r_var + data_sum))
    return p_var

def neg_bin_fit(vec, init=0.0001):
    ''' Function to fit negative binomial to data
    @param vec: The data vector used to fit the negative binomial distribution
    @param init: Set init to a number close to 0, and you will always converge
    '''
    if not vec:
        raise ValueError("Data must be specified")

    est_r = newton(r_derv, init, args=(vec,))
    est_p = p_equa(est_r, vec)
    return est_r, est_p
    

# r = 40
# p = 0.5
# size = 100000
#
# import matplotlib.pyplot as plt
# random_nb_data = np.random.negative_binomial(r, p, size)
# print(random_nb_data)
# a=[]
# for i in random_nb_data:
#     a.append(i)
# _ = plt.hist(a, bins='auto')
# plt.title("Histogram with 'auto' bins")
# plt.show()
# print(neg_bin_fit(a))