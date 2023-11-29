#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 09:22:26 2023

@author: kburdge
"""

import ultranest
import numpy as np
import matplotlib.pyplot as plt
from ultranest.plot import cornerplot
def create_data(x,a,b,sigma):
    return a*x+b+np.random.normal(0,sigma,N)

def linear_model(x,a,b):
    return a*x+b

N=100
x=np.random.uniform(0,10,N)
a=1.3
b=1.1
sigma=1
dy=np.ones(len(x))*sigma
y=create_data(x,a,b,sigma)
plt.errorbar(x,y,dy,ls=' ')

def gauss_prob(par,mu,std):

    return - 0.5 * (((par - mu) ** 2) / std ** 2)


param_names = ["a","b"]


def my_prior_transform(cube):
        params = cube.copy()

        # params[0] = cube[0]*10 - 5
        # params[1] = cube[1]*10 -5

        hi = -5
        lo = -10

        # uniform prior
        params[0] = cube[0]*(hi-lo) + lo
        params[1] = cube[1]*(hi-lo) + lo

        # # normal prior
        # import scipy
        # params[0] = scipy.stats.norm.ppf(cube[0],0,hi)
        # params[1] = scipy.stats.norm.ppf(cube[1],0,hi)
        
        return params
    
def my_likelihood(cube):
        a=cube[0]
        b=cube[1]
        
        predicted_values=linear_model(x,a,b)
        
        prob=np.sum(gauss_prob(predicted_values,y,dy))
        
        print(prob)
        return prob

# log directory: working_directory/results/
sampler = ultranest.ReactiveNestedSampler(param_names, my_likelihood, my_prior_transform,log_dir='results',resume=True)
result = sampler.run()
sampler.print_results()
sampler.plot_trace()
sampler.plot_run()
sampler.plot_corner()
