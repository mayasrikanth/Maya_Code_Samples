#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 22:48:06 2018

@author: mayasrikanth
"""

import numpy as np 
import math

def computeGrad(weights):
    u = weights[0]
    v = weights[1]
    
    partial_u = 2 * (math.exp(v) + 2*v*math.exp(-1 * u))
    partial_u *= (u*math.exp(v) - 2*v*math.exp(-1 * u))
    partial_v = 2* (u*math.exp(v) - 2*v*math.exp(-1*u))*(u*math.exp(v)- 2*math.exp(-1*u))
    return np.array([partial_u, partial_v])

def endDescent(weights):
    u = weights[0]
    v = weights[1]
    error = (u * math.exp(v) - 2 * v * math.exp(- 1 * u))**2
    print("Error: " + str(error))
    if(error <= (10) **-14):
        return True
    return False 



def gradDescent():
    # initialize (u, v) as [1, 1]
    weights = np.array([1, 1])
    runs = 0
    stop = False
    
    while(stop == False):
        grad_uv = computeGrad(weights)
        grad_uv = np.multiply(grad_uv, 0.1)
        weights = np.subtract(weights, grad_uv)
        runs += 1
        stop = endDescent(weights)
    
    print("Runs: " + str(runs))
    print(weights)

