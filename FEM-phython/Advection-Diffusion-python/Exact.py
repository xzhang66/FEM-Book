#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot exact solutions for the ddvection-diffusion Problem

Created on Wed Apr 29 10:26:04 2020

@author: xzhang@tsinghua.edu.cn, thurcni@163.com
"""

import numpy as np
from math import exp

import FEData as model
import matplotlib.pyplot as plt

    
def ExactSolution():
    """ 
    Plots the exact displacement distribution.
    """

    xa = np.arange(0,10,0.01) #1000*1
    ka = model.k[0]
    lene = 0.5
    v = -2.0*ka*model.PN/lene;
    if abs(v) > 300 :
        A = 0.0
    else :
        A = 1.0/(exp(-10.0*v/ka)-1.0)

    ya = np.zeros((1000,1))
    # exact displacement for xa
    for i in range(1000):
        if abs(v) > 300 :
            ya[i] = -A
        else :
            ya[i] = A*exp(-v*xa[i]/ka)-A #1000*1
    
    # plot displacement 
    line2, = plt.plot(xa,ya, '--r', label='Exact')