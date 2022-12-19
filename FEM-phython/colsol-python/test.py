#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test colsol function

Created on Fri Jun 19 14:06:31 2020

@author: xzhang@tsinghua.edu.cn
"""

import json
import numpy as np
from colsol import colsol

# Read in equations from JSON file
with open('Example_n_5.json') as f_obj:
        Equations = json.load(f_obj)

n = Equations['n']    # Number of equations
m = np.array(Equations['m'])    # The row number of the first nonzero element in each column
K = np.array(Equations['K'], dtype='f')    # The stiffness matrix
R = np.array(Equations['R'], dtype='f')    # The Right-hand-side load vector

[IERR, K, R] = colsol(n, m, K, R)

if IERR==0:
    print("K = \n", K)
    print("R = \n", R)