#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provides utilities used by FE analysis.
  1. gauss: Gauss quadrature rules.
  2. assembly: Global stiffness matrix and nodal force vector assembly.
  3. solvedr: Solving the stiffness equations by the reduction approach.

Created on Sun Apr 24 18:56:57 2020

@author: xzhang@tsinghua.edu.cn
"""

import numpy as np
import FEData as model

def gauss(ngp):
    """
    Get Gauss points in the parent element domain [-1, 1] and 
    the corresponding weights.
    
    Args: 
        ngp : (int) number of Gauss points.
        
    Returns: w,gp
        w  : weights.
        gp : Gauss points in the parent element domain.
    """
    
    if ngp == 1:
        gp = [0]
        w = [2]
    elif ngp == 2:
        gp = [-np.sqrt(3)/3, np.sqrt(3)/3]
        w = [1, 1]
    elif ngp == 3:
        gp = [-np.sqrt(3/5), np.sqrt(3/5), 0.0]
        w = [5/9, 5/9, 8/9]
    elif ngp == 4:
        gp = [-np.sqrt(525+70*np.sqrt(30))/35, np.sqrt(525+70*np.sqrt(30))/35, 
              -np.sqrt(525-70*np.sqrt(30))/35, np.sqrt(525-70*np.sqrt(30))/35]
        w = [(18-np.sqrt(30))/36, (18-np.sqrt(30))/36, 
             (18+np.sqrt(30))/36, (18+np.sqrt(30))/36]
    elif ngp == 5:
        gp = [0, 
              -np.sqrt(245-14*np.sqrt(70))/21, np.sqrt(245-14*np.sqrt(70))/21,
              -np.sqrt(245+14*np.sqrt(70))/21, np.sqrt(245+14*np.sqrt(70))/21]
        w = [128/225, 
             (322+13*np.sqrt(70))/900, (322+13*np.sqrt(70))/900,
             (322-13*np.sqrt(70))/900, (322-13*np.sqrt(70))/900
             ]
    else:
        print('\n Invalid Gauss quadrature rule: ngp = %1d (must be 1~5).'
              %(ngp))
        raise SystemExit
    
    return w, gp


def assembly(e,ke,fe):
    """
    Assemble element stiffness matrix and nodal force vector.
    
    Args:
        e   : (int) Element number
        ke  : (numpy(nen,nen)) element stiffness matrix
        fe  : (numpy(nen,1)) element nodal force vector
    """
    for loop1 in range(model.nen):
       i = model.LM[loop1,e]-1
       model.f[i] += fe[loop1]   # assemble nodal force vector
       
       for loop2 in range(model.nen):
          j = model.LM[loop2,e]-1
          model.K[i,j] += ke[loop1,loop2]   # assemble stiffness matrix


def solvedr():
    """
    Partition and solve the system of equations
        
    Returns:
        f_E : (numpy.array(nd,1)) Reaction force vector
    """
    nd = model.nd; neq=model.neq
    K_E = model.K[0:nd, 0:nd]
    K_F = model.K[nd:neq, nd:neq]
    K_EF =model. K[0:nd, nd:neq]
    f_F = model.f[nd:neq]
    d_E = model.d[0:nd]
    
    # solve for d_F
    d_F = np.linalg.solve(K_F, f_F - K_EF.T @ d_E) 

    # reconstruct the global displacement d
    model.d = np.append(d_E,d_F)
    
    # compute the reaction r
    f_E = K_E@d_E + K_EF@d_F
    
    # write to the workspace
    print('\nsolution d');  print(model.d)
    print('\nreaction f =', f_E)
    
    return f_E
