#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot exact solutions for some problems
  ExactSolution_TaperedBar: Plot the exact solution of the tapered elastic bar
  given in Example 5.2 in Fish's textbook
  ExactSolution_CompressionmBar: Plot the exact solution of the bar under
  compression given in Figure 5.13 in Fish's textbook
  ExactSolution_ConcentratedForce: Plot the exact solution of the bar under
  concentrated force and body force

Created on Wed Apr 29 10:26:04 2020

@author: xzhang@tsinghua.edu.cn, thurcni@163.com
"""

import numpy as np
from math import log, sqrt

from utitls import gauss
from Bar1DElem import Nmatrix1D, Bmatrix1D
import FEData as model

    
def ExactSolution_TaperedBar(ax1, ax2):
    """ 
    Plots the exact displacement and stress of the tapered elastic bar
    in axes ax1 and ax2, respectively.
    
    Args:
        ax1 : axis to draw displacement distribution
        ax2 : axis to draw stress distribution
    """
 
    # divide the problem domain into two regions 
    xa = np.arange(2,5,0.01)
    xb = np.arange(5,6,0.01)
         
    # exact displacement for xa 
    c1 = 72;  c2 = 1 - (c1/16)*log(2)
    u1 = -.5*xa + (c1/16)*np.log(xa) + c2
    
    # exact displacement for xb 
    c3 = 48;  c4 = log(5)/16*(c1-c3) + c2
    u2 = -.5*xb + (c3/16)*np.log(xb) + c4
    
    # plot displacement 
    ax1.plot(np.append(xa,xb),np.append(u1,u2), '--r', label='Exact')
     
    # exact stress for xa 
    ya = (36-4*xa)/xa
    
    # exact stress for xb 
    yb = (24-4*xb)/xb
    
    # plot stress 
    ax2.plot(np.append(xa,xb),np.append(ya,yb), '--r', label='Exact')
    
    
def ExactSolution_CompressionBar(ax1, ax2):
    """ 
    Plots the exact displacement and stress of a elastic bar under compression
    in axes ax1 and ax2, respectively.
    
    Args:
        ax1 : axis to draw displacement distribution
        ax2 : axis to draw stress distribution
    """
    xx = np.arange(0, 2, 0.01)

    # exact displacement for a bar under compression
    Ee = 10000
    ue = (-xx**3/6 + xx)/Ee 
    
    # plot displacement 
    ax1.plot(xx, ue, '--r',  label='Exact')
    
    # exact stress
    stre = (-xx**2/2 + 1)
    
    # plot stress 
    ax2.plot(xx,stre, '--r', label='Exact')


def ExactSolution_ConcentratedForce(ax1, ax2):
    """
    Plots the exact displacement and stress of a elastic bar under
    concentrated force and body force in axes ax1 and ax2, respectively.

    Args:
        ax1 : axis to draw displacement distribution
        ax2 : axis to draw stress distribution
    """
    E = 10000.0; A = 1.0
    l = 16.0   # the length of the bar is 2*l
    c = 1.0    # body force
    P1 = 100.0 # The concentrated force at the right end of the bar
    P2 = 100.0 # The concentrated force at x0 = 15l/16
    # divide the problem domain into two regions
    xa = np.arange(0, 15*l/16, 0.01)
    xb = np.arange(15*l/16, 2.0*l, 0.1)

    # exact displacement for xa
    u1 = c*(-xa**3/6 + 2*l**2*xa)/A/E + (P1 + P2)*xa/A/E
    # exact displacement for xb
    u2 = c*(-xb**3/6 + 2*l**2*xb)/A/E + P1*xb/A/E + 15*P2*l/16/A/E

    # plot displacement
    ax1.plot(np.append(xa, xb), np.append(u1, u2), '--r', label='Exact')

    # exact stress for xa
    ya = c*(-xa**2/2 + 2*l**2)/A + (P1 + P2)/A
    # exact stress for xb
    yb = c*(-xb**2/2 + 2*l**2)/A + P1/A

    # plot stress
    ax2.plot(np.append(xa, xb), np.append(ya, yb), '--r', label='Exact')


def ErrorNorm_CompressionBar():
    """ 
    Calculate and print the error norm (L2 and energy norm) of the elastic 
    bar under compression for convergence study
    """
    
    ngp = 3
    [w, gp] = gauss(ngp)    # extract Gauss points and weights
    
    L2Norm = 0
    EnNorm = 0
    
    L2NormEx = 0
    EnNormEx = 0
    
    for e in range(model.nel):
        
        de = model.d[model.LM[:,e]-1] # extract element nodal displacements
        IENe = model.IEN[:,e]-1       # extract local connectivity information
        xe = model.x[IENe]            # extract element x coordinates
        J = (xe[-1] - xe[0])/2        # compute Jacobian
        
        for i in range(ngp):
            xt = 0.5*(xe[0]+xe[-1])+J*gp[i]  # Gauss points in physical coordinates
            
            N = Nmatrix1D(xt,xe)     # shape functions matrix
            B = Bmatrix1D(xt,xe)     # derivative of shape functions matrix
            
            Ee = N@model.E[IENe]     # Young's modulus at element gauss points
            
            uh  = N@de               # displacement at gauss point
            uex = (-xt**3/6 + xt)/Ee # Exact displacement
            L2Norm += J*w[i]*(uex - uh)**2
            L2NormEx += J*w[i]*(uex)**2
            
            sh  = B@de               # strain at Gauss points
            sex = (-xt**2/2 + 1)/Ee  # Exact strain
            EnNorm += 0.5*J*w[i]*Ee*(sex-sh)**2
            EnNormEx += 0.5*J*w[i]*Ee*(sex)**2
    
    L2Norm = sqrt(L2Norm)
    L2NormEx = sqrt(L2NormEx)
    
    EnNorm = sqrt(EnNorm)
    EnNormEx = sqrt(EnNormEx)
    
    # print stresses at element gauss points
    print('\nError norms')
    print('%13s %13s %13s %13s %13s'
          %('h','L2Norm','L2NormRel','EnNorm','EnNormRel'))
    print('%13.6E %13.6E %13.6E %13.6E %13.6E\n'
          %(2/model.nel, L2Norm, L2Norm/L2NormEx, EnNorm, EnNorm/EnNormEx))
    
    return 2/model.nel, L2Norm, EnNorm


def ErrorNorm_ConcentratedForce(flag):
    """
    Calculate and print the error norm (L2 and energy norm) of the elastic
    bar under concentrated force and body force for convergence study

    Args:
        flag: (bool) whether x0 = 15*l/16 is the node or not
    """
    ngp = 3
    [w, gp] = gauss(ngp)  # extract Gauss points and weights

    E = 10000.0
    A = 1.0
    l = 16.0  # the length of the bar is 2*l
    c = 1.0  # body force
    P1 = 100.0  # The concentrated force at the right end of the bar
    P2 = 100.0  # The concentrated force at x0 = 15l/16
    x0 = 15*l/16

    L2Norm = 0
    EnNorm = 0

    L2NormEx = 0
    EnNormEx = 0

    P_element_index = int(x0 / (2*l/model.nel))

    for e in range(model.nel):

        if not flag:
            # The element which contains the concentrated force needs to be
            # handled individually when the force is not on a node
            if e == P_element_index:
                continue

        de = model.d[model.LM[:, e] - 1]  # extract element nodal displacements
        IENe = model.IEN[:, e] - 1  # extract local connectivity information
        xe = model.x[IENe]  # extract element x coordinates
        J = (xe[-1] - xe[0]) / 2  # compute Jacobian

        for i in range(ngp):
            xt = 0.5 * (xe[0] + xe[-1]) + J * gp[i]  # Gauss points in physical coordinates

            N = Nmatrix1D(xt, xe)  # shape functions matrix
            B = Bmatrix1D(xt, xe)  # derivative of shape functions matrix

            Ee = N @ model.E[IENe]  # Young's modulus at element gauss points

            uh = N @ de  # displacement at gauss point
            # Exact displacement
            if xt < x0:
                uex = c*(-xt**3/6 + 2*l**2*xt)/A/E + (P1 + P2)*xt/A/E
            else:
                uex = c*(-xt**3/6 + 2*l**2*xt)/A/E + P1*xt/A/E + 15*P2*l/16/A/E
            L2Norm += J * w[i] * (uex - uh) ** 2
            L2NormEx += J * w[i] * (uex) ** 2

            sh = B @ de  # strain at Gauss points
            # Exact strain
            if xt < x0:
                sex = c*(-xt**2/2 + 2*l**2)/A/E + (P1 + P2)/A/E
            else:
                sex = c*(-xt**2/2 + 2*l**2)/A/E + P1/A/E
            EnNorm += 0.5 * J * w[i] * Ee * (sex - sh) ** 2
            EnNormEx += 0.5 * J * w[i] * Ee * (sex) ** 2

    if not flag:
        # handle the element of P_element_index
        de = model.d[model.LM[:, P_element_index] - 1]  # extract element nodal displacements
        IENe = model.IEN[:, P_element_index] - 1  # extract local connectivity information
        xe = model.x[IENe]  # extract element x coordinates

        for i in range(ngp):
            J = (x0 - xe[0]) / 2  # compute Jacobian
            xt = 0.5 * (xe[0] + x0) + J * gp[i]  # Gauss points in physical coordinates

            N = Nmatrix1D(xt, xe)  # shape functions matrix
            B = Bmatrix1D(xt, xe)  # derivative of shape functions matrix

            Ee = N @ model.E[IENe]  # Young's modulus at element gauss points

            uh = N @ de  # displacement at gauss point
            uex = c * (-xt ** 3 / 6 + 2 * l ** 2 * xt) / A / E + (P1 + P2) * xt / A / E
            L2Norm += J * w[i] * (uex - uh) ** 2
            L2NormEx += J * w[i] * (uex) ** 2

            sh = B @ de  # strain at Gauss points
            sex = c * (-xt ** 2 / 2 + 2 * l ** 2) / A / E+ (P1 + P2) / A / E
            EnNorm += 0.5 * J * w[i] * Ee * (sex - sh) ** 2
            EnNormEx += 0.5 * J * w[i] * Ee * (sex) ** 2

        for i in range(ngp):
            J = (xe[-1] - x0) / 2  # compute Jacobian
            xt = 0.5 * (x0+ xe[-1]) + J * gp[i]  # Gauss points in physical coordinates

            N = Nmatrix1D(xt, xe)  # shape functions matrix
            B = Bmatrix1D(xt, xe)  # derivative of shape functions matrix

            Ee = N @ model.E[IENe]  # Young's modulus at element gauss points

            uh = N @ de  # displacement at gauss point
            uex = c * (-xt ** 3 / 6 + 2 * l ** 2 * xt) / A / E + P1 * xt / A / E + 15 * P2 * l / 16 / A / E
            L2Norm += J * w[i] * (uex - uh) ** 2
            L2NormEx += J * w[i] * (uex) ** 2

            sh = B @ de  # strain at Gauss points
            sex = c * (-xt ** 2 / 2 + 2 * l ** 2) / A / E + P1 / A / E
            EnNorm += 0.5 * J * w[i] * Ee * (sex - sh) ** 2
            EnNormEx += 0.5 * J * w[i] * Ee * (sex) ** 2

    L2Norm = sqrt(L2Norm)
    L2NormEx = sqrt(L2NormEx)

    EnNorm = sqrt(EnNorm)
    EnNormEx = sqrt(EnNormEx)

    # print stresses at element gauss points
    print('\nError norms')
    print('%13s %13s %13s %13s %13s'
          % ('h', 'L2Norm', 'L2NormRel', 'EnNorm', 'EnNormRel'))
    print('%13.6E %13.6E %13.6E %13.6E %13.6E\n'
          % (2 / model.nel, L2Norm, L2Norm / L2NormEx, EnNorm, EnNorm / EnNormEx))

    return 2*l / model.nel, L2Norm, EnNorm
