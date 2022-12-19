#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
'''
Global variables defining the FEM model
  Title: (str) Title of the problem to be solved.
  nsd  : (int) Number of space dimensions.
  ndof : (int) Number of degrees-of-freedom per node.
  nnp  : (int) Number of nodal points.
  nel  : (int) Number of elements.
  nen  : (int) Number of element nodes.
  neq  : (int) Number of equations (D.O.F)
  nd   : (int) Number of nodes on the essential boundary.

  CArea: (numpy.array(nel)) Element values of cross-sectional area.
  E    : (numpy.array(nel)) Element values of Young's modulus.
  leng : (numpy.array(nel)) Element values of length
  stress:(numpy.array(nel)) Element values of stress

  x    : (numpy.array(nnp))x coordinate.
  y    : (numpy.array(nnp))y coordinates.
  IEN  : (numpy.array(nen,nel)) Element connectivity array.
  LM   : (numpy.array(nen,nel)) Location matrix.
  K    : (numpy.array(neq,neq)) Global stiffness matrix
  f    : (numpy.array(neq,1)) Global nodal force vector            
  d    : (numpy.array(neq,1)) Solution vector

  plot_truss: (bool) Plot truss ?
  plot_node : plot node number ?
  plot_tex  : Convert figures into PGFPlots figures in LaTex file ?

Created on Sat May 9 18:34:00 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
'''

Title = ""
nsd   = 0
ndof  = 0
nnp   = 0
nel   = 0
nen   = 0
neq   = 0
nd    = 0

CArea = np.array([])
E     = np.array([])
leng  = np.array([])
stress= np.array([])

x     = np.array([])
y     = np.array([])
IEN   = np.array([[]])
LM    = np.array([[]])
K     = np.array([[]])
f     = np.array([[]])
d     = np.array([[]])

plot_truss = False
plot_node  = False
plot_tex   = False