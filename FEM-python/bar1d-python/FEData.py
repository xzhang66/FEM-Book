#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" 
Global variables defining the FEM model
  Title: (str) Title of the problem to be solved.
  nsd  : (int) Number of space dimensions.
  ndof : (int) Number of degrees-of-freedom per node.
  nnp  : (int) Number of nodal points.
  nel  : (int) Number of elements.
  nen  : (int) Number of element nodes.
  neq  : (int) Number of equations (D.O.F)

  E    : (numpy.array(nnp)) Nodal values Young's modulus.
  body : (numpy.arraay(nnp)) Nodal values body forces.
  CArea: (numpy.array(nnp)) Nodal values of cross-sectional area.
  ngp  : (int) Number of gauss points.

  flags: (numpy.array(neq))  Nodal DOF boundary condition flag:
         2 - located on the essential boundary;
         1 - located on the natural boundary.
  nd   : (int) Number of nodes on the essential boundary.
  e_bc : (numpy.array(neq)) Value of essential B.C.
  n_bc : (numpy.array(neq)) Value of natural B.C.

  np   : (int) Number of point forces.
  xp   : ((numpy.array(np))) Array of coordinates where point forces are applied.
  P    : (numpy.array(np)) Array of point forcess.

  plot_bar: (bool) Plot bar ?
  plot_nod: plot node number ?
  plot_tex : Convert figures into PGFPlots figures in LaTex file ?
  nplot: (int) Number of points in a element used to plot displacements 
         and stresses (10*nen).

  x    : (numpy.array(nnp))x coordinate.
  y    : (numpy.array(nnp))y coordinates, used only for the bar plot.
  IEN  : (numpy.array(nen,nel)) Element connectivity array.
  ID   : (numpy.array(neq) Identification matrix.
  LM   : (numpy.array(nen,nel)) Location matrix.

  K    : (numpy.array(neq,neq)) Global stiffness matrix
  f    : (numpy.array(neq,1)) Global nodal force vector            
  d    : (numpy.array(neq,1)) Solution vector

  Exact: Problem type whose exact solution is available  
         TaperedBar       : The tapered elastic bar in Example 5.2 (Fish's book)
         CompressionBar   : The bar under compression in Figure 5.13 (Fish's book)
         ConcentratedForce: The bar under concentrated force and body force

Important notices
  
  All indices in numpy.array are zero-based, but user provided FE model 
  is one-based numbered. For example, the x-coordinate of node 1 is stored
  in x[0]. Therefore, to obtain the nodal coordinates of element e, we use:
    xe = x[IEN[:,e]-1]

Created on Sun Apr 24 18:56:57 2020

@author: xzhang@tsinghua.edu.cn
"""

Title= None
nsd  = 2  
ndof = 1
nnp  = 0
nel  = 0
nen  = 2

neq  = 0
f = None            
d = None        
K = None   

# boundary conditions
flags= None        
e_bc = None        
n_bc = None

# element and material data (given at the element nodes)
E     = None
body  = None
CArea = None

# gauss integration
ngp = 0

# boundary conditions
flags = None
e_bc = None 
n_bc = None  
nd = 0

# point forces
np = None
xp = None
P  = None   

# output plots
plot_bar = None
plot_nod = None
nplot = 0
plot_tex = None

# define the mesh
x = None
y = None  
IEN = None

ID  = None
LM  = None   

Exact = None