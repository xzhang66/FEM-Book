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
  ngp  : (int) Number of gauss points.
  
  lx : (float) Length of plate in the x-direction.
  ly : (float) Length of plate element in the y-direction.
  nelx : (int) Number of elements in the x-direction.
  nely : (int) Number of elements in the y-direction.
  nenx : (int) Number of nodes in the x-direction (= nelx + 1).
  neny : (int) Number of nodes in the y-direction (= nely + 1).
  ae   : (float) Half of the length of element in the x-direction.
  be   : (float) Half of the length of element in the y-direction.
  h    : (float) thickness of the plate.
  q    : (float) The uniform load.
  
  flags: (numpy.array(neq))  Nodal boundary condition flag:
         2 - located on the essential boundary;
         1 - located on the natural boundary.
  nd   : (int) Number of nodes on the essential boundary.
  e_bc : (numpy.array(neq)) Value of essential B.C.

  nbe  : (int) Number of prescribed traction edges.
  n_bc : (numpy.array(nnp)) Value of natural B.C.
  
  P    : (numpy.array(neq)) Array of nodal external forces.
  b    : (numpy.array(nen*ndof, nel)) Element nodal forces.
  D	   : (numpy.array(3, 3)) elasticity matrix
  E	   : (float) Young's modulus
  nu   : (float) Poisson's ratio
  
  IEN  : (numpy.array(nen,nel)) Element connectivity array.
  ID   : (numpy.array(neq) Identification matrix.
  LM   : (numpy.array(nen,nel)) Location matrix.
  x    : (numpy.array(nnp))x coordinate.
  y    : (numpy.array(nnp))y coordinates, used only for the bar plot.

  K    : (numpy.array(neq,neq)) Global stiffness matrix
  f    : (numpy.array(neq,1)) Global nodal force vector            
  d    : (numpy.array(neq,1)) Solution vector

  plot_mesh       : plot mesh ?
  plot_nod        : plot node number ?
  plot_tex        : plot in latex tikz format ?
  nplot           : (int) Number of points in a element used to plot displacements
                    and moment (20).
"""

Title = None
nsd = 0
ndof = 0
nnp = 0
nel = 0
nen = 0
neq = 0
ngp = 0
nd = 0
nbe = 0

f = None
d = None
K = None

# geometric data
lx = 0.0
ly = 0.0
nelx = 0
nely = 0
nenx = 0
neny = 0
h = 0.0
ae = 0.0
be = 0.0

# boundary conditions
flags = None
e_bc = None
n_bc = None

# force conditions
P = None
b = None
q = 0.0

# material
D = None
E = 0.0
nu = 0.0

# define the mesh
x = None
y = None
IEN = None

ID = None
LM = None

# parameter for postprocess
plot_mesh = None
plot_nod = None
plot_tex = None
nplot = 20
