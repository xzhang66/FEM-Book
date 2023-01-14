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
  D    : (numpy.array(5, 5)) elasticity matrix
  E    : (float) Young's modulus
  ne   : (float) Poisson's ratio
  G	   : (float) shear modulus
  
  IEN  : (numpy.array(nen,nel)) Element connectivity array.
  ID   : (numpy.array(neq) Identification matrix.
  LM   : (numpy.array(nen,nel)) Location matrix.
  xb   : (numpy.array(nnp)) x coordinates of nodes at the bottom surface.
  yb   : (numpy.array(nnp)) y coordinates of nodes at the bottom surface.
  zb   : (numpy.array(nnp)) z coordinates of nodes at the bottom surface.
  xt   : (numpy.array(nnp)) x coordinates of nodes at the top surface.
  yt   : (numpy.array(nnp)) y coordinates of nodes at the top surface.
  zt   : (numpy.array(nnp)) z coordinates of nodes at the top surface.
  xI   : (numpy.array(nnp)) x coordinates of nodes on the midsurface.
  yI   : (numpy.array(nnp)) y coordinates of nodes on the midsurface.
  zI   : (numpy.array(nnp)) z coordinates of nodes on the midsurface.
  V3   : (numpy.array(3, nnp)) Vetors pointing from the bottom nodes to the top nodes.
  t    : (numpy.array(nnp)) Thickness of nodes on the midsurface.
  v3   : (numpy.array(3, nnp)) Unit vector normal to the midsurface.
  v1   : (numpy.array(3, nnp)) Unit vector normal to v3.
  v2   : (numpy.array(3, nnp)) Unit vector normal to v3.

  K    : (numpy.array(neq,neq)) Global stiffness matrix
  f    : (numpy.array(neq,1)) Global nodal force vector
  d    : (numpy.array(neq,1)) Solution vector
  w_c  : (float) Center deflection

  plot_mesh       : plot mesh ?
  plot_nod        : plot node number ?
  plot_centerline : plot deflection and moment Mx distributions along centerline ?
  plot_tex        : plot in latex tikz format ?
  nplot           : (int) Number of points in a element used to plot displacements(20).
"""

Title = None
nsd = 0
ndof = 0
nnp = 0
nel = 0
nen = 0
neq = 0
ngp = 1
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

xb = None
yb = None
zb = None
xt = None
yt = None
zt = None
xI = None
yI = None
zI = None
V3 = None
t  = None
v3 = None
v1 = None
v2 = None

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
ne = 0.0
G = 0.0

# define the mesh
IEN = None

ID = None
LM = None

# parameter for postprocess
plot_mesh = None
plot_nod = None
plot_tex = None
plot_centerline = None
nplot = 20

# center deflection
w_c = 0.0
