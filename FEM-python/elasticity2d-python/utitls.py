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
	gp = None
	w = None
	if ngp == 1:
		gp = [0]
		w = [2]
	elif ngp == 2:
		gp = [-0.57735027, 0.57735027]
		w = [1, 1]
	elif ngp == 3:
		gp = [-0.7745966692, 0.7745966692, 0.0]
		w = [0.5555555556, 0.5555555556, 0.8888888889]
	else:
		raise ValueError("The given number (ngp = {}) of Gauss points is too large and not implemented".format(ngp))
	return w, gp


def assembly(e, ke, fe):
	"""
	Assemble element stiffness matrix and nodal force vector.

	Args:
		e   : (int) Element number
		ke  : (numpy(nen,nen)) element stiffness matrix
		fe  : (numpy(nen,1)) element nodal force vector
	"""
	for loop1 in range(model.nen*model.ndof):
		i = model.LM[loop1, e]-1
		model.f[i] += fe[loop1]   # assemble nodal force vector

		for loop2 in range(model.nen*model.ndof):
			j = model.LM[loop2, e]-1
			model.K[i, j] += ke[loop1, loop2]   # assemble stiffness matrix


def solvedr():
	"""
	Partition and solve the system of equations

	Returns:
		f_E : (numpy.array(nd,1)) Reaction force vector
	"""
	nd = model.nd
	neq = model.neq
	K_E = model.K[0:nd, 0:nd]
	K_F = model.K[nd:neq, nd:neq]
	K_EF = model. K[0:nd, nd:neq]
	f_F = model.f[nd:neq]
	d_E = model.d[0:nd]

	if (neq > nd):
		print('\nCondition number of stiffness matrix: ', np.linalg.cond(K_F))

	# solve for d_F
	d_F = np.linalg.solve(K_F, f_F - K_EF.T @ d_E)

	# reconstruct the global displacement d
	model.d = np.append(d_E,d_F)

	# compute the reaction r
	f_E = K_E@d_E + K_EF@d_F - model.f[:model.nd]

	# write to the workspace
	print('\nsolution d')
	print(model.d)
	print('\nreaction f = \n', f_E)

	return f_E
