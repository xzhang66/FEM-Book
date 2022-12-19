#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provides methods to calculate element stiffness matrix

Created on Fri Jun 19 18:56:57 2020

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""

import FEData as model
from utitls import gauss
import numpy as np


def MindlinPlateElem(e):
	"""
	Calculate element stiffness matrix and element nodal body force vector

	Args:
		e : (int) element number

	Returns: ke, fe
		ke : (numpy(nen*ndof,nen*ndof)) element stiffness matrix
		fe : (numpy(nen*ndof,1)) element nodal force vector
	"""
	ke = np.zeros((model.nen*model.ndof, model.nen*model.ndof))
	fe = np.zeros((model.nen*model.ndof, 1))

	# get coordinates of element nodes
	je = model.IEN[:, e] - 1
	C = np.array([model.x[je], model.y[je]]).T

	ngpf = 2              # No. of gauss points for force calculation
	if model.ngp == 2:    # Full integration
		ngpb = model.ngp      # No. of gauss points for bending stiffness calculation
		ngps = model.ngp      # No. of gauss points for shear stiffness calculation
	elif model.ngp == 1:    # Selective reduced integration
		ngpb = model.ngp + 1  # No. of gauss points for bending stiffness calculation
		ngps = model.ngp      # No. of gauss points for shear stiffness calculation
	else:
		print("Error : Invalid value of ngp ({}) !".format(model.ngp))

	# get gauss points and weights
	wb, gpb = gauss(ngpb)
	ws, gps = gauss(ngps)
	wf, gpf = gauss(ngpf)

	# compute element bending stiffness matrix
	for i in range(ngpb):
		for j in range(ngpb):
			eta = gpb[i]
			psi = gpb[j]

			# derivative of the shape functions
			Bb, Bs, detJ = BmatMindlinPlate(eta, psi, C)

			# element stiffness matrix
			ke = ke + wb[i]*wb[j]*detJ*Bb.T@model.Db@Bb
	
	# compute element bending stiffness matrix
	for i in range(ngps):
		for j in range(ngps):
			eta = gps[i]
			psi = gps[j]

			# derivative of the shape functions
			Bb, Bs, detJ = BmatMindlinPlate(eta, psi, C)

			# element stiffness matrix
			ke = ke + ws[i]*ws[j]*detJ*Bs.T@model.Ds@Bs
	
	# compute element nodal force vector
	for i in range(ngpf):
		for j in range(ngpf):
			eta = gpf[i]
			psi = gpf[j]

			# shape functions matrix
			N = NmatMindlinPlate(eta, psi)
			# derivative of the shape functions
			Bb, Bs, detJ = BmatMindlinPlate(eta, psi, C)

			# element nodal force vector
			be = N@(model.b[:, e].reshape((-1, 1)))
			fe = fe + wf[i]*wf[j]*detJ*(N.T@be)

	return ke, fe


def NmatMindlinPlate(eta, psi):
	"""
	Calculate element shape function matrix N at coordinate xt

	Args:
		eta : The first parent coordinate
		psi : The second parent coordinate

	Returns:
		Element shape function matrix N
	"""
	N1 = 0.25*(1-eta)*(1-psi)
	N2 = 0.25*(1+eta)*(1-psi)
	N3 = 0.25*(1+eta)*(1+psi)
	N4 = 0.25*(1-eta)*(1+psi)

	return np.array([[N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0, 0],
					 [0, N1, 0, 0, N2, 0, 0, N3, 0, 0, N4, 0],
					 [0, 0, N1, 0, 0, N2, 0, 0, N3, 0, 0, N4]])


def BmatMindlinPlate(eta, psi, C):
	"""
	Calcualte bending kinematic matrix Bb and shear kinematic matrix Bs at coordinate xt

	Args:
		eta : The first parent coordinate
		psi : The second parent coordinate
		C   : The physical coordinates

	Returns:
		bending kinematic matrix Bb, shear kinematic matrix Bs and Jacobian determination
	"""
	#Calculate the Grad(N) matrix
	GN = 0.25*np.array([[psi-1, 1-psi, 1+psi, -psi-1],
						[eta-1, -eta-1, 1+eta, 1-eta]])

	# Compute Jacobian matrix
	J = GN@C
	detJ = np.linalg.det(J)

	DN = np.linalg.solve(J, GN)
	DN1x = DN[0, 0]
	DN2x = DN[0, 1]
	DN3x = DN[0, 2]
	DN4x = DN[0, 3]
	DN1y = DN[1, 0]
	DN2y = DN[1, 1]
	DN3y = DN[1, 2]
	DN4y = DN[1, 3]

	Bb = np.array([[DN1x, 0, 0, DN2x, 0, 0, DN3x, 0, 0, DN4x, 0, 0],
				   [0, DN1y, 0, 0, DN2y, 0, 0, DN3y, 0, 0, DN4y, 0],
				   [DN1y, DN1x, 0, DN2y, DN2x, 0, DN3y, DN3x, 0, DN4y, DN4x, 0]])
	
	# shape functions matrix
	N1 = 0.25*(1-eta)*(1-psi)
	N2 = 0.25*(1+eta)*(1-psi)
	N3 = 0.25*(1+eta)*(1+psi)
	N4 = 0.25*(1-eta)*(1+psi)
	
	Bs = np.array([[-N1, 0, DN1x, -N2, 0, DN2x, -N3, 0, DN3x, -N4, 0, DN4x],
				   [0, -N1, DN1y, 0, -N2, DN2y, 0, -N3, DN3y, 0, -N4, DN4y]])

	return Bb, Bs, detJ