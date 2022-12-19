#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provides methods to calculate element stiffness matrix

Created on Fri Jun 19 18:56:57 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""
import FEData as model
from utitls import gauss
import numpy as np


def PlateElem(e):
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

	# get gauss points and weights
	w, gp = gauss(model.ngp)

	# compute element stiffness matrix and element nodal force vector
	for i in range(model.ngp):
		for j in range(model.ngp):
			eta = gp[i]
			psi = gp[j]
			# shape functions matrix
			N = NmatPlate(eta, psi)
			# derivative of the shape functions
			B, detJ = BmatPlate(eta, psi)

			# element stiffness matrix
			ke = ke + w[i]*w[j]*detJ*(B.T@model.D@B)
			be = N@(model.b[:, e].reshape((-1, 1)))
			fe = fe + w[i]*w[j]*detJ*(N.T@be)
	return ke, fe


def NmatPlate(eta, psi):
	"""
	Calculate element shape function matrix N at coordinate xt

	Args:
		eta : The first parent coordinate
		psi : The second parent coordinate

	Returns:
		Element shape function matrix N
	"""
	# parent coordinates at nodes
	eta_I = np.array([-1, 1, 1, -1])
	psi_I = np.array([-1, -1, 1, 1])
	
	N = np.zeros((1, 12))
	
	for i in range(model.nen):
		N[0,3*i] = 0.125*( 1 +eta_I[i]*eta)*(1 + psi_I[i]*psi) * \
				(2 + eta_I[i]*eta + psi_I[i]*psi - eta**2 - psi**2)
				
		N[0,3*i+1] = 0.125*( 1 +eta_I[i]*eta)*(1 + psi_I[i]*psi) * \
				(-model.be * psi_I[i] * (1 - psi**2))
				
		N[0,3*i+2] = 0.125*( 1 +eta_I[i]*eta)*(1 + psi_I[i]*psi) * \
				(model.ae * eta_I[i] * (1 - eta**2))

	return N


def BmatPlate(eta, psi):
	"""
	Calcualte derivative of element shape function matrix B at coordinate xt by explicit expression

	Args:
		eta : The first parent coordinate
		psi : The second parent coordinate

	Returns:
		Derivative of element shape function matrix B and Jacobian determination
	"""
	# parent coordinates at nodes
	eta_val = np.array([-1, 1, 1, -1])
	psi_val = np.array([-1, -1, 1, 1])
	
	#Calculate the B_M matrix
	B = np.zeros((3, 12))
	for i in range(model.nen):
		B[:, 3*i:3*i+3] = 1.0 / (4 * model.ae * model.be) * np.array([[ \
								-3*model.be/model.ae*eta_val[i]*eta*(1+psi_val[i]*psi), \
								0, \
								-model.be*eta_val[i]*(1+3*eta_val[i]*eta)*(1+psi_val[i]*psi)], \
								[-3*model.be/model.ae*psi_val[i]*psi*(1+eta_val[i]*eta), \
								model.ae*psi_val[i]*(1+3*psi_val[i]*psi)*(1+eta_val[i]*eta), \
								0], \
								[eta_val[i]*psi_val[i]*(4-3*eta**2-3*psi**2), \
								model.be*eta_val[i]*(3*psi**2+2*psi_val[i]*psi-1), \
								model.ae*psi_val[i]*(1-2*eta_val[i]*eta-3*eta**2)]])

	# Compute Jacobian determination
	detJ = model.ae * model.be

	return B, detJ
