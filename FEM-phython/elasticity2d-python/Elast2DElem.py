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


def Elast2DElem(e):
	"""
	Calculate element stiffness matrix and element nodal body force vector

	Args:
		e : (int) element number

	Returns: ke, fe
		ke : (numpy(nen,nen)) element stiffness matrix
		fe : (numpy(nen,1)) element nodal force vector
	"""
	ke = np.zeros((model.nen*model.ndof, model.nen*model.ndof))
	fe = np.zeros((model.nen*model.ndof, 1))

	# get coordinates of element nodes
	je = model.IEN[:, e] - 1
	C = np.array([model.x[je], model.y[je]]).T

	ngpf = 2  # No. of gauss points for force calculation
	if model.ngp == 2:	# Full integration
		ngpk = 2   # No. of gauss points for stiffness calculation
	elif abs(model.ngp) == 1:	# Reduced integration with/without hourglass control
		ngpk = 1
	else:
		print("Error : Invalid value of ngp ({}) !".format(model.ngp))

	# get gauss points and weights
	wk, gpk = gauss(abs(ngpk))
	wf, gpf = gauss(abs(ngpf))

	# compute element stiffness matrix
	for i in range(ngpk):
		for j in range(ngpk):
			eta = gpk[i]
			psi = gpk[j]

			# derivative of the shape functions
			B, detJ = BmatElast2D(eta, psi, C)

			# element stiffness matrix
			ke = ke + wk[i]*wk[j]*detJ*(B.T@model.D@B)

	# Reduced integration with hourglass control
	if model.ngp == 1:		
		B, detJ = BmatElast2D(0.0, 0.0, C)		
		Ae = 4.0 * detJ  # Area of element

		hs_vectors = np.zeros((model.nen, 1))	# Hourglass shape vector (4*1)
		hs_vectors[0] = (C[1,0]*(C[2,1]-C[3,1]) + C[2,0]*(C[3,1]-C[1,1]) + C[3,0]*(C[1,1]-C[2,1])) / Ae
		hs_vectors[1] = (C[2,0]*(C[0,1]-C[3,1]) + C[3,0]*(C[2,1]-C[0,1]) + C[0,0]*(C[3,1]-C[2,1])) / Ae
		hs_vectors[2] = (C[3,0]*(C[0,1]-C[1,1]) + C[0,0]*(C[1,1]-C[3,1]) + C[1,0]*(C[3,1]-C[0,1])) / Ae
		hs_vectors[3] = (C[0,0]*(C[2,1]-C[1,1]) + C[1,0]*(C[0,1]-C[2,1]) + C[2,0]*(C[1,1]-C[0,1])) / Ae
		
		DNNorm = 0.0
		DN = DNmatElast2D(0.0, 0.0, C)
		for i in range(model.ndof):
			for I in range(model.nen):
				DNNorm = DNNorm + DN[i,I] * DN[i,I]
		
		hourglass_k = 0.01 * model.G * DNNorm

		for I in range(model.nen):
			for i in range(model.ndof):
				a = model.ndof * I + i
				for J in range(model.nen):
						b = model.ndof * J + i
						ke[a,b] = ke[a,b] + hs_vectors[I] * hs_vectors[J] * hourglass_k * Ae    

	# compute element nodal force vector
	for i in range(ngpf):
		for j in range(ngpf):
			eta = gpf[i]
			psi = gpf[j]

			# shape functions matrix
			N = NmatElast2D(eta, psi)
			# derivative of the shape functions
			B, detJ = BmatElast2D(eta, psi, C)

			# element nodal force vector
			be = N@(model.b[:, e].reshape((-1, 1)))
			fe = fe + wf[i]*wf[j]*detJ*(N.T@be)

	return ke, fe


def NmatElast2D(eta, psi):
	"""
	Calculate element shape function matrix N at coordinate xt

	Args:
		eta : The first parent coordinate
		psi : The second parent coordinate

	Returns:
		Element shape function matrix N
	"""
	N1 = 0.25*(1-psi)*(1-eta)
	N2 = 0.25*(1+psi)*(1-eta)
	N3 = 0.25*(1+psi)*(1+eta)
	N4 = 0.25*(1-psi)*(1+eta)

	return np.array([[N1, 0, N2, 0, N3, 0, N4, 0],
					 [0, N1, 0, N2, 0, N3, 0, N4]])


def BmatElast2D(eta, psi, C):
	"""
	Calcualte derivative of element shape function matrix B at coordinate xt

	Args:
		eta : The first parent coordinate
		psi : The second parent coordinate
		C   : The physical coordinates

	Returns:
		Derivative of element shape function matrix B and Jacobian determination
	"""
	#Calculate the Grad(N) matrix
	GN = 0.25*np.array([[eta-1, 1-eta, 1+eta, -eta-1],
						[psi-1, -psi-1, 1+psi, 1-psi]])

	# Compute Jacobian matrix
	J = GN@C
	detJ = np.linalg.det(J)

	BB = np.linalg.solve(J, GN)
	B1x = BB[0, 0]
	B2x = BB[0, 1]
	B3x = BB[0, 2]
	B4x = BB[0, 3]
	B1y = BB[1, 0]
	B2y = BB[1, 1]
	B3y = BB[1, 2]
	B4y = BB[1, 3]

	B = np.array([[B1x, 0, B2x, 0, B3x, 0, B4x, 0],
				  [0, B1y, 0, B2y, 0, B3y, 0, B4y],
				  [B1y, B1x, B2y, B2x, B3y, B3x, B4y, B4x]])

	return B, detJ

def DNmatElast2D(eta, psi, C):
	"""
	Calcualte derivative of element shape function matrix DN at coordinate xt

	Args:
		eta : The first parent coordinate
		psi : The second parent coordinate
		C   : The physical coordinates

	Return:
		Derivative of element shape function matrix BB
	"""

	#Calculate the Grad(N) matrix in parent coordinate
	GN = 0.25*np.array([[eta-1, 1-eta, 1+eta, -eta-1],
						[psi-1, -psi-1, 1+psi, 1-psi]])

	# Compute Jacobian matrix
	J = GN@C

	DN = np.linalg.solve(J, GN)

	return DN
