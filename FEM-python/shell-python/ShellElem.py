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


def ShellElem(e):
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
	C = np.array([model.xI[je], model.yI[je], model.zI[je]]).T # 8*3
	V3 = model.V3[:, je] # 3*8

	ngpx = model.ngp       # No. of gauss points in the xi direction
	ngpy = model.ngp       # No. of gauss points in the eta direction
	ngpz = 2               # No. of gauss points in the zeta direction

	# get gauss points and weights
	wx, gpx = gauss(ngpx)
	wy, gpy = gauss(ngpy)
	wz, gpz = gauss(ngpz)

	# compute element stiffness matrix
	for i in range(ngpx):
		for j in range(ngpy):
			for m in range(ngpz):
				xi = gpx[i]
				eta = gpy[j]
				zeta = gpz[m]

				# derivative of the shape functions
				B, detJ = BmatShell(xi, eta, zeta, C, V3)

				# element stiffness matrix
				ke = ke + wx[i]*wy[j]*wz[m]*detJ*B.T@model.D@B
	
	# compute element nodal force vector
	for i in range(ngpx):
		for j in range(ngpy):
			for m in range(ngpz):
				xi = gpx[i]
				eta = gpy[j]
				zeta = gpz[m]

				# shape functions matrix
				N = NmatShell(xi, eta, zeta)
				# derivative of the shape functions
				B, detJ = BmatShell(xi, eta, zeta, C, V3)

				# element nodal force vector
				be = N@(model.b[:, e].reshape((-1, 1)))
				fe = fe + wx[i]*wy[j]*wz[m]*detJ*(N.T@be)

	return ke, fe


def NmatShell(xi, eta, zeta):
	"""
	Calculate element shape function matrix N at coordinate xt

	Args:
		xi : The first parent coordinate
		eta : The second parent coordinate
		zeta : The third parent coordinate

	Returns:
		Element shape function matrix N 
	"""
	NI = np.zeros((model.nen, 1))
	NI[0, 0] = -0.25*(1-xi)*(1-eta)*(xi+eta+1)
	NI[1, 0] = 0.25*(1+xi)*(1-eta)*(xi-eta-1)
	NI[2, 0] = 0.25*(1+xi)*(1+eta)*(xi+eta-1)
	NI[3, 0] = -0.25*(1-xi)*(1+eta)*(xi-eta+1)
	NI[4, 0] = 0.5*(1-xi**2)*(1-eta)
	NI[5, 0] = 0.5*(1+xi)*(1-eta**2)
	NI[6, 0] = 0.5*(1-xi**2)*(1+eta)
	NI[7, 0] = 0.5*(1-xi)*(1-eta**2)

	N = np.zeros((model.nsd, model.nen*model.ndof)) # 3*40
	for i in range(model.nen):
		N[:, model.ndof*i : model.ndof*(i+1)] = np.array([ \
					[NI[i,0], 0, 0, \
					0.5*zeta*model.t[i,0]*NI[i,0]*model.v1[0,i], \
					-0.5*zeta*model.t[i,0]*NI[i,0]*model.v2[0,i]],
					[0, NI[i,0], 0, \
					0.5*zeta*model.t[i,0]*NI[i,0]*model.v1[1,i], \
					-0.5*zeta*model.t[i,0]*NI[i,0]*model.v2[1,i]],
					[0, 0, NI[i,0], \
					0.5*zeta*model.t[i,0]*NI[i,0]*model.v1[2,i], \
					-0.5*zeta*model.t[i,0]*NI[i,0]*model.v2[2,i]]])
	
	return N

def BmatShell(xi, eta, zeta, C, V3):
	"""
	Calcualte derivative of element shape function matrix B at coordinate xt

	Args:
		xi : The first parent coordinate
		eta : The second parent coordinate
		zeta : The third parent coordinate
		C   : The physical coordinates
		V3  : The vetors pointing from the bottom nodes to the top nodes

	Returns:
		Derivative of element shape function matrix B and Jacobian determination
	"""
	# shape function NI
	N = np.zeros((model.nen, 1))
	N[0, 0] = -0.25*(1-xi)*(1-eta)*(xi+eta+1)
	N[1, 0] = 0.25*(1+xi)*(1-eta)*(xi-eta-1)
	N[2, 0] = 0.25*(1+xi)*(1+eta)*(xi+eta-1)
	N[3, 0] = -0.25*(1-xi)*(1+eta)*(xi-eta+1)
	N[4, 0] = 0.5*(1-xi**2)*(1-eta)
	N[5, 0] = 0.5*(1+xi)*(1-eta**2)
	N[6, 0] = 0.5*(1-xi**2)*(1+eta)
	N[7, 0] = 0.5*(1-xi)*(1-eta**2)

	# derivative of shape function NI in parent coordinates
	DNDxi = np.zeros((model.nen, 1))
	DNDxi[0, 0] = 0.25*(2*xi+eta)*(1-eta)
	DNDxi[1, 0] = 0.25*(2*xi-eta)*(1-eta)
	DNDxi[2, 0] = 0.25*(2*xi+eta)*(1+eta)
	DNDxi[3, 0] = 0.25*(2*xi-eta)*(1+eta)
	DNDxi[4, 0] = -xi*(1-eta)
	DNDxi[5, 0] = 0.5*(1-eta**2)
	DNDxi[6, 0] = -xi*(1+eta)
	DNDxi[7, 0] = -0.5*(1-eta**2)
	
	DNDeta = np.zeros((model.nen, 1))
	DNDeta[0, 0] = 0.25*(xi+2*eta)*(1-xi)
	DNDeta[1, 0] = 0.25*(-xi+2*eta)*(1+xi)
	DNDeta[2, 0] = 0.25*(xi+2*eta)*(1+xi)
	DNDeta[3, 0] = 0.25*(-xi+2*eta)*(1-xi)
	DNDeta[4, 0] = -0.5*(1-xi**2)
	DNDeta[5, 0] = (1+xi)*(-eta)
	DNDeta[6, 0] = 0.5*(1-xi**2)
	DNDeta[7, 0] = (1-xi)*(-eta)
	
	# compute Jacobian matrix J
	DXDxi = np.zeros((model.nsd, 1))
	DXDeta = np.zeros((model.nsd, 1))
	DXDzeta = np.zeros((model.nsd, 1))
	for i in range(model.nen):
		DXDxi[0, 0] += DNDxi[i, 0] * (C[i, 0] + 0.5*zeta*V3[0, i])
		DXDxi[1, 0] += DNDxi[i, 0] * (C[i, 1] + 0.5*zeta*V3[1, i])
		DXDxi[2, 0] += DNDxi[i, 0] * (C[i, 2] + 0.5*zeta*V3[2, i])
		DXDeta[0, 0] += DNDeta[i, 0] * (C[i, 0] + 0.5*zeta*V3[0, i])
		DXDeta[1, 0] += DNDeta[i, 0] * (C[i, 1] + 0.5*zeta*V3[1, i])
		DXDeta[2, 0] += DNDeta[i, 0] * (C[i, 2] + 0.5*zeta*V3[2, i])
		DXDzeta[0, 0] += N[i, 0] * 0.5*V3[0, i]
		DXDzeta[1, 0] += N[i, 0] * 0.5*V3[1, i]
		DXDzeta[2, 0] += N[i, 0] * 0.5*V3[2, i]
	
	J = np.array(([DXDxi[0,0], DXDxi[1,0], DXDxi[2,0]], \
				 [DXDeta[0,0], DXDeta[1,0], DXDeta[2,0]], \
				 [DXDzeta[0,0], DXDzeta[1,0], DXDzeta[2,0]])) # 3*3
	detJ = np.linalg.det(J)
	J_inv = np.linalg.inv(J)

	# derivative of shape function matrix N in parent coordinates
	DuDxi = np.zeros((model.nsd*model.nsd, model.nen*model.ndof)) # 9*40
	for i in range(model.nen):
		DuDxi[0, model.ndof*i:model.ndof*(i+1)] = np.array([ \
					DNDxi[i,0], 0, 0, \
					0.5*zeta*model.t[i,0]*DNDxi[i,0]*model.v1[0,i], \
					-0.5*zeta*model.t[i,0]*DNDxi[i,0]*model.v2[0,i]])
		DuDxi[1, model.ndof*i:model.ndof*(i+1)] = np.array([ \
					DNDeta[i,0], 0, 0, \
					0.5*zeta*model.t[i,0]*DNDeta[i,0]*model.v1[0,i], \
					-0.5*zeta*model.t[i,0]*DNDeta[i,0]*model.v2[0,i]])
		DuDxi[2, model.ndof*i:model.ndof*(i+1)] = np.array([ \
					0, 0, 0, \
					0.5*model.t[i,0]*N[i,0]*model.v1[0,i], \
					-0.5*model.t[i,0]*N[i,0]*model.v2[0,i]])
		
		DuDxi[3, model.ndof*i:model.ndof*(i+1)] = np.array([ \
					0, DNDxi[i,0], 0, \
					0.5*zeta*model.t[i,0]*DNDxi[i,0]*model.v1[1,i], \
					-0.5*zeta*model.t[i,0]*DNDxi[i,0]*model.v2[1,i]])
		DuDxi[4, model.ndof*i:model.ndof*(i+1)] = np.array([ \
					0, DNDeta[i,0], 0, \
					0.5*zeta*model.t[i,0]*DNDeta[i,0]*model.v1[1,i], \
					-0.5*zeta*model.t[i,0]*DNDeta[i,0]*model.v2[1,i]])
		DuDxi[5, model.ndof*i:model.ndof*(i+1)] = np.array([ \
					0, 0, 0, \
					0.5*model.t[i,0]*N[i,0]*model.v1[1,i], \
					-0.5*model.t[i,0]*N[i,0]*model.v2[1,i]])
		
		DuDxi[6, model.ndof*i:model.ndof*(i+1)] = np.array([ \
					0, 0, DNDxi[i,0], \
					0.5*zeta*model.t[i,0]*DNDxi[i,0]*model.v1[2,i], \
					-0.5*zeta*model.t[i,0]*DNDxi[i,0]*model.v2[2,i]])
		DuDxi[7, model.ndof*i:model.ndof*(i+1)] = np.array([ \
					0, 0, DNDeta[i,0], \
					0.5*zeta*model.t[i,0]*DNDeta[i,0]*model.v1[2,i], \
					-0.5*zeta*model.t[i,0]*DNDeta[i,0]*model.v2[2,i]])
		DuDxi[8, model.ndof*i:model.ndof*(i+1)] = np.array([ \
					0, 0, 0, \
					0.5*model.t[i,0]*N[i,0]*model.v1[2,i], \
					-0.5*model.t[i,0]*N[i,0]*model.v2[2,i]])

	# transformation matrix J_9 between global and parent coordinates
	J_9 = np.zeros((model.nsd*model.nsd, model.nsd*model.nsd)) # 9*9
	for i in range(model.nsd):
		J_9[model.nsd*i:model.nsd*(i+1), model.nsd*i:model.nsd*(i+1)] = J_inv
	
	# derivative of shape function matrix N in global coordinates
	DuDx = J_9 @ DuDxi # 9*40

	# transformation matrix H between strain-displacement matrix B and displacement derivative matrix
	H = np.array(([1, 0, 0, 0, 0, 0, 0, 0, 0], \
				  [0, 0, 0, 0, 1, 0, 0, 0, 0], \
				  [0, 0, 0, 0, 0, 0, 0, 0, 1], \
				  [0, 1, 0, 1, 0, 0, 0, 0, 0], \
				  [0, 0, 0, 0, 0, 1, 0, 1, 0], \
				  [0, 0, 1, 0, 0, 0, 1, 0, 0])) # 6*9
	
	# strain-displacement matrix in global coordinates
	B1 = H @ DuDx # 6*40

	# Compute transformation matrix T between global and local coordinates
	temp = np.cross(DXDxi.T, DXDeta.T) # 1*3
	v3_ = temp.T / (temp[0, 0]**2 + temp[0, 1]**2 + temp[0, 2]**2)**0.5 # 3*1
	v1_ = (np.cross(np.array([1,0,0]), v3_.T)).T # 3*1
	v2_ = (np.cross(v3_.T, v1_.T)).T # 3*1
	
	T = np.array(([v1_[0,0]**2, v1_[1,0]**2, v1_[2,0]**2, v1_[0,0]*v1_[1,0], v1_[1,0]*v1_[2,0], v1_[2,0]*v1_[0,0]], \
				[v2_[0,0]**2, v2_[1,0]**2, v2_[2,0]**2, v2_[0,0]*v2_[1,0], v2_[1,0]*v2_[2,0], v2_[2,0]*v2_[0,0]], \
				[v3_[0,0]**2, v3_[1,0]**2, v3_[2,0]**2, v3_[0,0]*v3_[1,0], v3_[1,0]*v3_[2,0], v3_[2,0]*v3_[0,0]], \
				[2*v1_[0,0]*v2_[0,0], 2*v1_[1,0]*v2_[1,0], 2*v1_[2,0]*v2_[2,0], v1_[0,0]*v2_[1,0]+v1_[1,0]*v2_[0,0], \
				v1_[1,0]*v2_[2,0]+v1_[2,0]*v2_[1,0], v1_[2,0]*v2_[0,0]+v1_[0,0]*v2_[2,0]], \
				[2*v2_[0,0]*v3_[0,0], 2*v2_[1,0]*v3_[1,0], 2*v2_[2,0]*v3_[2,0], v2_[0,0]*v3_[1,0]+v2_[1,0]*v3_[0,0], \
				v2_[1,0]*v3_[2,0]+v2_[2,0]*v3_[1,0], v2_[2,0]*v3_[0,0]+v2_[0,0]*v3_[2,0]], \
				[2*v3_[0,0]*v1_[0,0], 2*v3_[1,0]*v1_[1,0], 2*v3_[2,0]*v1_[2,0], v3_[0,0]*v1_[1,0]+v3_[1,0]*v1_[0,0], \
				v3_[1,0]*v1_[2,0]+v3_[2,0]*v1_[1,0], v3_[2,0]*v1_[0,0]+v3_[0,0]*v1_[2,0]])) # 6*6

	# strain-displacement matrix in local coordinates
	B2 = T @ B1 # 6*40
	
	# ignore normal strain in z_ direction
	B = np.vstack((B2[0:2, :], B2[3:6, :])) # 5*9

	return B, detJ

def Areatop(xi, eta, Ct):
	"""
	Calcualte detJ at coordinate (xi, eta, 1)

	Args:
		xi : The first parent coordinate
		eta : The second parent coordinate
		Ct  : The physical coordinates at the top surface

	Returns:
		Jacobian determination
	"""
	# derivative of shape function NI in parent coordinates
	DNDxi = np.zeros((model.nen, 1))
	DNDxi[0, 0] = 0.25*(2*xi+eta)*(1-eta)
	DNDxi[1, 0] = 0.25*(2*xi-eta)*(1-eta)
	DNDxi[2, 0] = 0.25*(2*xi+eta)*(1+eta)
	DNDxi[3, 0] = 0.25*(2*xi-eta)*(1+eta)
	DNDxi[4, 0] = -xi*(1-eta)
	DNDxi[5, 0] = 0.5*(1-eta**2)
	DNDxi[6, 0] = -xi*(1+eta)
	DNDxi[7, 0] = -0.5*(1-eta**2)
	
	DNDeta = np.zeros((model.nen, 1))
	DNDeta[0, 0] = 0.25*(xi+2*eta)*(1-xi)
	DNDeta[1, 0] = 0.25*(-xi+2*eta)*(1+xi)
	DNDeta[2, 0] = 0.25*(xi+2*eta)*(1+xi)
	DNDeta[3, 0] = 0.25*(-xi+2*eta)*(1-xi)
	DNDeta[4, 0] = -0.5*(1-xi**2)
	DNDeta[5, 0] = (1+xi)*(-eta)
	DNDeta[6, 0] = 0.5*(1-xi**2)
	DNDeta[7, 0] = (1-xi)*(-eta)
	
	GN = np.vstack((DNDxi.T, DNDeta.T)) # 2*8
	
	# Compute Jacobian matrix
	J = GN@Ct
	detJ = np.linalg.det(J)
	
	return detJ