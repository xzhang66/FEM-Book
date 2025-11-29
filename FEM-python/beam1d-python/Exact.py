#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot exact solutions for some problems
	ExactSolution_Cantilever: Plot the exact solution of the cantilever given
	in Example 10.1 in Fish's textbook

Created on Aug. 15 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn, jsli@163.com
"""

import numpy as np
from math import log, sqrt

from utils import gauss
from Beam1DElem import Nmatrix1D, Bmatrix1D, Smatrix1D
import FEData as model

def ExactSolution_Fish_10_1(ax1, ax2, ax3):
	"""
	Plot the exact displacement, moment and shear force of the cantilever 
    beam (Example 10.1, Fish's book) in	ax1, ax2 and ax3, respectively.

	Args:
		ax1 : axis to draw displacement distribution
		ax2 : axis to draw moment distribution
		ax3 : axis to draw shear force distribution
	"""
	L = 12; a1 = 4; p1 = 10; a2 = 8; p2 = -5; a3 = 8; m = -20; p4 = 20
	p3 = 1; E = 1e4; I = 1.0

	x = np.arange(0, 12, 0.05)
	w = np.zeros(240, float)
	M = np.zeros(240, float)
	S = np.zeros(240, float)

	for index, xi in enumerate(x):
		if xi < a1:
			w1 = -p1*xi**2/(6*E*I)*(3*a1-xi)
			M1 = p1 * a1 * (1 - xi / a1)
		else:
			w1 = -p1*a1**2/(6*E*I)*(3*xi-a1)
			M1 = 0

		if xi < a2:
			w2 = -p2*xi**2/(6*E*I)*(3*a2-xi)
			M2 = p2 * a2 * (1 - xi / a2)
		else:
			w2 = -p2*a2**2/(6*E*I)*(3*xi-a2)
			M2 = 0

		w3 = -m*(xi**2)/(2*E*I)
		M3 = m

		w4 = -p4*xi**2*(3*L-xi)/(6*E*I)
		M4 = p4*(L-xi)

		if xi < a3:
			w5 = -p3*xi**2*(6*a3**2-4*xi*a3+xi**2)/(24*E*I)
			M5 = 1/2*p3*(a3-xi)**2
		else:
			w5 = -p3*a3*(0.5*a3)**2/(6*E*I)*(3*xi-a3/2)
			M5 = 0

		if xi <= 12:
			S1 = 20
		else:
			S1 = 0

		if xi <= 8:
			S2 = -5
		else:
			S2 = 0

		if xi < 8:
			S3 = (8-xi)*1
		else:
			S3 = 0

		if xi < 4:
			S4 = 10
		else:
			S4 = 0

		w[index] = w1+w2+w3+w4+w5
		M[index] = -(M1+M2+M3+M4+M5)
		S[index] = S1+S2+S3+S4

	ax1.plot(x, w, '--r', label='Exact')
	ax2.plot(x, M, '--r', label='Exact')
	ax3.plot(x, S, '--r', label='Exact')
    

def ExactSolution_Ex_6_1(ax1, ax2, ax3):
	"""
	Plot the exact displacement, moment and shear force of the cantilever 
    beam (Example 6-1) in ax1, ax2 and ax3, respectively. 

	Args:
		ax1 : axis to draw displacement distribution
		ax2 : axis to draw moment distribution
		ax3 : axis to draw shear force distribution
	"""
    
	E = 1e4; I = 1.0

	x = np.arange(0, 8, 0.01)
	w = np.zeros(800, float)
	M = np.zeros(800, float)
	S = np.zeros(800, float)

	for index, xi in enumerate(x):
		if xi < 4:
			w[index] = (-xi**4/24 + 14*xi**3/3 - 71*xi**2)/(E*I)
			M[index] = -xi**2/2 + 28*xi - 142
			S[index] = 28 - xi
		else:
			w[index] = (-xi**4/24 + 3*xi**3 - 51*xi**2 - 80*xi + 320/3)/(E*I)
			M[index] = -xi**2/2 + 18*xi - 102
			S[index] = 18 - xi

	ax1.plot(x, w, '--r', label='Exact')
	ax2.plot(x, M, '--r', label='Exact')
	ax3.plot(x, S, '--r', label='Exact')
	
def ErrorNorm_Ex_6_1():
	"""
	Calculate and print the error norms (L2, L∞ and L1 norm) for convergence study

	"""
	EI = 10000.0
	x0 = 4
	
	ngp = 3
	[w, gp] = gauss(ngp)  # extract Gauss points and weights for L2 and L1 norm
	
	num_in_node = 20  # resolution for L∞ norm
	
	# L2 norm
	L2Norm_d = 0
	L2Norm_m = 0
	L2Norm_s = 0

	L2NormEx_d = 0
	L2NormEx_m = 0
	L2NormEx_s = 0
	
	# L1 norm 
	L1Norm_d = 0
	L1Norm_m = 0
	L1Norm_s = 0

	L1NormEx_d = 0
	L1NormEx_m = 0
	L1NormEx_s = 0
	
	# L∞ norm
	LinfNorm_d = 0
	LinfNorm_m = 0
	LinfNorm_s = 0

	LinfNorm_d_Norm = 0
	LinfNormEx_m_Norm = 0
	LinfNormEx_s_Norm = 0
	
	# L2 and L1 norm
	for e in range(model.nel):
		de = model.d[model.LM[:, e] - 1]  # extract element nodal displacements
		IENe = model.IEN[:, e] - 1  # extract local connectivity information
		xe = model.x[IENe]  # extract element x coordinates
		J = (xe[model.nen - 1] - xe[0]) / 2  # compute Jacobian
		
		for i in range(ngp):
			xt = 0.5 * (xe[0] + xe[model.nen - 1]) + J * gp[i]  # Gauss points in physical coordinates
			
			N = Nmatrix1D(gp[i], xe)  # shape functions matrix
			Ee = model.EI[e]  # Young's modulus at element gauss points
			
			uh = N @ de  # displacement at gauss point
			# Exact displacement
			if xt < x0:
				uex = 1.0/EI*(-1.0/24.0*xt**4 + 14.0/3.0*xt**3 - 71.0*xt**2)
			else:
				uex = 1.0/EI*(-1.0/24.0*xt**4 + 3.0*xt**3 - 51.0*xt**2 - 80.0*xt + 320.0/3.0)
			L2Norm_d += J * w[i] * (uex - uh) ** 2
			L2NormEx_d += J * w[i] * (uex) ** 2
			
			L1Norm_d += J * w[i] * abs(uex - uh)
			L1NormEx_d += J * w[i] * abs(uex)
			
			B = Bmatrix1D(gp[i], xe)*1/J**2
			
			mh = Ee*B@de  # moments at gauss point
			# Exact moments
			if xt < x0:
				mex = -0.5*xt**2 + 28.0*xt - 142
			else:
				mex = -0.5*xt**2 + 18.0*xt - 102
			L2Norm_m += J * w[i] * (mex - mh) ** 2
			L2NormEx_m += J * w[i] * (mex) ** 2
			
			L1Norm_m += J * w[i] * abs(mex - mh)
			L1NormEx_m += J * w[i] * abs(mex)
			
			S = Smatrix1D(gp[i], xe)*1/J**3
			
			sh = Ee*S@de  # shear forces at gauss point
			# Exact shear forces
			if xt < x0:
				sex = 28.0 - xt
			else:
				sex = 18.0 - xt
			L2Norm_s += J * w[i] * (sex - sh) ** 2
			L2NormEx_s += J * w[i] * (sex) ** 2
			
			L1Norm_s += J * w[i] * abs(sex - sh)
			L1NormEx_s += J * w[i] * abs(sex)

	L2Norm_d = sqrt(L2Norm_d)
	L2NormEx_d = sqrt(L2NormEx_d)

	L2Norm_m = sqrt(L2Norm_m)
	L2NormEx_m = sqrt(L2NormEx_m)

	L2Norm_s = sqrt(L2Norm_s)
	L2NormEx_s = sqrt(L2NormEx_s)


	# L∞ norm
	for e in range(model.nel):
		de = model.d[model.LM[:, e] - 1]  # extract element nodal displacements
		IENe = model.IEN[:, e] - 1  # extract local connectivity information
		xe = model.x[IENe]  # extract element x coordinates
		J = (xe[model.nen - 1] - xe[0]) / 2  # compute Jacobian
		
		for i in range(num_in_node):
			s = 2 / num_in_node * i - 1   # parent coordinate
			xt = 0.5 * (xe[0] + xe[model.nen - 1]) + J * s   # Gauss points in physical coordinates
			
			N = Nmatrix1D(s, xe)  # shape functions matrix
			Ee = model.EI[e]  # Young's modulus at element gauss points
			
			uh = N @ de  # displacement at gauss point
			# Exact displacement
			if xt < x0:
				uex = 1.0/EI*(-1.0/24.0*xt**4 + 14.0/3.0*xt**3 - 71.0*xt**2)
			else:
				uex = 1.0/EI*(-1.0/24.0*xt**4 + 3.0*xt**3 - 51.0*xt**2 - 80.0*xt + 320.0/3.0)
			LinfNorm_d = max(LinfNorm_d, abs(uex - uh))
			LinfNorm_d_Norm = max(LinfNorm_d_Norm, abs((uex - uh)/uex))
			
			B = Bmatrix1D(s, xe)*1/J**2
			
			mh = Ee*B@de  # moments at gauss point
			# Exact moments
			if xt < x0:
				mex = -0.5*xt**2 + 28.0*xt - 142
			else:
				mex = -0.5*xt**2 + 18.0*xt - 102
			LinfNorm_m = max(LinfNorm_m, abs(mex - mh))
			LinfNormEx_m_Norm = max(LinfNormEx_m_Norm, abs((mex - mh)/mex))
			
			S = Smatrix1D(s, xe)*1/J**3
			
			sh = Ee*S@de  # shear forces at gauss point
			# Exact shear forces
			if xt < x0:
				sex = 28.0 - xt
			else:
				sex = 18.0 - xt
			LinfNorm_s = max(LinfNorm_s, abs(sex - sh))
			LinfNormEx_s_Norm = max(LinfNormEx_s_Norm, abs((sex - sh)/sex))
	
	# print error norms
	print('\nError norms')
	print('\nError L2 norms')
	print('%13s %13s %13s %13s %13s %13s %13s'
		% ('h', 'L2Norm_d', 'L2NormRel_d', 'L2Norm_m', 'L2NormRel_m', 'L2Norm_s', 'L2NormRel_s'))
	print('%13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n'
		% (8 / model.nel, L2Norm_d, L2Norm_d / L2NormEx_d, L2Norm_m, L2Norm_m / L2NormEx_m, L2Norm_s, L2Norm_s / L2NormEx_s))
	
	print('\nError L∞ norms')
	print('%13s %13s %13s %13s %13s %13s %13s'
		% ('h', 'LinfNorm_d', 'LinfNormRel_d', 'LinfNorm_m', 'LinfNormRel_m', 'LinfNorm_s', 'LinfNormRel_s'))
	print('%13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n'
		% (8 / model.nel, LinfNorm_d, LinfNorm_d_Norm, LinfNorm_m, LinfNormEx_m_Norm, LinfNorm_s, LinfNormEx_s_Norm))
	
	print('\nError L1 norms')
	print('%13s %13s %13s %13s %13s %13s %13s'
		% ('h', 'L1Norm_d', 'L1NormRel_d', 'L1Norm_m', 'L1NormRel_m', 'L1Norm_s', 'L1NormRel_s'))
	print('%13.6E %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n'
		% (8 / model.nel, L1Norm_d, L1Norm_d / L1NormEx_d, L1Norm_m, L1Norm_m / L1NormEx_m, L1Norm_s, L1Norm_s / L1NormEx_s))

	return 8 / model.nel, L2Norm_d, L2Norm_m, L2Norm_s, LinfNorm_d, LinfNorm_m, LinfNorm_s, L1Norm_d, L1Norm_m, L1Norm_s
