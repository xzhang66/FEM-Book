#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot exact solution for Mindlin plate

Created on Aug. 15 2020

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""

import FEData as model
import numpy as np
import math

def ExactSolution(ax1, ax2):
	"""
	Plot the exact deflection and moment Mx along the centerline
	of the plate in ax1 and ax2, respectively.

	Args:
		ax1 : axis to draw deflection distribution
		ax2 : axis to draw moment Mx distribution
	"""
	D = model.Db[0,0]
	q = model.q
	a = model.lx
	b = model.ly
	
	dx = 0.1
	nx  = math.ceil(a / dx)
	x = np.arange(-a/2, a/2, dx)
	y = 0.0
	w = np.zeros(nx, np.float)
	Mx = np.zeros(nx, np.float)

	pi = math.pi
	K = -4*q*a**2/pi**3
	m_all = np.array([1, 3, 5, 7])
	E = np.zeros(8, np.float)
	E[m_all[0]] = 0.3722*K
	E[m_all[1]] = -0.0380*K
	E[m_all[2]] = -0.0178*K
	E[m_all[3]] = -0.0085*K

	for index, xi in enumerate(x):
		w1 = 0.0
		w2 = 0.0
		w3 = 0.0
		w1_xx = 0.0
		w2_xx = 0.0
		w3_xx = 0.0
		w1_yy = 0.0
		w2_yy = 0.0
		w3_yy = 0.0
		for m_index, m in enumerate(m_all):
			a_m = m*pi*b/(2*a)
			b_m = m*pi*a/(2*b)
			
			# parameters of w1, w2 and w3 in superposition method
			A1 = 4*q*a**4/(pi**5*D) * (-1)**((m-1)/2.0)/m**5
			B1 = (a_m*math.tanh(a_m)+2)/(2*math.cosh(a_m))
			C1 = 1/(2*math.cosh(a_m))
			D1 = m*pi/a
			
			A2 = -a**2/(2*pi**2*D) * E[m]*(-1)**((m-1)/2.0)/(m**2*math.cosh(a_m))
			B2 = a_m*math.tanh(a_m)
			D2 = m*pi/a
			
			A3 = -b**2/(2*pi**2*D) * E[m]*(-1)**((m-1)/2.0)/(m**2*math.cosh(b_m))
			B3 = b_m*math.tanh(b_m)
			D3 = m*pi/b
			
			# calulate components of deflection
			w1 += A1 * math.cos(D1*xi) * (1-B1*math.cosh(D1*y)+C1*D1*y*math.sinh(D1*y))
			w2 += A2 * math.cos(D2*xi) * (D2*y*math.sinh(D2*y)-B2*math.cosh(D2*y))
			w3 += A3 * math.cos(D3*y) * (D3*xi*math.sinh(D3*xi)-B3*math.cosh(D3*xi))
			
			# calulate components of curvatures
			w1_xx += A1 * (-D1**2*math.cos(D1*xi)) * (1-B1*math.cosh(D1*y)+C1*D1*y*math.sinh(D1*y))
			w2_xx += A2 * (-D2**2*math.cos(D2*xi)) * (D2*y*math.sinh(D2*y)-B2*math.cosh(D2*y))
			w3_xx += A3 * math.cos(D3*y) * ((2-B3)*D3**2*math.cosh(D3*xi)+D3**3*xi*math.sinh(D3*xi))
			
			w1_yy += A1 * math.cos(D1*xi) * ((2*C1-B1)*D1**2*math.cosh(D1*y)+C1*D1**3*y*math.sinh(D1*y))
			w2_yy += A2 * math.cos(D2*xi) * ((2-B2)*D2**2*math.cosh(D2*y)+D2**3*y*math.sinh(D2*y))
			w3_yy += A3 * (-D3**2*math.cos(D3*y)) * (D3*xi*math.sinh(D3*xi)-B3*math.cosh(D3*xi))
			
		w[index] = w1 + w2 + w3
		w_xx = w1_xx + w2_xx + w3_xx
		w_yy = w1_yy + w2_yy + w3_yy
		Mx[index] = -D * (w_xx + model.nu * w_yy)

	xplot = np.arange(0, a, dx)
	line4, = ax1.plot(xplot, w, '--r', label='Exact')
	line5, = ax2.plot(xplot, Mx, '--r', label='Exact')