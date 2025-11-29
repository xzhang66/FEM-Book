#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot exact solution for the midsurface of the simply supported shell

Created on Aug. 15 2020

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""

import FEData as model
import numpy as np
import matplotlib.pyplot as plt
import math

def ExactSolution():
	"""
	Plot the exact deflection along the centerline.

	"""
	D = model.E / (12 * (1 - model.ne**2)) * model.t[0]**3
	q = model.q
	a = model.lx
	b = model.ly
	
	dx = 0.1
	nx  = math.ceil(a / dx)
	x = np.arange(-a/2, a/2, dx)
	y = 0.0
	w = np.zeros(nx, float)

	pi = math.pi
	m_all = np.array([1, 3, 5, 7])

	for index, xi in enumerate(x):
		w1 = 0.0
		
		for m_index, m in enumerate(m_all):
			a_m = m*pi*b/(2*a)
			
			# parameters of w1, w2 and w3 in superposition method
			A1 = 4*q*a**4/(pi**5*D) * (-1)**((m-1)/2.0)/m**5
			B1 = (a_m*math.tanh(a_m)+2)/(2*math.cosh(a_m))
			C1 = 1/(2*math.cosh(a_m))
			D1 = m*pi/a
			
			# calulate components of deflection
			w1 += A1 * math.cos(D1*xi) * (1-B1*math.cosh(D1*y)+C1*D1*y*math.sinh(D1*y))
			
		w[index] = w1

	xplot = np.arange(0, a, dx)
	line3, = plt.plot(xplot, w, '--r', label='Exact')