#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot exact solutions for some problems
	ExactSolution_Cantilever: Plot the exact solution of the cantilever given
	in Example 10.1 in Fish's textbook

Created on Aug. 15 2022

@author: jsli@163.com, xzhang@tsinghua.edu.cn
"""

import numpy as np
import math 
import matplotlib.pyplot as plt
import FEData as model

from Elast2DElem import NmatElast2D, BmatElast2D


def Exact(ax1):
	"""
	Plot the exact deflection and moment Mx along the centerline
	of the plate(Example 6-3) in ax1 and ax2, respectively.

	Args:
		ax1 : axis to draw deflection distribution
		ax2 : axis to draw moment Mx distribution
	"""
	a=1
	b=4
	pa=1
	pb=0
	E=13
	
	dx = 0.1
	nx  = math.ceil((b-a) / dx)
	x = np.arange(a, b, dx)
	sigma_rr = np.zeros(nx, np.float)

	for index, xi in enumerate(x):
		sigma_rr[index] = (pa*a**2-pb*b**2)/(b**2-a**2) - a**2*b**2/(b**2-a**2)/xi**2*(pa-pb)

	xplot = np.arange(a, b, dx)
	line5, = ax1.plot(xplot, sigma_rr, 'r', label='Exact')

def sigmarr():
	"""
	Plot deflection and moment Mx distributions along the radius
	
	"""
	nplot=21
	
	xplot = np.zeros(4*nplot)
	sigma_rr = np.zeros(4*nplot)
	
	xplot_osp = np.zeros(4)
	sigma_osp = np.zeros(4)
	
	e_all = np.array([7, 15, 23, 31])
	
	for index, e in enumerate(e_all):
		
		# get coordinate and deflection of element nodes
		je = model.IEN[:, e] - 1
		C = np.array([model.x[je], model.y[je]]).T
		de = model.d[model.LM[:,e]-1]
		
		# equally distributed coordinates on the psi = 1 line of an element
		
		xplot_e = np.linspace(C[1,0], C[2,0], nplot)
		etaplot = 0.0
		psiplot = (2*xplot_e - C[1,0] - C[2,0])/(C[2,0] - C[1,0])
		
		sigma_rr_e = np.zeros(nplot)
		sigma_all = np.zeros(3)
		n = np.array([0.980785280403230, 0.195090322016128])
		
		for i in range(nplot):
			psi = psiplot[i]
			B, detJ = BmatElast2D(etaplot, psi, C)
			sigma_all = model.D@B@de
			sigma_rr_e[i] = (sigma_all[0]+sigma_all[1])/2.0 + (sigma_all[0]-sigma_all[1])/2.0*n[0] - sigma_all[2]*n[1]
		
		xplot[index*nplot:(index+1)*nplot] = xplot_e[:]
		sigma_rr[index*nplot:(index+1)*nplot] = sigma_rr_e[:]
		
		xplot_osp[index] = xplot_e[10]
		sigma_osp[index] = sigma_rr_e[10]
		
	return xplot, sigma_rr, xplot_osp, sigma_osp