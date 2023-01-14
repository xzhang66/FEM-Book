#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot exact solution for displacement along the computational domain

Created on Thu 5 16:56:57 2023

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""

import FDData as model
import numpy as np
import matplotlib.pyplot as plt

def ExactSolution():
	"""
	Plot the exact displacement along the computational domain.

	"""
	np_exact = 201
	dx = (model.x_up - model.x_low) / (np_exact - 1)
	x_exact = np.zeros((np_exact, 1))
	u_exact = np.zeros((np_exact, 1))
	
	for i in range(np_exact):
		x_exact[i] = model.x_low + i * dx

	for i in range(np_exact):
		if x_exact[i] < (0.5 + model.t_end):
			u_exact[i] = 1
		else:
			u_exact[i] = 0

	line2, = plt.plot(x_exact, u_exact, '--r', label='Exact')