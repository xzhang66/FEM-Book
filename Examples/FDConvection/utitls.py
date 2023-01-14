#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provides utilities used by FD analysis.
  1. Apply_initial_condition: Apply initial condition.
  2. Apply_boundary_condition: Apply boundary condition.
  3. solve: Solving the discretized governing equations.

Created on Thu 5 16:56:57 2023

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""

import FDData as model


def Apply_initial_condition():
	"""
	Apply initial condition.

	"""
	for i in range(model.nx):
		if model.x[i] < 0.5:
			model.u_now[i] = 1
		else:
			model.u_now[i] = 0


def Apply_boundary_condition():
	"""
	Apply boundary condition.

	"""
	model.u_next[0] = 1


def solve():
	"""
	Solving the discretized governing equations.

	"""
	
	for it in range(model.nt):
		
		# solve the discretized governing equation
		for ix in range(1,model.nx):
			model.u_next[ix] = (1.0 - model.ratio) * model.u_now[ix] + model.ratio * model.u_now[ix-1]
		
		# calculate the values of the lower boundary node according to the boundary condition
		Apply_boundary_condition()
		
		for ix in range(model.nx):
			model.u_now[ix] = model.u_next[ix]
		model.t = model.t + model.dt
