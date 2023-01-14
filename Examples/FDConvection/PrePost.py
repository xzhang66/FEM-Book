#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Provide methods to create FD model from a json file, to plot the solution.

Created on Thu 5 16:56:57 2023

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""

import json
import FDData as model
import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib
from Exact import ExactSolution


def create_model_json(DataFile):
	"""
	Initialize the FEM model from file DataFile (in json format)
	"""

	with open(DataFile) as f_obj:
		FEData = json.load(f_obj)

	model.Title = FEData['Title']
	model.x_low = FEData['x_low']
	model.x_up = FEData['x_up']
	model.dx = FEData['dx']
	model.ratio = FEData['ratio']
	model.dt = model.dx * model.ratio
	model.t_end = FEData['t_end']
	
	model.nx = round((model.x_up - model.x_low) / model.dx) + 1
	model.nt = round((model.t_end - 0) / model.dt)
	model.x = np.zeros((model.nx, 1))
	for i in range(model.nx):
		model.x[i] = model.x_low + i * model.dx
	model.u_now = np.zeros((model.nx, 1))
	model.u_next = np.zeros((model.nx, 1))
	
	model.plot_curve = FEData['plot_curve']
	model.plot_tex = FEData['plot_tex']
	model.plot_region = FEData['plot_region']


def postprocess():
	"""
	Plot deflection distributions along centerline obtained by FE analysis.
	
	"""
	if model.plot_curve == 'yes':
		# plot displacement along the computational domain
		plt.title('FD analysis of the problem')
		plt.ylabel('u')
		plt.xlabel('x')

		line1, = plt.plot(model.x, model.u_now, 'o')
		line1.set_label('FD')

		# plot the exact displacement along the computational domain
		ExactSolution()

		plt.legend()
		plt.xlim(model.plot_region[0],model.plot_region[1])

		# Convert matplotlib figures into PGFPlots figures
		if model.plot_tex == "yes":
			tikzplotlib.save("convection-curve.pgf")

		plt.show()
