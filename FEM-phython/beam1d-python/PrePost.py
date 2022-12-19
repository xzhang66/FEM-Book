#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provide methods to create FE model for a beam from a json file, to plot the
beam, to print stresses at Gauss points, to plot displacement and stress
distributions obtained by FE analysis and exact solution.

Created on Aug. 11 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""

import numpy as np
import json
import matplotlib.pyplot as plt
import tikzplotlib

from utils import gauss
from Beam1DElem import Nmatrix1D, Bmatrix1D, Smatrix1D
from Exact import ExactSolution_Fish_10_1, ExactSolution_Ex_6_1
import FEData as model


def create_model_json(DataFile):
	"""
	Initialize the FEM model from file DataFile (in json format)
	"""

	# input data from json file
	with open(DataFile) as f_obj:
		FEData = json.load(f_obj)

	model.Title = FEData['Title']
	model.nsd = FEData['nsd']
	model.ndof = FEData['ndof']
	model.nnp = FEData['nnp']
	model.nel = FEData['nel']
	model.nen = FEData['nen']
	model.neq = model.ndof*model.nnp
	model.neqe = model.ndof*model.nen

	# initialize K, d and f
	model.f = np.zeros((model.neq, 1))
	model.d = np.zeros((model.neq, 1))
	model.K = np.zeros((model.neq, model.neq))

	# element and material data (given at the element nodes)
	model.E = np.array(FEData['E'])
	model.body = np.array(FEData['body'])
	model.CArea = np.array(FEData['CArea'])

	# gauss integration
	model.ngp = FEData['ngp']

	# boundary conditions
	model.flags = np.array(FEData['flags'])
	model.e_bc = np.array(FEData['e_bc'])
	model.n_bc = np.array(FEData['n_bc'])
	model.nd = FEData['nd']

	# point forces
	model.np = FEData['np']
	if model.np > 0:
		model.xp = np.array(FEData['xp'])
		model.P = np.array(FEData['P'])

	# output plots
	model.plot_beam = FEData['plot_beam']
	model.plot_nod = FEData['plot_nod']
	model.nplot = model.nen * 10
	model.plot_tex = FEData['plot_tex']
	model.Exact = FEData['Exact']

	# define the mesh
	model.x = np.array(FEData['x'])
	model.y = np.array(FEData['y'])
	model.IEN = np.array(FEData['IEN'], np.int)
	model.leng = model.x[model.IEN[1,:]-1] - \
                 model.x[model.IEN[0,:]-1]

	model.ID = np.zeros(model.neq, np.int)
	model.LM = np.zeros((model.neqe, model.nel), np.int)

	# generate LM and ID arrays
	setup_ID_LM()


def setup_ID_LM():
	""" Setup ID and LM arrays """
	count = 0
	count1 = 0

	# Reorder the D.O.F. to make the essential B.C. numbered first
	for i in range(model.neq):
		if model.flags[i] == 2:		# Essential boundary node
			count += 1
			model.ID[i] = count		# The reordered number of essential B.C
			model.d[count] = model.e_bc[i]
		else:
			count1 += 1
			model.ID[i] = model.nd + count1

	for i in range(model.nel):
		for j in range(model.nen):
			for k in range(model.ndof):
				ind = j*model.ndof + k
				model.LM[ind, i] = model.ID[model.ndof*(model.IEN[j,i] - 1) + k]


def naturalBC():
	""" Compute and assemble nodal boundary force vector """
	for i in range(model.neq):
		if model.flags[i] == 1:
			dof = model.ID[i] - 1
			model.f[dof] += model.n_bc[dof]


def plotbeam():
	""" Plot the beam """
	if model.plot_beam == 'yes':
		for i in range(model.nel):
			XX = np.array([model.x[model.IEN[0,i]-1], model.x[model.IEN[1,i]-1]])
			YY = np.array([model.y[model.IEN[0,i]-1], model.y[model.IEN[1,i]-1]])
			plt.figure(1)
			plt.plot(XX, YY)
			plt.plot(XX, -YY)
			plt.plot(XX, [0,0], '+r')

			if model.plot_nod == 'yes':
				plt.text(XX[0], 0, str(model.IEN[0, i]))
				plt.text(XX[1], 0, str(model.IEN[1, i]))

		plt.plot([ model.x[model.IEN[0,0] - 1], model.x[model.IEN[0,0] - 1]],
				 [-model.y[model.IEN[0,0] - 1], model.y[model.IEN[0,0] - 1]])
		plt.plot([ model.x[model.IEN[-1,-1] - 1], model.x[model.IEN[-1,-1] - 1]],
				 [-model.y[model.IEN[-1,-1] - 1], model.y[model.IEN[-1,-1] - 1]])
		plt.title('Beam Plot')
		plt.axis('equal')
		plt.show()

		# print some mesh parameters
		print('\n  Beam Params')
		print('No. of Elements  ' + str(model.nel))
		print('No. of Nodes     ' + str(model.nnp))
		print('No. of Equations ' + str(model.neq))

def disp_moment_and_shear(e, ax1, ax2, ax3):
	"""
	Print the moments and shear forces at the Gauss points, plot displacements,
	moments and shear forces distributions obtained by FE analysis.

	Args:
		e: (int) element number
		ax1 : axis to draw displacement distribution
		ax2 : axis to draw moment distribution
		ax3 : axis to draw shear force distribution
	"""
	de = model.d[model.LM[:,e]-1]
	IENe = model.IEN[:,e]-1
	xe = model.x[IENe]
	J = (xe[-1] - xe[0])/2
	w, gp = gauss(model.ngp)

	gauss_pt = np.zeros(model.ngp)
	moment_gauss = np.zeros(model.ngp)
	shear_gauss = np.zeros(model.ngp)

	for i in range(model.ngp):
		gauss_pt[i] = 0.5*(xe[0]+xe[-1])+J*gp[i]
		N = Nmatrix1D(gp[i], xe)
		B = Bmatrix1D(gp[i], xe)*1/J**2
		S = Smatrix1D(gp[i], xe)*1/J**3
		Ee = model.E[e]

		moment_gauss[i] = Ee*B@de
		shear_gauss[i] = Ee*S@de

	print("%8d %12.6f %12.6f %16.6f %16.6f %16.6f %16.6f"%
		  (e, gauss_pt[0], gauss_pt[1], moment_gauss[0], moment_gauss[1], shear_gauss[0], shear_gauss[1]))

	# equally distributed coordinate within an element
	xplot = np.linspace(xe[0], xe[-1], model.nplot)
	xplotgauss = (2*xplot - xe[0] - xe[-1])/(xe[-1] - xe[0])

	displacement = np.zeros(model.nplot)
	moment = np.zeros(model.nplot)
	shear = np.zeros(model.nplot)

	for i in range(model.nplot):
		xi = xplotgauss[i]
		N = Nmatrix1D(xi, xe)
		B = Bmatrix1D(xi, xe)*1/J**2
		S = Smatrix1D(xi, xe)*1/J**3
		Ee = model.E[e]
		displacement[i] = N@de
		moment[i] = Ee*B@de
		shear[i] = Ee*S@de

	# plot displacements and moments and shear forces
	line1, = ax1.plot(xplot, displacement)
	line2, = ax2.plot(xplot, moment)
	line3, = ax3.plot(xplot, shear)
	if e == 0:
		line1.set_label('FE')
		line2.set_label('FE')
		line3.set_label('FE')

def postprocessor():
	"""
	Loop over elements to print&plot displacements, moments and shear forces

	Args：
		BeamType: the beam type of exact solution
	"""

	print()
	print('Print stresses at the Gauss points \n')
	print('%8s %12s %12s %16s %16s %16s %16s'
		  % ("Element", "x(gauss1)", "x(gauss2)", "moment(gauss1)", "moment(gauss2)",
			 "shear force(gauss1)", "shear force(gauss2)"))
	print('-------------------------------------------------------------------------------------------------------')

	fig, (ax1, ax2, ax3) = plt.subplots(3,1)
	plt.tight_layout()

	ax1.set_title('FE analysis of 1D beam')
	ax1.set_ylabel('displacement')

	ax2.set_ylabel('moment')

	ax3.set_xlabel('x')
	ax3.set_ylabel('shear force')

	for e in range(model.nel):
		disp_moment_and_shear(e, ax1, ax2, ax3)

	if model.Exact == "Fish-10.1":
		ExactSolution_Fish_10_1(ax1, ax2, ax3)
	elif model.Exact == "Ex-6-1":
		ExactSolution_Ex_6_1(ax1, ax2, ax3)
	else:
		print('Exact solution for %s is not available'%(model.Exact))

	ax1.legend()
	ax2.legend()
	ax3.legend()
	plt.show()

	# Convert matplotlib figures into PGFPlots figures stored in a Tikz file,
	# which can be added into your LaTex source code by "\input{fe_plot.tex}"
	if model.plot_tex == "yes":
		tikzplotlib.save("fe_plot.tex")