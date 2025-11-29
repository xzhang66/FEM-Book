#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Provide methods to create FE model from a json file, to plot the mesh,
to plot deflection and moment Mx distributions along centerline obtained by FE analysis.

Created on Fri Jun 19 18:56:57 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""
import json
import FEData as model
import numpy as np
import matplotlib.pyplot as plt
from utitls import gauss
from PlateElem import NmatPlate, BmatPlate
import tikzplotlib
from Exact import ExactSolution_Ex_6_3


def create_model_json(DataFile):
	"""
	Initialize the FEM model from file DataFile (in json format)
	"""

	with open(DataFile) as f_obj:
		FEData = json.load(f_obj)

	model.Title = FEData['Title']
	model.nsd = FEData['nsd']
	model.ndof = FEData['ndof']
	model.nnp = FEData['nnp']
	model.nel = FEData['nel']
	model.nen = FEData['nen']
	model.nbe = FEData['nbe']
	model.neq = model.ndof * model.nnp

	# initialize K, d and f
	model.f = np.zeros((model.neq, 1))
	model.d = np.zeros((model.neq, 1))
	model.K = np.zeros((model.neq, model.neq))

	# geometric data
	model.h = FEData['h']
	model.lx = FEData['lx']
	model.ly = FEData['ly']
	model.nelx = FEData['nelx']
	model.nely = FEData['nely']
	model.nenx = model.nelx + 1
	model.neny = model.nely + 1
	model.ae = model.lx / (2 * model.nelx)
	model.be = model.ly / (2 * model.nely)
	if model.nelx % 2 != 0:
		print('No. of Elements  {}  is not even, can not get the centerline of x-axis'.format(model.nelx))
	if model.nely % 2 != 0:
		print('No. of Elements  {}  is not even, can not get the centerline of y-axis'.format(model.nely))
	
	# material properties
	model.E = FEData['E']
	model.nu = FEData['nu']
	model.D = model.E * model.h ** 3 / (12.0 * (1 - model.nu ** 2)) * \
				np.array([[1, model.nu, 0],
						[model.nu, 1, 0],
						[0, 0, (1-model.nu)/2]])

	# gauss integration
	model.ngp = FEData['ngp']

	# boundary conditions
	model.flags = np.array(FEData['flags'])
	model.nd = FEData['nd']
	if model.nbe > 0:
		model.n_bc = np.array(FEData['n_bc'])

	# The Essential B.C. is set to zero by default
	try:
		model.e_bc = np.array(FEData['e_bc'])
	except KeyError:
		model.e_bc = np.zeros((model.neq, 1))

	# force conditions
	# The F.C. is set to zero by default
	try:
		model.P = np.array(FEData['P'])
	except KeyError:
		model.P = np.zeros((model.neq, 1))

	try:
		model.b = np.array(FEData['b'])
	except KeyError:
		model.b = np.zeros((model.nen*model.ndof, model.nel))
	
	try:
		model.q = FEData['q']
	except KeyError:
		model.q = 0.0

	# define the mesh
	model.x = np.array(FEData['x'])
	model.y = np.array(FEData['y'])
	model.IEN = np.array(FEData['IEN'], dtype=int)

	# parameter for postprocess
	model.plot_mesh = FEData['plot_mesh']
	model.plot_nod = FEData['plot_nod']
	model.plot_tex = FEData['plot_tex']

	plot_mesh()

	model.ID = np.zeros(model.neq, dtype=int)
	model.LM = np.zeros((model.nen*model.ndof, model.nel), dtype=int)
	setup_ID_LM()


def point_and_trac():
	"""
	Add nodal forces and natural B.C. to the global force vector.
	"""
	# Assemble uniform load
	P_Q = np.zeros((model.neq, 1))
	for e in range(model.nel):
		P_Q[(model.IEN[0, e]-1)*3] += model.q * model.ae * model.be
		P_Q[(model.IEN[0, e]-1)*3+1] += model.q * model.ae * model.be * model.be / 3.0
		P_Q[(model.IEN[0, e]-1)*3+2] += model.q * model.ae * model.be * (-model.ae) / 3.0
		
		P_Q[(model.IEN[1, e]-1)*3] += model.q * model.ae * model.be
		P_Q[(model.IEN[1, e]-1)*3+1] += model.q * model.ae * model.be * (-model.be) / 3.0
		P_Q[(model.IEN[1, e]-1)*3+2] += model.q * model.ae * model.be * (-model.ae) / 3.0
		
		P_Q[(model.IEN[2, e]-1)*3] += model.q * model.ae * model.be
		P_Q[(model.IEN[2, e]-1)*3+1] += model.q * model.ae * model.be * model.be / 3.0
		P_Q[(model.IEN[2, e]-1)*3+2] += model.q * model.ae * model.be * model.ae / 3.0
		
		P_Q[(model.IEN[3, e]-1)*3] += model.q * model.ae * model.be
		P_Q[(model.IEN[3, e]-1)*3+1] += model.q * model.ae * model.be * (-model.be) / 3.0
		P_Q[(model.IEN[3, e]-1)*3+2] += model.q * model.ae * model.be * model.ae / 3.0
	
	model.f[model.ID - 1] = model.f[model.ID - 1] + P_Q

	# Assemble point forces
	model.f[model.ID - 1] = model.f[model.ID - 1] + model.P

	# Compute nodal boundary force vector
	for e in range(model.nbe):
		ft = np.zeros((4, 1))							# initialize nodal boundary force vector
		node1 = int(model.n_bc[0, e])					# first node
		node2 = int(model.n_bc[1, e])					# second node
		n_bce = model.n_bc[2:, e].reshape((-1, 1))		# traction value at node1

		# coordinates
		x1 = model.x[node1 - 1]
		y1 = model.y[node1 - 1]
		x2 = model.x[node2 - 1]
		y2 = model.y[node2 - 1]

		# edge length
		leng = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
		J = leng/2.0

		w, gp = gauss(model.ngp)

		for i in range(model.ngp):
			psi = gp[i]
			N = 0.5*np.array([[1-psi, 0, 1+psi, 0],
							  [0, 1-psi, 0, 1+psi]])

			traction = N@n_bce
			ft = ft + w[i]*J*(N.T@traction)

		# Assemble nodal boundary force vector
		ind1 = model.ndof*(node1 - 1)
		ind2 = model.ndof*(node2 - 1)

		model.f[model.ID[ind1] - 1, 0] += ft[0]
		model.f[model.ID[ind1 + 1] - 1, 0] += ft[1]
		model.f[model.ID[ind2] - 1, 0] += ft[2]
		model.f[model.ID[ind2 + 1] - 1, 0] += ft[3]


def setup_ID_LM():
	"""
	Calculate the ID and LM matrix according to model.flags and model.IEN matrix.
	"""
	count = 0
	count1 = 0
	for i in range(model.neq):
		if model.flags[i] == 2:
			# check if a node on essential boundary
			count += 1
			model.ID[i] = count
			model.d[count-1] = model.e_bc[i]
		else:
			count1 += 1
			model.ID[i] = model.nd + count1

	for i in range(model.nel):
		n = 0
		for j in range(model.nen):
			blk = model.ndof * (model.IEN[j, i] - 1)
			for k in range(model.ndof):
				model.LM[n, i] = model.ID[blk + k]
				n += 1

def plot_mesh():
	"""
	Plot the initial mesh and print the mesh parameters.
	"""
	if model.plot_mesh == 'yes':
		for i in range(model.nbe):
			# plot Natural B.C. in red lines
			node1 = model.n_bc[0, i]
			node2 = model.n_bc[1, i]

			# coordinates
			x1 = model.x[node1 - 1]
			y1 = model.y[node1 - 1]
			x2 = model.x[node2 - 1]
			y2 = model.y[node2 - 1]

			plt.plot([x1, x2], [y1, y2], color='r', linewidth=4)

		for i in range(model.nel):
			XX = [model.x[model.IEN[0, i] - 1], model.x[model.IEN[1, i] - 1], model.x[model.IEN[2, i] - 1],
				  model.x[model.IEN[3, i] - 1], model.x[model.IEN[0, i] - 1]]
			YY = [model.y[model.IEN[0, i] - 1], model.y[model.IEN[1, i] - 1], model.y[model.IEN[2, i] - 1],
				  model.y[model.IEN[3, i] - 1], model.y[model.IEN[0, i] - 1]]
			plt.plot(XX, YY, color='b')

			if model.plot_nod == 'yes':
				plt.text(XX[0], YY[0], str(model.IEN[0, i]))
				plt.text(XX[1], YY[1], str(model.IEN[1, i]))
				plt.text(XX[2], YY[2], str(model.IEN[2, i]))
				plt.text(XX[3], YY[3], str(model.IEN[3, i]))

		plt.title('Meshed plate')
		plt.xlabel(r'$X$')
		plt.ylabel(r'$Y$')

	print('  Mesh Params ')
	print('No. of Elements  {}'.format(model.nel))
	print('No. of Nodes     {}'.format(model.nnp))
	print('No. of Equations {}'.format(model.neq))


def postprocess():
	"""
	Plot deflection and moment Mx distributions along centerline obtained by FE analysis.
	
	"""
	if model.plot_mesh:
		# Convert matplotlib figures into PGFPlots figures
		if model.plot_tex == "yes":
			tikzplotlib.save("plate-mesh.pgf")
            
#		plt.savefig("plate-mesh.pdf")
		plt.show()

	# plot deflection and moment Mx distributions along centerline
	fig, (ax1, ax2) = plt.subplots(2,1)
	plt.tight_layout()

	ax1.set_title('FE analysis of centerline')
	ax1.set_ylabel('deflection')

	ax2.set_xlabel('x')
	ax2.set_ylabel('moment Mx')

	for e in range(model.nelx):
		n_e = int(model.nelx * (model.nely / 2 - 1)) + e
		centerline_deflection_Mx(n_e, ax1, ax2)

	# plot the exact deflection and moment Mx distributions along centerline
	ExactSolution_Ex_6_3(ax1, ax2)

	ax1.legend()
	ax2.legend()

	# Convert matplotlib figures into PGFPlots figures
	if model.plot_tex == "yes":
		tikzplotlib.save("plate-centerline.pgf")

#	plt.savefig("plate-centerline.pdf")
	plt.show()

def centerline_deflection_Mx(e, ax1, ax2):
	"""
	Plot deflection and moment Mx distributions along the psi = 1 line of an element
	
	"""
	# get coordinate and deflection of element nodes
	je = model.IEN[:, e] - 1
	C = np.array([model.x[je], model.y[je]]).T
	de = model.d[model.LM[:,e]-1]
	
	# equally distributed coordinates on the psi = 1 line of an element
	xplot = np.linspace(C[0,0], C[1,0], model.nplot)
	etaplot = (2*xplot - C[0,0] - C[1,0])/(C[1,0] - C[0,0])
	psiplot = 1.0

	deflection = np.zeros(model.nplot)
	moment_x = np.zeros(model.nplot)
	moment_all = np.zeros(3)
	
	for i in range(model.nplot):
		eta = etaplot[i]
		N = NmatPlate(eta, psiplot)
		B, detJ = BmatPlate(eta, psiplot)
		deflection[i] = N@de
		moment_all = -model.D@B@de
		moment_x[i] = moment_all[0]
	
	c0 = np.zeros(model.nplot)
	
	# plot deflection and moment Mx
	line1, = ax1.plot(xplot, c0, 'k')
	line2, = ax1.plot(xplot, deflection, 'b')
	line3, = ax2.plot(xplot, moment_x, 'b')
	if e - int(model.nelx * (model.nely / 2 - 1)) == 0:
		line2.set_label('FE')
		line3.set_label('FE')