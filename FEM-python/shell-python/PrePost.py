#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Provide methods to create FE model from a json file, to plot the mesh.

Created on Fri Jun 19 18:56:57 2020

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""

import json
import FEData as model
import numpy as np
import matplotlib.pyplot as plt
from utitls import gauss
from ShellElem import NmatShell, Areatop
import tikzplotlib
from Exact import ExactSolution


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
	model.xb = np.array(FEData['xb'])
	model.yb = np.array(FEData['yb'])
	model.zb = np.array(FEData['zb'])
	model.xt = np.array(FEData['xt'])
	model.yt = np.array(FEData['yt'])
	model.zt = np.array(FEData['zt'])
	model.xI = (model.xt + model.xb) / 2.0
	model.yI = (model.yt + model.yb) / 2.0
	model.zI = (model.zt + model.zb) / 2.0
	model.V3 = np.array([model.xt-model.xb, model.yt-model.yb, model.zt-model.zb]) #3*21
	model.t = np.zeros((model.nnp, 1))
	for i in range(model.nnp):
		model.t[i, 0] = (model.V3[0, i]**2 + model.V3[1, i]**2 + model.V3[2, i]**2)**0.5
	
	model.v3 = np.zeros((3, model.nnp))
	for i in range(model.nnp):
		model.v3[:, i] = model.V3[:, i] / model.t[i, 0]
	model.v1 = np.zeros((3, model.nnp))
	model.v2 = np.zeros((3, model.nnp))
	for i in range(model.nnp):
		model.v1[:, i] = (np.cross(np.array([1,0,0]), model.v3[:, i].T)).T
		model.v2[:, i] = (np.cross(model.v3[:, i].T, model.v1[:, i].T)).T
	
	model.lx = FEData['lx']
	model.ly = FEData['ly']
	model.nelx = FEData['nelx']
	model.nely = FEData['nely']
	if model.nelx % 2 != 0:
		print('No. of Elements  {}  is not even, can not get the center deflection'.format(model.nelx))
	if model.nely % 2 != 0:
		print('No. of Elements  {}  is not even, can not get the center deflection'.format(model.nely))

	# material properties
	model.E = FEData['E']
	model.ne = FEData['nu']
	model.G = model.E / (2.0 * (1.0 + model.ne))
	shcof = 5/6.0 #shear correction factor
	model.D = model.E/(1-model.ne**2) * \
					np.array([[1, model.ne, 0, 0, 0],
							[model.ne, 1, 0, 0, 0],
							[0, 0, (1-model.ne)/2, 0, 0],
							[0, 0, 0, shcof*(1-model.ne)/2, 0],
							[0, 0, 0, 0, shcof*(1-model.ne)/2]])

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
	
	# Add the uniform loada
	if model.nbe > 0:
		temp = np.zeros((model.nen*model.ndof+1, model.nel))
		for e in range(model.nel):
			for i in range(model.nen):
				temp[0, e] = e + 1
				temp[1+i*model.ndof:1+(i+1)*model.ndof, e] += np.array([0, 0, model.q, 0, 0])
		model.n_bc = np.hstack((model.n_bc, temp))
	elif model.nbe == 0:
		model.n_bc = np.zeros((model.nen*model.ndof+1, model.nel))
		for e in range(model.nel):
			for i in range(model.nen):
				model.n_bc[0, e] = e + 1
				model.n_bc[1+i*model.ndof:1+(i+1)*model.ndof, e] += np.array([0, 0, model.q, 0, 0])
	
	model.nbe = model.nbe + model.nel

	# define the mesh
	model.IEN = np.array(FEData['IEN'], dtype=int)

	# parameter for postprocess
	model.plot_mesh = FEData['plot_mesh']
	model.plot_nod = FEData['plot_nod']
	model.plot_centerline = FEData['plot_centerline']
	model.plot_tex = FEData['plot_tex']

	model.ID = np.zeros(model.neq, dtype=int)
	model.LM = np.zeros((model.nen*model.ndof, model.nel), dtype=int)
	setup_ID_LM()


def point_and_trac():
	"""
	Add the nodal forces and natural B.C. at the top surface to the global force vector.
	"""
	# Assemble point forces
	model.f[model.ID - 1] = model.f[model.ID - 1] + model.P

	# Compute nodal boundary force vector
	for i in range(model.nbe):
		ft = np.zeros((model.nen*model.ndof, 1))		# initialize nodal surface force vector
		e = int(model.n_bc[0, i]) - 1					# element number
		n_bce = model.n_bc[1:, i].reshape((-1, 1))		# traction value at node1

		# get coordinates of element nodes
		je = model.IEN[:, e] - 1
		Ct = np.array([model.xt[je], model.yt[je]]).T # 8*2
		
		ngp = model.ngp
		w, gp = gauss(ngp)
		
		for j in range(ngp):
			for m in range(ngp):
				xi = gp[j]
				eta = gp[m]
			
				# shape functions matrix at the top
				N = NmatShell(xi, eta, 1)
				# surface area
				detJ = Areatop(xi, eta, Ct)

				traction = N@n_bce
				ft = ft + w[j]*w[m]*detJ*(N.T@traction) # 40*1
		
		for loop in range(model.nen*model.ndof):
			n = model.LM[loop, e]-1
			model.f[n] += ft[loop]   # assemble nodal force vector

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
	Plot the initial mesh on the midsurface and print the mesh parameters.
	"""
	if model.plot_mesh == 'yes':
		for i in range(model.nel):
			XX = [model.xI[model.IEN[0, i] - 1], model.xI[model.IEN[4, i] - 1], model.xI[model.IEN[1, i] - 1],
				  model.xI[model.IEN[5, i] - 1], model.xI[model.IEN[2, i] - 1], model.xI[model.IEN[6, i] - 1],
				  model.xI[model.IEN[3, i] - 1], model.xI[model.IEN[7, i] - 1], model.xI[model.IEN[0, i] - 1]]
			YY = [model.yI[model.IEN[0, i] - 1], model.yI[model.IEN[4, i] - 1], model.yI[model.IEN[1, i] - 1],
				  model.yI[model.IEN[5, i] - 1], model.yI[model.IEN[2, i] - 1], model.yI[model.IEN[6, i] - 1],
				  model.yI[model.IEN[3, i] - 1], model.yI[model.IEN[7, i] - 1], model.yI[model.IEN[0, i] - 1]]
			plt.plot(XX, YY, color='b')

			if model.plot_nod == 'yes':
				plt.text(XX[0], YY[0], str(model.IEN[0, i]))
				plt.text(XX[1], YY[1], str(model.IEN[4, i]))
				plt.text(XX[2], YY[2], str(model.IEN[1, i]))
				plt.text(XX[3], YY[3], str(model.IEN[5, i]))
				plt.text(XX[4], YY[4], str(model.IEN[2, i]))
				plt.text(XX[5], YY[5], str(model.IEN[6, i]))
				plt.text(XX[6], YY[6], str(model.IEN[3, i]))
				plt.text(XX[7], YY[7], str(model.IEN[7, i]))

		plt.title('Meshed midsurface')
		plt.xlabel(r'$X$')
		plt.ylabel(r'$Y$')

	print('  Mesh Params ')
	print('No. of Elements  {}'.format(model.nel))
	print('No. of Nodes     {}'.format(model.nnp))
	print('No. of Equations {}'.format(model.neq))


def postprocess():
	"""
	Plot deflection distributions along centerline obtained by FE analysis.
	
	"""
	if model.plot_mesh == 'yes':
		plot_mesh()
		# Convert matplotlib figures into PGFPlots figures
		if model.plot_tex == "yes":
			tikzplotlib.save("shell-mesh.pgf")
            
#		plt.savefig("shell-mesh.pdf")
		plt.show()

	if model.plot_centerline == 'yes':
		# plot deflection distributions along centerline on the middle plane
		plt.title('FE analysis of centerline on the midsurface')
		plt.ylabel('deflection')
		plt.xlabel('x')

		for e in range(model.nelx):
			n_e = int(model.nelx * (model.nely / 2 - 1)) + e
			centerline_deflection(n_e)

		# plot the exact deflection and moment Mx distributions along centerline
		ExactSolution()

		plt.legend()

		# Convert matplotlib figures into PGFPlots figures
		if model.plot_tex == "yes":
			tikzplotlib.save("shell-centerline.pgf")

	#	plt.savefig("shell-centerline.pdf")
		plt.show()

def centerline_deflection(e):
	"""
	Plot deflection distributions along the eta = 1, zeta = 0 line of an element
	
	"""
	# get coordinate and deflection of element nodes
	je = model.IEN[:, e] - 1
	C = np.array([model.xI[je], model.yI[je], model.zI[je]]).T # 8*3
	de = model.d[model.LM[:,e]-1]
	
	# equally distributed coordinates on the psi = 1 line of an element
	xplot = np.linspace(C[0,0], C[1,0], model.nplot)
	xiplot = (2*xplot - C[0,0] - C[1,0])/(C[1,0] - C[0,0])
	etaplot = 1.0
	zetaplot = 0.0

	deflection = np.zeros(model.nplot)
	deflection_all = np.zeros(3)
	
	for i in range(model.nplot):
		xi = xiplot[i]
		N = NmatShell(xi, etaplot, zetaplot)
		deflection_all = N@de
		deflection[i] = deflection_all[2]
	c0 = np.zeros(model.nplot)
	
	# plot deflection
	line1, = plt.plot(xplot, c0, 'k')
	line2, = plt.plot(xplot, deflection, 'b')
	if e - int(model.nelx * (model.nely / 2 - 1)) == 0:
		line2.set_label('FE')