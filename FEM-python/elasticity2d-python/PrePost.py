#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Provide methods to create FE model from a json file, to plot the mesh,
to print stresses at nodes, to plot displacement and stress distributions
obtained by FE analysis.

Created on Fri Jun 19 18:56:57 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""
import json
import FEData as model
import numpy as np
import matplotlib.pyplot as plt
from utitls import gauss
from Elast2DElem import NmatElast2D, BmatElast2D
import matplotlib.colors
import matplotlib.cm
import tikzplotlib


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

	# material properties
	E = FEData['E']
	ne = FEData['nu']

	model.plane_strain = FEData['plane_strain']
	if model.plane_strain == 1:   # Plane strain
		E = E / (1.0 - ne * ne)
		ne = ne / (1.0 - ne)

	model.D = np.array([[1, ne, 0],
						[ne, 1, 0],
						[0, 0, (1-ne)/2]])*E/(1 - ne**2)
	model.G = E / (2.0 * (1.0 + ne))

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

	# define the mesh
	model.x = np.array(FEData['x'])
	model.y = np.array(FEData['y'])
	model.IEN = np.array(FEData['IEN'], dtype=np.int)

	# parameter for postprocess
	model.counter = np.zeros((model.nnp, 1))
	model.nodestress = np.zeros((model.nnp, 3))

	model.plot_mesh = FEData['plot_mesh']
	model.plot_nod = FEData['plot_nod']
	model.plot_disp = FEData['plot_disp']
	model.print_disp = FEData['print_disp']
	model.compute_stress = FEData['compute_stress']
	model.plot_stress_xx = FEData['plot_stress_xx']
	model.plot_mises = FEData['plot_mises']
	model.plot_tex = FEData['plot_tex']
	model.fact = np.double(FEData['fact'])

	plot_mesh()

	model.ID = np.zeros(model.neq, dtype=np.int)
	model.LM = np.zeros((model.nen*model.ndof, model.nel), dtype=np.int)
	setup_ID_LM()


def point_and_trac():
	"""
	Add the nodal forces and natural B.C. to the global force vector.
	"""
	# Assemble point forces
	model.f[model.ID - 1] = model.f[model.ID - 1] + model.P

	# Compute nodal boundary force vector
	for i in range(model.nbe):
		ft = np.zeros((4, 1))							# initialize nodal boundary force vector
		node1 = int(model.n_bc[0, i])					# first node
		node2 = int(model.n_bc[1, i])					# second node
		n_bce = model.n_bc[2:, i].reshape((-1, 1))		# traction value at node1

		# coordinates
		x1 = model.x[node1 - 1]
		y1 = model.y[node1 - 1]
		x2 = model.x[node2 - 1]
		y2 = model.y[node2 - 1]

		# edge length
		leng = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
		J = leng/2.0

		ngp = abs(model.ngp)
		w, gp = gauss(ngp)

		for j in range(ngp):
			psi = gp[j]
			N = 0.5*np.array([[1-psi, 0, 1+psi, 0],
							  [0, 1-psi, 0, 1+psi]])

			traction = N@n_bce
			ft = ft + w[j]*J*(N.T@traction)

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
			node1 = int(model.n_bc[0, i])
			node2 = int(model.n_bc[1, i])

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

		plt.title('Initial structure')
		plt.xlabel(r'$X$')
		plt.ylabel(r'$Y$')

	print('  Mesh Params ')
	print('No. of Elements  {}'.format(model.nel))
	print('No. of Nodes     {}'.format(model.nnp))
	print('No. of Equations {}'.format(model.neq))


def postprocess():
	"""
	1. Calculate the coordinates of deformed configuration.
	2. Plot the initial and deformed configuration in one figure.
	3. Print the element stress on Gauss Point.
	4. Calculate the nodal stress and plot the stress contours.
	"""
	# print the nodal displacement
	print_displacement()

	# plot the deformed configuration
	displacement()

	if model.plot_mesh or model.plot_disp:
        	# Convert matplotlib figures into PGFPlots figures
		if model.plot_tex == "yes":
			tikzplotlib.save("elasticity-mesh.pgf")
            
		plt.savefig("elasticity-mesh.pdf")
		plt.show()

	# Compute strains and stresses at gauss points
	if model.compute_stress == 'yes':
		print("\n                     Stress at Gauss Points")
		print("-----------------------------------------------------------------------------")
		for e in range(model.nel):
			print("Element  {}".format(e))
			print("-------------")
			get_stress(e)
			nodal_stress(e)

		stress_contours()


def displacement():
	"""
	Calculate the coordinates of deformed configuration.

	Plot the deformed configuration.
	"""
	if model.plot_disp == 'yes':
		dis = model.d[model.ID - 1]
		dis *= model.fact

		#  Deformed coordinates
		j = 0
		xnew = np.zeros(model.nnp)
		ynew = np.zeros(model.nnp)
		for i in range(0, model.nnp*model.ndof, model.ndof):
			xnew[j] = model.x[j] + dis[i]
			ynew[j] = model.y[j] + dis[i+1]
			j += 1

		xnew = xnew.T.squeeze()
		ynew = ynew.T.squeeze()
		# plot deformed shape over the initial configuration
		for i in range(model.nel):
			XXnew = [xnew[model.IEN[0, i] - 1], xnew[model.IEN[1, i] - 1], xnew[model.IEN[2, i] - 1],
				  xnew[model.IEN[3, i] - 1], xnew[model.IEN[0, i] - 1]]
			YYnew = [ynew[model.IEN[0, i] - 1], ynew[model.IEN[1, i] - 1], ynew[model.IEN[2, i] - 1],
				  ynew[model.IEN[3, i] - 1], ynew[model.IEN[0, i] - 1]]
			plt.plot(XXnew, YYnew, color='k')

		plt.title('Initial and deformed structure')


def print_displacement():
	"""
	Print the displacement of all nodes.
	"""
	if model.print_disp == 'yes':
		dis = model.d[model.ID - 1]

		print("\n                           Nodal displacement")
		print("-------------------------------------------------------------------------------")
		print("\tnode\tx\ty\t\t\tu_x\t\t\tu_y")

		for i in range(model.nnp):
			print("\t{}\t{}\t{}\t{:.15e}\t{:.15e}".format(i+1, model.x[i], model.y[i], dis[2*i], dis[2*i+1]))


def get_stress(e):
	"""
	Print the element stress on Gauss Point.

	Args:
		e   : The element number
	"""
	de = model.d[model.LM[:, e] - 1]		# extract element nodal displacements

	# get coordinates of element nodes
	je = model.IEN[:, e] - 1
	C = np.array([model.x[je], model.y[je]]).T

	ngp = abs(model.ngp)
	# get gauss points and weights
	w, gp = gauss(ngp)

	# compute strains and stress at the gauss points
	ind = 0
	number_gp = ngp*ngp
	X = np.zeros((number_gp, 2))
	strain = np.zeros((3, number_gp))
	stress = np.zeros((3, number_gp))
	for i in range(ngp):
		for j in range(ngp):
			eta = gp[i]
			psi = gp[j]
			# shape functions matrix
			N = NmatElast2D(eta, psi)
			# derivative of the shape functions
			B, detJ = BmatElast2D(eta, psi, C)

			Na = np.array([N[0, 0], N[0, 2], N[0, 4], N[0, 6]])
			X[ind, :] = Na@C
			strain[:, ind] = (B@de).T.squeeze()
			stress[:, ind] = (model.D@(strain[:, ind].reshape((-1, 1)))).T.squeeze()

			ind += 1

	print("\tx-coord\t\t\ty-coord\t\t\ts_xx\t\t\ts_yy\t\t\ts_xy")
	for i in range(number_gp):
		print("\t{}\t{}\t{}\t{}\t{}".format(X[i, 0], X[i, 1], stress[0, i], stress[1, i], stress[2, i]))


def nodal_stress(e):
	"""
	Calculate the nodal stress

	Args:
		e   : The element number
	"""
	de = model.d[model.LM[:, e] - 1]  # extract element nodal displacements

	# get coordinates of element nodes
	je = model.IEN[:, e] - 1
	C = np.array([model.x[je], model.y[je]]).T

	# parent coordinates at nodes
	psi_val = np.array([-1, 1, 1, -1])
	eta_val = np.array([-1, -1, 1, 1])

	# compute strains and stresses at the element nodes
	ind = 0
	strain = np.zeros((3, 4))
	stress = np.zeros((3, 4))
	for i in range(model.nen):
		eta = eta_val[i]
		psi = psi_val[i]

		B, detJ = BmatElast2D(eta, psi, C)

		strain[:, ind] = (B @ de).T.squeeze()
		stress[:, ind] = (model.D @ (strain[:, ind].reshape((-1, 1)))).T.squeeze()

		ind += 1

	model.counter[je] += np.ones((model.nen, 1))
	model.nodestress[je, :] += stress.T


def stress_contours():
	"""
	plot the stress contours
	"""
	if model.plot_stress_xx == 'yes':
		# sigma_xx contour
		# The vmin and vmax of colorbar is different in different problems
		# Please determine it by yourself for different cases.
		vmin = -150
		vmax = 200
		norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
		for i in range(model.nel):
			XX = [[model.x[model.IEN[0, i] - 1], model.x[model.IEN[1, i] - 1]],
				  [model.x[model.IEN[3, i] - 1], model.x[model.IEN[2, i] - 1]]]
			YY = [[model.y[model.IEN[0, i] - 1], model.y[model.IEN[1, i] - 1]],
				  [model.y[model.IEN[3, i] - 1], model.y[model.IEN[2, i] - 1]]]

			sxx = model.nodestress[model.IEN[:, i] - 1, 0] / \
                 (model.counter[model.IEN[:, i] - 1].T.squeeze())

			dd = [[sxx[0], sxx[1]], [sxx[3], sxx[2]]]

			plt.contourf(XX, YY, dd, cmap='jet', norm=norm)

		# outline
		for i in range(model.nel):
			XX = [model.x[model.IEN[0, i] - 1], model.x[model.IEN[1, i] - 1], model.x[model.IEN[2, i] - 1],
				  model.x[model.IEN[3, i] - 1], model.x[model.IEN[0, i] - 1]]
			YY = [model.y[model.IEN[0, i] - 1], model.y[model.IEN[1, i] - 1], model.y[model.IEN[2, i] - 1],
				  model.y[model.IEN[3, i] - 1], model.y[model.IEN[0, i] - 1]]
			plt.plot(XX, YY, color='k')

		plt.title(r'$\sigma_{xx}$ contours')
		plt.xlabel(r'$X$')
		plt.ylabel(r'$Y$')
		plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap='jet'))

        	# Convert matplotlib figures into PGFPlots figures
		if model.plot_tex == "yes":
			tikzplotlib.save("elasticity-sxx.pgf")

		plt.savefig("elasticity-sxx.pdf")
		plt.show()

	if model.plot_mises == 'yes':
		# Von Mises stress contour
		# The vmin and vmax of colorbar is different in different problems
		# Please determine it by yourself for different cases.
		vmin = 0
		vmax = 250
		norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
		for i in range(model.nel):
			XX = [[model.x[model.IEN[0, i] - 1], model.x[model.IEN[1, i] - 1]],
				  [model.x[model.IEN[3, i] - 1], model.x[model.IEN[2, i] - 1]]]
			YY = [[model.y[model.IEN[0, i] - 1], model.y[model.IEN[1, i] - 1]],
				  [model.y[model.IEN[3, i] - 1], model.y[model.IEN[2, i] - 1]]]

			sxx = model.nodestress[model.IEN[:, i] - 1, 0] / (model.counter[model.IEN[:, i] - 1].T.squeeze())
			syy = model.nodestress[model.IEN[:, i] - 1, 1] / (model.counter[model.IEN[:, i] - 1].T.squeeze())
			sxy = model.nodestress[model.IEN[:, i] - 1, 2] / (model.counter[model.IEN[:, i] - 1].T.squeeze())

			S1 = 0.5*(sxx + syy) + np.sqrt((0.5*(sxx - syy))**2 + sxy**2)
			S2 = 0.5*(sxx + syy) - np.sqrt((0.5*(sxx - syy))**2 + sxy**2)
			mises = np.sqrt(S1**2 + S2**2 - S1*S2)

			dd = [[mises[0], mises[1]], [mises[3], mises[2]]]

			plt.contourf(XX, YY, dd, cmap='jet', norm=norm)

		# outline
		for i in range(model.nel):
			XX = [model.x[model.IEN[0, i] - 1], model.x[model.IEN[1, i] - 1], model.x[model.IEN[2, i] - 1],
				  model.x[model.IEN[3, i] - 1], model.x[model.IEN[0, i] - 1]]
			YY = [model.y[model.IEN[0, i] - 1], model.y[model.IEN[1, i] - 1], model.y[model.IEN[2, i] - 1],
				  model.y[model.IEN[3, i] - 1], model.y[model.IEN[0, i] - 1]]
			plt.plot(XX, YY, color='k')

		plt.title(r'Von Mises $\sigma$ contours')
		plt.xlabel(r'$X$')
		plt.ylabel(r'$Y$')
		plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap='jet'))
        
		# Convert matplotlib figures into PGFPlots figures
		if model.plot_tex == "yes":
			tikzplotlib.save("elasticity-mises.pgf")

		plt.savefig("elasticity-mises.pdf")
		plt.show()
            