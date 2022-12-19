#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convergence analysis for degenerated shell element.

Plot w_c*D/q/L^4 - L/h cuvers when using full integration and reduced integration.

Created on Thu Apr 30 21:05:47 2020

@author: xzhang
"""

from PrePost import create_model_json
from Shell import FERun
import numpy as np
import FEData as model
import matplotlib.pyplot as plt
import tikzplotlib

# create FE model from DataFile in json format
#DataFile = "shell_4.json"
DataFile = "shell_16.json"
create_model_json(DataFile)

# calculate w_c*D/q/L^4
ratio = np.arange(5, 1000, 5)
nh = len(ratio)
wc_e = np.ones(nh) * 0.004062
wc_f = np.zeros(nh)
wc_r = np.zeros(nh)
for index, ri in enumerate(ratio):
	
	# the plate gets thinner as L/h increases
	h = model.lx / ri
	D = model.E * h**3 /(12 * (1 - model.ne**2))
	model.zt = np.ones(model.nnp) * h
	model.zI = (model.zt + model.zb) / 2.0
	model.V3 = np.array([model.xt-model.xb, model.yt-model.yb, model.zt-model.zb])
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
	
	# full integration
	model.ngp = 3

	# initialize K, d and f
	model.K = np.zeros((model.neq, model.neq))
	model.d = np.zeros((model.neq, 1))
	model.f = np.zeros((model.neq, 1))

	FERun(DataFile)

	wc_f[index] = model.w_c * D / (model.q * model.lx**4)

	# reduced integration
	model.ngp = 2

	# initialize K, d and f
	model.K = np.zeros((model.neq, model.neq))
	model.d = np.zeros((model.neq, 1))
	model.f = np.zeros((model.neq, 1))

	FERun(DataFile)

	wc_r[index] = model.w_c * D / (model.q * model.lx**4)

# plot the curve of w_c*D/q/L^4 vs. L/h
line1, = plt.plot(ratio, wc_r, color='b', label='Reduced integration')
line2, = plt.plot(ratio, wc_f, color='k', label='Full integration')
line3, = plt.plot(ratio, wc_e, color='r', label='Exact thin plate solution')

plt.xlim(5, 1000)
plt.xticks([5,50,100,1000])

plt.xlabel(r'$L/h$')
plt.ylabel(r'$w_cD/qL^4$')

plt.legend()
plt.grid()

tikzplotlib.save("wc_convergence.pgf")

#plt.savefig("wc_convergence.pdf")
plt.show()
