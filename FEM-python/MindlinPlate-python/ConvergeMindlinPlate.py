#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convergence analysis for Mindlin plate element. 

Plot w_c*D/q/L^4 - L/h cuvers when using full integration and selective reduced integration.

Created on Thu Apr 30 21:05:47 2020

@author: xzhang
"""

from sys import argv, exit
from PrePost import create_model_json
from MindlinPlate import FERun
import numpy as np
import FEData as model
import matplotlib.pyplot as plt
import tikzplotlib

nargs = len(argv)
if nargs == 2:
	DataFile = argv[1]
else:
	print("Usage ： MindlinPlate file_name")
	exit()

# create FE model from DataFile in json format
create_model_json(DataFile)

# calculate w_c*D/q/L^4
ratio = np.arange(5, 1000, 5)
nh = len(ratio)
wc_e = np.ones(nh) * 0.00126
wc_f = np.zeros(nh)
wc_s = np.zeros(nh)

for index, ri in enumerate(ratio):
	# the plate gets thinner as L/h increases
	model.h = model.lx / ri
	model.Db = np.array([[1, model.nu, 0],
						[model.nu, 1, 0],
						[0, 0, (1-model.nu)/2]])*model.E*model.h**3/(12.0*(1-model.nu**2))

	shcof = 5/6.0 #shear correction factor
	model.Ds = np.array([[1, 0],
						[0, 1]])*shcof*model.G*model.h

	# full integration
	model.ngp = 2

	# initialize K, d and f
	model.K = np.zeros((model.neq, model.neq))
	model.d = np.zeros((model.neq, 1))
	model.f = np.zeros((model.neq, 1))

	FERun("plate_64.json")

	wc_f[index] = model.wc

	# selective reduced integration
	model.ngp = 1

	# initialize K, d and f
	model.K = np.zeros((model.neq, model.neq))
	model.d = np.zeros((model.neq, 1))
	model.f = np.zeros((model.neq, 1))

	FERun("plate_64.json")

	wc_s[index] = model.wc

# plot the curve of w_c*D/q/L^4 vs. L/h
line1, = plt.plot(ratio, wc_s, color='b', label='Selective reduced integration')
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
