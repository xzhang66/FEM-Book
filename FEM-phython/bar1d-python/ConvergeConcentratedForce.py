#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convergence analysis for the bar under concentrated force and body force.

Plot the element length - L2/energy norm cuvers in logarithm scale for both
linear and quadratic elements, and obtain their convergence rates and the
intercepts by linear regression.

Created on Thu Apr 30 21:05:47 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""
from Bar1D import FERun
from Exact import ErrorNorm_ConcentratedForce

import numpy as np
import matplotlib.pyplot as plt

# Json data files for the concentrated force in an element
files_even = ("2_element_even.json", "4_element_even.json", 
              "8_element_even.json", "16_element_even.json")

# Json data files for the concentrated force on a node
files_node = ("2_element_node.json", "4_element_node.json", 
              "8_element_node.json", "16_element_node.json")

# Run FE analysis for all files of the concentrated force in an element
n_even = len(files_even)
h_even = np.zeros(n_even)
L2Norm_even = np.zeros(n_even)
EnNorm_even = np.zeros(n_even)
for i in range(n_even):
	FERun("Convergence/ConcentratedForce/" + files_even[i])

	# Calculate error norms for convergence study
	h_even[i], L2Norm_even[i], EnNorm_even[i] = ErrorNorm_ConcentratedForce(False)

# Run FE analysis for all files of the concentrated force on a node
n_node = len(files_node)
h_node = np.zeros(n_node)
L2Norm_node = np.zeros(n_node)
EnNorm_node = np.zeros(n_node)
for i in range(n_node):
	FERun("Convergence/ConcentratedForce/" + files_node[i])

	# Calculate error norms for convergence study
	h_node[i], L2Norm_node[i], EnNorm_node[i] = ErrorNorm_ConcentratedForce(True)

# plot the element length - error norm curve in logarithmic scale
fig,(axs) = plt.subplots(2,2)
plt.tight_layout()

axs[0,0].set_title('Concentrated force in an element', fontsize=9);
axs[0,0].set_ylabel('L_2 error', fontsize=8)
axs[0,0].xaxis.set_tick_params(labelsize=7)
axs[0,0].yaxis.set_tick_params(labelsize=7)
axs[0,0].set_xscale('log')
axs[0,0].set_yscale('log')
axs[0,0].plot(h_even,L2Norm_even)

axs[0,1].set_title('Concentrated force on a node', fontsize=9);
axs[0,1].set_ylabel('L_2 error', fontsize=8)
axs[0,1].xaxis.set_tick_params(labelsize=7)
axs[0,1].yaxis.set_tick_params(labelsize=7)
axs[0,1].set_xscale('log')
axs[0,1].set_yscale('log')
axs[0,1].plot(h_node,L2Norm_node)

axs[1,0].set_xlabel('Element length (m)', fontsize=8);
axs[1,0].set_ylabel('Energy error', fontsize=8)
axs[1,0].xaxis.set_tick_params(labelsize=7)
axs[1,0].yaxis.set_tick_params(labelsize=7)
axs[1,0].set_xscale('log')
axs[1,0].set_yscale('log')
axs[1,0].plot(h_even,EnNorm_even)

axs[1,1].set_xlabel('Element length (m)', fontsize=8)
axs[1,1].set_ylabel('Energy error', fontsize=8)
axs[1,1].xaxis.set_tick_params(labelsize=7)
axs[1,1].yaxis.set_tick_params(labelsize=7)
axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')
axs[1,1].plot(h_node,EnNorm_node)

# Linear regression
print("The L2/energy error norms are ")
a, C = np.polyfit(np.log(h_even),np.log(L2Norm_even),1)
print("    Concentrated force in an element    : ||e||_L2 = %e h^%g" %(np.e**C, a))
a, C = np.polyfit(np.log(h_node),np.log(L2Norm_node),1)
print("    Concentrated force on a node        : ||e||_L2 = %e h^%g" %(np.e**C, a))
a, C = np.polyfit(np.log(h_even),np.log(EnNorm_even),1)
print("    Concentrated force in an element    : ||e||_en = %e h^%g" %(np.e**C, a))
a, C = np.polyfit(np.log(h_node),np.log(EnNorm_node),1)
print("    Concentrated force on a node        : ||e||_en = %e h^%g\n" %(np.e**C, a))

# Convert matplotlib figures into PGFPlots figures stored in a Tikz file,
# which can be added into your LaTex source code by "\input{fe_plot.tex}"
import tikzplotlib
tikzplotlib.save("fe_convergence.tex")