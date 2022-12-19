#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convergence analysis for the bar under concentrated force and body force.

Plot the element length - L2/energy norm cuvers in logarithm scale for both
linear and quadratic elements, and obtain their convergence rates and the
intercepts by linear regression.

Created on Fri Aug 5 21:05:47 2022

@author: jsli@163.com, xzhang@tsinghua.edu.cn
"""
from Beam1D import FERun
from Exact import ErrorNorm_Ex_6_1

import numpy as np
import matplotlib.pyplot as plt

# Json data files for 2L element
files_2L = ("2-elements.json", "4-elements.json",
            "8-elements.json", "16-elements.json")

# Run FE analysis for all files using 2L element
n2L = len(files_2L)
h2 = np.zeros(n2L)
L2Norm_d = np.zeros(n2L)
L2Norm_m = np.zeros(n2L)
L2Norm_s = np.zeros(n2L)

LinfNorm_d = np.zeros(n2L)
LinfNorm_m = np.zeros(n2L)
LinfNorm_s = np.zeros(n2L)

L1Norm_d = np.zeros(n2L)
L1Norm_m = np.zeros(n2L)
L1Norm_s = np.zeros(n2L)
for i in range(n2L):
    FERun("Convergence/"+files_2L[i])

    # Calculate error norms for convergence study
    h2[i], L2Norm_d[i], L2Norm_m[i], L2Norm_s[i], LinfNorm_d[i], LinfNorm_m[i], LinfNorm_s[i], L1Norm_d[i], L1Norm_m[i], L1Norm_s[i] = ErrorNorm_Ex_6_1()

# plot the element length - error norm curve in logarithmic scale
fig,(axs) = plt.subplots(2,2)
plt.tight_layout()

axs[0,0].set_title('error', fontsize=9); 
axs[0,0].set_ylabel('displacement', fontsize=8)
axs[0,0].xaxis.set_tick_params(labelsize=7)
axs[0,0].yaxis.set_tick_params(labelsize=7)
axs[0,0].set_xscale('log')
axs[0,0].set_yscale('log')
line1, = axs[0,0].plot(h2,L2Norm_d,label='L_2')
line2, = axs[0,0].plot(h2,LinfNorm_d,label='L_∞')
line3, = axs[0,0].plot(h2,L1Norm_d,label='L_1')
plt.legend(handles=[line1,line2,line3],labels=['L_2','L_∞','L_1'] ,loc='best')

axs[0,1].set_title('error', fontsize=9); 
axs[0,1].set_ylabel('moment', fontsize=8)
axs[0,1].xaxis.set_tick_params(labelsize=7)
axs[0,1].yaxis.set_tick_params(labelsize=7)
axs[0,1].set_xscale('log')
axs[0,1].set_yscale('log')
line1, = axs[0,1].plot(h2,L2Norm_m,label='L_2')
line2, = axs[0,1].plot(h2,LinfNorm_m,label='L_∞')
line3, = axs[0,1].plot(h2,L1Norm_m,label='L_1')

axs[1,0].set_xlabel('Element length (m)', fontsize=8); 
axs[1,0].set_ylabel('shear force', fontsize=8)
axs[1,0].xaxis.set_tick_params(labelsize=7)
axs[1,0].yaxis.set_tick_params(labelsize=7)
axs[1,0].set_xscale('log')
axs[1,0].set_yscale('log')
line1, = axs[1,0].plot(h2,L2Norm_s,label='L_2')
line2, = axs[1,0].plot(h2,LinfNorm_s,label='L_∞')
line3, = axs[1,0].plot(h2,L1Norm_s,label='L_1')

axs[1,1].set_xlabel('Element length (m)', fontsize=8); 

# Convert matplotlib figures into PGFPlots figures stored in a Tikz file, 
# which can be added into your LaTex source code by "\input{fe_plot.tex}"
import tikzplotlib
tikzplotlib.save("fe_convergence.tex")

plt.savefig("convergence.pdf")
plt.show()

# Linear regression 
print("The L2 error norms are ")

a, C = np.polyfit(np.log(h2),np.log(L2Norm_d),1)
print("    Displacement   : ||e||_L2 = %e h^%g" %(np.e**C, a))
a, C = np.polyfit(np.log(h2),np.log(L2Norm_m),1)
print("    Moment         : ||e||_L2 = %e h^%g" %(np.e**C, a))
a, C = np.polyfit(np.log(h2),np.log(L2Norm_s),1)
print("    Shear force    : ||e||_L2 = %e h^%g" %(np.e**C, a))


print("The L∞ error norms are ")

a, C = np.polyfit(np.log(h2),np.log(LinfNorm_d),1)
print("    Displacement   : ||e||_L∞ = %e h^%g" %(np.e**C, a))
a, C = np.polyfit(np.log(h2),np.log(LinfNorm_m),1)
print("    Moment         : ||e||_L∞ = %e h^%g" %(np.e**C, a))
a, C = np.polyfit(np.log(h2),np.log(LinfNorm_s),1)
print("    Shear force    : ||e||_L∞ = %e h^%g" %(np.e**C, a))


print("The L1 error norms are ")

a, C = np.polyfit(np.log(h2),np.log(L1Norm_d),1)
print("    Displacement   : ||e||_L1 = %e h^%g" %(np.e**C, a))
a, C = np.polyfit(np.log(h2),np.log(L1Norm_m),1)
print("    Moment         : ||e||_L1 = %e h^%g" %(np.e**C, a))
a, C = np.polyfit(np.log(h2),np.log(L1Norm_s),1)
print("    Shear force    : ||e||_L1 = %e h^%g" %(np.e**C, a))


