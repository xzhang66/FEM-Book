#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot stress cuvers along radius, and obtain optimal stress point.

Created on Fri Aug 5 21:05:47 2022

@author: jsli@163.com, xzhang@tsinghua.edu.cn
"""
from Elasticity2D import FERun
from Exact import Exact, sigmarr

import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib

# Json data files for reduced integration
files_1GP = ("CylindricalCavity_1_03.json", "CylindricalCavity_1_04.json",
            "CylindricalCavity_1_045.json", "CylindricalCavity_1_0499.json")

# Json data files for full integration
files_2GP = ("CylindricalCavity_2_03.json", "CylindricalCavity_2_04.json",
            "CylindricalCavity_2_045.json", "CylindricalCavity_2_0499.json")

# Run FE analysis for all files using 2L element
nplot=21

n1GP = len(files_1GP)
n2GP = len(files_2GP)

xplot1GP = np.zeros((4*nplot, n1GP))
xplot2GP = np.zeros((4*nplot, n2GP))

xplot1GP_osp = np.zeros((4, n1GP))
xplot2GP_osp = np.zeros((4, n2GP))

sigma_rr_1GP = np.zeros((4*nplot, n1GP))
sigma_rr_2GP = np.zeros((4*nplot, n1GP))

sigma_rr_1GP_osp = np.zeros((4, n1GP))
sigma_rr_2GP_osp = np.zeros((4, n1GP))

for i in range(n1GP):
    FERun("OptimalStressPoint/"+files_1GP[i])

    # Calculate error norms for convergence study
    xplot1GP[:,i], sigma_rr_1GP[:,i], xplot1GP_osp[:,i], sigma_rr_1GP_osp[:,i] = sigmarr()

for i in range(n2GP):
    FERun("OptimalStressPoint/"+files_2GP[i])

    # Calculate error norms for convergence study
    xplot2GP[:,i], sigma_rr_2GP[:,i], xplot2GP_osp[:,i], sigma_rr_2GP_osp[:,i] = sigmarr()

#Plot deflection and moment Mx distributions along the radius
#Reduced integration
fig, (ax1) = plt.subplots(1,1)
ax1.set_title('FE analysis of centerline(1 Gauss point)')
ax1.set_ylabel('sigma_rr')
ax1.set_ylim(-1.0,0.4)

line1, = ax1.plot(xplot1GP[:,0], sigma_rr_1GP[:,0],'--',label='μ=0.3')
line2, = ax1.plot(xplot1GP[:,1], sigma_rr_1GP[:,1],'--',label='μ=0.4')
line3, = ax1.plot(xplot1GP[:,2], sigma_rr_1GP[:,2],'--',label='μ=0.45')
line4, = ax1.plot(xplot1GP[:,3], sigma_rr_1GP[:,3],'--',label='μ=0.499')

#Optimal stress point
ax1.scatter(xplot1GP_osp,sigma_rr_1GP_osp,edgecolors='k')

Exact(ax1)

ax1.legend()

tikzplotlib.save("OptimalStressPoint_1GP.tex")

plt.savefig("OptimalStressPoint_1GP.pdf")
plt.show()

#Full integration
fig, (ax1) = plt.subplots(1,1)
ax1.set_title('FE analysis of centerline(2 Gauss point)')
ax1.set_ylabel('sigma_rr')
ax1.set_ylim(-1.0,0.4)

line1, = ax1.plot(xplot2GP[:,0], sigma_rr_2GP[:,0],'--',label='μ=0.3')
line2, = ax1.plot(xplot2GP[:,1], sigma_rr_2GP[:,1],'--',label='μ=0.4')
line3, = ax1.plot(xplot2GP[:,2], sigma_rr_2GP[:,2],'--',label='μ=0.45')
line4, = ax1.plot(xplot2GP[:,3], sigma_rr_2GP[:,3],'--',label='μ=0.499')

#Optimal stress point
ax1.scatter(xplot1GP_osp,sigma_rr_2GP_osp,edgecolors='k')

Exact(ax1)

ax1.legend()

tikzplotlib.save("OptimalStressPoint_2GP.tex")

plt.savefig("OptimalStressPoint_2GP.pdf")
plt.show()

