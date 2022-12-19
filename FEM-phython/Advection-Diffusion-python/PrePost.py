#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provide methods to create FE model for a bar from a json file, to plot 
displacement distributions obtained by FE analysis and exact solution.

Created on Sun Apr 24 18:56:57 2020

@author: xzhang@tsinghua.edu.cn
"""

import numpy as np
import json
import matplotlib.pyplot as plt
import tikzplotlib

import FEData as model
from Exact import ExactSolution

def create_model_json(DataFile):
    """ 
    Initialize the FEM model from file DataFile (in json format)
    """

    # input data from json file
    with open(DataFile) as f_obj:
        FEData = json.load(f_obj)

    model.Title= FEData['Title']
    model.nsd  = FEData['nsd']
    model.ndof = FEData['ndof']
    model.nnp  = FEData['nnp']
    model.nel  = FEData['nel']
    model.nen  = FEData['nen']    
    model.neq  = model.ndof*model.nnp

    # initialize K, d and f 
    model.f = np.zeros((model.neq,1))            
    model.d = np.zeros((model.neq,1))        
    model.K = np.zeros((model.neq,model.neq))    

    # element and material data (given at the element nodes)
    model.body  = np.array(FEData['body'])
    model.CArea = np.array(FEData['CArea'])
    model.k     = np.array(FEData['k'])
    model.PN    = np.array(FEData['PN'])
    model.alpha = np.array(FEData['alpha'])

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
        model.P  = np.array(FEData['P'])

    # output plots
    model.nplot = model.nen*10
    model.plot_tex = FEData['plot_tex']

    # define the mesh
    model.x = np.array(FEData['x'])
    model.IEN = np.array(FEData['IEN'])

    model.ID  = np.zeros(model.neq,np.int)
    model.LM  = np.zeros((model.nen,model.nel),np.int)   
        
    if 'Exact' in FEData:
        model.Exact = FEData['Exact']
    else:
        model.Exact = None

    # generate LM and ID arrays
    setup_ID_LM()


def setup_ID_LM():
    """ Setup ID and LM arrays """
    count  = 0
    count1 = 0
    
    # Reorder the D.O.F. to make the essential B.C. numbered first
    for i in range(model.neq):
        if model.flags[i] == 2:   # Essential boundary node
            count = count + 1
            model.ID[i] = count   # The reordered number of essential B.C
            model.d[count-1] = model.e_bc[i] 
        else:
            count1 = count1 + 1
            model.ID[i] = model.nd + count1

    for i in range(model.nel):
        for j in range(model.nen):
            model.LM[j,i] = model.ID[model.IEN[j,i]-1]


def naturalBC():
    """ Compute and assemble nodal boundary force vector """
    for i in range(model.neq):
       if model.flags[i] == 1:
          node = model.ID[i]-1
          model.f[node] += model.CArea[node]*model.n_bc[node]

def postprocessor():
    """ 
    Plot displacement distributions obtained by FE analysis and by exact solution.
    """    
    
    for i in range(model.nel):
        if i == 0:
            plt.plot(model.x[model.IEN[:,i]-1], model.d[model.IEN[:,i]-1], '-*k', label='FE')
        else:
            plt.plot(model.x[model.IEN[:,i]-1], model.d[model.IEN[:,i]-1], '-*k')
    
    # plot the exact solution
    ExactSolution()

    plt.legend()
    plt.show()

    # Convert matplotlib figures into PGFPlots figures stored in a Tikz file, 
    # which can be added into your LaTex source code by "\input{fe_plot.tex}"
    if model.plot_tex == "yes":
        tikzplotlib.save("displacement.tex")