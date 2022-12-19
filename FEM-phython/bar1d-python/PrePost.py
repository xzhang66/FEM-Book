#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provide methods to create FE model for a bar from a json file, to plot the 
bar, to print stresses at Gauss points, to plot displacement and stress 
distributions obtained by FE analysis and exact solution.

Created on Sun Apr 24 18:56:57 2020

@author: xzhang@tsinghua.edu.cn
"""

import numpy as np
import json
import matplotlib.pyplot as plt
import tikzplotlib

import FEData as model
from utitls import gauss
from Exact import ExactSolution_TaperedBar, ExactSolution_CompressionBar, \
     ExactSolution_ConcentratedForce
from Bar1DElem import Nmatrix1D, Bmatrix1D


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
    model.E     = np.array(FEData['E'])
    model.body  = np.array(FEData['body'])
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
        model.P  = np.array(FEData['P'])

    # output plots
    model.plot_bar = FEData['plot_bar']
    model.plot_nod = FEData['plot_nod']
    model.nplot = model.nen*10
    model.plot_tex = FEData['plot_tex']

    # define the mesh
    model.x = np.array(FEData['x'])
    model.y = np.array(FEData['y'])  
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
            model.d[count] = model.e_bc[i] 
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

    
def plotbar():
    """ Plot the bar  """
    if model.plot_bar == 'yes' and model.nen == 3:
        for i in range(model.nel): 
            XX = np.array([
                model.x[model.IEN[0,i]-1], model.x[model.IEN[1,i]-1], 
                model.x[model.IEN[2,i]-1] ])
            YY = np.array([
                model.y[model.IEN[0,i]-1], model.y[model.IEN[1,i]-1],
                model.y[model.IEN[2,i]-1] ]) 
            plt.figure(1)
            plt.plot(XX,YY)
            plt.plot(XX,-YY) 
            plt.plot(XX, [0,0,0], '+r')   
     
            # check if user defined the plots of the global node numbering  
            if model.plot_nod == 'yes':    
                plt.text(XX[0],-1.5,str(model.IEN[0,i])) 
                plt.text(XX[1],-1.5,str(model.IEN[1,i])) 
                plt.text(XX[2],-1.5,str(model.IEN[2,i])) 
    
        plt.plot([ model.x[model.IEN[0,0]-1], model.x[model.IEN[0,0]-1]],
                 [-model.y[model.IEN[0,0]-1], model.y[model.IEN[0,0]-1]])
        plt.plot([ model.x[model.IEN[-1,-1]-1], model.x[model.IEN[-1,-1]-1]],
                 [-model.y[model.IEN[-1,-1]-1], model.y[model.IEN[-1,-1]-1]])
        plt.title('Bar Plot') 
        plt.show()
     
        # print some mesh parameters 
        print('\n  Bar Params') 
        print('No. of Elements  ' + str(model.nel)) 
        print('No. of Nodes     ' + str(model.nnp)) 
        print('No. of Equations ' + str(model.neq))
        
        
def disp_and_stress(e, d, ax1, ax2):
    """
    Print stresses at Gauss points, plot displacements and stress 
    distributions obtained by FE analysiss

    Args:
        e : (int) element number
        d : (numnp.array(neq,1)) solution vector
        ax1 : axis to draw displacement distribution
        ax2 : axis to draw stress distribution
    """
    de = model.d[model.LM[:,e]-1]  # extract element nodal displacements 
    IENe = model.IEN[:,e]-1       # extract element connectivity information 
    xe = model.x[IENe]            # extract element coordinates 
    J = (xe[-1] - xe[0])/2       # Jacobian  
    w, gp = gauss(model.ngp)      # Gauss points and weights
    
    # compute stresses at Gauss points 
    gauss_pt = np.zeros(model.ngp)
    stress_gauss = np.zeros(model.ngp)
    for i in range(model.ngp):
        xt  = 0.5*(xe[0]+xe[-1])+J*gp[i]  # Gauss point in the physical coordinates 
        gauss_pt[i] = xt         # store gauss point information   
             
        N  = Nmatrix1D(xt,xe)    # extract shape functions  
        B  = Bmatrix1D(xt,xe)    # extract derivative of shape functions  
      
        Ee = N@model.E[IENe]             # Young's modulus at element Gauss points    
        stress_gauss[i] = Ee*B@de  # compute stresses at Gauss points 
      
    # print stresses at element gauss points    
    print("%8d %12.6f %12.6f %16.6f %16.6f" 
          %(e,gauss_pt[0],gauss_pt[1],stress_gauss[0],stress_gauss[1]))

    # equally distributed coordinate within an element
    xplot = np.linspace(xe[0],xe[-1],model.nplot)               

    # compute displacements and stresses  
    displacement = np.zeros(model.nplot)
    stress = np.zeros(model.nplot)
    for i in range(model.nplot):
        xi = xplot[i]              # current coordinate   
        N  = Nmatrix1D(xi,xe)      # shape functions  
        B  = Bmatrix1D(xi,xe)      # derivative of shape functions    
         
        Ee = N@model.E[IENe]      # Young's modulus   
        displacement[i] = N@de     # displacement output 
        stress[i]       = Ee*B@de  # stress output s
         
    # plot displacements and stresses  
    line1, = ax1.plot(xplot,displacement)
    line2, = ax2.plot(xplot,stress) 
    if e == 0:
        line1.set_label('FE')
        line2.set_label('FE')


def postprocessor():
    """ 
    Print stresses at Gauss points, plot displacement and stress
    distributions obtained by FE analysis calling disp_and_stress and
    by exact solution calling ExactSolution.
    """    
    
    print()
    print('Print stresses at the Gauss points \n')
    print('%8s %12s %12s %16s %16s'
          %("Element","x(gauss1)","x(gauss2)","stress(gauss1)","stress(gauss2)"))
    print('--------------------------------------------------------------------')
    
    fig,(ax1,ax2) = plt.subplots(2,1)
    plt.tight_layout()

    ax1.set_title('FE analysis of 1D bar'); 
    ax1.set_ylabel('displacement')  

    ax2.set_xlabel('x')
    ax2.set_ylabel('stress')
    
    # loop over elements to compute the stresses
    for e in range(model.nel):
        # compute stresses and displacements for the current element
        disp_and_stress(e,model.d,ax1,ax2)
        
    # plot the exact solution
    if model.Exact == "TaperedBar":
        ExactSolution_TaperedBar(ax1,ax2)
    elif model.Exact == "CompressionBar":
        ExactSolution_CompressionBar(ax1,ax2)
    elif model.Exact == "ConcentratedForce":
        ExactSolution_ConcentratedForce(ax1, ax2)
    elif model.Exact != None:
        print('Error: Exact solution for %s is not available'%(model.Exact))

    ax1.legend()
    ax2.legend()
    plt.show()

    # Convert matplotlib figures into PGFPlots figures stored in a Tikz file, 
    # which can be added into your LaTex source code by "\input{fe_plot.tex}"
    if model.plot_tex == "yes":
        tikzplotlib.save("fe_plot.tex")