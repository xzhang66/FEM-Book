#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Advection-Diffusion - 1D FEM Program for Advection-Diffusion Problem with Galerkin and Petrov-Galerkin Method
  Element types  : 2-node linear element.
  Problem solved : 1D bar whose diffusion coefficient, cross-sectinal area and body 
  force are lineary interpolated within each element.

Usage:
   >>> Advection_Diffusion file_name
   
Command line arguments:
  file_name: File name in which the FE model is stored in json format

Created on Sun Apr 24 18:56:57 2020

@author: xzhang@tsinghua.edu.cn
"""

from sys import argv,exit

import numpy as np
import FEData as model
from Advection_DiffusionElem import Advection_DiffusionElem
from utitls import solvedr, assembly
from PrePost import create_model_json, naturalBC, postprocessor

def FERun(DataFile): 
    # create FE model from DataFile in json format
    create_model_json(DataFile)
    
    dd = np.zeros((model.neq,3))
    xx = np.zeros((model.neq,1))
    
    for i in range(model.alpha.__len__()):
        model.K = np.zeros((model.neq,model.neq))
        model.d = np.zeros((model.neq,1))
        count  = 0
        for j in range(model.neq):
            if model.flags[j] == 2:   # Essential boundary node
                count = count + 1
                model.ID[j] = count   # The reordered number of essential B.C
                model.d[count-1] = model.e_bc[j]
        
        # Element matrix computations and assembly
        for e in range(model.nel): 
            ke,fe = Advection_DiffusionElem(e,model.alpha[i]) 
            assembly(e, ke, fe)
            
        # Add nodal boundary force vector
        naturalBC()

        # Partition and solution
        solvedr()

        # Postprocessing
        postprocessor()
        
        for j in range(model.nel):
            dd[j,i] = model.d[model.IEN[0,j]-1]
        dd[model.neq-1,i] = model.d[model.IEN[1,model.nel-1]-1]
    
    for j in range(model.nel):
        xx[j] = model.x[model.IEN[0,j]-1]
    xx[model.neq-1] = model.x[model.IEN[1,model.nel-1]-1]
    
    fname = 'AD-Peclet-' + str(model.PN) + '.dat'
    fid = open(fname,'w')
    for i in range(model.neq):
        fid.write(str("%f" %xx[i,0]))
        fid.write('  ')
        fid.write(str("%f" %dd[i,0]))
        fid.write('  ')
        fid.write(str("%f" %dd[i,1]))
        fid.write('  ')
        fid.write(str("%f" %dd[i,2]))
        fid.write('\n')
    fid.close()

if __name__ == "__main__":    
    nargs = len(argv)
    if nargs == 2:
        DataFile = argv[1]
    else:
        print("Usage ： Advection_Diffusion file_name")
        exit()
    
    FERun(DataFile)