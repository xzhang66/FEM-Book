#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bar1D - 1D bar FE analysis program.
  Element types  : 2-node linear and 3-node quadratic element.
  Problem solved : 1D bar whose Young's modulus, cross-sectinal area and body 
  force are lineary interpolated within each element.

Usage:
   >>> Bar1D file_name [BarType]
   
Command line arguments:
  file_name: File name in which the FE model is stored in json format
  BarType  : Optinal. Plot the exact solution of the problem defined by BarType
    TaperedBar: The tapered elastic bar given in Example 5.2 in Fish's textbook
    CompressionBar: The bar under compression in Figure 5.13 in Fish's textbook

Created on Sun Apr 24 18:56:57 2020

@author: xzhang@tsinghua.edu.cn
"""

from sys import argv,exit

import FEData as model
from Bar1DElem import BarElem
from utitls import solvedr, assembly
from PrePost import create_model_json, naturalBC, plotbar, postprocessor

def FERun(DataFile): 
    # create FE model from DataFile in json format
    # create_model_json('Convergence/16-elements-3Q.json')
    create_model_json(DataFile)
    
    # plot the models
    plotbar()
    
    # Element matrix computations and assembly
    for e in range(model.nel): 
        ke,fe = BarElem(e) 
        assembly(e, ke, fe)
    
    # Add nodal boundary force vector
    naturalBC()
    
    # Partition and solution
    solvedr()
    
    # Postprocessing
    postprocessor()

    
if __name__ == "__main__":    
    nargs = len(argv)
    if nargs == 2:
        DataFile = argv[1]
    else:
        print("Usage ： Bar1D file_name")
        exit()
    
    FERun(DataFile)