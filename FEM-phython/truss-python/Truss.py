#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Truss - truss FE analysis program.
  Element types  : 1-, 2- or 3-dimension linear bar element.
  Problem solved : truss whose Young's modulus, cross-sectional area are known within each element
    is deformed at nodal forces.

Usage:
   >>> Truss file_name
   
Command line arguments:
  file_name: File name in which the FE model is stored in json format

Created on Sat May 9 18:34:00 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""
from sys import argv,exit
import FEData as model
from TrussElem import TrussElem
from PrePost import create_model_json, print_stress
from utitls import assembly, solvedr

def FERun(DataFile):
    # create FE model from DataFile in json format
    create_model_json(DataFile)

    # Element matrix computations and assembly
    for e in range(model.nel):
        ke = TrussElem(e)
        assembly(e, ke)
    
    # Partition and solution
    solvedr()

    # Postprocessing
    print_stress()


if __name__ == "__main__":
    nargs = len(argv)
    if nargs == 2:
        DataFile = argv[1]
    else:
        print("Usage ： Truss file_name")
        exit()

    FERun(DataFile)