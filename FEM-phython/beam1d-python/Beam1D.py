#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Beam1D - 1D beam FE analysis program.
  Element types  : 2-node beam element.
  Problem solved : 1D beam whose Young's modulus, cross-sectinal area and body
  force are lineary interpolated within each element.

Usage:
   >>> Beam1D file_name [BeamType]

Command line arguments:
  file_name: File name in which the FE model is stored in json format
  BeamType  : Optinal. Plot the exact solution of the problem defined by BeamType
    ExactSolution_Cantilever: the cantilever given in Example 10.1 in Fish's textbook

Created on Aug. 15, 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""

from sys import argv,exit

from PrePost import create_model_json, naturalBC, plotbeam, postprocessor
from Beam1DElem import BeamElem
from utils import assembly, solvedr
import FEData as model

def FERun(DataFile):
	# create FE model from DataFile in json format
	create_model_json(DataFile)

	# plot the models
	plotbeam()

	# Element matrix computations and assembly
	for e in range(model.nel):
		ke, fe = BeamElem(e)
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
		print("Usage ： Beam1D file_name")
		exit()

	FERun(DataFile)
