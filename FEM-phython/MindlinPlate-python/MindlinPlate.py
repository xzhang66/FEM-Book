#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MindlinPlate - Elastic plate bending FE analysis program.
  Element types  : 4-node quadratic Mindlin plate element.
  Problem solved : plate whose Young's modulus, Poisson's rate and thickness are known
    is deformed at surface pressure, boundary forces, nodal forces and body force.

Usage:
   >>> MindlinPlate file_name
or
   $ python MindlinPlate.py file_name

Command line arguments:
  file_name: File name in which the FE model is stored in json format

Created on Fri August 15 18:56:57 2021

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""

from sys import argv, exit
from PrePost import create_model_json, point_and_trac, postprocess
from MindlinPlateElem import MindlinPlateElem
import FEData as model
from utitls import assembly, solvedr

def FERun(DataFile):
	# Calculation and assembly of element matrices
	for e in range(model.nel):
		ke, fe = MindlinPlateElem(e)
		assembly(e, ke, fe)

	# Compute and assemble nodal boundary force vector and point forces
	point_and_trac()

	# Solution phase
	f_E = solvedr()


if __name__ == "__main__":
	nargs = len(argv)
	if nargs == 2:
		DataFile = argv[1]
	else:
		print("Usage ： MindlinPlate file_name")
		exit()

	# create FE model from DataFile in json format
	create_model_json(DataFile)

	# DataFile = "./plate_64.json"
	FERun(DataFile)

	# Post process
	postprocess()

