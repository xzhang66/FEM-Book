#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plate - Elastic thin plate bending FE analysis program.
  Element types  : 4-node quadratic plate element.
  Problem solved : thin plate whose Young's modulus, Poisson's rate and thickness are known
    is deformed at surface pressure, boundary forces, nodal forces and body force.

Usage:
   >>> Plate file_name
or
   $ python Plate.py file_name

Command line arguments:
  file_name: File name in which the FE model is stored in json format

Created on Fri Jun 19 18:56:57 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""
from sys import argv, exit
from PrePost import create_model_json, point_and_trac, postprocess
from PlateElem import PlateElem
import FEData as model
from utitls import assembly, solvedr


def FERun(DataFile):
	# create FE model from DataFile in json format
	create_model_json(DataFile)

	# Calculation and assembly of element matrices
	for e in range(model.nel):
		ke, fe = PlateElem(e)
		assembly(e, ke, fe)

	# Compute and assemble nodal boundary force vector and point forces
	point_and_trac()

	# Solution phase
	f_E = solvedr()

	# Post process
	postprocess()


if __name__ == "__main__":
	nargs = len(argv)
	if nargs == 2:
		DataFile = argv[1]
	else:
		print("Usage ： Elasticity2D file_name")
		exit()

	# DataFile = "./plate_16.json"
	FERun(DataFile)
