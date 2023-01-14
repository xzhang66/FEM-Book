#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FDConvection - 1D convection equation FD analysis program.
  Problem solved : 1D convective problem with the initial shape of step.

Usage:
   >>> FDConvection file_name
or
   $ python FDConvection.py file_name

Command line arguments:
  file_name: File name in which the FD model is stored in json format

Created on Thu 5 16:56:57 2023

@author: thujsli@163.com, xzhang@tsinghua.edu.cn
"""

from sys import argv, exit
from PrePost import create_model_json, postprocess
from utitls import Apply_initial_condition, solve
import FDData as model

def FERun(DataFile):
	# Apply initial condition
	Apply_initial_condition()

	# Solution phase
	solve()


if __name__ == "__main__":
	nargs = len(argv)
	if nargs == 2:
		DataFile = argv[1]
	else:
		print("Usage ： FDConvection file_name")
		exit()

	# create FD model from DataFile in json format
	create_model_json(DataFile)

	# DataFile = "./convection_1.0.json"
	FERun(DataFile)

	# Post process
	postprocess()

