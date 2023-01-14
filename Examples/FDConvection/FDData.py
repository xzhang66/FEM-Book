#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Global variables defining the FDM model
  Title : (str) Title of the problem to be solved.
  x_low : (float) Lower boundary of the computational domain.
  x_up  : (float) Upper boundary of the computational domain.
  dx    : (float) Space step length.
  dt    : (float) Time step length.
  ratio : (float) The ratio of time step length and space step length.
  t_end : (float) Total computational time for the problem.
  
  nx    : (int) Number of nodes in the computational domain.
  nt    : (int) Number of time steps.
  x     : (numpy.array(nx)) x coordinates of nodes in the computational domain.
  u_now : (numpy.array(nx)) Solution vector of current time step i.
  u_next: (numpy.array(nx)) Solution vector of next time step i+1.
  t     : (float) current time.

  plot_curve : plot displacement along the computational domain?
  plot_tex        : plot in latex tikz format ?
  plot_region     : (float) (numpy.array(2)) Postprocess domain of the problem.
"""
# Initialize time and space data
Title = None
x_low = 0.0
x_up = 0.0
dx = 0.0
dt = 0.0
ratio = 0.0
t_end = 0.0

nx = 0
nt = 0
x = None
u_now = None
u_next = None
t = 0.0

# Initialize parameters for postprocessing
plot_curve = None
plot_tex = None
plot_region = None
