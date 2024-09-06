#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: Oliver Rice

Python wrapper for running multiple instances of lare2d, varying certain parameters

"""


import os
import shutil
import numpy as np
import sys
from numpy import random
import time
from scipy.io import netcdf_file

#os.system('killall mf2d')

run = 0
ncores = 1

if os.uname()[1][-14:] == 'ham8.dur.ac.uk':
    hamilton_flag = 1
else:
    hamilton_flag = 0

nx = 64
ny = 64
nz = 64

x0 = -0.5; x1 = 0.5
y0 = -0.5; y1 = 0.5
z0 = -1.0/ny; z1 = 1.0

shearfact = 1.0
bfact = 1.0

nplots = 10
ndiags = 100
tmax = 1.0

eta = 1e-6

nu0_decay = 0.0

variables = np.zeros((30))

variables[0] = run
variables[1] = nx
variables[2] = tmax
variables[3] = nplots
variables[4] = ndiags
variables[5] = bfact
variables[6] = shearfact
variables[7] = x0
variables[8] = x1
variables[9] = y0
variables[10] = y1
variables[11] = z0
variables[12] = z1
variables[13] = ny
variables[14] = nx
variables[15] = nz

variables[16] = hamilton_flag


os.system('make')

np.savetxt('parameters/variables%03d.txt' % run, variables)   #variables numbered based on run number (up to 1000)

#Create initial condition using new init.py (potential field with arbitrary lower boundary and domain dimensions)
#compute_initial_condition(Grid(), lbound_fn, run)

if hamilton_flag < 0.5:
    os.system('/usr/lib64/openmpi/bin/mpiexec -np %d ./bin/lare3d %d' % (ncores, run))
