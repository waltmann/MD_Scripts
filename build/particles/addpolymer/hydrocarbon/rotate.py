# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 09:15:24 2016

@author: waltmann
"""

import sys
import os
import numpy as np

import rotate_polymer as rp

#############################
## Number of repeating units
#############################

n=10

##############
### v must be [x,y,z]
################
v=[-1,-1,-1]

#################
## File to read
#################

inputfile='HC_chain_n' + str(n) +'.xyz'

#################
## Call Function
#################

rp.align_vector(inputfile,v)



