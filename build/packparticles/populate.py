# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:30:49 2016

@author: nathan"""


import pack_particles as pp
import points

##############################
#  Particle file for import  #
##############################

particle_name='fccnp_cutbcc_fullpoly.xyz'

###########################################################
#  Number of particles to be initialized (must be a cube) #
###########################################################

A=64

################################################
#  Pack particles writing with cube root of A  #
################################################

pack = pp.read_particle([A],[particle_name])
pack.write_particles(4*4,4,4)
