# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 13:33:35 2016

@author: waltmann
"""
import os
import sys
from math import sqrt
 
 
 
 
###########################################################################
### This function simply copies an xyz file, and renames the particle types
###########################################################################
 
def rename_xyz(save,read,name):
     
     fid = open(save,'w')
     inputfile = open(read,'r')
     fdata=inputfile.read()
     data=fdata.splitlines()
     num=data[0]
     data2=data[2:]
     fid.write(('%s \n \n')%(num))
     for line in data2:
         s = line.split()
         fid.write(('%s %s %s %s \n')%(name,s[1],s[2],s[3])) 
        
##############################################################################
### This function scales the positions of xyz files in the x y or z directions
##############################################################################
 
def scale_xyz(save,read,scalex,scaley,scalez):
 
 
    #Open input file for reading in xyz format
 
    fin=open(read,'r')
    fdata=fin.read()
    datain=fdata.splitlines()
    num=datain[0]
    data=datain[2:]
 
    #Open output file for writing in xyz format
 
    fout=open(save,'w')
    fout.write(('%s \n \n')%(num))
 
    for line in data:
        s=line.split()
        x=float(s[1])*scalex
        y=float(s[2])*scaley
        z=float(s[3])*scalez
        fout.write(("%s %f %f %f\n")%(s[0],x,y,z))
 
##############################################################################
### This function shifts the positions of xyz files in the x y or z directions
##############################################################################
 
def shift_xyz(save,read,shiftx,shifty,shiftz):
 
 
    #Open input file for reading in xyz format
 
    fin=open(read,'r')
    fdata=fin.read()
    datain=fdata.splitlines()
    num=datain[0]
    data=datain[2:]
 
    #Open output file for writing in xyz format
 
    fout=open(save,'w')
    fout.write(('%s \n \n')%(num))
 
    for line in data:
        s=line.split()
        x=float(s[1])+shiftx
        y=float(s[2])+shifty
        z=float(s[3])+shiftz
        fout.write(("%s %f %f %f\n")%(s[0],x,y,z))
 
 
### Example run script ###
 
"""        
import os
import sys
 
import manipulate_xyz as m
 
save='newname.xyz'
read='tetra.xyz'
name='Au'
 
m.rename_xyz(save,read,name)
""" 