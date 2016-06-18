# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 09:16:36 2016

@author: nathan
"""

from __future__ import division
 
import numpy as np

import nanoparticle_core as npc


uc_size=2.5
uc_bcc = npc.NanoBcc(uc_size)
#uc_fcc = npc.NanoFcc(uc_size)

searchquery1='C'
vec_all=[]

with open('basefcc.xyz') as f1:
    lines=f1.readlines()
    for i, line in enumerate(lines):
        if line.startswith(searchquery1):
            split=str(line).split()
            vec_all.append(np.array([split[1],split[2],split[3]]).astype(np.float))
type = 'C'
with open("fccnp_cutbcc.xyz",'w') as f3:    
    for ind, vec in enumerate(vec_all):
        if uc_bcc.check_point(np.array(vec)):
            print(vec)
            f3.write('%s %1.6f %1.6f %1.6f'% (type, vec[0], vec[1], vec[2])+'\n')
f3.close()
            
