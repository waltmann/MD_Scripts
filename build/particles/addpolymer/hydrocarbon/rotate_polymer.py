# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 08:52:46 2016

@author: waltmann
"""
from __future__ import division
import os
import sys

import numpy as np
from numpy import linalg as LA



def angle_difference(opp1,ad1,opp2,ad2):
    theta1=0
    theta2=0
    notplane=0    
    if(ad1==0):
        #print('loop1 entry')
        if(opp1>0):
            theta1=np.pi/2.0
        elif(opp1<0):
            theta1= -1*np.pi/2.0
        else:
            notplane=1
    elif(opp1==0):
        if(ad1<0):
            theta1=np.pi
    else:
        theta1=np.arctan(opp1/ad1)
    if(ad2==0):
        #print('loop 2 entry')
        if(opp2>0):
            theta2=np.pi/2.0
            #print('theta2 = pi/2')
        elif(opp2<0):
            theta2= -1*np.pi/2.0
        else:
            notplane=1
    elif(opp2==0):
        if(ad2<0):
            theta2=np.pi
    else:
        theta2=np.arctan(opp2/ad2)        
    if(notplane==0):
        return theta1-theta2
    else:
        return 0

######################
##aligns the poltymer chain to v and then returns the polymer chain as a list
########################

def align_vector(readpoly,v):
    fin1 = open(readpoly,'r')
    fout = open('testchain.xyz','w')
    f1data=fin1.read()
    data1=f1data.splitlines()
    data1=data1[2:]
    P=np.array([[0.0,0.0,0.0]])
    count=0          
    for line in data1:
        d=line.split()
        pos= [[float(d[1]),float(d[2]),float(d[3])]]
        P=np.append(P,pos,axis=0)
        count+=1;
        #fout.write(("%s %f %f %f\n")%(d[0],pos[0],pos[1],pos[2]))
    first=P[1]
    last=P[3]
    negs=0
    adjust=0
    for k in range(0,3):
        if(v[k]<0):
            negs+=1
    if(negs%2==1):
        adjust=np.pi
    pv=[float(last[0])-float(first[0]),float(last[1])-float(first[1]),float(last[2])-float(first[2])]
    thetaz=angle_difference(v[0],v[1],pv[0],pv[1])+adjust
    print('theta z is ' + str(thetaz*180/np.pi))
    thetay=angle_difference(v[2],v[0],pv[2],pv[0])+adjust
    print('theta y is ' + str(thetay*180/np.pi))
    thetax=angle_difference(v[2],v[1],pv[2],pv[1]) + adjust
    print('theta x is ' + str(thetax*180/np.pi))
    Rx=[[1,0,0],[0,np.cos(thetax),-1*np.sin(thetax)],[0,np.sin(thetax),np.cos(thetax)]]
    Ry=[[np.cos(thetay),0,np.sin(thetay)],[0,1,0],[-1*np.sin(thetay),0,np.cos(thetay)]]
    Rz=[[np.cos(thetaz),-1*np.sin(thetaz),0],[np.sin(thetaz),np.cos(thetaz),0],[0,0,1]]
    Pz=np.dot(P,Rz)
    Py=np.dot(Pz,Ry)
    Px=np.dot(Py,Rx)
    P=Px
    #print(P)
    fout.write(str(count)+'\n\n')
    for i in range(1,count+1):
        fout.write(('C %s %s %s\n')%(P[i][0],P[i][1],P[i][2]))
    return P
    