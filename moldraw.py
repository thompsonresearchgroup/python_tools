#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 17:46:32 2018

@author: lee
"""
#  MOLDRAW:     A package of functions to draw molecules
#
#  FUNCTIONS:   xyzplot  - plot molecule from an xyz file
#               molplot  - plot molecule from molecule specification variables 

import numpy as np
import matplotlib.pyplot as plt
import shapes
from set_axes_equal import set_axes_equal
from mpl_toolkits.mplot3d import axes3d, Axes3D

#  XYZPLOT:     A function to draw a molecule from an xyz file. It is based on 
#               a MatLab script by Lee Thompson translated to python.
#
#  Usage:       xyzplot(file,scale,cscale)
#
#  Arguments:  file   - Filename containing coordinates in xyz form. 
#              scale  - Diameter of spheres representing atoms.
#              cscale - Diameter of cylinders representing bonds. 
#
#  Returns:    ax     - Axis handle for figure
#
#  Author:     Lee Thompson
#  Date:       17 July 2018

def xyzplot(file=None,scale=None,cscale=None):
    fid = open(file,'r')
    out = fid.readlines()
    out = out[2:len(out)]
    xyz = []
    for i in out:
        xyz.append(i.split())
    fid.close()

    AtNam = np.array(xyz)[:,0]
    AtX = (np.array(xyz)[:,1]).astype(np.float)
    AtY = (np.array(xyz)[:,2]).astype(np.float)
    AtZ = (np.array(xyz)[:,3]).astype(np.float)
    NAtoms = len(AtNam)
    
    fig,ax = molplot(AtNam,AtX,AtY,AtZ,NAtoms,scale,cscale)
    return(fig,ax)
    
#  MOLPLOT:     A function to draw a molecule from molecule specification data
#
#  Usage:       molplot(AtNam,AtX,AtY,AtZ,NAtoms,scale,cscale)
#
#  Arguments:  AtNam  - List of atom names.
#              AtX    - List of X coordinates.
#              AtY    - List of Y coordinates.
#              AtZ    - List of Z coordinates.
#              NAtoms - Number of Atoms.
#              scale  - Diameter of spheres representing atoms.
#              cscale - Diameter of cylinders representing bonds. 
#
#  Returns:    ax     - Axis handle for figure
#
#  Author:     Lee Thompson
#  Date:       17 July 2018
    
def molplot(AtNam,AtX,AtY,AtZ,NAtoms,scale,cscale):
        
    params = {'H':([0.80,0.80,0.80],0.25*scale),
              'C':([0.55,0.55,0.55],0.4*scale),
              'N':([0.01,0.01,0.90],0.4*scale),
              'O':([1.00,0.00,0.00],0.4*scale),
              'F':([0.70,1.00,1.00],0.4*scale),
              'S':([1.00,0.78,0.16],0.5*scale),
              'Pd':([0.00,0.41,0.52],0.7*scale)}
    
    fig = plt.figure()
#    ax = fig.gca(projection='3d')
    ax = Axes3D(fig)
    ax.set_aspect('equal')
    
    for i in range(NAtoms):
        color, r = params.get(AtNam[i],([0.00,1.00,1.00],0.4*scale))
        x,y,z = shapes.sphere(50)
        ax.plot_surface(r*x+AtX[i],r*y+AtY[i],r*z+AtZ[i],color=color,facecolor=color)
        for j in range(i+1,NAtoms):
            dist = np.sqrt((AtX[i]-AtX[j])**2 + (AtY[i]-AtY[j])**2 + (AtZ[i]-AtZ[j])**2)
            if dist < 1.6:
                x,y,z = shapes.cylinder([0.06*cscale],20,2,np.array([AtX[i],AtY[i],AtZ[i]]),np.array([AtX[j],AtY[j],AtZ[j]]))
#                ax.plot_surface(x,y,z,color=[0.50,0.50,0.50])
                ax.plot_surface(x,y,z,color=[0.00,0.00,0.00])
    
    set_axes_equal(ax)
    return(fig,ax)
   
