#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 11:59:47 2018

@author: lee
"""
#  CUBEREAD:    A function to import a cube file into variables in python
#
#  Usage:       cuberead(file)
#
#  Arguments:   file   - Filename of cube file. 
#
#  Returns:     cubemat - Voxel data, cubemat(NX,NY,NZ)
#               orig - origin of voxel data, orig(x,y,z)
#               Nstep - cartesian vector step of voxel, Nstep(x,y,z)
#               Atom - Atoms name vector based on proton mass, Atom(i)
#               AtMass - Proton mass vector, AtMass(i)
#               AtN - N Coordinate of ith atom, AtN(i)
#               NAtoms - Number of Atoms, scalar
#               NX - number of data points along x basis vector, scalar
#               NY - number of data points along y basis vector, scalar
#               NZ - number of data points along z basis vector, scalar
#
#  Author:     Lee Thompson
#  Date:       17 July 2018

import numpy as np

def cuberead(file=None):
    fid = open(file,'r')
    out = fid.readlines()

    head = out[2:6]
    header = []
    for i in head:
        header.append(i.split())
    
    NAtoms = int(header[0][0])
    NX = int(header[1][0])
    NY = int(header[2][0])
    NZ = int(header[3][0])
    orig = np.array([header[0][1],header[0][2],header[0][3]]).astype(np.float)
    Xstep = np.array([header[1][1],header[1][2],header[1][3]]).astype(np.float)
    Ystep = np.array([header[2][1],header[2][2],header[2][3]]).astype(np.float)
    Zstep = np.array([header[3][1],header[3][2],header[3][3]]).astype(np.float)
    
    Atom = np.zeros(NAtoms);
    AtMass = np.zeros(NAtoms);
    AtX = np.zeros(NAtoms);
    AtY = np.zeros(NAtoms);
    AtZ = np.zeros(NAtoms);

    atoms = out[6:6+NAtoms]
    molspec = []
    for i in atoms:
        molspec.append(i.split())
    
    if NAtoms != 0:
        for i in range(NAtoms):
            Atom[i] = int(molspec[i][0])
            AtMass[i] = float(molspec[i][1])
            AtX[i] = float(molspec[i][2])
            AtY[i] = float(molspec[i][3])
            AtZ[i] = float(molspec[i][4])          
        
    out = out[6+NAtoms:len(out)]
    cube = []
    for i in out:
        cube.append(i.split())
    fid.close() 
    
    cube = np.reshape(np.transpose(np.array(cube).astype(np.float)),(-1,1))    
#    voxels = np.reshape(cube,(-1,1))    
#    voxels[np.any(np.isnan),axis=2,:] = None

#    cubemat = np.zeros([NY,NX,NZ])
    cubemat = np.zeros([NX,NY,NZ])
    i = 0
    for ix in range(NX):
        for iy in range (NY):
            for iz in range(NZ):
#                cubemat[iy,ix,iz] = cube[i]
                cubemat[ix,iy,iz] = cube[i]
                i = i + 1
   
    return(cubemat,orig,Xstep,Ystep,Zstep,Atom,AtMass,AtX,AtY,AtZ,NAtoms,NX,NY,NZ)
