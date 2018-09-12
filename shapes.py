#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 12:38:37 2018

@author: lee
"""
#  SHAPES:      A package of functions to output surfaces of geometric forms 
#
#  FUNCTIONS:   sphere   - output x,y,z coordinates for 2D sphere
#               cylinder - output x,y,z coordinates of cylinder 

import numpy as np

#  SPHERE:      A function to draw a 2-sphere.
#
#  Usage:       x,y,z = sphere(N)
#
#  Arguments:  N - Number of vertices in rendering, default is 20
#
#  Returns:    X - The x-coordinates of each vertex in the sphere surface.
#              Y - The y-coordinates of each vertex in the sphere surface.
#              Z - The z-coordinates of each vertex in the sphere surface.
#
#  Author:     Lee Thompson
#  Date:       17 July 2018

def sphere(n=20):
    [phi,theta] = np.mgrid[0:2*np.pi:complex(n),0:np.pi:complex(n)]
    x = np.cos(phi)*np.sin(theta)
    y = np.sin(phi)*np.sin(theta)
    z = np.cos(theta)
    return(x,y,z)

#  CYLINDER:    A function to draw a N-sided cylinder based on the
#               generator curve in the vector R. It is based on a MatLab script
#               by Luigi Barone, Per Sundqvist and Wei Pan that itself is an update of 
#               http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=5468&objectType=file 
#               translated to python.
#
#  Usage:       X,Y,Z = cylinder([R], N, M, [r1], [r2])
#
#  Arguments:  R - List of radii at each distance along cylinder. If len(R)<M
#                  remaining points are radius of last value given 
#              N - The number of points around the circumference.
#              M - The number of divisions along the length of 
#                  the cylinder, should be at least equal to len(R).
#              r1 - List of start coordinates [x,y,z]
#              r2 - List of end coordinatea [x,y,z] 
#
#  Returns:    X - The x-coordinates of each vertex in the cylinder surface.
#              Y - The y-coordinates of each vertex in the cylinder surface.
#              Z - The z-coordinates of each vertex in the cylinder surface.
#
#  Author:     Lee Thompson
#  Date:       17 July 2018

def cylinder(R=None,N=None,M=None,r1=None,r2=None):
    theta = np.linspace(0,2*np.pi,N)
    m = M
    lr = len(R)
    assert m >= lr, 'Wrong input'
    
    if m < 2:                       # Only one radius value supplied.
        R = np.array([R,R])         # Add a duplicate radius to make
        m = 2                       # a cylinder. Should error above.
 
    if m > lr:
        R.extend([R[lr-1]]*(m-lr))

    X = np.zeros((m,N))
    Y = np.zeros((m,N))
    Z = np.zeros((m,N))
    
    v = (r2-r1)/np.sqrt(np.dot((r2-r1),(r2-r1)))
    R2 = np.random.rand(3)
    x2 = v - R2/np.dot(R2,np.transpose(v))
    x2 = x2/np.sqrt(np.dot(x2,np.transpose(x2)))
    x3 = np.cross(v,x2)
    x3 = x3/np.sqrt(np.dot(x3,np.transpose(x3)))
    
    r1x = r1[0]
    r1y = r1[1]
    r1z = r1[2]
  
    r2x = r2[0]
    r2y = r2[1]
    r2z = r2[2]

    x2x = x2[0]
    x2y = x2[1]
    x2z = x2[2]    
    
    x3x = x3[0]
    x3y = x3[1]
    x3z = x3[2]  
    
    time = np.linspace(0,1,m)
    for j in np.arange(0,m):
        t = time[j]
        X[j,:] = r1x + (r2x-r1x)*t + R[j]*np.cos(theta)*x2x + R[j]*np.sin(theta)*x3x
        Y[j,:] = r1y + (r2y-r1y)*t + R[j]*np.cos(theta)*x2y + R[j]*np.sin(theta)*x3y
        Z[j,:] = r1z + (r2z-r1z)*t + R[j]*np.cos(theta)*x2z + R[j]*np.sin(theta)*x3z
        
    return(X,Y,Z)

