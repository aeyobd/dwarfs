#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 10:33:26 2021

@author: asya

Centers all particles and removes outliers and returns new centered positions to datafile
Syntax: inputFile outputFile radius
"""

import numpy as np
import h5py
import sys

inputFile = sys.argv[1]
outputFile = sys.argv[2]
clipRadius = sys.argv[3]

def shrinkSphere(inpDM,r0, frac, Npart, inpGC):
    """
    used in makeCentredHDF5()
    centres using shrinking spheres method

    Parameters
    ----------
    boundParticles : np.array
        array of snapshots - each snapshot containing:
            R,x,y,z,vx,vy,vz,[mostBoundX,mostBoundY,mostBoundZ],pot,ids,mdm
    r0 : float
        radius of largest sphere.
    frac : float
        fraction to decrease radius every iteration.
    Npart : int
        min number of particles inside radius.

    Returns
    -------
    newData : TYPE
        centred version of input containing:
            R,x,y,z,vx,vy,vz,x_centre,y_centre,z_centre.
    """

    boundParticles = inpDM[0]
    pos0 = inpDM[1]

    x0 = pos0[0]
    y0 = pos0[1]
    z0 = pos0[2]

    x = boundParticles[1]
    y = boundParticles[2]
    z = boundParticles[3]

    vx = boundParticles[4]
    vy = boundParticles[5]
    vz = boundParticles[6]

    while True:
        r = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)

        okIdx = np.where(r < r0)[0]
        if len(okIdx) < Npart: break

        x = x[okIdx]
        y = y[okIdx]
        z = z[okIdx]

        vx = vx[okIdx]
        vy = vy[okIdx]
        vz = vz[okIdx]

        x0 = np.average(x)
        y0 = np.average(y)
        z0 = np.average(z)

        r0 *= frac

    vx0 = np.average(vx)
    vy0 = np.average(vy)
    vz0 = np.average(vz)

    # centre DM particles
    newX = boundParticles[1]-x0
    newY = boundParticles[2]-y0
    newZ = boundParticles[3]-z0
    newR = np.sqrt(newX**2+newY**2+newZ**2)

    newVX = boundParticles[4]-vx0
    newVY = boundParticles[5]-vy0
    newVZ = boundParticles[6]-vz0

    # if GC exist, centre them as well
    try:
        newGCX = inpGC[1]-x0
        newGCY = inpGC[2]-y0
        newGCZ = inpGC[3]-z0
        newGCR = np.sqrt(newGCX**2+newGCY**2+newGCZ**2)

        newGCVX = inpGC[4]-vx0
        newGCVY = inpGC[5]-vy0
        newGCVZ = inpGC[6]-vz0

        nGC = np.array([newGCR,newGCX,newGCY,newGCZ,newGCVX,newGCVY,newGCVZ])
    except TypeError:
        nGC = None

    newData=(np.array([newR,newX,newY,newZ,newVX,newVY,newVZ]),[x0,y0,z0,vx0,vy0,vz0],nGC)
    return newData

def makeCentredHDF5(filename,clipRadius):
    with h5py.File(filename, 'a') as f:
        pot = f[('PartType1/Potential')][:]
        ids = f[('PartType1/ParticleIDs')][:]

        pos = f[('PartType1/Coordinates')][()]
        x = pos[:,0]; y = pos[:,1]; z = pos[:,2]
        R = np.sqrt(x**2 + y**2 + z**2)

        vel = f[('PartType1/Velocities')][()]
        vx = vel[:,0]; vy = vel[:,1]; vz = vel[:,2]


        mostBoundX = np.mean(x)
        mostBoundY = np.mean(y)
        mostBoundZ = np.mean(z)

        centredData = shrinkSphere((np.array([R,x,y,z,vx,vy,vz,ids]),[mostBoundX,mostBoundY,mostBoundZ]),1,0.95,1000,None)

        x = centredData[0][1] #containing all x relative to centre
        y = centredData[0][2] #containing all y relative to centre
        z = centredData[0][3] #containing all z relative to centre
        vx = centredData[0][4]
        vy = centredData[0][5]
        vz = centredData[0][6]

        pos[:,0] = x
        pos[:,1] = y
        pos[:,2] = z

        vel[:,0] = vx
        vel[:,1] = vy
        vel[:,2] = vz

        R = np.sqrt(x**2 + y**2 + z**2)
        numInitial = len(R)
        clipRadius = float(clipRadius)
        mask = R < clipRadius

        R = R[mask]
        pos = pos[mask]
        ids = ids[mask]
        vel = vel[mask]

        numFinal = len(R)
        numDeleted = numInitial - numFinal


        del f['PartType1/Coordinates']
        del f['PartType1/ParticleIDs']
        del f['PartType1/Potential']

        #create new datasets
        f.create_dataset("PartType1/Coordinates", data=pos)
        f.create_dataset("PartType1/ParticleIDs", data=ids)
        f.create_dataset("PartType1/Velocities", data=vel)

        #change number of particles in file

        attrs = f['Header'].attrs
        NumPart_ThisFile = attrs['NumPart_ThisFile']
        NumPart_Total = attrs['NumPart_Total']

        # remove deleted
        NumPart_Total[1] = NumPart_Total[1]-numDeleted
        NumPart_ThisFile[1] = NumPart_ThisFile[1]-numDeleted

        # update header
        attrs.modify('NumPart_ThisFile', NumPart_ThisFile)
        attrs.modify('NumPart_Total', NumPart_Total)

        print('Trimmed ' + str(numDeleted)+ ' particles from file')

    return

f = h5py.File(outputFile, 'w')
snap = h5py.File(inputFile, 'r')
for j in snap:
    print('Copying ' + j)
    snap.copy(j, f)

snap.close()
f.close()

