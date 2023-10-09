## Rescales positions, velocities, and mass from natural units to GADGET simulation units
## Then returns new rescaled data to datafile
## Syntax: inputFile outputFile physicalScaleRadius massInPhysicalScaleRadius

import numpy as np
import h5py
import sys
import math
import random
import matplotlib.pyplot as plt

inputFile = sys.argv[1]
outputFile = sys.argv[2]
scaleRadius = float(sys.argv[3]) # scale radius (2.11 kpc for Fornax)
scaleMass = float(sys.argv[4]) # mass eclosed in scale radius (~2*10^8 M_sol for Fornax)

## G = 1, M_0 = 10^10 solar masses, R_0 = 1 kpc, T_0 = 4.718*10^6 yrs, V_0 = 207.4 km/s
G = 1
M_0 = 10**10
R_0 = 1
T_0 = 4.718*10**6
V_0 = 207.4

def rescale(filename,radius,mass):
    with h5py.File(filename, 'a') as f:
        pos = f['PartType1/Coordinates'][()]
        ids = f['PartType1/ParticleIDs'][()]
        pot = f['PartType1/Potential'][()]
        vel = f['PartType1/Velocities'][()]
        acc = f['PartType1/Acceleration'][()]
        mdm = f['Header'].attrs['MassTable'][1]

        pos = np.array(pos)
        ids = np.array(ids)
        pot = np.array(pot)
        vel = np.array(vel)
        acc = np.array(acc)

        ## Natural Units: (scaleR = 1, M(<scaleR) = 1)
        x = pos[:,0]
        y = pos[:,1]
        z = pos[:,2]

        vx = vel[:,0]
        vy = vel[:,1]
        vz = vel[:,2]

        ax = acc[:,0]
        ay = acc[:,1]
        az = acc[:,2]

        ## Rescaling
        x_s = x*radius/R_0
        y_s = y*radius/R_0
        z_s = z*radius/R_0

        r_ss = 1*radius/R_0
        M_s = 1*mass/M_0

        vx_s = vx*np.sqrt(G*M_s/r_ss)
        vy_s = vy*np.sqrt(G*M_s/r_ss)
        vz_s = vz*np.sqrt(G*M_s/r_ss)
        
        ax_s = np.zeros(len(ax))
        ay_s = np.zeros(len(ay))
        az_s = np.zeros(len(az))

        mdm_s = mdm*mass/M_0

        ## Update arrays with simulation units
        pos[:,0] = x_s
        pos[:,1] = y_s
        pos[:,2] = z_s

        vel[:,0] = vx_s
        vel[:,1] = vy_s
        vel[:,2] = vz_s

        acc[:,0] = ax_s
        acc[:,1] = ay_s
        acc[:,2] = az_s

        ## Delete old data
        del f['PartType1/Coordinates']
        del f['PartType1/Velocities']
        del f['PartType1/Acceleration']

        ## Add rescaled data
        f.create_dataset("PartType1/Coordinates", data=pos)
        f.create_dataset("PartType1/Velocities", data=vel)
        f.create_dataset("PartType1/Acceleration", data=acc)
        
        attrs = f['Header'].attrs
        MassTable = attrs['MassTable']
        MassTable[1] = mdm_s
        attrs.modify('MassTable', MassTable)
    return

f = h5py.File(outputFile, 'w')
snap = h5py.File(inputFile, 'r')
for j in snap:
    print('Copying ' + j)
    snap.copy(j, f)

snap.close()
f.close()

rescale(outputFile,scaleRadius,scaleMass)
