## Adds a galactocentric radius and velocity to all particles
## Syntax: inputFile outputFile [x,y,z] [vx,vy,vz]

import numpy as np
import h5py
import sys
import math

inputFile = sys.argv[1]
outputFile = sys.argv[2]
oPosition = sys.argv[3] #IN CODE UNITS
oVelocity = sys.argv[4] #IN CODE UNITS

def setup(filename,positionVector,velocityVector):
    with h5py.File(filename, 'a') as f:
        positionVector=positionVector.strip('[]').split(',')
        velocityVector=velocityVector.strip('[]').split(',')

        radius = math.sqrt(float(positionVector[0])**2+float(positionVector[1])**2+float(positionVector[2])**2)
        print('Setting particles at ',radius)
        
        pos = f['PartType1/Coordinates'][()]
        vel = f['PartType1/Velocities'][()]
        pos = np.array(pos)
        vel = np.array(vel)
       
        try:
            pos4 = f['PartType4/Coordinates'][()]
            vel4 = f['PartType4/Velocities'][()]
            pos4 = np.array(pos4)
            vel4 = np.array(vel4)

            x4 = pos4[:,0];y4 = pos4[:,1];z4 = pos4[:,2]
            vx4 = vel4[:,0];vy4 = vel4[:,1];vz4 = vel4[:,2]

            x_c4 = x4 + float(positionVector[0])
            y_c4 = y4 + float(positionVector[1])
            z_c4 = z4 + float(positionVector[2])

            vx_c4 = vx4 + float(velocityVector[0])
            vy_c4 = vy4 + float(velocityVector[1])
            vz_c4 = vz4 + float(velocityVector[2])

            pos4[:,0] = x_c4; pos4[:,1] = y_c4; pos4[:,2] = z_c4
            vel4[:,0] = vx_c4; vel4[:,1] = vy_c4; vel4[:,2] = vz_c4

            del f['PartType4/Coordinates']
            del f['PartType4/Velocities']

            f.create_dataset("PartType4/Coordinates", data=pos4)
            f.create_dataset("PartType4/Velocities", data=vel4)
        except Exception:
            pass
        
        x1 = pos[:, 0];y1 = pos[:, 1];z1 = pos[:, 2]
        vx1 = vel[:,0];vy1 = vel[:,1];vz1 = vel[:,2]
        
        # the new corrected position arrays
        x_c1 = x1 + float(positionVector[0])
        y_c1 = y1 + float(positionVector[1])
        z_c1 = z1 + float(positionVector[2])

        velocity=math.sqrt(float(velocityVector[0])**2+float(velocityVector[1])**2+float(velocityVector[2])**2)
        print('Nudging particles with ',velocity)
        
        # the new corrected velocity arrays
        vx_c1 = vx1 + float(velocityVector[0])
        vy_c1 = vy1 + float(velocityVector[1])
        vz_c1 = vz1 + float(velocityVector[2])

        pos[:,0] = x_c1; pos[:,1] = y_c1; pos[:,2] = z_c1
        vel[:,0] = vx_c1; vel[:,1] = vy_c1; vel[:,2] = vz_c1

        #delete old datasets
        del f['PartType1/Coordinates']
        del f['PartType1/Velocities']

        #create new datasets
        f.create_dataset("PartType1/Coordinates", data=pos)
        f.create_dataset("PartType1/Velocities", data=vel)
        
    return

f = h5py.File(outputFile, 'w')
snap = h5py.File(inputFile, 'r')
for j in snap:
    print('Copying ' + j)
    snap.copy(j, f)

snap.close()
f.close()

setup(outputFile,oPosition,oVelocity)
