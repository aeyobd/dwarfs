import numpy as np
import astropy.units as u
from astropy.coordinates import (Galactocentric,ICRS)
import h5py
import os

# GADGET Units
G = 1
M_u = 10**10 #M_solar
R_u = 1 #kpc
T_u = 4.714*10**6 #years
V_u = 207.4 #km/s


def getCoords(r,d,D,pm_r,pm_d,rv):
    """Returns positions in kpc and velocities in km/s"""
    v_sun = [11.1,240.3+12.24, 7.25] * (u.km / u.s)
    d_sun = 8.29 *u.kpc

    c = ICRS(ra=r*u.degree,dec=d*u.degree,distance=D*u.kpc,pm_ra_cosdec=pm_r*u.mas/u.yr,pm_dec=pm_d*u.mas/u.yr,radial_velocity=rv*u.km/u.s)
    g = c.transform_to(Galactocentric(galcen_v_sun=v_sun,galcen_distance=d_sun,z_sun=0*u.pc))
    return [(g.x.value,g.y.value,g.z.value),(g.v_x.value,g.v_y.value,g.v_z.value)]

def readData(path):
    """
    run this on the output of GADGET in a MW potential
    it will check which particles in each snapshot are bound and make a new hdf5 with all simulation snapshots,
    each one containing only bound particles, and centred on the centre of your halo
    """
    try:
        #if centred hdf5 file already exists, do nothing
        centredDataFile = h5py.File(path+'centredSnapshotData.hdf5', 'r')
        return path
    except OSError:

        print('creating new centred datafile in ',path)

        #define data
        files = np.sort([path+f for f in os.listdir(path) if f.endswith('.hdf5')])
        if(path[-11:-1]=='10GyrOrbit'):
            # GADGET Milky Way potential version for some reason creates extra snapshot
            # at the end which is a duplicate?
            files = files[:-1]

        #process data
        for i,file in enumerate(files):
            print('snapshot: ',i,end='\r')
            with h5py.File(file, 'r') as f:
                ids = f[('PartType1/ParticleIDs')][:]
                pot = f[('PartType1/Potential')][:]
                mdm = f['Header'].attrs['MassTable'][1] # mass of one particle
                time = f['Header'].attrs['Time'] # time of snapshot in code units

                x = f[('PartType1/Coordinates')][:,0]
                y = f[('PartType1/Coordinates')][:,1]
                z = f[('PartType1/Coordinates')][:,2]

                vx = f[('PartType1/Velocities')][:,0]
                vy = f[('PartType1/Velocities')][:,1]
                vz = f[('PartType1/Velocities')][:,2]

                R = np.sqrt(x**2 + y**2 + z**2)

                mostBoundIndex = np.where(pot==np.min(pot))[0][0]
                mostBoundX = x[mostBoundIndex]; mostBoundY = y[mostBoundIndex]; mostBoundZ = z[mostBoundIndex]

                type1 = np.array([R,x,y,z,vx,vy,vz,ids])

                if('PartType4' in [k for k in f.keys()]):
                    x = f[('PartType4/Coordinates')][:,0]
                    y = f[('PartType4/Coordinates')][:,1]
                    z = f[('PartType4/Coordinates')][:,2]

                    vx = f[('PartType4/Velocities')][:,0]
                    vy = f[('PartType4/Velocities')][:,1]
                    vz = f[('PartType4/Velocities')][:,2]

                    R = np.sqrt(x**2 + y**2 + z**2)

                    type4 = np.array([R,x,y,z,vx,vy,vz,f[('PartType4/ParticleIDs')][:]])
                else:
                    type4 = None

            #centre snapshot using shrinking spheres
            makeCentredHDF5(type1,1.0,0.95,1000,path,i,time,type4,[mostBoundX,mostBoundY,mostBoundZ],mdm) #defined in centring.py

        return path

def getApocentresPericentres(self):
        # compute radii of the halo over the course of the simulation
        with h5py.File(self.dataPath+'centredSnapshotData.hdf5','r') as f:
            centreR = np.array([np.sqrt(f[str(snap)].attrs['centre'][0]**2+
                                        f[str(snap)].attrs['centre'][1]**2+
                                        f[str(snap)].attrs['centre'][2]**2) \
                                for snap in range(self.numSnaps)])
        apoS = []; apoR = []; apoT = []
        periS = []; periR = []; periT = []
        for i in range(1,len(centreR)-1):
            r = centreR[i]
            if(r>centreR[i-1] and r>centreR[i+1]):
                apoS.append(i)
                apoR.append(r)
                apoT.append(self.sTimes[i])
            elif(r<centreR[i-1] and r<centreR[i+1]):
                periS.append(i)
                periR.append(r)
                periT.append(self.sTimes[i])
        apoS = np.array(apoS); apoR = np.array(apoR); apoT = np.array(apoT);
        periS = np.array(periS); periR = np.array(periR); periT = np.array(periT);
        return {'apocentres':{'snapshots':apoS,'radii':apoR,'times':apoT},
                'pericentres':{'snapshots':periS,'radii':periR,'times':periT}}

def getVc(self,snapshot,rax):
        """
        Returns circular velocity profile at a given snapshot
        Parameters
        ----------
        snapshot : int
            Snapshot for which profile should be obtained.
        rax : np.ndarray
            Radii at which to compute Vc in kpc

        Returns
        -------
        vs : numpy.array
            Array of velocities for each radius in km/s.
        """
        R = self.getRadii(snapshot) #replace this with galactocentric radii
        vc = []
        for r in rax:
            n = len(np.where(R<r)[0])
            vc.append(np.sqrt(n*self.mdm/r)*V_u) #replace self.mdm with the mass of each simulation particle
        return np.array(vc)
