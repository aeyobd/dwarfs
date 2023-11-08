# Integrator for a fully analytical McMillan2017-like potential from Rapha
# January 2020

# ALL FUNCTION INPUTS ARE IN PHYSICAL UNITS; INTERNAL CALCULATIONS ARE PERFORMED USING 'GADGET' UNITS DEFINED BELOW

import numpy as np
from numpy.linalg import norm

from scipy.special import erf
from scipy.special import spence
from scipy import integrate
import math

""" GADGET UNITS"""
G_u = 1
M_u = 1e10       # M_solar
R_u = 1          # kpc
T_u = 4.714e6    # years
V_u = 207.4      # km/s
A_u = 1.394e-12  # km/s^2

""" UNIVERSAL CONSTANTS """
G = 4.30091e-3 # pc (km/s)^2 / Msol

""" POTENTIAL CONSTANTS"""
# NFW halo
r_s = 20.2                                        # scale radius of MW [kpc]
M_r_s = 1.53e11/M_u                               # mass within the scale radius of MW
rho_0 = (M_r_s)/(4*np.pi*r_s**3*(np.log(2)-0.5))

# Hernquist stellar bulge
M_bulge = 2.1e10/M_u                              # mass of bulge [Msol]
a_bulge = 1.3                                     # scale length [kpc]

# thin disk component
a_thin = 3.9                                      # scale length [kpc]
b_thin = 0.31                                     # scale height [kpc]
M_thin = 5.9e10/M_u                               # mass

# thick disk component
a_thick = 4.4                                    # scale length [kpc]
b_thick = 0.92                                   # scale height [kpc]
M_thick = 2.0e10/M_u                             # mass

## --------------------------------------------------------- POTENTIALS -------------------------------------------------------------------------------

def bulgePotential(x,y,z):
    r = math.sqrt(x**2+y**2+z**2)
    pot = -G_u*M_bulge/(r+a_bulge)
    return pot

def diskPotential(x,y,z):
    R = np.sqrt(x**2+y**2)
    pot1 = -G_u*M_thin/math.sqrt(R**2+(a_thin+math.sqrt(z**2+b_thin**2))**2)
    pot2 = -G_u*M_thick/math.sqrt(R**2+(a_thick+math.sqrt(z**2+b_thick**2))**2)
    pot = pot1+pot2
    return pot

def haloPotential(x,y,z):
    r = math.sqrt(x**2+y**2+z**2)
    pot = -4*np.pi*G_u*rho_0*r_s**3/r*np.log(1+r/r_s)
    return pot

def totPot(x,y,z):
    return bulgePotential(x,y,z)+diskPotential(x,y,z)+haloPotential(x,y,z)

def totPotS(r):
    hP = -4*np.pi*G_u*rho_0*r_s**3/r*np.log(1+r/r_s)
    bP = -G_u*M_bulge/(r+a_bulge)
    dP = diskPotential(r,0,0)
    return hP+bP+dP

## --------------------------------------------------------- STUFF -----------------------------------------------------------------

def getMiyamotoNagai(x,y,z,physicalUnits=None):
    """Takes position [kpc] and returns acceleration [GADGET] or [km/s^2] in x, y, z from [total disk, thin disk, thick disk]"""
    #thin component
    xi_thin = np.sqrt(z**2+b_thin**2)
    den_thin = np.sqrt(x**2+y**2+(a_thin+xi_thin)**2)
    ax_thin = -G_u*M_thin*x/den_thin**3
    ay_thin = -G_u*M_thin*y/den_thin**3
    az_thin = -G_u*M_thin*z/den_thin**3*(1.+a_thin/xi_thin)

    #thick component
    xi_thick = np.sqrt(z**2+b_thick**2)
    den_thick = np.sqrt(x**2+y**2+(a_thick+xi_thick)**2)
    ax_thick = -G_u*M_thick*x/den_thick**3
    ay_thick = -G_u*M_thick*y/den_thick**3
    az_thick = -G_u*M_thick*z/den_thick**3*(1.+a_thick/xi_thick)

    axTot = ax_thin + ax_thick
    ayTot = ay_thin + ay_thick
    azTot = az_thin + az_thick
    if physicalUnits:
        return [(axTot,ayTot,azTot)*A_u,(ax_thin,ay_thin,az_thin)*A_u,(ax_thick,ay_thick,az_thick)*A_u]
    else:
        return [(axTot,ayTot,azTot),(ax_thin,ay_thin,az_thin),(ax_thick,ay_thick,az_thick)]

"""Dynamical friction stuff from Rapha"""
def Li2(x): return spence(1-x)

def rhoNFW(r,z):
    Msz = 1.53*10**11/M_u
    Rsz = 20.2
    return 1./(4.*np.pi*r*(r+Rsz)**2.)*Msz/(np.log(2.)-0.5)

def sigmaNFW(r,z):
    Msz = 1.53*10**11/M_u
    Rsz = 20.2
    x= r/Rsz
    sigma2adim = x*(x+1)/(np.log(4.)-1.)*(np.pi**2.-np.log(x)+3*(np.log(x+1))**2.+np.log(x+1)*(1+x**-2.-4./x-2./(x+1))-1./x-1./(x+1)**2.-6./(x+1)+6.*Li2(-x))
    return np.sqrt(G_u*Msz/Rsz*sigma2adim)

def adynmag(m, rho, sig, v, lnL):
    Y = norm(v)/(np.sqrt(2.)*sig)
    return -4.*np.pi*lnL*G_u*G_u*m*rho*(erf(Y)-2.*Y/np.sqrt(np.pi)*np.exp(-Y**2.))/(norm(v)**2)

def adynNFW(r,v,m,z):
    """Dynamical friction components"""
    rho = rhoNFW(norm(r),z)
    sig = sigmaNFW(norm(r),z)
    lnL = 2.2
    return adynmag(m,rho,sig,v,lnL)*v/norm(v)

def getMassNFW_H(r):
    """Takes radius [kpc] and returns bulge & halo mass [GADGET] enclosed within that radius"""
    M_r_NFW = 4*np.pi*rho_0*r_s**3*(np.log((r_s+r)/r_s)-r/(r_s+r))
    M_r_H = M_bulge*r**2/(r+a_bulge)**2
    return M_r_NFW,M_r_H

## ---------------------------------------------------------- CIRCULAR VELOCITIES ---------------------------------------------------------------------

def totVc(r,physicalUnits=None):
    """Returns total circular velocity [GADGET] or [km/s] on galactic plane"""
    a = np.sqrt(acceleration(r,0,0)[0]**2+acceleration(r,0,0)[1]**2+acceleration(r,0,0)[2]**2)
    if physicalUnits:
        return np.sqrt(a*r)*V_u
    return np.sqrt(a*r)

def haloVc(r,physicalUnits=None):
    """Returns halo circular velocity [km/s] on galactic plane"""
    m = getMassNFW_H(r)[0]
    vc = np.sqrt(G_u*m/r)
    if physicalUnits:
        return vc*V_u
    return vc

def bulgeVc(r,physicalUnits=None):
    """Returns bulge circular velocity [km/s] on galactic plane"""
    m = getMassNFW_H(r)[1]
    vc = np.sqrt(G_u*m/r)
    if physicalUnits:
        return vc*V_u
    return vc

def thinDiskVc(r,physicalUnits=None):
    MNx,MNy,MNz = getMiyamotoNagai(r,0,0)[1]
    a = np.sqrt(MNx**2+MNy**2+MNz**2)
    vc = np.sqrt(a*r)
    if physicalUnits:
        return vc*V_u
    return vc

def thickDiskVc(r,physicalUnits=None):
    MNx,MNy,MNz = getMiyamotoNagai(r,0,0)[2]
    a = np.sqrt(MNx**2+MNy**2+MNz**2)
    vc = np.sqrt(a*r)
    if physicalUnits:
        return vc*V_u
    return vc

def totalDiskVc(r,physicalUnits=None):
    MNx,MNy,MNz = getMiyamotoNagai(r,0,0)[0]
    a = np.sqrt(MNx**2+MNy**2+MNz**2)
    vc = np.sqrt(a*r)
    if physicalUnits:
        return vc*V_u
    return vc

## -------------------------------------------------------- INTEGRATOR ----------------------------------------------------------------------------------
def acceleration(x,y,z,physicalUnits=None):
    """Takes position [kpc] and returns acceleration at given radius [GADGET] or [km/s^2]"""
    r = np.sqrt(x**2+y**2+z**2)

    # spherically symmetric components
    M_r_NFW,M_r_H = getMassNFW_H(r) #mass within radius r for halo and bulge
    aSpherical = -G_u*(M_r_NFW+M_r_H)/r**2

    # disk
    MNx,MNy,MNz = getMiyamotoNagai(x,y,z)[0] #analytical acceleration due to disk

    ax = aSpherical*x/r+MNx #+adynNFW(np.array([x,y,z]),np.array([vx,vy,vz]),m,z)[0]
    ay = aSpherical*y/r+MNy #+adynNFW(np.array([x,y,z]),np.array([vx,vy,vz]),m,z)[1]
    az = aSpherical*z/r+MNz #+adynNFW(np.array([x,y,z]),np.array([vx,vy,vz]),m,z)[2]

    if physicalUnits:
        return ax*A_u,ay_*A_u,az*A_u
    return ax,ay,az

def leapfrog(time, timestep, ri, vi, physicalUnits=None):
    """
    Leapfrog integrator

    time: duration of orbit [years]
    timestep: delta time [years]
    ri: np.array([initialXPosition, initialYPosition, initialZPosition]) [kpc]
    vi: np.array([initialXVelocity, initialYVelocity, initialZVelocity]) [km/s]
    """
    time /= T_u
    timestep /= T_u

    # Initial Position
    x = ri[0]
    y = ri[1]
    z = ri[2]
    r = np.sqrt(x**2+y**2+z**2)

    #Arrays that will store positions over time, initialise with initial position
    xs = [x]
    ys = [y]
    zs = [z]
    rs=[r]

    # Initial Velocity (converted to GADGET)
    vx = vi[0]/V_u
    vy = vi[1]/V_u
    vz = vi[2]/V_u

    # Initial Acceleration
    ax,ay,az = acceleration(x,y,z)

    # Single half-step of Euler to obtain v_1/2
    vx = 0.5*ax*timestep+vx
    vy = 0.5*ay*timestep+vy
    vz = 0.5*az*timestep+vz

    t=0
    def boolT(t):
        if(time>0):
            return t<time
        else:
            return t>time
    while(boolT(t)):
        # Update position
        x = vx*timestep+x
        y = vy*timestep+y
        z = vz*timestep+z
        r = np.sqrt(x**2+y**2+z**2)

        # Update acceleration
        ax,ay,az = acceleration(x,y,z)

        # Update velocity
        vx = ax*timestep+vx
        vy = ay*timestep+vy
        vz = az*timestep+vz

        # Add new point to list
        xs.append(x)
        ys.append(y)
        zs.append(z)
        rs.append(r)

        t+=timestep

    if physicalUnits:
        return([xs,ys,zs,rs],[vx*V_u,vy*V_u,vz*V_u])

    return([xs,ys,zs,rs],[vx,vy,vz]) # position arrays, final velocity

def leapfrog_detailed(time, timestep, ri, vi, physicalUnits=None):
    """
    Leapfrog integrator

    time: duration of orbit [years]
    timestep: delta time [years]
    ri: np.array([initialXPosition, initialYPosition, initialZPosition]) [kpc]
    vi: np.array([initialXVelocity, initialYVelocity, initialZVelocity]) [km/s]
    """
    time /= T_u
    timestep /= T_u

    # Initial Position
    x = ri[0]
    y = ri[1]
    z = ri[2]
    r = np.sqrt(x**2+y**2+z**2)

    #Arrays that will store positions over time, initialise with initial position
    xs = [x]
    ys = [y]
    zs = [z]
    rs=[r]

    # Initial Velocity (converted to GADGET)
    vx = vi[0]/V_u
    vy = vi[1]/V_u
    vz = vi[2]/V_u

    vxs = [vx]
    vys = [vy]
    vzs = [vz]

    # Initial Acceleration
    ax,ay,az = acceleration(x,y,z)

    # Single half-step of Euler to obtain v_1/2
    vx = 0.5*ax*timestep+vx
    vy = 0.5*ay*timestep+vy
    vz = 0.5*az*timestep+vz

    t=0
    def boolT(t):
        if(time>0):
            return t<time
        else:
            return t>time
    while(boolT(t)):
        # Update position
        x = vx*timestep+x
        y = vy*timestep+y
        z = vz*timestep+z
        r = np.sqrt(x**2+y**2+z**2)

        # Update acceleration
        ax,ay,az = acceleration(x,y,z)

        # Update velocity
        vx = ax*timestep+vx
        vy = ay*timestep+vy
        vz = az*timestep+vz

        # Add new point to list
        xs.append(x); vxs.append(vx)
        ys.append(y); vys.append(vy)
        zs.append(z); vzs.append(vz)
        rs.append(r)

        t+=timestep

    if physicalUnits:
        return([xs,ys,zs,rs],[np.array(vxs)*V_u,np.array(vys)*V_u,np.array(vzs)*V_u])

    return([xs,ys,zs,rs],[vxs,vys,vzs]) # position arrays, final velocity

## ------------------------------------------------------------ USEFUL FUNCTIONS ----------------------------------------------------------------------

def initConditions(timeBack,timestep,p,v, physicalUnits=None):
    [xR,yR,zR,rsR],[vxi,vyi,vzi] = leapfrog(timeBack,timestep,np.array(p),np.array(v))
    if physicalUnits:
        return([xR[-1],yR[-1],zR[-1]],[vxi*V_u,vyi*V_u,vzi*V_u])
    return([xR[-1],yR[-1],zR[-1]],[vxi,vyi,vzi])

def getPeriApoTime(time,timestep,p,v):
    """Returns coordinates of apocentre closest to [t_approx] Gyr ago
    p: [x,y,z] coordinates in kpc
    v: [vx,vy,vz] coordinates in km/s
    time: approximate integrated time in years
    timestep: time step in years

    Returns pericentre [kpc], apocentre [kpc], and orbital period [Gyr]
    """
    apo = None
    peri = None
    period = None
    [x,y,z,r],[vx,vy,vz] = leapfrog(time,timestep,np.array(p),np.array(v),True)
    peaks = []
    for i in range(1,len(r)-1):
        if(r[i]<r[i-1] and r[i]<r[i+1]): # get pericentre closest to today
            peaks.append(i)
            if(peri==None):
                peri = r[i]
        elif(r[i]>r[i-1] and r[i]>r[i+1]): # get apocentre closest to today
            peaks.append(i)
            if(apo==None):
                apo = r[i]
    times = []
    if(len(peaks)>=2): #at least one apocentre and pericentre
        for i in range(1,len(peaks)):
            times.append(peaks[i]-peaks[i-1]) # difference between current peak and last peak
        period = abs(2*np.mean(times)*timestep/1e9)
    return (peri,apo,period)

def getInjectionCoords(p,v, t_approx, physicalUnits=None):
    """Returns coordinates of apocentre closest to [t_approx] Gyr ago
    p: [x,y,z] coordinates in kpc
    v: [vx,vy,vz] coordinates in km/s
    t_approx: approximate injection time in positive Gyr
    """
    t = -(t_approx*1.5)*1e9
    ts = -1e5
    apoTimes = np.array([])
    [x,y,z,r],[vx,vy,vz] = leapfrog(t,ts,np.array(p),np.array(v))
    for i in range(1,len(r)-1):
        if(r[i]>r[i-1] and r[i]>r[i+1]):
            apoTimes = np.append(apoTimes,i*ts)
    if(len(apoTimes)<1):
        apoTime = -(t_approx)*1e9 # If unbound, set injection to point [t_approx] Gyr ago
    else:
        apoTime = apoTimes[np.where(abs(apoTimes+(t_approx)*1e9)==min(abs(apoTimes+(t_approx)*1e9)))[0][0]] #Find apocentre closest to [t_approx] Gyr ago
    [x,y,z],[vx,vy,vz]=initConditions(apoTime,ts,np.array(p),np.array(v), physicalUnits)
    return [x,y,z],[vx,vy,vz],apoTime

def getRecentestPericentre(p,v, physicalUnits=None):
    """Returns coordinates of most recent apocentre
    p: [x,y,z] coordinates in kpc
    v: [vx,vy,vz] coordinates in km/s
    t_approx: approximate injection time in positive Gyr
    """
    t = -10*1e9
    ts = -1e5
    periTimes = np.array([])
    [x,y,z,r],[vx,vy,vz] = leapfrog(t,ts,np.array(p),np.array(v))
    for i in range(1,len(r)-1):
        if(r[i]<r[i-1] and r[i]<r[i+1]):
            periTimes = np.append(periTimes,i*ts)
    if(len(periTimes)<1):
        return [None,None,None],[None,None,None],None
    else:
        periTime = periTimes[np.argmax(periTimes)] #Find pericentre closest to [t_approx] Gyr ago
    [x,y,z],[vx,vy,vz]=initConditions(periTime,ts,np.array(p),np.array(v), physicalUnits)
    return [x,y,z],[vx,vy,vz],periTime
