#!/usr/bin/env python3
from math import pi
import numpy as np
from scipy.integrate import quad
import pandas as pd
from scipy.special import gammaincinv

import h5py
import sys

import warnings
warnings.filterwarnings("ignore")
# NB: when using interpolate, make sure X-array is increasing funct.


import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from matplotlib import rc
#rc('text', usetex=True)
plt.rcParams['font.size'] = 16
plt.rcParams['font.family'] = "serif"
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['figure.dpi'] = 150
plt.rcParams['figure.autolayout'] = True


G = 1
np.random.seed(667408)



print ("    * Reads model.hdf5. Model must be centred on 0,0,0. Units: G=1, Mtot=1 )")
filename = sys.argv[1]
f        = h5py.File(filename, 'r')
partmass = f["/PartType1/Masses"][0]
PartIDs  = np.copy(   f['/PartType1/ParticleIDs'][:]   )
x        = np.copy(   f['/PartType1/Coordinates'][:,0] )
y        = np.copy(   f['/PartType1/Coordinates'][:,1] )
z        = np.copy(   f['/PartType1/Coordinates'][:,2] )
vx       = np.copy(   f['/PartType1/Velocities'][:,0]  )
vy       = np.copy(   f['/PartType1/Velocities'][:,1]  )
vz       = np.copy(   f['/PartType1/Velocities'][:,2]  )
N        = len(x)
f.close()

print("partmass" , partmass)
print("N" , N)
print("x", x.shape)
print("v", vx.shape)
print("id", PartIDs.shape)

print("      *  Computing DM energies from model" )
#

r         = np.sqrt(x**2 + y**2 + z**2)


Mr        = partmass * np.arange(N)
#Mr = 1./float(N) * np.arange(N)


Ekin      = 0.5 * ( (vx)**2. +  (vy)**2. +  (vz)**2. )
#Ekin      = 0.5 * ( (vx)**2. +  (vy)**2. +  (vz)**2. )  /partmass * 1./float(N)


sortedr   = np.sort(r)
reverser  = sortedr[::-1]
SortIdx   = np.argsort(r)
UnsortIdx = np.argsort(SortIdx)
Phi       =  - G*  ( Mr / sortedr + np.cumsum(partmass / reverser )[::-1] )
Etot      = Phi[UnsortIdx] + Ekin


print("      *  Defining stellar distribution" )
#
def rhoS(r) :   return np.exp(-(r/rsS)**(1/n))
n  = 1
rsS = float(sys.argv[2]) # stellar scale radius in units of HDF5 file;
MS  = 4 * pi * quad( lambda r: r*r * rhoS(r) , 0., np.inf)[0] # non-normalized total stellar mass
MSr = np.vectorize(lambda x: 4*pi *  quad( lambda r: r*r * rhoS(r)/MS , 0., x  )[0]  )




print("      *  Binning DM model to compute DF" )
#
EN     = 100
DFNr   = 100
DFrmax = np.max(r)
DFrmin = np.min(r)
EPSREL = 1e-7

DFr    = np.logspace(np.log10(DFrmin),np.log10(DFrmax),num=DFNr)

NoutS =   1. - 4. * pi * quad( lambda r: r*r * rhoS(r)/MS , 0., DFrmax)[0]
NinS  =        4. * pi * quad( lambda r: r*r * rhoS(r)/MS , 0., DFrmin)[0]
print("     (!) Stars: there are %.2f per cent < DRrmin "%NinS )
print("     (!) and %.2f per cent > DFrmax "%NoutS )




print("      *  Building Psi and Nu arrays" )
psi  = - np.interp(DFr, sortedr, Phi)
nuDM, edges = np.histogram(r, bins=np.logspace(np.log10(DFrmin),np.log10(DFrmax),num=DFNr+1) )
nuDM        = nuDM * partmass  / ( 4./3 * pi * ((edges[1:])**3 - (edges[:-1])**3 )  )
nuS  =   rhoS(DFr) / MS


Mcum  = np.interp(DFr, sortedr, Mr)  # DM cumulative mass on array
MScum = MSr(DFr)               # stellar cumulative mass on array

plt.clf()
plt.plot(np.log10(DFr), np.log10(nuDM), color="red", label = "DM")
plt.plot(np.log10(DFr), np.log10(nuS),  color="blue", label = "stars")
plt.xlabel ("$\\log_{10} r$")
plt.ylabel ("$\\log_{10} \\nu_\\star$")
plt.legend()
plt.savefig("nu.png")


plt.clf()
plt.plot(np.log10(DFr), -psi, color="red", label = "DM")
plt.xlabel ("$\\log_{10} r$")
plt.ylabel ("$\\Phi(r)$")
plt.savefig("phi.png")


print("      *  Calculating gradients" )
# starts are collisionless tracers, i.e. psi doesnt give a f about nuS
dndpDM   = np.gradient(nuDM,   psi)
d2nd2pDM = np.gradient(dndpDM, psi)
dndpS   = np.gradient(nuS,   psi)
d2nd2pS = np.gradient(dndpS, psi)


less_than_one = np.where (Mcum < 1.0 - float(1./N) )[0]
less_than_one_idx = np.argmax(Mcum >= 1.0)


print("      *  Evaluating DFs (stars + DM)" )
fS  = np.vectorize( lambda e: 1./(np.sqrt(8)*pi*pi) * (quad( lambda p:  np.interp(p, psi[::-1], d2nd2pS[::-1]) / np.sqrt(e-p) , 0., e,  epsrel=EPSREL)[0]  ) ) # + np.interp(0., psi, dndp) / np.sqrt(e)   == 0 due to B.C.
fDM = np.vectorize( lambda e: 1./(np.sqrt(8)*pi*pi) * (quad( lambda p:  np.interp(p, psi[::-1], d2nd2pDM[::-1]) / np.sqrt(e-p) , 0., e,  epsrel=EPSREL)[0]  ) ) # + np.interp(0., psi, dndp) / np.sqrt(e)   == 0 due to B.C.


maxE = psi[0]
minE = maxE/float(EN)
E = np.linspace(minE,maxE,num=EN)

print("      *    ...stars..." )
DFS = fS(E)
if np.any (DFS < 0) :
  print("      *  Exit. star DF < 0, see df.dat. NO DM DF computed yet." )
  np.savetxt("df.dat", np.column_stack(( E, DFS)))
  sys.exit(0)


print("      *    ...DM..." )
DFDM = fDM(E)
if np.any (DFDM < 0) :
  print("      *  Exit. DM DF < 0, see df.dat" )
  np.savetxt("df.dat", np.column_stack(( E, DFS, DFDM)))
  #sys.exit(0)

print("      *  star + dm DF >= 0, all good" )


pd.DataFrame({"E": E, 'DFDM': DFDM, "DFs": DFS}).to_csv('DFDM.csv')

print("      *  Computing probabilities" )
probs = np.interp(-Etot, E, DFS)/ np.interp(-Etot, E, DFDM)
probs = probs/np.sum(probs[~np.isnan(probs)])

print(probs)
print("      *  Writing probability file stars.npy -- to be used with specific model.dat" )
np.save("stars.npy", np.column_stack( (PartIDs, probs) ) )


print("      *  Sanity check: density plot using new probabilities" )

nuS_Nbody, edges = np.histogram(r, weights=probs, bins=np.logspace(np.log10(DFrmin),np.log10(DFrmax),num=DFNr+1) )
nuS_Nbody        = nuS_Nbody  / ( 4./3 * pi * ((edges[1:])**3 - (edges[:-1])**3 )  )

plt.clf()
plt.plot(np.log10(DFr), np.log10(nuS), color="blue", linewidth=2, linestyle="dashed", label = "stars (analyt.)")
plt.plot(np.log10(DFr), np.log10(nuS_Nbody),  color="fuchsia", label = "stars (N-body)")
plt.ylim(np.max(np.log10(nuS))-3,0)
plt.xlabel ("$\\log_{10} r$")
plt.ylabel ("$\\log_{10} \\nu_\\star$")
plt.legend()
plt.savefig("nu_probs.png")

print("      *  All done :o)" )
