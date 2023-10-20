"""
Created on Mon Oct 18 17:37:17 2021

@author: asya

Samples observable parameters for an object and returns
orbits which lie at the desired percentiles of pericentre

"""

# ------------------------------ Imports --------------------------------------------------
import numpy as np
import pickle
import multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib
#from lib.util import getCoords
from integrator import getPeriApoTime,getInjectionCoords,leapfrog
from astropy.coordinates import SkyCoord,Galactocentric
from astropy import units as u


################################# INPUTS ##################################################
nCores = 1                                                 # number of cores used for multiprocessing
percentiles = [16,84]                                       # desired percentiles of pericentre
observables = {'ra': 177.3, 'dec': -18.4, 'D': 117.5,       # observables of object
               'sig_D': 1.1, 'pmra': -0.246,'sig_pmra': 0.052,
               'pmdec': -0.227, 'sig_pmdec': 0.026, 'rv': 87.5, 'sig_rv': 0.4}
N_coords=100                                              # number of samples to draw from data
peri_bins = np.linspace(0,60,50)
apo_bins = np.linspace(120,140,50)
title = 'Crater II'

###########################################################################################

# matplotlib.use('Agg')
# -------------------------------- Other Global Parameters --------------------------------

# Solar Parameters from Schonrich 2010
v_sun = [11.1,240.3+12.24, 7.25] * (u.km / u.s)
d_sun = 8.29 *u.kpc
gc_frame = Galactocentric(galcen_distance=d_sun)
# gc_frame = Galactocentric(galcen_distance=d_sun,galcen_v_sun=v_sun,z_sun=0*u.pc)

# -------------------------------- Plotting Settings --------------------------------------
#plt.style.use('mnras')

plt.rcParams['figure.figsize'] = [3.15, 3.15]
plt.rcParams['figure.facecolor'] = 'w'
plt.rcParams['figure.dpi'] = 150

def set_aspect(ax,log_x,log_y):
    if(log_x and log_y):
        ax.set_aspect(np.log10(ax.get_xlim()[1]/ax.get_xlim()[0])/np.log10(ax.get_ylim()[1]/ax.get_ylim()[0]))
    elif(log_x):
        ax.set_aspect(np.log10(ax.get_xlim()[1]/ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0]))
    elif(log_y):
        ax.set_aspect((ax.get_xlim()[1]-ax.get_xlim()[0])/np.log10(ax.get_ylim()[1]/ax.get_ylim()[0]))
    else:
        ax.set_aspect(np.diff(ax.get_xlim())/np.diff(ax.get_ylim()))

# ------------------------------ Functions ------------------------------------------------
def drawCoords(ra,dec,D,sig_D,pmra,sig_pmra,pmdec,sig_pmdec,rv,sig_rv,N):
    np.random.seed(1)
    dSet = np.random.normal(D,sig_D,N)
    np.random.seed(1)
    pmraSet = np.random.normal(pmra,sig_pmra,N)
    np.random.seed(1)
    pmdecSet = np.random.normal(pmdec,sig_pmdec,N)
    np.random.seed(1)
    rvSet = np.random.normal(rv,sig_rv,N)

    coordSet = np.array([np.full(N,ra),np.full(N,dec),dSet,pmraSet,pmdecSet,rvSet]).T

    return coordSet


def getCoords(ra, dec, D, pmra, pmdec, rv):
    sc = SkyCoord(ra=ra*u.degree, 
            dec=dec*u.degree, 
            distance=D*u.kpc, 
            radial_velocity=rv * u.km/u.s, 
            pm_ra_cosdec=pmra * u.mas/u.yr, 
            pm_dec=pmdec * u.mas / u.yr)
    x, y, z = sc.cartesian.get_xyz() / u.kpc
    vx = sc.cartesian.differentials.get("s").d_x / (u.km/u.s)
    vy = sc.cartesian.differentials.get("s").d_y / (u.km/u.s)
    vz = sc.cartesian.differentials.get("s").d_z / (u.km/u.s)

    return [(x, y, z), (vx, vy, vz)]



def helper(c):
        [(x,y,z),(vx,vy,vz)] = getCoords(c[0],c[1],c[2],c[3],c[4],c[5])
        peri,apo,t = getPeriApoTime(-4e9,-1e6,[x,y,z],[vx,vy,vz])
        return (peri,apo)

def periApoHist(coords):
    periSet = np.array([])
    apoSet = np.array([])
    coordSet = []

    print('multiprocessing...',end='\r')
    pool = mp.Pool(nCores)
    results = np.array(pool.map(helper, coords))
    pool.close()
    print('done multiprocessing')

    # select results where pericentre and apocentre exist and are nonzero
    ok_idx = (results[:,0]!=None) & (results[:,1]!=None)
    ok_idx2 = (results[:,0][ok_idx]>0) & (results[:,1][ok_idx]>0)
    periSet = results[:,0][ok_idx][ok_idx2]
    apoSet = results[:,1][ok_idx][ok_idx2]
    coordSet = coords[ok_idx][ok_idx2]

    return periSet,apoSet,coordSet

def dndr(rSet,bs):
    N, edges = np.histogram(rSet,bins=bs,density=True)
    dndr = N/(edges[1:]-edges[:-1])
    return edges[1:],dndr#(edges[1:]+edges[:-1])/2,dndr

def get_percentile(rSet,percentile):
    #percentile can be for ex. 16 or 84
    sortedIdx = np.argsort(rSet)
    sortedSet = rSet[sortedIdx]
    nSet = len(rSet)
    idx_insorted = round(percentile/100*nSet)
    p = sortedSet[idx_insorted]
    idx = sortedIdx[idx_insorted]
    return p,idx

# ------------------------------ Main -------------------------------------------
coordinates_present = {}; coordinates_past = {}
[(x_mean,y_mean,z_mean),(vx_mean,vy_mean,vz_mean)] = getCoords(observables['ra'],observables['dec'],
                                                               observables['D'],observables['pmra'],
                                                               observables['pmdec'],observables['rv'])
coordinates_present['mean'] = [x_mean,y_mean,z_mean,vx_mean,vy_mean,vz_mean]

peri_mean, apo_mean = helper([observables['ra'],observables['dec'],
                              observables['D'],observables['pmra'],
                              observables['pmdec'],observables['rv']])

coords = drawCoords(observables['ra'],observables['dec'],observables['D'],observables['sig_D'],
                    observables['pmra'],observables['sig_pmra'],observables['pmdec'],
                    observables['sig_pmdec'],observables['rv'],observables['sig_rv'],N_coords)
periSet,apoSet,coordSet = periApoHist(coords)


peris = []; apos = []
for per in percentiles:
    p,idx = get_percentile(periSet,per) #pericentre at this percentile
    peris.append(p); apos.append(apoSet[idx])
    c = coords[idx]
    [(x,y,z),(vx,vy,vz)] = getCoords(c[0],c[1],c[2],c[3],c[4],c[5])
    coordinates_present[str(per)] = [x,y,z,vx,vy,vz]

apoTimes = {}
for k,v in coordinates_present.items():
    [x,y,z,vx,vy,vz] = v
    [x2,y2,z2],[vx2,vy2,vz2],apoTime = getInjectionCoords([x,y,z],[vx,vy,vz], 10)
    coordinates_past[k] = [x2,y2,z2,vx2,vy2,vz2]
    apoTimes[k] = apoTime


f0 = open('coordinates-past.pkl','wb')
pickle.dump(coordinates_past, f0)
f0.close()

f1 = open('coordinates-present.pkl','wb')
pickle.dump(coordinates_present, f1)
f1.close()

f2 = open('orbits.out','w')
f2.write('Orbital x,y,z,vx,vy,vz [kpc and km/s] are saved in coordinates.pkl\n'+
         'Using gc_frame as defined in orbits.py, these can be converted back to sky coordinates like this:\n'+
         "c_sample = SkyCoord(x=dict['16'][0]*u.kpc, y=dict['16'][1]*u.kpc,z=dict['16'][2]*u.kpc,v_x=dict['16'][3]*u.km/u.s, v_y=dict['16'][4]*u.km/u.s, v_z=dict['16'][5]*u.km/u.s,frame=gc_frame)\n"+
         "c_sample.transform_to(coord.ICRS)\n\n"+
         'The percentile pericentres are: {} [kpc]\n'.format(peris)+
         'The mean pericentre is: {} [kpc]\n\n'.format(peri_mean)+
         'The percentile apocentres are: {} [kpc]\n'.format(apos)+
         'The mean apocentre is: {} [kpc]\n'.format(apo_mean))
f2.close()

#%% figure 1
fig, (ax0,ax1) = plt.subplots(1,2,figsize=(8,4))

peri_dist = dndr(periSet,peri_bins)
ax0.step(peri_dist[0],peri_dist[1],c='b')
ax0.fill_between(peri_dist[0],peri_dist[1],facecolor='b',alpha=0.2,step='pre')
for per in [16,50,84]:
    p,idx = get_percentile(periSet,per)
    idx2 = np.argmin(abs(peri_dist[0]-p))
    ax0.fill_between(peri_dist[0][idx2:], peri_dist[1][idx2:],facecolor='b',alpha=0.3,step='pre')
ax0.set_xlabel('$r_\mathrm{peri}$ [kpc]')

apo_dist = dndr(apoSet,apo_bins)
ax1.step(apo_dist[0],apo_dist[1],c='b')
ax1.fill_between(apo_dist[0],apo_dist[1],facecolor='b',alpha=0.2,step='pre')
for per in [16,50,84]:
    a,idx = get_percentile(apoSet,per)
    idx2 = np.argmin(abs(apo_dist[0]-a))
    ax1.fill_between(apo_dist[0][idx2:], apo_dist[1][idx2:],facecolor='b',alpha=0.3,step='pre')
ax1.set_xlabel('$r_\mathrm{apo}$ [kpc]')

for ax in [ax0,ax1]:
    ax.set_ylabel('d$N$ / d$r$')
    set_aspect(ax,False,False)
    ax.tick_params(axis='y',which='both',labelleft=False,length=0.)
    ax.tick_params(axis='x',which='minor',length=2,labeltop=False); ax.tick_params(axis='x',which='major',length=3,labeltop=False)

ax_t0 = ax0.secondary_xaxis('top'); ax_t0.tick_params(axis='x', direction='inout')
ax_t0.set_xticks([peri_mean]+peris); ax_t0.set_xticklabels(['orb{}'.format(i) for i in range(len(peris)+1)],rotation=45)

ax_t1 = ax1.secondary_xaxis('top'); ax_t1.tick_params(axis='x', direction='inout')
ax_t1.set_xticks([apo_mean]+apos); ax_t1.set_xticklabels(['orb{}'.format(i) for i in range(len(apos)+1)],rotation=45)

ax0.margins(x=0,y=0)
ax1.margins(x=0,y=0)
fig.suptitle(title)
plt.savefig('histogram',bbox_inches='tight',pad_inches=0)

#%% figure 2
fig, axes = plt.subplots(2,1+len(percentiles),figsize=(4*(1+len(percentiles)),8),sharex=True)

i = 0
for k,v in coordinates_present.items():
    [x,y,z,vx,vy,vz] = v
    apoTime = apoTimes[k]
    [x_list,y_list,z_list,r_list],[vx_tmp,vy_tmp,vz_tmp] = leapfrog(apoTime,-1e6,[x,y,z],[vx,vy,vz])
    axes[0][i].plot(x_list,y_list,c='b')
    axes[1][i].plot(x_list,z_list,c='b')
    axes[0][i].set_title('orb{}'.format(i))
    i+=1

maxlim = max(max([ax.get_ylim()[1] for ax in fig.axes]),axes[0][0].get_xlim()[1])
minlim = min(min([ax.get_ylim()[0] for ax in fig.axes]),axes[0][0].get_xlim()[0])

for ax in fig.axes:
    ax.set(xlim=(minlim,maxlim),ylim=(minlim,maxlim))
    set_aspect(ax,False,False)
    ax.set(xlabel='x [kpc]')
axes[0][0].set(ylabel='y [kpc]'); axes[1][0].set(ylabel='z [kpc]')


plt.savefig('orbits.pdf',bbox_inches='tight',pad_inches=0.2)

