In this section, we discuss our simulation methods from orbital estimation and initial conditions to post-processing and software.

We use Agama [@agama] for potential specification and initial conditions, and Gadget-4 [@gadget4] for orbit integration and dark matter simulations (with a custom patch to integrate with Agama). 

# Milky Way Potential 

@fig:v_circ_potential plots the circular velocity profiles of each component of our fiducially potential. We adopt the potential described in @EP2020, an analytic approximation of @mcmillan2011. The potential consisting of a stellar Bulge, a thin and thick disk, and a dark matter NFW halo. Especially in the regime where Sculptor and Ursa Minor orbit, the dark matter halo is the dominant component of the potential. The details of the stellar disk and bulge are less important. 

The galactic bulge is described by a @hernquist1990 potential. The potential is

$$
\Phi(r) = - \frac{GM}{r + a}
$$
where  $a=1.3\,{\rm kpc}$ is the scale radius and $M=2.1 \times 10^{10}\,\Mo$ is the total mass.

The thin and thick disks are represented with the @miyamoto+nagai1975  cylindrical potential:
$$
\Phi(R, z) = \frac{-GM}{\left(R^2 + \left[a + \sqrt{z^2 + b^2}\right]^{2}\right)^{1/2}}
$$

where $a$ is the disc radial scale length, $b$ is the scale height, and $M$ is the total mass of the disk. For the thin disk,  $a=3.944\,$kpc, $b=0.311\,$kpc, $M=3.944\times10^{10}\,$M$_\odot$. For the thick disk, $a=4.4\,$kpc, $b=0.92\,$kpc, and $M=2\times10^{10}\,$M$_\odot$.

The halo is a NFW dark matter halo (REF) with $M_s=79.5\times10^{10}\,$M$_\odot$.  and $r_s = 20.2\,$kpc.

Variations to the potential of the inner disk (exclusion of a bar) should minimally affect our results as no orbit we consider reaches less than ~15 kpc of the MW centre. We exclude the mass evolution of the halo from this analysis. Over $10\,$Gyr, this would be fairly significant (factor of $\sim 2$in MW mass, REF) but since we want to determine the upper limit of tidal effects, it is safe to neglect this. 

Finally, in chapter REF, we consider the influence of the large Milky Cloud on the orbits and evolution of Scl. We adopt the @vasiliev2024 multipole approximation of an N-body simulation of the LMC and MW. Their initial conditions are

- MW halo:
- MW bulge (static):
- MW disk (static):
- LMC halo:

![Circular velocity of potential](figures/v_circ_potential.png){#fig:v_circ_potential} 

Figure: Circular velocity profile of @EP2020 potential. We also show the MW and LMC potential from @vasiliev2024 for their L3M11 model.

# Orbital Estimation

To estimate the orbits of a dwarf galaxy, we perform a Monte Carlo sampling of the present-day observables. We integrate each sampled observable back in time for 10 Gyr in Gadget as massless point particles outputting positions every 5Myr (with otherwise similar parameters to n-body runs below.)

To convert from Gaia to galactocentric coordinates, we use the @astropy v4 galactocentric frame. This frame assumes the galactic centre (Sagittarius A*) appears:

- J2000 'ra' = 17h45m37.224s, 'dec' = -28°56'10.23", $\mu_{\alpha*}=-3.151\pm0.018$ mas/yr, $\mu_\delta=-5.547\pm0.026$ from the appendix and Table 2 of @reid+brunthaler2004.
- distance = $8.122\pm0.033\,$kpc, solar radial velocity = $11 + 1.9 \pm 3$ km/s.  from  @gravitycollaboration+2018.

Finally, adding that the sun is $20.8\pm0.3\,$pc above the disk from bennett+bovy2019,  and using  the procedure described in @drimmel+poggio2018, the solar velocity relative to the galactic rest frame is 

- 'v_sun' = $[-12.9 \pm 3.0, 245.6 \pm 1.4, 7.78 \pm 0.08]$ km/s



We select the initial conditions from the location of the first apocentre of a selected orbit.

While we assume that the galaxy can be described as a point particle, not subject to (internal) dynamical friction, this approximation works well if the galaxy is not strongly affected by tides (like Scl). Analytic approximations of dynamical friction tend to be inadequate to more accurately model the N-body trajectory.

# N-Body Modelling

## Initial conditions

We use Agama [@agama] to generate initial conditions. We initially assume galaxies are described by an NFW dark matter potential (REF) and the stars are merely collisionless tracers embedded in this potential (added on in post-processing). The density is a cubic-exponentially truncated with a profile
$$
\rho_{\rm tNFW} = e^{-(r/r_t)^3}\ \rho_{\rm NFW}(t)
$$ {#eq:trunc_nfw}
where we adopt $r_t = 20 r_s$ or approximately $r_{200}$ for our Sculptor-like fiducial halo. (Using $r_{200}$ depends on the chosen scale of the halo, we chose to be scale free here.) This truncation causes the density to drop quickly past $20 r_s$ but retails > 90 $\%$ agreement with an NFW within $18r_s$. This choice is purely for numerical convenience, however it is unlikely that the outer density profile of loosely bound particles past $20$ kpc affects the evolution of the halo. A more detailed cosmological consideration would be ideal.

## Isolation runs and simulation parameters

To ensure that the initial conditionss of the simulation are dynamically relaxed and well-converged, we run the simulation in isolation (no external potential) for 5 Gyr (or about 3 times the crossing timescale at the virial radius). Our fiducial isolation halo uses $r_s=2.76$ kpc and $M_s = 0.29 \times 10^{10}$ Msun, but can be easily rescaled for any length or mass scale. 

We adopt a softening length of 
$$
h_{\rm grav} = 0.014 \left(\frac{r_s}{2.76\,{\rm kpc}}\right)\left(\frac{N}{10^7}\right)^{-1/2}.
$$
See appendix REF for a discussion of this choice, which is similar to the @power+2003  suggested softening. 

We use the relative tree opening criterion with the accuracy parameter set to 0.005, and adaptive time stepping with integration accuracy set to 0.01. 

### Numerical convergence

As discussed in @power+2003, the region where density is converged is related to the region which becomes collisionally relaxed over the time of the universe (so about 5-10 Gyr). (i.e. where the collisionless assumption breaks down). The relaxation timescale is given by $t_{\rm relax} = t_{\rm circ} N(r) / 8\ln N(r)$ or otherwise, 
$$
t_{\rm relax}(r) := t_{\rm circ}(r) \frac{N(r)}{8\,\ln N(r)}
= {t_{\rm circ}(r_{200})} \frac{\sqrt{200}}{8} \frac{N(r)}{\ln N(r)} \left(\frac{\bar \rho (r)}{\rho_{\textrm crit}}\right)^{-1/2}
$$
This works out to be about 6-10 times  (increasing with particle number) our adopted softening length for NFW halos given our assumptions. As such, at full resolution, we can only trust density profiles down to $\sim10\epsilon$, just enough to resolve stellar density profiles. 

Note that by decreasing the truncation radius, the effective particle number is improved, so future experiments could be less generous with this parameter to improve computational performance (order 30% better resolution).

![Numerical halo convergence](figures/iso_converg_num.png){#fig:numerical_convergance}

Figure: Numerical convergence test for circular velocity as a function of log radius for simulations with different total numbers of particles in isolation. Residuals in lower panel are relative to NFW. The initial conditions are dotted and the converged radius is marked by arrows (REF). 

## Orbital runs

To perform the simulations of a given galaxy in a given potential, we centre the isolation run's final snapshot and place the dwarf galaxy in the specified orbit in the given potential. We typically run the simulation for 10 Gyr, which allows us to orbit slightly past the expected initial conditions. 

# Post Processing

## Centring

Shrinking spheres centres inspired by @power+2003

- Recursively shrink radius by 0.975 quantile and recalculating centroid until radius is less than ~1kpc or fewer than 0.1% of particles remain.
- Remove bound particles (using instantanious potential).
- Use previous snapshot centre and acceleration to predict new centre for next snapshot and only include particles included in previous centring timestep.

The statistical centring uncertainty for the 1e7 particle isolation run is of order 0.003 kpc, however oscillations in the centre are of order 0.03 kpc. This is about three times the softening length but is less than the numerically converged radius scale. 

## Velocity profiles



- Circular velocity is computed assuming spherical symmetry and only shown for every 200th particle (ranked from the centre outwards)

$$
v_{\rm circ}(r) = \sqrt{\frac{GM(r)}{r}}
$$



- $v_{\rm circ\ max}$ and $r_{\rm circ,\ max}$ is found by least-squares fitting the NFW functional velocity form to the points of the velocity profile that have the 10% highest velocities. This is not a necessarily good fit, especially as the halo becomes stripped, but accurate enough to find a reliable maximum. 

## Stellar Probabilities

We assign stellar probabilities (assuming spherical, isotropic symmetry) via Eddington inversion. 

Let $\Psi$ be the potential (normalized to vanish at infinity) and  ${\cal E}$ is the binding energy ${\cal E} = \Psi - 1/2 v^2$. If we know $f({\cal E})$, the distribution function (phase-space density in energy), then we assign the stellar weight for a given particle with energy ${\cal E}$ is 
$$
P_\star({\cal E}) = \frac{f_\star({\cal E})}{f_{\rm halo}({\cal E})}.
$$
While $f({\cal E})$ is a phase-space density, the differential energy distribution includes an additional $g({\cal E})$ occupation term (BTXXX).

We use Eddington inversion to find the distribution function, (eq. 4-140b in BT87)
$$
f({\cal E}) = \frac{1}{\sqrt{8}\, \pi^2}\left( \int_0^{\cal E} \frac{d^2\rho}{d\Psi^2} \frac{1}{\sqrt{{\cal E} - \Psi}}\ d\Psi + \frac{1}{\sqrt{\cal E}} \left(\frac{d\rho}{d\Psi}\right)_{\Psi=0} \right).
$$

In practice the right, boundary term is zero as $\Psi \to 0$ as $r\to\infty$, and if $\rho \propto r^{-n}$ at large $r$ and $\Psi \sim r^{-1}$ then $d\rho / d\Psi \sim r^{-n+1}$ which goes to zero provided that $n > 1$.



$\Psi$ is taken as known from the underlying assumed analytic dark matter potential. $\rho_\star$ is pre-specified or calculated via Abel integration / deprojection



We consider both a Plummer and 2-dimensional exponential (Exp2D) stellar models:

### Exp2D

Perhaps the simplest in form, the 2D exponential is commonly used to describe dwarf galaxy's density profiles as well as stellar disks. The stellar profile is given by 
$$
\Sigma_\star(R) = A e^{-R / R_s}
$$
where $A$ is the normalization and $R_s$ is the exponential scale radius. The 3D deprojected profile is found through Abel inversion [e.g. @errani+2024]:
$$
\rho_\star (r) =- \frac{1}{\pi}\int_r^\infty \frac{d\Sigma}{dR} \frac{1}{\sqrt{R^2 - r^2}} dR  = \frac{\Sigma_0}{\pi R_s^2}\,K_0(r/R_s)
$$
where $K_0$ is the 0th order modified Bessel function of the second type. 

### Plummer

The Plummer density profile is generally considered to be a reasonable empirical fit to Globular Clusters, where it was first proposed [@sadfjk]. Plummer profiles are also sometimes used to describe dwarf spheroidal galaxies REFS and galactic bulges REF. A Plummer profile is defined by
$$
\Sigma(R) = \frac{M}{4\pi R_s^2} \left(1 + \frac{R^2}{R_s^2}\right)^{-2} ,
$$
where $M$ is the total mass and $R_s$ is the characteristic scale radius. The 3D density is 
$$
\rho(r) = \frac{3M}{4\pi\,R_s^3} \left(1 + \frac{r^2}{R_s^2}\right)^{-5/2}.
$$



## Stellar density profiles and velocity dispersion

- Stellar velocity dispersion is calculated for all stars within 1kpc (3D) of the centre by assuming $\sigma_{\rm los} = \sqrt{\sigma_x^2 + \sigma_y^2 + \sigma_z^2}/ \sqrt 3$, i.e. isotropic velocity dispersion
- We calculate density profiles similar to stars and assume a constant bin width in $\log R$ of $2 {\rm IQR} / \sqrt[3]{N}$ (Freedman-Diaconis prescription). $R$ may be either the cylendrical radius in $x-y$ or the on-sky project $\xi, \eta$ tangent plane coordinates. 
- Because our dwarfs are assumed to be spherical/isotropic, we retain this assumption when calculating predicted density profiles. 







# Appendix: Numerical convergence and parameters

Here, we describe some convergence tests to ensure our methods and results are minimally impacted by numerical limitations and assumptions. See @power+2003 for a detailed discussion of various assumptions and parameters used in N-body simulations. 

### Softening

@power+2003 suggest the empirical rule that the ideal softening (balancing integration time and only compromising resolution in collisional regime) is 
$$
h_{grav} = 4 \frac{R_{200}}{\sqrt{N_{200}}}
$$

For our isolation halo ($M_s=2.7$, $r_s=2.76$) and with $10^7$ particles, this works out to be $0.044\,{\rm kpc}$.We adpoted the slightly smaller softening which was reduced by a factor of $\sqrt{10}$ which appears to improve agreement slightly in the innermost regions. However, this choice likely unnecessarily increases computation time for the relative gain in accuracy.

![Softening convergence](/Users/daniel/thesis/figures/iso_converg_softening.png){#fig:softening_convergence}

### Timestepping and force accuracy

In general, we use adaptive timestepping and relative opening criteria for gravitational force computations. To verify that these choices and associated accuracy parameters minimally impact convergence or speed, we show a few more isolation runs (using only 1e5 particles)

- constant timestep (...), approximantly half of minimum timestep with adaptive timestepping
- geometric opening, with $\theta = 0.5$.
- strict integration accuracy, (facc = ....)

### Alternative methods

- FMM
- PMM-tree
- Gadget2
- etc.

![Isolation method convergence](/Users/daniel/thesis/figures/iso_converg_methods.png){#fig:methods_convergence}

### Fiducial Parameters

Note that we use code units which assume that $G=1$ for convenience and numerical stability. The conversion between code units to physical units is (for our convention):

- 1 length = 1 kpc
- 1 mass unit = $10^{10}$ Msun
- 1 velocity unit = 207.4 km/s
- 1 time unit = 4.715 Myr

Most parameters below are not too relevant or have been discussed or are merely dealing with cpu and IO details. The changes between simulation runs primarily affect the integration time, output frequency, and softening. Otherwise, we leave all other parameters fixed. 

 ```
#======IO parameters======

#---Filenames
InitCondFile                initial
OutputDir                   ./out
SnapshotFileBase            snapshot
OutputListFilename          outputs.txt

#---File formats 
ICFormat                    3       # use HDF5
SnapFormat                  3 

#---Mem & CPU limits
TimeLimitCPU                86400
CpuTimeBetRestartFile       7200
MaxMemSize                  2400

#---Time
TimeBegin                   0
TimeMax	                   2120      # 10 Gyr

#---Output frequency
OutputListOn                0
TimeBetSnapshot             10
TimeOfFirstSnapshot         0 
TimeBetStatistics           10
NumFilesPerSnapshot         1
MaxFilesWithConcurrentIO    1 


#=======Gravity======


#---Timestep accuracy
ErrTolIntAccuracy           0.01
CourantFac                  0.1     # ignored; for SPH
MaxSizeTimestep             0.5
MinSizeTimestep             0.0 

#---Tree algorithm
TypeOfOpeningCriterion      1       # Relative
ErrTolTheta                 0.5     # mostly used for Barnes-Hut
ErrTolThetaMax              1.0     # (used only for relative)
ErrTolForceAcc              0.005   # (used only for relative)

#---Domain decomposition: should only affect performance
TopNodeFactor                       3.0
ActivePartFracForNewDomainDecomp    0.02

#---Gravitational Softening
SofteningComovingClass0      0.044  # HALO dependent
SofteningMaxPhysClass0       0      # ignored; for cosmological
SofteningClassOfPartType0    0
SofteningClassOfPartType1    0


#=======Miscellanius=======

# probably do not need to change the options below 

#---Unit System
UnitLength_in_cm            1 
UnitMass_in_g               1
UnitVelocity_in_cm_per_s    1 
GravityConstantInternal     1

#---Cosmological Parameters 
ComovingIntegrationOn	   0 # no cosmology
Omega0	                   0
OmegaLambda                 0 
OmegaBaryon                 0
HubbleParam                 1
Hubble                      100
BoxSize                     0

#---SPH
ArtBulkViscConst             0.8
MinEgySpec                   0
InitGasTemp                  100

#---Initial density estimate (SPH)
DesNumNgb                   64
MaxNumNgbDeviation          1
 ```



