# Simulation Methods

In this section, we discuss our simulation methods from orbital estimation to initial conditions to simulation parameters and software.

We use Agama [@agama] for potential specification and initial conditions, and Gadget-4 [@gadget4] for orbit integration and dark matter simulations (with a custom patch to integrate with Agama). 

## Potential 

Our fiducial potential is as described in @EP2020. This potential is an analytic approximation of @mcmillan2011, consisting of a stellar Bulge, disk, and a dark matter NFW halo. @fig:v_circ_potential plots the circular velocity profiles. 

The @hernquist1990 density profile for the galactic bulge is parameterized in terms of a characteristic mass, $M$, and radius $a$. 

$$
\rho(r) = \frac{M}{2\pi} \frac{a}{r} \frac{1}{(r+a)^3}
$$
Hernquist potential with $r_s=1.3$, mass = 2.1. 

- Thin disk. miyamoto+nagai disk with parameters $a=3.944$, $b=0.311$, $M=3.944$,
- Thick disk. miyamoto+nagai disk with parameters $a=4.4$, $b=0.92$, and $M=2$

$$
\Phi(R, z) = \frac{-GM}{\left(R^2 + \left[a + \sqrt{z^2 + b^2}\right]^{2}\right)^{1/2}}
$$



- NFW Dark matter halo. $M_s=79.5$  (double check). and $r_s = 20.2$.

Variations to the potential of the inner disk (exclusion of a bar) should minimally affect our results as no orbit we consider reaches less than ~15 kpc of the MW centre. We exclude the mass evolution of the halo from this analysis. Over 10 Gyr, this would be fairly significant (~2x in MW mass) but since we want to determine the upper limit of tidal effects, it is safe to neglect this. 

![Circular velocity of potential](figures/v_circ_potential.png){#fig:v_circ_potential} 

**Figure:** Circular velocity profile of @EP2020 potential. We also show the MW and LMC potential from @vasiliev2024 for their L3M11 model.

## Orbital Estimation

To estimate the orbits of a dwarf galaxy, we perform a Monte Carlo sampling of the present-day observables. We integrate each sampled observable back in time for 10 Gyr in Gadget as massless point particles outputing every 0.5Myr (with otherwise similar parameters to n-body runs below.)

To convert from Gaia to galactocentric coordinates, we use the astropy v4 galactocentric frame. This frame assumes an 

- J2000 'ra' = 17h45m37.224s, 'dec' = -28°56'10.23" (appendix), 'pmracosdec=-3.151±0.018 mas/yr', 'pmdec=-5.547±0.026' (Table 2). from @reid+brunthaler2004
- distance = 8.122 ± 0.033 kpc, solar radial velocity = 11 + 1.9 ± 3 km/s.  from  @gravitycollaboration+2018.
- 'z_sun' = 20.8±0.3pc from bennett+bovy2019
- 'v_sun' is slightly more complicated, relying on parameters from above and using  the procedure described in @drimmel+poggio2018. This results in 
  - 'v_sun' = [-12.9 ± 3.0, 245.6 ± 1.4, 7.78 ± 0.08] km/s

## Initial conditions

We use Agama [@agama] to generate initial conditions. We initially assume galaxies are described by an NFW dark matter potential (REF) and the stars are merely collisionless tracers embedded in this potential. We adopt a truncation radius of 100 times the scale radius (likely larger than necessary, but should not affect results.)

### Stellar Probabilities

We assign stellar probabilities via Eddington inversion using the following method:

- $\Psi$ is taken as known from the underlying assumed analytic dark matter potential.
- $\rho_\star$ is pre-specified or calculated via abel integration / deprojection
- $f(\epsilon)$ (the energy distribution function) is determined by Edington inversion

To find the distribution function, we use Eddington inversion (eq. 4-140b in BT87)
$$
f({\cal E}) = \frac{1}{\sqrt{8}\, \pi^2}\left( \int_0^{\cal E} \frac{d^2\rho}{d\Psi^2} \frac{1}{\sqrt{{\cal E} - \Psi}}\ d\Psi + \frac{1}{\sqrt{\cal E}} \left(\frac{d\rho}{d\Psi}\right)_{\Psi=0} \right)
$$


where ${\cal E}$ is the binding energy, and $\Psi$ is the potential normalized to go to zero at infinity. In practice the right, boundary term is zero as $\Psi \to 0$ as $r\to\infty$, and if $\rho \propto r^{-n}$ at large $r$ and $\Psi \sim r^{-1}$ then $d\rho / d\Psi \sim r^{-n+1}$ which goes to zero provided that $n > 1$. 

We consider both a Plummer and 2-dimensional exponential (Exp2D) stellar profile:

#### Exp2D

The 2D exponential stellar profile is given  by 
$$
\Sigma_\star(R) = A e^{-R / R_s}
$$
where $A$ is some normalization and $R_s$ is the exponential scale radius. The 3D deprojected profile is found through abel inversion (e.g. Rapha...), i.e.
$$
\rho_\star (r) =- \frac{1}{\pi}\int_r^\infty \frac{d\Sigma}{dR} \frac{1}{\sqrt{R^2 - r^2}} dR  = \frac{\Sigma_0}{\pi R_s^2}\,K_0(r/R_s)
$$
where $K_0$ is the 0th order modified Bessel function of the second type. 

#### Plummer

A Plummer profile is defined by
$$
\Sigma(R) = \frac{M}{4\pi R_s^2} \left(1 + \frac{R^2}{R_s^2}\right)^{-2} ,
$$
where $M$ is the total mass and $R_s$ is the characteristic scale radius. The density is 
$$
\rho(r) = \frac{3M}{4\pi\,R_s^3} \left(1 + \frac{r^2}{R_s^2}\right)^{-5/2}.
$$


## Dark Matter simulations

We adopt a softening length of 
$$
h_{\rm grav} = 0.014 \left(\frac{r_s}{2.76\,{\rm kpc}}\right)\left(\frac{N}{10^7}\right)^{-1/2}
$$
see appendix for a discussion of this choice, which is similar to the @power+2003  suggested softening. 



### Analysis

Shrinking spheres centres inspired by @power+2003

- Recursively shrink radius by 0.975 quantile and recalculating centroid until radius is less than ~1kpc
- Remove bound particles (using instantanious potential).
- Use previous snapshot centre and acceleration to predict new centre for next snapshot. 

# Appendix: Convergence tests and Simulation Parameters

Here, we describe some convergence tests to ensure our methods and results are minimally impacted by numerical limitations and assumptions. See @power+2003 for a detailed discussion of various assumptions and parameters used in N-body simulations. 

### Softening

@power+2003 suggest the empirical rule that the ideal softening (balancing integration time and only compromising resolution in collisional regime) is 
$$
h_{grav} = 4 \frac{R_{200}}{\sqrt{N_{200}}}
$$

For our isolation halo ($M_s=2.7$, $r_s=2.76$) and with $10^7$ particles, this works out to be $0.044\,{\rm kpc}$.We adpoted the slightly smaller softening which was reduced by a factor of $\sqrt{10}$ based on initial numerical experiments. However, this choice is likely too small and a larger softening would minimally impact our results but improve computation speed. 

### Particle Number

As discussed in @power+2003, the region where density is converged is related to the region which becomes collisionally relaxed over the time of the universe (so about 5-10 Gyr). (i.e. where the collisionless assumption breaks down). The relaxation timescale is given by $t_{\rm relax} = t_{\rm circ} N(r) / 8\ln N(r)$ or otherwise, 
$$
t_{\rm relax}(r) := t_{\rm circ}(r) \frac{N(r)}{8\,\ln N(r)}
= {t_{\rm circ}(r_{200})} \frac{\sqrt{200}}{8} \frac{N(r)}{\ln N(r)} \left(\frac{\bar \rho (r)}{\rho_{\textrm crit}}\right)^{-1/2}
$$
This works out to be about 6-10 times  (increasing with particle number) our adopted softening length for NFW halos given our assumptions. As such, at full resolution, we can only trust density profiles down to $10\epsilon$ or about $0.14\,{\rm kpc} (r_s/2.76\,{\rm kpc})$, just enough to resolve stellar density profiles. 

Note that by decreasing the truncation radius, the effective particle number is improved, so future experiments could be less generous with this parameter to improve computational performance (order 30% better resolution).

![num convergence](figures/iso_converg_num.png)

Figure: Numerical convergence test for circular velocity as a function of log radius for simulations with different total numbers of particles in isolation. Residuals in lower panel are relative to NFW. The initial conditions are dotted and the converged radius is marked by arrows (REF). 



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





