# Simulation Methods

In this section, we discuss our simulation methods from orbital estimation to initial conditions to simulation parameters and software.

We use Agama [@vasiliev] for potential specification and initial conditions, and Gadget-4 [@springle2021] for orbit integration and dark matter simulations (with a custom patch to integrate with Agama). 

## Potential 

Our fiducial potential is as described in @EP2020. This potential is an analytic approximation of @McMillan2011, consisting of a stellar Bulge, disk, and a dark matter NFW halo. 

The @hernquist1990 density profile for the galactic bulge is parameterized in terms of a characteristic mass, $M$, and radius $a$. 

Density Profile
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

**Figure:** Circular velocity profile of @EP2020 potential.

## Orbital Estimation

To estimate the orbits of a dwarf galaxy, we perform a Monte Carlo sampling of the present-day observables. We integrate each sampled observable back in time for 10 Gyr in Gadget as massless point particles outputing every 0.5Myr (with otherwise similar parameters to n-body runs below.)

To convert from Gaia to galactocentric coordinates, we use the astropy v4 galactocentric frame. This frame assumes an 

- J2000 'ra' = 17h45m37.224s, 'dec' = -28°56'10.23" (appendix), 'pmracosdec=-3.151±0.018 mas/yr', 'pmdec=-5.547±0.026' (Table 2). from @reid+brunthaler2004
- distance = 8.122 ± 0.033 kpc, solar radial velocity = 11 + 1.9 ± 3 km/s.  from  @gravitycollaboration+2018.
- 'z_sun' = 20.8±0.3pc from bennett+bovy2019
- 'v_sun' is slightly more complicated, relying on parameters from above and using  the procedure described in @drimmel+poggio2018. This results in 
  - 'v_sun' = [-12.9 ± 3.0, 245.6 ± 1.4, 7.78 ± 0.08] km/s

## Initial conditions

We use Agama [@vasiliev2019] to generate initial conditions. We initially assume galaxies are described by an NFW dark matter potential and the stars are merely collisionless tracers embedded in this potential

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



### Stellar profiles

We consider several different stellar profiles in this work. 

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



See @power+2003 for a discussion of this.

The gravitational softening is based on:
$$
h_{grav} = 4 \frac{R_{200}}{\sqrt{N_{200}}}
$$

I in fact find that dividing the softening by a factor of sqrt(10) gives slightly more precise results with a negligible increase in compute time.

# Appendix: Convergence tests and Simulation Parameters

Here, we describe some convergence tests to ensure our methods and results are minimally impacted by numerical limitations and assumptions. See @power2003 for a detailed discussion of various assumptions and parameters used in N-body simulations. 

### Particle Number. 

### Timestepping

### Gravitational methods

### Gravitational accuracy

### Fiducial Parameters

 ```
 # -----------------------------------IO---------------------------------------
 
 #---Filenames
 InitCondFile                initial
 OutputDir                   ./out
 SnapshotFileBase            snapshot
 OutputListFilename          outputs.txt
 
 #---File formats 
 ICFormat                    3       # use HDF5
 SnapFormat                  3 
 
 #---Mem & CPU limits
 TimeLimitCPU                86400   # TUNE
 CpuTimeBetRestartFile       7200    # TUNE 
 MaxMemSize                  2400    # TUNE to system (MB)
 
 #---Time
 TimeBegin                   0
 TimeMax	                   2120      # 10 Gyr
 
 #---Output frequency
 OutputListOn                0
 TimeBetSnapshot             10      # change as needed
 TimeOfFirstSnapshot         0       # MATCH to TimeBegin
 TimeBetStatistics           10      # change as needed 
 NumFilesPerSnapshot         1
 MaxFilesWithConcurrentIO    1 
 
 
 # ----------------------------------Gravity-----------------------------------
 
 # MATCH all the below if continuing run
 
 #---Timestep accuracy
 ErrTolIntAccuracy           0.01    # MATCH
 CourantFac                  0.1     # ignored; for SPH
 MaxSizeTimestep             0.5     # MATCH
 MinSizeTimestep             0.0 
 
 #---Tree algorithm
 TypeOfOpeningCriterion      1       # 0: Barnes-Hut, 1: Relative
 ErrTolTheta                 0.5     # mostly used for Barnes-Hut
 ErrTolThetaMax              1.0     # (used only for relative)
 ErrTolForceAcc              0.005   # (used only for relative)
 
 #---Domain decomposition: should only affect performance
 TopNodeFactor                       3.0
 ActivePartFracForNewDomainDecomp    0.02
 
 #---Gravitational Softening
 SofteningComovingClass0      # FILL IN
 SofteningMaxPhysClass0       0      # ignored; for cosmological
 SofteningClassOfPartType0    0
 SofteningClassOfPartType1    0
 
 
 #--------------------------------Misc-------------------------------------
 
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





