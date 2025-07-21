

# Numerical convergence and parameters {#sec:extra_convergence}

Here, we describe some convergence tests to ensure our methods and results are minimally impacted by numerical limitations and assumptions. See @power+2003 for a detailed discussion of various assumptions and parameters used in N-body simulations. 

## Softening

One challenge of N-body integration is close, collisional encounters, assuming point particles, cause divergences in the local force. However, this should not occur in a *collisionless* simulation. As a result, most collisionless N-body codes adopt a gravitational softening, a length scale below which the force of gravity begins to decrease between point particles.

@power+2003 empirically suggest that the ideal softening is
$$
h_{\rm grav} = 4 \frac{R_{200}}{\sqrt{N_{200}}},
$$

where $h_{\rm grav}$ is the softening length, and $N_{200}$ is the number of particles within $R_{200}$. This choice balances integration time and only compromises resolution in collisional regime.

For our isolation halo ($M_s=2.7$, $r_s=2.76$) and with $10^7$ particles, this works out to be $0.044\,{\rm kpc}$.We adpoted the slightly smaller softening which was reduced by a factor of $\sqrt{10}$ which appears to improve agreement slightly in the innermost regions. 

![Softening convergence](figures/iso_converg_softening.png){#fig:softening_convergence}

## Time stepping and force accuracy

In general, we use adaptive timestepping and relative opening criteria for gravitational force computations. To verify that these choices and associated accuracy parameters minimally impact convergence or speed, we show a few more isolation runs (using only 1e5 particles)

- constant timestep (...), approximantly half of minimum timestep with adaptive timestepping
- geometric opening, with $\theta = 0.5$.
- strict integration accuracy, (facc = ....)

![Isolation method convergence](figures/iso_converg_methods.png){#fig:methods_convergence}

## Fiducial parameters

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
TimeMax	                    2120     # 10 Gyr

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
SofteningComovingClass0     0.044  # HALO dependent
SofteningMaxPhysClass0      0      # ignored; for cosmological
SofteningClassOfPartType0   0
SofteningClassOfPartType1   0


#=======Miscellanius=======

# probably do not need to change the options below 

#---Unit System
UnitLength_in_cm            1 
UnitMass_in_g               1
UnitVelocity_in_cm_per_s    1 
GravityConstantInternal     1

#---Cosmological Parameters 
ComovingIntegrationOn	      0 # no cosmology
Omega0	                    0
OmegaLambda                 0 
OmegaBaryon                 0
HubbleParam                 1
Hubble                      100
BoxSize                     0

#---SPH
ArtBulkViscConst            0.8
MinEgySpec                  0
InitGasTemp                 100

#---Initial density estimate (SPH)
DesNumNgb                   64
MaxNumNgbDeviation          1
 ```

