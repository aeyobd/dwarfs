# Appendix / Extra Notes

## Additional density profile tests



![Density profiles](figures/scl_density_methods_extra.png){#fig:sculptor_observed_profiles}

Figure: Density profiles for various assumptions for Sculptor. PSAT is our fiducial 2-component J+24 sample, circ is a 2-component bayesian model assuming circular radii, simple is the series of simple cuts described in Appendix ?, bright is the sample of the brightest half of stars (scaled by 2), DELVE is a sample of  RGB  stars (background subtracted and rescaled to match).

Note that a full rigorous statistical analysis would require a simulation study of injecting dwarfs into Gaia and assessing the reliability of various methods of membership and density profiles. This is beyond the scope of this thesis. 



```
SELECT TOP 1000
       *
FROM delve_dr2.objects
WHERE 11 < ra
and ra < 19
and -37.7 < dec
and dec < -29.7

```

## Velocity fit parameters.

Savage-Dickey calculated bayes factor

*Is it interesting that the velocity dispersion of Scl seems to increase significantly with Rell?* 

| study      | mean              | sigma           | $\partial \log\sigma / \partial \log R$ | $\partial v_z / \partial x$ (km/s/deg) | $\theta_{\rm grad} / ^{\circ}$ | $\log bf$ |
| ---------- | ----------------- | --------------- | --------------------------------------- | -------------------------------------- | ------------------------------ | --------- |
| all        |                   |                 |                                         |                                        |                                |           |
|            | $111.17\pm0.22$   | $9.68\pm0.17$   | -                                       | -                                      | -                              | 0         |
|            | $111.18 \pm 0.23$ | $9.64\pm0.16$   | -                                       | $4.8\pm1.3$                            | $-147_{-12}^{+15}$             | -3.7*     |
|            | $111.15\pm0.22$   | $9.70\pm0.16$   | $0.055\pm0.021$                         | -                                      | -                              | -0.8      |
| tolstoy+23 |                   |                 |                                         |                                        |                                |           |
|            | $111.2 \pm 0.3$   | $9.77 \pm 0.18$ | -                                       | -                                      | -                              | 0         |
|            | $111.2\pm0.3$     | $9.73\pm0.18$   | --                                      | $4.8\pm1.4$                            | $-154_{-12}^{+16}$             | -2.5      |
|            | $111.2 \pm 0.3$   | $9.71\pm0.18$   | $0.081 \pm 0.023$                       | --                                     | --                             | -3.3      |
| walker+09  |                   |                 |                                         |                                        |                                |           |
|            | $111.0\pm0.3$     | $9.57\pm0.21$   | --                                      | --                                     | --                             | 0         |
|            | $111.1\pm0.3$     | $9.54\pm0.21$   | -                                       | $5.3_{-1.6}^{+1.8}$                    | $-134_{-16}^{+22}$             | -2.0      |
|            | $111.0\pm0.3$     | $9.61\pm0.21$   | $0.03\pm0.03$                           | --                                     | --                             | +1.6      |
| apogee     |                   |                 |                                         |                                        |                                |           |
|            | $109.9\pm0.8$     | $8.3\pm0.6$     | --                                      | --                                     | --                             | --        |
|            | $109.9\pm0.8$     | $8.3\pm0.6$     | --                                      | $6\pm3$                                | $-151_{-36}^{+44}$             | +0.3      |
|            | $109.9\pm0.8$     | $8.3\pm0.7$     | $0.05\pm0.08$                           | --                                     | --                             | +1.1      |

Table: MCMC fits for  different RV datasets for Scl amoung 3 different models.





| study   | mean           | sigma               | $\log bf_{\rm sigma}$ | $\log bf_{\rm grad}$ |
| ------- | -------------- | ------------------- | --------------------- | -------------------- |
| all     | $-245.9\pm0.3$ | $8.76\pm0.24$       | +1.5                  | +2.2                 |
| pace    | $-244.5\pm0.4$ | $9.1\pm0.3$         | +0.2                  | +1.1                 |
| spencer | $-246.9\pm0.4$ | $8.8\pm0.3$         | +1.8                  | -0.3                 |
| apogee  | $-248.2\pm1.6$ | $9.0_{-1.1}^{+1.3}$ | +0.8                  | +0.8                 |

Table: MCMC fits for UMi velocity dispersion.



|      | Study       | Instrument | Nspec | Nstar | Ngood | Nmemb | $\delta v_{\rm med}$ |
| ---- | ----------- | ---------- | ----- | ----- | ----- | ----- | -------------------- |
| Scl  | combined    |            | 8945  | 2280  | 2034  | 1920  | 0.9                  |
|      | tolstoy+23  | FLAMES     | 3311  | 1701  | 1522  | 1481  | 0.65                 |
|      | sestito+23a | GMOS       | 2     | 2     | 2     | 2     | 13                   |
|      | walker+09   | MMFS       | 1818  | 1522  | 1417  | 1330  | 1.8                  |
|      | APOGEE      | APOGEE     | 5082  | 253   | 102   | 98    | 0.6                  |
| UMi  | combined    |            | 4714  | 1225  | 1148  | 831   | 2.3                  |
|      | sestito+23b | GRACES     | 5     | 5     | 5     | 5     | 1.8                  |
|      | pace+20     | DEIMOS     | 1716  | 1538  | 829   | 682   | 2.5                  |
|      | spencer+18  | Hectoshell | 1407  | 970   | 596   | 406   | 0.9                  |
|      | APOGEE      | APOGEE     | 9500  | 279   | 37    | 32    | 0.9                  |

Table: Summary of velocity measurements and derived properties. sestito+2023a number of members depends on spatial model used. 

# Numerical Convergence

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

