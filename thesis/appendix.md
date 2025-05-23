# Density profile tests





![Density profiles](figures/scl_density_methods_extra.png){#fig:sculptor_observed_profiles}

Figure: Density profiles for various assumptions for Sculptor. PSAT is our fiducial 2-component J+24 sample, circ is a 2-component bayesian model assuming circular radii, simple is the series of simple cuts described, bright is the sample of the brightest half of stars (scaled by 2), DELVE is a sample of  RGB  stars (background subtracted and rescaled to match).

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

# Velocity modelling and comparisons

Here, we describe in additional detail, our methods and comparisons for RV modelling between studies.

Savage-Dickey calculated Bayes factor using Silverman-bandwidth KDE smoothed samples from posterior/prior.



|      | Study       | Instrument | Nspec | Nstar | Ngood | Nmemb | $\delta v_{\rm med}$ | $R_{\rm xmatch}$/arcmin |
| ---- | ----------- | ---------- | ----- | ----- | ----- | ----- | -------------------- | ----------------------- |
| Scl  | combined    |            | 8945  | 2280  | 2034  | 1918  | 0.9                  |                         |
|      | tolstoy+23  | FLAMES     | 3311  | 1701  | 1522  | 1480  | 0.65                 | --                      |
|      | sestito+23a | GMOS       | 2     | 2     | 2     | 2     | 13                   | --                      |
|      | walker+09   | MMFS       | 1818  | 1522  | 1417  | 1329  | 1.8                  | 3                       |
|      | APOGEE      | APOGEE     | 5082  | 253   | 102   | 98    | 0.6                  | --                      |
| UMi  | combined    |            | 4714  | 1225  | 1148  | 831   | 2.3                  |                         |
|      | sestito+23b | GRACES     | 5     | 5     | 5     | 5     | 1.8                  | --                      |
|      | pace+20     | DEIMOS     | 1716  | 1538  | 829   | 682   | 2.5                  | 1                       |
|      | spencer+18  | Hectoshell | 1407  | 970   | 596   | 406   | 0.9                  | ?                       |
|      | APOGEE      | APOGEE     | 9500  | 279   | 37    | 32    | 0.9                  | --                      |

Table: Summary of velocity measurements and derived properties. 



| study      | mean              | sigma           | $\partial \log\sigma / \partial \log R$ | $\partial v_z / \partial x$ (km/s/deg) | $\theta_{\rm grad} / ^{\circ}$ | $\hat R$ | $n_{\rm eff}$ | $\log B_2/B_1$ | WAIC | LOO  |
| ---------- | ----------------- | --------------- | --------------------------------------- | -------------------------------------- | ------------------------------ | -------- | ------------- | -------------- | ---- | ---- |
| all        |                   |                 |                                         |                                        |                                |          |               |                |      |      |
|            | $111.18\pm0.23$   | $9.71\pm0.17$   | -                                       | -                                      | -                              |          |               | 0              |      |      |
|            | $111.19 \pm 0.23$ | $9.68\pm0.17$   | -                                       | $4.3\pm1.3$                            | $-149_{-13}^{+17}$             |          |               | -1.7           |      |      |
|            | $111.16\pm0.23$   | $9.73\pm0.17$   | $0.056\pm0.021$                         | -                                      | -                              |          |               | -0.9           |      |      |
| tolstoy+23 |                   |                 |                                         |                                        |                                |          |               |                |      |      |
|            | $111.3 \pm 0.3$   | $9.80 \pm 0.18$ | -                                       | -                                      | -                              |          |               | 0              |      |      |
|            | $111.3\pm0.3$     | $9.78\pm0.18$   | --                                      | $4.3\pm1.4$                            | $-154_{-13}^{+19}$             |          |               | -1.4           |      |      |
|            | $111.2 \pm 0.3$   | $9.74\pm0.19$   | $0.085 \pm 0.023$                       | --                                     | --                             |          |               | -4.5           |      |      |
| walker+09  |                   |                 |                                         |                                        |                                |          |               |                |      |      |
|            | $111.0\pm0.3$     | $9.57\pm0.21$   | --                                      | --                                     | --                             |          |               | 0              |      |      |
|            | $111.1\pm0.3$     | $9.53\pm0.21$   | -                                       | $5.2_{-1.6}^{+1.8}$                    | $-134_{-16}^{+23}$             |          |               | -2.4           |      |      |
|            | $111.0\pm0.3$     | $9.61\pm0.21$   | $0.03\pm0.03$                           | --                                     | --                             |          |               | +1.6           |      |      |
| apogee     |                   |                 |                                         |                                        |                                |          |               |                |      |      |
|            | $109.9\pm0.8$     | $8.3\pm0.6$     | --                                      | --                                     | --                             |          |               | --             |      |      |
|            | $109.9\pm0.8$     | $8.3\pm0.6$     | --                                      | $6\pm3$                                | $-151_{-36}^{+44}$             |          |               | +0.4           |      |      |
|            | $109.9\pm0.8$     | $8.3\pm0.7$     | $0.05\pm0.08$                           | --                                     | --                             |          |               | +1.1           |      |      |

Table: MCMC fits for  different RV datasets for Sculptor among 3 different models. {#tbl:scl_rv_mcmc short="Sculptor RV fits"}





| study   | mean           | sigma               | $\hat R$ | $n_{\rm eff}$ | $\log bf_{\rm sigma}$ | $\log bf_{\rm grad}$ |
| ------- | -------------- | ------------------- | -------- | ------------- | --------------------- | -------------------- |
| all     | $-245.9\pm0.3$ | $8.76\pm0.24$       |          |               | +1.5                  | +2.2                 |
| pace    | $-244.5\pm0.4$ | $9.1\pm0.3$         |          |               | +0.2                  | +1.1                 |
| spencer | $-246.9\pm0.4$ | $8.8\pm0.3$         |          |               | +1.8                  | -0.3                 |
| apogee  | $-248.2\pm1.6$ | $9.0_{-1.1}^{+1.3}$ |          |               | +0.8                  | +0.8                 |

Table: MCMC fits for UMi velocity dispersion. {#tbl:umi_rv_mcmc short="Ursa Minor RV fits"}



![Scl velocity sample](/Users/daniel/thesis/figures/scl_rv_scatter.pdf)

Figure: RV members of Sculptor plotted in the tangent plane coloured by corrected velocity difference from mean $v_z - \bar v_z$ . The black ellipse marks the half-light radius in @fig:scl_selection. The black and green arrows mark the proper motion (PM, GSR frame) and derived velocity gradient (rot) vectors (to scale).



![Scl velocity gradient](/Users/daniel/thesis/figures/scl_vel_gradient_scatter.pdf) 

Figure: The corrected LOS velocity along the best fit rotational axis. RV members are black points, the systematic $v_z$ is the horizontal grey line, blue lines represent the (projected) gradient from MCMC samples, and the orange line is a rolling median (with a window size of 50).

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

