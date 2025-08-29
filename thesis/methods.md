In the previous chapters, we have seen that the Scl and UMi classical dSphs have outer density profiles that appear to deviate from the exponential law that approximates well all other classical dSphs. Our main intension is to assess whether such excesses result from the effects of Galactic tides. To this purpose, we intend to use N-body simulations of the evolution of CDM halos in a Galactic potential, constrained to have the orbital parameters consistent with a dwarf's present-day position and velocity. We shall assume that the Galactic potential is the static, analytic potential inferred by @mcmillan+2011 from observations of kinematic tracers. We also assume the potential of each dwarf may be initially approximated by a cuspy NFW profile. Since the dwarfs in question are heavily dark matter dominated, we shall use a carefully selected sample of dark matter particles to mimic and track the evolution of an embedded tracer stellar component. In this Chapter, I describe the Galactic potential used, the orbital estimation method, the initial conditions setup, and the N-body method used. 



# Orbital estimation

To choose the possible orbits of a dwarf galaxy, we perform a Monte Carlo sampling of the present-day observables. The present-day position, distance modulus, LOS velocity, and proper motions are each sampled from normal distributions given the reported uncertainties in [@tbl:scl_obs_props; @tbl:umi_obs_props]. We integrate each sampled observable back in time for 10 Gyr using \agama [@agama]. Dynamical is not expected to impact orbits substantially because of the low masses of the dwarfs, and because the dwarfs orbit mainly in the periphery of the Galaxy. When selecting the initial position and velocity of an N-body model, we select the position and velocity of the first apocentre occurring 10 Gyr ago.

To convert from Gaia to Galactocentric coordinates, we use the Astropy v4 Galactocentric frame [@astropycollaboration+2022]. This frame assumes the Galactic centre is at position $\alpha = {\rm 17h\,45m\,37.224s}$, $\delta = -28^\circ\,56'\,10.23''$ with proper motions $\mu_{\alpha*}=-3.151\pm0.018\ \masyr$ , $\mu_\delta=-5.547\pm0.026 \masyr$ (from the appendix and Table 2 of @reid+brunthaler2004). The Galactic centre distance from the Sun is $8.122\pm0.033\,$kpc with a radial velocity = $11 + 1.9 \pm 3\,\kms$,  from @gravitycollaboration+2018. Finally, adding that the Sun is $20.8\pm0.3\,$pc above the disk from @bennett+bovy2019,  and using  the procedure outlined in @drimmel+poggio2018, the Solar velocity relative to the Galactic rest frame is $\v_\odot = [-12.9 \pm 3.0, 245.6 \pm 1.4, 7.78 \pm 0.08]$ km/s. The uncertainties in the reference frame are typically smaller than the uncertainties on a dwarf galaxy's position and velocity. 

## Milky Way Potential



We adopt the Milky Way potential described in @EP2020, which is an analytic approximation to that proposed by @mcmillan2011. @fig:v_circ_potential plots the circular velocity profiles of each component and the total circular velocity profile for our fiducial profile. The potential includes a stellar bulge, a thin and thick disk, and a dark matter NFW halo. 

The galactic bulge is described by a @hernquist1990 potential,

$$
\Phi(r) = - \frac{GM}{r + a},
$$
where  $a=1.3\,{\rm kpc}$ is the scale radius and $M=2.1 \times 10^{10}\,\Mo$ is the total mass. The thin and thick disks are represented with the @miyamoto+nagai1975  cylindrical potential:

$$
\Phi(R, z) = \frac{-GM}{\left(R^2 + \left[a + \sqrt{z^2 + b^2}\right]^{2}\right)^{1/2}}
$$

where $a$ is the disc radial scale length, $b$ is the scale height, and $M$ is the total mass of the disk. For the thin disk,  $a=3.944\,$kpc, $b=0.311\,$kpc, $M=3.944\times10^{10}\,$M$_\odot$. For the thick disk, $a=4.4\,$kpc, $b=0.92\,$kpc, and $M=2\times10^{10}\,$M$_\odot$. The halo is a NFW dark matter halo ([@eq:nfw]) with $r_{\rm max} = 43.7\,$kpc and $\V_{\rm max} = 191\,\kms$. 



![Circular velocity of potential](figures/v_circ_potential.pdf){#fig:v_circ_potential} 

Figure: Circular velocity profile of @EP2020 potential. 



## Orbits of Sculptor



Next, we calculate the range of possible orbits of Sculptor in our fiducial Milky Way potential. @fig:scl_orbits illustrates point particle orbits for 100 samples of Sculptor's observed kinematics integrated backwards 5 Gyr in both galactocentric coordinate slices (x, y, z) and in galactocentric radius with time. In terms of the x, y, z planes, note that all sampled orbits of Sculptor have nearly the same morphology---they each orbit mostly in the y-z plane, passing through the same number of pericentres and maintaining similar orbits across 5Gyr. Sculptor's orbital history is relatively well-constrained. 

We select an orbit with the $\sim 3\sigma$ smallest pericentre among all possible orbits. We achieve this by taking the median parameters of all orbits with a pericentre less than the 0.0027th quantile pericentre, yielding a pericentre of 43 kpc.^[only slightly smaller than the 0.00135th quantile pericentre] It is unlikely given the current observations that Sculptor has a smaller pericentre. Since tidal effects depend most strongly on the pericentre of an orbit, this orbit should maximize the tidal force from the Milky Way. 

We take the first apocentre after a lookback time of 10 Gyr, or at 9.43 Gyr, as the initial conditions for our model of Sculptor. @tbl:scl_orbits notes the initial conditions for our simulation.





![Sculptor Orbits](/Users/daniel/thesis/figures/scl_xyzr_orbits.pdf){#fig:scl_orbits}

Figure: The orbits of Sculptor in a static Milky Way potential in galactocentric $x$, $y$, and $z$ coordinates. The Milky Way is at the centre with the disk lying in the $x$--$y$ plane. Our selected `smallperi` orbit is plotted in black and light blue transparent orbits represent the past 5Gyr orbits sampled over Sculptor observables in [@tbl:scl_obs_props]. The orbit of sculptor is well-constrained in this potential and it is unlikely to achieve a smaller pericentre than our selected orbit.





| Property             | Mean                    | SmallPeri                 | LMC  |
| -------------------- | ----------------------- | ------------------------- | ---- |
| distance             | 83.2                    | 82.6                      |      |
| pmra                 | 0.099                   | 0.134                     |      |
| pmdec                | -0.160                  | -0.198                    |      |
| Vlos                 | 111.2                   | 111.2                     |      |
| $t_i$                | -8.74                   | -9.43                     |      |
| $\hat{x}_{i}$        | [16.13, 92.47, 39.63]   | [-2.49, -42.78, 86.10]    |      |
| $\vec{v}_i$          | [-2.37, -54.70, 128.96] | [-20.56, -114.83, -57.29] |      |
| pericentre           | 53                      | 43                        |      |
| apocentre            | 102                     | 96                        |      |
| $t_{\rm last\ peri}$ | -0.45                   | -0.46                     |      |
| Numer of peris       |                         |                           |      |

Table: Properties of selected orbits for Sculptor. The mean orbit represents the observational mean from @tbl:scl_obs_props. The Smallperi represents instead the $3\sigma$ smallest pericentre, which we use to provide an upper limit on tidal effects. {#tbl:scl_orbits  short="Sculptor Selected Orbits"}



## Orbits of Ursa Minor



![Ursa Minor Orbits](/Users/daniel/thesis/figures/umi_xyzr_orbits.pdf)

## 

ra = 227.242
dec = 67.2221
distance = 64.6
pmra = -0.158
pmdec = 0.05
radial_velocity = -245.75

t_i = -9.53
pericentre = 29.64
apocentre = 74.88
t last peri = -0.80
x_i = [-16.48 69.92 21.05]
v_i = [16.32 39.86 -116.99]

# Initial conditions

We use \agama [@agama] to generate initial conditions. We initially assume galaxies are described by an NFW dark matter potential [@eq:nfw] and the stars are merely collisionless tracers embedded in this potential (added on in post-processing). The density is a cubic-exponentially truncated with a profile
$$
\rho_{\rm tNFW} = e^{-(r/r_t)^3}\ \rho_{\rm NFW}(t)
$$
where we adopt $r_t = 20 r_s$ or approximately $r_{200}$ for our Sculptor-like fiducial halo. Using $r_{200}$ for $r_t$ would depend on the chosen scale of the halo. So our adopted $r_t$ is an approximate upper limit of $r_{200}$ for typical dwarf galaxy halos [@ludlow+2016]. It is unlikely that the outer density profile of loosely bound particles past $20$ kpc affects the tidal evolution of a subhalo. 

## Halos for Sculptor

From the observed properties of Sculptor, we can infer reasonable $\Lambda$CMD dark matter halo hosts. [@tbl:scl_derived_props] reports our inferred halo and kinematic properties of Sculptor. First, taking the absolute magnitude from @munoz+2018 with the mass-to-light ratio from @woo+courteau+dekel2008 (1.7 with $\sim$ 0.17 dex uncertainty), the total current stellar mass of Sculptor is $M_\star \sim 3.1 \times 10^6 \Mo$. Based on the stellar mass-$v_{\rm max}$ relation [from @fattahi+2018; see also @fig:smhm], the halo should have $v_{\rm max} \sim 31 \pm 3 \,\kms$. Finally, using the @ludlow+2016 $z=0$ mass-concentration relation, this constraint translates into a $r_{\rm max} \sim 6 \pm 2\,{\rm kpc}$, or alternatively $M_{200} = 5\times10^9$ with $c=13$. 

@fig:scl_halos illustrates these estimates visually. The stellar-mass $v_{\rm max}$ constraint translates to an estimate of $v_{\rm max}$ only (horizontal band). The mass-concentration relation describes the relationship between $v_{\rm max}$ and $r_{\rm max}$, an approximately diagonal linear band. And finally, we include a curved line which illustrates where the initial velocity dispersion is approximantly larger than the present-day observed dispersion $v_{\rm circ}(R_h) / \sqrt{3} \sim \sigma_v \gtrsim 9\,\kms$. 

The intersection of the constraints on the halo leads us to select an initial halo, `compact`, which is more compact than the mean but with a observationally consistent initial velocity dispersion. We also have the `light` halo which occupies the intersection between 3$\sigma$ lowest mass and the velocity dispersion constraints, which we come back to in @sec:scl_lmc. 



While there is some range in the choice of initial halo, reasonable changes to the initial halo do not substantially affect the tidal evolution for Sculptor. The velocity dispersion, in particular, translates to an estimate of the mass within $r_h$. As a result, halos with the same velocity dispersion only differ in total mass, but the density slope and circular velocity near the stellar component is similar. These halos will have a similarly affected stellar component. 





| parameter            | value                                         | reference                                                    |
| -------------------- | --------------------------------------------- | ------------------------------------------------------------ |
| $L_\star$            | $1.8\pm0.2\times10^6\ L_\odot$                |                                                              |
| $M_\star$            | $3.1_{-1.0}^{+1.6} \times10^6\ {\rm M}_\odot$ |                                                              |
| $M_\star / L_\star$  | $1.7\times 10^{\pm 0.17}$                     | @woo+courteau+dekel2008                                      |
| $v_{\rm circ, max}$  | $31\pm 3\,\kms$                               | relation from @fattahi+2018                                  |
| $r_{\rm circ, max}$  | $6 \pm 2$ kpc                                 | $v_{\rm circ}$ with @ludlow+2016 mass-concentration relation |
| $M_{200}$            | $0.5 \pm 0.2\times10^{10}\ M_0$               |                                                              |
| $c_{\rm NFW}$        | $13.1_{-2.8}^{+3.6}$                          |                                                              |
| $r_{\rm break, obs}$ | $25 \pm 5$ arcmin                             | @sestito+2024a                                               |
| $t_{\rm break, obs}$ | $110\pm30$ Myr                                | $\sigma_v$, $r_{\rm break}$ with @                           |

Table: To derive total mass, we use the absolute magnitude from @munoz+2018 (see also [@tbl:scl_obs_props]) along with the stellar mass to light ratio (1.7 with 0.17 dex uncertainty) from @woo+courteau+dekel2008. Sculptor's total stellar mass is then $M_\star \sim 3.1_{-1.0}^{1.6} \times 10^6\,\Mo$. From the stellar mass to dark matter halo's characteristic velocity $v_{\rm circ}$ in @fattahi+2018, we expect Sculptor to have $v_{\rm circ} \approx 31 \pm 3 \kms$. {#tbl:scl_derived_props  short="Derived Properties of Sculptor"}



![Sculptor initial halos](/Users/daniel/thesis/figures/scl_initial_halos.pdf){#fig:scl_halos}

Figure: The suggested halos of Sculptor compared to cosmological predictions. **TODO: simplify $\sigma_v$ calculation**.



| Halo name | Rmax | Vmax | M200 | c    |
| --------- | ---- | ---- | ---- | ---- |
| compact   | 3.2  | 31   | 0.33 | 21   |
| lmc       | 4.2  | 31   | 0.39 | 17   |
| small     | 2.5  | 25   |      |      |

Table: The parameters for our initial Sculptor halos. {#tbl:scl_ini_halos  short="Initial halos of Sculptor"}

## Halos for Ursa Minor



| parameter            | value                              | reference               |
| -------------------- | ---------------------------------- | ----------------------- |
| $L_\star$            | $3.5 \pm 0.1 \times 10^5\,L_\odot$ |                         |
| $M_\star$            | $7_{-2}^{+3} \times 10^5\,\Mo$     |                         |
| $M_\star / L_\star$  | 1.9 (pm 0.17 dex)                  | @woo+courteau+dekel2008 |
| $M_{200}$            | $3_{-2}^{+4} \times 10^9\,\Mo$     |                         |
| $c$                  | 14?                                |                         |
| $v_{\rm circ, max}$  | $27_{-6}^{+7}\,\kms$               |                         |
| $r_{\rm circ, max}$  | $5_{-2}^{+1}$ kpc                  |                         |
| $r_{\rm break, obs}$ | $30 \pm 5$ arcmin                  |                         |
| $t_{\rm break, obs}$ | $120\pm30$ Myr                     |                         |

Table: Derived properties of Ursa Minor. {#tbl:umi_derived_props  short="Derived Properties of Ursa Minor"}



| Halo name | Rmax | Vmax | M200 | c    |
| --------- | ---- | ---- | ---- | ---- |
| fiducial  | 5    | 37   | 0.67 | 17   |
| compact   | 4    | 38   | 0.62 | 21   |

Table: Initial halos for Ursa Minor. {#tbl:umi_ini_halos  short="Ursa Minor Initial Halos"}

![Ursa Minor initial halos](/Users/daniel/thesis/figures/umi_initial_halos.pdf)

Figure: The suggested halos of Ursa Minor compared to cosmological predictions.

## 

# N-Body modelling

Modelling gravitational evolution for large systems requires special methods. Perhaps the simplest method to compute the evolution of dark matter is through *N-body simulations*. A dark matter halo is represented as a large number of dark matter particles (bodies). Each body is essentially a Monte Carlo sample of the underlying phase-space distribution. Note that dark matter (and galaxies) are often assumed to be *collisionless*—particles are not strongly affected by close, *collisional* gravitational encounters which substantially change the momenta of involved bodies. In contrast, star clusters are often collisional so neglecting these encounters may not be a reasonable approximation. While we use individual gravitating bodies in N-body simulations, the Newtonian gravitational force is softened to be a Plummer sphere as to limit strongly collisional encounters.

Naively, the Newtonian gravitational force requires adding together the forces from each particle on each particle, with a computational cost that scales quadratically with the number of particles, or $O(N^2)$. With this method, simulating a large number of particles, such as $10^6$, would require $10^{12}$ force evaluations at each time step, making cosmological and high-resolution studies unfeasible. However, only long-range gravitational interactions tend to be important for CDM, so we can utilize the *tree method* to compute the gravitational force vastly more efficiently.

The first gravitational tree code was introduced in @barnes+hut1986, and is still in use today. We utilize the massively parallel code *Gadget 4* [@gadget4].  Particles are spatially split into an *octotree*. The tree construction stars with one large node, a box containing all of the particles. If there is more than one particle in a box/node , the box is then divided into 8 more nodes (halving the side length in each dimension) and this step is repeated until each node only contains 1 particle. With this heirarchichal organization, if a particle is sufficiently far away from a node, then the force is well approximated by the force from the centre of mass of the node. As such, each force calculation only requires a walk through the tree, only descending farther into the tree as necessary to retain accuracy. The total force calculations reduce from $O(N^2)$ to $O(N\,\log N)$, representing orders of magnitude speedup. Modern codes such as *Gadget* utilize other performance tricks, such as splitting particles across many supercomputer nodes, efficient memory storage, adaptive time stepping, and parallel file writing to retain fast performance for large scale simulations, forming the foundation for many cosmological simulation codes. 

**could be 1 paragraph**



## Isolation runs and simulation parameters 

To ensure that the initial conditionss of the simulation are dynamically relaxed and well-converged, we run the simulation in isolation (no external potential) for 5 Gyr (or about 3 times the crossing timescale at the virial radius). Our fiducial isolation halo uses $r_s=2.76$ kpc and $M_s = 0.29 \times 10^{10}\,\Mo$ , but can be easily rescaled for any length or mass scale. 

For our simulation parameters, we adopt a softening length of 
$$
h_{\rm grav} = 0.014 \left(\frac{r_s}{2.76\,{\rm kpc}}\right)\left(\frac{N}{10^7}\right)^{-1/2}.
$$
See Appendix @sec:extra_convergence for a discussion of this choice, which is similar to the @power+2003  suggested softening. We use the relative tree opening criterion with the accuracy parameter set to 0.005, and adaptive time stepping with integration accuracy set to 0.01. 

## Numerical fidelity

@fig:numerical_convergence illustrates how well our numerical setup is able to reproduce the desired initial conditions, before and after running the model in isolation. Circular velocity is computed assuming spherical symmetry and only shown for every 200th particle (ranked from the centre outwards). With increasing particle number, the circular velocity profile maintains close agreement with the expected NFW velocity profile until $r_{\rm relax}$, as marked by arrows. The reduction of mass (and consequently circular velocity) interior to $r_{\rm relax}$ is likely due to collisional effects which would continue to reduce with higher resolution. Typically, $r_{\rm relax}(10\Gyr)$ is about 6-10 times our adopted softening length, increasing with particle number. As such, at full resolution, we can only trust density profiles down to $\sim10$ times the softening length, sufficient to resolve stellar density profiles. 

![Numerical halo convergence](figures/iso_converg_num.pdf){#fig:numerical_convergence}

Figure: Numerical convergence test for circular velocity as a function of log radius for simulations with different total numbers of particles in isolation. Residuals in bottom panel are relative to NFW. The initial conditions are dotted and the converged radius is marked by arrows ([@eq:t_relax]). Note that a slight reduction in density starting around $r = 30 $kpc is expected given our truncation choice. **Add softening**



## Orbital evolution



To perform the simulations of a given galaxy in a given potential, we centre the isolation run's final snapshot ([@sec:shrinking_spheres]) and place the dwarf galaxy in the specified orbit in the given potential. We typically run the simulation for 10 Gyr, which allows us to orbit slightly past the expected initial conditions. 

## Centring {#sec:shrinking_spheres}

Accurately determining the centre of the subhalo at each timestep is essential to most analysis. We use a shrinking-spheres centre method inspired by @power+2003. First, we start with an initial centre estimate from the last timestep. Then, we calculate the radius of all particles from the centre, remove particles with a radius beyond the 0.975 quantile of the centre, and recalculating centroid until radius is less than ~1kpc or fewer than 0.1% of particles remain. Finally, we remove all unbound particles based on the instantaneous N-body potential with all particles. For all future timesteps, we use the snapshot centre's position, velocity, and acceleration to predict the location of the next centre. We also consider only particles retained from the previous iteration.

The statistical centring uncertainty for the full resolution ($10^7$ particle) isolation run is of order 0.003 kpc, but oscillations in the centre are of order 0.03 kpc. This is about three times the softening length but is less than the numerically converged radius scale. 



## Stellar probabilities

We "paint" stars onto dark matter particles using the particle tagging method [e.g. @bullock+johnston2005], assuming spherical symmetry. Let $\Psi$ be the potential (normalized to vanish at infinity) and  ${\cal E}$ is the binding energy ${\cal E} = \Psi - 1/2 v^2$. If we know $f({\cal E})$, the distribution function (phase-space density in energy), then we assign the stellar weight for a given particle with energy ${\cal E}$ is 

$$
P_\star({\cal E}) = \frac{f_\star({\cal E})}{f_{\rm halo}({\cal E})}.
$$
While $f({\cal E})$ is a phase-space density, the differential energy distribution includes an additional $g({\cal E})$ occupation term (BT87). We use Eddington inversion to find the distribution function, (eq. 4-140b in BT87)

$$
f({\cal E}) = \frac{1}{\sqrt{8}\, \pi^2}\left( \int_0^{\cal E} \frac{d^2\rho}{d\Psi^2} \frac{1}{\sqrt{{\cal E} - \Psi}}\ d\Psi + \frac{1}{\sqrt{\cal E}} \left(\frac{d\rho}{d\Psi}\right)_{\Psi=0} \right).
$$

In practice the right, boundary term is zero as $\Psi \to 0$ as $r\to\infty$, and if $\rho \propto r^{-n}$ at large $r$ and $\Psi \sim r^{-1}$ then $d\rho / d\Psi \sim r^{-n+1}$ which goes to zero provided that $n > 1$. We take $\Psi$ from the underlying assumed analytic dark matter potential. $\rho_\star$ can be calculated from the surface density, $\Sigma_\star$, via the inverse Abel transform. 

We find the stellar profiles created in this manner are stable in the isolated systems and retail excellent agreement with the assumed stellar density profile.

