In the previous chapters, we have seen that the Scl and UMi classical dSphs have outer density profiles that appear to deviate from the exponential law that approximates well all other classical dSphs. Our main intention is to assess whether such excesses result from the effects of Galactic tides. To this purpose, we intend to use N-body simulations of the evolution of CDM halos in a Galactic potential, constrained to have the orbital parameters consistent with a dwarf's present-day position and velocity. We shall assume that the Galactic potential is the static, analytic potential inferred by @mcmillan2011 from observations of kinematic tracers. We also assume the potential of each dwarf may be initially approximated by a cuspy NFW profile. Since the dwarfs in question are heavily dark matter dominated, we shall use a carefully selected sample of dark matter particles to mimic and track the evolution of an embedded tracer stellar component. In this Chapter, I describe the Galactic potential used, the orbital estimation method, the initial conditions setup, and the N-body method used. 



# Orbital estimation

To explore the possible orbits of a dwarf galaxy, we perform a Monte Carlo sampling of the present-day observables. The present-day position, distance modulus, LOS velocity, and proper motions are each sampled from normal distributions given the reported uncertainties in [@tbl:scl_obs_props; @tbl:umi_obs_props]. We integrate each sampled position/velocity back in time for 10 Gyr using \agama{} [@agama]. Dynamical friction is not expected to impact orbits substantially because of the low masses and large pericentres of the dwarfs. 

## Galactocentric frame

To convert observed positions and velocities to Galactocentric coordinates, we use the Astropy v4 Galactocentric frame [@astropycollaboration+2022]. This frame assumes the Galactic centre is at position $\alpha = {\rm 17h\,45m\,37.224s}$, $\delta = -28^\circ\,56'\,10.23''$ with proper motions $\mu_{\alpha*}=-3.151\pm0.018\ \masyr$ , $\mu_\delta=-5.547\pm0.026 \masyr$ [from the appendix and Table 2 of @reid+brunthaler2004]. The Galactic centre is at a distance from the Sun of $8.122\pm0.033\,$kpc with a radial velocity = $11 + 1.9 \pm 3\,\kms$  [@gravitycollaboration+2018]. The Sun is assumed to be $20.8\pm0.3\,$pc above the disk [@bennett+bovy2019]. Using  the procedure outlined in @drimmel+poggio2018, the Solar velocity relative to the Galactic rest frame is then $\V_\odot = [-12.9 \pm 3.0, 245.6 \pm 1.4, 7.78 \pm 0.08]$ km/s. The uncertainties in the reference frame are typically smaller than the uncertainties on a dwarf galaxy's position and velocity. 

## Milky Way Potential

We adopt the Milky Way potential described in @EP2020, which is an analytic approximation to that proposed by @mcmillan2011. @fig:v_circ_potential plots the circular velocity profiles of each component and the total circular velocity profile for our fiducial profile. The potential includes a stellar bulge, a thin and thick disk, and a dark matter NFW halo. 

The Galactic bulge is described by a @hernquist1990 potential,

$$
\Phi(r) = - \frac{GM}{r + a},
$$
where  $a=1.3\,{\rm kpc}$ is the scale radius and $M=2.1 \times 10^{10}\,\Mo$ is the total mass. The thin and thick disks are represented with the @miyamoto+nagai1975  cylindrical potential,
$$
\Phi(R, z) = \frac{-GM}{\left(R^2 + \left[a + \sqrt{z^2 + b^2}\right]^{2}\right)^{1/2}},
$$
where $a$ is the disc radial scale length, $b$ is the scale height, and $M$ is the total mass of the disk. For the thin disk,  $a=3.944\,$kpc, $b=0.311\,$kpc, and $M=3.944\times10^{10}\,$M$_\odot$. For the thick disk, $a=4.4\,$kpc, $b=0.92\,$kpc, and $M=2\times10^{10}\,$M$_\odot$. The halo is an NFW dark matter halo ([@eq:nfw]) with $\rmax = 43.7\,$kpc and $\vmax = 191\,\kms$. 



![Circular velocity of potential](figures/v_circ_potential.pdf){#fig:v_circ_potential} 

Figure: Circular velocity profile of @EP2020 potential. The total circular velocity (thick black line) is composed of an NFW halo (green dashed line), a think and thick @miyamoto+nagai1975 disk (orange dash-dotted line), and a @hernquist1990 bulge (blue dotted line).

## Orbits of Sculptor

Sculptor's orbital history is relatively well-constrained. @fig:scl_orbits illustrates point particle orbits for 100 samples of Sculptor's observed kinematics integrated backwards 5 Gyr in both Galactocentric coordinate slices ($x$, $y$, $z$) and in Galactocentric radius with time. All sampled orbits of Sculptor have nearly the same morphology---the orbit primarily resides in the $y$–$z$ plane and has a similar number of periods and pericentre/apocentre. 

To maximize tidal effects, we select an orbit with the $\sim 3\sigma$ smallest pericentre among all possible orbits. We achieve this by taking the median parameters of all orbits with a pericentre less than the 0.0027th quantile pericentre, yielding a pericentre of 43 kpc. Given the current observations, it is unlikely that Sculptor has a significantly smaller pericentre than our selected orbit. 

We take the first apocentre after a look-back time of 10 Gyr, or at 9.43 Gyr, as the initial conditions for our model of Sculptor, noted in @tbl:scl_orbits.





![Sculptor Orbits](/Users/daniel/thesis/figures/scl_xyzr_orbits.pdf){#fig:scl_orbits}

Figure: The orbits of Sculptor in a static Milky Way potential in Galactocentric $x$, $y$, and $z$ coordinates (top) and in Galactocentric radius $r$ versus time (bottom). The Milky Way is at the centre with the disk lying in the $x$--$y$ plane. Our selected `smallperi` orbit is plotted in black and light blue transparent orbits represent the past 5Gyr orbits sampled over Sculptor observables in [@tbl:scl_obs_props]. The orbit of sculptor is well-constrained in this potential and it is unlikely to achieve a smaller pericentre than our selected orbit.



![Ursa Minor Orbits](/Users/daniel/thesis/figures/umi_xyzr_orbits.pdf){#fig:umi_orbits}

Figure: Similar to @fig:scl_orbits, the orbits of Ursa Minor in a static Milky Way potential in Galactocentric $x$, $y$, and $z$ coordinates. In the lower panel, we show the radius versus time for only three orbits of Ursa Minor.



| Property                         | Mean                    | SmallPeri                 |
| -------------------------------- | ----------------------- | ------------------------- |
| distance / kpc                   | 83.2                    | 82.6                      |
| $\pmra / \masyr$                 | 0.099                   | 0.134                     |
| $\pmdec / \masyr$                | -0.160                  | -0.198                    |
| LOS velocity / $\kms$            | 111.2                   | 111.2                     |
| $t_i / \Gyr$                     | -8.74                   | -9.43                     |
| $\vec{x}_{i} / \kpc$             | [16.13, 92.47, 39.63]   | [-2.49, -42.78, 86.10]    |
| $\vec{v}_i / \kms$               | [-2.37, -54.70, 128.96] | [-20.56, -114.83, -57.29] |
| pericentre / kpc                 | 53                      | 43                        |
| apocentre / kpc                  | 102                     | 96                        |
| $t_{\rm last\ peri} / {\rm Gyr}$ | -0.45                   | -0.46                     |
| Number of pericentres            | 6                       | 6                         |

Table: Properties of selected orbits for Sculptor. The mean orbit represents the observational mean from @tbl:scl_obs_props. The Smallperi represents instead the $3\sigma$ smallest pericentre, which we use to provide an upper limit on tidal effects. {#tbl:scl_orbits  short="Sculptor Selected Orbits"}



## Orbits of Ursa Minor

Similar to Sculptor, Ursa Minor has a well-constrained orbit. @fig:umi_orbits shows 100 random point-orbits of Ursa Minor. Initially, we select an orbit with approximately the $3\sigma$ smallest pericentre as in Sculptor. 

Ursa Minor's N-body orbit diverges from the point particle orbit. To ensure the final conditions of the N-body simulation are close to the intended final position, we iteratively adjust Ursa Minor's initial conditions. Initially, starting with low-resolution runs, we adjust the cylindrical actions of the initial orbit by the final difference in actions at the end of orbital evolution. After the initial actions have converged (2 iterations), we change the initial action angles by the final difference in action angles. This method converges within 4 iterations to an orbit agreeing with the observed kinematics of Ursa Minor.



| Property                    | Mean                    | SmallPeri.1             | SmallPeri.5             |
| --------------------------- | ----------------------- | ----------------------- | ----------------------- |
| distance / kpc              | $70.1$                  | 64.6                    |                         |
| $\pmra / \masyr$            | $-0.124 \pm 0.17$       | -0.158                  |                         |
| $\pmdec / \masyr$           | $0.078\pm0.17$          | 0.05                    |                         |
| $v_{\rm LOS} / \kms$        | $-245.9\pm 1$           | -245.75                 |                         |
| $t_i / \Gyr$                | -8.74                   | -9.53                   | -9.53                   |
| $\vec{x}_{i} / \kpc$        | [4.88, -65.11, 50.78]   | [-16.48 69.92 21.05]    | [17.40, 74.51, 21.34]   |
| $\vec{v}_i / \kms$          | [-34.28, 77.11, 101.38] | [16.32, 39.86, -116.99] | [14.27, 48.62, -114.08] |
| pericentre / kpc            | 37                      | 29                      | 28                      |
| apocentre / kpc             | 83                      | 75                      | 72                      |
| $t_{\rm last\ peri} / \Gyr$ | $-0.96$                 | $-0.80$                 | $-0.81$                 |
| Number of pericentres       | 9                       | 6                       | 6                       |

Table: Properties of selected orbits for Ursa Minor. The "smallperi" orbit is the initial point orbit and the "smallperi.5" is the initial orbit for the N-body simulation. {#tbl:umi_orbits short="Ursa Minor Selected Orbits"}

# Initial conditions

We use \agama{} [@agama] to generate the initial N-body dark matter halo. We assume galaxies are described by an NFW dark matter potential ([@eq:nfw]). We also assume the stars do not contribute to the potential. The dark matter density is truncated in the outer regions by
$$
\rho_{\rm tNFW} = e^{-(r/r_t)^3}\ \rho_{\rm NFW}(t),
$$
where we adopt $r_t = 20 r_s$. 

## Initial dark matter halos for Sculptor and Ursa Minor

From the observed properties of Sculptor and Ursa Minor, we infer reasonable \LCDM{} initial halo conditions. 

[@tbl:derived_props] reports our inferred halo and kinematic properties of Sculptor and Ursa Minor. First, taking the absolute magnitudes from @munoz+2018 with the mass-to-light ratio from @woo+courteau+dekel2008 (with $\sim$ 0.17 dex uncertainty), the total current stellar mass of Sculptor and Ursa Minor are $M_\star \sim 3.1 \times 10^6 \Mo$ and $M_\star \sim 7 \times 10^5 \Mo$. Based on the stellar mass-$\vmax$ relation [from @fattahi+2018; see also @fig:smhm], Sculptor and Ursa Minor's halos should have $\vmax \approx 31 \,\kms$ and $\vmax \approx 27\,\kms$. Finally, using the @ludlow+2016 $z=0$ mass-concentration relation, this constraint translates into a $\rmax \approx 6 {\rm kpc}$ and $\rmax \approx 5\,\kpc$ for each galaxy.

@fig:scl_halos illustrates these estimates visually for both galaxies. The stellar-mass $\vmax$ constraint translates to an estimate of $\vmax$ only (horizontal band). The mass-concentration relation describes the relationship between $\vmax$ and $\rmax$ (diagonal linear band). And finally, we include a curved line which illustrates the halo $\vmax$ which has a specified initial LOS velocity dispersion given $\rmax$ $, \vcirc(R_h) / \sqrt{3} \approx \sigma_v$. $R_h \approx 0.24\,\kpc$ for both galaxies. 

While there is some range in the choice of initial halo, reasonable changes to the initial halo do not substantially affect the tidal evolution for either galaxy. The observed velocity dispersion is directly related to the mass within a half-light radius [e.g., @wolf+2010]. As a result, halos with the same velocity dispersion may differ in total mass but should have similar mass within $R_h$. So, the tidal effects on stars should be similar for halos with similar velocity dispersions.

Our selected halos (in [@tbl:initial_halos]) for each simulation run are based on the cosmological constraints and the need to match the present-day velocity dispersion at the end of the simulation. 



| parameter           | Sculptor                                      | Ursa Minor                         |
| ------------------- | --------------------------------------------- | ---------------------------------- |
| $L_\star$           | $1.8\pm0.2\times10^6\ L_\odot$                | $3.5 \pm 0.1 \times 10^5\,L_\odot$ |
| $M_\star$           | $3.1_{-1.0}^{+1.6} \times10^6\ {\rm M}_\odot$ | $7_{-2}^{+3} \times 10^5\,\Mo$     |
| $M_\star / L_\star$ | $1.7\times 10^{\pm 0.17}$                     | $1.9 \times 10^{\pm 0.17}$         |
| $\vmax$             | $31\pm 3\,\kms$                               | $27_{-6}^{+7}\,\kms$               |
| $\rmax$             | $6 \pm 2$ kpc                                 | $5_{-2}^{+1}$ kpc                  |
| $M_{200}$           | $0.5 \pm 0.2\times10^{10}\ M_0$               | $3_{-2}^{+4} \times 10^9\,\Mo$     |
| $c_{\rm NFW}$       | $13_{-3}^{+4}$                                | 14?                                |

Table: Inferred properties of the stellar component and halo for Sculptor and Ursa Minor. We record the total luminosity, stellar mass, mass-to-light ratio, dark matter halo $\vmax$ and $\rmax$, and dark matter halo virial mass $M_{200}$ and concentration $c_{\rm NFW}$. {#tbl:derived_props  short="Derived Properties of Sculptor and Ursa Minor"}



![Sculptor initial halos](figures/initial_halo_selection.pdf){#fig:scl_halos}

Figure: Selection of initial halos for Sculptor and Ursa Minor. The grey line and pink line with shaded regions represent the @ludlow+2016 mass-concentration relation and @fattahi+2018 SMHM relation respectively. The curved lines represent the velocity dispersion of the initial halo given the present-day half-light radius [via the @wolf+2010 mass estimator].



| Halo name     | $\rmax$ | $\vmax$ | $M_{200}$ | $c_{\rm NFW}$ |
| ------------- | ------- | ------- | --------- | ------------- |
| Scl: fiducial | 3.2     | 31      | 0.33      | 21            |
| Scl: small    | 2.5     | 25      |           |               |
| UMi: fiducial | 4       | 38      | 0.62      | 21            |

Table: The initial conditions for our initial dark matter halos. {#tbl:initial_halos  short="Initial halos"}

# Numerical methods

To simulate the tidal evolution of galaxies, we use "N-body" simulations integrated with the parallel, gravitational-tree program \gadget{} [@gadget4]. The N-body method calculates and evolves the gravitational accelerations between a large number of collisionless particles to approximate the dynamical evolution of matter. To approximate a collisionless system (i.e., without strong particle-particle gravitational deflections), the force is softened below a "softening length." Regions which are smaller than the softening length, or have few particles, are not well-resolved. Tree codes organize particles into a spatial tree, enabling grouping of the gravitational forces from nearby particles. A gravitational tree code substantially reduces the required number of force computations, as compared to the exact, direct-summation method. 

## Isolation runs and simulation parameters 

To ensure that the initial conditions of the simulation are dynamically relaxed and well-converged, we run a halo first in isolation using \gadget{}. Since gravity is scale-free, we use the same isolation run for all halos. We adopt $\rmax = 6.0\,$kpc and $\vmax = 31\,\kms$  for the isolation halo based on Sculptor's mean properties. We run this model for 5 Gyr (about three times the free fall timescale $t_{\rm ff} = \vcirc / r \approx 1.5\,\Gyr$ at $r_{200}=36\,$kpc).

For our simulation parameters, we adopt a softening length of 
$$
h_{\rm grav} = 0.014\,{\rm kpc}\left(\frac{r_{\rm max}}{6.0\,{\rm kpc}}\right)\left(\frac{N}{10^7}\right)^{-1/2}.
$$ {#eq:softening_length}
See Appendix @sec:extra_convergence for a discussion of this choice, which is similar to the @power+2003 suggested softening. We use the relative tree opening criterion with the accuracy parameter set to 0.005, and adaptive time stepping with integration accuracy set to 0.01. 

## Numerical fidelity

@fig:numerical_convergence illustrates how well our numerical setup is able to reproduce the desired initial conditions, before and after running the model in isolation. This figure shows that our numerical methods are able to approximate well an NFW halo down to an innermost radius that strongly depends on resolution. The larger the number of particles, the smaller the radius that is effectively "resolved" in a given simulation. For the Sculptor halo shown in this figure (with $\rmax = 6.0\,$kpc and $\vmax = 31\,\kms$ ), a simulation with $10^7$ particles is needed to resolve the innermost 100 pc. For reference, the half-light radius of Sculptor is roughly 100 pc, which means that at least 10 million particles would be needed to follow faithfully its tidal evolution. Vertical arrows in @fig:numerical_convergence indicate the "convergence radius" defined by @power+2003 [eq.~13] for NFW halos formed in cosmological N-body simulations. This radius marks the region where collisional effects driven by the finite number of particles used to describe the innermost regions of a halo become important. The softening length (from @eq:softening_length) is typically a few times smaller than the converged length.

![Numerical halo convergence](figures/iso_converg_num.pdf){#fig:numerical_convergence}

Figure: Numerical convergence test for circular velocity as a function of log radius for simulations with different total numbers of particles in isolation. Residuals in bottom panel are relative to NFW. The initial conditions are dotted, the converged radius is marked by arrows [from @power+2003 eq. 13], and the softening length is marked by a vertical bar. Note that a slight reduction in density starting around $r = 30\,\kpc$ is expected given our truncation choice. 



## Orbital evolution

Next, we evolve the halo in the Galactic potential. We scale the relaxed snapshot and softening length to match the initial halo [in @tbl:initial_halos], and shift the snapshot to the initial conditions from the orbital analysis (see [@tbl:scl_orbits; @tbl:umi_orbits]). We then evolve the full N-body NFW model forward in time in the Galactic potential and follow it in time until the present time, when the halo is closest to the observed position of the galaxy. 

## Tidal mass losses {#sec:shrinking_spheres}

To accurately follow the evolution of a halo, it is important to determine the centre of the halo at each time chosen for analysis. We use a shrinking-spheres centre method inspired by @power+2003. First, we start with an initial centre estimate from the last timestep. Then, we calculate the radius of all particles from the centre, remove particles with a radius beyond the 0.975 quantile of the centre, and recalculating the centre of mass. The procedure is repeated until the selection radius is less than ~1kpc or fewer than 0.1% of particles remain. After a centre has been chosen, we remove all unbound particles based on the \gadget{} calculated potential of the halo. For all future timesteps, we consider only particles retained from the previous iteration.

The statistical centring uncertainty for the full resolution ($10^7$ particle) isolation run is of order 0.003 kpc, but fluctuations are observed of  $\sim 0.03\,\kpc$. This is about three times the softening length but is less than the numerically converged radius scale. 



## The stellar component {#sec:painting_stars}

We "paint" stars onto dark matter particles using the particle-tagging method [e.g. @bullock+johnston2005], assuming spherical symmetry. We initially assume exponential stars with $R_s = 0.10\,,\kpc$ for both galaxies.  

Let $\Psi$ be the potential (normalized to vanish at infinity) and  ${\cal E}$ the binding energy ${\cal E} = \Psi - 1/2 v^2$. If we know $f({\cal E})$, the distribution function (phase-space density in energy), then we assign a stellar weight to a given particle with energy ${\cal E}$ using
$$
P_\star({\cal E}) = \frac{f_\star({\cal E})}{f_{\rm halo}({\cal E})}.
$$
While $f({\cal E})$ is a phase-space density, the differential energy distribution includes an additional $g({\cal E})$ occupation term (BT87). We use Eddington inversion to find the distribution function, (eq. 4-140b in BT87)
$$
f({\cal E}) = \frac{1}{\sqrt{8}\, \pi^2}\left( \int_0^{\cal E} \frac{d^2\rho}{d\Psi^2} \frac{1}{\sqrt{{\cal E} - \Psi}}\ d\Psi + \frac{1}{\sqrt{\cal E}} \left(\frac{d\rho}{d\Psi}\right)_{\Psi=0} \right).
$$
In practice the right, boundary term is zero as $\Psi \to 0$ as $r\to\infty$, and if $\rho \propto r^{-n}$ at large $r$ and $\Psi \sim r^{-1}$ then $d\rho / d\Psi \sim r^{-n+1}$ which goes to zero provided that $n > 1$. We take $\Psi$ from the underlying assumed analytic dark matter potential. $\rho_\star$ can be calculated from the surface density, $\Sigma_\star$, via the inverse Abel transform. 

We find the stellar profiles created in this manner are stable in the isolated systems and show excellent agreement with the assumed stellar density profile.

