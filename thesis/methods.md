In the previous Chapters, we have seen that the Scl and UMi classical dSphs have outer density profiles that appear to deviate from the exponential law that approximates well all other classical dSphs. Our main intention is to assess whether such deviations result from the effects of Galactic tides. To this purpose, we use N-body simulations of the evolution of CDM halos in a Galactic potential, constrained to have the orbital parameters consistent with a dwarf's present-day position and velocity. We shall assume that, over the past 10 Gyr, the Galactic potential is the static, analytic potential inferred by @mcmillan2011 from observations of kinematic tracers. We also assume that the potential of each dwarf may be initially approximated by a cuspy NFW profile. Since the dwarfs in question are heavily dark-matter-dominated, we shall use a carefully selected sample of dark matter particles to emulate the evolution of an embedded tracer stellar component. In this Chapter, we describe our choice of Galactic potential, orbital estimation, initial conditions setup, and N-body methods. 



# Orbital estimation {#sec:orbital_estimation}

To explore the possible orbits of a dwarf galaxy, we perform a Monte Carlo sampling of the present-day observables. The present-day position, distance modulus, LOS velocity, and proper motions are each sampled from normal distributions given the reported uncertainties in [@tbl:scl_obs_props; @tbl:umi_obs_props]. We integrate each sampled position/velocity back in time for 10 Gyr using \agama{} [@agama]. Dynamical friction is not expected to impact orbits substantially because of the low masses and large pericentres of the dwarfs, so we assume a single point-mass particle for the backwards integration.

## Galactocentric frame

To convert observed positions and velocities to Galactocentric coordinates, we use the Astropy v4 Galactocentric frame [@astropycollaboration+2022]. Our Cartesian Galactocentric coordinates here assume the Galactic centre is at $[x, y, z] = [0,0,0]$, where $x$ is the direction from the sun to the Galactic centre, $y$ is the direction of the motion of the Local Standard of Rest, and $z$ is the direction perpendicular to the Galactic plane. The coordinate frame is also right-handed, such that the $z$-angular momentum of the sun is negative (since the sun is at $x<0$). In this frame, the solar position is $[-8.122 \pm 0.033, 0, 0.0208 \pm 0.003]\, \kpc$ [@gravitycollaboration+2018; @bennett+bovy2019] and the solar velocity is $\V_\odot = [-12.9 \pm 3.0, 245.6 \pm 1.4, 7.78 \pm 0.08]$ km/s [@reid+brunthaler2004; @drimmel+poggio2018; @gravitycollaboration+2018]. The uncertainties in this reference frame are typically smaller than the uncertainties on a dwarf galaxy's distance and tangential velocity. 

## Milky Way potential

We adopt the Milky Way potential described in @EP2020, which is an analytic approximation to that proposed by @mcmillan2011. @fig:v_circ_potential plots the circular velocity profiles of each component and the total circular velocity profile for this potential. The potential includes a stellar bulge, a thin and thick disk, and a dark matter NFW halo. 

The Galactic bulge is described by a @hernquist1990 potential,

$$
\Phi(r) = - \frac{GM}{r + a},
$$
where  $a=1.3\,{\rm kpc}$ is the scale radius and $M=2.1 \times 10^{10}\,\Mo$ is the total mass. The thin and thick disks are represented with the @miyamoto+nagai1975  cylindrical potential,
$$
\Phi(R, z) = \frac{-GM}{\left(R^2 + \left[a + \sqrt{z^2 + b^2}\right]^{2}\right)^{1/2}},
$$
where $a$ is the disc radial scale length, $b$ is the scale height, and $M$ is the total mass of the disk. For the thin disk,  $a=3.944\,$kpc, $b=0.311\,$kpc, and $M=3.944\times10^{10}\,$M$_\odot$. For the thick disk, $a=4.4\,$kpc, $b=0.92\,$kpc, and $M=2\times10^{10}\,$M$_\odot$. The halo is an NFW dark matter halo ([@eq:nfw]) with $\rmax = 43.7\,$kpc and $\vmax = 191\,\kms$, or $M_{200} = 126.6\times 10^{10}\,\Mo$ and $r_s=20.2\,\kpc$.



![Circular velocity of the Milky Way potential](figures/v_circ_potential.pdf){#fig:v_circ_potential} 

Figure: Circular velocity profile of @EP2020 potential. The total circular velocity (thick black line) is composed of an NFW halo (green dashed line), thin and thick @miyamoto+nagai1975 disks (orange dash-dotted line), and a @hernquist1990 bulge (blue dotted line).

## Sculptor's orbit {#sec:scl_smallperi}

Sculptor's orbit in the assumed potential is relatively well-constrained. [@fig:scl_orbits] illustrates point particle orbits for 100 samples of Sculptor's observed kinematics integrated backwards for $10\,\Gyr$ in both Galactocentric Cartesian coordinates ($x$, $y$, $z$) and in Galactocentric radius with time. These orbits sample the uncertainties in distance, proper motion, and radial velocity, as given in [@tbl:scl_obs_props]. All sampled orbits of Sculptor have nearly the same shape---the orbit primarily resides in the $y$–$z$ plane and completes a similar number of periods and pericentric and apocentric passages. 

To maximize tidal effects, we select an orbit with the $\sim 3\sigma$ smallest pericentre among all possible orbits integrated backwards for $10\,\Gyr$. We achieve this by taking the median parameters of all orbits with a pericentre less than the 0.0027th quantile pericentre, yielding a pericentre of 43 kpc. Given the current observations, it is unlikely that Sculptor has had a significantly smaller pericentre than our selected orbit, which we refer to as the \smallperi{} orbit. We take the first apocentre after a look-back time of 10 Gyr, or at $\sim9.1\,\Gyr$, as the initial conditions for our model of Sculptor, noted in @tbl:orbit_ics.





![Sculptor's possible orbits](/Users/daniel/thesis/figures/scl_xyzr_orbits.pdf){#fig:scl_orbits}

Figure: The orbits of Sculptor in a static Milky Way potential in Galactocentric $x$, $y$, and $z$ coordinates (top) and in Galactocentric radius $r$ versus time (bottom). The Milky Way is at the centre with the disk lying in the $x$--$y$ plane. Our selected \smallperi{} orbit is plotted in black, and light green transparent lines represent orbits sampled over Sculptor observables in [@tbl:scl_obs_props]. The orbit of Sculptor is well-constrained in this potential, and it is unlikely to achieve a smaller pericentre than the \smallperi{} orbit. 



![Ursa Minor's possible orbits](/Users/daniel/thesis/figures/umi_xyzr_orbits.pdf){#fig:umi_orbits}

Figure: Similar to @fig:scl_orbits, the orbits of Ursa Minor in a static Milky Way potential in Galactocentric $x$, $y$, and $z$ coordinates. In the lower panel, we show the radius versus time for only three orbits of Ursa Minor, representing the \smallperi{} point-mass orbit (black), the mean orbit, and the orbit with the $3\sigma$-largest pericentre. 



| Property              | Scl: \smallperi{} | Scl: `LMC-flyby` | Umi: \smallperi{} |
| --------------------- | ----------------- | ---------------- | ----------------- |
| distance / kpc        | 82.6              | 72.9             | 64.6              |
| $\pmra / \masyr$      | 0.134             | 0.137            | -0.158            |
| $\pmdec / \masyr$     | -0.198            | -0.157           | 0.050             |
| LOS velocity / $\kms$ | 111.2             | 111.2            | -245.75           |
| $t_i / \Gyr$          | -9.17             | -2.0             | -9.67             |
| ${x}_{i} / \kpc$      | -2.49             | 4.30             | -17.40            |
| ${y}_{i} / \kpc$      | -42.78            | 138.89           | 74.51             |
| ${z}_{i} / \kpc$      | 86.10             | 125.88           | 21.34             |
| $\V_{x\,i} / \kms$    | -20.56            | 6.89             | 14.27             |
| $\V_{y\,i} / \kms$    | -114.83           | -56.83           | 48.62             |
| $\V_{z\,i} / \kms$    | -57.29            | 52.09            | -114.08           |

Table: The orbital initial conditions for models presented in this work. The observables represent the medians from orbital integration used to derive the orbits. Instead, the initial position and velocity represent the initialization of the actual N-body model. The \smallperi{} represents instead the $3\sigma$ smallest pericentre, which we use to provide an upper limit on tidal effects. We describe the \texttt{LMC-flyby} orbit in @sec:scl_lmc. {#tbl:orbit_ics  short="Orbit initial conditions"}

## Ursa Minor's orbit {#sec:orbit_corrections}

Similar to Sculptor, Ursa Minor has a well-constrained orbit in the assumed MW potential. @fig:umi_orbits shows 100 random point-mass orbits of Ursa Minor. As for Sculptor, we select an orbit with approximately the $3\sigma$ smallest pericentre, the Ursa Minor $\smallperi{}$ orbit (see @fig:umi_orbits and @tbl:orbit_ics). 

We shall see later that the orbit of Ursa Minor's N-body model differs from the point particle orbit because of the effects of tidal mass loss. To ensure the final conditions of the N-body simulation are close to the intended final position, we iteratively adjust Ursa Minor's initial conditions. Initially, starting with low-resolution runs, we adjust the cylindrical actions of the initial orbit by the final difference in actions at the end of orbital evolution. After the initial actions have converged (2 iterations), we change the initial action angles by the final difference in action angles. This method converges within 4 iterations to an orbit agreeing with the observed kinematics of Ursa Minor. Since Sculptor's orbit is less strongly affected by tides, we do not carry out this correction for Sculptor.

# Initial conditions

We use \agama{} [@agama] to generate the initial N-body dark matter halo. We assume galaxies are described by a spherical, isotropic NFW dark matter potential ([@eq:nfw]). We also assume the stars do not contribute to the potential. The dark matter density is truncated in the outer regions by
$$
\rho_{\rm tNFW}(r) = e^{-(r/r_t)^3}\ \rho_{\rm NFW}(r),
$$ {#eq:trunc_nfw}
where we adopt $r_t = 20 r_s$. 

## Initial dark matter halos for Sculptor and Ursa Minor

The observed half-light radius and velocity dispersion of Sculptor and Ursa Minor, together with the mass-concentration relation of \LCDM{} halos, determine the properties of the N-body models adopted for each dwarf.

[@tbl:derived_props] lists our inferred halo and kinematic properties of Sculptor and Ursa Minor. First, taking the absolute magnitudes from @munoz+2018 with the mass-to-light ratio from @woo+courteau+dekel2008, the total current stellar mass of Sculptor and Ursa Minor are $M_\star \approx 3.1 \times 10^6 \Mo$ and $M_\star \approx 7 \times 10^5 \Mo$, respectively. Based on the stellar mass-$\vmax$ relation [from @fattahi+2018], Sculptor and Ursa Minor's halos should have $\vmax \approx 31 \,\kms$ and $\vmax \approx 27\,\kms$. Finally, using the @ludlow+2016 $z=0$ mass-concentration relation, this constraint translates into a $\rmax \approx 6 {\rm kpc}$ and $\rmax \approx 5\,\kpc$ for each galaxy. 

The observed velocity dispersion constrains the total mass within the stellar half-light radius.  From the @wolf+2010 mass estimator, the mass contained within the 3D half-light radius $r_h$ is
$$
M(r_h) \approx \frac{3}{G} \sigma_\V^2\,r_h
$$ {#eq:wolf_mass}
for a velocity dispersion $\sigma_\V$. For Scl and UMi, the cosmological mean parameters above would predict $\sigma_{\V, i}  \approx 8.5\,\kms$ and $8.0\,\kms$ assuming $r_h=240\,$pc. This is below the observed values of $9.7\,\kms$ and $8.7\,\kms$, and tidal evolution will further reduce the dispersion. To match the observed $\sigma_\V$ at the end of the simulation, we choose $\vmax = 31\,\kms$ and $\rmax=3.2\,\kpc$ for Scl and $\vmax=38\,\kms$ and $\rmax=4\,\kpc$ for UMi.

While there is some range in the choice of initial halo, reasonable changes to the initial halo do not substantially affect the tidal evolution for either galaxy.  So, the tidal effects on stars should be similar for halos with similar velocity dispersions. 



| parameter           | Sculptor                                      | Ursa Minor                              |
| ------------------- | --------------------------------------------- | --------------------------------------- |
| $L_\star$           | $1.8\pm0.2\times10^6\ L_\odot$                | $3.5 \pm 0.1 \times 10^5\,L_\odot$      |
| $M_\star$           | $3.1_{-1.0}^{+1.6} \times10^6\ {\rm M}_\odot$ | $7_{-2}^{+3} \times 10^5\,\Mo$          |
| $M_\star / L_\star$ | $1.7\times 10^{\pm 0.17}$                     | $1.9 \times 10^{\pm 0.17}$              |
| $\vmax$             | $31\pm 3\,\kms$                               | $27_{-2}^{+3}\,\kms$                    |
| $\rmax$             | $6 \pm 2$ kpc                                 | $5_{-1}^{+2}$ kpc                       |
| $M_{200}$           | $0.5 \pm 0.2\times10^{10}\ M_0$               | $0.3_{-0.1}^{+0.2} \times 10^{10}\,\Mo$ |
| $c_{\rm NFW}$       | $13_{-3}^{+4}$                                | $13.5_{-3}^{+4}$                        |

Table: Inferred properties of the stellar component and halo for Sculptor and Ursa Minor. We record the total luminosity, stellar mass, mass-to-light ratio, dark matter halo $\vmax$ and $\rmax$, and dark matter halo virial mass $M_{200}$ and concentration $c_{\rm NFW}$. Uncertainties are the 16-84th percentile derived using Monte-Carlo sampling, assuming 0.035 and 0.1 dex uncertainties in the $\vmax(M_\star)$ and $c(M_{200})$ relations.  {#tbl:derived_props  short="Derived Properties of Sculptor and Ursa Minor"}



![Initial halo parameter choice](figures/initial_halo_selection.pdf){#fig:scl_halos}

Figure: The blue circle indicates our choice of initial halo parameters for Sculptor and Ursa Minor. The grey line and green line with shaded regions represent the @ludlow+2016 mass-concentration relation and @fattahi+2018 SMHM relation, respectively. The curved orange lines represent the velocity dispersion of an exponential stellar component with the observed half-light radius of Sculptor and Ursa Minor.

# Numerical methods

## The N-body code: \gadget{}

To simulate the tidal evolution of galaxies, we use N-body simulations integrated with the parallel, gravitational-tree program \gadget{} [@gadget4]. The N-body method calculates and evolves the gravitational accelerations between a large number of collisionless particles to approximate the dynamical evolution of matter. To approximate a collisionless system (i.e., without strong particle-particle gravitational deflections), the force is tapered below a "softening length." Resolution is limited by the number of particles and the softening length. Tree codes organize particles into a spatial tree, enabling the grouping of the gravitational forces from nearby particles. A gravitational tree code substantially reduces the required number of force computations, as compared to the exact, direct-summation method. 

## Isolation runs and simulation parameters

To ensure that the initial conditions of the simulation are dynamically relaxed and well-converged, we run a halo first in isolation using \gadget{}. Since gravity is scale-free, we use the same isolation run for all halos and rescale the results to the desired values of size and mass. We adopt a fiducial value of  $\rmax = 6.0\,$kpc and $\vmax = 31\,\kms$  for the isolation halo based on Sculptor's mean properties. We run this model for 5 Gyr (about one-half crossing time $t_{\rm cross} = 2\pi\,r /\vcirc  \approx 9\,\Gyr$ at $r_{200}=36\,$kpc).

For our simulation parameters, we adopt a softening length of 
$$
h_{\rm grav} = 0.014\,{\rm kpc}\left(\frac{r_{\rm max}}{6.0\,{\rm kpc}}\right)\left(\frac{N}{10^7}\right)^{-1/2},
$$ {#eq:softening_length}
where the force is softened by a spline kernel [eqs. 70--71 in @springel+yoshida+white2001; but with the spline characteristic radius $2.8h_{\rm grav}$, see @gadget4]. See Appendix -@sec:extra_convergence for a discussion of this choice and our simulation parameters, which is similar to the @power+2003 suggested softening.

## Numerical fidelity

@fig:numerical_convergence illustrates how well our numerical setup is able to reproduce the desired initial conditions, before and after running the model in isolation. This figure shows that our numerical methods are able to approximate well an NFW halo down to a given innermost radius that strongly depends on resolution. The larger the number of particles, the smaller the radius that is effectively "resolved" in a given simulation. For the Sculptor-like halo shown in this figure (with $\rmax = 6.0\,$kpc and $\vmax = 31\,\kms$), a simulation with $10^7$ particles is needed to resolve the innermost 100 pc. For reference, the half-light radius of Sculptor is roughly 100 pc, which means that at least 10 million particles are needed to follow faithfully its tidal evolution. Vertical arrows in @fig:numerical_convergence indicate the "convergence radius" defined by @power+2003 [their eq. 13] for NFW halos formed in cosmological N-body simulations. This radius marks the region where collisional effects driven by the finite number of particles used to describe the innermost regions of a halo become important over a Hubble time in a cosmological simulation. The softening length (from @eq:softening_length) is typically a few times smaller than the converged length.

![Numerical convergence of the N-body simulation](figures/iso_converg_num.pdf){#fig:numerical_convergence}

Figure: Numerical convergence test for circular velocity as a function of log radius for simulations with different total numbers of particles in isolation. Residuals in the bottom panel are relative to NFW. The initial conditions are dotted, the converged radius is marked by arrows [eq. 13, @power+2003], and the softening length is marked by a vertical bar. Note that a slight reduction in density starting around $r = 30\,\kpc$ is expected given our outer truncation choice. 

## Orbital evolution

Next, we evolve the halo in the Galactic potential. We scale the relaxed snapshot and softening length to match the initial halo, and shift the snapshot to the initial conditions inferred from the orbital analysis (see [@tbl:orbit_ics]). We then evolve the full N-body NFW model forward in time in the Galactic potential until the present time, when the halo is closest to its present-day observed position in the MW halo. 

## Halo centring {#sec:shrinking_spheres}

To accurately follow the evolution of a halo, it is important to determine the centre of the self-bound halo remnant at each time chosen for analysis. We use a shrinking-spheres centre method inspired by @power+2003. First, we start with an initial centre estimate using all bound particles from the previous snapshot. Then, we calculate the distance of all particles from the centre, remove particles with a distance beyond the 0.975 quantile of the centre, and recalculate the centre of mass. The procedure is repeated until the selection radius is less than ~1kpc or fewer than 0.1% of particles remain. After a centre has been chosen, we remove all unbound particles based on the \gadget{} calculated potential of the halo. For all future timesteps, we consider only particles retained from the previous iteration.

The statistical centring uncertainty for the full resolution ($10^7$ particle) isolation run is of order 0.003 kpc, but fluctuations are observed of order  $\sim 0.03\,\kpc$. This is about three times the softening length but is less than the numerically converged radial scale. 



## Sculptor and Ursa Minor's initial stellar components {#sec:painting_stars}

We "paint" stars onto dark matter particles using the particle-tagging method [e.g., @bullock+johnston2005], assuming spherical symmetry. We initially assume stars follow a projected exponential law ([@eq:exponential_law]) with $R_s = 0.10\,\kpc$ for both galaxies.  The tagging method assigns a probability to each dark matter particle, which is proportional to the "light-to-mass" ratio required to match the assumed stellar light profile. We briefly describe the procedure next, but refer interested readers to @EP2020.

Let $\Psi$ be the potential (normalized to vanish at infinity) and ${\cal E}$ the binding energy, ${\cal E} = \Psi - 1/2 v^2$. If we know the distribution function $f({\cal E})$ (the phase-space density as a function of energy), then we assign a stellar weight to a given particle with energy ${\cal E}$ using
$$
P_\star({\cal E}) = \frac{f_\star({\cal E})}{f_{\rm halo}({\cal E})}.
$$
While $f({\cal E})$ is a phase-space density, the differential energy distribution includes an additional $g({\cal E})$ occupation term [@BT1987]. We use Eddington inversion to find the distribution function, [eq. 4-140b in @BT1987]
$$
f({\cal E}) = \frac{1}{\sqrt{8}\, \pi^2}\left( \int_0^{\cal E} \frac{d^2\rho}{d\Psi^2} \frac{1}{\sqrt{{\cal E} - \Psi}}\ d\Psi + \frac{1}{\sqrt{\cal E}} \left(\frac{d\rho}{d\Psi}\right)_{\Psi=0} \right).
$$
In practice the right, boundary term is zero.^[If, at large $r$, $\rho \propto r^{-n}$ with $n>1$, so $\Psi \sim r^{-1}$ then $d\rho / d\Psi \sim r^{-n+1}$.] We take $\Psi$ from the underlying assumed analytic dark matter potential. $\rho_\star$ can be calculated from the surface density, $\Sigma_\star$, via the inverse Abel transform. 

@fig:scl_umi_initial_isolation shows the initial circular velocity profiles of stars and dark matter, demonstrating that these galaxies are stable in isolation and that the stellar component contributes negligibly to the total mass ($M \sim \V_{\rm circ}^2 R$). Our cosmologically-motivated, observationally-based, and numerically-converged initial conditions enable us to accurately consider tidal effects in the next Chapter.



![Initial halo velocity profiles](figures/initial_velocity.pdf){#fig:scl_umi_initial_isolation}

Figure: The initial circular velocity profiles of dark matter (blue) and stars (orange) for Sculptor and Ursa Minor. The initial conditions are dotted, and the isolation-evolved profiles are solid. The green cross marks the present-day half-light radius and velocity dispersion, and the black band represents the 1-sigma mean density of the MW at pericentre across orbits. Initial conditions are stable in isolation, and mass is dominated by dark matter. 

<!-- JFN: indicate jacobi radius in above figure, maybe just intersection of lines?-->
