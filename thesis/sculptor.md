In this chapter, we present our results from our simulations of Sculptor and Ursa Minor. In each case, we first discuss the range of reasonable orbits and initial dark matter halos for each galaxy. Then, we describe how tides affect the evolution of each galaxy. We furthermore consider how the effects of the LMC changes the results, including using additional N-body simulations in evolving potentials. Finally, we review the predictions and characteristics of the present-day stellar populations and properties from our models.

# Sculptor

Before running N-body simulations, we need to understand both the typical halos which could host a sculptor like galaxy and the possible orbits Sculptor.





@fig:scl_halos presents approximate constraints based on 





### Derived properties and setup

| parameter            | value                                         | reference                                                    |
| -------------------- | --------------------------------------------- | ------------------------------------------------------------ |
| $L_\star$            | $1.8\pm0.2\times10^6\ L_\odot$                |                                                              |
| $M_\star$            | $3.1_{-1.0}^{+1.6} \times10^6\ {\rm M}_\odot$ |                                                              |
| $M_\star / L_\star$  | $1.7\times 10^{\pm 0.17}$                     | @woo+courteau+dekel2008                                      |
| $M_{200}$            | $0.5 \pm 0.2\times10^{10}\ M_0$               |                                                              |
| $c_{\rm NFW}$        | $13.1_{-2.8}^{+3.6}$                          |                                                              |
| $v_{\rm circ, max}$  | $31\pm 3\,\kms$                               | relation from @fattahi+2018                                  |
| $r_{\rm circ, max}$  | $6 \pm 2$ kpc                                 | $v_{\rm circ}$ with @ludlow+2016 mass-concentration relation |
| $r_{\rm break, obs}$ | $25 \pm 5$ arcmin                             | @sestito+2024a                                               |
| $t_{\rm break, obs}$ | $110\pm30$ Myr                                | $\sigma_v$, $r_{\rm break}$ with @                           |

Table: To derive total mass, we use the absolute magnitude from @munoz+2018 (see also [@tbl:scl_obs_props]) along with the stellar mass to light ratio (1.7 with 0.17 dex uncertainty) from @woo+courteau+dekel2008. Sculptor's total stellar mass is then $M_\star \sim 3.1_{-1.0}^{1.6} \times 10^6\,\Mo$. From the stellar mass to dark matter halo's characteristic velocity $v_{\rm circ}$ in @fattahi+2018, we expect Sculptor to have $v_{\rm circ} \approx 31 \pm 3 \kms$. {#tbl:scl_derived_props  short="Derived Properties of Sculptor"}

For the mass-to-light ratio, we use the values from @woo+courteau+dekel2008. Note that we assume a scatter of 0.17 dex for the mass-to-light ratios as reported for mass uncertainties. @delosreyes+2024 also indicate that other methods for stellar mass-to-light ratios result in similar scatters and systematic biases. However these values are just approximate. 

![Sculptor initial halos](figures/scl_initial_halos.pdf){#fig:scl_halos}

Figure: The suggested halos of Sculptor compared to cosmological predictions.



| Halo name | Rmax | Vmax | M200 | c    |
| --------- | ---- | ---- | ---- | ---- |
| compact   | 3.2  | 31   | 0.33 | 21   |
| lmc       | 4.2  | 31   | 0.39 | 17   |
| small     | 2.5  | 25   |      |      |

Table: The parameters for our initial Sculptor halos. {#tbl:scl_ini_halos  short="Initial halos of Sculptor"}

## Milky Way tides

### Orbital properties

![Sculptor Orbits](figures/scl_xyzr_orbits.pdf){#fig:scl_orbits}

Figure: The orbits of Sculptor in a static Milky Way potential in galactocentric $x$, $y$, and $z$ coordinates. The Milky Way is at the centre with the disk lying in the $x$--$y$ plane. Our selected `smallperi` orbit is plotted in black and light blue transparent orbits represent the past 5Gyr orbits sampled over Sculptor observables in [@tbl:scl_obs_props]. The orbit of sculptor is well-constrained in this potential and it is unlikely to achieve a smaller pericentre than our selected orbit.  \the\textwidth



To select an orbit with ~ maximum possible observationally-consistent tidal forces, we take the median parameters for orbits with a pericentre less than 2x the 3$\sigma$ minimum pericentre and derive the following orbit. The initial conditions are from the 



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

Table: Properties of selected orbits for Sculptor. The mean orbit represents the observational mean from @tbl:scl_obs_props. The Smallperi represents instead the $3\sigma$ smallest pericentre, which we use to provide an upper limit on tidal effects. {#tbl:scl_obrits  short="Sculptor Selected Orbits"}



## Tidal effects

- dark matter loss
- little influence to stars
- smallest pericentre does not help

![Sculptor simulation snapshots](figures/scl_sim_images.png)

Figure: Images of the dark matter evolution over a selection of past apocentres and the present day position. Limits range from -150 to 150 kpc in the $y$-$z$ (approximately orbital) plane and the colourscale is logarithmic spanning 5 orders of magnitude between the maximum and minimum values. In this image, stars occupy only ever a few pixels so are not plotted. 



![Sculptor Tidal Tracks](figures/scl_tidal_track.pdf)

Figure: The tidal tracks for the smallperi orbit. Todo: add velocity dispersion plot to RHS



![Sculptor velocity dispersion evolution](figures/scl_sigma_v_time.pdf)

Figure: Evolution of stellar velocity dispersion within 1 kpc for different Scl models. In all cases, the evolution is mild. Note that binarity may reduce the inflate the observed velocity dispersion by  ~ 1 km/s, so the conservative lower limit is around 8 km/s.





![Sculptor initial and final density profiles](figures/scl_smallperi_i_f.pdf)

Figure: Effects on exponential initial stars. TODO: plot 2D sky proj. stars





![Sculptor Plummer initial and final density profiles](figures/scl_plummer_i_f.pdf)

Figure: effects on Plummer initial stars.

## Effects of the LMC



### Orbital effects

As discussed in @battaglia+2022, Sculptor's orbit is strongly influenced by the presence of an LMC. Figure @fig:scl_lmc_orbit_effect, the addition of an LMC reduces Scl's pericentre with the MW and implies that Scl may be on its first infall with the MW, in contrast to the discussion above. Thus, the LMC has a critical impact on the evolution of Scl. Additionally, Scl passes rapidly but relatively close to the LMC in the past 100 Myr. 



Finally, we also consider the influence of the large Milky Cloud on the orbits and evolution of Scl. We adopt the @vasiliev2024 multipole approximation of an N-body simulation of the LMC and MW. Their initial conditions are

- MW halo:
- MW bulge (static):
- MW disk (static):
- LMC halo:

We now consider 3 orbits under the L3M11 potential model from @vasiliev2024. We focus on this LMC potential as a lighter LMC or MW should only reduce tidal impacts and the recent orbit of Scl is minimally affected by the LMC potential. 

![Sculptor Orbits with LMC](figures/scl_lmc_xyzr_orbits.pdf){#fig:scl_lmc_orbit_effect}

Figure: This figure is similar to @fig:scl_orbits except that we are showing the orbits with and without an LMC. In the bottom row, the distance from Sculptor (or the LMC) to the MW is plotted (left), and the Sculptor - LMC distance (right.)



**tinyperilmc**

- ra = 15.0183
  dec = -33.7186
  distance = 73.1
  pmra = 0.137
  pmdec = -0.156
  radial_velocity = 111.2

  t_i = -2.00
  pericentre = 38.82
  apocentre = 187.50
  t last peri = -0.33
  x_i = [4.30 138.89 125.88]
  v_i = [6.89 -56.83 52.09]

Smallperi

- ra = 15.0183
  dec = -33.7186
  distance = 84.0
  pmra = 0.166
  pmdec = -0.237
  radial_velocity = 111.2

  t_i = -2.00
  pericentre = 28.00
  apocentre = 190.45
  t last peri = -0.39
  x_i = [12.79 168.76 87.33]
  v_i = [0.37 -57.69 61.40]



### Tidal effects

![Sculptor Simulation Snapshots with LMC](figures/scl_lmc_sim_images.pdf)





![Sculptor initial and final density with LMC](figures/scl_lmc_i_f.pdf)

Figure: The tidal effects on the stellar surface density due to the LMC today.

# Ursa Minor

## Derived properties

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

![Ursa Minor initial halos](figures/umi_initial_halos.pdf)

Figure: The suggested halos of Ursa Minor compared to cosmological predictions.

## Orbital properties

![Ursa Minor Orbits](figures/umi_xyzr_orbits.pdf)

## Tidal Effects

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





![Ursa Minor simulation snapshots](figures/umi_sim_images.png)

Figure: Ursa Minor simulation images. 



Figure: Velocity dispersion evolution of Ursa Minor



![Ursa Minor simulated density profiles](figures/umi_smallperi_i_f.pdf)

Figure: The tidal effects on the stellar surface density.

## Effects of the LMC



![Ursa Minor orbits with LMC](figures/umi_lmc_xyzr_orbits.pdf)

Figure: Orbits of Ursa Minor with (orange) and without (green) an LMC. The final positions of Ursa Minor and the LMC are plotted as scatter points and the solid blue line represents the LMC trajectory. Note that the LMC only increases Ursa Minor's peri and apo-centres, weakening any tidal effect. Interestingly, there is a change that Ursa Minor may have once been bound to the LMC (diverging orange lines at top left of middle panel.)

# Summary