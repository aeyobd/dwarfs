As discussed in Chapters [-@sec:introduction; -@sec:data], Sculptor and Ursa Minor have unusually extended stellar light profiles. Using our setup introduced in Chapter [-@sec:methods], in this Chapter, we present the results of N-body simulations tailored to each galaxy. We find that tides do induce moderate dark matter loss. However, the compact stellar distribution for each galaxy is nearly unaffected. We also show that the LMC strongly perturbs Scl's orbit. But, the tidal effects when considering an LMC are yet weaker than in the Milky Way only case. We conclude that tides are unlikely to affect the stellar component of either galaxy.



# Tidal effects on Sculptor



The Milky Way tides drive significant mass loss for Sculptor. @fig:scl_sim_images shows projected dark matter distribution of our N-body simulation for Sculptor in the $y$-$z$ plane for 5 previous pericentres. This simulation uses the fiducial initial halo on the \smallperi{} orbit. The dark matter forms large streams which orbit the Milky Way several times. Note that, in contrast to the outer regions, the central dark matter cusp is similar in appearance across all simulation snapshots. Despite the mass loss, Sculptor's orbit only deviates slightly from a point particle orbit, reaching a similar position to the observed location in @fig:scl_sim_images.

@fig:scl_tidal_track shows the initial and final circular velocity profile of the dark matter, and the evolution of the maximum circular velocity. Sculptor's final dark matter halo evolves from $\vmax=31\,\kms$ at $\rmax=3.2\,\kms$ to $\vmax = 22\,\kms$ at $\rmax = 1.3\,\kpc$. In terms of bound mass, the mass evolves from  $0.40 \times 10^{10}\,\Mo$ to $0.026\times10^{10}\,\Mo$. The circular velocity maximum evolves along the tidal track predicted by @EN2021. 



![Sculptor simulation snapshots](figures/scl_sim_images.png){#fig:scl_sim_images}

Figure: Images of the dark matter evolution over a selection of past apocentres and the present day position. Limits range from -150 to 150 kpc in the $y$-$z$ (approximately orbital) plane and the colourscale is logarithmic spanning 5 orders of magnitude between the maximum and minimum values. The green dot represents the final (observed) position of the galaxy and the grey curve represents the orbit over one radial oscillation.





| Property                  | smallperi             | smallperi-Plummer | LMC     |
| ------------------------- | --------------------- | ----------------- | ------- |
| final heliocen distance   | 81.6 kpc              | --                |         |
| $\vmax_f / \vmax_i$       | 0.695                 | --                | 1.0     |
| $\rmax_f / \rmax_i$       | 0.406                 | --                | 1.0     |
| fractional dm mass loss   | 0.0893                | --                | 1.00    |
| $\sigma_v i$              | 9.8                   | 11?               |         |
| $\sigma_vf$               | 8.8                   | ?8                |         |
| fractional stellar mass   | $1 -2.1\times10^{-6}$ | $1 - x$           | $1 - y$ |
| r_h_stars_f / r_h_stars_i | 1.12                  |                   |         |
| Jacobi radius             | 148 arcmin, 3.5 kpc   |                   |         |
| Break radius              | 110 arcmin, 2.6 kpc   |                   |         |

Table: The present-day properties for the three simulations of Sculptor we consider here.



![Sculptor Tidal Tracks](figures/scl_tidal_track.png){#fig:scl_tidal_track}

Figure: The tidal tracks for the smallperi orbit. **Todo**: add velocity dispersion plot to RHS?

While mass loss leads to expansion of the stellar component, no features of tidal disruption or disequilibrium would be currently observable. @fig:scl_smallperi_i_f shows the final projected density of stars on the sky, and the initial and final radially averaged density profiles. No non-spherical density features are apparent, even at 5 decades fainter in surface density than the central density. In addition, the initial and final density profiles look identical up to some scale, with a tidally-induced excess of stars only appearing $\sim 3$ orders of magnitude fainter than the faintest our density profile measures. 

The analytically motivated break and tidal radii also support a weak tidal effect on stars. In both cases, these radii work out to be $\gtrsim 100$ arcminutes, outside where we can measure the density profile. In particular, because we have chosen the most extreme observationally permissible orbit, it is unlikely, for similar $\Lambda$CDM motivated initial conditions, that these radii would ever approach where the break in Sculptor's density profile appears. 



![Sculptor initial and final density profiles](figures/scl_smallperi_i_f.pdf){#fig:scl_smallperi_i_f}

Figure: The tidal effects on Scl's stellar component, for the \smallperi{} orbit with the fiducial halo and exponential stars with $R_s=0.10\,\kpc$. **Left:** the final 2D projected density of stars on the sky, with colours representing logarithmically increasing density. The circle marks the break radius. **Right:** The initial (dotted) and final (solid) stellar density profiles as compared to the observed stellar density profile. Arrows mark the break and Jacobi ([@eq:jacobi_radius]) radii. **sculptor in figure title, orbit trace opaque, monotonic colours?**





![Sculptor Plummer initial and final density profiles](figures/scl_plummer_i_f.pdf){#fig:scl_smallperi_plummer_i_f}

Figure: Similar to [@fig:scl_smallperi_i_f] except for Plummer initial stars with $R_h = 0.20\,\kpc$. While a faint stream may be visible with deeper observations, effects on the stellar profile are minimal.







## Effects of the LMC {#sec:scl_lmc}

The Milky Way isn't the only galaxy in town. Recently, work has shown that the infall of the LMC may substantially affects the Milky Way system [e.g., @erkal+2019; @cautun+2019; @garavito-camargo+2021; @vasiliev2023]. The mass of the LMC may be as high as 1/5 of the Milky Way mass. Including an LMC in the Milky Way potential may change conclusions about properties and orbits of the Milky Way satellite system, including Scl [e.g., @patel+2020; @battaglia+2022].

### How the LMC changes Sculptor's orbit

To explore the effects of the LMC on Sculptor, we use the publicly available L3M11 model from @vasiliev2024. While @vasiliev2024 present other models, the L3M11 model has the highest MW and LMC masses, and the (recent) orbital history appears to be independent of the potential details. This model stars with a lighter Milky Way than our fiducial, @mcmillan2011-like model. The L3M11 potential is an evolving multipole approximation of an N-body simulation with both a live MW and LMC dark matter halo. The potential additionally includes the reflex motion of the MW due to the LMC. The MW was initialized as an NFW dark matter halo with $r_s=16.5\,$kpc and mass $M_{\rm 100}= 11\times10^{11}\Mo$, and the LMC was a NFW halo with $r_s=11.7$ and $M_{100} = 2.76 \times 10^{11} \Mo$. This model had a previous LMC pericentre at about 6 Gyr ago.

For our orbit, we select the orbit with the $3\sigma$-smallest LMC pericentre. Results are similar when selecting for a small MW pericentre instead. We also use, instead, the `small` halo for Sculptor, as the tidal effects do not reduce the velocity dispersion, unlike for the MW-only model.

@fig:scl_lmc_orbits_effect displays the effect of including an LMC in the potential.   The green samples are in the initial MW-only potential in the `L3M11` model, and the orange samples are integrated in the evolving MW and LMC `L3M11` model. The past 1 Gyr is similar in both cases, but the orbits diverge significantly afterwards. The recent passage of Sculptor with the LMC around 0.1 Gyr ago allows for Sculptor to begin as far as 300 kpc from the Milky Way centre. The evolving potential also adds significant long term variability in the possible orbits of Sculptor. Sculptor, however, is orbiting in the opposite direction of the LMC so is likely not associated with the LMC system. 

Given the large uncertainties of the LMC model, we conservatively double all of the observational parameters of Sculptor. This has a similar effect to including LMC mass and orbital uncertainties but is considerably simpler. From this larger range of orbits, we once again select the orbit with the median final observables of all orbits with pericentres less than the 0.0027th quantile. This orbit, the `smallperilmc` orbit is plotted in black and is only integrated up to 2 Gyr ago, to isolate recent tidal effects. 



![Sculptor Orbits with LMC](figures/scl_lmc_xyzr_orbits.png){#fig:scl_lmc_orbits_effect}

Figure: This figure is similar to @fig:scl_orbits except that we are showing the orbits with and without an LMC. In the bottom row, the distance from Sculptor (or the LMC) to the MW is plotted (left), and the Sculptor - LMC distance (right.) 





| Property                         | smallperi-LMC |
| -------------------------------- | ------------- |
| distance / kpc                   | 73.1          |
| $\pmra / \masyr$                 | 0.137         |
| $\pmdec / \masyr$                | --0.156       |
| LOS velocity / $\kms$            | 111.2         |
| $t_i / \Gyr$                     | -2            |
| ${x}_{i} / \kpc$                 | 4.30          |
| ${y}_{i} / \kpc$                 | 138.89        |
| ${z}_{i} / \kpc$                 | 125.88        |
| $\V_{x\,i} / \kms$               | 6.89          |
| $\V_{y\,i} / \kms$               | -56.83        |
| $\V_{z\,i} / \kms$               | 52.09         |
| pericentre / kpc                 | 38.82         |
| peri-LMC / kpc                   | ?             |
| apocentre / kpc                  | 187.50        |
| $t_{\rm last\ peri} / {\rm Gyr}$ | -0.46         |
| $t_{\rm last\ peri-LMC} / \Gyr$  | -0.1?         |
| Number of pericentres            | 1-2           |

Table: Properties of selected orbits for Sculptor. The mean orbit represents the observational mean from @tbl:scl_obs_props. The \smallperi{} represents instead the $3\sigma$ smallest pericentre, which we use to provide an upper limit on tidal effects. {#tbl:scl_orbits  short="Sculptor Selected Orbits"}



### Tidal effects from the LMC

Unlike the previous mdoel, this model only has one pericentre with the Milky Way. @fig:scl_lmc_sim_images shows snapshots of Sculptor over the past 2 Gyr while marking the position of the LMC. With only one pericentre, and a larger one than the `smallperi` orbit, Sculptor's dark matter is substantially less disrupted. And while, based on the tidal tensor values, the LMC induces a greater instantaneous tidal effect than the Milky Way, Sculptor's dark matter component does not show strong effects due to the LMC. 

Finally, @fig:scl_lmc_i_f shows the projected on-sky final stellar distribution and the initial and final stellar density profile for this model with exponential stars. The stellar component does not change at all. While the break radius set by the time since the last LMC pericentre would agree with the location of the break in the observed density profile, no stellar effect would be observable. 



The LMC flyby encounter is an approximately impulsive encounter, in contrast with more adiabadic mass loss due to the Milky Way. Impulsive encounters tend to inject energy into the stellar and dark matter distribution, and can initially cause the galaxy to contract. In addition, the tidal force is required to be far larger than for slower, adiabadic mass loss, because the galaxy experiences this tidal field for far less time. Even in models where Sculptor passes through the the LMC, tidal tails do not form immediately after this encounter. 



![Sculptor Simulation Snapshots with LMC](figures/scl_lmc_sim_images.png){#fig:scl_lmc_sim_images}





![Sculptor initial and final density with LMC](figures/scl_lmc_i_f.pdf){#fig:scl_lmc_i_f}

Figure: The tidal effects on the stellar surface density due to the LMC today.

### Sculptor's long-term orbital history is unknown



![Sculptor Orbits with LMC](figures/scl_lmc_orbits_mass_effect.png){#fig:scl_lmc_orbits_mass}

Figure: The long-term orbital history of Sculptor is uncertain. In light, transparant lines, the orbits of Sculptor in three different LMC/MW mass models from @vasiliev2024 are shown for random samples of the observed properties. The LMC orbits are in solid, thick lines of the corresponding colour. The L2M11 has a lighter LMC mass, and the L3M10 model has a lighter MW mass than our fiducial L3M11 LMC model. 

A final note is that this model only considers the recent tidal effects .As illustrated in the range of possible orbits in @fig:scl_lmc_orbits_effect, there is a chance that Sculptor experienced an extremely small pericentre with the Milky Way. This pericentre has the potential to substantially rearrange the stellar component, drive large tidal mass loss, and create a Plummer-like stellar density profile. However, this encounter is highly dependent on the choice of the MW-LMC potential model, and may not occur at all. Long time integration in dynamic potentials amplifies uncertainties in the inputs and may not be reliable. 

To illustrate the possible effects of this previous pericentre, we create a model with a goal of matching the $3\sigma$ smallest pericentre with the MW. Because of strong tidal shocking, the final conditions depend strongly on the initial conditions. In addition, actions are not conserved in this evolving potential, so we use the @vasiliev2024 method of iteratively updating the initial conditions by empirically estimating the Jacobian. Our final position reaches adequate, but not perfect, agreement with the intended final position. 

As a result of this proceedure, while the point particle reaches a pericentre of $\sim 4\,\kpc$, the adjusted N-body orbit instead only reaches a pericentre of $\sim 12\,\kpc$. 

# Ursa Minor



| parameter                     | Scl               | UMi               |                                |
| ----------------------------- | ----------------- | ----------------- | ------------------------------ |
| $r_{\rm break}$ observed      | $25 \pm 5$ arcmin | $30 \pm 5$ arcmin | @sestito+2024a; @sestito+2024b |
| $t_{\rm last\ peri}$ required | $110\pm30$ Myr    | $120\pm30$ Myr    | @eq:r_break                    |
| $r_{\rm peri}$ required       | 16kpc             | 14kpc             | @eq:r_jacobi                   |
|                               |                   |                   |                                |
| $r_J$ from orbits             |                   |                   |                                |
| $r_{\rm break}$ from orbits   |                   |                   |                                |

Table: Given the observed break radii, the required time of last pericentre and pericentre to produce the observed tidal effect.





![Ursa Minor simulation snapshots](figures/umi_sim_images.png)

Figure: Ursa Minor simulation images. 



Figure: Velocity dispersion evolution of Ursa Minor



![Ursa Minor simulated density profiles](figures/umi_smallperi_i_f.pdf)

Figure: The tidal effects on the stellar surface density.

## Effects of the LMC



![Ursa Minor orbits with LMC](figures/umi_lmc_xyzr_orbits.png)

Figure: Orbits of Ursa Minor with (orange) and without (green) an LMC. The final positions of Ursa Minor and the LMC are plotted as scatter points and the solid blue line represents the LMC trajectory. Note that the LMC only increases Ursa Minor's peri and apo-centres, weakening any tidal effect. Interestingly, there is a change that Ursa Minor may have once been bound to the LMC (diverging orange lines at top left of middle panel.)

