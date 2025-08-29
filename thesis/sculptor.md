In this chapter, we present our results from our simulations of Sculptor and Ursa Minor. In each case, we first discuss the range of reasonable orbits and initial dark matter halos for each galaxy. Then, we describe how tides affect the evolution of each galaxy. We furthermore consider how the effects of the LMC changes the results, including using additional N-body simulations in evolving potentials. Finally, we review the predictions and characteristics of the present-day stellar populations and properties from our models.

# Sculptor

## Milky Way tides



## Tidal effects

The Milky Way tides indeed affect Sculptor's dark matter halo. @fig:scl_sim_images shows projected dark matter distribution of our N-body simulation for Sculptor in the $y$-$z$ plane for 5 previous pericentres. This simulation uses the `compact` initial halo on the `smallperi` orbit. The dark matter forms large streams which orbit the Milky Way several times. However, note that the central dark matter cusp is similar in appearance across all simulation snapshots. 

More quantitatively, the halo loses about 99% of its initial dark matter mass. This corresponds in a reduction from the initial $v_{\rm max} = 31\,\kms$ to $v_{\rm max} = 22\,\kms$.  However, the stellar component of Sculptor is relatively unaffected. 

While mass loss leads to adiabatic expansion of the stellar component, no features of tidal disruption or disequilibrium are currently observable. @fig:scl_smallperi_i_f shows the final projected density of stars on the sky, and the initial and final radially averaged density profiles. No non-spherical density features are apparent, even at 5 decades fainter in surface density than the central density. In addition, the initial and final density profiles look identical up to some scale, with a tidally-induced excess of stars only appearing $\sim 3$ orders of magnitude fainter than the faintest our density profile measures. 

The analytically motivated break and tidal radii corroborate a weak tidal effect on stars. In both cases, these radii work out to be $\gtrsim 100$ arcminutes, outside where we can measure the density profile. In particular, because we have chosen the most extreme observationally permissible orbit, it is unlikely, for similar $\Lambda$CDM motivated initial conditions, that these radii would ever approach where the break in Sculptor's density profile appears. 

![Sculptor simulation snapshots](figures/scl_sim_images.png){#fig:scl_sim_images}

Figure: Images of the dark matter evolution over a selection of past apocentres and the present day position. Limits range from -150 to 150 kpc in the $y$-$z$ (approximately orbital) plane and the colourscale is logarithmic spanning 5 orders of magnitude between the maximum and minimum values. In this image, stars occupy only ever a few pixels so are not plotted. 



![Sculptor Tidal Tracks](figures/scl_tidal_track.pdf){#fig:scl_tidal_track}

Figure: The tidal tracks for the smallperi orbit. Todo: add velocity dispersion plot to RHS



![Sculptor velocity dispersion evolution](figures/scl_sigma_v_time.pdf)

Figure: Evolution of stellar velocity dispersion within 1 kpc for different Scl models. In all cases, the evolution is mild. Note that binarity may reduce the inflate the observed velocity dispersion by  ~ 1 km/s, so the conservative lower limit is around 8 km/s. **TODO**, add bound mass with time. Maybe combine with tidal tracks and radius / time orbit.





![Sculptor initial and final density profiles](figures/scl_smallperi_i_f.pdf){#fig:scl_smallperi_i_f}

Figure: Effects on exponential initial stars. TODO: plot 2D sky proj. stars





![Sculptor Plummer initial and final density profiles](figures/scl_plummer_i_f.pdf){#fig:scl_smallperi_plummer_i_f}

Figure: effects on Plummer initial stars.

## Effects of the LMC {#sec:scl_lmc}

The Milky Way isn't the only galaxy in town. Becoming more apparent, the infall of the LMC substantially affects the Milky Way system. The mass of the LMC may be as high as 1/5 of the Milky Way mass, so including an LMC in the Milky Way potential may change conclusions about properties and orbits of the Milky Way satellite system. 

### Orbital effects

As discussed in @battaglia+2022, Sculptor's orbit is strongly influenced by the presence of an LMC. 

To explore the effects of the LMC on Sculptor, we use the publically available L3M11 model from @vasiliev2024. This model stars with a lighter Milky Way than our fiducial @EP2020, @mcmillan2011 like model. The L3M11 potential is an evolving multipole approximation of an N-body simulation with both a live MW and LMC dark matter halo. The MW intital conditions were an NFW with $r_s=16.5\,$kpc and mass $M_{\rm 100}= 11\times10^{11}\Mo$, and the LMC was a NFW halo with $r_s=11.7$ and $M_{100} = 2.76 \times 10^{11} \Mo$. This model had a previous LMC pericentre at about 6 Gyr ago.

@fig:scl_lmc_orbits_effect displays the effect of including an LMC in the potential.   The green samples are in the initial MW-only potential in the `L3M11` model, and the orange samples are integrated in the evolving MW and LMC `L3M11` model. The past 1 Gyr is similar in both cases, but the orbits diverge significantly afterwards. The recent passage of Sculptor with the LMC around 0.1 Gyr ago allows for Sculptor to begin as far as 300 kpc from the Milky Way centre. The evolving potential also adds significant long term variability in the possible orbits of Sculptor. Sculptor, however, is orbiting in the opposite direction of the LMC so is likely not associated with the LMC system. 

Given the large uncertainties of the LMC model, we conservatively double all of the observational parameters of Sculptor. This has a similar effect to including LMC mass and orbital uncertainties but is considerably simpler. From this larger range of orbits, we once again select the orbit with the median final observables of all orbits with pericentres less than the 0.0027th quantile. This orbit, the `smallperilmc` orbit is plotted in black and is only integrated up to 2 Gyr ago, to isolate recent tidal effects. 



![Sculptor Orbits with LMC](figures/scl_lmc_xyzr_orbits.pdf){#fig:scl_lmc_orbits_effect}

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



### Tidal effects

Unlike the previous mdoel, this model only has one pericentre with the Milky Way. @fig:scl_lmc_sim_images shows snapshots of Sculptor over the past 2 Gyr while marking the position of the LMC. With only one pericentre, and a larger one than the `smallperi` orbit, Sculptor's dark matter is substantially less disrupted. And while, based on the tidal tensor values, the LMC induces a greater instantaneous tidal effect than the Milky Way, Sculptor's dark matter component does not show strong effects due to the LMC. 

Finally, @fig:scl_lmc_i_f shows the projected on-sky final stellar distribution and the initial and final stellar density profile for this model with exponential stars. The stellar component does not change at all. While the break radius set by the time since the last LMC pericentre would agree with the location of the break in the observed density profile, no stellar effect would be observable. 



The LMC flyby encounter is an approximately impulsive encounter, in contrast with more adiabadic mass loss due to the Milky Way. Impulsive encounters tend to inject energy into the stellar and dark matter distribution, and can initially cause the galaxy to contract. In addition, the tidal force is required to be far larger than for slower, adiabadic mass loss, because the galaxy experiences this tidal field for far less time. Even in models where Sculptor passes through the the LMC, tidal tails do not form immediately after this encounter. 

A final note is that this model only considers the recent tidal effects .As illustrated in the range of possible orbits in @fig:scl_lmc_orbits_effect, there is a chance that Sculptor experienced an extremely small pericentre with the Milky Way. This pericentre has the potential to substantially rearrange the stellar component, drive large tidal mass loss, and create a Plummer-like stellar density profile. However, this encounter is highly dependent on the choice of the MW-LMC potential model, and may not occur at all. Long time integration in dynamic potentials amplifies uncertainties in the inputs and may not be reliable. **TODO: I would love to get this model to work...**

![Sculptor Simulation Snapshots with LMC](figures/scl_lmc_sim_images.pdf){#fig:scl_lmc_sim_images}





![Sculptor initial and final density with LMC](figures/scl_lmc_i_f.pdf){#fig:scl_lmc_i_f}

Figure: The tidal effects on the stellar surface density due to the LMC today.

# Ursa Minor







![Ursa Minor simulation snapshots](figures/umi_sim_images.png)

Figure: Ursa Minor simulation images. 



Figure: Velocity dispersion evolution of Ursa Minor



![Ursa Minor simulated density profiles](figures/umi_smallperi_i_f.pdf)

Figure: The tidal effects on the stellar surface density.

## Effects of the LMC



![Ursa Minor orbits with LMC](figures/umi_lmc_xyzr_orbits.pdf)

Figure: Orbits of Ursa Minor with (orange) and without (green) an LMC. The final positions of Ursa Minor and the LMC are plotted as scatter points and the solid blue line represents the LMC trajectory. Note that the LMC only increases Ursa Minor's peri and apo-centres, weakening any tidal effect. Interestingly, there is a change that Ursa Minor may have once been bound to the LMC (diverging orange lines at top left of middle panel.)

# Summary