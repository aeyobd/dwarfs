# Gaia Membership Selection

![Sculptor selection criteria](/Users/daniel/thesis/figures/scl_selection.png){#fig:sculptor_selection}

Figure: The selection criteria for Scl members. Members are red and all field stars (satisfying quality criteria) are in light grey. **Top left:** tangent plane. **Top Right:** Colour magnitude diagram.

We use J+24 data. J+24 select members using a multi-component Baysian algorithm:

- Remove stars with poor astrometry or photometry, no colour excess (@lindegren+2018 equation C.2), 3$\sigma$ consistency of measured parallax with dwarf distance (near zero with @lindegren+2018 zero-point correction), and absolute RA and Dec proper motions less than 10$\,{\rm mas\ yr^{-1}}$.

- Spatial likelihood based on a double exponential $\Sigma_\star \propto e^{-r/r_s} + B\,e^{-r/r_{\rm outer}}$ where the inner scale radius is fixed. 

- Stars are assigned a likelihood based on the location on the CMD (using Padova isochrones including an intrinsic 0.1 CMD width in colour convolved with colour and distance modulus)

- Background KDE density maps for the CMD and PM are constructed using the other quality-selected stars outside of $5R_h$, where the satellite density would be orders of magnitude less than the background (even in the presence of extended tidal features).

- Likelihoods normalized to unity to represnt a PDF

  

  Membership probabilities are then given by

$$
P_{\rm sat} = \frac{f_{\rm sat}{\cal L}_{\rm sat}}{f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}} = \frac{1}{1 + \frac{(1-f_{\rm sat}){\cal L}_{\rm bg}}{f_{\rm sat}{\cal L}_{\rm sat} }}
$$

where $f_{\rm sat}$ is the fraction of stars belonging to the system inside the given field, ${\cal L}_{\rm sat}$ is the likelihood of a star belonging to the satellite, and ${\cal L}_{\rm bg}$ is the likelihood of the star belonging to the background. Each likelihood is calculated as a product of the CMD, PM, and spatial likelihoods:
$$
{\cal L} = {\cal L}_{\rm space}\ {\cal L}_{\rm PM}\ {\cal L}_{\rm CMD}
$$

The above formula suggests that a cut in $P_{\rm sat}$ is equivalent to the cut in likelihoods
$$
\frac{{\cal L}_{\rm sat}}{{\cal L}_{\rm bg}} > \frac{(1-f_{\rm sat})/f_{\rm sat}}{1/P_{\rm sat}- 1}
$$

Note that if we remove the spatial component of the likelihood, then $f_{\rm sat}$ represents a global normalization.



Not shown here, we explore simple cuts of the stars, using absolute cuts in parallax, proper motions, and the CMD. The results are similar to the nospace model.



### Searches for tidal tails

![Tidal tails](/Users/daniel/thesis/figures/scl_tidal_tails.png){#fig:sculptor_tidal_tails}

Figure: The distribution of member stars (orange), PM & CMD only selected stars (blue) and all stars (passing quality cuts, black).





- There are no apparent overdensities in the PM & CMD only selected stars to suggest the presence of a tidal tail

- This means that at least at the level of where the background density dominates, we can exclude models which produce tidal tails brighter than a density of $\Sigma_\star \approx 10^{-2}\,\text{Gaia-stars\ arcmin}^{-2} \approx 10^{-6} \, {\rm M_\odot\ kpc^{-2}}$ (TODO assuming a distance  of ... and stellar mass of ...). 

  


## Density Profile Reliability and Uncertainties

- How well do we know the density profiles? 
- What uncertainties affect derived density profiles? 
- Can we determine if Gaia, structural, or algorithmic systematics introduce important errors in derived density profiles?



J+24's algorithm takes spatial position into account, assuming either a one or two component exponential density profile. When deriving a density profile, this assumption may influence the derived density profile, especially when the galaxy density is fainter than the background of similar appearing stars. To remedy this and estimate where the background begins to take over, we also explore a cut based on the likelihood ratio of only the CMD and PM components. This is in essence assuming that the spatial position of a star contains no information on it's membership probability (a uniform distribution like the background)

To incorporate the structural uncertainties and robustly model the sampling uncertainty, we construct the following bootstrap model

- Centre is varied by a centring error, estimated from the standard normal error of the positions plus the systematic shift of the mean

- Position angle and ellipticity are sampled from a normal distribution given the reported uncertainties 

- $f_{\rm sat}$ is sampled from ...

  

  **TODO**: Look into normalization of Likelihoods and check how $f_{\rm sat}$ matters. If $f_{\rm sat}$ is related to the normalizations of fg / bg densities, and other likelihoods are area-normalized to 1, then this makes life much easier. 

  Test if psat weighted density profiles are similar

  Save MC density profile outputs

  

  

  



![Density profiles](/Users/daniel/thesis/figures/scl_density_methods.png){#fig:sculptor_observed_profiles}

Figure: 



# Comparison of the Classical dwarfs

- Using J+24 data, we validate
  - Check that PSAT, magnitude, no-space do not affect density profile shape too significantly
- Our "high quality" members all have > 50 member stars and do not depend too highly on the spatial component, mostly corresponding to the classical dwarfs
- We fit Sérsic profiles to each galaxy
  - The Sérsic index, $n$, is a measure of the deviation from an exponential. Exponentials have $n=1$, whereas more extended dwarf galaxies will have higher $n$
- To better estimate the uncertainties due to unknown galaxy properties and flexibility in the likelihood cut, we can 

![Sculptor and UMi versus classical dwarfs](/Users/daniel/thesis/figures/classical_dwarf_profiles.png)

Figure: Density profiles for each dwarf galaxy.



**TODO**: Use updated density profiles (nospace) with uncertainties included to MCMC fit Sérsic profiles to every dwarf galaxy. 

- Fornax
- Leo I
- Sculptor
- Leo II
- Carina 
- Sextans I
- Ursa Minor
- Draco
- Antlia II?

## Radial Velocity Measurements





# Simulation-based motivations

To motivate why a tidal interaction may give rise to the observed density profiles, we create a toy simulation following @PNM2008. 

- NFW initial conditions (sculptor like, vcm, rcm)

- Evolved in x-y plane using @EP2020 potential for ~ 5Gyr with pericentre of 15 kpc and apocentre of 100 kpc. 

- Exponential initial stellar profile.

  

As a dark matter halo is perturbed on a pericentric passage with the milky way,

- Tidal stress heats halo slightly
- Mass loss, particularly of loosely bound particles

The stellar component tracers will similarly follow the behaviour of the dark matter. 

An emperical estimate of where the simulation's stars are becoming unbound is, as stated in @PNM2008, the break radius
$$
R_b = C\,\sigma_{v}\,\Delta t
$$
where $\sigma_v$ is the present line of sight velocity dispersion , $\delta t$ is the time since pericentre, and $C \approx 0.55$ is a fit. The idea motivating this equation is stars in the inner regions will have dynamically equilibriated to the new potential (phase mixed), however the outer regions are no longer in steady state, so we have to wait until the crossing time reaches them as well.



As illustrated in @fig:toy_profiles, the density profile initially stars off exponential. At increasing times since the first pericentric passage, the break radius, appearing as an apparent separation between the slopes of the inner and outer profile, increases. 

![Idealized simulations match Scl and UMi](/Users/daniel/thesis/figures/scl_umi_vs_idealized.png){#fig:toy_profiles}

Figure: Sculptor and UMi's profiles are well-matched to an idealized simulation



![Break radius validation](/Users/daniel/thesis/figures/idealized_break_radius.png){#fig:idealized_break_radius}

Figure: The break radius of the simulations is set by 



From this argument, we note that the following properties must be approximately true for tides to occur:

- Close enough pericentre. The other break radius $r_J$ implies that if the host density is 3x the satellite, stars will be lost
- Corresponding time since last pericentre: If the time since last pericentre is not $\sim$ consistent with an observed break in the density profile, then tides 
- Halo evolution. As found in @EN2021, galaxies evolve along well defined tidal tracks (assuming spherical, isotropic, NFW halo, which may not be true, see ...). These tracks tend to "puff up" the stellar component while also removing dark matter mass, leaving a smaller, compacter DM halo with a more extended stellar component.
  - This information is mostly related to the statistical initial distribution of satellites from cosmology [ludlow+2016; @fattahi+2018]

# Summary

- Of the classical dwarfs, UMi & Scl stand out statistically, with high $n$ given their luminosity
- Including fainter dwarf galaxies, Boo 3 and Boo 1 appear to also have extended density profiles
  - Deeper data would be required to robustly measure this



