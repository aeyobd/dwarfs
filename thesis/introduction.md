# Background & Past Work

- What is dark matter? Why do we look at dwarfs?
- Forms of dark matter, lambda-CDM, and dwarf galaxies
- How does gravity affect dwarfs, theory of tidal perturbations
  - @EN2021, @PNM2008, etc.
- Instances of dwarfs undergoing weird processess
- Alternative processes and uncertainties in the evolution of dwarfs

The classical dwarfs are some of the earliest discovered systems, begining with @shapley1938



- @fattahi+2013, @fattahi+2018
- @sanchez-salcedo+hernandez2007: mond in dsph
- @mayer+2001 theory of tidal stripping
- @IH1995 structural parameters
- @mateo1998, @simon2019 ARAA on ultrafaints.
- @walker+penarrubia2011: measuring cores and cusps
- @battaglia+nipoti2022: dynamics and dwarf galaxies review 
- @dsouza+bell2022, @santistevan+2024. challenges to backwards time integration



![Picture of Sculptor](/Users/daniel/thesis/figures/scl_des_dr2.png){#fig:scl_image width=390pt height=390pt}

Figure: Image of the Sculptor dwarf spheroidal galaxy from Dark Energy Survey Data Release 2 [@abbott+2021; image created with HiPS2FITS]. Sculptor appears as a fairly prominent, extended over density of predominantly faint, red stars. 0.5 degree field of view centred on Sculptor. **TODO: combine with UMi, maybe add LMC and Fornax?**



![Picture of Ursa Minor](figures/umi_DSS2_0.75deg.png){#fig:umi_image width=390pt height=390pt}

Figure: Image of Ursa Minor dwarf galaxy from the Digitized Sky Survey 2 (0.75 deg field of view tangent plane). UMi appears as a diagonal/elliptical haze of faint, reddish stars from the top left to the bottom right.  Even as classical dwarf, Ursa minor is fairly diffuse and does not stand obviously out from the background.





![Dwarf galaxies sky position](figures/mw_satellites_onsky_scl_umi_for_1.png){#fig:mw_satellite_system width=390pt weight=219.375pt}

Figure: The location of MW dwarf galaxies on the sky. Dwarf galaxies are taken as confirmed dwarfs with MW or LMC hosts from the @pace2024 catalogue (version 1.0.3. We label all classical dwarfs and the Magellanic clouds, changing symbols based on the relevance to this work (classical dwarfs, Magellanic clouds REF, and the systems in question or comparison). The background image is from ESA/Gaia/DPAC (https://www.esa.int/ESA_Multimedia/Images/2018/04/Gaia_s_sky_in_colour2). 



## Motivation 

- Sculptor and Ursa Minor are possibly unusual.



![Idealized simulations match Scl and UMi](figures/scl_umi_vs_penarrubia.png){#fig:toy_profiles}

Figure: Sculptor and UMi's profiles are well-matched to @PNM2008.  





# Cosmological context

![Cosmological Power Spectrum](figures/power_spectrum.png){#fig:cosmological_power_spectrum width=100% }

Figure: The matter power spectrum under different assumptions for dark matter. Dwarf galaxies occupy the middle and low end of the blue region (10^10 - 10^8 solar masses), enabling a unique window into properties of dark matter on small scales. The smaller scales we can understand dark matter, the better we are able to test different models of dark matter. figure 1 from @bechtol+2022. 



Figure: Density profiles of comological simulated halos, matching approximantly the NFW formula. 

Origin of NFW profiles



# Theoretical Background

- Cosmology foundations
  - power spectrum plot

- NFW plot (density or energy), maybe borrow from paper?
  - Explain origin of NFW
- Dwarf galaxy formation, halos

N-body DM simulations

- Collisionless Boltzmann equation and meaning of such simulations
- Assumptions & context & past work
- Evolution under tidal field



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

## Introduction to Dark Matter simulations

In this section, we will cover

- How are dark matter simulations conducted
- Interpretations and uncertainties 
- Methods

### Background

### The N-body method

- Nbody as a solution to the collisionless boltzman equation
- limitations and interpretation (MC sampling of phase space)
- Tree approximations and the BH method





## Tidal effects on dark matter halos

- EN2020 tidal tracks
- peñarrubia2009/etc.
- assumtions and limitations



![Break radius validation](/Users/daniel/thesis/figures/idealized_break_radius.png){#fig:idealized_break_radius}

Figure: The break radius of the simulations is set by the time since pericentre. 



From this argument, we note that the following properties must be approximately true for tides to occur:

- Close enough pericentre. The other break radius $r_J$ implies that if the host density is 3x the satellite, stars will be lost
- Corresponding time since last pericentre: If the time since last pericentre is not $\sim$ consistent with an observed break in the density profile, then tides 
- Halo evolution. As found in @EN2021, galaxies evolve along well defined tidal tracks (assuming spherical, isotropic, NFW halo, which may not be true, see ...). These tracks tend to "puff up" the stellar component while also removing dark matter mass, leaving a smaller, compacter DM halo with a more extended stellar component.
  - This information is mostly related to the statistical initial distribution of satellites from cosmology [ludlow+2016; @fattahi+2018]

## Potentials

### NFW

From @nfw1996, eqns. 3 & 4
$$
\frac{\rho}{\rho_c} = \frac{\delta_c}{(r/r_s)(1+r/r_s)^2}
$$
where 
$$
\delta_c = \frac{200}{3}\frac{c^3}{[\ln(1+c)-c/(1+c)]}
$$
$c$ is concentration parameter, $\rho_c$ is the critical density of the universe, and $r_s$ is the characteristic scale length of the halo.

The NFW halo is sometimes described by $M_{200}$. $r_{200}$ is the radius at which the mean density of the halo interior to $r_{200}$ is 200 times the critical density of the universe, and $M_{200}$ is the mass contained inside $r_{200}$. As equations:
$$
\rho_{200} = 200\rho_{c} = \frac{M_{200}}{(4\pi/3) r_{200}^3}
$$

$$
r_{200} = \sqrt[3]{\frac{1}{200}\frac{3M_{200}}{4\pi \rho_{\rm c
}}}
$$


$$
M_{200} = \frac{4\pi}{3} r_{200}^3\ \rho_{200}
$$

$M_{200}$ is also sometimes called the virial mass of the halo. $r_{200}$ is directly related to $r_s$ by
$$
r_{200} = c\,r_s
$$

The circular velocity in terms of $v_{200} = \sqrt{G M_{200} / R_{200}}$ is
$$
\left(v_{\rm circ}/v_{200}\right)^2 = \frac{A(x)/x}{A(c)/c},
$$

or in terms of $M_s$ and $r_s$, 
$$
v_{\rm circ}^2 = \frac{G M(r)}{r} = \frac{G M_s A(r/r_s)}{r}.
$$
Another parameterization of the NFW profile is in terms of the maximum circular velocity $v_{\rm circ}^{\rm max}$ and the radius at which it is reached $r_{\rm circ}^{\rm max}$. Given the scale radius, 
$$
r_{\rm circ}^{\rm max} = \alpha\ r_s
$$

where $\alpha\approx2.16258$, and $v_{\rm circ}^{\rm max}$ can be found from either of the equations above.

