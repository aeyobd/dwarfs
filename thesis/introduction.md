Dwarf galaxies host, in many ways, the most extreme galactic environments in the universe. These little galaxies are typically defined to fainter than the LMC or SMC [$M_V \gtrsim -18$ , e.g.,@hodge1971; @mcconnachie2012]. As the smallest class of galaxies, dwarfs are the most numerous galaxies. Dwarf galaxies are highly *dark-matter dominated*, with mass to light ratios sometimes exceeding 100 to 1000 $M_\odot/ L_\odot$.  Because of their small total masses, many dwarf galaxies are *quenched*, with little to no recent star formation. Their stellar populations are *relics* from the early universe, consisting of many of the oldest and most metal poor stars. Understanding the properties of dwarf galaxies thus has implications across astronomy, from cosmological structure formation on the smallest scale to extremely metal-poor stellar populations to chemical evolution to galactic dynamics. Of the classical dwarfs, the Sculptor and Ursa Minor dwarf spheroidal galaxies stand out with a more extended density profile relative to an exponential, hinting at tidal effects. As a case study in tidal effects and our understanding of dwarf galaxy formation and evolution, we aim to understand the dynamical history of these galaxies. 

In this section, we first describe cosmological structure formation and the role of dwarf galaxies. Then, we explore the observational history of dwarf galaxies, the *Gaia* telescope, open questions, and why Sculptor and Ursa Minor stand out. We  discuss idealized simulations, tidal effects, break radii, and recent developments in tidal simulations. Finally, we provide a roadmap to the thesis as a whole at the end of this section.

# Cosmological context

We only understand a tiny fraction of the universe's composition. The leading theory of cosmology, $\Lambda$CDM (cold dark matter), posits that the universe is composed of about 68% dark energy ($\Lambda$), 27% dark matter (DM), and 5% regular baryons[^baryons]  [@planckcollaboration+2020]. While the composition of dark matter and dark energy remains elusive, we know their general properties. Dark energy causes the acceleration of the expansion of the universe on large scales. We do not discuss dark energy here—it does not substantially affect the local group today. Dark matter, instead, makes up the vast majority of mass in galaxies. Typically, galaxies have baryonic to dark matter ratios of between 1:5 to beyond 1:1000 for faint dwarf galaxies. In $\Lambda$CMD, dark matter is assumed to only interact gravitationally, passing through matter without effect (transparent to light, or *dark*). Dark matter is also *cold*, i.e. typical velocities much smaller than the speed of light in the early universe. Implications of dark matter properties range from cosmological structural formation, galaxy structure, and galaxy interactions. 

[^baryons]: Astronomers like to change definitions of words. *Baryons* here means baryons+leptons, i.e. any standard model massive fermion. The photon energy density is negligible today

## Structure formation and dwarf galaxies

The very early universe was almost featureless. Our earliest observations of the universe stem from the cosmic microwave background (CMB)---revealing a uniform, isotropic, near-perfect blackbody emission. But tiny perturbations in the CMB, temperature fluctuations of 1 part in 10,000, became the seeds of large-scale cosmological structure. Governed by cosmological expansion, gravitational collapse, and baryonic physics, each perturbation grows into larger structures. Initially, baryonic matter was coupled to radiation and resisted collapse. Dark matter, only influenced by gravity instead, freely collapsed into the first structures. Each self-gravitating overdensity of dark matter is known as a *halo*. After recombination, where electrons combined with atomic nuclei to form atoms, baryons decoupled from radiation and fell into the dark matter halos. The densest pockets of baryons later formed the first stars and galaxies.

Dark matter halos, and their associated galaxies, rarely evolve in isolation. Instead, structure formation is *hierarchical*. Small dark matter halos collapse first and hierarchically merge into progressively larger halos [e.g.,@blumenthal+1984; @blumenthal+1986; @white+rees1978; @white+frenk1991; @somerville+dave2015]. Structure formation happens on a wide range of scales. The largest structures become the cosmic web—composed of voids, filaments, and clusters—and the smallest structures directly observable host dwarf galaxies. Hierarchical assembly is evident through the large scale structure of the universe, remnants of past mergers within the milky way, and tidal disruption of dwarf galaxies and their streams around nearby galaxies.

Small-scale structure formation is sensitive to deviations from $\Lambda$CDM cosmology [e.g.,@bechtol+2022]. Many alternative models, such as warm dark matter or self-interacting dark matter, may smooth out small-scale features and reduce the abundance of small halos or change their structure [e.g.,@lovell+2014]. Dwarf galaxies, occupying the smallest dark matter halos, are promising windows into small-scale cosmological features. Alternative dark matter properties can influence dwarf galaxy structure, formation, and tidal evolution. The nearby dwarf galaxies of the local group, with detailed observations, present a promising opportunity to understand the evolution of these objects and test our understanding of cosmology. To understand the evolution of these objects, we need to understand the general predictions from $\Lambda$CDM for the properties of dwarf galaxies. 

## Typical structure of dark matter halos

In $\Lambda$CDM cosmological simulations, dark matter halos are remarkably self-similar. In \citet{NFW1996, NFW1997}, hereafter NFW, the authors observe that the spherical radial density profiles $\rho(r)$ are well described by a two-parameter law: 
$$
\rho/\rho_s= \frac{1}{(r/r_s)(1+r/r_s)^2},
$$
where $r_s$ is a scale radius and  $\rho_s$ a scale density . This profile has shown remarkable success in describing $\Lambda$CMD halos across several orders of magnitude in mass. NFW profiles are *cuspy*, where the density continuously rises like $\rho \sim 1/r$ at small radii $r \ll r_s$. The steepness of the density profile increases at $r \sim r_s$ and at large radii, the density  falls off like $\rho \sim 1/r^3$. @fig:nfw_density shows an example of a NFW halo and a cored dark matter halo. 

The total mass of an NFW profile diverges, so halos are commonly characterized by an overdensity criterion. The virial mass, $M_{200}$, is defined as the mass within a region $r_{200}$ containing a mean enclosed density 200 times[^200] the critical density of the universe:
$$
M_{200} =(4\pi/3) \ r_{200}^3\ 200\rho_{\rm crit}, \qquad {\rm where} \quad \rho_{\rm crit}(z) = 3H(z) / 8\pi G
$$
A standard second parameter is the halo concentration, $c=r_{200} / r_s$ describing how the characteristic size of the halo compares to its virial radius.  In this case, the scale density is a function of $c$ alone, $\rho_s = (200/3)\,\rho_{\rm crit} c^3 / [\log(1+c) - c/(1+c)]$ [@NFW1996]. However, we characterize halos by their circular velocity profiles. The circular velocity, $v_{\rm circ}(r) = \sqrt{G M(r) / r}$, reaches a maximum of $v_{\rm max}$ at radius  $r_{\rm max} = \alpha\,r_s$ where $\alpha\approx2.16258$. $v_{\rm max}$ is related to both the total halo mass and the observed line of sight velocity disperion, and $r_{\rm max}$ relates to the scale radius. 

While an NFW profile has two free parameters ($M_{200}$ and $c$), these are not independent. Smaller dark matter halos often collapse earlier, when the universe was denser. As a result, small subhalos tend to be more concentrated [e.g.,@NFW1997]. The relationship between $M_{200}$ and c, or the mass-concentration relation, how more massive halos become less concentrated in a predictable way [e.g., @bullock+2001; @ludlow+2016]. While there is a general trend in concentration with mass, the concentration values still have substantial scatter. Other parameters such as the halo spin or shape or slight changes in the halo density profile add additional variability to the present-day distribution of halos [see e.g., @navarro+2010; @dicintio+2013; @dutton+maccio2014, darkEXP ]. 



[^200]: 200 smells like an arbitrary value. From approximate analytic arguments about the mean density of a self-gravitating spherical collapse halo in an expanding universe, the expected virial overdensity is XXX, rounding up to 200 (REFS).

![Example density profiles](figures/example_density_profiles.pdf){#fig:nfw_density}

Figure: Example stellar and dark matter density profiles for a Sculptor-like galaxy. The dark matter is more extended and massive than the star across the entire galaxy.

## Connecting dark matter to stars

Both cosmology and observations find fundamental relationships between properties of galaxies, particularly for the total mass and luminosity. A key prediction of $\Lambda$CDM is the stellar mass halo mass relation, describing the amount of stellar mass forms in a (sub)halo of a given size. In particular, the SMHM relation grows especially steep in the dwarf galaxy regime---many dwarf galaxies are formed in halos of similar masses. @Fig:smhm shows the stellar mass versus $v_{\rm max}$ maximum circular velocity (proxy for halo mass) for apostle dwarfs from @fattahi+2018. While there is some scatter, the range of possible $v_{\rm max}$ is fairly narrow across $\sim 5$ decades in stellar mass. Thus, if we know the initial stellar mass of a dwarf galaxy, we also know, at some level, its dark matter mass.

Several challenges complicate a simple SMHM trend including environment, assembly history, and tidal effects. By being closer to a massive host, most dwarfs quench earlier, resulting in lower stellar masses [e.g., @christensen+2024]. Additionally, the time of formation  (relative to reionization) can influence the resulting stellar content [@kim+2024]. Finally, tides influence both the dark matter and stellar mass but in different amounts [e.g., @PNM2008]. Consequently, tides may reduce the halo mass more than the stellar mass, adding additional scatter to the SMHM trend, particularly for Milky Way satellites [e.g.,@fattahi+2018]. Understanding the effects of tides on local group dwarf galaxies is essential to determining where and how these galaxies formed in a cosmological context. 

![Stellar-mass halo-mass relation](/Users/daniel/Library/Application Support/typora-user-images/image-20250715100059332.png){#fig:smhm width=100%}

Figure: The stellar mass halo mass relation in terms of $v_{\rm max}$, halo maximum circular velocity related to halo mass, and $M_{\rm str}$ stellar mass. SRight: the mass-concentration relation for NFW halos, but parameterized in terms of $r_{\rm max}$ and $v_{\rm max}$. Together, these plots represent the cosmologically expected properties of underlying dark matter halos in any dwarf galaxy. The velocity dispersion directly constrains the mass  contained within a half-light radius, so the underlying halo in $\Lambda$CDM is well-constrained. Adapted from figure 1 of @fattahi+2018, stellar mass and $v_{\rm max}$ of APOSTLE satellites at peak $v_{\rm max}$.  

# Observational context

Dwarf galaxies have long questioned theories of cosmology and galaxy formation. The discovery of Fornax and Sculptor in 1938 [@shapley1938][^lmc_discovery], with no known analogues, already presented a conundrum. Shapley presented these dwarfs as a new type of *stellar system* resembling the Magellanic clouds and globular clusters but did not attempt to speculate on the exact nature. While generally understood to be galaxies based on the inferred luminosities and sizes, the nature and formation remained unclear [e.g.,@hodge1971; @gallagher+wyse1994]. 

The earliest spectroscopic work hinted that dwarf galaxies may contain substantial dark matter. From early determinations of the velocity dispersion for Sculptor and Ursa Minor dSph  [e.g., @aaronson1983, @aaronson+olszewski1987],  inferred mass-to-light ratios were at least 10 times the expectation for globular clusters. While only determined by a few stars initially, these values have solidified with time with larger and more precise samples [e.g., @hargreaves+1994]. Additionally, as we expect the dark matter component to be much more extended than the stars, the total mass-to-light ratios are much higher, reaching about 1000 $M/L$ for the galaxies considered here, requiring high-density concentrations of dark matter. Subsequently, several theories attempting to understand the formation and observed properties of these objects were proposed. Examples include: dwarf galaxies are undergoing tidal dissolution resulting in extreme mass-to-light ratios [e.g., @kuhn+miller1989], presence of massive central black holes [e.g., @strobel+lake1994], the formation of dark matter free out-of-equilibrium "tidal dwarfs" from past mergers [e.g., @lynden-bell1982, @kroupa1997], or modified theories of gravity [@milgrom1995].  However, consistency with ordinary galaxy formation was also not out of the question [e.g., @dekel+silk1986]. Since then, we have known that dwarf galaxies are among the darkest objects in the universe, and understanding their properties is critical to understanding the universe. 

Dwarf galaxies span a large range of physical sizes, luminosities, and morphologies. Broadly, there are three classes of dwarf galaxies based on luminosity, as illustrated by @fig:galaxy_images. Local **bright dwarfs galaxies** with magnitudes $-18 \lesssim M_V \lesssim  -14$, often exhibit irregular morphologies and recent star formation.  @fig:galaxy_images shows the LMC as an example of an irregular (and slightly bright) dwarf galaxy displaying a bar. **Classical dwarf** galaxies occupy intermediate luminosities ( $-14 \lesssim M_V  \lesssim -7.7$). Typically, these systems are old, gas-poor, and spheroidal. All Milky Way satellites discovered before digital sky surveys are classicals.  The 12 classical dwarfs satellites of our galaxy are Sagittarius, Fornax, Leo I, Sculptor, Antlia II, Leo II, Carina, Draco, Ursa Minor, Canes Venatici I, Sextans I, and Crater II.[^dsph_suffix]   The **ultra-faint**s occupy the very faintest magnitudes ($-7.7 \lesssim M_V$). These galaxies have minuscule stellar masses, tend to be more compact, and are the darkest known galaxies. An example is Boötes V as shown in [@fig:galaxy_images]. Altogether, dwarf galaxies span about 7 orders of magnitude in stellar mass.



[^lmc_discovery]: technically, the LMC and SMC may be classified as dwarf galaxies, but these were likely always known to humans at southern latitudes. 
[^dsph_suffix]: While formally the dwarf galaxy names we discuss contain" dwarf spheroidal" (dSph), e.g. Sculptor dSph, we omit this suffix for brevity.





![Dwarf Galaxy Pictures](figures/galaxy_pictures.pdf){#fig:galaxy_images width=390pt height=390pt}

Figure: Images of the LMC (DSS2), Sculptor (DES 2), Ursa Minor (UNWISE with Gaia point sources over-plotted), and Bootes V (SDSS). While the LMC is very prominent in the sky, even the classical dwarf galaxies are not obvious except in terms of star counts (e.g.,@fig:scl_selection). All (non-LMC) dwarfs occupy the central third to sixth of the image. **TODO**: replace Boötes V with For





![Dwarf galaxies sky position](figures/mw_satellites_1.jpg){#fig:mw_satellite_system width=390pt}

Figure: The location of MW dwarf galaxies on the sky. We label the classical dwarf galaxies (green diamonds), fainter dwarfs (blue squares), globular clusters (orange circles), and ambiguous systems (pink open hexagons). Globular clusters are more centrally concentrated, but dwarf galaxies are preferentially found away from the MW disk. Sculptor and Ursa Minor are highlighted as two dwarfs we study later. The background image is from ESA/Gaia/DPAC (https://www.esa.int/ESA_Multimedia/Images/2018/04/Gaia_s_sky_in_colour2). Dwarf galaxies (confirmed), globular clusters, and ambiguous systems are from the @pace2024 catalogue (version 1.0.3). 

## The era of *Gaia* 

A star's or galaxy's position and velocity are fundamental quantities for understanding its orbit and origin. Unfortunately, determining distances and full 3D velocities is challenging. The most direct measurement of distance, parallax, requires precise tracking of a star's sky position across a year. And while line-of-sight (LOS) velocities are easily determined from spectroscopy, tangental velocities, derived from proper motions, are much more challenging. Typically, measuring proper motions requires accurate (much less than arcsecond) determinations of small changes in a star's position over baselines of several years to decades. The full 6D position and velocity information for stars was historically only known for a handful of stars.

*Gaia* has redefined astrometry, providing photometry, positions, proper motions, and parallaxes for over 1 billion stars [@gaialcollaboration+2021]. Launched in 2013, *Gaia* is a space-based, all-sky survey telescope situated at the Sun-Earth L2 Lagrange point [@gaiacollaboration+2016]. While *Gaia* completed its space-based mission in 2025, two forthcoming data releases remain. Determining absolute parallax measurement is facilitated by the observation that stars in different regions of the sky are affected by parallax motion with different phases. By imagining two regions separated by 106.5 degrees on the same focal plane, *Gaia* measures tiny changes in relative positions of stars across small and large angles. Combining measurements from multiple epochs across several years, an absolute all-sky reference frame is derived from which parallax and proper motions are derived [@gaiacollaboration+2016]. In addition to astrometry, *Gaia* measures photometry in the wide *G* band (330-1050nm) and colours from the blue photometer (BP, 330-680 nm) and red photometer (RP, 640-1050 nm). *Gaia* additionally provides low resolution BP-RP spectra and radial velocity measurements of bright stars (*Gaia* radial velocity magnitudes <16) [@gaiacollaboration+2016]. For this work, *Gaia*'s most relevant measurement are $G$ magnitude, $G_{\rm BP} - G_{\rm RP}$ colour, $(\alpha, \delta)$ position, and $(\mu_{\alpha*}, \mu_\delta)$ proper motions.[^pmra_cosdec]

*Gaia* has revolutionized many fields of astronomy including the Milky Way and the Local Group. For example, 6D kinematic classification of Milky Way stars has lead to the discovery of past mergers or Milky Way building blocks like *Gaia*-sausage Enceladus [e.g., @helmi+2018], out-of-equilibrium structures like the *Gaia* snail [e.g., @antoja+2018l, and effects of the Milky Way's bulge and bar [@hunt+vasiliev2025]. Moving to the Milky Way halo, numerous streams  [@bonaca+price-whelan2025]. *Gaia* has revealed the hierarchical formation and complex, evolving structure of our own galaxy.

For Milky Way satellites, *Gaia* has enabled well-constrained orbital analysis and facilitated precise stellar membership determination. Before, proper motions of dwarf galaxies were spase or undetermined and measured by the Hubble Space Telescope or other methods for only a few systems. *Gaia* allowed for the first systematic determinations of Milky Way satellite proper motions, and of unprecedented precision [@MV2020a; @pace+li2019]. While the proper motion uncertainty on a typical dwarf member star is often with large uncertainties, by combining the proper motions for 100s to 1000s of stars in *Gaia*, incredibly precise proper motion measurements can be determined, sometimes only limited by *Gaia*'s systematic error floor  [e.g.,@MV2020a]. Proper motions have ushered in a new dynamical era for MW satellites studies, where we often have precise determinations of the dwarf galaxy's recent orbit under a given MW potential. In addition, *Gaia* helps seperate out contaiminating MW foreground stars. By measuring parallaxes or proper motions significantly different from the dwarf galaxy, many stars can be immediately removed as non-members [e.g., @battaglia+2022, @jensen2024].



[^pmra_cosdec]: The proper motions $\mu_\alpha$ and $\mu_\delta$ are the apparent rates of change in right ascension $\alpha$ and declination $\delta$, typically in units of mili-arcsecond (mas) per year. $\mu_{\alpha*}$ allows where $\mu_{\alpha*} = \mu_\alpha \cos \delta$ corrects for projection effects in $\alpha$.



## Dwarf galaxies today

Today, we know the Milky Way system is teeming with satellites. [@fig:mw_satellite_system] shows the MW satellite system, including dwarf galaxies, ambiguous systems, and globular clusters. Advances in telescopes and observations has accelerated progress on the nature of dwarf galaxy and introduced new questions. Deep digital photometric surveys have more than quadrupled the number of known Milky Way dwarf galaxies [@simon2019]. Upcoming and ongoing surveys, like the Vera Rubin Observatory's Legacy Survey for Space and Time, will continue to probe fainter and fainter dwarf galaxies. In addition, large aperture multi-object spectrographs have revealed the detailed and complex inner chemodynamical structure of dwarf galaxies. Beyond precise structural and kinematic properties of dwarf galaxies, modern observations allow for the separation of multiple stellar populations, detailed constraints on the dark matter density profiles, and hints of tidal disruption or stellar halos.

Of local dwarfs, the classical systems still remain the best studied and with most strongly constrained parameters. While extending 1-2 degrees across the sky, making deriving properties trickier [e.g., @mateo1998], these dwarfs contain large numbers of bright (giant) stars, allowing thousands of stars to be observed with deep photometry and spectroscopy [e.g., @tolstoy+2023; @pace+2020]. As a result, the determination of fundamental properties such as the position, size, orientation, distances, proper motions, line-of-sight (LOS) velocities, and velocity dispersions $\sigma_v$, are all relatively well constrained and converged among different studies today. However, ongoing research continues to redefine our understanding of the detailed structure of Milky Way satellites, hinting that these objects may be more complex than at first glance. 

Observations of dwarf galaxies have been the origin of several disputes or *small-scale* problems in cosmology [see review @bullock+boylan-kolchin2017]. For example, the mismatch between the expected number of dwarf galaxies from simulations and the observed abundance was known as the *missing satellites problem*. Additionally, a number of observations showed that many dwarf galaxies, but not all, possess dark matter "cores" as measured from circular velocity curves [e.g., @moore1994; @adams+2014; @oh+2015] or dwarf spheroidal [albeit less constraining, e.g., @walker+penarrubia2011; @read+walker+steger2019]. As a result, numerous alternative forms of dark matter have been advocated as solutions. *Warm* dark matter, relativistic in the early universe but cooler now, smooths out small-scale features and softens the cusps of dark matter halos. *Self-interacting* dark matter instead can form cores but also the cores can collapse into a density peak.

However, these tensions have eased as a result of better consideration of baryonic physics. For instance, understanding observational limitations and completeness [e.g., @kim+peter+hargis2018]. Recent hydrodynamic simulations in particular have shown that strong feedback can produce dark matter cores [e.g., @tollet+2016; @fitts+2017; @orkney+2021]. Altogether, the challenges and discussion around if $\Lambda$CDM correctly predicts the properties illustrates how dwarf galaxies represent ongoing 



## Stellar density profiles

Surface density profiles efficiently characterize the shape of a galaxy. At the most basic, fitting stellar densities provide properties such as the shape, location, size, and orientation of a dwarf galaxy. However, the details of a stellar density profile are essential for interpreting the total mass of dwarf galaxies, their assembly and dynamical history, and understanding correlations between dwarf galaxies [e.g., @herrmann+hunter+elmegreen2013, @querci+2025, @lee+2018]. 

Typically, four different surface density profiles are applied to dwarf galaxies: Exponential, Plummer, King, or Sérsic profiles [e.g., @munoz+18]. The exponential profile is perhaps the simplest, as a 1-parameter profile
$$
\Sigma_{\rm exp} = \Sigma_0\exp(-R / R_s)
$$

For a long time, exponential surface density profiles have been used as a description for spiral galaxies [@mateo1998; discussion below]. 

To fit globular cluster density profiles, @plummer1911 proposed a 1-parameter solution for a density polytrope^[where density and pressure are assumed to be correlated],
$$
\Sigma_{\rm Plummer} = \frac{M}{\pi R_h^2}\frac{1}{(1 + (R/R_h)^2)^2},
$$

where $M$ is the total mass and $R_h$ is the 2D half-light radius. Now mostly superseded by the King profile for globular clusters, the Plummer model is still a good fit to dwarf spheroidal [e.g., moskowitz+walker2020].

The @king1962 profile, also an empirical fit to globular clusters, is also often applied to dwarf galaxies, especially in older literature. Using three parameters, a core radius $R_c$, a truncation radius $R_t$, and a characteristic density, $\Sigma_0$, the  King profile is 
$$
\Sigma_{\rm King} = \Sigma_0\left(\frac{1}{\sqrt{1 + (R/R_c)^2}} - \frac{1}{\sqrt{1+(R_t/R_c)^2}}\right).
$$

In much of the older literature, $R_t$ was often called and interpreted as "tidal radius" as well, after the similar interpretation for globular clusters [e.g., @IH1995, @hodge1961]. 

Finally, the @sersic1963 profile represents a generalization of an exponential profile, and describes most galaxy light profiles well. Typically parameterized in terms of a half-light radius $R_h$, the density at half-light radius $\Sigma_h$ and a Sérsic index $n$, the profile's equation is
$$
\Sigma_{\rm S\acute ersic} = \Sigma_0 \exp\left[-b_n \,  \left((R/R_h)^{1/n} - 1\right)\right]
$$
where $b_n$ is a constant depending on $n$. $n=1$ provides an exponential profile and $n=4$ recovers @devaucouleurs1948's profile for elliptical galaxies. While a Sérsic profile is not always used for dwarf galaxies, @munoz+2018 advocate for using the Sérsic profile since the added flexibility allows for more profiles to be fit. 

While there are not clear theoretical explanations why any profile is best, the Exponential density profile is commonly assumed for dwarf spheroidal galaxies. @faber+lin1983 was one of the first to demonstrate that an Exponential density profile is a reasonable empirical  fit to dwarf galaxies, theorizing that dwarf spheroidal evolve from exponential disk galaxies, maintaining a similar light profile. Many later photometric works for dwarf spheroidal galaxies apply exponential fits, finding that exponential and king profiles perform similarly [@mateo1998; @IH1995; @mcconnachie+irwin2006; @cicuendez+2018]. As a result, many studies assume an exponential density profile in theoretical or observational modelling [e.g., @kormendy1985; @martin+2016 @MV2020a; @battaglia+2022; @kowalczyk+2013]. On the other hand, some authors fitting Sérsic profiles (often in addition to Exponential), often concluding that the added flexibility of a Sérsic produces better fits [@vanzee+barton+skillman2004; @munoz+2018; @wang+2019]. Note that the Sérsic indicies are typically $n \lesssim 1$, implying that most galaxies tend to be an exponential or slightly steeper. 

Beyond local group dwarf spheroidals, exponential densities appear to be common, but sometimes with modifications. In particular, many extragalactic dwarf Elliptical and blue compact / Irr galaxies galaxies are fit well with an exponential + central cusp / nuclear region [@caldwell+bothun1987,  @noeske+2003].





However, an exponential profile was often noted to not fit every dwarf galaxy, particular for more distant dwarf galaxies. Even starting from @aparicio1997,  and @graham+guzman2003  [for coma cluster dE], some people have noted deviations from this simple empirical rule. Additionally, @hunter+elmegreen; @herrmann+hunter+elmegreen2013; @herrmann+hunter+elmegreen2016; @lee+2018 note that at least in a photometric sample of more distant dwarf disky/Irr/blue compact dwarfs, that many or most dwarfs show two-part exponential profiles. Irr do represent a substantially different class though so it is unclear if these conclusions apply. Note that while diversity is also commented on here, dwarf galaxies with outer flattening profiles are often a minority in these studies. 



Altogether, while there is some natural variation in the density profiles of dwarf galaxies, an exponential is an excellent first-order approximation. Typically, deviations from exponentials are in a direction of steeper outer cutoff or changes to the inner slope of a dwarf galaxy (due to a nuclear region/etc.). As such, deviations in directions of density profiles becoming steeper in the very outer regions are unusual. Explaining the origin, similarity, and diversity of dwarf galaxy density profiles is a pressing question in our theories of dwarf galaxy formation and evolution.



- @caldwell+1992 mostly well fit outer exponential, inner deviations photometric M31

- @makarov+2012 isolated local volume dwarf galaxy with central deprivation but otherwise exponential.

- @moskowitz+walker2020, generalized plummer fits. Typically steeper and normal models similar, so not significant differences. 

- Figure out how to relate steeper ^^ model, King, and Exponential...

  

Irregardless of if all dwarf galaxies are indeed exponential, there is nevertheless diversity requiring interpretation in the density profiles of dwarf galaxies. 

## Sculptor and Ursa Minor: Hints of tidal signatures?

Sculptor and Ursa Minor may appear to be typical dwarf galaxies at first glance. [@tbl:scl_obs_props; @tbl:umi_obs_props] describe the present-day properties of each galaxy. Indeed, Sculptor as the first discovered dwarf galaxy, is often described as a "prototypical" dSph. However, both galaxies have a long history of speculation that they may be influenced by the Milky Way's gravity (see discussion). With *Gaia* data and spectroscopic followup, the detection of a density excess has become more unmistakable. Using @jensen+2024's algorithm (described briefly below), @sestito+2023a; @sestito+2023b report that a "kink" in the density profile, beginning around 30 arcmin for both Sculptor and Ursa Minor.  They spectroscopically followup some of the most distant stars, finding multiple members between 6-12 half-light radii from the centre of each dwarf (REF fig. XXX). If dwarfs initially begin with exponential profiles like Fornax, then these stars should be much rarer. 

Perhaps the most straighforward explanation for the density profiles of Sculptor and Ursa Minor are tidal interactions. @sestito+2023a; @sestito+2023b conclude that tides are likely the explanation. @Fig:scl_umi_fnx_vs_penarrubia reproduces their comparison, showing that Sculptor and Ursa Minor's density profiles deviate substantially from fornax's in the outskirts, and are well-described by the tidal model from @PNM2008. Our goal is to determine, assuming $\Lambda$CMD, if the effects of the Milky Way (or other satellites) may indeed create observable tidal signatures in these galaxies. If tides cannot explain these features, Sculptor and Ursa Minor may instead contain a extended stellar "halo" or second component of the galaxy—illustrating evidence of complex history or formation in each galaxy and forcing formation models to confront a diversity of density profiles in the local group.



![Sculptor and Ursa Minor match tidal models](./figures/scl_umi_vs_penarrubia.pdf){#fig:scl_umi_vs_penarrubia}

Figure: A plot of the surface density profiles of Sculptor, Ursa Minor, and Fornax scaled to their half-light radius and the density at half-light radius (data described in @sec:data). The lines represent the initial and final models from @PNM2008's model (Fig. 4 top left, stellar segregation of  0.20 moving on a highly eccentric orbit peri:apo = 1 : 100)).



| parameter        | value                                                        | Source |
| ---------------- | ------------------------------------------------------------ | ------ |
| $\alpha$         | $15.0183 \pm 0.0012^\circ$                                   | 1      |
| $\delta$         | $-33.7186 \pm 0.0007^\circ$                                  | "      |
| distance modulus | $19.60 \pm 0.05$                                             | 2      |
| distance         | $83.2 \pm 2$ kpc                                             | "      |
| $\mu_{\alpha*}$  | $0.099 \pm 0.002 \pm 0.017$ mas yr$^{-1}$                    | 3      |
| $\mu_\delta$     | $-0.160 \pm 0.002_{\rm stat} \pm 0.017_{\rm sys}$ mas yr$^{-1}$ | "      |
| LOS velocity     | $111.2 \pm 0.3\ {\rm km\,s^{-1}}$                            | 4      |
| $\sigma_v$       | $9.7\pm0.2\ {\rm km\,s^{-1}}$                                | "      |
| $R_h$            | $9.79 \pm 0.04$ arcmin                                       | 1      |
| ellipticity      | $0.37 \pm 0.01$                                              | "      |
| position angle   | $94\pm1^\circ$                                               | "      |
| $M_V$            | $-10.82\pm0.14$                                              | "      |

Table: Observed properties of Sculptor. References are: 1. @munoz+2018 Sérsic fits, 2. @tran+2022 RR lyrae distance, 3. @MV2020b, 4. @arroyo-polonio+2024. {#tbl:scl_obs_props  short="Observed Properties of Sculptor"}



| parameter        | value                                                        | Source |
| ---------------- | ------------------------------------------------------------ | ------ |
| $\alpha$         | $ 227.2420 \pm 0.0045$˚                                      | 1      |
| $\delta$         | $67.2221 \pm 0.0016$˚                                        | "      |
| distance modulus | $19.23 \pm 0.11$                                             | 2      |
| distance         | $70.1 \pm 3.6$ kpc                                           | "      |
| $\mu_\alpha*$    | $-0.124 \pm 0.004 \pm 0.017$ mas yr$^{-1}$                   | 3      |
| $\mu_\delta$     | $0.078 \pm 0.004_{\rm stat} \pm 0.017_{\rm sys}$ mas yr$^{-1}$ | "      |
| LOS velocity     | $-245.9 \pm 0.3_{\rm stat} \pm 1_{\rm sys}$ km s$^{-1}$      | 4      |
| $\sigma_v$       | $8.6 \pm 0.3$                                                | "      |
| $R_h$            | $11.62 \pm 0.1$ arcmin                                       | 1      |
| ellipticity      | $0.55 \pm 0.01$                                              | "      |
| position angle   | $50 \pm 1^\circ$                                             | "      |
| $M_V$            | $-9.03 \pm 0.05$                                             | "      |

Table: Observed properties of Ursa Minor. References are: (1) @munoz+2018 Sérsic fits, (2) @garofalo+2025 RR lyrae distance, (3) @MV2020a, (4) @pace+2020, average of MMT and Keck results.  {#tbl:umi_obs_props  short="Observed Properties of Ursa Minor"}



# Simulating tidal effects

Since the discovery of dwarf galaxies, a large body of work has speculated, considered, or simulated tidal effects. While pre-*Gaia* work often did not know the specific orbits of dwarf galaxies, the theory for satellites serves as an excellent framework for understanding specific systems. 

Simulating dwarf galaxies accurately in a cosmological context remains a substantial challenge. Cosmological simulations can now predict the overall abundance of larger mass dwarf galaxies [e.g., ] and predict generally effects of tides [e.g., @riley+2024]. But, dwarf galaxies are often barely resolved leading to numerical disruption of dwarf galaxies [@santos-santos+2024]. The highest resolution simulation of a Milky Way dark matter halo, the Aquarius project [@springel+2008], acchieved a DM resolution of $1.712 \times10^3 \Mo$, enough to barely resolve Sculptor like halos (see methods).  On the other hand, idealized simulations are able to reach high resolution and numerical convergence for a single dwarf galaxy. For instance, our simulations later are about 3x higher resolution than Aquarius but with a fraction of the computational cost (400x less particles). Idealization does make numerous simplifying assumptions: neglecting mergers, cosmological context and evolution, mass assembly, and often baryonic physics. We use idealized simulations here which are more powerful in accurately assessing tidal effects after infall, given that the idealizations do not impact these galaxies's recent history too much. As such, a long history of work has resolved to understand the tidal evolution of dwarf galaxies using these *idealized* simulations. 

Some of the earliest theory work on tidal stripping of dwarf galaxies originate from @oh+lin+aarseth1995; @piatek+pryor1995; @moore+davis1994; @johnston+spergel+hernquist1995. Already, these works used similar techniques as we continue to use today and layed the general foundation for interpreting tidal mass loss of galaxies. Unfortunately, most models at this point ar King models (@king1966, different from king density profiles discussed earlier) and assume smaller masses than we would today since the mass-to-light ratios of dwarf galaxies were in question. However, the qualitative results of what happens during tidal mass loss remains. Mass is predominantly lost through the lagrange points at pericentres. This is essentially an energy truncation. Additionally, the inner dwarf density profile may remain in tact, but the outer profile can begin to flatten with tidal effects. Additionally, the velocity dispersion in the centre decreases. Tidal arms are formed from particles leaving in either L1 or L2, which enter into slightly lower energy, interior, leading orbits, or higher energy, exterior, trailing orbits. 

More recent works include @read+2006; @bullock+johnston2005; @PNM2008; @penarrubia+2009; @klimentowski+2009; @errani+2023a; @fattahi+2018; @stucker+2023; @wang+2017.

With precise orbits and a better understanding of the Milky Way potential and system, more recent work began to directly probe the dynamical histories of individual dwarf galaxies (although early work began this for  Sagittarius / etc). Examples include @iorio+2019 for Sculptor, @borukhovetskaya+2022; @dicintio+2024 for Fornax; @borukhovetskaya+2022a for Antlia II. Our goal is to apply a similar framework to Sculptor and Ursa Minor.

## Simulating large gravitational systems: The N-body method

Modelling gravitational evolution is essential for understanding properties of galaxies. Perhaps the simplest method to compute the evolution of dark matter is through *N-body simulations*. A dark matter halo is represented as a large number of dark matter particles (bodies). Each body represents a monte-carlo sample from the underlying matter distribution. However, galaxies are often assumed to be *collisionless*—particles are not strongly affected by small-scale gravitational encounters. In contrast, star clusters evolve differently as star-star gravitational collisions changes the dynamics of the system. While we use individual gravitating bodies in N-body simulations, the Newtonian gravitational force is softened to be a Plummer sphere, so encounters closer than a softening length do not substantially impact the dynamical evolution.

Naively, the newtonian gravitational force requires adding together the forces from each particle at each particle, causing a computational cost that scales quadratically with the number of particles, or $O(N^2)$. With this method, simulating a large number of particles, such as 10^6, would require about 10^12 force evaluations at each time step, making cosmological and high-resolution studies unfeasible. However, only long-range gravitational interactions tend to be important for CMD, so we can utilize the *tree method* to compute the gravitational force vastly more efficiently.

The first gravitational tree code was introduced in @barnes+hut1986, and is still in use today. We utilize the massively parallel code *Gadget 4* [@gadget4].  Particles are spatially split into an *octotree*. The tree construction stars with one large node, a box containing all of the particles. If there is more than one particle in a box/node , the box is then divided into 8 more nodes (halving the side length in each dimension) and this step is repeated until each node only contains 1 particle. With this heirarchichal organization, if a particle is sufficiently far away from a node, then the force is well approximated by the force from the centre of mass of the node. As such, each force calculation only requires a walk through the tree, only descending farther into the tree as necessary to retain accuracy. The total force calculations reduce from $O(N^2)$ to $O(N\,\log N)$, representing orders of magnitude speedup. Modern codes such as *Gadget* utilize other performance tricks, such as splitting particles across many supercomputer nodes, efficient memory storage, adaptive time stepping, and parallel file writing to retain fast performance for large scale simulations, forming the foundation for many cosmological simulation codes. 

## Why do particles leave?



## Break and tidal radii

Analytic approximations are powerful tools to generally determine the influence of tides given a dwarf galaxy and a host. A variety of "tidal radii" exist in the literature, we focus on two simple forms, the Jacobi radius and the break radius. 

The **Jacobi radius** represents the approximate radius where stars become unbound for a galaxy in a circular orbit around a host galaxy. Calculated from an approximation of the location of the L1 and L2 lagrange points, the Jacobi radius is where the mean density of the dwarf galaxy is three times the mean interior density of the host galaxy at pericentre, or
$$
3\bar \rho_{\rm MW}(r_{\rm peri}) = \bar \rho_{\rm dwarf}(r_J).
$$
If $r_J$ occurs within the visible extent of a galaxy, we should expect to find unbound, *extratidal* stars. While most valid for circular orbits, assuming $r_{\rm peri}$ for the host-dwarf distance works as most stars are lost near pericentre. 

We use the **break radius** as defined in @penarrubia+2009, marking where the galaxy is still in disequilibrium. The break radius $r_{\rm break}$ is proportional to the velocity dispersion $\sigma_v$ and time since pericentre $\Delta t$, 
$$
r_{\rm break} = C\,\sigma_{v}\,\Delta t
$$
where the scaling constant $C \approx 0.55$ is a fit. $r_{\rm break}$ describes where the dynamical timescale is longer than the time since the perturbation, i.e. the radius within which the galaxy should have dynamically relaxed.  As illustrated in @fig:toy_profiles, for an idealized model with exponential stars in a NFW halo, after some time past pericentre, the stellar component is smooth but contains a change in slope around $r_{\rm break}$, where the radial velocities of the stars becomes positive as they still readjust to a new equilibrium.



![Break radius validation](figures/idealized_break_radius.pdf){#fig:idealized_break_radius}

Figure: The break radius of the simulations is set by the time since pericentre.  The initial conditions are Sculptor-like, exponential stars embedded in NFW, evolved in @EP2020 potential with a pericentre of 15 kpc and apocentre of 100 kpc. In this model, the jacobi radius is close to the break radius. See @sec:methods for a description of our simulation setup.

## Tidal evolution

While the break and jacobi radii help with understanding when tidal effects matter, we still need to know what these effects do and how a galaxy evolves. 

As a galaxy evolves in a tidal field, several changes happen:

1. *Mass is lost*. In particular, particles and stars on weakly bound orbits are most likely to be removed by tides. 
2. *Steams form*. If tides are strong enough, unbound mass becomes part of tidal tails, evolving along a similar orbit but leading or trailing the galaxy.
3. *Mass redistributes*. Bound mass of the galaxy redistributes. As mentioned in the last section, this is visible as a wave of outward moving material, with the outermost material reaching equilibrium last. 
4. *A new equilibrium*. With mass loss, the gravitational potential decreases, resulting in a more compact dark matter halo and stars which adiobadically expand to a larger scale radius.

Assuming that galaxies are spherical, isotropic, and evolve in a constant tidal field, most dwarf galaxies should evolve along similar tidal tracks. @EN2021 derive tidal tracks for dwarf galaxies, showing that NFW halos all evolve along a similar trajectory in terms of $r_{\rm max}$ and $v_{\rm max}$. 

$$
\frac{v_{\rm max}}{v_{\rm max, 0}} = 
2^\alpha 
\left(\frac{r_{\rm max}}{r_{\rm max, 0}}\right)^{\beta}\left[1 + \left(r_{\rm max} / r_{\rm max, 0}\right)^2\right]^{-\beta}
$$
where $\alpha=0.4$, $\beta=0.65$ are empirical fits. As illustrated in @fig:tidal_tracks, this formula works for both circular and elliptical orbits and is independent of the initial subhalo size or distance to the host. 

![Tidal tracks of dwarf galaxies](/Users/daniel/Library/Application Support/typora-user-images/image-20250715095615423.png){#fig:tidal_tracks width=80%}

Figure: Tidal tracks of dwarf galaxies, the logarithm of maximum circular velocity and radius of relative to the initial conditions for satellites on a variety of orbits. Almost all satellites follow the tidal track suggested by @EN2021. Adapted from fig. 6 of @EN2021.

# Thesis outline

In this thesis, our goal is to review the evidence for an extended density profile in Ursa Minor and Sculptor, to assess the impact of tidal effects on each galaxy, and to discuss possible interpretations for the structure of these galaxies. 

In chapter 2, we describe how we compute observational density profiles from @jensen+2024. In chapter 3, we review our simulation methods. Next, we present our results for the tidal stripping of Sculptor and Ursa Minor in Chapter 4. We discuss our results, limitations, and implications in Chapter 5. Finally, in Chapter 6, we summarize this work and discuss future directions for similar work and the field of dwarf galaxies. 



