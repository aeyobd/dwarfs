Dwarf galaxies host, in many ways, the most extreme galactic environments in the universe. These galaxies are typically defined to be fainter than the Large Magellanic Cloud (LMC),  with $M_V \gtrsim -18$ or similarly $M_\star \lesssim 10^9 M_\odot$ [e.g.,  @mcconnachie2012; @bullock+boylan-kolchin2017]. Because the galaxy luminosity function increases towards fainter objects, dwarfs are the most numerous of galaxies [e.g., @blanton+2005; @mao+2021]. Dwarf galaxies are also highly *dark-matter dominated*, with mass to light ratios which may exceed 1000 $M_\odot/ L_\odot$ [implying $\sim1000$ times more dark matter than stellar mass, e.g., @simon+geha2007; @hayashi+2023].  

Except for the Magellanic Clouds, most dwarf galaxy satellites of the Milky Way (MW) are *quenched*, with little to no recent star formation [e.g., @weisz+2014]. Indeed, most faint MW satellites contain stellar populations which are *relics* from the early universe, consisting of many of the oldest and most metal-poor stars known [@simon2019]. Understanding the properties of dwarf galaxies thus has implications across astronomy, from small-scale cosmological structure formation to the origins of the first stars. 

In this Chapter, we first describe the general observed properties of local dwarf galaxies. Next, we summarize our understanding of the cosmological origin of dwarf galaxies. We later review recent advancements and pending questions concerning dwarf galaxies, and introduce the puzzle posed by the extended stellar density profiles of Sculptor and Ursa Minor. Then, we discuss the theory of tidal evolution. We end with a brief roadmap to the remainder of this dissertation. 

# Observations of dwarf galaxies

Dwarf galaxies have long raised conundrums for theories of galaxy formation. The discovery of Fornax and Sculptor in 1938 [@shapley1938][^lmc_discovery], with no known analogues at the time, already presented an enigma. H. Shapley presented these dwarfs as a new type of *stellar system* resembling the Magellanic Clouds and globular clusters but did not attempt to speculate on their nature. While dwarf galaxies were soon understood to be galaxies based on the inferred luminosities and sizes, their exact nature remained unclear for decades [e.g.,@hodge1971; @gallagher+wyse1994]. 

The earliest spectroscopic work hinted that dwarf galaxies may contain substantial amounts of dark matter. From velocity dispersion measurements for dwarf spheroidal (dSph) galaxies,  inferred mass-to-light ratios were at least 10 times larger than for globular clusters [e.g., @aaronson1983; @aaronson+olszewski1987]. While uncertain initially, these values were later corroborated with larger and more precise samples [e.g., @hargreaves+1994]. At the time, several theories were proposed to explain these unusually high mass-to-light ratios. Examples include: ongoing tidal disruption inflating inferred velocity dispersions [e.g., @kuhn+miller1989], the presence of massive central black holes [e.g., @strobel+lake1994], or modified theories of gravity [@milgrom1995].  Over time, a consensus developed where the high mass-to-light ratios of dwarf galaxies were due to the presence of a massive dark matter halo [e.g., @dekel+silk1986; @wechsler+tinker2018]. Since then, the properties of dwarf galaxies have played an increasingly important role in our understanding of the clustering of dark matter on small scales [e.g., @bullock+boylan-kolchin2017; @sales+2022].  

Today, a common definition for a (dwarf) galaxy is that of a gravitationally bound stellar system with dark matter.^[Or, more generally, systems inconsistent with Newtonian dynamics of visible matter alone [@willman+strader2012].] In contrast, star clusters (like globular clusters) have no clear evidence for dark matter. The boundary between these two classes blurs for faint, compact stellar associations. Systems with characteristics of both globular clusters and dwarf galaxies are known as "ambiguous" systems [e.g., @smith+2024]. 

Dwarf galaxies span a large range of sizes, luminosities, and morphologies. Broadly, there are three classes of dwarf galaxies based on luminosity. Local **bright dwarfs** with magnitudes $-14 \gtrsim M_V \gtrsim  -18$ or stellar masses^[Assuming stellar mass-to-light ratio of $1 \Mo / L_{\odot}$, may be $\sim 2 \Mo/L_\odot$ for older populations.] $3\times10^7\,\Mo \lesssim M_\star \lesssim 10^9\,\Mo$, often exhibit irregular morphologies and recent star formation.  @fig:galaxy_images shows the Large Magellanic Cloud (LMC) as an example of an irregular, bright dwarf galaxy, where most stars are in a rotationally-supported disk (seen nearly face-on) with a prominent bar. **Classical dwarf spheroidals**[^dsph_suffix]  occupy intermediate luminosities ( $-7.7 \gtrsim M_V  \gtrsim -14$ or $10^5\,\Mo \lesssim M_\star \lesssim 3\times10^7\,\Mo$). Typically, these systems are old, non-star-forming, gas-poor, and elliptical. All Milky Way satellites discovered before digital sky surveys are classicals, and these systems remain among the best studied dwarfs. The **ultra-faint dwarfs** occupy the very faintest magnitude range ($M_V \gtrsim -7.7$ or $M_\star \lesssim 10^5\,\Mo$). These galaxies have minuscule stellar masses, typically compact sizes, and very metal-poor stellar populations [see the review by @simon2019]. Altogether, known dwarf galaxies span more than 15 magnitudes, or over 6 decades in stellar mass.

Most well-studied dwarf galaxies lie in the vicinity of the Milky Way, the *Local Group* of galaxies. The Local Group is defined as the group consisting of galaxies within $\sim 1$ Mpc from the MW-Andromeda centre [e.g., @mcconnachie2012 and references therein]. Today, we know that the Milky Way system is teeming with dwarfs, many of which are satellites of either the MW or Andromeda (M31). [@fig:mw_satellite_system] shows the MW satellite system, including dwarf galaxies, globular clusters, and ambiguous systems. This nearby population of dwarf galaxies is amenable to resolved studies aimed at investigating their detailed history and structure. 

[^lmc_discovery]: Technically, the Large and Small Magellanic Clouds (LMC, SMC) are also classified as dwarf galaxies, but these were likely always known to humans at southern latitudes. 
[^dsph_suffix]: While formally the dwarf galaxy names we discuss contain "dwarf spheroidal" (dSph), e.g., Sculptor dSph, we omit this suffix for brevity. Additionally, the 12 classical dwarf satellites of our Galaxy are (in order of decreasing luminosity) Sagittarius, Fornax, Leo I, Sculptor, Antlia II, Leo II, Carina, Draco, Ursa Minor, Canes Venatici I, Sextans I, and Crater II. Antlia II, Crater II, and Canes Venatici I are the only post-digital sky surveys additions.



![Images of dwarf galaxies](figures/galaxy_pictures.png){#fig:galaxy_images width=100%}

Figure: Images^[Created with hips2fits (https://alasky.cds.unistra.fr/hips-image-services/hips2fits), a service provided by CDS.] of the LMC [Digitized Sky Survey II, @lasker+1996], Fornax [DES DR2, @abbott+2021], Sculptor (DES DR2),  and Ursa Minor [UNWISE, @lang2014; @meisner+lang+schlegel2017; @meisner+lang+schlegel2017a, with \textit{Gaia} point sources overplotted]. The grey ellipse represents the half-light radius for the three dwarf spheroidals, and the luminosity is derived from the absolute V-band magnitude of each galaxy.



![The on-sky distribution of Milky Way satellites](figures/mw_satellites_1.jpg){#fig:mw_satellite_system width=100%}

Figure: The location of MW satellites on the sky. We label the classical dwarf galaxies (green diamonds), fainter dwarfs (blue squares), globular clusters (orange circles), and ambiguous systems (pink open hexagons). Globular clusters are more centrally concentrated, but dwarf galaxies are preferentially found away from the MW disk. Sculptor and Ursa Minor are highlighted as two dwarfs we study later. The background image is from ESA/Gaia/DPAC.^[https://www.esa.int/ESA_Multimedia/Images/2018/04/Gaia_s_sky_in_colour2] Dwarf galaxies (confirmed), globular clusters, and ambiguous systems are from the @pace2024 catalogue (version 1.0.3). 

# Dwarf galaxies in a cosmological context

We only understand a fraction of the universe's composition. The leading theory of cosmology,  Lambda Cold Dark Matter (\LCDM{}), posits that the universe is composed of about 68% dark energy ($\Lambda$), 27% dark matter (DM), and 5% baryons[^baryons]  [@planckcollaboration+2020]. While the composition of dark matter and dark energy remains elusive, we know their general properties. Dark energy drives the acceleration of the expansion of the universe on large scales. We do not discuss dark energy here—it does not substantially affect the Local Group today. Dark matter, instead, makes up the vast majority of mass in galaxies. Typically, galaxies have baryonic-to-dark matter ratios of between 1:5 to beyond 1:1000 for faint dwarf galaxies [e.g., @hayashi+2023].

[^baryons]:  In a classic astronomer's corruption of jargon, *Baryons* here means baryons and leptons, i.e., protons, neutrons, and electrons.



In \LCDM{}, dark matter is assumed to interact only gravitationally. Light and matter pass through dark matter unimpeded---in this sense, dark matter is transparent. Dark matter is also assumed to be *cold*, i.e., with typical velocities much smaller than the speed of light in the early universe. If dark matter is cold, then it should condense on all scales, from the size of galaxy clusters to smaller than the faintest dwarf galaxies. Implications of dark matter properties include cosmological structure, galaxy formation, and galaxy interactions. 

## Structure formation in \LCDM{}

The very early universe was almost featureless. Our earliest observations of the universe stem from the cosmic microwave background (CMB)---displaying a nearly uniform, isotropic blackbody emission [e.g., @ryden2016]. But tiny perturbations in the CMB, temperature fluctuations of 1 part in 100,000, reveal the underlying seeds of large-scale cosmological structure. In an expanding universe, gravitational instability makes CDM overdensities grow and collapse hierarchically onto larger structures. Initially, baryonic matter was coupled to radiation and resisted collapse. Dark matter, only influenced by gravity, freely collapsed into the first structures. Mass perturbations sufficiently small and overdense become self-gravitating structures, known as *halos* [e.g., @galaxiesbook]. After recombination, where electrons combined with atomic nuclei to form atoms, baryons decoupled from radiation and fell into dark matter halos, where they condensed at the centre through radiative energy losses. The densest pockets of baryons later formed the first stars and galaxies.

Dark matter halos and their associated galaxies rarely evolve in isolation. Instead, \LCDM{} structure formation is *hierarchical*. Small dark matter halos collapse first and hierarchically merge into progressively larger halos [e.g., @white+rees1978; @blumenthal+1984; @white+frenk1991]. Hierarchical assembly is evident through the large-scale structure of the universe, remnants of past mergers within the Milky Way, and tidal disruption of dwarf galaxies and their streams around nearby galaxies.

Small-scale structure formation is sensitive to deviations from \LCDM{} [e.g., @bechtol+2022]. One key prediction of \LCDM{} is that mass perturbations are expected to exist on all scales, and are largest on the smallest scales, so we would expect the formation of halos on all scales. Many alternative models, such as warm dark matter, may smooth out small-scale features and reduce the abundance of small halos or change their structure [e.g., @lovell+2014]. Dwarf galaxies, which occupy the smallest dark matter halos, are promising probes into the behaviour of dark matter on small scales.

## The structure of cold dark matter halos {#sec:NFW}

In \LCDM{} cosmological simulations, dark matter halos are remarkably self-similar. \citet{NFW1996, NFW1997} observe that the spherically-averaged density profiles $\rho(r)$ are universally well described by a two-parameter law,
$$
\rho/\rho_s= \frac{1}{(r/r_s)(1+r/r_s)^2},
$$ {#eq:nfw}
where $r_s$ is a scale radius and  $\rho_s$ a scale density. This profile, known by the author's initials NFW, has shown remarkable success at describing \LCDM{} halos across several orders of magnitude in mass. NFW profiles are *cuspy*, where the density rises like $\rho \sim 1/r$ at small radii $r \ll r_s$. The steepness of the density profile increases gradually with radius, and at large radii the density falls off like $\rho \sim 1/r^3$. The solid blue curve in @fig:nfw_density shows an example NFW halo.

The total mass of an NFW profile formally diverges, so halo masses are conventionally defined using an overdensity criterion. The virial mass, $M_{200}$, is defined as the mass within a radius, $r_{200}$, containing a mean enclosed density 200 times[^200] the critical density of the universe:
$$
M_{200} =200\,\frac{4\pi}{3} \ r_{200}^3\ \rho_{\rm crit}, \qquad {\rm where} \quad \rho_{\rm crit}(z) = 3H(z)^2 / 8\pi G,
$$
and $H(z)$ is the Hubble constant, which depends on redshift. Another way of characterizing NFW halos is through the concentration parameter, $c=r_{200} / r_s$, which describes how the characteristic radial scale of the halo compares to the virial radius.  Using this parameter, the scale density is a function of $c$ alone, $\rho_s = (200/3)\,\rho_{\rm crit} c^3 / [\log(1+c) - c/(1+c)]$ [@NFW1996]. 

An equivalent, alternative characterization of NFW halos uses their circular velocity profiles. The circular velocity, $\vcirc(r) = \sqrt{G M(r) / r}$, reaches a maximum $\vmax$ at radius  $\rmax \approx 2.16258\,r_s$. $\vmax$ and $r_{\rm max}$,  like $M_{200}$ and $c$, fully specify an NFW halo.



The two parameters of an NFW profile are not independent. Lower-mass dark matter halos typically collapse earlier, when the universe was denser. As a result, low mass subhalos tend to be more concentrated [e.g., @NFW1997]. The relationship between $M_{200}$ and c, or the mass-concentration relation, describes the mean trend of concentration with mass or, equivalently, the dependence of $\vmax$ on $\rmax$ [e.g., @bullock+2001; @ludlow+2014]. The left panel of @fig:smhm illustrates the present-day mass-concentration from @ludlow+2016. While concentration tends to decrease with increasing mass, the relation has substantial scatter. Other parameters, such as the halo spin or shape, may affect the scatter of the mass-concentration relation, but their effect is typically expected to be small [@navarro+2010; @dicintio+2013; @dutton+maccio2014]. 



[^200]: For the collapse of a uniform spherical density, the virialized overdensity would be $\Delta = 18\pi^2\approx 178$ for a critical universe $\Omega_m = 1$. This is commonly rounded to $\Delta = 200$. While this parameter may be closer to $\Delta \approx 100$ for our universe, $\Delta$ also increases with redshift [see, e.g., eq. 6 from @bryan+norman1998]. 

![Example dark matter and stellar density profiles](figures/example_density_profiles.png){#fig:nfw_density width=252pt}

Figure: Density profiles in log 3D density versus log 3D radius for stars and dark matter in a Fornax-like galaxy. The dark matter is more extended and massive than the star across the entire galaxy. This galaxy has a stellar mass $M_\star \approx 2.5\times10^7\,\Mo$ with half-light radius $0.65\,\kpc$ and 2D exponential stars [based on, @munoz+2018; @woo+courteau+dekel2008]. The corresponding cosmological-mean halo has $\vmax=40\,\kpc$ and $\rmax=8\,\kpc$,  or $M_{200} = 1\times10^{10}\,\Mo$ and $c=12.5$. 

![Cosmological mass-concentration and stellar mass-halo mass relations](figures/cosmological_means.pdf){#fig:smhm}

Figure: **Left** The NFW halo concentration $c=r_{200} / r_s$ as a function of virial mass $M_{200}$. The solid line with 1-$\sigma$ shaded region is the mass-concentration relation from @ludlow+2016 for $z=0$. **Middle**: Equivalent to the left except in terms of the halo maximum circular velocity, $\vmax$, and radius where the velocity is maximized, $\rmax$.  **Right** Stellar mass (top) as a function of maximum circular velocity. The solid line with the 1-$\sigma$ shaded region is the relation from @fattahi+2018 with scatter points simulated central galaxies from \apostle{} in @fattahi+2018. The pink star illustrates the location of the Fornax galaxy, whose density profiles are shown in @fig:nfw_density. 

## Galaxy formation in \LCDM{} {#sec:galaxy_formation}

The observed abundance of galaxies may be compared with the abundance of \LCDM{} halos to derive constraints regarding which galaxies inhabit which halos. One simple technique, dubbed "abundance matching," assumes a tight relation between the stellar mass of a galaxy and the mass of the halo it inhabits [@li+white2009; @moster+naab+white2013]. 

The right-hand panel of @fig:smhm shows the stellar mass versus halo mass relation (SMHM, with halo mass represented by $\vmax$) predicted by \LCDM{} cosmological hydrodynamical simulations of Local Group analogues from the \apostle{} project [@sawala+2016].^[\apostle{} simulated Local Group analogues in a \LCDM{} cosmological context with the hydrodynamical setup from the \eagle{} simulations [@crain+2015; @schaye+2015].] While there is some scatter, the range of predicted $\vmax$ is fairly narrow across $\sim 4$ decades in stellar mass. This figure indicates that the SMHM relation becomes increasingly steep in the dwarf galaxy regime—many dwarf galaxies are formed in halos of similar mass. Because lower mass galaxies have shallower potential wells, the energetic output of evolving stars and possibly supermassive black holes (i.e., "feedback") becomes more effective at removing gas. Reionization additionally suppresses late star formation in the faintest galaxies. As a result, the resulting stellar mass of a dwarf galaxy is highly sensitive to the details of halo assembly and evolution.

In \LCDM{} galaxy formation, the majority of mass in a dwarf galaxy comes from the extended dark matter halo. @fig:nfw_density shows an example exponential stellar component for the Fornax dwarf galaxy with its surrounding dark matter NFW halo [with parameters matching the @ludlow+2016; @fattahi+2018; relations]. Where the stars are densest, the dark matter remains nearly an order of magnitude higher in density. Stars make a small contribution to the gravitational structure of dwarf galaxies---indeed, stars are reasonably approximated as tracer particles of the underlying dark matter halo. In addition, the stellar component is typically confined to the central regions of the dark matter halo. 

Several factors affect the SMHM trend, including environment, assembly history, tidal effects, and the details of galaxy formation. For example, effects like ram-pressure stripping (removal of gas in the dwarf galaxy due to pressure from the host's circumgalactic medium) and tidal removal of gas cause star formation to quench [e.g., @christensen+2024]. Additionally, the time of formation (relative to reionization) can influence the resulting stellar content [@kim+2024]. Finally, Galactic tides reduce both the dark matter and stellar mass but in different amounts, adding additional scatter to the SMHM trend for satellites [e.g., @PNM2008; @fattahi+2018]. Understanding the effects of tides on Local Group dwarf galaxies may help us understand where and how these galaxies formed in a cosmological context. 

## Challenges and questions concerning dwarf galaxies

Observations of dwarf galaxies have been the origin of several disputes or *small-scale* problems for \LCDM{} [see reviews by @bullock+boylan-kolchin2017; @sales+2022]. For example, the mismatch between the number of dwarf galaxies and the predicted abundance of \LCDM{} halos has been known as the *missing satellites problem*. Additionally, several observations suggest that some dwarf galaxies, although not all, possess dark matter "cores,"  [e.g., @moore1994; @adams+2014; @oh+2015; @walker+penarrubia2011; @read+walker+steger2019],^[In detail, gas-phase rotation curves are better able to differentiate between cores and cusps, whereas stellar kinematics is less constraining.] contrary to the expectation from \LCDM{} of "cuspy" inner dark matter profiles [@NFW1996; @NFW1997]. As a result, alternative forms of dark matter have been advocated as solutions, such as Warm or Self-Interacting Dark Matter.

However, some of these tensions have eased as a result of improved understanding of baryonic physics. For example, recent hydrodynamic simulations, in particular, have shown that strong feedback can produce dark matter cores [e.g., @navarro+eke+frenk1996, @tollet+2016; @fitts+2017; @benitez-llambay+2019; @orkney+2021]. Several open questions remain, concerning, e.g., the nearly planar distribution of luminous Milky Way satellites, the details of the sizes and rotation curves of dwarf galaxies, and the existence and nature of stellar halos in dwarf galaxies [e.g., @sales+2022]. Altogether, the numerous past and ongoing challenges for \LCDM{} in the dwarf galaxy regime illustrate the opportunity for dwarf galaxies to our the understanding of galaxy formation and dark matter physics.

# The structure of nearby dwarf galaxies

## The *Gaia* mission

Since Local Group dwarfs are nearby, they are resolved into individual stars, and therefore, we can study these galaxies on a star-by-star basis. As a result, it is possible to measure the 3D velocity and position of a star if we can measure its position, distance, line-of-sight (LOS) velocity, and proper motion. Unfortunately, determining distances and full 3D velocities is challenging. The most direct measurement of distance, the parallax, requires precise tracking of a star's sky position across a year. And while line-of-sight (LOS) velocities are relatively easily determined from spectroscopy, tangential velocities, derived from proper motions and distances, are much more challenging. Typically, measuring proper motions requires precise ($\ll$ arcsecond) determinations of small changes in a star's position over baselines of years to decades. The full 6D position and velocity information for stars has, until recently, been known for only a handful of stars.

Launched in 2013, *Gaia* is a space-based, all-sky survey telescope situated at the Sun-Earth L2 Lagrange point [@gaiacollaboration+2016]. *Gaia* has redefined astrometry, providing photometry, positions, proper motions, and parallaxes for over 1 billion stars [@gaiacollaboration+2021]. While *Gaia* completed its space-based mission in 2025, two further data releases are still expected. 

Determining absolute parallax measurement is facilitated by the observation that stars in different regions of the sky are affected by parallax motion with different phases. By imaging two regions separated by 106.5 degrees on the same focal plane, *Gaia* measures changes in the relative positions of stars across small and large angles. Combining measurements from multiple epochs across several years, an absolute all-sky reference frame is derived from which parallax and proper motions are calculated. In addition to astrometry, *Gaia* measures photometry in the wide *G* band (330--1050nm) and colours from the blue photometer (BP, 330--680 nm) and red photometer (RP, 640--1050 nm). *Gaia* additionally provides low-resolution BP-RP spectra and radial velocity measurements of bright stars [of magnitudes $G_{\rm RVS} < 16$, @gaiacollaboration+2016]. For our work, *Gaia*'s most relevant measurements are $G$ magnitude, $G_{\rm BP} - G_{\rm RP}$ colour, $(\alpha, \delta)$ position, and $(\mu_{\alpha*}, \mu_\delta)$ proper motions.[^pmra_cosdec]

## *Gaia*'s impact on Milky Way studies

*Gaia* has revolutionized our understanding of Milky Way structure. For example, the 6D dynamical measurements and metallicities of MW stars led to the (re)discovery of past mergers or Milky Way building blocks like *Gaia*-Sausage Enceladus [e.g., @helmi+2018; @belokurov+2018; but see also @meza+2005], out-of-equilibrium structures like the *Gaia* snail [e.g., @antoja+2018], and dynamical effects of the Milky Way's spiral arms and the bar in the solar neighbourhood [@hunt+vasiliev2025 and references therein]. In the Milky Way halo, *Gaia* has helped find and constrain numerous stellar streams  [@ibata+malhan+martin2019; @bonaca+price-whelan2025]. Altogether, *Gaia* has revealed the hierarchical formation and complex, evolving structure of our own Galaxy.

For Milky Way satellites, *Gaia* has improved orbital analysis and facilitated robust stellar membership determinations. Before *Gaia*, few galaxies had precisely measured proper motions [e.g., using the Hubble Space Telescope, @piatek+2005; @sohn+2017]. *Gaia* allowed for some of the first systematic and precise determinations of Milky Way satellite proper motions [@pace+li2019; @MV2020a]. While the proper motion uncertainty of a typical dwarf member star is often large, by combining the proper motions of 100s or 1000s of stars from *Gaia*, precise average proper motion measurements can be determined, sometimes only limited by *Gaia*'s systematic error floor  [e.g., @MV2020a]. Proper motions have thus ushered in a new era for MW satellite dynamical studies, where we can derive precise orbits for any satellite, assuming a given MW potential. In addition, *Gaia* helps establish membership by filtering contaminating MW foreground stars. By measuring parallaxes and/or proper motions, many more background and foreground stars can be classified as non-members [e.g., @battaglia+2022; @jensen+2024].



[^pmra_cosdec]: The proper motions $\mu_\alpha$ and $\mu_\delta$ are the apparent rates of change in right ascension, $\alpha$, and declination, $\delta$, typically in units of milli-arcsecond (mas) per year. $\mu_{\alpha*} = \mu_\alpha \cos \delta$ corrects for projection effects in $\alpha$.



## Dwarf galaxy light profiles {#sec:exponential_profiles}

Projected luminosity/stellar density profiles efficiently characterize the radial structure of a galaxy. At its most basic, light profiles synthesize properties such as the shape, size, and orientation of a dwarf galaxy. In addition, the details of a stellar density profile can help interpret a galaxy's assembly and dynamical history [e.g., @penarrubia+2009; @lee+2018; @querci+2025]. Note that for resolved galaxies, these profiles are expressed in stellar count densities instead of surface brightness.

Four different surface density laws are frequently used to parameterize dwarf galaxy profiles: Exponential, Plummer, King, or Sérsic profiles [e.g., @munoz+2018]. The exponential profile is perhaps the simplest, defined in terms of the central surface density, $\Sigma_0$, and projected scale radius, $R_s$:
$$
\Sigma_{\rm exp} = \Sigma_0\exp(-R / R_s).
$$ {#eq:exponential_law}
This profile is also often applied to the radial light distribution of galaxy disks [@devaucouleurs1959a; @freeman1970; @kent1985].

To fit globular cluster density profiles, @plummer1911 proposed a profile based on a self-gravitating polytrope,^[where density and pressure are assumed to be related by a power law]
$$
\Sigma_{\rm Pl} = \frac{\Sigma_0}{(1 + (R/R_h)^2)^2},
$$
where $\Sigma_0$ is the central surface density and $R_h$ is the projected half-light radius. Now mostly superseded by the King profile for globular clusters, the Plummer model is still a good fit to many dwarf spheroidals [e.g., @moskowitz+walker2020].

The @king1962 profile, also a fit to globular clusters, is also used to describe dwarf galaxies, more so in older literature. Using three parameters, a core radius $R_c$, a truncation radius $R_t$, and a characteristic density, $\Sigma_0$, the  King profile may be written as 
$$
\Sigma_{\rm K} = \Sigma_0\left(\frac{1}{\sqrt{1 + (R/R_c)^2}} - \frac{1}{\sqrt{1+(R_t/R_c)^2}}\right); \qquad \text{for } R<R_t
$$
and $\Sigma_{\rm K}=0$ for $R \geq R_t$. In much of the older literature, $R_t$ was interpreted as a "tidal radius," after an analogous interpretation for globular clusters [e.g., @hodge1961; @IH1995]. 

Finally, the @sersic1963 profile represents a generalization of an exponential profile and describes most dwarf galaxy light profiles well. Typically parameterized in terms of a half-light radius $R_h$, the density at half-light radius $\Sigma_h$, and a Sérsic index $n$, the profile's equation is
$$
\Sigma_{\rm S} = \Sigma_h \exp\left[-b_n \,  \left((R/R_h)^{1/n} - 1\right)\right]
$$
where $b_n$ the coefficient which solves $\Gamma(2n) = 2\gamma(2n, b_n)$ with $\Gamma$ the Gamma function and $\gamma$ the lower incomplete gamma function [@graham+driver2005]. A Sérsic profile with $n=1$ is equivalent to an exponential profile, while $n=4$ recovers \citepos{devaucouleurs1948} profile for elliptical galaxies. Although a Sérsic profile is less commonly applied to dwarf galaxies, @munoz+2018 advocate for the Sérsic profile since the added flexibility allows more profiles to be fit with a single law. 

While there are no clear theoretical preferences for any of these profiles, exponential density profiles have been commonly used for dwarf spheroidal galaxies. @faber+lin1983 were among the first to demonstrate that an exponential law is a reasonable empirical fit, theorizing that dwarf spheroidals may have evolved from exponential disk galaxies and maintained a similar light profile. Later, @read+gilmore2005 showed that exponential profiles may originate from mass loss during the evolution of dwarf galaxies. Tides are a possible mechanism for this transformation---the *tidal stirring hypothesis* [@mayer+2001a; @klimentowski+2009]. However, the theoretical origin of exponential disks is also unknown.[^exp_disk]

[^exp_disk]: While studied for far longer, the exponential origin of disk galaxies is no better understood. Ideas range from scattering of stars [@elmegreen+struck2013; @wu+2022] to angular momentum transport, to disk viscosity-driven radial gas flows [@lin+pringle1987; @wang+lilly2022] or spherical collapse [@freeman1970].

Many subsequent photometric studies of dwarf spheroidal galaxies have used exponential fits, finding that exponential and King profiles both provide good descriptions in many cases  [@binggeli+sandage+tarenghi1984; @mateo1998; @mcconnachie+irwin2006; @cicuendez+2018]. More recently, @moskowitz+walker2020 fit instead generalized Plummer profiles, but most of their fits would be consistent with a single-component exponential.[^moskowitz] As a result, it has become conventional to assume an exponential density profile to describe dwarf galaxies in theoretical or observational modelling [e.g., @kowalczyk+2013; @martin+2016; @MV2020a; @battaglia+2022]. 

[^moskowitz]: A single-component exponential profile is close to their "Steeper" profile in @moskowitz+walker2020 over the range of their data. However, most of their galaxies have little evidence to prefer the Steeper or Plummer density profile.

Dwarf galaxies outside the Local Group commonly follow exponential profiles, but sometimes with modifications. For example, many extragalactic dwarf elliptical, blue compact dwarf, and irregular dwarf galaxies are better described with an exponential profile to which a central cusp or nuclear region is added [@caldwell+bothun1987;  @noeske+2003]. On the other hand, some studies find an inner density decrement relative to exponentials [e.g., @caldwell+1992; @makarov+2012], or that dwarfs are better fit by two nested exponentials [e.g., @aparicio+1997; @graham+guzman2003; @hunter+elmegreen2006; @lee+2018]. It is unclear how these conclusions apply to the dwarf spheroidals of the Local Group.

Altogether, while there is some variation in the density profiles of dwarf galaxies, an exponential is an excellent first-order approximation. Typically, deviations from exponentials are in the direction of a steeper outer cutoff or changes to the inner slope of a dwarf galaxy (due, for example, to a nuclear star cluster). Flattened density profiles in the outer regions are more unusual. Explaining in detail the origin, similarity, and diversity of dwarf galaxy density profiles is an open question for theories of dwarf galaxy formation and evolution.

## The extended light profiles of Sculptor and Ursa Minor: Hints of tidal signatures? {#sec:scl_umi_obs_tides}

Sculptor (Scl) and Ursa Minor (UMi) appear to be typical dwarf spheroidal galaxies at first glance (see @fig:galaxy_images). [@tbl:scl_obs_props; @tbl:umi_obs_props] describe the structural parameters of each galaxy.  Sculptor, as the first discovered classical dSph, is even described as a "prototypical" dSph [e.g., @mcconnachie2012]. 

However, many studies have speculated that Scl and UMi have been influenced by the Milky Way's tidal field.  Already, @innanen+papp1979 found RR Lyrae candidate Scl members [from @vanagt1978] out to 180' in an elongated distribution, speculating this to be tidal disruption. Later density profile determinations noted Scl's elongation and apparent outer density excess [@eskridge1988; @IH1995; @walcher+2003;  @westfall+2006; but see also @coleman+dacosta+bland-hawthorn2005]. These studies often interpreted these features as evidence of tidal effects [e.g., @walcher+2003] or sometimes a "halo"of stars surrounding the dwarf [@westfall+2006]. 

UMi has attracted similar suspicions.  @hargreaves+1994 first detected a velocity gradient in UMi, suggestive of tidal disruption. Later, @martinez-delgado+2001 found stars far beyond the 
King-profile "tidal radius," aligned with UMi's elongation, consistent with tidally-stripping simulations in @gomez-flechoso+martinez-delgado2003. @palma+2003 furthermore detected S-shaped isophotes and an extended population of "extratidal" stars. 

Adding to the evidence, \citet{sestito+2023a, sestito+2023b} report a "kink" in the density profile of each galaxy.  They spectroscopically follow up distant stars, finding members as far as 6 and 12 half-light radii from the centre of each dwarf. If these dwarfs had exponential profiles, like Fornax, then these far-outlying stars should be much rarer. 

Sculptor and Ursa Minor are poorly described by an exponential profile. The left panel of @fig:scl_umi_vs_fornax shows the density profiles of Sculptor, Ursa Minor, and Fornax (see @sec:observations for details on how these profiles are measured). Compared to Fornax, both Sculptor and Ursa Minor show an excess of stars outside $\log R/R_h\approx 0.4$, implying densities which exceed 100 times the density of the exponential fit at large radii.

A goal of this work is to determine if tidal effects are indeed responsible for the extended outer light profiles of Sculptor and Ursa Minor. If tides cannot explain these features, these features may instead be due to an extended stellar "halo" or second component of the galaxy—suggestive of a complex star formation or assembly history.



![The extended stellar profiles of Sculptor and Ursa Minor](./figures/scl_umi_fornax_exp.pdf){#fig:scl_umi_vs_fornax}

Figure: Surface density profiles of Sculptor (orange squares), Ursa Minor (red triangles), and Fornax (green circles) scaled to their half-light radius and the density at half-light radius (data described in @sec:observations). The solid black line is an exponential profile ([@eq:exponential_law]). Scl and UMi show a clear excess over an exponential at large radii.



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

Table: Observed properties of Sculptor. References are: (1) @munoz+2018 [Sérsic fit], (2) @tran+2022 [RR lyrae distance], (3) @MV2020b, (4) @arroyo-polonio+2024. {#tbl:scl_obs_props  short="Observed properties of Sculptor"}



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

Table: Observed properties of Ursa Minor. References are: (1) @munoz+2018 [Sérsic fit], (2) @garofalo+2025 [RR lyrae distance], (3) @MV2020a, (4) @pace+2020, average of MMT and Keck results with systematic uncertainty from Appendix -@sec:extra_rv_models discussion.  {#tbl:umi_obs_props  short="Observed properties of Ursa Minor"}



# Interpreting tidal signatures {#sec:tidal_theory}

The Local Group hosts several examples of ongoing tidal disruption. The Magellanic stream, a massive, gas-rich feature emanating from the Magellanic clouds, is believed to arise partially from the MW's tides [@putman+1998; @diaz+bekki2012; @donghia+fox2016]. Other clear examples of tidal streams include the Sagittarius stream, the Andromeda Giant Southern stream, and the Tucana III stream [e.g., @ibata+gilmore+irwin1994; @ibata+2001;  @li+2018]. These examples illustrate that hierarchical accretion remains an active process. Interpreting such observations relies on simulations of tidal disruption. 

Cosmological simulations struggle to resolve tidal effects on dwarfs. Since many dwarfs are near the resolution limit, they are vulnerable to artificial disruption [e.g., @vandenbosch+2018; @santos-santos+2025]. To overcome numerical challenges,  idealized simulations model a single subhalo in an analytic host potential, achieving excellent numerical convergence. For example, the simulations we describe later reach three times higher resolution than Aquarius [@springel+2008] with 400 times fewer particles. Idealized simulations make numerous simplifications, neglecting mergers, cosmological context, mass assembly, and often baryonic physics [e.g., @hayashi+2003; @bullock+johnston2005; @klimentowski+2009; @ogiya+2019]. Cosmological simulations appear to predict that tidal disruption is more common than idealized simulations of MW dwarfs suggest, but the role of numerics and assumptions in this discrepancy are unclear [@panithanpaisal+2021; @shipp+2023; @riley+2024]. We shall use idealized simulations here to assess tidal effects after the satellite's infall into the MW halo. 

Idealized simulations predict clear properties of tidally disrupting dwarf spheroidal galaxies. Tidally stripped stars form *tidal streams*---stellar tails with a bulk velocity gradient [e.g., @moore+davis1994; @johnston+spergel+hernquist1995; @read+2006]. Most mass loss happens near pericentre, where tides are strongest. However, the central structure of a dwarf galaxy often remains undisturbed [@oh+lin+aarseth1995; @piatek+pryor1995]. For instance, NFW halos are also found to be resilient to full tidal disruption [@EP2020], but cored dark matter halos may disrupt fully and faster [e.g., @penarrubia+2010; @errani+2023a]. 

To first order, tidal mass loss peels away the outer layers of a dwarf galaxy in energy space.  \citet{drakos+taylor+benson2020, drakos+taylor+benson2022, amorisco2021} showed that tidal effects are nearly entirely described as the removal of particles above a truncation energy [see also @choi+weinberg+katz2009]. @stucker+2023 generalized this idea, creating a model for adiabatic tidal mass loss in an isotropic tidal field. Their model explains the resilience of NFW halos against full tidal disruption and the origin of well-defined "tidal tracks" [as observed in @PNM2008; @green+vandenbosch2019; @EN2021].

With precise orbital constraints and improved models of the Milky Way potential, recent studies have continued to probe the dynamical histories of individual dwarf galaxies. \citet{battaglia+sollima+nipoti2015, borukhovetskaya+2022, dicintio+2024} both ran simulations tuned to Fornax, showing that this galaxy's stellar component or globular clusters are likely not affected by tides. Similarly, @borukhovetskaya+2022a analyzed Crater II, showing that the present-day structure is challenging to reconcile with NFW initial conditions and Galactic tides. Most relevantly, @iorio+2019 also tailored simulations to Scl, finding weak Galactic tidal influence.

Building on this body of work, we will use idealized simulations to understand the severity of tidal effects on Sculptor and Ursa Minor.

## Tidal and "break" radii {#sec:break_radii}

For a given orbit in a given potential, there are characteristic radii which help gauge the effects of tides on a dwarf galaxy.

The **Jacobi radius** represents the approximate radius where stars become unbound for a galaxy in a circular orbit around a host galaxy.^[The Jacobi radius was derived at least as early as @laplace1798. This radius also bears other names, such as the Hill radius [from @hill1878]. Likely only named after Jacobi for the Jacobi integral [@jacobi1836].] Calculated from an approximation of the location of the $L_1$ and $L_2$ Lagrange points, the Jacobi radius is where the mean density of the dwarf galaxy is roughly three times the mean interior density of the host galaxy at pericentre, or
$$
3\bar \rho_{\rm MW}(r_{\rm peri}) \approx \bar \rho_{\rm dwarf}(r_J),
$$ {#eq:r_jacobi}
[@BT1987 eq. 7-84]. If $r_J$ is comparable to the visible extent of a galaxy, we should expect to find clear signs of tidal disturbance. While strictly valid for circular orbits, assuming $r_{\rm peri}$ for the host-dwarf distance works as most stars are lost during pericentric passages. 

We also use the **break radius** as defined in @penarrubia+2009, marking the outermost radius within which the dwarf has been able to achieve equilibrium after pericentric passage in a highly eccentric orbit. The break radius $r_{\rm break}$ is proportional to the velocity dispersion, $\sigma_v$, and time elapsed since pericentre, $\Delta t$, 
$$
r_{\rm break} = C\,\sigma_{v}\,\Delta t
$$ {#eq:r_break}
where the empirical constant is $C \approx 0.55$. $r_{\rm break}$ describes where the dynamical timescale becomes longer than the time since the perturbation, i.e., the radius within which the galaxy is dynamically relaxed. 

## A simple tidal simulation

To illustrate the effects of tides on an NFW halo due to a massive host, we consider a toy model. We evolve an N-body realization of an NFW subhalo with $\rmax=5\,\kpc$ and $\vmax = 27\,\kms$ orbiting a static NFW host halo with $\rmax=25\,\kpc$ and $\vmax = 207.4\,\kms$. These choices are motivated by the inferred structure and masses of dSphs and the Milky Way, respectively. We simulate the stellar component by assigning stellar weights to each particle based on the relative densities of the distribution functions in energy space (see @sec:painting_stars). The stars initially follow a 2D exponential profile with a scale radius $R_s=0.25\,\kpc$ and total mass $5\times10^5\,\Mo$, embedded within the inner dark matter halo. We start the model at apocentre on an orbit with a pericentre of $20\,\kpc$ and apocentre of $100\,\kpc$. See [@sec:methods] for a more detailed description of our numerical setup.

@fig:idealized_break_radius illustrates the properties of this idealized simulation at the second apocentre, as seen from the centre of the host. The projected density of stars is relatively undisturbed and spherical in the centre, but becomes non-isotropic outside the break radius and shows nascent tidal tails. These tidally disturbed stars appear as an extended, outer density excess relative to the initial conditions. This excess appears just outside the break radius. The break radius also marks where the mean 3D radial velocity of the stars (with respect to the galaxy's centre)  becomes non-zero — the galaxy is out of equilibrium outside $r_{\rm break}$. 

To first order, the final density profile of this toy simulation indeed resembles the outer excess in Sculptor and Ursa Minor, motivating this thesis.



![Example tidal simulation](figures/idealized_break_radius.pdf){#fig:idealized_break_radius}

Figure: Example density and velocity distributions of an idealized dwarf galaxy model undergoing tidal stripping. **Top**: The projected 2D stellar density in the $y'$--$z'$ plane for the initial (left) and final (right) simulation. $y'$ is the direction of tangential orbital motion, and $z'$ is the direction of orbital angular momentum, as measured from the centre of the host.  The dashed green circle represents the break radius (@eq:r_break) and the blue arrow points in the orbital direction.  **Bottom left**: The projected stellar density profile for the initial (dotted) and final (solid) simulation snapshot. **Bottom right**: the mean 3D radial velocity (the dot product of relative position and velocity relative to dwarf centre) as a function of projected radius. The green arrow marks the break radius in both lower panels.

# Thesis outline

The goal of this thesis is to review the evidence for an extended density profile in Ursa Minor and Sculptor, to assess the impact of tidal effects on each galaxy, and to discuss possible interpretations for the structure of these galaxies. 

In Chapter [-@sec:observations], we describe how we compute observational density profiles following @jensen+2024. In Chapter [-@sec:methods], we review our simulation methods. Next, we present our results for the tidal effects on Sculptor and Ursa Minor in Chapter [-@sec:results]. We discuss the implications of our results and end with a summary and outlook in Chapter [-@sec:discussion]. 



