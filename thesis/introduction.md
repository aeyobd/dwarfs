Dwarf galaxies host, in many ways, the most extreme galactic environments in the universe. These little galaxies are typically defined to be fainter than the LMC  [$M_V \gtrsim -18$ or similarly $L \lesssim 10^9 L_\odot$, e.g., @hodge1971; @mcconnachie2012]. As the smallest class of galaxies, dwarfs are the most cosmologically numerous galaxies. Dwarf galaxies are highly *dark-matter dominated*, with mass to light ratios exceeding 100 to 1000 $M_\odot/ L_\odot$.  Because of their small total masses, many dwarf galaxy satellites of the Milky Way are *quenched*, with little to no recent star formation. Milky Way dwarfs contain stellar populations which are *relics* from the early universe, consisting of many of the oldest and most metal poor stars. Understanding the properties of dwarf galaxies thus has implications across astronomy, from cosmological structure formation on the smallest scale to extremely metal-poor stellar populations to chemical evolution to galactic dynamics. Of the classical dwarfs, we find that the Sculptor and Ursa Minor dwarf spheroidal galaxies stand out with a more extended density profile relative to an exponential, hinting at tidal effects. As a case study in tidal effects and our understanding of dwarf galaxy formation and evolution, we aim to understand the dynamical history of these galaxies. 

In this section, we first describe cosmological structure formation and the role of dwarf galaxies. Then, we explore the observational history of dwarf galaxies, the *Gaia* telescope, open questions, and why Sculptor and Ursa Minor stand out. We  discuss idealized simulations, tidal effects, break radii, and recent developments in tidal simulations. Finally, we provide a roadmap to the thesis as a whole at the end of this section.

# Observational context

Dwarf galaxies have long raised conundrums in theories of galaxy formation. The discovery of Fornax and Sculptor in 1938 [@shapley1938][^lmc_discovery], with no known analogues, already presented such an enigma. Shapley presented these dwarfs as a new type of *stellar system* resembling the Magellanic Clouds and globular clusters but did not attempt to speculate on the exact nature. While dwarf galaxies were quickly understood to be galaxies based on the inferred luminosities and sizes, their nature and formation remained unclear for decades [e.g.,@hodge1971; @gallagher+wyse1994]. 

The earliest spectroscopic work hinted that dwarf galaxies may contain substantial amounts of dark matter. From early determinations of the velocity dispersion for the Sculptor and Ursa Minor dwarf spheroidal (dSph) galaxies  [e.g., @aaronson1983, @aaronson+olszewski1987],  inferred mass-to-light ratios were at least 10 times larger than the values from globular clusters scaled to the same sizes. While rather uncertain initially, these values have corroborated with larger and more precise samples [e.g., @hargreaves+1994]. Subsequently, several theories attempting to understand the formation and observed properties of these objects were proposed. Examples include: dwarf galaxies are undergoing tidal dissolution resulting in extreme mass-to-light ratios [e.g., @kuhn+miller1989], presence of massive central black holes [e.g., @strobel+lake1994], the formation of dark matter-free "tidal dwarfs" from past mergers [e.g., @lynden-bell1982, @kroupa1997], or modified theories of gravity [@milgrom1995].  However, consistency with CDM galaxy formation was also not out of the question [e.g., @dekel+silk1986]. Since then, we have known that dwarf galaxies are among the darkest objects in the universe, and understanding their properties is critical to understanding dark matter. 

Dwarf galaxies span a large range of physical sizes, luminosities, and morphologies. Broadly, there are three classes of dwarf galaxies based on luminosity, as illustrated by @fig:galaxy_images. Local **bright dwarfs** with magnitudes $-18 \lesssim M_V \lesssim  -14$, often exhibit irregular morphologies and recent star formation.  @fig:galaxy_images shows the Large Magellanic Cloud (LMC) as an example of an irregular, bright dwarf galaxy displaying a bar. **Classical dwarfs**  occupy intermediate luminosities ( $-14 \lesssim M_V  \lesssim -7.7$). Typically, these systems are old, gas-poor, and spheroidal. All Milky Way satellites discovered before digital sky surveys are classicals.  The 12 classical dwarfs satellites of our Galaxy are Sagittarius, Fornax, Leo I, Sculptor, Antlia II, Leo II, Carina, Draco, Ursa Minor, Canes Venatici I, Sextans I, and Crater II.[^dsph_suffix]   The **ultra-faint**s occupy the very faintest magnitudes ($-7.7 \lesssim M_V$). These galaxies have minuscule stellar masses, tend to be more compact, and are the darkest known galaxies. Altogether, dwarf galaxies span more than 15 orders in absolute magnitude.



Since Local Group dwarfs are nearby, we can study these galaxies on a star-by-star basis. As a result, it is possible to measure the 3D velocity and position of a star if we can measure the stars position, line-of-sight (LOS) velocity, and proper motions. Unfortunately, determining distances and full 3D velocities is challenging. The most direct measurement of distance, parallax, requires precise tracking of a star's sky position across a year. And while line-of-sight (LOS) velocities are easily determined from spectroscopy, tangential velocities, derived from proper motions and distances, are much more challenging. Typically, measuring proper motions requires accurate (much less than arcsecond) determinations of small changes in a star's position over baselines of several years to decades. The full 6D position and velocity information for stars has historically been known for only a handful of stars.



[^lmc_discovery]: Technically, the Large and Small Magellanic Clouds (LMC, SMC) may be classified as dwarf galaxies, but these were likely always known to humans at southern latitudes. 
[^dsph_suffix]: While formally the dwarf galaxy names we discuss contain "dwarf spheroidal" (dSph), e.g. Sculptor dSph, we omit this suffix for brevity.





![Dwarf Galaxy Pictures](figures/galaxy_pictures.png){#fig:galaxy_images width=390pt height=390pt}

Figure: Images of the LMC (DSS2), Fornax (DES DR2), Sculptor (DES DR2),  and Ursa Minor (UNWISE with Gaia point sources over-plotted). Each image includes the galaxy's luminosity and a 200 pc scale bar. *TODO: Mark half-light radius with circle?*





![Dwarf galaxies sky position](figures/mw_satellites_1.jpg){#fig:mw_satellite_system width=390pt}

Figure: The location of MW dwarf galaxies on the sky. We label the classical dwarf galaxies (green diamonds), fainter dwarfs (blue squares), globular clusters (orange circles), and ambiguous systems (pink open hexagons). Globular clusters are more centrally concentrated, but dwarf galaxies are preferentially found away from the MW disk. Sculptor and Ursa Minor are highlighted as two dwarfs we study later. The background image is from ESA/Gaia/DPAC (https://www.esa.int/ESA_Multimedia/Images/2018/04/Gaia_s_sky_in_colour2). Dwarf galaxies (confirmed), globular clusters, and ambiguous systems are from the @pace2024 catalogue (version 1.0.3). 

## The era of *Gaia* 

*(db: these paragraphs start the same...)*

*Gaia* has redefined astrometry, providing photometry, positions, proper motions, and parallaxes for over 1 billion stars [@gaiacollaboration+2021]. Launched in 2013, *Gaia* is a space-based, all-sky survey telescope situated at the Sun-Earth L2 Lagrange point [@gaiacollaboration+2016]. While *Gaia* completed its space-based mission in 2025, two forthcoming data releases remain. Determining absolute parallax measurement is facilitated by the observation that stars in different regions of the sky are affected by parallax motion with different phases. By imagining two regions separated by 106.5 degrees on the same focal plane, *Gaia* measures tiny changes in relative positions of stars across small and large angles. Combining measurements from multiple epochs across several years, an absolute all-sky reference frame is derived from which parallax and proper motions are derived. In addition to astrometry, *Gaia* measures photometry in the wide *G* band (330-1050nm) and colours from the blue photometer (BP, 330-680 nm) and red photometer (RP, 640-1050 nm). *Gaia* additionally provides low resolution BP-RP spectra and radial velocity measurements of bright stars (*Gaia* radial velocity magnitudes <16) [@gaiacollaboration+2016]. For this work, *Gaia*'s most relevant measurement are $G$ magnitude, $G_{\rm BP} - G_{\rm RP}$ colour, $(\alpha, \delta)$ position, and $(\mu_{\alpha*}, \mu_\delta)$ proper motions.[^pmra_cosdec]

*Gaia* has revolutionized many fields of astronomy including studies of the Milky Way and the Local Group. The Local Group is defined as galaxies which are approximately bound to the Milky Way-Andromeda system, or within about 1 Mpc from the Milky Way [e.g., @mcconnachie2012 and references therein]. For example, 6D kinematic classification of Milky Way stars led to the discovery of past mergers or Milky Way building blocks like *Gaia*-Sausage Enceladus [e.g., @helmi+2018], out-of-equilibrium structures like the *Gaia* snail [e.g., @antoja+2018, and dynamical effects of the Milky Way's spiral arms and the bar in the solar neighbourhood [@hunt+vasiliev2025]. Moving to the Milky Way halo, *Gaia* has helped find and constrain numerous stellar streams  [@bonaca+price-whelan2025]. Altogether, *Gaia* has revealed the hierarchical formation and complex, evolving structure of our own Galaxy.

For Milky Way satellites, *Gaia* has enabled well-constrained orbital analysis and facilitated precise stellar membership determination. Before, proper motions of dwarf galaxies were sparse or undetermined. Few galaxies had precisely measured proper motions, often from Hubble Space Telescope observations [e.g., @sohn+2017]. *Gaia* allowed for the first systematic determinations of Milky Way satellite proper motions of unprecedented precision [@MV2020a; @pace+li2019]. While the proper motion uncertainty on a typical dwarf member star is often large, by combining the proper motions for 100s to 1000s of stars in *Gaia*, precise proper motion measurements can be determined, sometimes only limited by *Gaia*'s systematic error floor  [e.g.,@MV2020a]. Proper motions have ushered in a new dynamical era for MW satellite studies, where we can derive precise orbits assuming a given MW potential. In addition, *Gaia* helps separate out contaminating MW foreground stars. By measuring parallaxes or proper motions, many background and foreground stars can be removed as non-members [e.g., @battaglia+2022, @jensen+2024].



[^pmra_cosdec]: The proper motions $\mu_\alpha$ and $\mu_\delta$ are the apparent rates of change in right ascension, $\alpha$, and declination, $\delta$, typically in units of mili-arcsecond (mas) per year. $\mu_{\alpha*} = \mu_\alpha \cos \delta$ corrects for projection effects in $\alpha$.



## Dwarf galaxies today

Today, we know that the Milky Way system is teeming with satellites. [@fig:mw_satellite_system] shows the MW satellite system, including dwarf galaxies, ambiguous systems (whose nature as a dwarf-like or cluster-like system remains uncertain), and globular clusters. Advances in observational techniques has accelerated progress in our understanding of dwarf galaxies and introduced new questions. Deep digital photometric surveys have more than quadrupled the number of known Milky Way dwarf galaxies [@simon2019]. Upcoming and ongoing surveys, like the Vera Rubin Observatory's Legacy Survey for Space and Time [@ivezic+2019], will likely identify even fainter dwarf galaxies. In addition, large aperture multi-object spectrographs have revealed the detailed and complex inner chemodynamical structure of dwarf galaxies. Beyond precise structural and kinematic properties of dwarf galaxies, modern observations allow for the separation of multiple stellar populations, detailed constraints on the dark matter density profiles, and hints of tidal disruption or stellar halos.

Of local dwarfs, the classical systems still remain the most studied galaxies with the best constrained parameters. Because many Milky Way dwarfs extend 1-2 degrees across the sky, measurements of some of their properties can be challenging [e.g., @mateo1998]. These dwarfs also contain large numbers of bright (giant) stars which can be followed up spectroscopically [e.g., @tolstoy+2023; @pace+2020]. As a result, the determination of fundamental properties such as the position, size, orientation, distances, proper motions, line-of-sight (LOS) velocities, and velocity dispersions $\sigma_v$, are all relatively well constrained today. However, ongoing research continues to redefine our understanding of the detailed structure of Milky Way satellites, hinting that these objects may be more complex than expected at first glance. 

Observations of dwarf galaxies have been the origin of several disputes or *small-scale* problems for $\Lambda$CDM [see review @bullock+boylan-kolchin2017]. For example, the mismatch between the expected number of dwarf galaxies from simulations and the observed abundance was known as the *missing satellites problem*. Additionally, a number of observations suggested that some dwarf galaxies, although not all, possess dark matter "cores" [e.g., @moore1994; @adams+2014; @oh+2015; @walker+penarrubia2011; @read+walker+steger2019].^[In detail, gas-phase rotation curves are better able to differentiate between cores and cusps, whereas stellar kinematics is less constraining.] As a result, alternative forms of dark matter have been advocated as solutions, such as Warm and Self-Interacting Dark Matter.

However, these tensions have eased as a result of improved understanding of baryonic physics. For example, recent hydrodynamic simulations in particular have shown that strong feedback can produce dark matter cores [e.g., @tollet+2016; @fitts+2017; @benitez-llambay+2019; @orkney+2021]. Altogether, the numerous past challenges for $\Lambda$CDM in the dwarf galaxy regime illustrates the opportunity for dwarf galaxies to test the understanding of galaxy formation and dark matter physics.



## Typical stellar density profiles

Projected density profiles efficiently characterize the shape of a galaxy. At its most basic, density profiles provide properties such as the shape, location, size, and orientation of a dwarf galaxy. In addition, the details of a stellar density profile can help interpret the total mass of dwarf galaxies, their assembly and dynamical history, and to understand scaling laws between dwarf galaxy properties [e.g., @herrmann+hunter+elmegreen2013, @penarrubia+2009, @querci+2025, @lee+2018]. We note that for resolved galaxies, these profiles are in terms of stellar count density instead of surface magnitudes.

Four different surface density laws are frequently used to parameterize dwarf galaxy profiles: Exponential, Plummer, King, or Sérsic profiles [e.g., @munoz+2018]. The exponential profile is perhaps the simplest, defined in terms of the central surface density, $\Sigma_0$, and scale radius, $R_s$:
$$
\Sigma_{\rm exp} = \Sigma_0\exp(-R / R_s).
$$

This profile is often applied to the radial light distribution of galaxy disks [@freeman1970; other classic paper?].

To fit globular cluster density profiles, @plummer1911 proposed a profile based on a polytropic solution^[where density and pressure are assumed to be related by a power law],
$$
\Sigma_{\rm Plummer} = \frac{\Sigma_0}{(1 + (R/R_h)^2)^2},
$$

where $\Sigma_0$ is the central surface density and $R_h$ is the 2D half-light radius. Now mostly superseded by the King profile for globular clusters, the Plummer model is still a good fit to many dwarf spheroidals [e.g., @moskowitz+walker2020].

The @king1962 profile, also an empirical fit to globular clusters, is also often applied to dwarf galaxies, especially in older literature. Using three parameters, a core radius $R_c$, a truncation radius $R_t$, and a characteristic density, $\Sigma_0$, the  King profile is 
$$
\Sigma_{\rm King} = \Sigma_0\left(\frac{1}{\sqrt{1 + (R/R_c)^2}} - \frac{1}{\sqrt{1+(R_t/R_c)^2}}\right).
$$

In much of the older literature, $R_t$ was often called and interpreted as a "tidal radius", after the similar interpretation for globular clusters [e.g., @IH1995, @hodge1961]. 

Finally, the @sersic1963 profile represents a generalization of an exponential profile, and describes most galaxy light profiles well. Typically parameterized in terms of a half-light radius $R_h$, the density at half-light radius $\Sigma_h$ and a Sérsic index $n$, the profile's equation is
$$
\Sigma_{\rm S\acute ersic} = \Sigma_h \exp\left[-b_n \,  \left((R/R_h)^{1/n} - 1\right)\right]
$$
where $b_n$ is a constant depending on $n$. $n=1$ provides an exponential profile and $n=4$ recovers @devaucouleurs1948's profile for elliptical galaxies. While a Sérsic profile is not always used for dwarf galaxies, @munoz+2018 advocate for the Sérsic profile since the added flexibility allows for more profiles to be fit. 

While there are no clear theoretical preferences for any of these profiles, exponential density profiles are commonly used for dwarf spheroidal galaxies. @faber+lin1983 were among the first to demonstrate that an exponential law is a reasonable empirical fit to dwarf galaxies, theorizing that dwarf spheroidal may have evolved from exponential disk galaxies and maintained a similar light profile. Later, @read+gilmore2005 showed that exponential profiles may originate from general mass loss during the evolution of dwarf galaxies. Many later photometric works for dwarf spheroidal galaxies have used exponential fits, finding that exponential and King profiles both provide good descriptions in many cases  [@binggeli+sandage+tarenghi1984, @mateo1998; @mcconnachie+irwin2006; @cicuendez+2018]. As a result, many studies assume an exponential density profile for dwarf galaxies in theoretical or observational modelling [e.g., @martin+2016 @MV2020a; @battaglia+2022; @kowalczyk+2013 but is for disk]. 

An exponential is far from the only density profile advocated for. Some authors fitting Sérsic or King profiles (often in addition to Exponential), find that the additional parameter somewhat improves fits relative to an exponential [@IH1995; @vanzee+barton+skillman2004; @munoz+2018; @wang+2019]. Note that the Sérsic indices are typically $n \lesssim 1$, implying that most galaxies tend to be an exponential or slightly steeper in the outer regions. 

Beyond Local Group dwarf spheroidals, exponential profiles appear to be common, but sometimes with significant modifications. In particular, many extragalactic dwarf elliptical and blue compact / Irr galaxies are fit well with an exponential + central cusp / nuclear region [@caldwell+bothun1987,  @noeske+2003]. In constrast, some studies find an inner decrement in density compared to exponentials [e.g., @caldwell+1992; @makarov+2012].

However, exponential profiles do not fit every dwarf galaxy. Even starting from @aparicio+1997,  and @graham+guzman2003  [for coma cluster dE], some people have noted deviations from this simple empirical rule. Additionally, @hunter+elmegreen2006; @herrmann+hunter+elmegreen2013; @herrmann+hunter+elmegreen2016; @lee+2018 note that at least in a photometric sample of more distant dwarf disky/Irr/blue compact dwarfs, that many or most dwarfs show two-part exponential profiles. Irr do represent a substantially different class though so it is unclear if these conclusions apply. @caldwell+1992's photometric study of M31 dwarfs show that the outer density profile typically fits well to an exponential but with inner deviations. Note that while diversity is also commented on here, dwarf galaxies with outer flattening profiles are often a minority in these studies. 

Alternatively, @moskowitz+walker2020, use a Plummer profile and modified Plummer profiles with steeper outer slopes ($\Gamma = d\log \Sigma / d \log R = -8$ compared to $\Gamma = -4$). For most galaxies, they do not have convincing statistical evidence to prefer one model over another, but show that for 8 / 15 of dwarfs with sufficient data, a plummer or steeper profile is preferred, with only the stepper one ever strongly preferred.  

Altogether, while there is some natural variation in the density profiles of dwarf galaxies, an exponential is an excellent first-order approximation. Typically, deviations from exponentials are in a direction of steeper outer cutoff or changes to the inner slope of a dwarf galaxy (due to a nuclear region/etc.). Flattened density profiles in the outer regions are more unusual. Explaining the origin, similarity, and diversity of dwarf galaxy density profiles is a pressing question for theories of dwarf galaxy formation and evolution. Whether dwarf galaxies are indeed exponential or not, the diversity of density profile shapes of dwarf galaxies needs to be explained. 

## Sculptor and Ursa Minor: Hints of tidal signatures?

Sculptor and Ursa Minor may appear to be typical dwarf galaxies at first glance. [@tbl:scl_obs_props; @tbl:umi_obs_props] describe the present-day properties of each galaxy. Sculptor, as the first discovered classical dSph, is often described as a "prototypical" dSph. There has long been speculation that both Sculptor and Ursa Minor may have been influenced by the Milky Way's tidal field (see @sec:discussion).  @sestito+2023a; @sestito+2023b report that a "kink" in the density profile, beginning around 30 arcmin for both Sculptor and Ursa Minor.  They spectroscopically follow up some of the most distant stars, finding multiple members between 6-12 half-light radii from the centre of each dwarf. If dwarfs had exponential profiles like Fornax, then these stars should be much rarer. 



Sculptor and Ursa Minor are not well-described by an exponential profile. The left panel of @fig:scl_umi_vs_penarrubia shows our density profiles of Sculptor, Ursa Minor, and Fornax (see @sec:observations for our methods). Compared to Fornax, both Sculptor and Ursa Minor show an excess of stars beginning around $\log R/R_h\approx 0.4$. The right panel of @fig:scl_umi_vs_penarrubia compares the same density profiles to an initial and final density profile from a simple tidal model of an exponential dwarf spheroidal galaxy embedded in a NFW halo in the Milky Way field, described in @sec:tidal_theory. Fornax agrees well with the exponential initial conditions. However, Sculptor and Ursa Minor are better described by the more extended, tidally evolved density profile.



A goal of this work is to determine, assuming $\Lambda$CMD, if the effects of the Milky Way (or other satellites) may indeed create observable tidal signatures in Sculptor and Ursa Minor. If tides cannot explain these features, these features may instead be a extended stellar "halo" or second component of the galaxy—suggestive of a complex history or formation in each galaxy, which needs to be explained by galaxy formation models of Local Group dwarf galaxies.



![Sculptor and Ursa Minor match tidal models](./figures/scl_umi_vs_penarrubia.pdf){#fig:scl_umi_vs_penarrubia}

Figure: A plot of the surface density profiles of Sculptor, Ursa Minor, and Fornax scaled to their half-light radius and the density at half-light radius (data described in @sec:observations). The right panel adds an idealized model with initial, final, and break radius as a dotted line, solid line and an arrow respectively [see @sec:break_radii].



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





# Cosmological context

We only understand a tiny fraction of the universe's composition. The leading theory of cosmology, $\Lambda$CDM (Lambda Cold Dark Matter), posits that the universe is composed of about 68% dark energy ($\Lambda$), 27% dark matter (DM), and 5% regular baryons[^baryons]  [@planckcollaboration+2020]. While the composition of dark matter and dark energy remains elusive, we know their general properties. Dark energy causes the acceleration of the expansion of the universe on large scales. We do not discuss dark energy here—it does not substantially affect the Local Group today. Dark matter, instead, makes up the vast majority of mass in galaxies. Typically, galaxies have baryonic to dark matter ratios of between 1:5 to beyond 1:1000 for faint dwarf galaxies. In $\Lambda$CMD, dark matter is assumed to interact only gravitationally. Light passes through dark matter unimpeded---in this sense dark matter is transparent. Dark matter is also commonly assumed to be *cold*, i.e. typical velocities much smaller than the speed of light in the early universe. Implications of dark matter properties range from cosmological structural formation, galaxy structure, and galaxy interactions. 

[^baryons]: Astronomers like to change definitions of words. *Baryons* here means baryons+leptons, i.e. any standard model massive fermion.

## Structure formation and dwarf galaxies

The very early universe was almost featureless. Our earliest observations of the universe stem from the cosmic microwave background (CMB)---revealing a uniform, isotropic, near-perfect blackbody emission. But tiny perturbations in the CMB, temperature fluctuations of 1 part in 10,000, reveal the underlying seeds of large-scale cosmological structure.^[The power spectrum from the CMB is complicated by baryons, which lead to oscillations in power based on the sound-crossing timescale.] Governed by cosmological expansion, gravitational collapse, and baryonic physics, each overdensity grows into larger structures. Initially, baryonic matter was coupled to radiation and resisted collapse. Dark matter, only influenced by gravity instead, freely collapsed into the first structures. For mass perturbations both smaller than the horizon and overdense enough to collapse, these overdensities of dark matter become self-gravitating, known as *halos*. After recombination, where electrons combined with atomic nuclei to form atoms, baryons decoupled from radiation and fell into the dark matter halos. The densest pockets of baryons later formed the first stars and galaxies.

Dark matter halos, and their associated galaxies, rarely evolve in isolation. Instead, structure formation is *hierarchical*. Small dark matter halos collapse first and hierarchically merge into progressively larger halos [e.g.,@blumenthal+1984; @white+rees1978; @white+frenk1991]. Structure formation happens on a wide range of scales. The largest structures become the cosmic web—composed of voids, filaments, and clusters—and the smallest structures directly observable host dwarf galaxies. Hierarchical assembly is evident through the large scale structure of the universe, remnants of past mergers within the Milky Way, and tidal disruption of dwarf galaxies and their streams around nearby galaxies.

Small-scale structure formation is sensitive to deviations from $\Lambda$CDM cosmology [e.g., @bechtol+2022]. One key prediction of $\Lambda$CDM is that mass perturbations are expected to exist on all scales, and are largest on the smallest scales, so we would expect the formation of halos on all scales. Many alternative models, such as warm dark matter or self-interacting dark matter, may smooth out small-scale features and reduce the abundance of small halos or change their structure [e.g., @lovell+2014]. Dwarf galaxies, which occupy the smallest dark matter halos, are promising windows into small-scale cosmological features. Alternative dark matter properties can influence dwarf galaxy abundance, formation, structure, and tidal evolution. The nearby dwarf galaxies of the Local Group, with detailed observations, present a promising opportunity to understand the evolution of these objects and test our understanding of cosmology. To understand the evolution of these objects, we need to understand the general predictions of $\Lambda$CDM for the properties of dwarf galaxies. 



*Warm* dark matter, relativistic in the early universe but cooler now, smooths out small-scale features and softens the cusps of dark matter halos. *Self-interacting* dark matter instead can form cores but also the cores can collapse into a density peak. **TODO: incorporate this here.**

## Typical structure of dark matter halos

In $\Lambda$CDM cosmological simulations, dark matter halos are remarkably self-similar. In \citet{NFW1996, NFW1997}, hereafter NFW, the authors observe that the spherically-averaged density profiles $\rho(r)$ are universally well described by a two-parameter law: 
$$
\rho/\rho_s= \frac{1}{(r/r_s)(1+r/r_s)^2},
$$ {#eq:nfw}
where $r_s$ is a scale radius and  $\rho_s$ a scale density . This profile has shown remarkable success in describing $\Lambda$CMD halos across several orders of magnitude in mass. NFW profiles are *cuspy*, where the density continuously rises like $\rho \sim 1/r$ at small radii $r \ll r_s$. The steepness of the density profile increases gradually with radius, and at large radii the density  falls off like $\rho \sim 1/r^3$. @fig:nfw_density shows an example of a NFW halo. 

The total mass of an NFW profile diverges, so halos are commonly defined using an overdensity criterion. The virial mass, $M_{200}$, is defined as the mass within a radius $r_{200}$ containing a mean enclosed density 200 times[^200] the critical density of the universe:
$$
M_{200} =200\,\frac{4\pi}{3} \ r_{200}^3\ \rho_{\rm crit}, \qquad {\rm where} \quad \rho_{\rm crit}(z) = 3H(z)^2 / 8\pi G
$$
A standard second parameter is the halo concentration, $c=r_{200} / r_s$ describing how the characteristic size of the halo compares to its virial radius.  In this case, the scale density is a function of $c$ alone, $\rho_s = (200/3)\,\rho_{\rm crit} c^3 / [\log(1+c) - c/(1+c)]$ [@NFW1996]. However, we characterize halos by their circular velocity profiles. The circular velocity, $v_{\rm circ}(r) = \sqrt{G M(r) / r}$, reaches a maximum of $v_{\rm max}$ at radius  $r_{\rm max} \approx 2.16258\,r_s$. $v_{\rm max}$ is related to both the total halo mass and the observed line of sight velocity dispersion, and $r_{\rm max}$ relates to the scale radius. 

While an NFW profile has two free parameters ($M_{200}$ and $c$), these are not independent. Lower mass dark matter halos often collapse earlier, when the universe was denser. As a result, low mass subhalos tend to be more concentrated [e.g.,@NFW1997]. The relationship between $M_{200}$ and c, or the mass-concentration relation, describes the mean trend of concentration with mass [e.g., @bullock+2001; @ludlow+2016]. The left panel of @fig:smhm illustrates the present-day mass-concentration relation in terms of $r_{\rm max}$ and $v_{\rm max}$. While there is a general trend for concentration to increase with mass, the concentration values still have substantial scatter. Other parameters such as the halo spin or shape may affect the scatter of the mass-concentration relation, but their effect is typically expected to be small [@navarro+2010; @dicintio+2013; @dutton+maccio2014]. 



[^200]: For the collapse of a uniform spherical density, the virialized overdensity  would be $\Delta = 18\pi^2\approx 178$ for a critical universe $\Omega_m = 1$. This is commonly rounded to $\Delta = 200$. While this parameter may be closer to $\Delta \approx 100$ for our universe, $\Delta$ also increases with redshift [using eq. 6 from @bryan+norman1998]. 

![Example density profiles](figures/example_density_profiles.pdf){#fig:nfw_density}

Figure: Density profiles in log density versus log radius for the stars and dark matter of a Sculptor-like galaxy. The dark matter is more extended and massive than the star across the entire galaxy---dynamical evolution is driven by the dark matter mass distribution. 

![Stellar-mass halo-mass relation](figures/cosmological_means.pdf){#fig:smhm}

Figure: **Left** Radius of maximum circular velocity $r_{\rm max}$ as a function of maximum circular velocity. The solid line with 1-$\sigma$ shaded region and dashed line are the relations from @ludlow+2016 for $z=0$ and $z=2$ respectively. **Right** Stellar mass (top) as a function of maximum circular velocity. The solid line with the 1-$\sigma$ shaded region is the relation from @fattahi+2018 with scatter points simulated central galaxies from APOSTLE in @fattahi+2018. 

## Connecting dark matter to stars

Unfortunately, directly observing dark matter remains infeasible (lest dark matter wouldn't be dark). Fortunately, dwarf spheroidal galaxies come with a convenient population of visible stars. Stars typically represent a small fraction of the total mass of a dwarf galaxy. For example, @fig:nfw_density shows the typical density profile for the stars and dark matter for a Sculptor-like galaxy. At all radii, dark matter are at least ~10 times more dense, and the total virial mass is often ~1000 times higher than the stellar mass alone. As a result, the dynamics and formation of dwarf galaxies is governed by the underlying dark matter halo. 

Both cosmology and observations find fundamental relationships between the dynamical or total mass and stellar masses of galaxies.**Improve last sentence.** A prediction of cosmological simulations is the Stellar mass Halo Mass relation (SMHM), describing the mass of stars which form in a halo of a given size. In particular, the SMHM relation grows especially steep in the dwarf galaxy regime---many dwarf galaxies are formed in halos of similar mass. @Fig:smhm shows the stellar mass versus $v_{\rm max}$ maximum circular velocity (proxy for halo mass) for APOSTLE dwarfs from @fattahi+2018. The APOSTLE project [@sawala+2016] uses the hydrodynamical setup from the EAGLE project [@crain+2015; @schaye+2015] to simulate Local Group analogues in a $\Lambda$CDM cosmological context (db: reorder these sentences). While there is some scatter, the range of predicted $v_{\rm max}$ is fairly narrow across $\sim 4$ decades in stellar mass. Because lower mass galaxies have increasingly shallow potential wells, feedback becomes more effective at removing gas. Re-ionization additionally suppresses late star formation in the faintest galaxies. As a result, the resulting stellar mass of a galaxy is highly sensitive to the dark matter mass, especially for faint dwarf galaxies. If we know the initial stellar mass of a dwarf galaxy, we also know, at some level, the properties of its dark matter halo.

Several challenges complicate a simple SMHM trend including environment, assembly history, and tidal effects. Additionally, the detailed physics of galaxy formation can impact the SMHM relationship. By being closer to a massive host, many dwarfs quench earlier, resulting in lower stellar masses [e.g., @christensen+2024]. In particular, effects like ram-pressure stripping (removal of gas in the dwarf galaxy due to pressure from host's circumgalactic medium) and tidal removal of gas cause star formation to quench [REFS]. Additionally, the time of formation (relative to reionization) can influence the resulting stellar content [@kim+2024]. Finally, tides influence both the dark matter and stellar mass but in different amounts [e.g., @PNM2008]. Consequently, tides may reduce the halo mass more than the stellar mass, adding additional scatter to the SMHM trend, particularly for low-mass satellites [e.g., @fattahi+2018]. Understanding the effects of tides on Local Group dwarf galaxies is essential to determine where and how these galaxies formed in a cosmological context. 







# Simulating tidal effects {#sec:tidal_theory}

Since the discovery of dwarf galaxies, a large body of work has either considered or simulated tidal effects. While pre-*Gaia* work often had few constraints on orbits of dwarf galaxies, the theory of tidal mass loss remains largely the same. 

Simulating dwarf galaxies accurately in a cosmological context remains a substantial challenge. Currently, cosmological simulations can predict the overall abundance of larger mass? dwarf galaxies [e.g., ] and broadly predict effects of tides [e.g., @riley+2024]. But, dwarf galaxies are often near the resolution limit. Insufficient resolution can lead to artificial disruption of dwarf galaxies and overestimation of tidal effects [e.g., @santos-santos+2025]. For example, the highest resolution cosmological simulation of a Milky Way-size dark matter halo, the Aquarius project [@springel+2008], achieved a DM resolution of $1.712 \times10^3 \Mo$, enough to barely resolve the stellar components of Sculptor like halos (see chapter [@-sec:methods]).  To address this challenge,  idealized simulations only simulate a single subhalo in an approximate host potential. As a result, idealized simulations can reach excellent numerical convergence at the cost of realism. For instance, our simulations later are about 3x higher resolution than Aquarius at a fraction of the computational cost (400x less particles). Idealized simulations make numerous simplifications: neglecting mergers, cosmological context and evolution, mass assembly, and often baryonic physics. We use idealized simulations here which are more powerful in accurately assessing tidal effects after infall, given that the idealizations do not impact these galaxies's recent history too much. 

Some of the earliest theory work on tidal mass loss of dwarf galaxies originate from @oh+lin+aarseth1995; @piatek+pryor1995; @moore+davis1994; @johnston+spergel+hernquist1995. Already, these works used techniques similar to what we continue to use today and laid the foundation for understanding tidal effects. While the detailed assumptions have evolved (e.g. no longer assuming mass-follows-light), many of their conclusions still hold true. More recent work expanding this theory include @read+2006; @bullock+johnston2005; @PNM2008; @penarrubia+2009; @klimentowski+2009; @errani+2023a; @fattahi+2018; @stucker+2023; @wang+2017.

With precise orbits and a better understanding of the Milky Way potential and system, more recent work began to directly probe the dynamical histories of individual dwarf galaxies (although early work began this for  Sagittarius / etc). Examples include @iorio+2019 for Sculptor, @borukhovetskaya+2022; @dicintio+2024 for Fornax; @borukhovetskaya+2022a for Antlia II. Our goal is to apply a similar framework to Sculptor and Ursa Minor.

***These past three paragraphs need revisions...***

## Simulating large gravitational systems: The N-body method

Modelling gravitational evolution for large systems requires special methods. Perhaps the simplest method to compute the evolution of dark matter is through *N-body simulations*. A dark matter halo is represented as a large number of dark matter particles (bodies). Each body is essentially a Monte Carlo sample of the underlying phase-space distribution. Note that dark matter (and galaxies) are often assumed to be *collisionless*—particles are not strongly affected by close, *collisional* gravitational encounters which substantially change the momenta of involved bodies. In contrast, star clusters are often collisional so neglecting these encounters may not be a reasonable approximation. While we use individual gravitating bodies in N-body simulations, the Newtonian gravitational force is softened to be a Plummer sphere as to limit strongly collisional encounters.

Naively, the Newtonian gravitational force requires adding together the forces from each particle on each particle, with a computational cost that scales quadratically with the number of particles, or $O(N^2)$. With this method, simulating a large number of particles, such as $10^6$, would require $10^{12}$ force evaluations at each time step, making cosmological and high-resolution studies unfeasible. However, only long-range gravitational interactions tend to be important for CMD, so we can utilize the *tree method* to compute the gravitational force vastly more efficiently.

The first gravitational tree code was introduced in @barnes+hut1986, and is still in use today. We utilize the massively parallel code *Gadget 4* [@gadget4].  Particles are spatially split into an *octotree*. The tree construction stars with one large node, a box containing all of the particles. If there is more than one particle in a box/node , the box is then divided into 8 more nodes (halving the side length in each dimension) and this step is repeated until each node only contains 1 particle. With this heirarchichal organization, if a particle is sufficiently far away from a node, then the force is well approximated by the force from the centre of mass of the node. As such, each force calculation only requires a walk through the tree, only descending farther into the tree as necessary to retain accuracy. The total force calculations reduce from $O(N^2)$ to $O(N\,\log N)$, representing orders of magnitude speedup. Modern codes such as *Gadget* utilize other performance tricks, such as splitting particles across many supercomputer nodes, efficient memory storage, adaptive time stepping, and parallel file writing to retain fast performance for large scale simulations, forming the foundation for many cosmological simulation codes. 



## Tidal evolution

A galaxy in equilibrium will remain in equilibrium, unless acted upon by an external force. As an example, @fig:lagrange_points illustrates the effective potential around a Sculptor-like NFW galaxy in a circular orbit around a Milky Way-like galaxy. The effective potential accounts for the centrifugal force in the rotating frame. There are two key saddle points in @fig:lagrange_points, labeled $L_1$ and $L_2$ (as Lagrange points, there are another 3 not as relevant). Since a particles trajectory is constrained to be within the equipotential surface equal to the particle's total energy, the easiest (requiring the lowest energy) path for a particle to escape the satellites gravitational influence is through the windows around $L_1$ and $L_2$. Otherwise, the potential steeply increases. $L_1$ and $L_2$ are co-linear with the satellite and host-origin. In general, particles which have higher energy than the energy at $L_1$ and $L_2$ are most commonly unbound from the galaxy. 

The full evolution of a dwarf galaxy in a tidal follows an approximate sequence (although somewhat concurrent as well.)

1. *Mass is lost*. In particular, particles and stars on weakly bound orbits are most likely to be removed by tides. Particles escape through $L_1$ and $L_2$
2. *Steams form*. Unbound mass becomes part of a stream or tidal tails. These particles follow similar orbits but are slightly higher or lower energy (depending on which side the particle escaped from).
3. *Mass redistributes*. Bound mass of the galaxy redistributes. As mentioned in the last section, this is visible as a wave of outward moving material, with the outermost material reaching equilibrium last. 
4. *A new equilibrium*. With mass loss, the gravitational potential decreases, resulting in a more compact dark matter halo and stars which adiabadically expand to a larger scale radius.

A tidally affected galaxy contains predictable observational clues to this process. Nearby to the galaxy, the density profile is flattened as a result of newly unbound material. Further away from the galaxy, this would be visible as a stream. Additionally, a satellite stream contains a velocity gradient along the path. Finally, due to mass loss, a satellite will have a reduced velocity dispersion, larger stellar size, and 

Dwarf galaxies evolve along *tidal tracks*. From N-body simulations across a numerous initial conditions and orbits, NFW halos evolve along a narrow track in terms of maximum circular velocity and radius $v_{\rm max}$ and $r_{\rm max}$ . @EN2021 derive an empirical fit to these tidal tracks, finding that 

$$
\frac{v_{\rm max}}{v_{\rm max, 0}} = 
2^\alpha 
\left(\frac{r_{\rm max}}{r_{\rm max, 0}}\right)^{\beta}\left[1 + \left(r_{\rm max} / r_{\rm max, 0}\right)^2\right]^{-\beta},
$$
where $\alpha=0.4$ and $\beta=0.65$. As illustrated in @fig:tidal_tracks, this formula works for both circular and elliptical orbits and is independent of the initial subhalo size or distance to the host.  ***Rewrite this paragrah***

![Lagrange points](./figures/lagrange_points.pdf){#fig:lagrange_points}

Figure: The potential contours for a 2-body galaxy system in the rotating frame. $L_1$ and  and $L_2$ is most efficient as the required energy (the effective potential) is the smallest. 



![Tidal tracks of dwarf galaxies](/Users/daniel/Library/Application Support/typora-user-images/image-20250715095615423.png){#fig:tidal_tracks width=80%}

Figure: Tidal tracks of dwarf galaxies, the logarithm of maximum circular velocity and radius of relative to the initial conditions for satellites on a variety of orbits. Almost all satellites follow the tidal track suggested by @EN2021. Adapted from fig. 6 of @EN2021.

## Break and tidal radii {#sec:break_radii}

Following the detailed discussion of tidal evolution, we can interpret two simple analytic approximations of where tidal effects should be visible.

The **Jacobi radius** represents the approximate radius where stars become unbound for a galaxy in a circular orbit around a host galaxy. Calculated from an approximation of the location of the $L_1$ and $L_2$ Lagrange points, the Jacobi radius is where the mean density of the dwarf galaxy is roughly three times the mean interior density of the host galaxy at pericentre, or
$$
3\bar \rho_{\rm MW}(r_{\rm peri}) \approx \bar \rho_{\rm dwarf}(r_J).
$$
If $r_J$ occurs within the visible extent of a galaxy, we should expect to find relatively clear signs of tidal disturbance. While strictly valid for circular orbits, assuming $r_{\rm peri}$ for the host-dwarf distance works as most stars are lost near pericentre. 

We also use the **break radius** as defined in @penarrubia+2009, marking where the galaxy is still in disequilibrium. The break radius $r_{\rm break}$ is proportional to the velocity dispersion ,$\sigma_v$, and time elapsed since pericentre ,$\Delta t$, 
$$
r_{\rm break} = C\,\sigma_{v}\,\Delta t
$$
where the scaling constant $C \approx 0.55$ is derived empirically. $r_{\rm break}$ describes where the dynamical timescale is longer than the time since the perturbation, i.e. the radius within which the galaxy should have dynamically relaxed.  As illustrated in @fig:idealized_break_radius, for an idealized model with exponential stars in a NFW halo, shortly after pericentre, the stellar component is smooth but contains a change in slope around $r_{\rm break}$. This radius is visible in the stellar distribution as non-spherical S-shaped overdensities of stars. Also, $r_{\rm break}$ is where where the mean radial velocities of the stars becomes positive, i.e. the system is out of equilibrium and adjusting to the new density profile.





![Break radius validation](figures/idealized_break_radius.pdf){#fig:idealized_break_radius}

Figure: Example density and velocity distributions of an idealized simulation shortly after pericentre. **Top left**: The 2D density profile for the initial and final simulation with the break radius marked.  The break radius of the simulations is set by the time since pericentre.  **Top right**: The projected 2D stellar density in the $x$-$y$ plane. The green circle represents the break radius and the grey arrow points towards the host centre. **Bottom left**: the mean radial velocity (dot product of relative position and velocity relative to dwarf centre) as a function of 2D radius. **Bottom right**: The mean radial velocity as projected into 2D bins. The initial conditions are Sculptor-like, exponential stars embedded in NFW, evolved in @EP2020 potential with a pericentre of 10 kpc and apocentre of 100 kpc shortly after the 3rd pericentre. See @sec:methods for a description of our simulation setup. *Note to self:* This figure may be too complex, not sure the best way to present everything yet.

# Thesis outline

In this thesis, our goal is to review the evidence for an extended density profile in Ursa Minor and Sculptor, to assess the impact of tidal effects on each galaxy, and to discuss possible interpretations for the structure of these galaxies. 

In chapter [@-sec:observations], we describe how we compute observational density profiles from @jensen+2024. In chapter [@-sec:methods], we review our simulation methods. Next, we present our results for the tidal effects on Sculptor and Ursa Minor in chapter [@-sec:results], We discuss our results, limitations, and implications in chapter [@-sec:discussion]. Finally, chapter [@-sec:summary] summarizes this work and discuss future directions for similar work and the field of dwarf galaxies. 



