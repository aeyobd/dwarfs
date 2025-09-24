

# Reliability and comparison of observed density profiles

## Caveats

While the J+24 method is excellent and has been verified, there are several possible limitations we discuss here. In general, while these limitations are likely real, they do not substantially affect our conclusions up to where we suggest the density profiles are unreliable.

**Spatial likelihood.** J+24 method was designed in particular to detect the presence of a density excess and find individual stars at large radii to be followed up. We are more interested in accurately quantifying the density profile and size of any perturbations. One potential problem with using J+24's candidate members is that the algorithm assumes the density is either described by a single or double exponential. If this model does not accurately match the actual density profile of the dwarf galaxy, we want to understand the impact of this assumption.In particular, in @Fig:scl_observed_profiles, notice that the $P_{\rm sat}$ selection method produces small errorbars, even when the density is more than 1 dex below the local background. These stars are likely selecting stars from the statistical MW background consistent with UMi PM / CMD, recovering the assumed density profile. As a result, the reliability of these density profiles below the CMD+PM background may be questionable. A more robust analysis, removing this particular density assumption, would be required to more appropriately represent the knowledge of the density profile as the background begins to dominate. 

**Absolute cuts in Bayesian inference and uncertainty misrepresentation.** A more self consistent model would fit the density profile to the entire field at once, eliminating possible misrepresentation of the uncertainties. In the Appendix to this section, we *will* also discuss an alternative method which runs a MCMC model using the likelihoods above to solve for the density in each elliptical bin. From these tests, we note that the density profile and uncertainties derived from the J+24 sample are reliable insofar as the dwarf density is above the background from MW CMD+PM-consistent interlopers. We estimate that this effect comprimises the density profiles past $\log R / {\rm arcmin} = 1.8$ for Sculptor and Ursa Minor, but agreement is good before then. 

***Gaia* systematics**. While *Gaia* has shown excellent performance, some notable limitations may introduce problems in our interpretation and reliability of density profiles. Gaia systematics in proper motions and parallaxes are typically smaller than the values for sources of magnitudes $G\in[18,20]$. Since we use proper motions and parallaxes as general consistency with the dwarf, and factor in systematic uncertainties in each case, these effects should not be too significant. However, the systematic proper motion uncertainties becomes the dominant source of uncertainty in the derived systemic proper motions of each galaxy (see [@tbl:scl_obs_props; @tbl:umi_obs_props]).

**Completeness**. *Gaia* shows high but imperfect completeness, particularly showing limitations in crowded fields and for faint sources ($G\gtrapprox20$). As discussed in @fabricius+2021, for the high stellar densities in globular clusters, the completeness relative to HST varies significantly with the stellar density. However, the typical stellar densities of dwarf galaxies are much lower, at about 20 stars/arcmin = 90,000 stars / degree, lower than the lowest globular cluster densities and safely below the crowding limit of 750,000 objects/degree for BP/RP photometry. In @fabricius+2021, for the lowest density globular clusters, the completeness down to $G\approx 20$ is $\sim 80\%$. Closely separated stars pose problems for Gaia's on-board processing, as the pixel size is 59x177 mas on the sky. This results in a reduction of stars separated by less than 1.5" and especially for stars separated by less than 0.6 arcseconds. The astrometric parameters of closely separated stars furthermore tends to be of lower quality [@fabricius+2021]. However, even for the denser field of Fornax, only about 3% of stars have a neighbour within 2 arc seconds, so multiplicity should not affect completeness too much (except for unresolved binaries). One potential issue is that the previous analyses do not account for our cuts on quality and number of astrometric parameters. These could worsen completeness, particularly since the BP-RP spectra are more sensitive to dense fields. In the appendix to this section, we test if magnitude cuts impact the resulting density profiles, finding that this is likely not an issue. 

**Structural uncertainties**. J+24 do not account for structural uncertainties in dwarfs for the two component case. We assume constant ellipticity and position angle. Dwarf galaxies, in reality, are not necessarily smooth or have constant ellipticity. J+24 test an alternative method using circular radii for the extended density component, and we find these density profiles are very similar to the fully elliptical case, even when assuming circular bins for the circular outer component. As such, even reducing the assumed ellipticity from $0.37-0.55$ to 0 does not substantially impact the density profiles. 

**Comparison to literature density profiles.** In Appendix XX, we compare our derived Sculptor and Ursa Minor profiles to various previous determinations. We find good agreement, and an excess in the density profiles has long been acknowledged (as discussed below.)

**Outlook.** Finally, an excellent test of systematics and methods in *Gaia* is to compare against another survey. In the appendix, we show that profiles derived from J+24 agree with DELVE or UNIONS data within uncertainties. Systematics, completion effects, and selection methods are unlikely to substantially change the density profiles presented here. 



# Limitations of simulations

Variations to the potential of the inner disk (exclusion of a bar, spiral arms) should minimally affect our results as no orbit we consider reaches less than ~15 kpc of the MW centre. We exclude the mass evolution of the halo from this analysis. Over $10\,$Gyr, this would be fairly significant (factor of $\sim 2$in MW mass, REF) but since we want to determine the upper limit of tidal effects, it is safe to neglect this. 



## Challenges to backwards-time integration {#sec:scl_orbit_variability}

The long term orbits of satellites are highly uncertain. Limitations of the assumed potential such as neglecting triaxiality, halo twisting, the location of subhalos / halo lumpiness, and the evolution of the potential make point orbits in a static potential a poor approximation [@dsouza+bell2022]. However, mass growth of the Milky Way implies that satellites would be less bound in the past. But, energy and angular momentum are not necessarily conserved, and pericentric distances similarly overestimated [@santistevan+2023]. 

As an illustration of the challenges, even in an idealized potential, of backwards time integration, @fig:scl_orbit_lmc_uncert shows orbits of Scl in different magellanic cloud potentials. Before a lookback time of 4 Gyr, Scl may be at almost any distance between 4 and 300 kpc from the MW at a given time. Indeed, Scl may have even experienced a more extreme earlier MW pericentre. 

However, a single extreme pericentre may be insufficient to create the observed tidal signatures. As an example, we take an orbit with a MW pericentre of just 4 kpc. The Jacobi radius should create observed tidal effects. However, because of the extreme radial nature of this orbit (also reaching an apocentre of 300 kpc), the pericentric passage causes more stellar mass loss but the observed stellar component would still appear smooth. 

Ursa Minor's long term orbital history is an interesting  case. A possibility is that UMi was previously a satellite of the LMC in the two passage LMC scenario [@vasiliev2024], and may have even been "tidally-processed" by the MW. But, once again, if UMi was stripped in the first LMC pericentre, than its orbit cannot be too tightly bound. We do briefly investigate the impact of an extreme LMC pericentre with UMi, but this model still runs into similar issues as the MW pericentric passage---UMi is unlikely to experience this tidal field enough to change the inner stellar distribution. 

Finally, we note that the bias in long-term time integration is towards the underestimation of the pericentre. Effects like the mass growth of the MW and dynamical friction cause long-term orbits to shrink in peri and apocentric distances. Modeling a frictionless MW system would 

Unfortunantly, even at their $2\sigma$ deviation of an increase of a factor of 2 in previous tidal forces [@santistevan+2023], this would only bring the Jacobi radii down to $ \gtrsim 2\,\kpc$. While a possibility, an extrodinary miscalculation of the orbital history of Scl and UMi having multiple pericentres with $4\times $ the tidal force seems an unlikely explanation. Furthermore, the time of last pericentre is recent enough orbital analysis should not deviate substantially enough to predict a consistent break radius. 







![Sculptor Orbits with LMC](/Users/daniel/thesis/figures/scl_lmc_orbits_mass_effect.png){#fig:scl_orbit_lmc_uncert}

Figure: The long-term orbital history of Sculptor (**top**) and UMi (**bottom**) are uncertain. In both panels,  light, transparent lines represent randomly-sampled  orbit of the satellites (alla ref) in three different LMC/MW mass models from @vasiliev2024. The LMC orbits are in solid, thick lines of the corresponding colour. The L2M11 has a lighter LMC mass, and the L3M10 model has a lighter MW mass than our fiducial L3M11 LMC model. **todo: add UMi**

Our LMC model only simulates recent tidal effects---the long-term orbit of Scl is essentially unknown. @fig:scl_lmc_orbits_effect illustrates the orbits of Scl for three different LMC-MW models in @vasiliev2024: L2M11 has a lighter LMC and L3M10 has a lighter MW than the L3M11 potential we used. All orbits diverge significantly after 2 Gyr ago. In addition, other affects such as the lumpiness of the MW halo, the growth of the MW halo, dynamical friction, and limitations of this potential make long-term backwards time integration yet more unreliable [e.g., @dsouza+bell2022]. 

In @fig:scl_lmc_orbits_effect, some orbits have a small $\lesssim 10\,\kpc$ pericentre with the Milky Way. In the Appendix, we consider a model with a previous MW pericentre of 4kpc. While this encounter should affect the inner stellar distribution of Sculptor based on the Jacobi radius, the observed light profile remains unaffected.  Because this model experiences just one tidal shock, stronger tidal streams form, but the tidal impact is not strong enough to reshape the observed stellar distribution. We conclude that 



As discussed in @vasiliev2024, if the LMC has had a past, second pericentre with the Milky Way, there is a chance that UMi was previously bound to the LMC. We find that this occurs about 15% of the time for this choice of potential. However, backwards integration may not appropriately represent the true probability of a Magellanic origin, since tidal stripping is not time reversible. Because of the possibility of this past encounter, and the high uncertainties associated with backwards integration, we note this as a possibility but defer to the discussion for our interpretations.

# Comparison and review of previous work

For both Sculptor and Ursa minor, these galaxies have been studied extensively in both a theoretical and observational context. While many works considered mass modelling, star formation histories, chemistry, and other observational properties, we focus on the observations and models concerning tides and dynamical evolution here.



## Sculptor

A number of works have previously speculated that dwarf galaxies such as Sculptor and Ursa Minor may have signatures of tidal perturbations, however not without contention. While @hodge1961 and @demers+krautter+kunkel1980 were some of the earliest work to derive the density profile of Sculptor, @innanen+papp1979 was perhaps the first to speculate that Sculptor may harbour a substantial population of "extratidal stars" (stars beyond the tidal radius), finding candidate members out to 180' in an elongated distribution from @vanagt1978's catalogue of variable stars.[^uranus] Additionally, @eskridge1988 showed a possible excess of stars in Sculptor relative to a @king1962 and Exponential density profile beginning around 50', but suggest that this excess may not be unusual. Later work by @IH1995, @walcher+2003, and @westfall+2006 also showed evidence of an extended component of Sculptor's density profile (among other dwarfs), interpreting these stars as either evidence of  tidal debris or a dwarf galaxy halo. On the other hand, @coleman+dacosta+bland-hawthorn2005 show that Sculptor is well described by a two component density profile, they additionally mention that it is unlikely many of Sculptor's stars (less than 2%) are extratidal. While @munoz+2018 is the most recent and deepest photometric study of Sculptor, they only cover an area of out to radii of about 30', so they are unable to study the possible extended components of Sculptor detected in past works. In addition to the presence of "extratidal stars", many studies note that Sculptor's stars appear to become more elliptical with radius, consistent with tidal effects (@IH1995; @westfall+2006).

[^uranus]:  @innanen+papp1979 also speculate about Uranus's satellite distribution, covering gravitational dynamics at scales of 100,000km and scales of 10,000,000,000,000,000 km.



As mentioned earlier, a number of works have found two, distinct, chemodynamical populations in Sculptor [@tolstoy+2004;  @battaglia2008], revealed to be also age gradients [@deboer+2011], in the distribution functions[@breddels+helmi2014], chemistry [@kirby+2009] and later followed up with larger samples [@tolstoy+2023; @arroyo-polonio+2024].

While motivated by @sestito+2023a, we arrive at a different conclusion as to the tidal nature of the excess. 

@iorio+2019 is the most similar work to ours on Sculptor. They find a similar story, showing that Scl may be affected by large DM mass loss but the stellar component is likely unchanged. Our models and orbits are in agreement with theirs, only our updated observables provide tighter constraints. However, they do not consider the LMC, which we show is a critical ingredient in the tidal history of Scl. Otherwise, few works have modelled tidal effects in detail. For example, @tchiorniy+genina2025 consider tidal effects on the inner density and mass estimation of dwarf galaxies, but find limited effects and only consider mean orbits in static potentials without an LMC.

## Ursa Minor

Many works have claimed the detection of substructure and signatures of tidal disruption in UMi. Starting with the first density profile from @hodge1964, Ursa Minor was shown to be highly elliptical---not unexpected for tidal effects. Studies by @olszewski+aaronson1985, @demers+1995, @IH1995, @kleyna+1998, @battinelli+demers1999, and @bellazzini+2002  note that Ursa Minor appears to contain substructure along its major axis, but without strong interpretation for the causes. In particular, there appears to be a kinematically cold clump [e.g., @pace+2014]. If a star cluster, than this only survives for several Gyrs provided the UMi's halo is not cored [@kleyna+2003] and without too much substructure/halos [@iora+2012]. @wilkinson+2004 interpret a drop in the velocity dispersion in the outskirts as an additional kinematically-cold popultion. 

@hargreaves+1994 first detected a velocity gradient in UMi similar to those expected for tidal disruption. One of the first direct suggestions of tidal features came from @martinez-delgado+2001, who find that stars extend far beyond the nominal tidal radius for Ursa Minor in the direction of the galaxies elongation, interpreting this as evidence for tidal effects, and with corresponding simulations by @gomez-flechoso+martinez-delgado2003. @palma+2003 further added evidence for Ursa Minor's peculiar morphology including S-shaped contours, a possible extended extratidal population of stars, and a failure for a King density profile to adequately capture the actual stellar distribution. Each of these characteristics appear to indicate strong tidal disruption. 

However, @munoz+2018's more modern photometry show a more regular distribution of stars, similar to our results. 



## Tides in General

Our conclusions about weak tidal effects are not alone: a number of works have claimed either common or rare tidal effects on local dwarf galaxies. 

Some works have theorized that tides are key to forming dwarf spheroidals as we know it [@tsujimoto+shigeyama2002; @mayer+2001a; ], @lynden-bell1982a and a magellanic origin for all dwarf galaxies. 

However, some cosmological simulations show a large number of disrupting satellites. For example, @riley+2024 find that most satellites in Auriga appear to be stripping (losing more than 3% stellar mass, similar to our UMi extreme orbit with a plummer initial condition). 



@penarrubia+2009 [also defining @eq:r_break] use their break radius to show that, based on orbits then known, that the break radii were not consistent with earlier interpretations of tides in Scl and UMi. In addition, @pace+erkal+li2022 use the criterion that the Jacobi radius should be larger than the half-light radius (based on the observed half-light mean density), finding that the only local dwarfs, with or without the LMC, flagging a number of likely disrupting dwarfs but excluding Scl and UMi. @read+2006 furthermore consider the signatures of tides on the velocity dispersion, showing that a lack of a rising, outer velocity dispersion profile is an indication that dwarfs are unlikely to be disrupting (which is not seen in many dwarfs).

Tides are also known to generally affect the maximum circular velocity. Our assumptions and dynamical effects are similar to the @robles+bullock2021 comparison of dark matter-only simulations in a MW potential.



# Understanding the extended density profiles of dwarf galaxies

A number of recent work has also confirmed the presence of extended, likely non-tidal, density profiles. For example, @chiti+2021; @chiti+2023 spectroscopically confirm members out to $~9 R_h$ in Tucana II. 

While some galaxies are speculated to be undergoing tidal effects and thus show extended profiles [Boo I, @roderick+2016], etc. However, (recent) tides cannot always be the explanation, as this thesis shows.x

Besides extended density profiles, many dwarf galaxies also show evidence for multiple chemodynamical stellar populations. Typically, the older population is more extended, kinematically hotter, and lower metallicity, whereas the younger population is more compact, colder, and higher metallicity. Examples include Sculptor [@tolstoy+2004, @arroyo-polonio+2024], Ursa Minor [@pace+2020], Carina [@battaglia+2012, @fabrizio+2016, @kordopatis+2016], Fornax [@battaglia+2006, @amorisco+evans2012, @delpino+aparicio+hidalgo2015], Sextans [@battaglia+2011; @cicuendez+battaglia2018; @roederer+2023 ], and Andromeda II [@mcconnachie+arimoto+irwin2007; @ho+2012; @delpino+2017]. Dwarf spheroidal galaxies are likely more complex than a single component stellar population.

We discuss scenarios ranging from more mundane to extraordinary (in our opinion) below. We review the plausability of each scenario later.

### Dynamical heating of old stars

For a variety of reasons, older stars in a galaxy may be preferentially "hotter" (i.e., higher typical random velocities) than younger stars. 

In the case of dwarf galaxies, a few mechanisms have been proposed which may heat the stellar components: heating by strong feedback, heating by dark sub-subhalos, and heating 

Because dwarf galaxies have extended star formation histories, like the stellar-feedback scattering method, 

@el-badry+2016. Feedback drives fluctuations in size + radial migration => population gradients.

A number of other explanations below involve dynamical heating in some aspect. 





### Multi-epoch star formation

Many simulations predict that a population gradient should form naturally in dwarf spheroidals. In some cases, this may appear as distinct chemodynamical populations [@kawata+2006]. In this continuous star formation case, halos may additionally be supplied with stars migrating out to the extended halo through initial supernovae velocity kicks or dynamical heating of preferentially old stars [@stinson+2009; @revaz+jablonka2018].

The star formation history of dwarf galaxies is typically thought to be *bursty*: a series of discrete episodes of intense star formation separated by periods of quiescence [e.g., @salvadori+ferrara+schneider2008; @valcke+derijcke+dejonghe2008;  @wheeler+2019]. Quiessence may also be driven by reionization causing a temporary pause in star formation [@benitez-llambay+2015].

Star formation may also be reignited by varios environmental encounters. For example, tidal compression [@mayer+2001a], encounters with gaseous streams [surry nature paper], shock with the MW halo [@wang+2024a?].

Besides the creation of multiple stellar populations, feedback can also drive potential fluctuations which scatter stars into the halo [@maxwell+2012].



- @wright+2019

- @azartash-namin+2024

- @dong+lin+murray2003

- @bermejo-climent+2018

- @mercado+2021: formation of metallicity gradients from SFH 

  

One challenge for these models is why Sculptor and Ursa Minor would form a different stellar component than other dwarf spheroidals. 

@starkenburg+helmi+sales2016: disk disruption from the accretion of a dark stellar halo, driving morphological changes and a starburst. 

### Stellar accretion

@ricotti+polisensky+cleland2022, others...

### Past major mergers 

When galaxies merge, traces of the effects of the merger can remain for long afterwards. For dwarf galaxies, mergers may appear as multiple populations which may have different sizes, chemistry, and kinematics. 

A number of studies has investigated mergers in dwarf galaxies, showing that this is a feasible scenario for halo formation. Mergers can be categorized according to if they have gas, but the dynamics are similar in either case. Gas only provides fuel for post-merger star formation of a new population of stars. 

@stierwalt+2015: observational evidence for dwarf-dwarf interactions.

A few galaxies have been suspected to have formed from a merger. @lokas+2014, @fouquet+2017, @amorisco+evans+vandeven2014: And II from merger, 

@tarumi+yoshida+frebel2021: Shows cosmological simulations in context of Tuc II, merger of two galaxies in the early universe produces a more extended density profile and metallicity gradient.



- @benitez-llambay+2016: mergers formting bimodal stellar populations through heating old stars and forming new population. Gradients form naturally here but mergers excasperate them. 
- @deason+2014
- @deason+2022: Mergers may be common in cosmological simulation. In particular, intermediate mass mergers of 1:5 may create the most visible stellar halo. However, their predicted densities of the halos are fainter than ours
- @querci+2025
- rotation from merger: @cardona-barrero+2021
- @genina+2019: \apostle{} galaxies: Mergers with gas can heat old population and form new population subsequently. Additionally, ram-pressure stripping and interactions with cosmic filimants can remove outer gas and trigger inner star formation. 
- @stierwalt+2015: tiny titans.

### Preprocessing

Another alternative is that the dwarf galaxies have been "preprocessed" by another dwarf/LMC like galaxy earlier on and this has resulted in their unique properties [e.g., @santistevan+2023; @riley+2024]. From the orbit integrations above, it is possible that UMi was a satellite of the LMC and preprocessed by this system. However, it is unclear how many halos in the Milky Way would be more massive than Scl / UMi and have had a chance to interact with these galaxies. 

### Tidal dwarfs

On e possibility is that Sculptor and Ursa Minor are not formed from dark matter halos but instead due to the collision of gas. This possibility is not a convincing explanation because we know the star formations histories are extended and relatively old, and the galaxies are close to apocentre, making achieving the observed velocity dispersion challenging in this framework. A key approximation here is that the dwarf galaxy should be almost completely radially anisotropic, leading to a radial velocity gradient. [@kroupa1997]

 We do not detect any strong gradient in velocity or dispersion with radius in Appendix -@rv_models. 

### Alternatives to $\Lambda$CDM

Alternatives to the prevailing cosmology we assume may help explain the nature of these dwarf galaxies or complicate our interpretations. 

SIDM significantly complicates tidal evolution. SIDM halos in isolation are not static---by e-transfering heat through collisions, these halos first undergoe core formation and then core collapse. Besides structural changes, SIDM adds that pressure from the host DM halo can change the structure (e.g. aiding core collapse) or remove material from the inner galaxy (analogous to ram pressure stripping). While we have shown a core alone is unlikely to change our conclusions, there may be a scenario where tidal effects through this SIDM-evaporation may help heat or induce mass loss.  Also, note that SIDM has many of its own challenges in explaining the properties of dwarf galaxies [e.g., @zavala+vogelsberger+walker2013], so 

Fuzzy dark matter can also heat up stars owing to density perturbations in the form of interference fringes. Similar to other intrinsic heating methods, this model would likely affect all dwarfs similarly, so we should detect similar halos around similar dwarfs. 

MOND: @sanchez-salcedo+hernandez2007: mond in dsph. 



## Disentangling the origin of a stellar halo.



In this thesis, we have investigated tidal effects on the dwarf galaxies Sculptor and Ursa Minor. 

We first showed, using data from J+24, that each galaxy contains a substantial, extended population when compared to an exponential density profile. We show that this population is reasonably independent of the details of the data selection or methodology, implying that this "density excess" is likely a real feature of each galaxy. Sculptor and Ursa Minor are more extended in the outer regions than other nearby dwarf classical, given current observations.

We then investigated if tides were a permissable explanation. By modelling each galaxy based on cosmological initial conditions, we showed that tides do not strongly affect either galaxy. However, the presence of an LMC in the MW possibly complicates our conclusions, but we find that the effect of the LMC means Scl could be on first infal and has experienced a significantlly weaker tidal history, and Ursa Minor may only have a more extended orbit. In all of our investigation, we show that tides are very unlikely to affect the currently-observed density profiles. A stream is possibly detectable, depending on the initial stellar extent of each galaxy and how much deeper future observations are able to reach.

As a result, we conclude that the observed extended density of each galaxy is inconsistent with a tidal origin, and requires a different explanation. In this (final) chapter, we have summarized possible other explanations in the literature,  including the effects of mergrs, episodic star formation histories, stellar halos, 

While we have discussed a number of different avenues through which stellar halos form, these need not be exclusive.

As an example, consider +frenk+carlos+++: when two gas-rich dwarf galaxies collide, the merger can introduce new star formation, possibly forming a new, distint population of stars. The formation of stars may drive a burst which heats up old stars. The merger may scatter old stars into a stellar halo. And the merger may distribute the two populations of the old galaxies into separate components. As a result, many different processes can contribute to the formation of chemically and dynamically distinct populations. 

Is there hope for determining the nature of a stellar halo? Perhaps the most promising tool will be chemistry. In particular, if we reach a point where we can chemically distinguish the populations of different galaxies, then chemistry provides an opportunity to test the merger history of galaxies. 

We propose a categorization of the (\LCDM{}) scenarios above 

- Gas-poor mergers or accretion: more than one *chemodynamically independent* population (e.g. GSE in the MW halo)
- Gas-rich merger: more than one *chemodynamically independent* population plus an additional, younger population with possibly shared heritage.
- Continuous dynamical heating (e.g. subsubhalos) : chemically contiguous population with an old to young population gradient. 
- Strong stellar feedback: Associated bursty star formation history in chemical abundance space, but reasonably chemically-contiguous population. 

In addition, the first two scenarios could show significant variation among dwarfs (especially for (near) major-mergers), whereas the last two should be relatively similar among similar dwarfs. 

Key observations will be more precise and larger samples of radial velocities and proper motions (e.g. HST-Promo collaboration, or successors to Gaia, thirty-metre class telescopes, @evslin2016) and detailed chemistry for many stars. Combining these informations to reveal a detailed, chemodynamical picture of local dwarfs will provide some of the most favorable opportunities to understand the detailed origin and assembly of our Galactic neightbours, and by extension, our own galaxy.
