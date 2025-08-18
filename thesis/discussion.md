

# Reliability and comparison of observed density profiles

## Caveats

While the J+24 method is excellent and has been verified, there are several possible limitations we discuss here. In general, while these limitations are likely real, they do not substantially affect our conclusions up to where we suggest the density profiles are unreliable.

**Spatial likelihood.** J+24 method was designed in particular to detect the presence of a density excess and find individual stars at large radii to be followed up. We are more interested in accurately quantifying the density profile and size of any perturbations. One potential problem with using J+24's candidate members is that the algorithm assumes the density is either described by a single or double exponential. If this model does not accurately match the actual density profile of the dwarf galaxy, we want to understand the impact of this assumption.In particular, in @Fig:scl_observed_profiles, notice that the $P_{\rm sat}$ selection method produces small errorbars, even when the density is more than 1 dex below the local background. These stars are likely selecting stars from the statistical MW background consistent with UMi PM / CMD, recovering the assumed density profile. As a result, the reliability of these density profiles below the CMD+PM background may be questionable. A more robust analysis, removing this particular density assumption, would be required to more appropriately represent the knowledge of the density profile as the background begins to dominate. 

**Uncertainty misrepresentation.** A more self consistent model would fit the density profile to the entire field at once, eliminating possible misrepresentation of the uncertainties. In the Appendix to this section, we *will* also discuss an alternative method which runs a MCMC model using the likelihoods above to solve for the density in each elliptical bin. From these tests, we note that the density profile and uncertainties derived from the J+24 sample are reliable insofar as the dwarf density is above the background from MW CMD+PM-consistent interlopers. We estimate that this effect comprimises the density profiles past $\log R / {\rm arcmin} = 1.8$ for Sculptor and Ursa Minor, but agreement is good before then. 

***Gaia* systematics**. While *Gaia* has shown excellent performance, some notable limitations may introduce problems in our interpretation and reliability of density profiles. Gaia systematics in proper motions and parallaxes are typically smaller than the values for sources of magnitudes $G\in[18,20]$. Since we use proper motions and parallaxes as general consistency with the dwarf, and factor in systematic uncertainties in each case, these effects should not be too significant. However, the systematic proper motion uncertainties becomes the dominant source of uncertainty in the derived systemic proper motions of each galaxy (see [@tbl:scl_obs_props; @tbl:umi_obs_props]).

**Completeness**. *Gaia* shows high but imperfect completeness, particularly showing limitations in crowded fields and for faint sources ($G\gtrapprox20$). As discussed in @fabricius+2021, for the high stellar densities in globular clusters, the completeness relative to HST varies significantly with the stellar density. However, the typical stellar densities of dwarf galaxies are much lower, at about 20 stars/arcmin = 90,000 stars / degree, lower than the lowest globular cluster densities and safely below the crowding limit of 750,000 objects/degree for BP/RP photometry. In @fabricius+2021, for the lowest density globular clusters, the completeness down to $G\approx 20$ is $\sim 80\%$. Closely separated stars pose problems for Gaia's on-board processing, as the pixel size is 59x177 mas on the sky. This results in a reduction of stars separated by less than 1.5" and especially for stars separated by less than 0.6 arcseconds. The astrometric parameters of closely separated stars furthermore tends to be of lower quality [@fabricius+2021]. However, even for the denser field of Fornax, only about 3% of stars have a neighbour within 2 arc seconds, so multiplicity should not affect completeness too much (except for unresolved binaries). One potential issue is that the previous analyses do not account for our cuts on quality and number of astrometric parameters. These could worsen completeness, particularly since the BP-RP spectra are more sensitive to dense fields. In the appendix to this section, we test if magnitude cuts impact the resulting density profiles, finding that this is likely not an issue. 

**Structural uncertainties**. J+24 do not account for structural uncertainties in dwarfs for the two component case. We assume constant ellipticity and position angle. Dwarf galaxies, in reality, are not necessarily smooth or have constant ellipticity. J+24 test an alternative method using circular radii for the extended density component, and we find these density profiles are very similar to the fully elliptical case, even when assuming circular bins for the circular outer component. As such, even reducing the assumed ellipticity from $0.37-0.55$ to 0 does not substantially impact the density profiles. 

**Outlook.** Finally, an excellent test of systematics and methods in *Gaia* is to compare against another survey. In the appendix, we show that profiles derived from J+24 agree with DELVE or UNIONS data within uncertainties. Systematics, completion effects, and selection methods are unlikely to substantially change the density profiles presented here. 

## Comparison to literature



A number of works have previously speculated that dwarf galaxies such as Sculptor and Ursa Minor may have signatures of tidal perturbations, however not without contention. While @hodge1961 and @demers+krautter+kunkel1980 were some of the earliest work to derive the density profile of Sculptor, @innanen+papp1979 was perhaps the first to speculate that Sculptor may harbour a substantial population of "extratidal stars" (stars beyond the tidal radius), finding candidate members out to 180' in an elongated distribution from @vanagt1978's catalogue of variable stars.[^uranus] Additionally, @eskridge1988 showed a possible excess of stars in Sculptor relative to a @king1962 and Exponential density profile beginning around 50', but suggest that this excess may not be unusual. Later work by @IH1995, @walcher+2003, and @westfall+2006 also showed evidence of an extended component of Sculptor's density profile (among other dwarfs), interpreting these stars as either evidence of  tidal debris or a dwarf galaxy halo. On the other hand, @coleman+dacosta+bland-hawthorn2005 show that Sculptor is well described by a two component density profile, they additionally mention that it is unlikely many of Sculptor's stars (less than 2%) are extratidal. While @munoz+2018 is the most recent and deepest photometric study of Sculptor, they only cover an area of out to radii of about 30', so they are unable to study the possible extended components of Sculptor detected in past works. In addition to the presence of "extratidal stars", many studies note that Sculptor's stars appear to become more elliptical with radius, consistent with tidal effects (@IH1995; @westfall+2006).

[^uranus]:  Interestingly, @innanen+papp1979 also speculate about Uranus's satellite distribution, covering gravitational dynamics at two very different scales.

Ursa Minor was often noted not for a density excess but for unusual features. Starting with the first density profile from @hodge1964, Ursa Minor is known to be highly elliptical. Studies by @olszewski+aaronson1985, @demers+1995, @IH1995, @kleyna+1998, and @bellazzini+2002  note that Ursa Minor appears to contain substructure along its major axis, but without strong interpretation for the causes. As one interpretation @kleyna+2003 suggest that this feature is a long-lived star cluster residing in a cored dark matter halo to be dynamically stable. One of the first direct suggestions of tidal features came from @martinez-delgado+2001, who find that stars extend far beyond the nominal tidal radius for Ursa Minor in the direction of the galaxies elongation, interpreting this as evidence for tidal effects. @palma+2003 corroborate many of these earlier findings, showing further evidence for Ursa Minor's peculiar morphology including S-shaped contours, a possible extended extratidal population of stars, and a failure for a King density profile to adequately capture the actual stellar distribution. Each of these characteristics appear to indicate strong tidal disruption. Using velocity observations, @pace+2014 additionally show that there are multiple components in spatial-velocity information. However, @munoz+2018's more modern photometry show a more regular distribution of stars. Irregardless, Ursa Minor has had enough work suggesting peculiarities that a deeper investigation into the possibility of tidal effects is worthwhile. 



However, everything discussed so far is before *Gaia*. Thus, the knowledge of the orbits of these systems is largely unknown and most surveys could only filter out stars using photometry. As a result, the tidal radius could only be guessed based on either the density profile or current distance of the dwarf galaxy. With *Gaia* (discussed more below), recent work has used Bayesian frameworks to derive systematic proper motions of many dwarf galaxies (e.g. @MV2020a) and filter away foreground stars using proper motions and parallax to detect distant members, and study the 6D internal kinematics of these galaxies (e.g. @tolstoy+2023). Most relevant here, @jensen+2024 present a bayesian algorithm to determine likely members in *Gaia* (described below), and detect extended secondary components for 9 dwarf galaxies, including Sculptor and Ursa Minor. @sestito+2023a and @sestito+2023b followed up a few of the most distant stars detected in @jensen+2024 confirming that members exist in each galaxy as far ass 9-12 $R_h$ from the centres. In this chapter, we discuss the origin and reliability of these detections, confirming that the density profiles are indeed robust. These density profiles then provide the scientific motivation to explore the tidal interpretation in the next chapters. 



# Limitations of simulations

Variations to the potential of the inner disk (exclusion of a bar, spiral arms) should minimally affect our results as no orbit we consider reaches less than ~15 kpc of the MW centre. We exclude the mass evolution of the halo from this analysis. Over $10\,$Gyr, this would be fairly significant (factor of $\sim 2$in MW mass, REF) but since we want to determine the upper limit of tidal effects, it is safe to neglect this. 



## Long term orbital history

@dsouza+bell2022, @santistevan+2024. challenges to backwards time integration

## Deviations from the NfW

- @dicintio+2013

# Comparison to other work

## Sculptor

Theoretical work on Sculptor

- @battaglia+2008
- @iorio+2019
- @amorisco+zavala+deboer2014
- @battaglia+2008
- @breddels+2013
- @breddels+helmi2013
- @richardson+fairbairn2014
- @SFW2017
- @innanen+papp1979
- @wilkinson+2002
- @yang+2025: chemical evolution in Scl.
- @skuladottir+2024; @skuladottir+2021; @lee+jeon+bromm2024, pop III high res spectro?
- @wang+2024a hydrosim of dwarf galaxies like sculptor
- @tang+2023: apogee modeling of Scl.
- @pandey+west2022 chem evo OMEGA and isotopes.
- @an+koposov2022: distance gradients?
- @kawata+2006, 2pop origins?
- @grcevich+putman2009, @carignan+1998 ambiguous HI clouds near Scl ?
- @agnello+evans2012, claim that Sculptor cannot hold two populations in NFW halo?

Observational work on Scl

- @sestito+2023a
- @westfall+2006 wide degree survey for extended structure.
- @tolstoy+2023, @arroyo-polonio+2023, @arroyo-polonio+2024
- @eskridge1988, @eskridge1988a, @eskridge1988b
- @coleman+dacosta+bland-hawthorn2005
- @DQ1994
- @WMO2009
- @IH1995, 
- @munoz+2018: Using Megacam to derive density profiles and structural properties of many dwarf spheroidal galaxies.
- @kirby+2009
- @martinez-vazquez+2015, @pietrzynski+2008
- @grebel1996
- @barbosa+2025: Using DECam to derive narrowband photometric metallicity gradient and search for metal poor stars in Scl. 

Future ideas:

- @evslin2016: measuring ellipticities of halos w TMT.



## Ursa Minor

- @sestito+2023b
- @pace+2020
- @pace+2014 substructure
- @bellazzini+2002
- @hargreaves+1994
- @martinez-delgado+2001
- @munoz+2005: velocity dispersion profile
- @palma+2003
- @spencer+2018
- @vitral+2023
- @piatek+2005: old HST Pm meas.
- @gallagher+2003: no ionizaed gas.
- @shetrone+cote+stetson2001: RGBs
- @wilkinson+2004
- @kleyna+2003: dynamical fossil.
- @pryor+kormendy1990, @lake1990 MLR / velocity dispersion
- @armandroff+olszewski+pryor1995; @olszewski+pryor+armandroff1996 MLR
- @aaronson+olszewski+hodge1983; @aaronson1983 carbon stars and velocities
- @stetson1984 early spectroscopy.
- @tsujimoto+shigeyama2002 chem evo.



Theoretical work

- robles+bullock2021
- @caproni+lanfranchi2021, @caproni+2015: simulation of gass loss.
- @bajkova+bobylev2017, local orbits
- @gomez-flechoso+martinez-delgado2003: MLR of UMi
- @lynden-bell1976 early discussion of orbits and LMC plane.

# Are Sculptor and Ursa Minor typical?

## The formation of Exponential profiles

As mentioned earlier, many dwarf galaxies appear to have exponential stellar density profiles. 

Disk galaxies have long been known to be exponential-like across a wide range of scales. As such, several theoretical works have aimed to undertsand the formation of exponential stellar disks. One prevailing theory posits that scattering of stars in a disk naturally forms an exponential [@elmegreen+struck2013, @wu+2022], or that the disk is due to conservation of angular momentum during spherical collapse. 

Given that dwarf galaxies were noticed to be exponential like above, @faber+lin1983 proposed that dwarf spheroidals formed from disky galaxies maintaining the typical exponential density profile

@read+gilmore2005 propsed that exponential dSph form through impusive mass loss, redistributing the stellar component into an exponential like profile. 

- disk galaxy formation @fall+efstathiou1980, @mestel1963
- valcke+derijcke+dejonghe2008 with early idealized hydro simulations show Sérsic values from simulations between 0.8 and 1 for dSph, indicating a natural formation of Exponential or slightly steeper density profiles. 

One possible explanation is that @mayer+2001a

Other explanations range from the effects of mass loss, feedback, angular momentum,  

@klimentowski+2007, @klimentowski+2009. 

In summary, while there is some theoretical reasons why exponential spheroidal stellar profiles may form, the emperical sucess of an exponential law still requires explanations. Deviations from this exponential thus test wether this exponential law is indeed near-universal or if other ingredients in the formation and history of dwarf galaxies may be responsible for changes from the empirical expectation.

## Alternatives to exponential profiles

- elmegreen+hunter2006: Discuss the formation of double exponentials in disk galaxies. While a simpler idealized / Semi-analytic model of star formation, show that assuming exponentail gas and KS-relation with SF threshhold, that reduction of turbulance (driven by a number of processes) in outer disk leads to a steeper decline in SFR as crossing the critical SF threshold becomes less likely. 

# Understanding the extended density profiles of dwarf galaxies

Are dwarf galaxies indeed expected to be one-component exponential-like density profiles? The suggestion of exponential density profiles dates back to @faber+lin1983. Expand...

A number of recent work has also confirmed the presence of extended, likely non-tidal, density profiles. For example, @chiti+2021; @chiti+2023 spectroscopically confirm members out to $~9 R_h$ in Tucana II.  

- @revaz+jablonka2018
- @tau+vivas+martinez-vazquez2024 (RRL stars )
- @roderick+2016 Boötes I
- @mcconnachie+penarrubia+navarro2007 (velocity space & multiple populations).
- @cicuendez+battaglia2018

## Multi-epoch star formation and feedback

The star formation history of dwarf galaxies is typically thought to be *bursty*: a series of discrete episodes of intense star formation separated by periods of quiescence. 

For example, in the simulations of @wheeler+2019, etc., dwarf galaxies exhibit strongly bursty SFHs. 

- @maxwell+2012
- @wright+2019
- @azartash-namin+2024

### Multi-component density profiles

A related note is that several dwarf galaxies show evidence for multiple chemodynamical components in their stars. While spectra tend to focus on the inner regions (high probability members), the existence of multiple-component stellar populations in the inner regions of dwarf galaxies hints at a complex star formation history capable of creating the observed profiles.

- @benitez-llambay+2016
- @mercado+2021: formation of metallicity gradients from SFH 
- @revaz+jablonka2018 natural formation of multi-component and gradients from cosmological simulations of dwarf galaxies. Originate from dynamical heating of metal rich population and outside in star formation?
- @el-badry+2016. Feedback drives fluctuations in size + radial migration => population gradients.
- 

Observational evidence:

- @arroyo-polonio+2024
- @pace+2020
- @fabrizio+2016, @kordopatis+2016 (multi-component in Carina)
- @battaglia+2006, @amorisco+evans2012 for fornax

## Past mergers and accretion

- @deason+2014
- @deason+2022
- @ricotti+polisensky+cleland2022
- @querci+2025
- @amorisco+evans+vandeven2014
- @tarumi+yoshida+frebel2021: Shows cosmological simulations in context of Tuc II, merger of two galaxies in the early universe produces a more extended density profile and metallicity gradient.
- @lokas+2014: And II from merger, 
- rotation from merger: @cardona-barrero+2021



## Preprocessing

Another alternative is that the dwarf galaxies have been "preprocessed" by another dwarf/LMC like galaxy earlier on and this has resulted in their unique properties

## Tidal dwarfs

On e possibility is that Sculptor and Ursa Minor are not formed from dark matter halos but instead due to the collision of gas. This possibility is not a convincing explanation because we know the star formations histories are extended and relatively old, and the galaxies are close to apocentre, making achieving the observed velocity dispersion challenging in this framework.

## Alternatives to $\Lambda$CDM

Alternatives to the prevailing cosmology we assume may help explain the nature of these dwarf galaxies. t

- WDM: forms cores, doesn't help
- SIDM: can cause faster disruption, doesn't fix orbital path or change tidal effects?
- MOND: @sanchez-salcedo+hernandez2007: mond in dsph
