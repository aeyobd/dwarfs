# Abstract

To provide an empirical benchmark for simulations in the next chapters, we aim to understand the observed properties of Sculptor and Ursa Minor. In this section, we first compile past measurements of each galaxy. We describe the @jensen+2024 Bayesian methodology to select high probability member stars from *Gaia*. The derived density profiles are robust to changes in assumptions. Fornax is well represented by a 2D exponential profile. Instead, Sculptor and Ursa Minor show an excess of stars relative to an 2D exponential in the outer regions, better described my a Plummer profile. Finally, we combine various radial velocity samples with @jensen+2024's members to derive systemic velocities and dispersions for Sculptor and Ursa Minor.



# Introduction

Since the discovery of dwarf galaxies around the Milky Way, observational work has attempted to measure and refine the basic properties of these objects. While the Milky Way's satellites are close (by extragalactic standards), their low numbers of stars and large areas on the sky present challenges for observational work. The classical systems discussed here each extend about 1-2 degrees across the sky ([@fig:scl_selection; @fig:umi_selection; @fig:fornax_selection]). 



In particular, *Gaia* has provided a wealth of high quality data, including the distances to nearby stars and proper motions. This allows for some of the first measurements on the motions for many local group dwarf galaxies [e.g.; @MV2020]. Accurate distances and velocities are invaluable for understanding the orbital history of satellites.  



High-quality large spectroscopic samples combined with *Gaia* has opened new windows into understanding the internal kinematics and history of dwarf galaxies. <Discuss battaglia, pace, etc. and conclusions about dwarf galaxies.>



In [@tbl:scl_obs_props; @tbl:umi_obs_props] we present  adopted observed properties for each galaxy.



| parameter        | value                                                        | Source                                |
| ---------------- | ------------------------------------------------------------ | ------------------------------------- |
| $\alpha$         | $15.0183 \pm 0.0012^\circ$                                   | 1                                     |
| $\delta$         | $-33.7186 \pm 0.0007^\circ$                                  | 1                                     |
| distance modulus | $19.60 \pm 0.05$ (RR lyrae)                                  | 2                                     |
| distance         | $83.2 \pm 2$ kpc                                             | 2                                     |
| $\mu_{\alpha*}$  | $0.099 \pm 0.002 \pm 0.017$ mas yr$^{-1}$                    | 3                                     |
| $\mu_\delta$     | $-0.160 \pm 0.002_{\rm stat} \pm 0.017_{\rm sys}$ mas yr$^{-1}$ | 3                                     |
| radial velocity  | $111.2 \pm 0.2\ {\rm km\,s^{-1}}$                            | @arroyo-polonio+2024; @sec:rv_results |
| $\sigma_v$       | $9.7\pm0.2\ {\rm km\,s^{-1}}$                                | @arroyo-polonio+2024;                 |
| $R_h$            | $9.79 \pm 0.04$ arcmin                                       | 4                                     |
| $R_{h,inner}$    | $9.46 \pm 0.26$ arcmin                                       | XREF                                  |
| ellipticity      | $0.36 \pm 0.01$                                              | 1                                     |
| position angle   | $92\pm1^\circ$                                               | 1                                     |
| $M_V$            | $-10.82\pm0.14$                                              | 1                                     |

Table: Observed properties of Sculptor. References are: 1. @munoz+2018, 2. @tran+2022, 3. @MV2020b, 4 @MV2020a. {#tbl:scl_obs_props  short="Observed Properties of Sculptor"}



| parameter          | value                                                        | Source                          |
| ------------------ | ------------------------------------------------------------ | ------------------------------- |
| $\alpha$           | $ 227.2420 \pm 0.0045$˚                                      | 1                               |
| $\delta$           | $67.2221 \pm 0.0016$˚                                        | 1                               |
| distance modulus   | $19.23 \pm 0.11$ (RR lyrae)                                  | 2                               |
| distance           | $70.1 \pm 3.6$ kpc                                           | 2                               |
| $\mu_\alpha*$      | $-0.124 \pm 0.004 \pm 0.017$ mas yr$^{-1}$                   | 3                               |
| $\mu_\delta$       | $0.078 \pm 0.004_{\rm stat} \pm 0.017_{\rm sys}$ mas yr$^{-1}$ | 3                               |
| radial velocity    | $-245.9 \pm 0.3_{\rm stat} \pm 1_{\rm sys}$ km s$^{-1}$      | @pace+2020; see @sec:rv_results |
| $\sigma_v$         | $8.7 \pm 0.3$                                                | @pace+2020                      |
| $R_h$              | $11.62 \pm 0.1$ arcmin                                       | 1                               |
| $R_{h, \rm inner}$ | $11.7\pm0.5$                                                 |                                 |
| ellipticity        | $0.55 \pm 0.01$                                              | 1                               |
| position angle     | $50 \pm 1^\circ$                                             | 1                               |
| $M_V$              | $-9.03 \pm 0.05$                                             | 1                               |

Table: Observed properties of Ursa Minor. References are: (1) @munoz+2018, (2) @garofalo+2025, (3) @MV2020a, (4) average of @pace+2020 and @spencer+2018.  {#tbl:umi_obs_props  short="Observed Properties of Ursa Minor"}

## The *Gaia* mission

Some of the most fundamental properties of astronomical objects are their position and velocity. Unfortunately, determining distances to stars is nontrivial. Additionally, while line-of-sight or radial velocities (RVs) are easily determined from spectroscopy, the tangental velocities, perpendicular to RV, are only measurable through proper motions, typically requiring precise astrometry as well. *Gaia*'s mission is to produce extremely precise astrometry enabling measurements of unprecedented accuracy and scale for proper motions, parallaxes, and magnitudes. 

*Gaia* was designed to revolutionize proper motion and parallax measurements. *Gaia* is a space-based, all-sky survey telescope with two primary 1.45x0.5m mirrors situated at the Sun-Earth L2 lagrange point [@gaiacollaboration+2016]. *Gaia* was launched in XXX, completing its mission in 2025 (but with two more data releases planned). By imagining two patches of sky on the same focal plane, separated by a fixed angle of 106.5degrees, *Gaia* is able to measure absolute proper motions by comparing the apparent shifts of stars in different regions of the sky. In addition to precise astrometric information, *Gaia* measures the magnitude of stars in the very wide *G* band (330-1050nm), blue and red colours using the blue and red photometers (BP and RP, 330-680, 640-1050 respectively), and takes low resolution BP-RP spectra and radial velocity measurements of bright stars (magnitudes <16?).

*Gaia* has revolutionized many astronomical disciplines, the least of which is local group and Milky Way science. While proper motions of a dwarf galaxies has been measured in a case by case basis by the Hubble Space Telescope, a full systematic determination of proper motions for most dwarf galaxies was unavailable until *Gaia* [@MV2020a; @battaglia+2022]. *Gaia* has furthermore allowed for the detection of many substructures of the halo, streams, 

# Gaia Membership Selection

Here, we briefly describe J+24's membership estimation method. J+24 use a Bayesian framework incorporating proper motion (PM), colour-magnitude diagram (CMD), and spatial information to determine the probability that a given star belongs to the satellite. By accounting for PM in particular, J+24 produces low contamination samples of candidate member stars. J+24 extends the algorithm presented in @MV2020a; @MV2020b but additionally includes an optional secondary, extended spatial component to find possible members as far as ~10 half-light radii $R_h$ from some dwarf galaxies. See also similar work by @pace+li2019; @battaglia+2022; @pace+erkal+li2022; @qi+2022.

To create a high-quality sample, J+24 select stars initially from Gaia within a 2--4 degree circular region centred on the dwarf satisfying: 

- Solved astrometry, magnitude, and colour.
- Renormalized unit weight error, ${\rm ruwe} \leq 1.3$, ensuring high quality astrometry. `ruwe` is a measure of the excess astrometric noise on fitting a consistent parallax-proper motion solution [see @lindegrenXXX].
- 3$\sigma$ consistency of measured parallax with dwarf's distance (dwarf parallax is very small; with @lindegren+2021 zero-point correction).
-  Absolute proper motions, $\mu_{\alpha*}$, $\mu_\delta$, less than 10$\,{\rm mas\ yr^{-1}}$. (Corresponds to tangental velocities of $\gtrsim 500$ km/s at distances larger than 10 kpc.)
- Corrected colour excess is within expectations: $|C^*| \leq 3\,\sigma_{C^*}(G)$, with $C^*$ and $\sigma_{C^*}$ from @riello+2021. Removes stars with unreliable photometry.
- De-reddened $G$ magnitude is between $22 > G > G_{\rm TRGB} - 5\delta\rm {DM}$. Removes very faint stars and stars significantly brighter than the tip of the red giant branch (TRGB) magnitude plus the distance modulus uncertainty $\delta {\rm DM}$. 
- Colour is between $-0.5 < {\rm BP - RP} <  2.5$ (dereddened).

Photometry is dereddened with @schlegel+finkbeiner+davis1998 extinction maps.

J+24 define likelihoods ${\cal L}$ representing the probability density that a star is consistent with either the MW stellar background (${\cal L}_{\rm bg}$) or the satellite galaxy (${\cal L}_{\rm sat}$). In either case, the likelihoods are the product of a spatial, PM, and CMD term:
$$
{\cal L} = {\cal L}_{\rm space}\ {\cal L}_{\rm PM}\ {\cal L}_{\rm CMD}.
$$

Each likelihoods is normalized over their respective 2D parameter space for both the satellite. Since the likelihoods are normalized, $f_{\rm sat}$, representing the fraction of member stars in the field, controls the relative scaling between the satellite and background likelihoods. The total likelihood for any star in this model is then
$$
{\cal L}_{\rm tot} = f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}
$$ {#eq:Ltot}
The probability that any star belongs to the satellite is then given by
$$
P_{\rm sat} = \frac{f_{\rm sat}{\cal L}_{\rm sat}}{f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}}.
$$

For the satellite's spatial likelihood, J+24 consider both one-component and a two-component cases. The one component model is constructed as a single exponential profile  ($\Sigma \propto e^{R_{\rm ell} / R_s}$), with $R_s$ fixed to the value in the table in @MV2020a. Additionally, structural uncertainties (for position angle, ellipticity, and scale radius) are sampled over to construct the final likelihood map. The two-component model instead adds a second exponential, $\Sigma_\star \propto e^{-R/R_s} + B\,e^{-R/R_{\rm outer}}$. The inner scale radius is fixed, and the outer scale radius and magnitude of the second component $R_{\rm outer}$, $B$  are solved for. Structural properties are not accounted for.

The PM likelihood is a bivariate gaussian with variance and covariance equal to each star's proper motions. J+24 assume the stellar PM errors are the main source of uncertainty.

 The satellite's CMD likelihood is based on a Padova isochrone [@girardi+2002]. The isochrone has a matching metallicity and 12 Gyr age. The (gaussian) colour width is assumed to be 0.1 mag plus the Gaia colour uncertainty at each magnitude. The HB is modelled as a constant magnitude extending blue of the CMD (mean magnitude of -2.2, 12 Gyr HB stars and a 0.1 mag width plus the mean colour error). A likelihood map is constructed by sampling the distance modulus in addition to the CMD width, taking the maximum of RGB and HB likelihoods.

The background likelihoods are instead empirically constructed. Stars stars outside of 5$R_h$ passing the quality cuts estimate the background density in PM and CMD space. The density is a sum of bivariate gaussians with variances based on Gaia uncertainties (and covariance for proper motions).  The spatial background likelihood is assumed to be constant. 

J+24 derive $\mu_{\alpha*}$, $\mu_\delta$, $f_{\rm sat}$ (and $B$, $R_{\rm outer}$ for two-component) through an MCMC simulation with likelihood from [@Eq:Ltot]. Priors are weakly informative or uniform. The proper motion single component prior is same as @MV2020a: a normal distribution with mean 0 and standard deviation $100\,\kms$. If 2-component spatial, instead is a uniform distribution spanning 5$\sigma$ of single component case w/ systematic uncertainties. $f_{\rm sat}$  (and $B$) has a uniform prior  0--1. $R_{\rm outer}$ has a uniform prior only restricting $R_{\rm outer} > R_s$. The mode of each parameter from the MCMC are then reported and used to calculate the final $P_{\rm sat}$ values. 

We adopt a probability cut of $P_{\rm sat} = 0.2$ as our fiducial sample. Most stars are assigned probabilities close to either 0 or 1, so the choice of probability threshold is only marginally significant. Additionally, even for a probability cut of 0.2, the purity of the resulting sample with RV measurements is very high (~90%, J+24). Note there is likely a systematic bias in using stars with RV measurements to measure purity. Fainter stars are less likely to have been targeted and have poorer astrometry and colour information. We find no difference in the resulting density distributions when restricting stars to be brighter than a specific magnitude.

## Selected samples

In [@Fig:sculptor_selection; @Fig:umi_selection; @Fig:fornax_selection], we illustrate the resulting samples from the algorithm in the tangent plane ($\xi$, $\eta$, *does this need defined?*), CMD, and proper motion space. For all galaxies, each criteria plays a commensurate role in sifting out nonmembers. The CMD is well defined and probable members only extend a few times the colour uncertainty from the CMD. In proper motion space, selected stars are within $\masyr$ from the systematic proper motion. The stars furthest away in proper motion space typically have large uncertainties $\sim 1\masyr$, so most members are reasonably consistent. In any case, the algorithm does enable removing many stars which do not have consistent proper motions with the dwarf. Finally, the spatial likelihood reduces the probabilities of stars distant from the dwarfs centres. Member stars are rare outside of $\sim5R_h$, resulting from both the likelihood specification and the lower density of background stars consistent with the CMD and PM of the satellite. 



One limitation of the J+24 method is that it assumes a spatial form of the dwarf galaxy, possibly limiting the detection of very distant stars or features. [@Fig:sculptor_selection; @Fig:umi_selection] also illustrates the distribution of stars selected without a spatial criterion. We define the *CMD+PM* selection as
$$
{\cal L}_{\rm CMD,\ sat}\ {\cal L}_{\rm PM,\ sat} > {\cal L}_{\rm CMD,\ bg}\ {\cal L}_{\rm PM,\ bg}
$$ {#eq:sel_cmd_pm}
with the likelihoods from J+24 as described above. These stars are distributed similar to the fiducial (probable members) sample but with a fainter uniform distribution across the entire field. This illustrates the approximate background of stars which may be confused as members. Additionally, since there is no clear spatial structure in the *CMD+PM* (or *all*) sample outside several $R_h$, it is unlikely that there are additional faint, tidal features detectable with *Gaia*. Not shown here, we also try a variety of simpler, absolute cuts and thresholds, finding no evidence of other structure in both *Gaia* and DELVE/UNIONS data. 

Finally, we illustrate the location of RV-confirmed members from @sec:rv_obs in each panel. Because RV targets tend to be brighter than a typical *Gaia* candidate, RV members typically have more precise PMs. However, the RV members fill about the same area on the CMD down to the magnitude limit. Additionally, while the number of stars with RV measurements decreases at large radii, there are still confirmed RV members as far as $\sim 10R_h$ for both galaxies. These stars illustrate the extended profiles of each galaxy. Note that for an exponential profile, 99.95% of stars should be within $6R_h$, so the discovery of any stars beyond $6R_h$ hints at a deviation from an exponential profile. These *extratidal* stars hint at something tidal in nature, and our goal is to explore possible interpretations. 





![Sculptor sample selection](figures/scl_selection.pdf){#fig:sculptor_selection}

Figure: The selection criteria for Scl members. Light grey points represent field stars (satisfying quality criteria), turquoise points CMD+PM selected stars ([@eq:sel_cmd_pm]), blue crosses probable members (2-component), and RV confirmed members indigo diamonds. We mark the two stars from @sestito+2023a with rust-outlined indigo stars. **Top:** Tangent plane $\xi, \eta$. The orange ellipse represents 3 half-light radii. **Bottom left:** Colour magnitude diagram in Gaia $G$ magnitude versus $G_{\rm BP} - G_{\rm RP}$ colour. We plot the Padova 12Gyr MH=-1.68 isochone in orange. The black bar in the top left represents the median colour error. **Bottom right:** Proper motion in declination $\mu_\delta$ vs RA $\mu_{\alpha*}$ (corrected) the orange circle represents the @MV2020b proper motion. The black cross represents the median proper motion error.



![Ursa Minor sample selection](figures/umi_selection.pdf){#fig:umi_selection}

Figure: Similar to @Fig:sculptor_selection except for Ursa Minor. We outline RV members outside of $6R_h$ in black stars (from @sestito+2023b, @pace+2020 and @spencer+2018). 



![Fornax sample selection](figures/fornax_selection.pdf){#fig:fornax_selection}

Figure: Similar to @Fig:sculptor_selection except for Fornax. RV members are from APOGEE and @WMO2009 with same methodology as Sculptor and Ursa Minor. The isochrone is instead 2Gyr old -0.99 metallicity.

## Density Profiles

Density profiles are an essential observational constraint for our later simulations. To derive density profiles, we use 0.05 dex bins in log radius. We remove bins at smaller (larger) radii than the first bin to contain no stars, working outwards. We use symmetric poisson uncertainties. As discussed below, these uncertainties are straightforward but are likely under-representative. 

[@Fig:scl_observed_profiles; @Fig:umi_observed_profiles] show the derived density profiles for Scl and UMi. We calculate density profiles for different selections of stars from above: all quality stars, CMD+PM only, and the probable member samples. In each case, all samples are the same towards the inner regions of the satellite. Classical dwarfs dominate the stellar density in their cores. However, the density profile *all stars* plateaus at the total background in the field at radii of 30-60 arcminutes. By restricting stars to being most likely satellite members by CMD + PM, the background is reduced by 1-2 dex. The CMD+PM background plateau likely represents the density of background stars which could be mistaken as members.  Finally, the background (BG)-subtracted profile results from subtracting the apparent background in the *all* profile. 

Note that the probable members (fiducial) density profile continues to confidently estimate the density profile below the CMD+PM background. These points are likely unreliable (see discussion below). However, before this point, both the BG subtracted and probable members density profiles are strikingly similar. Assumptions about the details of the likelihood and spatial dependence have marginal influence on the resulting density profile when the satellite is higher density than the background. Thus the detection of excesses of stars in Sculptor (past 20 arcmin) and UMi (past 15-20 arcmin) are robust.

Nearby to UMi, there is a small, likely unassociated, ultrafaint star cluster, Muñoz 1 [@munoz+2012]. The cluster is at a relative position of $(\xi, \eta) \approx(-42, -15)$  arcminutes, corresponding to an elliptical radius of 36 arcminutes. However this cluster does not have a bright RGB, so has few stars brighter than a $G$ mag of 22. While Muñoz 1 may contribute $\sim$ 5 stars to the density profile of UMi, this would possibly only affect 1-2 bins due to its compact size.

![Sculptor density profiles](figures/scl_umi_fnx_density_methods.pdf){#fig:scl_observed_profiles}

Figure: The density profile of Sculptor, Ursa Minor, and Fornax for different selection criteria, ploted as log surface dencity versus log elliptical radius. Residuals are with respect to th interpolated *probable members* density. *Probable members* selects stars with PSAT > 0.2 considering PM, CMD, and spatial, *CMD+PM* select stars more likely to be members according to CMD and PM only, *all* selects any high quality star, and *BG subtracted* is the background-subtracted density derived from high-quality stars. We mark the half-light radius ( vertical dashed line) and the break radius (black arrow, REF). 

## Caveats

The J+24 method was designed in particular to detect the presence of a density excess and find individual stars at large radii to be followed up. We are more interested in accurately quantifying the density profile and size of any perturbations. One potential problem with using J+24's candidate members is that the algorithm assumes the density is either described by a single or double exponential. If this model does not accurately match the actual density profile of the dwarf galaxy, we would like to understand how strongly influenced the density profile is by this assumption. 

In particular, in @Fig:umi_observed_profiles, notice that the PSAT method produces small errorbars, even when the density is $>1$dex below the local background. These stars are likely selecting stars from the statistical MW background consistent with UMi PM / CMD, recovering the assumed density profile. As a result, the reliability of these density profiles below the CMD+PM background may be questionable. A more robust analysis, removing this particular density assumption, would be required to more appropriately represent the knowledge of the density profile as the background begins to dominate. 

J+24 do not account for structural uncertainties in dwarfs for the two component case. We assume constant ellipticity and position angle. Dwarf galaxies, in reality, are not necessarily smooth and constant. They do test an alternative method using circular radii for the extended density component, and we find these density profiles are very similar. We also assume a constant ellipticity and position angle here.

Finally, our approximation to the poisson uncertainties Additionally, a more self consistent model would fit the density profile to the entire field at once, eliminating possible misrepresentation of the uncertainties.

While *Gaia* has shown excellent performance, some notable limitations may introduce problems in our interpretation and reliability of density profiles.

Gaia systematics in proper motions and parallaxes are typically smaller than the values for sources of magnitudes $G\in[18,20]$. Since we use proper motions and parallaxes as general consistency with the dwarf, and factor in systematic uncertainties in each case, these effects should not be too significant. However, the systematic proper motion uncertainties becomes the dominant source of uncertainty in the derived systemic proper motions of each galaxy (see [@tbl:scl_obs_props; @tbl:umi_obs_props]).

*Gaia* shows high but imperfect completeness, particularly showing limitations in crowded fields and for faint sources ($G\gtrapprox20$). As discussed in @fabricius+2021, fir the high stellar densities in globular clusters, the completeness relative to HST varies significantly with the stellar density. However, the typical stellar densities of dwarf galaxies are much lower, at about 20 stars/arcmin = 90,000 stars / degree, lower than the lowest globular cluster densities and safely below the crowding limit of 750,000 objects/degree for BP/RP photometry. , where the completeness down to $G\approx 20$ is $\sim 80\%$. Closely separated stars pose problems for Gaia's on-board processing, as the pixel size is 59x177 mas on the sky. This results in a reduction of stars separated by less than 1.5" and especially for stars separated by less than 0.6 arcseconds. The astrometric parameters of closely separated stars furthermore tends to be of lower quality [@fabricius+2021]. However, even for the denser field of Fornax, only about 3% of stars have a neighbour within 2 arc seconds, so multiplicity should not affect completeness too much (except for unresolved binaries). One potential issue is that the previous analysis do not account for our cuts on quality and number of astrometric parameters. These could worsen completeness, particularly since the BP-RP spectra are more sensitive to dense fields. In **REF**, we test if magnitude cuts impact the resulting density profiles, finding that this is likely not an issue. 



# Comparison and conclusions

To illustrate the differences between each dwarf galaxy, in @Fig:classical_dwarfs_densities, we compare Sculptor, Ursa Minor, and Fornax against 2D-exponential and Plummer density profiles (REF). While all dwarfs appear similar in the inner regions, each dwarf diverges in the outer regions relative to an exponential. Relative to an exponential, Fornax is underdense but Sculptor and Ursa Minor are both overdense. A Plummer profile instead provides a more reasonable fit to Scl and UMi.

In summary, we have used J+24 data to derive the density profiles for Fornax, Sculptor, and Ursa Minor. In each case, the density profile is robust against different selection criteria. Both Sculptor and Ursa Minor show strong evidence for deviations from an exponential profile. We also compile velocity measurements to derive the systemic motions and velocity dispersions of each galaxy. We find evidence for a velocity gradient in Sculptor of $4.3\pm1.3\,\kmsdeg$. We find no evidence of additional (velocity or stellar) substructure in either galaxy. Our goal in the following chapters is to test if tides provide a viable explanation for the observed properties of Sculptor and Ursa Minor.



![Classical dwarf density profiles](figures/scl_umi_fornax_exp_fit.pdf){#fig:classical_dwarfs_densities}

Figure: The density profiles of Sculptor, Ursa Minor, and Fornax compared to Exp2D and Plummer density profiles. Dwarf galaxies are scaled to the same half-light radius and density at half-light radius (fit from the inner 3 scale radii exponential recursively. ). Sculptor and Ursa Minor have an excess of stars in the outer regions (past $\log R/R_h \sim 0.3$) compared with Sculptor. Sculptor and Ursa Minor's density profiles are created with the 2-component J+24 model, but the excess does not change significantly for the 1-component model.



# Appendix

## Density profile tests





![Density profiles](/Users/daniel/thesis/figures/scl_density_methods_extra.pdf){#fig:sculptor_observed_profiles}

Figure: Density profiles for various assumptions for Sculptor. PSAT is our fiducial 2-component J+24 sample, circ is a 2-component bayesian model assuming circular radii, simple is the series of simple cuts described, bright is the sample of the brightest half of stars (scaled by 2), DELVE is a sample of  RGB  stars (background subtracted and rescaled to match).



![Density profiles](figures/scl_density_methods_j24.pdf){#fig:sculptor_observed_profiles_j24} 

Figure: Comparison of density profiles for each J+24 method. The circ

![UMi Density profiles](figures/umi_density_methods_extra.pdf){#fig:umi_observed_profiles}

Figure: Similar to [@fig:sculptor_observed_profiles] except for Ursa Minor

![UMi density methods](figures/umi_density_methods_j24.pdf){#fig:umi_observed_profiles_j24} 

Note that a full rigorous statistical analysis would require a simulation study of injecting dwarfs into Gaia and assessing the reliability of various methods of membership and density profiles. This is beyond the scope of this thesis. 



```
SELECT TOP 1000
       *
FROM delve_dr2.objects
WHERE 11 < ra
and ra < 19
and -37.7 < dec
and dec < -29.7

```

