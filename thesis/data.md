An ongoing research question is whether local dwarf spheroidal galaxies contain substructure such as stellar halos. Some systems appear to be simple in structure, yet recent work has reported detection of stellar halos in several classical dwarfs. To provide an empirical benchmark for simulations in the next chapters, we examine the observed stellar structure of Sculptor and Ursa Minor. We describe the @jensen+2024 Bayesian methodology to select high-probability member stars from *Gaia*. The derived density profiles are robust to alternative selection criteria. While Fornax is well represented by an exponential density profile, we find that Sculptor and Ursa Minor exhibit an excess of stars in their outer regions relative to an exponential, with a Plummer profile providing a better fit. In following chapters,  we consider possible physical explanations for these extended profiles.



# *Gaia* sample selection

The detection of faint features in resolved stellar systems requires careful consideration of whether a given star belongs to the system. Without removing contamination from foreground/background features, faint features may be lost in the noise or be of uncertain association. Historically, background separation was performed by the colour-magnitude diagram alone (for example, matched filter methods like REF). Now with *Gaia* data, stellar parallax and proper motion are available as additional membership indicators. 

Here, we use the @jensen+2024's (hereafter J+24) membership probabilities for *Gaia* observations. J+24 use a Bayesian framework incorporating proper motion (PM), colour-magnitude diagram (CMD), and spatial information to determine the probability that a given star belongs to the satellite or background. By accounting for PM in particular, J+24 produces low contamination samples of candidate member stars out to large distances from a dwarf galaxy. J+24 extends the algorithm presented in @MV2020a; @MV2020b by additionally includes a secondary, extended spatial component. J+24 are able to find candidate members out to ~10 half-light radii ($R_h$). Similar recent work has also included @pace+li2019; @battaglia+2022; @pace+erkal+li2022; @qi+2022.

J+24 construct likelihoods, ${\cal L}$, representing the probability density that a star is consistent with either the foreground/background, ${\cal L}_{\rm bg}$, or the satellite galaxy, ${\cal L}_{\rm sat}$. In either case, the likelihoods are the product of a spatial, PM, and CMD component,
$$
{\cal L} = {\cal L}_{\rm space}\ {\cal L}_{\rm PM}\ {\cal L}_{\rm CMD}.
$$

For the satellite, the *spatial likelihood* is specified either as a one- or two-component elliptical exponential profile, including a second component only when statistically favoured. The *proper-motion likelihood* quantifies the agreement of a star's motion with the dwarf galaxy's systemic motion, accounting for observational uncertainties. The *CMD likelihood* measures the consistency of a star's *Gaia* $G$ and $G_{\rm BP}- G_{\rm RP}$ with theoretical isochrones for the galaxy. For the background likelihoods, the spatial likelihood is a uniform distribution over the field, and the background CMD and PM likelihoods are constructed empirically from the stars in an annulus surrounding the satellite. Each likelihood is normalized as a probability density over the respective parameter space.

The total likelihood is a mixture model of the satellite and background, weighted by the fraction of stars in the field belonging to the satellite, $f_{\rm sat}$: 
$$
{\cal L}_{\rm tot} = f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}.
$$
The probability that a given star belongs to the satellite is then
$$
P_{\rm sat} = 
\frac{f_{\rm sat}\,{\cal L}_{\rm sat}}{{\cal L}_{\rm tot}}
= \frac{f_{\rm sat}{\cal L}_{\rm sat}}{f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}}.
$$

J+24 fit this model using Monte Carlo Markov cain simulations, solving for the proper motions, satellite membership fraction ($f_{\rm sat}$), and structural properties of the second exponential density profile (if included). The median parameters from the samples are then used to calculate the $P_{\rm sat}$ we use for sample selection.

For our fiducial sample, we adopt a minimum probability of $P_{\rm sat} = 0.2$. Most stars have $P_{\rm  sat}$ values which are nearly 0 or 1. So, the exact choice of probability threshold has little effect on the resulting sample. Even at our relatively generous probability threshold of 0.2, the purity remains high when validated against spectroscopic line-of-sight (LOS) velocities (~90%, J+24). We note, however, in these purity estimates bay be biased. Stars with LOS velocities are typically brighter, and resultantly have more precise *Gaia* measurements. However, we find that our conclusions are unchanged when limiting samples to only the brightest stars. Altogether, the J+24 method provides a high-quality, low-contamination sample of stars for analysis into the stellar distributions of dwarf galaxies, which we will now investigate. 

# Observed stellar distributions



We use the tangent plane and elliptical radius to analyze stellar distributions while avoiding distortions due to projection effects and dwarf galaxy's elliptical shape. The tangent plane coordinates $\xi$ and $\eta$ are offsets in RA and declination as measured on the plane tangent to the galaxy centre. To account for the elliptical shape of the galaxy, we use $R_{\rm ell}$, which we define as the circularized elliptical radius,
$$
R_{\rm ell}^2 = a\,b\,({\xi'}^2 / a^2 + {\eta'}^2 / b^2),
$$
 where $\xi'$ and $\eta'$ are the tangent plane coordinate rotated to align with dwarf galaxy's major and minor axis, and $a$ and $b$ are the semi-major and semi-minor axis of the galaxy. $R_{\rm ell}$ more meaningfully represents the distance of a feature or star from the galaxy centre than circular radii, especially for highly elliptical systems like Ursa Minor. 

We illustrate the likelihood-based selection through a progressive tightening of criteria. First, we consider a minimally-refined sample, **all**, only excluding stars with poor astrometry, unreliable photometry, or inconsistent parallaxes. Next, incorporating CMD and PM information, we use the **CMD+PM** to test a selection method agnostic to position. This sample includes stars where the CMD and PM combined likelihood favours satellite membership:
$$
{\cal L}_{\rm CMD,\ sat}\ {\cal L}_{\rm PM,\ sat} > {\cal L}_{\rm CMD,\ bg}\ {\cal L}_{\rm PM,\ bg}.
$$
Our fiducial sample, **$P_{\rm sat} > 0.2$**, includes a spatial likelihood as described in J+24. Finally, the **RV members** sample adds line-of-sight velocity information  (see Appendix for details). While the most confident members and a useful validation sample, we do not use this sample for stellar density analysis due to its complex selection function. 

[@Fig:scl_selection; @Fig:umi_selection; @Fig:fornax_selection] shows our fiducial sample for each galaxy, and the impact of alternate samples, in tangent-plane coordinates, the CMD, and in proper-motion space. The *all* sample extends uniformly across the tangent plane, but includes a substantial background population. The *CMD+PM* sample has a much lower background density in the tangent plane, revealing the satellite more clearly.  The fiducial sample appears similar to the *CMD+PM* sample in the CMD and PM planes, but excludes improbably distant stars. The *RV member* sample also traces out similar distributions in CMD and PM space as the fiducial sample, with less dispersion likely reflecting the typically brighter magnitudes of stars with spectroscopic followup. Each selection component plays an important role in ensuring a high-quality membership criterion. 

Based on the fiducial sample's distribution in [@fig:scl_selection; @fig:umi_selection; @fig:fornax_selection], each galaxy appears to be more similar than different. Ursa Minor does stand out with its higher ellipticity. Sculptor and Ursa Minor show a population of members past $6R_h$, but Fornax does not, even though Fornax has many more stars. Fornax also has a more blue CMD, indicative of ongoing star formation. Otherwise, there are no standout morphological features between these galaxies. 

One limitation of the J+24 method is that it assumes a spatial form of the dwarf galaxy, which may impact the detection of distant stars or density deviations. However, there are no clear extensions or overdensities in the *CMD+PM* selected sample, which does not include a spatial likelihood. If even fainter, more extended features exist, they are likely not detectable with *Gaia*.

To illustrate the extent of Sculptor and Ursa Minor, we highlight the most distant, spectroscopically confirmed members with red stars in [@fig:scl_selection; @fig:umi_selection]. In particular, @sestito+2023a; @sestito+2023b targeted the most distant stars using the J+24 candidate membership list. These stars illustrate an extended stellar distribution for each galaxy. For an exponential profile, 99.95% of stars fall within $6R_h$, so the presence of very distant stars provides direct evidence for deviations from a simple, exponential profile. 



![Sculptor sample selection](figures/scl_selection.png){#fig:scl_selection width=70%}

Figure: The distributions of various samples of *Gaia* stars for Sculptor. We plot light grey points for field stars (with consistent parallaxes and reliable photometry and astrometry), turquoise points for *CMD+PM* selected stars ([@eq:sel_cmd_pm]), blue `x`s for the fiducial sample ($P_{\rm sat} > 0.2$), and indigo diamonds for the RV confirmed members. We mark the two stars from @sestito+2023a with rust-outlined indigo stars. **Top:** Tangent plane $\xi, \eta$. The orange ellipse represents 3 half-light radii. **Bottom left:** Colour magnitude diagram in Gaia $G$ magnitude versus $G_{\rm BP} - G_{\rm RP}$ colour. We plot the Padova 12Gyr [Fe/H]=-1.68 isochone in orange. The black bar in the top left represents the median colour error. **Bottom right:** Proper motion in declination $\mu_\delta$ vs RA $\mu_{\alpha*}$ (corrected) the orange circle represents the @MV2020b proper motion. The black cross represents the median proper motion error.



![Ursa Minor sample selection](figures/umi_selection.png){#fig:umi_selection width=70%}

Figure: Similar to @fig:scl_selection except for Ursa Minor. We outline RV members outside of $6R_h$ with red stars (from @sestito+2023b, @pace+2020 and @spencer+2018). 

![Fornax sample selection](figures/fornax_selection.png){#fig:fornax_selection width=70%}

Figure: Similar to @fig:scl_selection except for Fornax. RV members are from @WMO2009. While possibly a limitation of RV sample selection, Fornax does not show the same extended outer halo of probable members as Sculptor or Ursa Minor despite having many more stars.

# Density profiles

Density profiles summarize the radial distribution of stars and provide key constraints for our simulations in future chapters. If density profiles deviate from an exponential, then the very distant stars we noted are likely part of larger structure of the galaxy.Density profiles also are a useful to differentiate galaxies.

We derive density profiles using a simple histogram method. Stars are binned in constant width bins in $\log R_{\rm ell}$ of 0.05 dex. We ignore bins interior or exterior to the first empty bin in either direction. We use symmetric poisson uncertainties. For the reference $R_h$, we use the values from @munoz+2018 based on deeper photometry and maximum likelihood fits of Sérsic profiles. 

[@Fig:scl_observed_profiles] show the derived density profiles for Sculptor, Ursa Minor, and Fornax. We calculate density profiles for different selections of stars from above: *all* stars, CMD+PM only, and the fiducial sample. In each case, all samples are the same towards the inner regions of the satellite. Classical dwarfs dominate the stellar density in their cores. However, the density profile *all stars* plateaus at the total background in the field at radii of 30-60 arcminutes. By restricting stars to being most likely satellite members by CMD + PM, the background is reduced by 1-2 dex. The CMD+PM background plateau likely represents the density of background stars which could be mistaken as members.  Finally, the background (BG)-subtracted profile results from subtracting the apparent background in the *all* profile. 

While the density profiles show excellent agreement, the reliability past $\log R / R_h\approx 1.8$ is uncertain. After this radius, the fiducial profile continues to estimate a continuous density falloff for 1-2 magnitudes below the CMD+PM selection background. Since the primary difference between the *CMD+PM* and fiducial ($P_{\rm sat}$ ) selection is the inclusion of the spatial likelihood, we propose that the divergence here reveals the density profile likelihood specification rather than a real feature of the background. More sophisticated models in the Appendix corroborate this interpretation---there is not sufficient information to estimate the dwarf galaxy density in this regime. Therefore, we mark the density past this point for all dwarf galaxies with open circles (**TODO**). 

Nearby to UMi, there is a small, likely unassociated, ultrafaint star cluster, Muñoz 1 [@munoz+2012]. The cluster is at a relative position of $(\xi, \eta) \approx(-42, -15)$  arcminutes, corresponding to an elliptical radius of 36 arcminutes. However this cluster does not have a bright RGB, so has few stars brighter than a $G$ mag of 22. The cluster is not apparent in the *Gaia* star density. While Muñoz 1 may contribute $\sim$ 5 stars to the density profile of UMi, this would possibly only affect 1-2 bins due to its compact size.

To illustrate the differences between each dwarf galaxy, in @Fig:classical_dwarfs_densities, we compare Sculptor, Ursa Minor, and other classical dwarfs against exponential and Plummer density profiles. All dwarf galaxies appear to be well described by an exponential profile in the inner regions. However, while an exponential continues to be a good description for other dwarf galaxies, Sculptor and Ursa Minor diverge, showing a substantial excess over an exponential, better fit by a Plummer instead. In fact, Sculptor and Ursa Minor show excesses of between 2-3 orders of magnitude beyond what is expected of an exponential at the very outer regions! 



![Sculptor density profiles](figures/scl_umi_fnx_density_methods.pdf){#fig:scl_observed_profiles}

Figure: The density profile of Sculptor, Ursa Minor, and Fornax for different selection criteria, ploted as log surface dencity versus log elliptical radius.  *CMD+PM* select stars more likely to be members according to CMD and PM only, *all* selects any high quality star, and *BG subtracted* is the background-subtracted density derived from *all* stars. We mark the half-light radius (vertical dashed line) and the break radius (black arrow, REF). 





![Classical dwarf density profiles](figures/classical_dwarf_profiles.pdf){#fig:classical_dwarfs_densities}

Figure: The density profiles of Sculptor, Ursa Minor, and other classical dwarfs compared to Exp2D and Plummer density profiles. Dwarf galaxies are scaled to the same half-light radius and density at half-light radius. Sculptor and Ursa Minor have an excess of stars in the outer regions (past $\log R/R_h \sim 0.3$) compared with other classical dwarfs. 

