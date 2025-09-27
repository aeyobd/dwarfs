As discussed in [@sec:exponential_profiles], the projected luminosity/stellar mass density profile of dwarf galaxies is generally well-described by an exponential law ([@eq:exponential_law]). A prototypical example is the Fornax dSph — well-fit by an exponential over 4.5 decades in surface density. On the other hand, dSphs like Sculptor or Ursa Minor have profiles which deviate significantly from an exponential profile fitted to the inner regions, with a clear excess of stars/light in the outer regions. In this chapter, we critically review the density profile of the MW classical dwarfs.  This thesis concerns the origin of this excess, with a tidal interpretation explored in subsequent chapters.

# Satellite stellar membership with *Gaia*

Measuring the light profile of a resolved galaxy requires careful consideration of whether any given star belongs to the system or not. Without removing contamination from foreground/background sources, faint features may be lost in the noise or be of uncertain association. When only photometric data was available, the membership of stars was ascertained using the colour-magnitude diagram alone [e.g., matched filter methods like those used by @rockosi+2002]. Now that *Gaia* data are available, stellar parallax and proper motion are also available to improve membership assignment. 

Here, we use the @jensen+2024's (hereafter J+24) membership probabilities from *Gaia* data. J+24 used a Bayesian framework incorporating proper motion (PM), colour-magnitude diagram (CMD), and spatial information to determine the probability that a given star belongs to the satellite or foreground/background. By accounting for PM in particular, J+24 produced low contamination samples of candidate member stars out to large distances from a dwarf galaxy. J+24 extended the algorithm presented in @MV2020a; @MV2020b by additionally including a secondary, extended spatial component. J+24 detected candidate members out to ~10 half-light radii from the centres of some galaxies ($R_h$). Similar recent work has also included @pace+li2019; @battaglia+2022; @pace+erkal+li2022; @qi+2022.

J+24's method uses likelihoods, ${\cal L}$, representing the probability density that a star is consistent with either the foreground/background, ${\cal L}_{\rm bg}$, or the satellite galaxy, ${\cal L}_{\rm sat}$. In either case, the likelihoods are the product of a spatial, PM, and CMD component,
$$
{\cal L} = {\cal L}_{\rm space}\ {\cal L}_{\rm PM}\ {\cal L}_{\rm CMD}.
$$

For a satellite, the *spatial likelihood* is specified either as a one- or two-component elliptical exponential profile. A second component is included only when the preferred amplitude of the second component is non-zero. The *proper-motion likelihood* quantifies the agreement of a star's motion with the dwarf galaxy's systemic motion, accounting for observational uncertainties. The *CMD likelihood* measures the consistency of a star's *Gaia* $G$ magnitude and $G_{\rm BP}- G_{\rm RP}$ colour with theoretical isochrones for the galaxy. For the background likelihoods, the spatial likelihood is a uniform distribution over the field, and the background CMD and PM likelihoods are constructed empirically from the stars in an annulus far from the satellite. Each likelihood is normalized as a probability density over the respective parameter space. Appendix [-@sec:density_extra] discusses the details of quality cuts and the likelihood calculation.

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

J+24 fit this model using Monte Carlo Markov chain simulations, solving for the proper motions, satellite membership fraction ($f_{\rm sat}$), and structural properties of the second exponential density profile (if included). The median parameters from the samples are then used to calculate the $P_{\rm sat}$ we use for sample selection. 

For our fiducial sample, we adopt a minimum probability of $P_{\rm sat} = 0.2$. We do not filter on magnitudes explicitly, but J+24's quality cuts typically only include stars with $G < 21$. We use the $P_{\rm sat}$ values from the elliptical 2-component runs if a galaxy shows evidence for an outer component, the 1-component run otherwise. Most stars have $P_{\rm  sat}$ values which are nearly 0 or 1, so the exact choice of probability threshold has little effect on the resulting sample. Even at our relatively generous probability threshold of 0.2, the purity remains high when validated against spectroscopic line-of-sight (LOS) velocities (~90%, J+24). We note, however, that these purity estimates may be biased. Stars with LOS velocities are typically brighter and, as a result, have more precise *Gaia* measurements. However, we find that our conclusions are unchanged when limiting samples to only the brightest stars. Altogether, the J+24 method provides a high-quality, low-contamination sample of dwarf galaxy candidate member stars, which we will now investigate. 

# The effects of membership criteria 

We  analyze stellar distributions in the tangent plane and consider the projected shape of a galaxy. The tangent plane coordinates $\xi$ and $\eta$ are offsets in right ascension (RA) and declination as measured on the plane tangent relative to the galaxy centre. To account for the elliptical shape of the galaxy, we use $R_{\rm ell}$, which we define as the circularized elliptical radius,
$$
R_{\rm ell}^2 = a\,b\,\left(\frac{{\xi'}^2}{a^2} + \frac{{\eta'}^2}{b^2} \right),
$$
 where $\xi'$ and $\eta'$ are the tangent plane coordinate rotated to align with a dwarf galaxy's major and minor axis, and $a$ and $b$ are the semi-major and semi-minor axis of the galaxy.

We illustrate the likelihood-based membership selection criteria through a progressive tightening of criteria. First, we consider a minimally-refined sample, **all**, which only excludes stars with poor astrometry, unreliable photometry, or inconsistent parallaxes. Next, incorporating CMD and PM information, we use the **CMD+PM** to test a selection method agnostic to the spatial position. This sample includes stars where the CMD and PM combined likelihood favours satellite membership:
$$
{\cal L}_{\rm CMD,\ sat}\ {\cal L}_{\rm PM,\ sat} > {\cal L}_{\rm CMD,\ bg}\ {\cal L}_{\rm PM,\ bg}.
$$ {#eq:sel_cmd_pm}
Next, our **fiducial** sample, $P_{\rm sat} > 0.2$, includes a spatial likelihood as described in J+24. Finally, the **RV members** (for radial velocity) sample adds line-of-sight velocity information (see Appendix [-@sec:rv_obs] for details). While the latter is the best motivated membership sample, we do not use this sample for stellar density analysis due to its incompleteness and complex selection function. 

[@Fig:scl_selection; @Fig:umi_selection; @Fig:fornax_selection] show our fiducial sample for three dSphs (Sculptor, Ursa Minor, and Fornax), as well as impact of different samples, in tangent-plane coordinates, the CMD, and proper-motion space. The "all" sample extends uniformly across the tangent plane, but includes a substantial background population. The "CMD+PM" sample has a much lower background density in the tangent plane, revealing the satellite more clearly.  The fiducial sample appears similar to the "CMD+PM" sample in the CMD and PM planes, but excludes improbably distant stars. The "RV member" sample also traces out similar distributions in CMD and PM space as the fiducial sample, with less dispersion, likely reflecting the brighter magnitudes of stars with spectroscopic follow-up. Each selection criteria—spatial, CMD, and PM—contributes towards a high-quality membership assignment. 

Based on the fiducial sample's distribution in [@fig:scl_selection; @fig:umi_selection; @fig:fornax_selection], Ursa Minor stands out because of its higher ellipticity. While Sculptor and Ursa Minor show an extended population of members past $6R_h$ in radius Fornax does not despite having many more stars. Fornax also has a bluer CMD, indicative of more recent star formation. Otherwise, there are no major morphological differences between these galaxies. 

One limitation of the J+24 method is the assumption of a specific density profile for the dwarf galaxy (a one- or two-component exponential), which may impact the membership probability of distant stars. However, there are no clear extensions or overdensities in the "CMD+PM" selected sample, which does not include a spatial likelihood. If even fainter, more extended features exist, they are not obviously detectable with *Gaia*.

To illustrate the extent of Sculptor and Ursa Minor, we highlight distant, spectroscopically confirmed members with red-outlined stars in [@fig:scl_selection; @fig:umi_selection]. In particular, @sestito+2023a; @sestito+2023b targeted distant, bright stars from the J+24 candidate membership list. The most distant stars confirmed in their work lie at radii $\gtrsim 10 R_h$ from the centres of each dwarf.  For an exponential profile, 0.05% of stars fall outside a radius of $6R_h$ and less than $1/10^6$ fall outside $10R_h$---the mere presence of very distant member stars provides evidence for deviations from an exponential profile. 



![Sculptor sample selection](figures/scl_selection.png){#fig:scl_selection}

Figure: The distributions of various samples of *Gaia* stars for Sculptor. In each panel, we plot light grey points for the "all" sample of stars, turquoise points for "CMD+PM" selected stars ([@eq:sel_cmd_pm]), blue squares for the "fiducial" sample ($P_{\rm sat} > 0.2$), and indigo diamonds for the "RV members" sample. We mark the two far-outlier stars from @sestito+2023a with rust-outlined indigo stars. **Left:** Tangent plane $\xi, \eta$. The orange ellipses represents 3 and 6 half-light radii. **Top right:** Colour magnitude diagram in Gaia $G$ magnitude versus $G_{\rm BP} - G_{\rm RP}$ colour. We plot a Padova $12\,$Gyr, ${\rm [Fe/H] }=-1.68$ isochone in orange. The black bar in the top left represents the median colour error. **Bottom right:** Proper motion in declination $\mu_\delta$ vs RA $\mu_{\alpha*}$ (corrected). The orange point marks the systemic @MV2020b proper motion. The black cross represents the median proper motion error. 



![Ursa Minor sample selection](figures/umi_selection.png){#fig:umi_selection}

Figure: Similar to @fig:scl_selection except for Ursa Minor. We outline "velocity confirmed" members outside a radius of $6R_h$ with red stars [from @sestito+2023b; @pace+2020; @spencer+2018]. We also mark the location of Muñoz 1 with a pink circle.

![Fornax sample selection](figures/fornax_selection.png){#fig:fornax_selection}

Figure: Similar to @fig:scl_selection except for Fornax. The isocrone is instead for a 2Gyr stellar population. RV measurements are from @WMO2009. Fornax does not show the same extended outer halo of probable members as Sculptor or Ursa Minor despite having many more stars.

# Density profiles

We derive density profiles by binning member stars in constant width bins in $\log R_{\rm ell}$ of 0.05 dex. We ignore bins interior or exterior to the first empty bin in either direction. We use symmetric Poisson uncertainties as error bars in the density estimate at each bin. 

[@Fig:scl_observed_profiles] show the derived density profiles for Sculptor, Ursa Minor, and Fornax.  We calculate density profiles for three different samples from above: "all", "CMD+PM", and the "fiducial" sample. In each case, all samples coincide towards the inner regions of the satellite. However, the density profile from "all" plateaus at the total background in the field at radii of 30-60 arcminutes, depending on the galaxy. By restricting stellar membership with CMD + PM, the background is reduced by 1-2 dex. The "CMD+PM" background plateau represents the density of background stars which could be mistaken as members because of their coincident colours and proper motions.  Finally, the "all – background" profile results from subtracting the apparent background in the "all" profile. 

While the density profiles of all samples agree in the inner regions, they begin to deviate outside $6R_h$. Outside this radius, the fiducial profile extends  1–2 magnitudes below the "CMD+PM" selection background. Because the fiducial profile in this region depends critically on the spatial likelihood, its shape is vulnerable to the assumed density profile (an exponential in J+24). We explore this in more detail in Appendix [-@sec:density_extra].

Near UMi, there is a small ($R_h\sim 0.5'$), likely unassociated, ultrafaint star cluster, Muñoz 1 [@munoz+2012]. The cluster is at a relative position of $(\xi, \eta) \approx(-42, -15)$  arcminutes, corresponding to an elliptical radius of 37 arcminutes. Muñoz 1 does not have a bright RGB, so has few stars above the *Gaia* magnitude limit. The cluster has little effect on the elliptically-averaged density profile (see location on [@fig:scl_observed_profiles]).

Our density profiles are robust to alternative methodologies and magnitude cuts. Density profiles based on J+24 candidates may be biased by their assumed spatial likelihoods. To address this, we consider a Bayesian model with a non-parametric spatial likelihood. We also explore selections using absolute cuts, and the effect of structural assumptions of the density profile (Appendix -@sec:extra_density). These methodologies are consistent except when the derived density drops below the background of satellite-like stars (from the CMD+PM sample). We also do not find evidence of magnitude biases in the density profiles, indicating inhomogeneous magnitude-completeness is small. Our density profiles furthermore agree with photometric surveys and similar literature (see @sec:extra_density). We conclude that selection criteria do not influence our results.



@Fig:classical_dwarfs_densities compares the fiducial density profiles of  Sculptor, Ursa Minor, Fornax, and other classical dwarf galaxies. Of the classicals, we exclude Antlia II, due to the extremely high background, and Sagittarius, which was not included in J+24. The density profiles are scaled to match at the half-light radius, taken from @munoz+2018. All of the classical dwarfs appear to be well described by an exponential profile in the inner regions.[^differing_depth] In the outer regions, however, Sculptor and Ursa Minor deviate, and show a clear outer excess over an exponential law (solid black line). These galaxies are better fit by a Plummer law (dashed black line). The deviation from an exponential grows outwards, and at $\sim 8 R_h$, may reach 2 orders of magnitude. The remainder of this thesis will be devoted to assessing whether the outer excess shown by Scl and UMi are due to Galactic tides. 





[^differing_depth]: Comparing density profiles of dwarf galaxies is complicated by variations in the effective depth between galaxies. However, the deviations between Scl, UMi and other galaxies is apparent even when data covers all dwarfs. 



![Density profiles for different Gaia samples](figures/scl_umi_fnx_density_methods.pdf){#fig:scl_observed_profiles}

Figure: The density profile of Sculptor, Ursa Minor, and Fornax for different selection criteria, ploted as log surface dencity versus log elliptical radius.  The samples are: "all" selects any high quality and parallax-consistent star (grey circles),"all--background" subtracts a uniform background density from the "all" profile (orange pentagons), "CMD+PM" select stars according to CMD and PM only (turquoise triangles), and "fiducial" also includes spatial information (blue squares). We mark the half-light radius with a vertical dashed line and the background density with the horizontal grey line. For Ursa Minor, we show the expected location of Muñoz 1 stars as a horizontal bar (ranging from $\pm3$ times the cluster's half-light radii). 





![Classical dwarf density profiles](figures/classical_dwarf_profiles.pdf){#fig:classical_dwarfs_densities}

Figure: The density profiles of Sculptor, Ursa Minor, and other classical dwarfs compared to Exp2D and Plummer density profiles. Dwarf galaxies are scaled to the same half-light radius and density at half-light radius. Sculptor and Ursa Minor have an excess of stars in the outer regions (past $\log R/R_h \sim 0.3$) compared with other classical dwarfs. 





