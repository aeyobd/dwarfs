To provide an empirical benchmark for simulations in the next chapters, we aim to understand the observed properties and density profiles of Sculptor and Ursa Minor. In this section, we first compile past measurements of each galaxy. We then describe the @jensen+2024 Bayesian methodology to select high probability member stars from *Gaia*. We show the derived density profiles are robust to changes in assumptions. While Fornax is well represented by a 2D exponential profile, Sculptor and Ursa Minor show an excess of stars relative to an 2D exponential in the outer regions, better described by a Plummer profile. Future chapters assess possible explanations for these excess.





# *Gaia* data

Here, we briefly describe @jensen+2024's (hereafter J+24) membership estimation method. J+24 use a Bayesian framework incorporating proper motion (PM), colour-magnitude diagram (CMD), and spatial information to determine the probability that a given star belongs to the satellite. By accounting for PM in particular, J+24 produces low contamination samples of candidate member stars out to large distances from a dwarf galaxy. J+24 extends the algorithm presented in @MV2020a; @MV2020b but additionally includes an optional secondary, extended spatial component to find possible members as far as ~10 half-light radii $R_h$ from several dwarf galaxies. See also similar work by @pace+li2019; @battaglia+2022; @pace+erkal+li2022; @qi+2022.

J+24 define likelihoods ${\cal L}$ representing the probability density that a star is consistent with either the MW stellar background (${\cal L}_{\rm bg}$) or the satellite galaxy (${\cal L}_{\rm sat}$). In either case, the likelihoods are the product of a spatial, PM, and CMD term:
$$
{\cal L} = {\cal L}_{\rm space}\ {\cal L}_{\rm PM}\ {\cal L}_{\rm CMD}.
$$

Each likelihood is normalized over their respective 2D parameter space for both the satellite. To control the relative frequency of member and background stars, $f_{\rm sat}$ representing the fraction of member stars in the field. The total likelihood for any star in this model is the sum of the satellite and background likelihoods, weighted by their relative frequencies,
$$
{\cal L}_{\rm tot} = f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}.
$$
The probability that any star belongs to the satellite is then given by 
$$
P_{\rm sat} = 
\frac{f_{\rm sat}\,{\cal L}_{\rm sat}}{{\cal L}_{\rm tot}}
= \frac{f_{\rm sat}{\cal L}_{\rm sat}}{f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}}.
$$



# Selected stars

In [@Fig:scl_selection; @Fig:umi_selection; @Fig:fornax_selection], we illustrate the resulting samples from the algorithm in the tangent plane, CMD, and proper motion space. The tangent plane coordinates $\xi$ and $\eta$ are the distances in RA and declination as measured on a plane tangent to the dwarf galaxies centre. $R_{\rm ell}$ is defined as the circularized elliptical radius, so $R_{\rm ell}^2 = a\,b\,({\xi'}^2 / a^2 + {\eta'}^2 / b^2)$ where $\xi'$ and $\eta'$ are rotated to align with the position angle direction $\theta$. 

For our fiducial sample, we adopt a probability cut of $P_{\rm sat} = 0.2$. Most stars are assigned probabilities close to either 0 or 1, so the choice of probability threshold is only marginally significant. Additionally, even for a probability cut of 0.2, the purity of the resulting sample with RV measurements is very high (~90%, J+24). Note there is likely a systematic bias in using stars with RV measurements to measure purity. Fainter stars are less likely to have been targeted and have poorer astrometry and colour information. We find no difference in the resulting density distributions when restricting stars to be brighter than a specific magnitude.

For all galaxies, each criteria is commensurate in sifting out nonmembers. The CMD is well defined, and probable members only extend a few times the colour uncertainty from the CMD. In proper motion space, selected stars are within $\masyr$ from the systematic proper motion. The stars furthest away in proper motion space typically have large uncertainties $\sim 1\masyr$, so most members are reasonably consistent. In any case, the algorithm removes stars with inconsistent PMs. Finally, the spatial likelihood reduces the probabilities of stars distant from the dwarfs centres. Member stars are rare outside of $\sim5R_h$, resulting from both the likelihood specification and the lower density of background stars consistent with the CMD and PM of the satellite. 



One limitation of the J+24 method is that it assumes a spatial form of the dwarf galaxy, possibly limiting the detection of very distant stars or features. [@Fig:scl_selection; @Fig:umi_selection] also illustrates the distribution of stars selected without a spatial criterion. We define the *CMD+PM* selection as
$$
{\cal L}_{\rm CMD,\ sat}\ {\cal L}_{\rm PM,\ sat} > {\cal L}_{\rm CMD,\ bg}\ {\cal L}_{\rm PM,\ bg}
$$ {#eq:sel_cmd_pm}
with the likelihoods from J+24 as described above. These stars are distributed similar to the fiducial (probable members) sample but with a background uniform distribution across the entire field. This illustrates the approximate background of stars which may be confused as members. Additionally, since there is no clear spatial structure in the *CMD+PM* (or *all*) sample outside several $R_h$, it is unlikely that there are additional faint, tidal features detectable with *Gaia*. Not shown here, we also try a variety of simpler, absolute cuts and thresholds, finding no evidence of other structure in both *Gaia* and DELVE/UNIONS data. 

Finally, we illustrate the location of RV-confirmed members from @sec:rv_obs in each panel. Because RV targets tend to be brighter than a typical *Gaia* candidate, RV members typically have more precise PMs. However, the RV members fill about the same area on the CMD down to the magnitude limit. Additionally, while the number of stars with RV measurements decreases at large radii, there are still confirmed RV members as far as $\sim 10R_h$ for both galaxies. These stars illustrate the extended profiles of each galaxy. Note that for an exponential profile, 99.95% of stars should be within $6R_h$, so the discovery of any stars beyond $6R_h$ hints at a deviation from an exponential profile. These *extratidal* stars hint at something tidal in nature, and our goal is to explore possible interpretations. 





![Sculptor sample selection](figures/scl_selection.png){#fig:scl_selection width=70%}

Figure: The selection criteria for Scl members. Light grey points represent field stars (satisfying quality criteria), turquoise points CMD+PM selected stars ([@eq:sel_cmd_pm]), blue xs probable members (2-component), and RV confirmed members indigo diamonds. We mark the two stars from @sestito+2023a with rust-outlined indigo stars. **Top:** Tangent plane $\xi, \eta$. The orange ellipse represents 3 half-light radii. **Bottom left:** Colour magnitude diagram in Gaia $G$ magnitude versus $G_{\rm BP} - G_{\rm RP}$ colour. We plot the Padova 12Gyr [Fe/H]=-1.68 isochone in orange. The black bar in the top left represents the median colour error. **Bottom right:** Proper motion in declination $\mu_\delta$ vs RA $\mu_{\alpha*}$ (corrected) the orange circle represents the @MV2020b proper motion. The black cross represents the median proper motion error.



![Ursa Minor sample selection](figures/umi_selection.png){#fig:umi_selection width=70%}

Figure: Similar to @fig:scl_selection except for Ursa Minor. We outline RV members outside of $6R_h$ in black stars (from @sestito+2023b, @pace+2020 and @spencer+2018). 

![Fornax sample selection](figures/fornax_selection.png){#fig:fornax_selection width=70%}

Figure: Similar to @fig:scl_selection except for Fornax. RV members are from @WMO2009. While possibly a limitation of RV sample selection, Fornax does not show the same extended outer halo of probable members as Sculptor or Ursa Minor despite having many more stars.

# Density profiles

Density profiles are an essential observational constraint for our later simulations. To derive density profiles, we use 0.05 dex bins in log radius. We remove bins at smaller (larger) radii than the first empty bin working outwards. We use symmetric poisson uncertainties. As discussed below, these uncertainties are straightforward but are likely under-representative. We use the values from @munoz+2018's Sérsic maximum likelihood fits for $R_h$, which are more precisely derived given their deeper photometry.

[@Fig:scl_observed_profiles] show the derived density profiles for Sculptor, Ursa Minor, and Fornax. We calculate density profiles for different selections of stars from above: all quality stars, CMD+PM only, and the probable member samples. In each case, all samples are the same towards the inner regions of the satellite. Classical dwarfs dominate the stellar density in their cores. However, the density profile *all stars* plateaus at the total background in the field at radii of 30-60 arcminutes. By restricting stars to being most likely satellite members by CMD + PM, the background is reduced by 1-2 dex. The CMD+PM background plateau likely represents the density of background stars which could be mistaken as members.  Finally, the background (BG)-subtracted profile results from subtracting the apparent background in the *all* profile. 

Note that the probable members (fiducial) density profile continues to confidently estimate the density profile below the CMD+PM background. These points are likely unreliable (see discussion below). However, before this point, both the BG subtracted and probable members density profiles are strikingly similar. Assumptions about the details of the likelihood and spatial dependence have marginal influence on the resulting density profile when the satellite is higher density than the background. Thus the detection of excesses of stars in Sculptor (past 20 arcmin) and UMi (past 15-20 arcmin) are robust.

Nearby to UMi, there is a small, likely unassociated, ultrafaint star cluster, Muñoz 1 [@munoz+2012]. The cluster is at a relative position of $(\xi, \eta) \approx(-42, -15)$  arcminutes, corresponding to an elliptical radius of 36 arcminutes. However this cluster does not have a bright RGB, so has few stars brighter than a $G$ mag of 22. While Muñoz 1 may contribute $\sim$ 5 stars to the density profile of UMi, this would possibly only affect 1-2 bins due to its compact size.

To illustrate the differences between each dwarf galaxy, in @Fig:classical_dwarfs_densities, we compare Sculptor, Ursa Minor, and other classical dwarfs against 2D-exponential and Plummer density profiles (REF). All dwarf galaxies appear to be well described by an exponential profile in the inner regions. However, while an exponential continues to be a good description for other dwarf galaxies, Sculptor and Ursa Minor diverge, showing a substantial excess over an exponential, better fit by a Plummer instead. In fact, Sculptor and Ursa Minor show excesses of between 2-3 orders of magnitude beyond what is expected of an exponential at the very outer regions! 



![Sculptor density profiles](figures/scl_umi_fnx_density_methods.pdf){#fig:scl_observed_profiles}

Figure: The density profile of Sculptor, Ursa Minor, and Fornax for different selection criteria, ploted as log surface dencity versus log elliptical radius.  *CMD+PM* select stars more likely to be members according to CMD and PM only, *all* selects any high quality star, and *BG subtracted* is the background-subtracted density derived from *all* stars. We mark the half-light radius (vertical dashed line) and the break radius (black arrow, REF). 





![Classical dwarf density profiles](figures/classical_dwarf_profiles.pdf){#fig:classical_dwarfs_densities}

Figure: The density profiles of Sculptor, Ursa Minor, and other classical dwarfs compared to Exp2D and Plummer density profiles. Dwarf galaxies are scaled to the same half-light radius and density at half-light radius. Sculptor and Ursa Minor have an excess of stars in the outer regions (past $\log R/R_h \sim 0.3$) compared with other classical dwarfs. 



# Summary

In summary, we have used J+24 data to derive the density profiles for the classical dwarf galaxies. In each case, the density profile is robust against different selection criteria. Both Sculptor and Ursa Minor show strong evidence for deviations from an exponential profile. We find no evidence of additional (velocity or stellar) substructure in either galaxy. Our density profiles are consistent with many previous works, all indicating divergence from a King or Exponential density profile with a break between 30 and 60 arcminutes. Our goal in the following chapters is to test if tides provide a viable explanation for the observed properties of Sculptor and Ursa Minor.





