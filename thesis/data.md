# Gaia Membership Selection

Gaia provides unprecedented accuracy in proper motions and magnitudes. As such, Gaia data is uniquely excellent to produce low-contamination samples of likely member stars belonging to satellites. Here, we breifly describe J+24's membership estimation and discuss how this informs our observational knoledge of each galaxies density profile. In general, J+24 use a Bayesian framework incorporating proper motion (PM, colour-magnitude diagram (CMD), and spatial information to determine the probability that a given star belongs to the satellite. J+24 extends @MV2020a (see also @pace+li2019, etc.).

J+24 select stars initially from Gaia satisfying: 

- High quality astrometry (`ruwe <= 1.3`)
- 3$\sigma$ consistency of measured parallax with dwarf distance + uncertainty (near zero with @lindegren+2018 zero-point correction)
-  Absolute RA and Dec proper motions less than 10$\,{\rm mas\ yr^{-1}}$
- Good photometry
- No colour excess (@lindegren+2018 equation C.2)
- G > 22 and less than 5$\sigma$ above TRGB, and between -0.5 and 2.5 in Bp - Rp.  

J+24 calculate the probability that any star belongs to either the satellite or the MW background as
$$
P_{\rm sat} = \frac{f_{\rm sat}{\cal L}_{\rm sat}}{f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}}
$$

where the satellite (sat) and background (bg) likelihoods are simply the product of the PM, CMD, and spatial components: 
$$
{\cal L} = {\cal L}_{\rm space}\ {\cal L}_{\rm PM}\ {\cal L}_{\rm CMD}
$$

The satellite likelihood is constructed as

- CMD: The CMD is from  Padova [@girardi+2002] isochrone at 12 Gyr with a colour width of 0.1 mag  plus the Gaia colour uncertanty and the observed metallicity of the dwarf. The HB is modelled as a constant magnitude extending blue of the CMD with a 0.1 mag width plus the mean colour error. A likelihood map is constructed by sampling the distance modulus in addition to the CMD width, taking the maximum of RGB and HB likelihoods.
- Spatial: A single exponential ($\Sigma \propto e^{R_{\rm ell} / R_s}$) using mean values of the structural parameters (**verify**). For Scl and UMi, this is instead a double exponential $\Sigma_\star \propto e^{-R/R_s} + B\,e^{-R/R_{\rm outer}}$ where the inner exponential remains fixed and structural parameters are not accounted for. Structural parameter uncertainties are not accounted for.
- PM. A bivariate gaussian with variance and covariance equal to each star's proper motions. Each star's proper motions uncertainty are assumed to be the dominant uncertainty. 

The background likelihood is constructed as:

- CMD : Constructed as a KDE using the other quality-selected stars outside of $5R_h$ in the catalogue. The width of each star is projected as a bivariate Gaussian.
- PM: same as CMD except in PM space.
- Spatial: a constant likelihood. 

Note that each likelihood map is normalized over the respective parameter space. In order to represent the difference in frequency of background and forground stars, $f_{\rm sat}$ represents the field fraction of member stars.  

In J+24, a MCMC simulation is ran using the above likelihood to solve for the following parameters

- Systemic proper motions $\mu_\alpha$, $\mu_\delta$. Prior: within 5$\sigma$ of single component case w/ systematic uncertainties. Single component prior is based on MV2020 (RV members / expected velocity dispersion of halo at distance?)
- $f_{\rm sat}$ density normalization (uniform between 0 and 1)
- Spatial component parameters $B$ and $R_{\rm outer}$ for extended profiles (Scl and UMi here.)

The median parameters from the simulation are then used to calculate the final $P_{\rm sat}$ values we use here. 

We adopt a probability cut of $P_{\rm sat} = 0.2$ as our fiducial sample. Most stars are assigned probabilities close to 0 or 1, so the choice of probability threshhold is not too significant. Additionally, even for a probability cut of 0.2, the purity of the resulting sample with RV measurements is very high (~90%, J+24). (Note: there is likely a high systematic bias in using stars with RV measurements to measure purity. Fainter stars have poorer measurements and distant stars are less likely to have been targeted. )

## Resulting Samples



![Sculptor selection criteria](figures/scl_selection.png){#fig:sculptor_selection}

Figure: The selection criteria for Scl members. Probable members (2-component) are orange, and all field stars (satisfying quality criteria) are in light grey. **Top:** Tangent plane. **Bottom left:** Colour magnitude diagram. **Bottom right:** Proper motion. 



![Ursa Minor Selection](figures/umi_selection.png){#fig:umi_selection}

Figure: Similar to @fig:sculptor_selection except for Ursa Minor. UMi features a very extended density profile with some stars ~ 6$R_h$ including a RV member. UMi is also highly elliptical compared to other classical dwarfs. 







Figure: Similar to @fig:sculptor_selection except for Fornax. Note that we only use one-component probabilities here. Fornax is the quintessential unperturbed dwarf galaxy. 

### Searches for tidal tails

- There are no apparent overdensities in the PM & CMD only selected stars to suggest the presence of a tidal tail

- We have tried selective matched filter 

- This means that at least at the level of where the background density dominates, we can exclude models which produce tidal tails brighter than a density of $\Sigma_\star \approx 10^{-2}\,\text{Gaia-stars\ arcmin}^{-2} \approx 10^{-6} \, {\rm M_\odot\ kpc^{-2}}$ (TODO assuming a distance  of ... and stellar mass of ...). 

  

## Density Profiles

Our primary observational constraint is the density profile of a dwarf galaxy. 

To derive density profiles, we use 0.05 dex bins in log radius (i.e. the bins are derived from 10^(minimum(logR):0.05:maximum(logR))). The density in each bin is then (from Poisson statistics)
$$
\Sigma_b = N_{\rm memb} / A_{\rm bin} \pm \sqrt{N_{\rm memb}} / A_{\rm bin}
$$
where $N_{\rm memb}$ is the number of members in the bin and $A_{\rm bin}$ is the area of the bin. As discussed below, these uncertainties underrepresent the true uncertainty on multiple accounts. We retain Poisson errors for simplicity here. 



In Figures @fig:scl_observed_profiles, @fig:umi_observed_profiles, we show the derived density profiles for each galaxy for samples similar to in the selection plots above. In each case, all samples are the same towards the inner regions of the satellite, illustrating that these density profiles are dominated by satellite stars in the centre. The sample containing all stars reaches a plateau at the total background in the field. However, restricting stars to being most likely satellite members by CMD + PM, the background is much lower. This plateau likely represents the real background of background stars which could be mistaken as members.  Finally, we have the probable members and bg-subtracted densities. BG subtracted is based on the `all` density profile, subtracting the mean background density for stars beyond the last point of the BG subtracted profile. 

Note that the probable members (fiducial) density profile continues to confidently estimate the density profile below the CMD+PM likely star background. These points are likely unreliable (see discussion below). However, before this point, both the BG subtracted and probable members density profiles are strikingly similar. Assumptions about the details of the likelihood and spatial dependence have marginal influence on the resulting density profile when the satellite is higher density than the background. 

![Sculptor density profiles](figures/scl_density_methods.png){#fig:scl_observed_profiles}

Figure: The density profile of Sculptor for different selection criteria. *probable members* selects stars with PSAT > 0.2 considering PM, CMD, and spatial, *CMD+PM* select stars more likely to be members according to CMD and PM only, *all* selects any high quality star, and *BG subtracted* is the background-subtracted density derived from high-quality stars. 



![Ursa Minor density profiles](figures/umi_density_methods.png){#fig:umi_observed_profiles}

Figure: Similar to @fig:scl_observed_profiles except for Ursa Minor.



![Fornax density profiles](figures/fornax_density_methods.png){#fig:fornax_observed_profiles}

Figure: Similar to @fig:scl_observed_profiles except for Fornax.



# Comparison of the Classical dwarfs

Classical dwarfs are often the brightest dwarfs in the sky.c The density profiles of classical dwarf galaxies is thus well measured, enabling detailed comparisons. 



| galaxy     | R_h (exp inner 3Rs) | num cand |
| ---------- | ------------------- | -------- |
| Fornax     | $17.8\pm0.6$        | 23,154   |
| Leo I      | $3.7\pm0.2$         | 1,242    |
| Sculptor   | $9.4\pm0.3$         | 6,888    |
| Leo II     | $2.4\pm0.3$         | 347      |
| Carina     | $8.7\pm0.4$         | 2,389    |
| Sextans I  | $20.2 \pm0.9$       | 1,830    |
| Ursa Minor | $11.7 \pm 0.5$      | 2,122    |
| Draco      | $7.3\pm0.3$         | 1,781    |

Using the same methods above, we select members from J+24's stellar probabilities. We use the one-component exponential density profiles with 



Figure: Density profiles for each dwarf galaxy. Here, we use the 1-component exponential stellar probabilities from J+24. Dwarf galaxies are scaled by our derived $R_h$ values.

## Density Profile Reliability and Uncertainties

- How well do we know the density profiles? 
- What uncertainties affect derived density profiles? 
- Can we determine if Gaia, structural, or algorithmic systematics introduce important errors in derived density profiles?
- Using J+24 data, we validate
  - Check that PSAT, magnitude, no-space do not affect density profile shape too significantly
- Our "high quality" members all have > 50 member stars and do not depend too highly on the spatial component, mostly corresponding to the classical dwarfs



J+24's algorithm takes spatial position into account, assuming either a one or two component exponential density profile. When deriving a density profile, this assumption may influence the derived density profile, especially when the galaxy density is fainter than the background of similar appearing stars. To remedy this and estimate where the background begins to take over, we also explore a cut based on the likelihood ratio of only the CMD and PM components. This is in essence assuming that the spatial position of a star contains no information on it's membership probability (a uniform distribution like the background)

## Caveats

The J+24 method was designed to determine high probability members for spectroscopic followup in particular. Note that we instead care about retrieving a reliable density profile.

In particular, in @fig:umi_observed_profiles, notice that the PSAT method produces artifically small errorbars even when the density is >1dex below the local background. These stars are likely selecting stars from the statistical MW background consistent with UMi PM / CMD, recovering the assumed density profile. As a result, the reliability of faint features in these density profiles is questionable and a more robust analysis, removing this particular density assumtion, would be required to more appropriately represent the knowledge of the density profile as the background begins to dominate. 

J+24 do not account for structural uncertainties in dwarfs. This is a not insignificant source of uncertainty in the derived density profile

We assume constant ellipticity and position angle. Dwarf galaxies, in reality, are not necessarily smooth and constant.



# Comparison and conclusions



To illustrate the differences between each dwarf galaxy, in @fig:classical_dwarfs_densities, we compare Scl, UMi, and Fnx against exponential and plummer density profiles (**TODO: state these somewhere**). While all dwarfs have marginal differences in the inner regions, each dwarf diverges in the outer regions relative to an exponential. In particular, while Fnx is underdense, Scl and UMi are both overdense, approximately fitting a Plummer density profile instead. 

In summary, we have used J+24 data to derive the density profiles for Fornax, Sculptor, and Ursa Minor. In each case, the density profile is robust against different selection criteria. Both Sculptor (Ursa Minor) show strong (weak) evidence for deviations from an exponential profile. We will explore a tidal explanation for these features in this work.



![Classical dwarf density profiles](figures/scl_umi_fornax_exp_fit.png){#fig:classical_dwarfs_densities}

Figure: The density profiles of Sculptor, Ursa Minor, and Fornax compared to Exp2D and Plummer density profiles. Dwarf galaxies are scaled to the same half-light radius and density at half-light radius (fit from the inner 3 scale radii exponential recursively. )

# Appendix / Extra Notes

## Additional density profile tests



![Density profiles](figures/scl_density_methods_extra.png){#fig:sculptor_observed_profiles}

Figure: Density profiles for various assumptions for Sculptor. PSAT is our fiducial 2-component J+24 sample, circ is a 2-component bayesian model assuming circular radii, simple is the series of simple cuts described in Appendix ?, bright is the sample of the brightest half of stars (scaled by 2), DELVE is a sample of  RGB  stars (background subtracted and rescaled to match).

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



