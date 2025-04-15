# Gaia Membership Selection

Gaia provides unprecedented accuracy in proper motions and magnitudes. As such, Gaia data is uniquely excellent to produce low-contamination samples of likely member stars belonging to satellites. Here, we breifly describe J+24's membership estimation and discuss how this informs our observational knoledge of each galaxies density profile. In general, J+24 use a Bayesian framework incorporating proper motion (PM), colour-magnitude diagram (CMD), and spatial information to determine the probability that a given star belongs to the satellite. J+24 extends @MV2020a (see also @pace+li2019, @battaglia+2022, @pace+erkal+li2022, etc.).

J+24 select stars initially from Gaia within a 2--4 degree circular region centred on the dwarf satisfying:

- Solved parallax, proper motions, colour, and magnitudes. 
- High quality astrometry (`ruwe <= 1.3`)
- 3$\sigma$ consistency of measured parallax with dwarf distance + uncertainty (typically near zero; with @lindegren+2018 zero-point correction).
-  Absolute RA and Dec proper motions less than 10$\,{\rm mas\ yr^{-1}}$ (corresponding to tangental velocities of $\gtrsim 500$ km/s at distances larger than 10 kpc.)
- No colour excess (@lindegren+2018 equation C.2)
- G < 22 and less than 5$\sigma$ above TRGB, and between -0.5 and 2.5 in Bp - Rp.  

Photometry is dereddened using @schlegel+1988 extenction. 

J+24 calculate the probability that any star belongs to either the satellite or the MW background as
$$
P_{\rm sat} = \frac{f_{\rm sat}{\cal L}_{\rm sat}}{f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}}
$$

where the satellite (sat) and background (bg) likelihoods are simply the product of the PM, CMD, and spatial components: 
$$
{\cal L} = {\cal L}_{\rm space}\ {\cal L}_{\rm PM}\ {\cal L}_{\rm CMD}
$$

The satellite likelihood is constructed as

- CMD: The CMD is the lowest metallicity isochrone from Padova [@girardi+2002] with age 12 Gyr with a colour width of 0.1 mag  plus the Gaia colour uncertainty at each magnitude. The HB is modelled as a constant magnitude extending blue of the CMD. The HB magnitude is the mean magnitude of HB stars from most metal poor isochrone with a 0.1 mag width plus the mean colour error. A likelihood map is constructed by sampling the distance modulus in addition to the CMD width, taking the maximum of RGB and HB likelihoods.
- Spatial: A single exponential ($\Sigma \propto e^{R_{\rm ell} / R_s}$) accounting for structural uncertainty (sampled over position angle, ellipticity, and half light radius).
- Alternative spatial: For Scl and UMi, this is instead a double exponential $\Sigma_\star \propto e^{-R/R_s} + B\,e^{-R/R_{\rm outer}}$ where the inner exponential remains fixed. Structural parameter uncertainties are not accounted for.
- PM. A bivariate gaussian with variance and covariance equal to each star's proper motions. Each star's proper motions uncertainty are assumed to be the dominant uncertainty. 

The background likelihood is constructed as:

- CMD : Constructed as a density map using the other quality-selected stars outside of $5R_h$ in the catalogue. The map is a sum of bivariate gaussians for each star with standard deviations based on the Gaia uncertainties. 
- PM: same as CMD except in PM space.
- Spatial: a constant likelihood. 

Note that each likelihood map is normalized over the respective 2D parameter space. In order to represent the difference in frequency of background and forground stars, $f_{\rm sat}$ represents the field fraction of member stars.  

In J+24, a MCMC simulation is ran using the above likelihood to solve for the following parameters

- Systemic proper motions $\mu_\alpha$, $\mu_\delta$.  Single component prior is same as @MV2020: a normal distribution with mean 0 and standard deviation 100 km/s. If 2-component spatial, instead is a uniform distribution  spanning 5$\sigma$ of single component case w/ systematic uncertainties.
- $f_{\rm sat}$ density normalization. Prior is a uniform distribution between 0 and 1.
- Spatial component parameters $B$ is uniform from 0-1 and $R_{\rm outer}$ is uniform and greater than $R_s$ for extended profiles (Scl and UMi here.)

The mode of each parameter from the MCMC are then used to calculate the final $P_{\rm sat}$ values we use here. 

We adopt a probability cut of $P_{\rm sat} = 0.2$ as our fiducial sample. Most stars are assigned probabilities close to 0 or 1, so the choice of probability threshhold is not too significant. Additionally, even for a probability cut of 0.2, the purity of the resulting sample with RV measurements is very high (~90%, J+24). (Note: there is likely a high systematic bias in using stars with RV measurements to measure purity. Fainter stars have poorer measurements and distant stars are less likely to have been targeted. )

## Resulting Samples

In figures @fig:sculptor_selection, @fig:umi_selction, we illustrate the resulting samples from the algorithm. In each case, each criteria plays a roll: proper motions are centred around the dwarf systemic motion, the CMD is well defined, and stars only within a few $R_h$ are included. We also plot the RV members found in general and in @sestito+2024, @sestito+2024b.

The tangent plane, $\xi$, $\eta$, is the projection that

We also illustrate the approximate result of removing the spatial component from the likelihood. We define the CMD+PM selection as stars satisfying
$$
{\cal L}_{\rm CMD,\ sat}\ {\cal L}_{\rm PM,\ sat} > {\cal L}_{\rm CMD,\ bg}\ {\cal L}_{\rm PM,\ bg}
$$


These stars appear similar to the fiducial (probable members) sample, but instead also appear as an approximately uniform distribution across the entire field. This illustrates the approximate background of stars which may be confused as members. Additionally, since there is no clear spatial structure in the CMD+PM sample, it is unlikely that there are tidal tails detectable with Gaia. Not shown here, we also try a variety of simpler, absolute cuts and thresholds, finding no extended structure beyond what is detected in J+24. 

This means that at least at the level of where the background density dominates, we can exclude models which produce tidal tails brighter than a density of $\Sigma_\star \approx 10^{-2}\,\text{Gaia-stars\ arcmin}^{-2} \approx 10^{-6} \, {\rm M_\odot\ kpc^{-2}}$ (TODO assuming a distance  of ... and stellar mass of ...). 



![Sculptor selection criteria](figures/scl_selection.png){#fig:sculptor_selection}

Figure: The selection criteria for Scl members. Probable members (2-component) are orange, and all field stars (satisfying quality criteria) are in light grey. **Top:** Tangent plane. **Bottom left:** Colour magnitude diagram. **Bottom right:** Proper motion. 



![Ursa Minor Selection](figures/umi_selection.png){#fig:umi_selection}

Figure: Similar to @fig:sculptor_selection except for Ursa Minor. UMi features a very extended density profile with some stars ~ 6$R_h$ including a RV member. UMi is also highly elliptical compared to other classical dwarfs. 



## Density Profiles

Our primary observational constraint is the density profile of a dwarf galaxy. 

To derive density profiles, we use 0.05 dex bins in log radius (i.e. the bins are derived from 10^(minimum(logR):0.05:maximum(logR))). The density in each bin is then (from Poisson statistics)
$$
\Sigma_b = N_{\rm memb} / A_{\rm bin} \pm \sqrt{N_{\rm memb}} / A_{\rm bin}
$$
where $N_{\rm memb}$ is the number of members in the bin and $A_{\rm bin}$ is the area of the bin's annulus in 2D. As discussed below, these uncertainties underrepresent the true uncertainty on multiple accounts. We retain Poisson errors for simplicity here. 

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



# Radial velocity modeling

For both Sculptor and Ursa Minor, we construct samples of radial velocity measurements and cross match each sample to produce estimates of the total velocity dispersion and line of sight velocity for each galaxy. 

## Methodology

Our RV methodology is

- Cross match all catalogues to J+24 Gaia members. If Gaia ID's are included, we use these. Otherwise, we match to the nearest star within 2 arcseconds. 
- We combine the mean RV measurement from each study using the inverse-variance weighted mean and standard uncertainty. 

$$
v_i = \frac{1}{\sum w}\sum_j w\ v_{i, j} \\
\delta v_i = \sqrt{\frac{1}{\sum w^2}}
$$

where $w_i = 1/s_i^2$, and we estimate the inter-study standard deviation
$$
s_{v}^2 = \frac{1}{\sum w} \sum w (v_{i, j} - v_i)^2
$$
Note that we do not include the standard deviation if a source does measure the velocity multiple times.

If the total standard deviation is greater than five times the velocity standard uncertainty (reduced), then we remove the star since either the star is a binary or contains erroneous measurements or other uncertainties. Note that since stars typically only have very few measurements, these choices are highly approximate and likely limited.

The combined RV likelihood is then
$$
{\cal L} = {\cal L}_{\rm space} {\cal L}_{\rm CMD} {\cal L}_{\rm PM} {\cal L}_{\rm RV}
$$
where
$$
{\cal L}_{\rm RV} = N(\mu_{v}, \sigma_{v}^2 + (\delta v_i)^2)
$$
where $\mu_v$ and $\sigma_v$ are the systemic velocity and dispersion, and $\delta v_i$ is the individual measurement uncertainty. Typically, the velocity dispersion will dominate the uncertainty budget here. 

Similar to above, we retain stars with the resulting membership probability of greater than 0.2.



Finally, we need to correct the coordinate frames. We need a proper motion estimate of each star. Since most stars have Gaia measured proper motions, but the uncertainties on each proper motion measurement are large, we combine the system and individual star's proper motion by weighting the systemic proper motion by the z-score of the consistency of the measurement with the systemic proper motion. 
$$
\mu_{*}' = \frac{1}{w_\star + w_{\rm sys}}\left(\mu_{\star}w_\star + \mu_{\rm sys} w_{\rm sys}\right)
$$
where we model 
$$
{w_{\rm sys}} = P_{\rm sat}(\delta \mu ^2 + \sigma_{v, tan}^2)^{-1}
$$
and the weight of the star is just the inverse variance weight plus the gaia systematics. 

The corrected radial velocities are then 
$$
v_z = v_{\rm los, gsr}\cos\theta  + v_{\rm tan, R} \sin\theta
$$
where $v_{\rm tan, R} = d(\mu_{\alpha*}\cos\theta + \mu_\delta \sin \theta)$ is the radial component of the proper motion with respect to the centre of the galaxy. Since $\theta \lesssim 1^\circ$, this correction is small ($\sim 1$km\,s$^{-1}$), but can prevent the introduction of a spurious gradient on the sky [@WMO2008]. The uncertainty is then the velocity uncertainty plus the distance uncertainties times the PM uncertainty from above. We then use the $v_z$ values for the following modelling, however repeating with plain, solar-frame velocities does not substantially affect the results too much . 

## Sculptor

For Sculptor, we combine radial velocity measurements

- APOGEE (crossmatched to stars meeting basic quality cuts in field)

- @sestito+2024a

- @tolstoy+2024

  

We adopt a very loose prior of 
$$
\mu_{v} = N(111, 10) \\
\log\sigma_{v} = N(1, 1)
$$
Our resulting velocity dispersions are
$$

$$
Comparing to past work, note that our values are reasonable. 



## Ursa Minor

- APOGEE
- @sestito+2024b
- @pace+2022
- @spencer+2018



## Discussion and limitations

Our model here is relatively simple. Some things which we do not account for include

- interstudy systematics and biases
- Inappropriate uncertainty reporting.
- Binarity. Binay stars are expected to induce an average velocity dispersion of xxx on an individual star. 
- Selection effects. RV studies each have their own targeted program and 

We partially address several of these by also performing the analysis on large RV samples for each galaxy. This allows for individual comparison illustrating the differences in selection effects. 
