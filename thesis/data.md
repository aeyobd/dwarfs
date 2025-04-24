# Gaia Membership Selection

Gaia provides unprecedented accuracy in proper motions and magnitudes. Gaia data is uniquely excellent to produce low-contamination samples of likely member stars belonging to satellites. Here, we breifly describe J+24's membership estimation and discuss how this informs our observational knoledge of each galaxies density profile. In general, J+24 use a Bayesian framework incorporating proper motion (PM), colour-magnitude diagram (CMD), and spatial information to determine the probability that a given star belongs to the satellite. J+24 extends @MV2020a (see also @pace+li2019, @battaglia+2022, @pace+erkal+li2022, etc.).

J+24 select stars initially from Gaia within a 2--4 degree circular region centred on the dwarf. J+24 only consider stars with:

- Solved parallax, proper motions, colour, and magnitudes.
- High quality astrometry (`ruwe <= 1.3`).
- 3$\sigma$ consistency of measured parallax with dwarf distance + uncertainty (typically near zero; with @lindegren+2018 zero-point correction).
-  Absolute RA and Dec proper motions less than 10$\,{\rm mas\ yr^{-1}}$ (corresponding to tangental velocities of $\gtrsim 500$ km/s at distances larger than 10 kpc.).
- No colour excess (@lindegren+2018 equation C.2).
- 22 > G > 5$\sigma$ brighter than TRGB (in practice astrometry ruwe => 21 > G).
- Between -0.5 and 2.5 in Bp - Rp.

Photometry is dereddened using @schlegel+1988 extinction maps.

J+24 calculate the probability that any star belongs to either the satellite or the MW background with
$$
P_{\rm sat} = \frac{f_{\rm sat}{\cal L}_{\rm sat}}{f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}}
$$

where $f_{\rm sat}$ is the fraction of candidate members in the field, and ${\cal L}_{\rm sat}$ and ${\cal L}$_{\rm bg}$ are the satellite (sat) and background (bg) likelihoods respectively. Each likelihood is the product of a spatial, proper motion, and CMD term:
$$
{\cal L} = {\cal L}_{\rm space}\ {\cal L}_{\rm PM}\ {\cal L}_{\rm CMD}.
$$

The satellite likelihoods are specified as

- CMD: The CMD is the lowest metallicity isochrone from Padova [@girardi+2002] with age 12 Gyr with a colour width of 0.1 mag  plus the Gaia colour uncertainty at each magnitude. The HB is modelled as a constant magnitude extending blue of the CMD. The HB magnitude is the mean magnitude of HB stars from most metal poor isochrone with a 0.1 mag width plus the mean colour error. A likelihood map is constructed by sampling the distance modulus in addition to the CMD width, taking the maximum of RGB and HB likelihoods.
- Spatial: A single exponential ($\Sigma \propto e^{R_{\rm ell} / R_s}$) accounting for structural uncertainty (sampled over position angle, ellipticity, and half light radius).
- Alternative spatial: For Scl and UMi, this is instead a double exponential $\Sigma_\star \propto e^{-R/R_s} + B\,e^{-R/R_{\rm outer}}$ where the inner exponential remains fixed. Structural parameter uncertainties are not accounted for.
- PM. A bivariate gaussian with variance and covariance equal to each star's proper motions. Each star's proper motions uncertainty are assumed to be the dominant uncertainty. 

The background likelihoods are constructed as:

- CMD : Constructed as a density map using the other quality-selected stars outside of $5R_h$ in the catalogue. The map is a sum of bivariate gaussians for each star with standard deviations based on the Gaia uncertainties. 
- PM: same as CMD except in PM space ($\mu_{\alpha*}$ and $\mu_\delta$ with respecive uncertainties and covariance.)
- Spatial: a constant likelihood. 

Each likelihood map is normalized over the respective 2D parameter space. $f_{\rm sat}$ thus is the only term controlling the total abundance of satellite stars relative to the background.

In J+24, a MCMC simulation is ran using the above total likelihood to solve for the following parameters:

- Systemic proper motions $\mu_\alpha$, $\mu_\delta$.  Single component prior is same as @MV2020: a normal distribution with mean 0 and standard deviation 100 km/s. If 2-component spatial, instead is a uniform distribution  spanning 5$\sigma$ of single component case w/ systematic uncertainties.
- $f_{\rm sat}$ density normalization. Prior is a uniform distribution between 0 and 1.
- Spatial component parameters $B$ is uniform from 0-1 and $R_{\rm outer}$ is uniform and greater than $R_s$ for extended profiles (Scl and UMi here.)

The mode of each parameter from the MCMC are then used to calculate the final $P_{\rm sat}$ values we use here. 

We adopt a probability cut of $P_{\rm sat} = 0.2$ as our fiducial sample. Most stars are assigned probabilities close to 0 or 1, so the choice of probability threshhold is not too significant. Additionally, even for a probability cut of 0.2, the purity of the resulting sample with RV measurements is very high (~90%, J+24). However,  there is likely a systematic bias in using stars with RV measurements to measure purity. Fainter stars have poorer astrometry and are less likely to have been targeted. We find no difference in the resulting density distributions when restricting stars to be brighter than a specific magnitude.

## Selected samples

In figures @fig:sculptor_selection, @fig:umi_selction, we illustrate the resulting samples from the algorithm in the tangent plane ($\xi$, $\eta$, *does this need defined?*), CMD, and proper motion space. For both galaxies, each criteria plays a commensurate role in sifting out members. Proper motions are clustered around the dwarf systemic motion, the CMD is well defined including the horizontal branch, and stars only within a few $R_h$ are included. We also mark stars with consistent radial velocities (RVs) (see below), which trace each feature albeit more cautiously. 

Figure REF also illustrats the distribution of stars selected without a spatial criterion. We define the CMD+PM selection as stars satisfying
$$
{\cal L}_{\rm CMD,\ sat}\ {\cal L}_{\rm PM,\ sat} > {\cal L}_{\rm CMD,\ bg}\ {\cal L}_{\rm PM,\ bg}
$$
These stars are distributed similar to the fiducial (probable members) sample but with a fainter uniform distribution across the entire field. This illustrates the approximate background of stars which may be confused as members. Additionally, since there is no clear spatial structure in the CMD+PM sample, it is unlikely that there are additional faint, tidal features detectable with Gaia. Not shown here, we also try a variety of simpler, absolute cuts and thresholds, finding no extended structure beyond what is detected in J+24. (*Is it worth using this to calculate an OOM upper limit on the density of tidal tails?*).



![Sculptor selection criteria](figures/scl_selection.png){#fig:sculptor_selection}

Figure: The selection criteria for Scl members. Probable members (2-component) are orange, and all field stars (satisfying quality criteria) are in light grey. **Top:** Tangent plane. We outline in star symbols the two stars from @sestito+2023a. **Bottom left:** Colour magnitude diagram. **Bottom right:** Proper motion. 



![Ursa Minor Selection](figures/umi_selection.png){#fig:umi_selection}

Figure: Similar to @fig:sculptor_selection except for Ursa Minor. We outline RV members outside of $3R_h$ in black stars (from @sestito+2023b and @pace+2020 and @spencer+2018).



## Density Profiles

Density profiles are an essential observational constraint for our later simulations. To derive density profiles, we use 0.05 dex bins in log radius (i.e. the bins are derived from 10^(minimum(logR):0.05:maximum(logR))). We remove bins at smaller (larger) radii than the first bin to contain no stars, working outwards. We use poisson uncertainties ($\delta \Sigma_i / \Sigma_i = 1/\sqrt{N_i}$ for a derived density $\Sigma_i$ in a bin with $N_i$ stars). As discussed below, these uncertainties are straightforward but are likely under-represented. 

@fig:scl_observed_profiles, @fig:umi_observed_profiles show the derived density profiles for Scl and UMi. We calculate density profiles for different selections of stars from above: all quality stars, CMD+PM only, and the probable member samples. In each case, all samples are the same towards the inner regions of the satellite. Classical dwarfs dominate the stellar density in their cores. However, the density profile *all stars* plateaus at the total background in the field at radii of 30-60 arcminutes. By restricting stars to being most likely satellite members by CMD + PM, the background is reduced by 1-2 dex. The CMD+PM background plateau likely represents the density of background stars which could be mistaken as members.  Finally, the background (BG)-subtracted profile results from subtracting the apparent background in the *all* profile. 

Note that the probable members (fiducial) density profile continues to confidently estimate the density profile below the CMD+PM background. These points are likely unreliable (see discussion below). However, before this point, both the BG subtracted and probable members density profiles are strikingly similar. Assumptions about the details of the likelihood and spatial dependence have marginal influence on the resulting density profile when the satellite is higher density than the background. Thus the detection of excesses of stars in Sculptor (past 20 arcmin) and UMi (past 15-20 arcmin) are robust.

One interesting note is for UMi, there is a small ultrafaint star cluster at a position of (-100, -15) arcmin @munoz+2012. However this cluster does not have a bright RGB so has few stars brighter than g mag of 22. While munoz 1 may contribute about 4 stars to the density profile of UMi, this would only apply to the bin at XXX because of the small half light radius. We do not subtract this out. 

![Sculptor density profiles](figures/scl_density_methods.png){#fig:scl_observed_profiles}

Figure: The density profile of Sculptor for different selection criteria. *probable members* selects stars with PSAT > 0.2 considering PM, CMD, and spatial, *CMD+PM* select stars more likely to be members according to CMD and PM only, *all* selects any high quality star, and *BG subtracted* is the background-subtracted density derived from high-quality stars. 



![Ursa Minor density profiles](figures/umi_density_methods.png){#fig:umi_observed_profiles}

Figure: Similar to @fig:scl_observed_profiles except for Ursa Minor.



![Fornax density profiles](figures/fornax_density_methods.png){#fig:fornax_observed_profiles}

Figure: Similar to @fig:scl_observed_profiles except for Fornax.

## Caveats

The J+24 method was designed in particular to detect the presence of a density excess and find individual stars at large radii to be followed up. We are more interested in accurately quantifying the density profile and size of any perturbations. One potential problem with using J+24's candidate members is that the algorithm assumes the density is either described by a single or double exponential. If this model does not accurately match the actual density profile of the dwarf galaxy, we would like to understand how strongly influenced the density profile is by this assumption.

In particular, in @fig:umi_observed_profiles, notice that the PSAT method produces artifically small errorbars even when the density is >1dex below the local background. These stars are likely selecting stars from the statistical MW background consistent with UMi PM / CMD, recovering the assumed density profile. As a result, the reliability of these density profiles below the CMD+PM background may be questionable. A more robust analysis, removing this particular density assumption, would be required to more appropriately represent the knowledge of the density profile as the background begins to dominate. 

J+24 do not account for structural uncertainties in dwarfs for the two component case. We assume constant ellipticity and position angle. Dwarf galaxies, in reality, are not necessarily smooth and constant. They do test an alternative method using circular radii for the extended density component, and we find these density profiles are very similar. We also assume a constant ellipticity and position angle here.



# Radial velocity modeling

For both Sculptor and Ursa Minor, we construct literature samples of radial velocity measurements. We combine these samples with J+24's members to produce RV consistent stars and to compute velocity dispersion, systematic velocities, and test for the appearance of velocity gradients. 

First, we crossmatch all catalogues to J+24 Gaia stars. If a study did not report GaiaDR3 source ID's, we match to the nearest star within 2 arcseconds. We exclude stars not matched to Gaia for simplicity.

We combine the mean RV measurement from each study using the inverse-variance weighted mean and standard uncertainty. 
$$
\bar v = \frac{1}{\sum w_i}\sum_i w_i\ v_{i} 
$$
$$
\delta v = \sqrt{\frac{1}{\sum_i w_{i}^2}}
$$

where $w_i = 1/s_i^2$, and we estimate the inter-study standard deviation
$$
s^2 = \frac{1}{\sum w_i} \sum_i w_i (v_{i} - \bar v)^2
$$
We remove stars with significant velocity dispersions as measured between or within a study:
$$
P\left(\chi2(N-1) < \frac{s^2}{\delta v^2}\right) > 0.001
$$
where $s, \delta v, N$ are the weighted standard deviation, weighted standard error, and number of observations. We apply this cut to both within a study and between studies, which typically removes stars with reduced $\chi^2 \gtrsim 2$. 

The combined RV likelihood is then
$$
{\cal L} = {\cal L}_{\rm space} {\cal L}_{\rm CMD} {\cal L}_{\rm PM} {\cal L}_{\rm RV}
$$
where
$$
{\cal L}_{\rm RV, sat} = N(\mu_{v}, \sigma_{v}^2 + (\delta v_i)^2)
$$
$$
{\cal L}_{\rm RV, bg} = N(0, \sigma_{\rm halo}^2)
$$

where $\mu_v$ and $\sigma_v$ are the systemic velocity and dispersion of the satellite, and $\delta v_i$ is the individual measurement uncertainty. Typically, the velocity dispersion will dominate the uncertainty budget here. We assume a halo/background velocity dispersion of a constant $\sigma_{\rm halo} = 100$  km/s [e.g. @brown+2010].

Similar to above, we retain stars with the resulting membership probability of greater than 0.2.

Finally, we need to correct the coordinate frames for the solar motion and on-sky size of the galaxy. The first step is to subtract out the solar motion from each radial velocity, corresponding to a typical gradient of ~3 km/s across the field. The next step is to account for the slight differences in the direction of each radial velocity. Define the $\hat z$ direction to point parallel to the direction from the sun to the centre of the dwarf galaxy. Then if $\phi$ is the angular distance between the centre of the galaxy and the individual star, the corrected radial velocity is then
$$
v_z = v_{\rm los, gsr}\cos\phi  - v_{\alpha}\cos\theta \sin\phi - v_\delta \sin\theta\sin\phi
$$
where $v_{\rm tan, R} = d(\mu_{\alpha*}\cos\theta + \mu_\delta \sin \theta)$ is the radial component of the proper motion with respect to the centre of the galaxy.This correction is of the order $v_{\rm tan}\theta$ so induces a gradient of about $1 km/s/degree$ for sculptor [see @WMO2008]. The uncertainty is then the velocity uncertainty plus the distance uncertainties times the PM uncertainty from above. We then use the $v_z$ values for the following modelling, however repeating with uncorrected, heliocentric velocities does not significantly affect the results . 

For the priors on the satellite velocity dispersion and systematic velocity, we use
$$
\mu_{v} = N(0, \sigma_{\rm halo}^2) \\
\sigma_{v} = U(0, 20\,{\rm km\,s^{-1}})
$$
where $\sigma_{\rm halo} = 100\,{\rm km\,s^{-1}}$ is the velocity dispersion of the MW halo adopted above, a reasonable assumption for dwarfs in orbit around the MW. 

## Sculptor



**Figure:** Velocity histogram of Scl and UMi.



![Scl velocity sample](figures/scl_rv_2dhist.pdf)

Figure: A plot of the corrected los velocities for Scl binned in tangent plane coordinates. We detect a slight rotational gradient towards the bottom right. **TODO**: scatter plot and gradient here



![Scl velocity gradient](figures/scl_vel_gradient_binned.pdf) 

Figure: A velocity gradient in Sculptor! The arrow marks the gradient induced by Scl's proper motion on the sky. **Todo: add running median and scatter points...**



For Sculptor, we combine radial velocity measurements from APOGEE, @sestito+2023a, @tolstoy+2023, and @WMO2009. @tolstoy+2023 and @WMO2009 provide the bulk of the measurements. We find that there is no significant velocity shift in crossmatched stars between catalogues. After crossmatching to Gaia and excluding significant inter-study dispersions, we have a sample of XXXX members.

We additionally add a velocity gradient for Sculptor (adding parameters A, B $\sim N(0, 6)$km/s/deg). We derive a systemic velocity of $111.2\pm0.2$ km/s with velocity dispersion $9.67\pm0.16$ km/s. Our values are very consistent with previous work [e.g. @WMO2009, @arroyo-polonio+2024, @tolstoy+2023]. See Appendix for a more detailed comparison and inter-study tests.



We detect a gradient of $4.8\pm1.3$ km/s/deg  at a position angle of $-147_{-12}^{+15}$ degrees. 

Compared to past work, @battaglia+2008, @arroyo-polonio+2021. 

## Ursa Minor

![UMi velocity sample](figures/umi_rv_hist.pdf)

Figure: UMi velocities. More or less gaussian.

For UMi, we collect radial velocities from, APOGEE, @sestito+2023b, @pace+2020, and @spencer+2018.

We shifted the velocities of @spencer+2018 ($-1.1$ km/s) and @pace+2020 ($+1.1$ km/s) to reach the same scale. We found 183 crossmatched common stars (passing 3$\sigma$ RV cut, velocity dispersion cut, and PSAT J+24 > 0.2 w/o velocities). Since the median difference in velocities in this crossmatch is about 2.2 km/s, we adopt 1 km/s as the approximate systematic error here. 



We detect a mean $-245.9\pm0.3$  and velocity dispersion of $8.76\pm0.24$.

## Discussion and limitations

Our model here is relatively simple. Some things which we note as systematics and are challenging to account fully for are

- Inter-study systematics and biases. While basic crossmatches and a simple velocity shift, combining data from multiple instruments is challenging. 
- Inappropriate uncertainty reporting. Inspection of the variances compared to the standard deviations within a study seems to imply that errors are accurately reported. APOGEE notes that their RV uncertainties are known to be underestimates but are a small proportion of our sample.
- Binarity. While not too large of a change for classical dwarfs, this could inflate velocity dispersions of about 9 km/s by about 1 km/s @spencer+2017. Thus, our measurement is likely slightly inflated given the high binarity fractions measured in classical dwarfs [@arroyo-polonio+2023, @spencer+2018]. 
- Selection effects. RV studies each have their own selection effects. I do not know how to correct for this.

Because the derived parameters are similar for the two different larger surveys we consider for UMi and Scl, we note that many of these effects are likely not too significant (with the exception of the systemic motion of UMi.)

# Comparison and conclusions



To illustrate the differences between each dwarf galaxy, in @fig:classical_dwarfs_densities, we compare Scl, UMi, and Fnx against exponential and plummer density profiles (**TODO: state these somewhere**). While all dwarfs have marginal differences in the inner regions, each dwarf diverges in the outer regions relative to an exponential. In particular, while Fnx is underdense, Scl and UMi are both overdense, approximately fitting a Plummer density profile instead. 

In summary, we have used J+24 data to derive the density profiles for Fornax, Sculptor, and Ursa Minor. In each case, the density profile is robust against different selection criteria. Both Sculptor (Ursa Minor) show strong (weak) evidence for deviations from an exponential profile. We also compile velocity measurements to derive the systemic motions and velocity dispersions of each galaxy. We find evidence for a velocity gradient in Scl of XXX. We find no evidence of additional (velocity or stellar) substructure in either galaxy. Our goal is thus to explain why Scl and UMi have an excess of stars in their outer regions and why Scl may have a velocity gradient.



![Classical dwarf density profiles](figures/scl_umi_fornax_exp_fit.png){#fig:classical_dwarfs_densities}

Figure: The density profiles of Sculptor, Ursa Minor, and Fornax compared to Exp2D and Plummer density profiles. Dwarf galaxies are scaled to the same half-light radius and density at half-light radius (fit from the inner 3 scale radii exponential recursively. )



## 

