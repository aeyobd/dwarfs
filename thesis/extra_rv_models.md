# Line-of-sight velocities: Sample selection and modelling  {#sec:extra_rv_models} 

In this section, we analyze the observed line-of-sight (LOS) velocity distributions for Sculptor (Scl) and Ursa Minor (UMi). Our aim is to test for kinematic tidal signatures by combining literature LOS velocities with J+24's membership framework. Our derived systemic velocities and dispersions are consistent with literature. We find weak evidence for a velocity gradient in Scl. As the gradient is misaligned with the orbit, the gradient more likely reflects intrinsic rotation than tidal disruption.  Scl's velocity dispersion also rises within the effective radius, but predominantly in the inner regions, contrary to tidal disruption. We find no evidence of a gradient in mean velocity or velocity dispersion in UMi. We conclude that Scl and UMi do not show clear features of tidal disruption given current velocity observations. 

## Data processing  and selection

We compile LOS velocities from several spectroscopic surveys. For Sculptor, we combine @tolstoy+2023; @walker+2009; @sestito+2023a; and \apogee{} [DR17, @abdurrouf+2022]. For Ursa Minor, we combine @spencer+2018; @pace+2020; @sestito+2023b; and \apogee{}. We then crossmatch all catalogues to J+24 Gaia stars. If a study did not report GaiaDR3 source ID's, we match to the nearest star within 1--3 arcseconds. We combine measurements of the same star using inverse-variance weighting. 

To reduce likely binaries, we remove stars with significant velocity dispersions.^[Specifically, using that $\chi^2=\frac{s^2}{\delta \bar v^2}$, we remove stars with a $\chi^2$ larger than the 99.9th percentile of the $\chi^2$ distribution with $N-1$ measurements.]

We build on J+24's likelihood by adding multiplicative terms in the total likelihood for the velocity consistency (see @sec:extra_density). We assume that the satellite and background $v_{\rm los, gsr}$ distributions are Gaussian.  We assume a halo/background velocity dispersion of a constant $\sigma_{\rm halo} = 100\,\kms$  with mean 0 [e.g. @brown+2010]. We select stars with velocity-informed satellite membership probabilities of greater than 0.2. For Scl, we find 1918 unique members and UMi, 831. 

For UMi, we shifted the velocities of @spencer+2018 ($-1.1\,\kms$) and @pace+2020 ($+1.1\,\kms$ ) to account for a systematic velocity offset. Otherwise all studies appear to agree.

We correct the velocities for the solar motion and on-sky size of the galaxy. We transform the velocities into the galactic standard of rest (GSR, i.e. same location as ICRS but velocities relative to the galactic centre), and then transform velocities to correct for the apparent gradient induced by the dwarf's proper motion.^[Specifically, sub] We define $v_{\rm gsr}'$ to be velocities in the The correction from both effects induces an apparent gradient of about $1.3\,\kmsdeg$ for Sculptor and less for Ursa Minor [see also @WMO2008; @strigari2010]. The uncertainty on this velocity correction are less than the individual star uncertainties and our derived velocity gradient.

## MCMC modeling

We fit Monte-Carlo Markov chain (MCMC) models to solve for the systemic velocity, velocity dispersion, and possible gradient. We assume that the galaxy follows a planar velocity gradient. The mean velocity at a given point on the sky, $\mu$, is assumed to be,

$$
\mu(\xi, \eta) = \mu_0 + a\,\xi + b\,\eta
$$

for tangent-plane coordinates $\xi$ and $\eta$, systemic velocity $\mu_0$, and velocity gradient slopes $a$ and $b$. The velocity dispersion at a given position, $\sigma$, is assumed to depend as a power-law on elliptical radius $R_{\rm ell}$ alone:
$$
\log \sigma = \log \sigma_0 + c\,\log(R_{\rm ell} / R_h)
$$
where $\sigma_0$ is the system's velocity dispersion at $R_h$, and $c$ is the velocity dispersion gradient slope. We use weakly-informative priors, as described in @tbl:scl_rv_mcmc, and with $a, b \sim N(0, 6^2)\,\kmsdeg$. 

## Results {#sec:rv_results}



![LOS velocity fit to Scl](figures/scl_umi_rv_fits.pdf)

Figure: Velocity histogram of Scl and UMi in terms of projected-correction GSR LOS velocity ([@eq:v_z]). Black points with error bars are from the crossmatched observed sample, and green lines represent MCMC fits to the velocity dispersion.

Figure: The observed velocity dispersion gradient (black) in 10 equal number bins and the derived slopes from the model fitting. Scl shows moderate evidence for an increasing velocity dispersion with radius. UMi's evidence for a radial gradient is weaker, and Scl's velocity dispersion appears to be flat outside of $R_h$. 

\input{rv_table.tex}



![Scl velocity gradient](figures/scl_rv_scatter_gradient.png){#fig:scl_velocity_gradient_scatter} 

Figure: **Top**  members of Sculptor plotted in the tangent plane coloured by corrected velocity difference from mean $v_z - \bar v_z$ . The black ellipse marks the half-light radius in @fig:scl_selection. The black and green arrows mark the proper motion (PM, GSR frame) and derived velocity gradient (rot) vectors (to scale). **Bottom**: The corrected LOS velocity along the best fit rotational axis. RV members are black points, the systematic $v_z$ is the horizontal grey line, blue lines represent the (projected) gradient from MCMC samples, and the orange line is a rolling median (with a window size of 50).



We derive a systemic velocity for Sculptor of $111.3\pm0.2\,\kms$with velocity dispersion $9.64\pm0.16\,\kms$. Our values are very consistent with previous work [e.g. @walker+2009, @arroyo-polonio+2024, @battaglia+2008].  We detect a moderately significant gradient of $4.3\pm1.3\,\kmsdeg$   at a position angle of $-149_{-13}^{+17}$ degrees,  similar to past detections [e.g., @arroyo-polonio+2024; @battaglia+2008; but see also @strigari2010; @martinez-garcia+2023].

@fig:scl_velocity_gradient_scatter plots the combined velocity sample and the MCMC samples for the velocity gradient. While some samples have a consistent gradient with 0, most samples have a positive velocity gradient, consistent with the rolling median trend. The velocity gradient in Scl is misaligned with the proper motion

We derive a mean $-245.8\pm0.3_{\rm stat}\,\kms$  and velocity dispersion of $8.8\pm0.2\,\kms$ for UMi. We do not find evidence for a velocity gradient.   This is consistent with previous work [@pace+2020; somewhat with @spencer+2018; @martinez-garcia+2023].

## Discussion and caveats

*Binarity*. Stars in binary systems can inflate the inferred system's velocity dispersion. For Scl and UMi with high measured velocity dispersion, binaries likely add $\sim 1\,\kms$ to the velocity dispersion [@spencer+2018; @gration+2025]. 

*Multiple populations*. Both Sculptor and Ursa Minor likely contain multiple populations with different kinematics[@arroyo-polonio+2024, @pace+2020, @tolstoy+2004]. We do not model these separately---it is unclear how to define a velocity dispersion in such a case.

*Inter-study biases*. While basic crossmatches and a simple velocity shift, combining data from multiple instruments is challenging. Studies may also have biased selection effects or misrepresentative uncertainties, biasing our results. However, besides the $\sim1\,\kms$ velocity offset in Ursa Minor, our results are consistent across samples. 

*Tangental motions.* Precise star-by-star proper motions would further test for kinematic disequilibrium. Presently, *Gaia* uncertainties are too large to permit such studies for Sculptor and Ursa Minor.^[Specifically, typical proper motion uncertntanties for faint, member stars of Scl and UMi are around $0.5\,\masyr$, corresponding to velocities $\sim 100\,\kms$.]

## Summary

In this section, we analyze literature samples of LOS velocity measurements for both Scl and UMi. In each case, we find systemic LOS velocities and velocity dispersions consistent with past work. We detect weak evidence for a velocity gradient in Scl, however the gradient is mis-aligned with the proper motion of Scl so is unlikely to be of tidal origin. In both galaxies, the velocity dispersion profile is consistent with a constant velocity dispersion with radius. The lack of observational evidence for ongoing tidal disruption in the velocity distribution of stars in Scl and UMi further supports our interpretations in the main text.  
