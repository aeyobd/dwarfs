# The reliability of derived density profiles {#sec:extra_density}

In this section, we detail the J+24 membership algorithm and test the resulting density profiles for Scl and UMi. We compare profiles derived with alternative methodologies and literature results. We also present a non-parametric Bayesian density profile, which reproduces J+24's results in the inner regions but diverges when the background dominates. In summary, we find the density profiles are robust up to the background-limited radius derived here.

## Bayesian membership probabilities

To create a high-quality sample, J+24 select stars initially from Gaia within a 2--4 degree region around each satellite with high-quality astrometry, reliable photometry, and consistent parallaxes and broadly consistent proper motions and colours. Stars are removed if they have excess astrometric noise [${\rm ruwe} > 1.3$, see @lindegren+2021], colour excess $3\sigma$ outside of expectations [from @riello+2021], proper motions magnitudes $>10\,\masyr$ in $\mu_{\alpha*}$ or $\mu_\delta$, magnitudes brighter than the tip of the red giant branch or fainter than $G=22$, or colours outside $-0.5 < G_{\rm BP} - G_{\rm RP} <  2.5$. Photometry is dereddened with @schlegel+finkbeiner+davis1998 extinction maps.

J+24 construct satellite likelihoods in spatial, proper motion, and CMD space following expected satellite properties.  J+24 model the spatial likelihood as either a one or two-component exponential (@eq:exponential_law). The structural parameters of the inner component are fixed, and marginalized over if one-component. The outer profile scale radius and normalization are free parameters. 

The PM likelihood is a bivariate gaussian with variance and covariance equal to each star's proper motions. 

J+24 model the CMD as a Padova isochrone [@girardi+2002], with matching metallicity, 12 Gyr age (or 2 Gyr for Fornax), and Gaussian width of 0.1 mag plus the *Gaia* colour uncertainty at each magnitude. The horizontal branch is modelled as a constant-magnitude^[Specifically, the mean magnitude of the 12 Gyr, ${\rm [Fe/H]}=-2.2$, Padova isochrone's horizontal branch.] sequence extending blue of the CMD with the same width as the RGB. The CMD likelihoods are marginalized over the distance modulus and take the maximum of the RGB and horizontal branch likelihoods.

The background likelihoods are determined empirically as a kernel density estimates from stars outside 5$R_h$ in PM and CMD space. The spatial background likelihood is uniform.  

J+24 derive the distributions of parameters (proper motions, satellite fraction $f_{\rm sat}$, and second spatial component if included) through Monte Carlo Markov chain sampling with broad or weakly-informative priors. The posterior modes are used to calculate the final $P_{\rm sat}$ values. 

## Possible biases in *Gaia*-derived density profiles {#sec:density_extra}

In this section, we test how assumptions in J+24's algorithm or biases in *Gaia* may affect density profiles. In all cases, density profiles converge until $R_{\rm ell} \approx 60\,{\rm arcmin}$, where background contamination likely dominates. 

**Spatial likelihood.** As J+24 include a spatial likelihood, our density profiles risks "double fitting"---the spatial distribution of stars informs the membership catalogue from which density profiles are derived. If the assumed spatial likelihood is misrepresentative, the derived density profiles may bias towards the assumed likelihood. @fig:scl_umi_density_extras illustrates that selecting the one- or two-component spatial models affects the outer extent of the dwarf galaxy. Promisingly, the density excess persists when assuming a one-component exponential profile. We revisit this concern using spatially independent methods below. 

**Structural uncertainties**. The assumed structural parameters of a dwarf (centre, position angle, ellipticity, scale radius) may be biased or vary radially. As a test, we use a membership sample J+24 derive assuming the outer component is circular. Even when binning in circular radius, the density profiles remain unchanged (the `circ` points in @fig:scl_umi_density_extras). Thus, our conclusions appear to be robust to an extreme shift in ellipticity. Reasonable changes in other structural parameters are unlikely to have a stronger effects.

**Completeness**. *Gaia* appears to have high but imperfect completeness, more limited for crowded fields, faint sources ($G\gtrsim20$), and BP-RP magnitudes. In @fabricius+2021, the completeness down to $G\approx 20$ is $\sim 80\%$ for low-density globular clusters. As dwarfs are even less dense, their completeness is likely higher. To test for magnitude-dependent biases, we show density profiles derived for only the brighter half of member stars in @fig:scl_umi_density_extras. We find no substantive change compared to our fiducial density profile. If *Gaia*'s incompleteness is inhomogeneous, then these variations likely do not affect our density profiles.





![Sculptor density methodology comparison](figures/density_methods_extra.pdf){#fig:scl_umi_density_extras}

Figure: Density profiles for various assumptions for Sculptor. `2-exp` is the fiducial double-exponential-likelihood J+24 sample, `1-exp` instead is a one-component exponential spatial likelihood, `simple` uses the simple position-independent selection criteria,  `circ` is a 2-component bayesian model assuming circular radii, `bright` only includes stars brighter than the median magnitude, and `DELVE` or `UNIONS` use photometry from external surveys with background subtraction. In all cases, the density profiles are nearly identical until the background-limited regime (grey shaded region). 



## Comparison against alternate samples {#sec:simple_selection}

Next, to compare against simpler methods, we create samples using absolute filters in the CMD, PM space. We also compare against similarly selected samples from deep photometry-only data. 

For the `simple` *Gaia* samples, we require high-quality astrometry (`ruwe` < 1.3), parallaxes $3\sigma$-consistent with the dwarf's distance,  proper motions with $1\masyr$ radius of the galaxy's value, and within the dwarf's empirical CMD (@fig:extra_cmd).

We also determine density profiles from the deeper photometric surveys. For Scl, we use DELVE DR2 survey [@drlica-wagner+2022]. We select sources within an ellipse of radius 150 arcminutes, categorized as likely stars, with reliable $g$ and $r$ magnitudes (associated flags $\leq 4$), and within Scl's observed CMD (@fig:extra_cmd). Ursa Minor is within the Ultraviolet Near-Infrared Optical Northern Survey [UNIONS, @gwyn+2025]. We select sources within an ellipse of radius 230 arminutes, with no `FLAGS_CFIS` set, the `s21` and `s31` both $<3$ (i.e. not extended sources), and within the CMD cut in @fig:extra_cmd. @fig:delve_unions_tangent shows the resulting spatial distribution for both samples. The only remarkable structure is Muñoz 1 nearby to Ursa Minor. 

Finally, @fig:scl_umi_density_extras compares the density profiles from each sample. All samples agree until the "limiting radius", where the `simple` and DELVE/UNIONS samples reach their respective backgrounds. The very outskirts of J+24's density profiles may also extend below the background of apparent member stars, complicating the reliability of the outer density profiles. 



![Colour-Magnitude sample selection](figures/extra_cmd_selection.png){#fig:extra_cmd width=80%}

Figure: Colour-magnitude diagram cuts used in samples in this section. The likely-member stars in the inner $1R_h$ of each dwarf (black points) are used as a guide to determine the CMD cuts (orange polygons). Light grey points instead show the full distribution of stars across the field.





![DELVE and UNIONS spatial distribution of stars](figures/delve_unions_tangent.pdf){#fig:delve_unions_tangent}

Figure: The distribution of DELVE and UNIONS selected stars in Scl and UMi as grey points overdrawn with isophotes. The clump to the lower right of UMi is Muñoz 1. 



## A Bayesian density profile {#sec:mcmc_hists}

To address the concerns discussed above, we consider here a non-parametric model to fit the density in each bin. We demonstrate that the non-parametric model fails to find evidence for stellar density where J+24 still detects members. Notably, our density profiles diverge substantially from J+24 for Antlia II. We suggest these differences arise from background contamination, motivating our "limiting radii" representing where densities may become unreliable.



### Methodology

As a non-parametric but similar model to J+24, we consider a piecewise constant spatial likelihood. In this model, the stars are divided into radial bins which are then considered independently. Each bin has only one free parameter, the fraction of satellite stars in that bin, $f_{\rm sat}$. This model thus directly derives the density profile from the data in a single step. If a bin contains insufficient information to estimate a precise satellite density, than the posterior $f_{\rm sat}$ should reflect this uncertainty. 

We parameterize $f_{\rm sat}$ in terms of a log-relative density, $\theta$
$$
\begin{split}
\theta &\equiv \log_{10}({\Sigma_{\rm sat}}/{\Sigma_{\rm bg}}) \\
f_{\rm sat} &= \frac{10^{\theta}}{1 + 10^{\theta}}.
\end{split}
$$
We then adopt a broad uniform prior on $\theta$ from -12 to 6. The CMD and PM likelihoods are unchanged from J+24, except the systemic PM is fixed for efficiency. 

To bin stars, we hold fixed J+24's structural parameters and create bins of the wider of the width $\Delta \log R=0.05$ or containing the next 20 stars. We then derive posterior $f_{\rm sat}$ distributions with MCMC (48 walkers of 1000 steps each, and using the No-U-Turn sampler as implemented in Turing.jl). The final density profile is directly derived from $f_{\rm sat}$ and the number of stars in each bin.

### Results



| Galaxy           | $R_{\rm limit} / R_h$ | $R_{\rm limit} / '$ |
| ---------------- | --------------------- | ------------------- |
| Fornax           | 5.25                  | 79.1                |
| Sculptor         | 6.39                  | 64.1                |
| Leo I            | 4.23                  | 13.5                |
| Ursa Minor       | 6.42                  | 86.4                |
| Leo II           | 3.63                  | 8.76                |
| Carina           | 4.16                  | 33.3                |
| Draco            | 3.59                  | 27.3                |
| Canes Venatici I | 1.95                  | 12.5                |
| Sextans I        | 3.42                  | 67.9                |
| Crater II        | 1.93                  | 39.0                |

Table: For each classical dwarf, the limiting radius $R_{\rm limit}$ in units of $R_h$ and arcminutes. $R_{\rm limit}$ represents where there no longer appears to be evidence of stars in *Gaia* using the nonparametric MCMC density profiles. {#tbl:mcmc_props short="The limiting radii of Gaia-derived density profiles"}

[@fig:mcmc_hists; @fig:mcmc_hists2] show the derived, MCMC density profiles as compared to J+24. In general, both methodologies are consistent. However, J+24 tend to systematically overestimate faint densities and confidently derive densities where the MCMC model fails to derive a density estimate. 

In the outskirts of satellites, more background / foreground stars may have consistent PM and CMD properties than satellite members. Improperly estimating the satellite's density in this "background-limited" regime likely affects the inferred density of members. If the satellite's density is severely overestimated, then a J+24-like sample may select many additional background stars. The density of candidates would then be biased towards the assumed density. As J+24 only use a one or two-component density profile across the entire dwarf, the density profile is likely biased where the CMD+PM background dominates. The nonparametric MCMC model instead does not assume a local density, so should better represent the underlying satellite density.

Antlia II represents an extreme example---J+24 and the non-parametric density model systematically disagree across the entire galaxy. Additionally, the piecewise model derives densities with larger uncertainties and over a smaller range than J+24. The methodological divergence here likely arises due to an extraordinary background/foreground of MW stars. *Gaia* data alone may be unable to properly constrain the density of such background-contaminated objects. 

To properly compare density profiles before background-limiting effects become important, we only calculate our profiles in the main text out to the radii in @tbl:mcmc_props. We derive these "background-limited" radii based on the outermost derived density in the MCMC non-parametric model with an uncertainty lower than 1 dex. This closely corresponds to the background density from the "CMD+PM" samples in the main text (see @fig:scl_observed_profiles).



![Probabilistic density profiles](figures/mcmc_histograms.png){#fig:mcmc_hists width=100%}

Figure: A comparison between the MCMC histogram method and J+24. The MCMC samples in each bin are black transparent dots (with added jitter), and the J+24 derived density profiles (i.e., $P_{\rm sat} < 0.2$) are orange solid dots.



![Probabilistic density profiles continued](figures/mcmc_histograms2.png){#fig:mcmc_hists2 width=100%}

Figure: @fig:mcmc_hists continued. 



## Comparison to literature

Finally, we compare our density profiles against a sample of literature-derived density profiles in [@fig:scl_lit_profiles; @fig:umi_lit_profiles]. Deviations between different profiles are slight, despite the range of methods across time, and all density profiles extending into the outskirts of Scl and UMi show a similar overdensity to the ones we find. The extended stellar profiles of Scl and UMi appear to be a robust result across the literature. 



![Sculptor literature density profiles](figures/scl_literatre_profiles.pdf){#fig:scl_lit_profiles}

Figure: A comparison of the Scl derived density profile and historical works. Unlike most profiles in this thesis, this density profile is plotted with respect to the semi-major elliptical radius ($a = R_{\rm ell} / \sqrt{1-{\rm ell}}$). The solid black line is a 2D exponential with corresponding scale radius to Scl's half-light radius, and the residuals in the bottom panel are with respect to this profile. References are (in order), @munoz+2018; @westfall+2006; @walcher+2003; @eskridge1988; @demers+krautter+kunkel1980; and @hodge1961. 

![Ursa Minor literature density profiles](figures/umi_literature_profiles.pdf){#fig:umi_lit_profiles}

Figure: Similar to @fig:scl_lit_profiles except for Ursa Minor. The references are (in order) @sato+2025 [derived using their minor axis profile]; @palma+2003; @martinez-delgado+2001; @kleyna+1998; @IH1995; and @Hodge1964. 





## Summary

In this Appendix, we discussed the methodological details of the J+24 sample selection algorithm. We consider possible biases due to the assumed spatial likelihood, *Gaia*'s completeness, and structural parameters. In all cases, we find these assumptions likely do not cause major biases in the derived density profiles. However, we show that J+24's density profiles may become unreliable when dropping below the background of satellite-like stars. Finally, comparing our density profiles against the literature, we find our results to be consistent. We conclude that the detection of an extended density profile in Scl and UMi is robust to incompleteness, methodology, alternative surveys, and across the literature.
