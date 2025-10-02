# The reliability of derived density profiles {#sec:extra_density}

In this section, we describe the details of J+24's likelihoods. We then test the reliability of the derived density profiles for Scl and UMi, comparing methodologies, samples, and to literature. Finally, we present a comparison of J+24's samples to a non-parametric Bayesian-derived density profile, finding similar results except in the very outer regions. Past a "limiting radius", the density profiles become susceptible to the likelihood specification and may select mostly background stars. In conclusion, the density profiles from J+24 converge among methods out to our limiting radii derived here. 

## Bayesian membership probabilities

To create a high-quality sample, J+24 select stars initially from Gaia within a 2--4 degree circular region centred on the dwarf satisfying: 

- Solved astrometry, magnitude, and colour. 
- Renormalized unit weight error, ${\rm ruwe} \leq 1.3$, ensuring high quality astrometry. `ruwe` is a measure of the excess astrometric noise on fitting a consistent parallax-proper motion solution [see @lindegren+2021].  
- 3$\sigma$ consistency of measured parallax with dwarf's distance, with @lindegren+2021's zero-point correction. Note that the typical satellite parallax is very small. 
- Absolute proper motions, $\mu_{\alpha*}$, $\mu_\delta$, less than 10$\,{\rm mas\ yr^{-1}}$. This corresponds to tangental velocities of $\gtrsim 500$ km/s at distances larger than 10 kpc.
- Corrected colour excess is within 3$\sigma$ of the expected distribution from @riello+2021. Removes stars with unreliable photometry. 
- De-reddened $G$ magnitude is between $22 > G > G_{\rm TRGB} - 5\sigma_{\rm DM}$. Removes very faint stars and stars significantly brighter than the tip of the red giant branch (TRGB) magnitude plus the distance modulus uncertainty $\sigma_{\rm DM}$. 
- De-reddened colour is between $-0.5 < {\rm BP - RP} <  2.5$. Removes stars substantially outside the expected CMD.

Photometry is dereddened with @schlegel+finkbeiner+davis1998 extinction maps.

J+24 define likelihoods ${\cal L}$ measuring consistency with the MW stellar background (${\cal L}_{\rm bg}$) or the satellite galaxy (${\cal L}_{\rm sat}$). The total likelihoods are the product of a spatial, PM, and CMD term:
$$
{\cal L} = {\cal L}_{\rm space}\ {\cal L}_{\rm PM}\ {\cal L}_{\rm CMD}.
$$

Each likelihood is normalized over their respective 2D parameter space for both the satellite. To control the relative frequency of member and background stars, $f_{\rm sat}$ representing the fraction of member stars in the field. The total likelihood for any star in this model is the $f_{\rm sat}$-weighted sum of the foreground and background, 
$$
{\cal L}_{\rm tot} = f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}.
$$ {#eq:Ltot}
Finally, the probability that any star belongs to the satellite is then given by 
$$
P_{\rm sat} = 
\frac{f_{\rm sat}\,{\cal L}_{\rm sat}}{{\cal L}_{\rm tot}}
= \frac{f_{\rm sat}{\cal L}_{\rm sat}}{f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}}.
$$



For the satellite's spatial likelihood, J+24 consider both one-component and a two-component density models. The one component model is constructed as a single exponential profile  ( surface density $\Sigma \propto e^{R_{\rm ell} / R_s}$), with scale radius $R_s$. $R_s$ is fixed to the equivalent $R_h$ value from @munoz+2018's Sérsic fit. Additionally, structural uncertainties (for position angle, ellipticity, and scale radius) are sampled over to construct the final likelihood map. The two-component model instead adds a second exponential, $\Sigma_\star \propto e^{-R/R_s} + B\,e^{-R/R_{\rm outer}}$. The inner scale radius is fixed, and the outer scale radius and magnitude of the second component $R_{\rm outer}$, $B$  are free parameters. Structural property uncertainties are not included in the two-component model.

The PM likelihood is a bivariate gaussian with variance and covariance equal to each star's proper motions. J+24 assume the stellar PM errors are the main source of uncertainty.

The satellite's CMD likelihood is based on a Padova isochrone [@girardi+2002]. The isochrone has a matching metallicity and 12 Gyr age (except 2 Gyr is used for Fornax). The isochrone is given a Gaussian colour width of 0.1 mag plus the *Gaia* colour uncertainty at a given magnitude. The horizontal branch is modelled as a constant magnitude extending blue of the CMD (mean magnitude of -2.2, 12 Gyr HB stars and a 0.1 mag width plus the mean colour error). A likelihood map is constructed by sampling the distance modulus in addition to the CMD width, taking the maximum of RGB and HB likelihoods.

The background likelihoods are instead empirically constructed. Stars stars outside of 5$R_h$ passing the quality cuts estimate the background density in PM and CMD space. The density is a sum of bivariate gaussians with variances based on Gaia uncertainties (and covariance for proper motions).  The spatial background likelihood is assumed to be constant. 

J+24 derive $\mu_{\alpha*}$, $\mu_\delta$, $f_{\rm sat}$ (and $B$, $R_{\rm outer}$ for two-component) through an MCMC simulation with likelihood from [@Eq:Ltot]. Priors are only weakly informative. The proper motion single component prior is same as @MV2020a: a normal distribution with mean 0 and standard deviation $100\ \kms$. If using the 2-component model, the prior is instead a uniform distribution spanning 5$\sigma$ of single component case w/ systematic uncertainties. $f_{\rm sat}$  (and $B$) has a uniform prior from 0 to 1. $R_{\rm outer}$ has a uniform prior only restricting $R_{\rm outer} > R_s$. The mode of each parameter from the MCMC are then reported and used to calculate the final $P_{\rm sat}$ values. 

## Alternative methodologies

### Simplistic cuts

Besides a global background subtraction (as used in @sec:observations), absolute cuts in CMD and PM represent perhaps the next simplest method to construct density profiles. We create cuts using only the following

- (BP-RP, G) within the vertices defined in @tbl:colour_cuts
- PM within a $1\masyr$ radius of the galaxies measured value. 
- Parallax within $3\sigma$ of the dwarf distance
- `ruwe` < 1.3

| Scl   |       |       |       | UMi   |       |        |       |
| ----- | ----- | ----- | ----- | ----- | ----- | ------ | ----- |
| Gaia  |       | DELVE |       | Gaia  |       | Unions |       |
| BP-RP | G     | g-r   | g     | BP-RP | G     | r-i    | r     |
| -0.21 | 20.58 | 0.32  | 22.58 | 0.49  | 20.83 | -0.485 | 24.13 |
| -0.1  | 20.1  | 0.51  | 23.07 | 1.24  | 20.9  | 0.145  | 23.11 |
| 0.28  | 19.77 | 0.56  | 21.49 | 1.18  | 19.65 | 0.201  | 21.41 |
| 0.7   | 19.64 | 0.71  | 19.7  | 1.23  | 18.12 | 0.21   | 20.05 |
| 1.02  | 18.5  | 0.96  | 18.49 | 1.39  | 17.09 | 0.288  | 17.23 |
| 1.14  | 17.59 | 0.84  | 18.35 | 1.32  | 16.94 | 0.226  | 17.21 |
| 1.41  | 16.73 | 0.6   | 19.65 | 0.96  | 18.58 | 0.151  | 19.3  |
| 1.68  | 15.94 | 0.42  | 21.42 | 0.88  | 19.22 | 0.145  | 19.75 |
| 1.94  | 15.94 |       |       | 0.4   | 19.52 | 0.026  | 19.86 |
| 1.91  | 16.8  |       |       | -0.03 | 19.78 | 0.011  | 19.52 |
| 1.68  | 17.01 |       |       | 0.06  | 20.16 | -0.155 | 19.95 |
| 1.45  | 17.87 |       |       | 0.83  | 19.72 | -0.22  | 20.34 |
| 1.32  | 19.16 |       |       | 0.21  | 20.83 | -0.176 | 20.47 |
| 1.26  | 20.1  |       |       |       |       | -0.098 | 20.33 |
| 1.35  | 21    |       |       |       |       | 0.008  | 19.84 |
| 0.46  | 21    |       |       |       |       | 0.142  | 19.75 |
| 0.82  | 20.27 |       |       |       |       | 0.098  | 21.45 |
| 0.39  | 20.33 |       |       |       |       | -0.052 | 22.38 |
| 0.16  | 20.76 |       |       |       |       | -0.388 | 23.42 |

Table: Colour cuts used in samples in this section. Each pair of columns contains the vertices of a polygon in the CMD within which we select stars to derive density profiles. {#tbl:colour_cuts short="Colour cuts for density profiles."}



### Alternative samples

For the DELVE sample [@drlica-wagner+2022], we select everything in data release 2 between RA of 11 and 19, and declination of -37.7 and -29.7.

We also consider samples from the Ultraviolet Near-Infrared Optical Northern Survey [UNIONS @gwyn+2025] for Ursa Minor. We select stars within a 5 degree tangent square of UMi with no `FLAGS_CFIS` set and `s21` and `s31` both $<3$ (not extended sources).

### A Bayesian histogram model

To address the concerns discussed above, we consider here a non-parametric model to fit the density in each bin. This model extends the J+24 framework with two notable differences. (1) The systemic proper motions are held fixed. This allows for a much more efficient calculation of the likelihoods. (2) The spatial likelihood for the satellite is a piecewise constant density function.While this model theoretically has many more parameters to fit (one for each bin), each parameter is independently estimated based on only the stars in that specific bin.

We segregate the data into radial bins and independently fit the model to only stars in each bin. This fits many small models instead of a global large model. In detail, this means that, in each bin, the normalized density of the background and satellite are *equal*. The only fit parameter is the $f_{\rm sat}$ in the particular bin. Therefore, the likelihood in each bin is the same as above except ${\cal L}_{\rm space} = 1$ for the foreground and background. 

In detail, we parameterize $f_{\rm sat}$ in terms of the log-relative satellite density. The model construction is then
$$
\theta_i \equiv \log_{10}({\Sigma_{\rm sat}}/{\Sigma_{\rm bg}})
$$

$$
f_{\rm sat} = \frac{10^{\theta_i}}{1 + 10^{\theta_i}}
$$

$$
\theta_i \sim {\rm Uniform}(-12, 6)
$$

$$
{\cal L}_{\rm tot} = f_{\rm sat}{\cal L}_{\rm CMD, PM} + (1-f_{\rm sat}){\cal L}_{\rm CMD, PM}
$$



For each satellite, we hold the structural parameters fixed to those in J+24, create bins which are the larger of the bin containing 20 stars or 0.05dex in $\log R_{\rm ell}$. We then sample each chain 1000 steps with 48 realizations using Turing.jl and the No-U Turns Sampler. The chains appear to converge excellently.

Based on the $f_{\rm sat}$ values in each bin, the number of members is just the number of stars in the bin times $f_{\rm sat}$ for the bin, allowing the direct derivation of the density profiles from this model.

We have also briefly explored a model which solves for $f_{\rm sat}$ for random realizations of the global structural uncertainties. These do not appear to influence the results substantially. Another possible inconsistency is that the half-light radii are estimated from a different sample. We do briefly rederive structural parameters for the sample using Sérsic fits (although with a fixed origin).



## Comparison of density profiles {#sec:density_extra}

In this section, we discuss additional tests and verification of the derived density profiles. In particular, we check that methodology (simpler cuts, circularized radii, algorithm version) do not substantially affect the density profile. We also compile density profiles presented in the literature as reference. In all cases, the density profiles appear to have excellent convergence out to $\log R_{\rm ell} / {\rm arcmin} \approx 1.8$, about the distance where the background dominates. 

**Spatial likelihood.** J+24 method was designed in particular to detect the presence of a density excess and find individual stars at large radii to be followed up. We are more interested in accurately quantifying the density profile and size of any perturbations. One potential problem with using J+24's candidate members is that the algorithm assumes the density is either described by a single or double exponential. If this model does not accurately match the actual density profile of the dwarf galaxy, we want to understand the impact of this assumption.In particular, in @Fig:scl_observed_profiles, notice that the $P_{\rm sat}$ selection method produces small errorbars, even when the density is more than 1 dex below the local background. These stars are likely selecting stars from the statistical MW background consistent with UMi PM+CMD, recovering the assumed density profile. As a result, the reliability of these density profiles below the CMD+PM background may be questionable. A more robust analysis, removing this particular density assumption, would be required to more appropriately represent the knowledge of the density profile as the background begins to dominate. 

***Gaia* systematics**. While *Gaia* has shown excellent performance, some notable limitations may introduce problems in our interpretation and reliability of density profiles. Gaia systematics in proper motions and parallaxes are typically smaller than the values for sources of magnitudes $G\in[18,20]$. Since we use proper motions and parallaxes as general consistency with the dwarf, and factor in systematic uncertainties in each case, these effects should not be too significant. However, the systematic proper motion uncertainties becomes the dominant source of uncertainty in the derived systemic proper motions of each galaxy (see [@tbl:scl_obs_props; @tbl:umi_obs_props]).

**Completeness**. *Gaia* shows high but imperfect completeness, particularly showing limitations in crowded fields and for faint sources ($G\gtrapprox20$). As discussed in @fabricius+2021, for the high stellar densities in globular clusters, the completeness relative to HST varies significantly with the stellar density. However, the typical stellar densities of dwarf galaxies are much lower, at about 20 stars/arcmin = 90,000 stars / degree, lower than the lowest globular cluster densities and safely below the crowding limit of 750,000 objects/degree for BP/RP photometry. In @fabricius+2021, for the lowest density globular clusters, the completeness down to $G\approx 20$ is $\sim 80\%$. Closely separated stars pose problems for Gaia's on-board processing, as the pixel size is 59x177 mas on the sky. This results in a reduction of stars separated by less than 1.5" and especially for stars separated by less than 0.6 arcseconds. The astrometric parameters of closely separated stars furthermore tends to be of lower quality [@fabricius+2021]. However, even for the denser field of Fornax, only about 3% of stars have a neighbour within 2 arc seconds, so multiplicity should not affect completeness too much (except for unresolved binaries). One potential issue is that the previous analyses do not account for our cuts on quality and number of astrometric parameters. These could worsen completeness, particularly since the BP-RP spectra are more sensitive to dense fields.

 In [@fig:scl_density_extras; @fig:umi_density_extras] we show the results of limiting the density profile to the brightest half of stars. We find no substantive differences. If *Gaia*'s incompleteness is inhomogeneous, then the variations in *Gaia*'s completeness do not strongly affect the density profiles.

**Structural uncertainties**. J+24 do not account for structural uncertainties in dwarfs for the two component case. We assume constant ellipticity and position angle. Dwarf galaxies, in reality, are not necessarily smooth or have constant ellipticity. J+24 test an alternative method using circular radii for the extended density component, and we find these density profiles are very similar to the fully elliptical case, even when assuming circular bins for the circular outer component. As such, even reducing the assumed ellipticity from $0.37-0.55$ to 0 does not substantially impact the density profiles. 



![Scl density comparison](figures/scl_density_methods_extra.pdf){#fig:scl_density_extras}

Figure: Density profiles for various assumptions for Sculptor. PSAT is our fiducial 2-component J+24 sample, circ is a 2-component bayesian model assuming circular radii, simple is the series of simple cuts described, bright is the sample of the brightest half of stars (scaled by 2), DELVE is a sample of  RGB  stars (background subtracted and rescaled to match).

Comparison of density profiles for each J+24 method. The fiducial is a 2-component elliptical model. However, the 1-component is still elliptical but only contains 1 component and the circular model assumes a circular outer density profile and bins in circular bins instead of elliptical bins. 



![Scl literature density profiles](figures/scl_literatre_profiles.pdf)

Figure: A comparison of the Scl derived density profile and historical works. Unlike most profiles in this thesis, this density profile is plotted with respect to the semi-major elliptical radius ($a = R_{\rm ell} / \sqrt{1-{\rm ell}}$). The solid black line is a 2D exponential with corresponding scale radius to Scl's half-light radius, and residuals in the lower panel are relative to this exponential. 





### Ursa Minor

![UMi density comparison](figures/umi_density_methods_extra.pdf){#fig:umi_density_extras}

Figure: Similar to [@fig:scl_observed_profiles] except for Ursa Minor.

![UMi density profiles in literatures](figures/umi_literature_profiles.pdf)

Figure: Similar to 



| Galaxy           | $R_{\rm limit} / R_h$ | $R_{\rm limit} / '$ | num stars | $n$                       |
| ---------------- | --------------------- | ------------------- | --------- | ------------------------- |
| Fornax           | 5.25                  | 79.1                | 23154     | $0.794^{+0.034}_{-0.018}$ |
| Sculptor         | 6.39                  | 64.1                | 7024      | $1.281^{+0.04}_{-0.04}$   |
| Leo I            | 4.23                  | 13.5                | 1242      | $0.719^{+0.059}_{-0.053}$ |
| Ursa Minor       | 6.42                  | 86.4                | 2314      | $1.425^{+0.10}_{-0.093}$  |
| Leo II           | 3.63                  | 8.76                | 347       | $0.707^{+0.11}_{-0.092}$  |
| Carina           | 4.16                  | 33.3                | 2389      | $1.025^{+0.08}_{-0.076}$  |
| Draco            | 3.59                  | 27.3                | 1781      | $1.106^{+0.084}_{-0.078}$ |
| Canes Venatici I | 1.95                  | 12.5                | 156       | $1.18^{+0.36}_{-0.28}$    |
| Sextans I        | 3.42                  | 67.9                | 1830      | $1.136^{+0.10}_{-0.09}$   |
| Crater II        | 1.93                  | 39.0                | 507       | $0.535^{+0.12}_{-0.09}$   |

Table: For each classical dwarf, we have: the BG-limited radius $R_{\rm limit}$, where the density of stars is no longer reliably derived, the approximate number of member stars in Gaia, and the derived Sérsic indices from the density profiles.  {#tbl:mcmc_props short="Properties of probabilistically-derived density profiles"}



## A comparison of density profiles

[@fig:mcmc_hists; @fig:mcmc_hists2] show the derived, MCMC histogram density profiles s compared to J+24. In general, both methodologies agree well. However, J+24 tend to systematically overestimate faint densities and confidently derive densities where the MCMC model fails. Because J+24 use only one or two components across the entire dwarf galaxy, the faint regions become sensitive to the spatial likelihood. Antlia II is an extreme case, where the satellite is hidden behind a large number of foreground stars. We derive a much more poorly-constrained and fainter total density profile than J+24, and the divergence between the MCMC method here and J+24's profile likely shows that *Gaia* observations are insufficient to properly constrain Antia II's density profile. This is likely due to uncertainties in the background density of satellite-like stars in Antlia II.  

. Once the true density of stars drops below the background of satellite-like stars, J+24's method will occasionally select consistent background stars with a density directly dependent on the spatial likelihood. The uncertainties in the outer regions are likely underestimated. 

To properly compare density profiles before background-limiting effects become important, we only calculate our profiles in the main text out to the radii in @tbl:mcmc_props. We derive these "background-limited" radii based on the outermost derived density in the MCMC histogram with an uncertainty lower than 1 dex. This typically corresponds to the empirical background from the CMD+PM samples in the main text. 



![Probabilistic density profiles](figures/mcmc_histograms.png){#fig:mcmc_hists}

Figure: A comparison between the MCMC histogram method and J+24. The MCMC samples in each bin are black transparent dots, and the J+24 derived density profiles (with their $P_{\rm sat} < 0.2$) are orange solid dots.



![Probabilistic density profiles continued](figures/mcmc_histograms2.png){#fig:mcmc_hists2}

Figure: A continuation of @fig:mcmc_hists.
