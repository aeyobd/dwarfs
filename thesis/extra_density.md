# The reliability of derived density profiles {#sec:extra_density}

In this section, we describe the details of J+24's likelihoods. We then test the reliability of the derived density profiles for Scl and UMi, comparing methodologies, samples, and to literature. Finally, we present a comparison of J+24's samples to a non-parametric Bayesian density profile, finding similar results in the inner regions. Past a "limiting radius", J+24's density profiles become susceptible to the likelihood specification. However, all methods agree up to the limiting radii. In conclusion, the density profiles from J+24 are robust except in the very outer regions. 

## Bayesian membership probabilities

To create a high-quality sample, J+24 select stars initially from Gaia within a 2--4 degree circular region centred on the dwarf satisfying: 

- Solved astrometry, magnitude, and colour. 
- Renormalized unit weight error, ${\rm ruwe} \leq 1.3$. Ensures high quality astrometry. `ruwe` is a measure of the excess astrometric noise on fitting a consistent parallax-proper motion solution [see @lindegren+2021]. 
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

The background likelihoods are instead empirically constructed. Stars stars outside of 5$R_h$ passing the quality cuts estimate the background density in PM and CMD space. The density is a sum of bivariate gaussians with variances based on Gaia uncertainties (and covariance for proper motions).  The spatial background likelihood is assumed to be uniform over the field. 

J+24 derive $\mu_{\alpha*}$, $\mu_\delta$, $f_{\rm sat}$ (and $B$, $R_{\rm outer}$ for two-component) through an MCMC simulation with likelihood from [@Eq:Ltot]. Priors are only weakly informative. The proper motion single component prior is same as @MV2020a: a normal distribution with mean 0 and standard deviation $100\ \kms$. If using the 2-component model, the prior is instead a uniform distribution spanning 5$\sigma$ of single component case w/ systematic uncertainties. $f_{\rm sat}$  (and $B$) has a uniform prior from 0 to 1. $R_{\rm outer}$ has a uniform prior only restricting $R_{\rm outer} > R_s$. The mode of each parameter from the MCMC are then reported and used to calculate the final $P_{\rm sat}$ values. 

## Alternative sample selection {#sec:simple_selection}

As J+24's probabilistic model may systematically bias the density profile (see discussion below), we consider samples selected both from *Gaia* and external surveys using simple, absolute selection cuts.

For *Gaia* data, we filter stars by four simple selection criteria. We require `ruwe` < 1.3 as above, selecting only stars with high-quality astrometry. We then require the star's parallax to be $3\sigma$-consistent with the dwarf's distance. Next, we select stars with proper motions with $1\masyr$​ radius of the galaxy's measured value. Finally, we filter stars to lie on the dwarf's CMD in ($G_{\rm BP} - G_{\rm RP}$,$G$), within the vertices defined in @tbl:colour_cuts. This polygon is derived based on the colour-magnitude diagram of stars in the centre of the dwarf galaxy, where the dwarf dominates the stellar density. 

In addition to *Gaia*, we select stars in the Scl field from the deeper, photometric DELVE survey [@drlica-wagner+2022]. We select sources DELVE data release 2 within an ellipse of radius 150 arcminutes. We then only include sources which are categorized as likely stars and sources with reliable (flags <=4) $g$ and $r$ magnitudes. We then filter by the $g$ and $r$ magnitudes, selecting the empirical CMD cut using the vertices in @tbl:colour_cuts. 

Ursa Minor is within the Ultraviolet Near-Infrared Optical Northern Survey [UNIONS @gwyn+2025]. We select sources within an ellipse of radius 230 arminutes. We keep sources with no `FLAGS_CFIS` set and the `s21` and `s31` both $<3$ (i.e. not extended sources). We then use the polygon bounded by the vertices in @tbl:colour_cuts to select the observed red giant branch, horizontal branch, and main sequence of the dwarf. 

![CMD cuts](figures/extra_cmd_selection.png)

Figure: Colour-magnitude diagram cuts used in samples in this section. Orange points are selected, black are non-selected stars within $1R_h$, and grey stars are for all stars. **todo, only plot black/grey and cmd outline...**





![DELVE and UNIONS spatial distribution of stars](figures/delve_unions_tangent.pdf)

Figure: The distribution of DELVE and UNIONS selected stars in Scl and UMI

## Uncertainties and possible biases in *Gaia*-derived density profiles {#sec:density_extra}

In this section, we discuss additional tests and verification of the derived density profiles. In particular, we check that alternative methodologies do not substantially affect the density profile. In all cases, the density profiles appear to have excellent convergence out to $\log R_{\rm ell} / {\rm arcmin} \approx 1.8$, about the distance where the background dominates. 

**Spatial likelihood.** J+24's algorithm was designed specifically to detect the presence of a density excess in satellites and identify individual stars for spectroscopic follow-up. We are instead interested in constraining a dwarf galaxy's stellar structure. As such, one potential problem with using J+24's candidate members is that the stellar density model is either a single or double exponential. If the (double) exponential model does not adequately describe a dwarf galaxy, then J+24's derived density profiles may be systematically biased towards the assumed (double) exponential profile. 

As an example, in the lower panels of @fig:scl_density_extras; @fig:umi_density_extras, the density profile derived from the fiducial (J+24) sample appears to be well-constrained at densities far below the CMD+PM or `simple` model's background density.  In addition, the fiducial and 1-component density profiles deviate, despite being derived identically besides the spatial likelihood assumption. 

These stars are likely selecting stars from the statistical MW background consistent either galaxy's CMD and PMs, recovering the assumed density profile. As a result, the reliability of these density profiles below the CMD+PM background may be questionable. A more robust analysis, removing this particular density assumption, would be required to more appropriately represent the knowledge of the density profile as the background begins to dominate. 

We revisit possible biases in spatial likelihood in @sec:mcmc_hists. 

**Completeness**. *Gaia* shows high but imperfect completeness, particularly showing limitations in crowded fields and for faint sources ($G\gtrapprox20$). In @fabricius+2021, for the lowest density globular clusters, the completeness down to $G\approx 20$ is $\sim 80\%$. As the lowest density globular clusters are still far denser than dwarf galaxies, so the completeness is likely higher for classical dwarfs. *Gaia* is likely more incomplete in terms of BP-RP photometry (as we use here) and for stars with small pairwise separations ($\lesssim 1.4''$). The influence of these details is unclear. Nevertheless, *Gaia* should be more complete for brighter stars.

We find little evidence for magnitude-dependent biases in the density profiles. In [@fig:scl_density_extras; @fig:umi_density_extras] we show density profiles derived using only stars brighter than the median magnitude. We find no substantive differences compared to our fiducial density profile. If *Gaia*'s incompleteness is inhomogeneous, then these variations do not appear to affect the density profiles.

If *Gaia* has systematic density biases, the density profile should vary among different facilities. We show the density profiles derived from deep, ground-based UNIONS and DEVLE photometry (lower panels of [@fig:scl_density_extras; @fig:umi_density_extras]). For both Scl and UMi, these alternative surveys result in near-identical density profiles. Variations between density profiles derived from these surveys appear to also be with uncertainties. 

**Structural uncertainties**. The assumed structural parameters of a dwarf may be misrepresentative or vary with radius. J+24 derive probabilities where the outer stellar component is circular instead of elliptical. Using the circular 2-component probabilities and binning in circular radius, we find little difference in the derived density profiles (the `circ` points in [@fig:scl_density_extras; @fig:umi_density_extras]). Thus, our conclusions are robust to extreme variations in the ellipticity. Reasonable changes in other structural parameters (centre, position angle) are unlikely to have a more substantial effect.



![Sculptor density methodology comparison](figures/density_methods_extra.pdf){#fig:scl_density_extras}

Figure: Density profiles for various assumptions for Sculptor. `2-exp` is the fiducial double-exponential-likelihood J+24 sample, `1-exp` instead is a one-component exponential spatial likelihood, and `simple` uses the simple position-independent selection criteria   `circ` is a 2-component bayesian model assuming circular radii, simple is the series of simple cuts described, bright is the sample of the brightest half of stars (scaled by 2).



## A Bayesian density profile {#sec:mcmc_hists}

To address the concerns discussed above, we consider here a non-parametric model to fit the density in each bin. 



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



### Methodology

This model extends the J+24 framework with two notable differences. (1) The systemic proper motions are held fixed. This allows for a much more efficient calculation of the likelihoods. (2) The spatial likelihood for the satellite is a piecewise constant density function.While this model theoretically has many more parameters to fit (one for each bin), each parameter is independently estimated based on only the stars in that specific bin.

We segregate the data into radial bins and independently fit the model to only stars in each bin. This fits many small models instead of a global large model. In detail, this means that, in each bin, the normalized density of the background and satellite are *equal*. The only fit parameter is the $f_{\rm sat}$ in the particular bin. Therefore, the likelihood in each bin is the same as above except ${\cal L}_{\rm space} = 1$ for the foreground and background. 

In detail, we parameterize $f_{\rm sat}$ in terms of the log-relative satellite density. The model construction is then
$$
\begin{align}
\theta_i &\equiv \log_{10}({\Sigma_{\rm sat}}/{\Sigma_{\rm bg}}) \\
f_{\rm sat} &= \frac{10^{\theta_i}}{1 + 10^{\theta_i}} \\
\theta_i &\sim {\rm Uniform}(-12, 6) \\
{\cal L}_{\rm tot} &= f_{\rm sat}{\cal L}_{\rm CMD, PM} + (1-f_{\rm sat}){\cal L}_{\rm CMD, PM} \\
\end{align}
$$


For each satellite, we hold the structural parameters fixed to those in J+24, create bins which are the larger of the bin containing 20 stars or 0.05dex in $\log R_{\rm ell}$. We then sample each chain 1000 steps for 48 walkers using Turing.jl and the No-U Turns Sampler. 

Based on the $f_{\rm sat}$ values in each bin, the number of members is just the number of stars in the bin times $f_{\rm sat}$ for the bin, allowing the direct derivation of the density profiles from this model.

We have also briefly explored a model which solves for $f_{\rm sat}$ for random realizations of the global structural uncertainties. These do not appear to influence the results substantially. Another possible inconsistency is that the half-light radii are estimated from a different sample. We do briefly rederive structural parameters for the sample using Sérsic fits (although with a fixed origin).

### Results

[@fig:mcmc_hists; @fig:mcmc_hists2] show the derived, MCMC histogram density profiles compared to J+24. In general, both methodologies agree well. However, J+24 tend to systematically overestimate faint densities and confidently derive densities where the MCMC model fails. Because J+24 use only one or two components across the entire dwarf galaxy, the faint regions become sensitive to the spatial likelihood. Antlia II is an extreme case, where the satellite is hidden behind a large number of foreground stars. We derive a much more poorly-constrained and fainter total density profile than J+24, and the divergence between the MCMC method here and J+24's profile likely shows that *Gaia* observations are insufficient to properly constrain Antia II's density profile. This is likely due to uncertainties in the background density of satellite-like stars in Antlia II.  

Once the true density of stars drops below the background of satellite-like stars, J+24's method will occasionally select consistent background stars with a density directly dependent on the spatial likelihood. The uncertainties in the outer regions are likely underestimated. 

To properly compare density profiles before background-limiting effects become important, we only calculate our profiles in the main text out to the radii in @tbl:mcmc_props. We derive these "background-limited" radii based on the outermost derived density in the MCMC histogram with an uncertainty lower than 1 dex. This typically corresponds to the empirical background from the CMD+PM samples in the main text. 



![Probabilistic density profiles](figures/mcmc_histograms.png){#fig:mcmc_hists}

Figure: A comparison between the MCMC histogram method and J+24. The MCMC samples in each bin are black transparent dots, and the J+24 derived density profiles (with their $P_{\rm sat} < 0.2$) are orange solid dots.



![Probabilistic density profiles continued](figures/mcmc_histograms2.png){#fig:mcmc_hists2}

Figure: A continuation of @fig:mcmc_hists.



## Comparison to literature

We finally compare our density profiles against a sample of literature-derived density profiles in [@fig:scl_lit_profiles; @fig:umi_lit_profiles]. The J+24-derived density profiles are representative of state-of-the-art density profiles. A few photometric surveys may extend to slightly fainter densities, but typically with large uncertanties. Deviations between sources are slight despite the range of methods across time. All density profiles extending into the outskirts of Scl and UMi show a similar overdensity to the ones we find. The extended stellar profiles of Scl and UMi appear to be a robust result across the literature. 



![Sculptor literature density profiles](figures/scl_literatre_profiles.pdf){#fig:scl_lit_profiles}

Figure: A comparison of the Scl derived density profile and historical works. Unlike most profiles in this thesis, this density profile is plotted with respect to the semi-major elliptical radius ($a = R_{\rm ell} / \sqrt{1-{\rm ell}}$). The solid black line is a 2D exponential with corresponding scale radius to Scl's half-light radius, and the residuals in the bottom panel are with respect to this profile. **remove delve, umi title...,, combine??**

![Ursa Minor literature density profiles](figures/umi_literature_profiles.pdf)

Figure: Similar to @fig:scl_lit_profiles except for Ursa Minor. The references are (in order) @sato+2025; @palma+2003; @martinez-delgado+2001; @kleyna+1998; @IH1995; and @Hodge1964. 





## Summary

In this Appendix, we discussed the methodological details of the J+24 sample selection algorithm. We consider possible biases due to the assumed spatial likelihood, *Gaia*'s completeness, and structural parameters. In all cases, we find these assumptions likely do not cause major biases in the derived density profiles. However, we show that J+24's density profiles become unreliable when dropping below the background of satellite-like stars. We address this concern through an algorithm rederiving the density in elliptical bins. We find excellent agreement with J+24 for most galaxies out to a background-limited radius. Finally, comparing our density profiles against the literature, we find our results to be consistent. We conclude that the detection of an extended density profile in Scl and UMi is robust to incompleteness, methodology, alternative surveys, and across the literature.
