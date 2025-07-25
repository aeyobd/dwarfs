# Data selection

## Describution of algorithm

To create a high-quality sample, J+24 select stars initially from Gaia within a 2--4 degree circular region centred on the dwarf satisfying: 

- Solved astrometry, magnitude, and colour. 
- Renormalized unit weight error, ${\rm ruwe} \leq 1.3$, ensuring high quality astrometry. `ruwe` is a measure of the excess astrometric noise on fitting a consistent parallax-proper motion solution [see @lindegren+2021].  
- 3$\sigma$ consistency of measured parallax with dwarf's distance (dwarf parallax is very small; with @lindegren+2021 zero-point correction). 
- Absolute proper motions, $\mu_{\alpha*}$, $\mu_\delta$, less than 10$\,{\rm mas\ yr^{-1}}$. (Corresponds to tangental velocities of $\gtrsim 500$ km/s at distances larger than 10 kpc.) 
- Corrected colour excess is within 3$\sigma$ of the expected distribution from @riello+2021. Removes stars with unreliable photometry. 
- De-reddened $G$ magnitude is between $22 > G > G_{\rm TRGB} - 5\sigma_{\rm DM}$. Removes very faint stars and stars significantly brighter than the tip of the red giant branch (TRGB) magnitude plus the distance modulus uncertainty $\sigma_{\rm DM}$. 
- Colour is between $-0.5 < {\rm BP - RP} <  2.5$ (dereddened). Removes stars substantially outside the expected CMD.

Photometry is dereddened with @schlegel+finkbeiner+davis1998 extinction maps.

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



For the satellite's spatial likelihood, J+24 consider both one-component and a two-component density models. The one component model is constructed as a single exponential profile  ( surface density $\Sigma \propto e^{R_{\rm ell} / R_s}$), with scale radius $R_s$ fixed to the value in table 1 of @MV2020a from @munoz+2018 (for a Sérsic fit). Additionally, structural uncertainties (for position angle, ellipticity, and scale radius) are sampled over to construct the final likelihood map. The two-component model instead adds a second exponential, $\Sigma_\star \propto e^{-R/R_s} + B\,e^{-R/R_{\rm outer}}$. The inner scale radius is fixed, and the outer scale radius and magnitude of the second component $R_{\rm outer}$, $B$  are free parameters. Structural property uncertainties are not included in the two-component model.

The PM likelihood is a bivariate gaussian with variance and covariance equal to each star's proper motions. J+24 assume the stellar PM errors are the main source of uncertainty.

 The satellite's CMD likelihood is based on a Padova isochrone [@girardi+2002]. The isochrone has a matching metallicity and 12 Gyr age (except 2 Gyr is used for Fornax). The (gaussian) colour width is assumed to be 0.1 mag plus the Gaia colour uncertainty at each magnitude. The horizontal branch is modelled as a constant magnitude extending blue of the CMD (mean magnitude of -2.2, 12 Gyr HB stars and a 0.1 mag width plus the mean colour error). A likelihood map is constructed by sampling the distance modulus in addition to the CMD width, taking the maximum of RGB and HB likelihoods.

The background likelihoods are instead empirically constructed. Stars stars outside of 5$R_h$ passing the quality cuts estimate the background density in PM and CMD space. The density is a sum of bivariate gaussians with variances based on Gaia uncertainties (and covariance for proper motions).  The spatial background likelihood is assumed to be constant. 

J+24 derive $\mu_{\alpha*}$, $\mu_\delta$, $f_{\rm sat}$ (and $B$, $R_{\rm outer}$ for two-component) through an MCMC simulation with likelihood from [@Eq:Ltot]. Priors are only weakly informative. The proper motion single component prior is same as @MV2020a: a normal distribution with mean 0 and standard deviation $100\ \kms$. If 2-component spatial, instead is a uniform distribution spanning 5$\sigma$ of single component case w/ systematic uncertainties. $f_{\rm sat}$  (and $B$) has a uniform prior  0--1. $R_{\rm outer}$ has a uniform prior only restricting $R_{\rm outer} > R_s$. The mode of each parameter from the MCMC are then reported and used to calculate the final $P_{\rm sat}$ values. 



## Additional density tests

In this section, we discuss additional tests and verification of the derived density profiles. In particular, we check that methodology (simpler cuts, circularized radii, algorithm version) do not substantially affect the density profile. We also compile density profiles presented in the literature as reference. In all cases, the density profiles appear to have excellent convergence out to $\log R_{\rm ell} / {\rm arcmin} \approx 1.8$, about the distance where the background dominates. 

Discuss selection criteria for DELVE and UNIONS samples, literature comparison, simple selection criteria, MCMC density profiles and when @jensen+2024 becomes background-limited.



![Scl density comparison](figures/scl_density_methods_extra.pdf){#fig:scl_density_extras}

Figure: Density profiles for various assumptions for Sculptor. PSAT is our fiducial 2-component J+24 sample, circ is a 2-component bayesian model assuming circular radii, simple is the series of simple cuts described, bright is the sample of the brightest half of stars (scaled by 2), DELVE is a sample of  RGB  stars (background subtracted and rescaled to match).



![Scl density methods](figures/scl_density_methods_j24.pdf){#fig:scl_density_j24_methods} 

Figure: Comparison of density profiles for each J+24 method. The fiducial is a 2-component elliptical model. However, the 1-component is still elliptical but only contains 1 component and the circular model assumes a circular outer density profile and bins in circular bins instead of elliptical bins. 



![UMi density comparison](figures/umi_density_methods_extra.pdf){#fig:umi_density_extras}

Figure: Similar to [@fig:scl_observed_profiles] except for Ursa Minor



![UMi density methods](figures/umi_density_methods_j24.pdf){#fig:umi_density_j24_methods} 

Figure: Similar to [@fig:scl_density_j24_methods] except for Ursa Minor. 





## Comparison to Literature

Here, we compare our density profiles against past derivations of density profiles for Sculptor and Ursa Minor
