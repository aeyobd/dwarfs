# Background & Past Work

- What is dark matter? Why do we look at dwarfs?
- Forms of dark matter, lambda-CDM, and dwarf galaxies
- How does gravity affect dwarfs, theory of tidal perturbations
  - EN21, peñarrubia+09, etc.
- Instances of dwarfs undergoing weird processess
- Alternative processes and uncertainties in the evolution of dwarfs

### Sculptor

- sestito+23a
- tolstoy+22
- battaglia+09
- arroyo-polonio+22



| parameter                | value                                                        | Source    |
| ------------------------ | ------------------------------------------------------------ | --------- |
| $\alpha$                 | $15.0183 \pm 0.0012$˚                                        | M+18      |
| $\delta$                 | $-33.7186 \pm 0.00072$˚                                      | M+18      |
| distance                 | $83.2 \pm 2$ kpc                                             | Tran+22   |
| $\mu_\alpha \cos \delta$ | $0.099 \pm 0.002 \pm 0.017$ mas yr$^{-1}$                    | MV20a     |
| $\mu_\delta$             | $-0.160 \pm 0.002_{\rm stat} \pm 0.017_{\rm sys}$ mas yr$^{-1}$ | MV20a     |
| RV                       | $111.03 \pm 0.23$                                            | This work |
| $\sigma_v$               | $9.61\pm0.16$                                                | This work |
| $r_h$                    | $12.33 \pm 0.05$ arcmin                                      | MV20*     |
| ell                      | $0.36 \pm 0.01$                                              | M+18      |
| PA                       | $92\pm1$                                                     | M+18      |
| $M_V$                    | $-10.82\pm0.14$                                              | M+18      |
| $\Upsilon_\star$         | $1.5 \pm 0.3$                                                | assumed   |

Table: Measured properties of Sculptor

## Gaia Membership Selection

![Sculptor selection criteria. This will be a super long caption to test the figure list.](figures/scl_selection.png){#fig:sculptor_selection}

We use J+24 data. J+24 select members using a multi-component Baysian algorithm:

- Remove stars with poor astrometry or photometry, no colour excess (@lindegren+2018 equation C.2), 3$\sigma$ consistency of measured parallax with dwarf distance (near zero with @lindegren+2018 zero-point correction), and absolute RA and Dec proper motions less than 10$\,{\rm mas\ yr^{-1}}$.

- Spatial likelihood based on a double exponential $\Sigma_\star \propto e^{-r/r_s} + B\,e^{-r/r_{\rm outer}}$ where the inner scale radius is fixed. 

- Stars are assigned a likelihood based on the location on the CMD (using Padova isochrones including an intrinsic 0.1 CMD width in colour convolved with colour and distance modulus)

- Background KDE density maps for the CMD and PM are constructed using the other quality-selected stars outside of $5R_h$, where the satellite density would be orders of magnitude less than the background (even in the presence of extended tidal features).

- Likelihoods normalized to unity to represnt a PDF

  

  Membership probabilities are then given by

$$
P_{\rm sat} = \frac{f_{\rm sat}{\cal L}_{\rm sat}}{f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}} = \frac{1}{1 + \frac{(1-f_{\rm sat}){\cal L}_{\rm bg}}{f_{\rm sat}{\cal L}_{\rm sat} }}
$$

where $f_{\rm sat}$ is the fraction of stars belonging to the system inside the given field, ${\cal L}_{\rm sat}$ is the likelihood of a star belonging to the satellite, and ${\cal L}_{\rm bg}$ is the likelihood of the star belonging to the background. Each likelihood is calculated as a product of the CMD, PM, and spatial likelihoods:
$$
{\cal L} = {\cal L}_{\rm space}\ {\cal L}_{\rm PM}\ {\cal L}_{\rm CMD}
$$



Not shown here, we explore simple cuts of the stars, using absolute cuts in parallax, proper motions, and the CMD. The results are similar to the nospace model.

The above formula suggests that a cut in $P_{\rm sat}$ is equivalent to the cut in likelihoods
$$
\frac{{\cal L}_{\rm sat}}{{\cal L}_{\rm bg}} > \frac{(1-f_{\rm sat})/f_{\rm sat}}{1/P_{\rm sat}- 1}
$$

Note that if we remove the spatial component of the likelihood, then $f_{\rm sat}$ represents a global normalization

### Searches for tidal tails

![Density profiles](figures/scl_tidal_tails.png){#fig:sculptor_tidal_tails}

- The above figure shows the distribution of member stars (orange), PM & CMD only selected stars (blue) and all stars (passing quality cuts, black)

- There are no apparent overdensities in the PM & CMD only selected stars to suggest the presence of a tidal tail

- This means that at least at the level of where the background density dominates, we can exclude models which produce tidal tails brighter than a density of $\Sigma_\star \approx 10^{-2}\,\text{Gaia-stars\ arcmin}^{-2} \approx 10^{-6} \, {\rm M_\odot\ kpc^{-2}}$ (TODO assuming a distance  of ... and stellar mass of ...). 

  


## Density Profile Reliability and Uncertainties

- How well do we know the density profiles? 
- What uncertainties affect derived density profiles? 
- Can we determine if Gaia, structural, or algorithmic systematics introduce important errors in derived density profiles?



J+24's algorithm takes spatial position into account, assuming either a one or two component exponential density profile. When deriving a density profile, this assumption may influence the derived density profile, especially when the galaxy density is fainter than the background of similar appearing stars. To remedy this and estimate where the background begins to take over, we also explore a cut based on the likelihood ratio of only the CMD and PM components. This is in essence assuming that the spatial position of a star contains no information on it's membership probability (a uniform distribution like the background)

To incorporate the structural uncertainties and robustly model the sampling uncertainty, we construct the following bootstrap model

- Centre is varied by a centring error, estimated from the standard normal error of the positions plus the systematic shift of the mean

- Position angle and ellipticity are sampled from a normal distribution given the reported uncertainties 

- $f_{\rm sat}$ is sampled from ...

  

  **TODO**: Look into normalization of Likelihoods and check how $f_{\rm sat}$ matters. If $f_{\rm sat}$ is related to the normalizations of fg / bg densities, and other likelihoods are area-normalized to 1, then this makes life much easier. 

  Test if psat weighted density profiles are similar

  Save MC density profile outputs

  

  

  







![Density profiles](figures/scl_density_methods.png){#fig:sculptor_observed_profiles}

Figure: add simple cuts & delve to convince more...



# Comparison of the Classical dwarfs

- Using J+24 data, we validate
  - Check that PSAT, magnitude, no-space do not affect density profile shape too significantly
- Our "high quality" members all have > 50 member stars and do not depend too highly on the spatial component, mostly corresponding to the classical dwarfs
- We fit Sérsic profiles to each galaxy
  - The Sérsic index, $n$, is a measure of the deviation from an exponential. Exponentials have $n=1$, whereas more extended dwarf galaxies will have higher $n$
- To better estimate the uncertainties due to unknown galaxy properties and flexibility in the likelihood cut, we can 

![image-20250313130050775](/Users/daniel/Library/Application Support/typora-user-images/image-20250313130050775.png)

![image-20250313130110550](/Users/daniel/Library/Application Support/typora-user-images/image-20250313130110550.png)

![image-20250313130043114](/Users/daniel/Library/Application Support/typora-user-images/image-20250313130043114.png)



## Radial Velocity Measurements



## Summary

- Of the classical dwarfs, UMi & Scl stand out statistically, with high $n$ given their luminosity
- Including fainter dwarf galaxies, Boo 3 and Boo 1 appear to also have extended density profiles
  - Deeper data would be required to robustly measure this



**TODO**: Use updated density profiles (nospace) with uncertainties included to MCMC fit Sérsic profiles to every dwarf galaxy. 
