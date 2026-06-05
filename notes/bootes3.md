Boötes III (Boo3) is an unusual system. Boo3 is one of few MW satellites wish

Boötes III (Boo3) is one of few suspected tidally-disrupting Milky Way satellites, and is among the largest and most diffuse faint dwarf galaxies known. While having a small Milky Way pericentre (<10 kpc), Boo3 remains spare in the literature and has no resolved stream stars reported. Here, we will characterize the tidal disruption of Boo3 with idealized N-body simulations. Specifically, we investigate the impact of Boo3's orbit, the dwarf's initial structure, and the inner MW potential. We then assess the observational characteristics of disruption in Boo3 under each scenario and predict the properties of a possible stellar stream. 





# Introduction

The Sloan Digital Sky Survey heralded a new era for dwarf galaxies. Photographic plate surveys (e.g., PSS) were able to discover the *classical* Milky Way dwarf galaxies. But with digital sky surveys, deeper imaging and more quantitatively analyzed, astronomers discovered a host of new, ultra-faint dwarf galaxies. 

Among Milky Way satellites, Boötes III (Boo III) is a peculiar galaxy. Boo III is significantly larger than most dwarf galaxies (half-light radius of ~300 parsecs) while also being relatively faint. Boo III also has an orbit which may possibly cause disruption of the galaxy (see Fig. 1). 

Boo III was originally discovered in @grillmair2009 through matched filter analysis of SDSS data (along with four other streams). In particular, @grillmair2009 find Boo III within a stream dubbed "Styx" which they claim is the ongoing disruption of the dwarf galaxy. They note that Boo III might be disturbed and . They findthat the surface density profile of Boo III is extended $\Sigma \sim r^{-1}$ but nothing is well measured in this paper. The derived distance is $\sim 46$ kpc. 



@carlin+sand2018 present a sample of 15 RV members, first measuring the velocity dispersion, 



## History of the BooIII

Boo III was discovered in @grillmair2009 through matched filter method. In their paper, the only derived properties are the approximate distance (46 kpc or $18.35\pm0.01$ for distance modulus), the centroid (ra, dec = 209.3˚, 26.8˚), a density profile ($\Sigma \sim r^{-1}$). They discover the Styx stream and claim that Boo III might be associated (in the middle of the stream and at a similar distance). 



@correnti+bellazzini+ferraro2009 follow up the initial detection and detect and analyze a possible population of Blue Horizontal Branch (BHB) stars. From this population,they derive a similar density profile, a centroid of ($209.7\pm1.4, 26.8\pm0.6$ degrees), a distance modulus $18.58\pm0.05 \pm 0.14$. By extrapolating from the number of BHB stars detected, they infer$M_V = -5.8\pm0.5$ and an approximate ellipticity of 0.5. 



@carlin+2009 follow up the initial detection with MMT spectroscopy. They find 20 identified members (missing many due to high foreground). Their main results are a high velocity dispersion $\sigma_\text{v} = 14. 0 \pm 3.2$km/s, and the derived systemic velocity $\text{v} = 197.5\pm3.8$km/s. They also derive [Fe/H]$ \approx -2.1\pm0.2$, and a high metallicity dispersion. With a high dispersion, they interpret this system as likely undergoing active tidal disruption.



@massari+helmi2018 use Gaia DR 2 to derive proper motions for seven ultra faint dwarf galaxies, including Boo III. Their results from 34 possible members are:$\pmra = −1.21 \pm0.13$, $\pmdec =  −0.92 \pm0.17$, covariance = 0.23 (not includingg the 0.035mas/yr systematic).

@carlin+sand2018 perform the first orbital modelling of Boo III. They derive a proper motion of (μα cos d, μδ)=(−1.14, −0.98)±(0.18, 0.20) mas yr−1 based on LOS velocity selected members. By excluding members from @carlin+2009 with inconsistent proper motions, they further re-derive the systematic velocity (197.1\pm3.6) and dispersion(10.7\pm3.5km/s). Based on their orbital analysis, the pericentre is 10-12\pm6 kpc and they assert that Boo III is likely disrupting and is consistent with the position of Styx.

@moskowitz+walker2020: Stellar density profiles and structural parameters for many dwarfs including Boo III. 

@vivas+martinez-vazquez+walker2020: Updated RRL census. derive new distance to Boo III ($18.34  \pm 0.19$) with 7 RRL stars. (One previously known from @seaser+2014). They find two RRL beyond the tidal radius, and mention that the disruption of the galaxy is a known fact. 



@tau+2024 find possible RRL stars nearby Bootes III out to very far distances, but some may have came from Sgr stream. They find 32 total RRL members, but not in line with Styx and with a gradient consistent with Sgr contaimination. 

@jensen+2024: finding Boo III extended structure

@Yang+2025, weird features/morphology in Boo III with new observations (photometry). 



# The observational status of Boötes III

While the distance and proper motions of Bootes III appear to reasonably well-characterized, we revisit the structural and line-of-sight velocity properties of the galaxy. 



| parameter                       | value                                     | Source         |
| ------------------------------- | ----------------------------------------- | -------------- |
| $\alpha$                        | $209.47\pm0.11$                           | @struct_params |
| $\delta$                        | $26.68\pm0.06$                            | @struct_params |
| distance modulus                | $18.34 \pm 0.19$                          | vivas+2020     |
| distance                        | $46.56 \pm 4$ kpc                         | vivas+2020     |
| $\mu_\alpha \cos \delta$        | $-1.16 \pm 0.02 \pm 0.017$ mas yr$^{-1}$  | battaglia+2022 |
| $\mu_\delta$                    | $-0.88  \pm 0.01 \pm 0.017$ mas yr$^{-1}$ | battaglia+2022 |
| RV                              | $188 \pm 2.2$ km/s                        | @kin_params    |
| $\sigma_v$                      | $7.7_{-1.5}^{+2.0}$ km/s                  | @kin_params    |
| $R_h$ (Plummer, sphericalized)* | $44_{-6}^{+7}$                            | @struct_params |
| ell                             | $0.42_{-0.14}^{+0.11}$z                   | @struct_params |
| PA                              | $89\pm9$                                  | @struct_params |
| $M_V$                           | $-5.1\pm0.24$                             | @luminosity    |
| $M_\star / \mathrm{M}_\odot$    | $(14\pm3) \times10^3$                     | @luminosity    |
| [Fe/H]                          | $-2.5\pm0.1$                              | @kin_params    |



## Survey Data {#gaia_data}

To derive observed properties in what follows, we use three data sources: *Gaia* photometry and astrometry, and the Geha+2026 catalogue of radial velocities and metallicities. 

Our *Gaia* observations are processed based on an expansion of Jensen+2026. While their algorithm is able to detect 

In the appendix, we rederive parameters using DELVE instead of *Gaia* with similar results. While DELVE is deeper, the background is also higher and DELVE does not have proper motions to reduce background. As a result, the foreground in *Gaia* is at least 10 times lower than for DELVE, and we consider the *Gaia* parameters more reliable and less subject to systemics on varying incompletenesses of ground-based surveys

## Structural properties

We derive structural parameters using a Monte-Carlo Markov Chain (MCMC) model to fit an assumed density structure of the dwarf galaxy. 

The total likelihood is a mixture of the satellite likelihood ( $\mathcal L_\text{sat}$) and background likelihood ($\mathcal L_\text{bg}$) for each star:
$$
\mathcal{L} =  f\,\mathcal{L}_\text{sat} + (1-f) \mathcal{L}_\text{bg}
$$
where $f$ represents the fraction of stars in the current field belonging to the satellite. 

Following J+24, each likelihood is a product of a CMD, PM, and Spatial likelihood. 
$$
\mathcal L = \mathcal L_\text{CMD}\ \mathcal L_\text{PM}\ \mathcal L_\text{S}
$$
The CMD and PM likelihoods are fixed from J+24. The spatial likelihood is instead one of a Plummer or Sérsic density profiles 
$$
\Sigma_P(R) = \frac{1}{\pi R_h^2} \left(1 + \frac{R^2}{R_h^2}\right)^{-2}
$$

$$
\Sigma_S(R) = \Sigma_h\exp(-b_n (R/R_h)^{1/n})
$$

We also consider both circular and elliptical profiles. The circular $R$ is the Euclidian distance, and elliptical $R$ is the sphericalized coordinates given an ellipticity $e$ and position angle $\theta$.  **verify these**

For the priors, we assume 
$$
\begin{cases}
\Delta \xi \sim N(0, 15)' \\
\Delta \eta \sim N(0, 15)' \\
\ln R_h  \sim N(\ln 30, 1)' \\
f_\text{sat} \sim U(0, 0.5)
\end{cases}
$$
We also consider extensions of the base model, including an elliptical model with 
$$
\begin{cases}
e \sim U(0, 0.99) \\
\theta \sim U(0, 180)^\circ
\end{cases}
$$
In addition, we consider a Sérsic model, where $n\sim U(0.4, 12)$. 

While two (much brighter) globular clusters are in the same field, excluding these regions has a negligible impact on our derived parameters, likely because our priors constrain the location of the fitted object. From the @harris2010 database, NGC 5272 is at 205.54842, +28.37728 and of size 2.31', NGC 5466 is at 211.36371, 28.53444 and of size 2.30'. These correspond to elliptical (circular) radii of about 200' (230) and 170' (150'), so radial bins between about 140-230 are contaiminated with cluster. 



| model       | $\Delta \xi/'$ | $\Delta \eta/'$ | $e$                    | $\theta / \deg$  | $R_h\ /\ '$      | $n_\text{S\'ersic}$ | $N_\text{sat}$ |
| ----------- | -------------- | --------------- | ---------------------- | ---------------- | ---------------- | ------------------- | -------------- |
| Plummer     | $9.0\pm5.7$    | $-7.2 \pm 3.4$  | $0.42_{-0.11}^{+0.14}$ | $88.5\pm9.5$     | $43.8_{-6}^{+7}$ | --                  | $112 \pm 15$   |
| Exponential | $9.4\pm6$      | $-6.1\pm5$      |                        |                  |                  |                     | $130 \pm20$    |
| Sérsic      | $12\pm4$       | $-9_{-3}^{+4}$  | $0.23\pm0.14$          | $93^{+20}_{-19}$ | $67_{-15}^{+23}$ | $2.9_{-0.7}^{+0.9}$ | $57_{-7}^{+8}$ |

## Density profiles

While the likelihoods above can provide a membership list of stars, using this to then derive a density profile of the satellite risks "double fitting" the surface density. Instead, we derive density profiles by taking a sample of stars selected based on the likelihood of the foreground to background ratios, excluding the spatial likelihood terms. The resulting density profiles, and comparison MCMC parameterized density profiles, are shown in Fig. X



![image-20260408203004227](/Users/daniel/Library/Application Support/typora-user-images/image-20260408203004227.png)

Figure: The density profiles of Boo III using a sample of *Gaia* stars as selected based on their combined CMD+PM likelihood. 

## Luminosity

The only existing estimate of Boo III's luminosity is from @corenneti+2018 which simply count the number of HB stars and determine that $M_V \sim -5.8 \pm 0.5$. 

With a characterization of the number of stars in Gaia data with $G < 21$, we can determine an approximate luminosity of Boo III as follows (see also @munoz+2018).

To sample stars from the IMF, we use the following procedure:

1. Draw a new number of observed stars and distance modulus. In practice, we use the MCMC samples from the Plummer fit for the observed number of stars.
2. Draw stars from Kroupa IMF until we reach the target of "observed" stars (i.e. magnitudes above the threshold), using an isochrone (PADOVA Fe/H=-2.19, age=12Gyr) to convert masses into magnitudes.
3. Sum the total mass of all drawn stars and their luminosity. 



Our results for the magnitude uncertainties are in @tab:luminosity. We find magnitudes consistent with uncertainties (with the exp and Sérsic profiles favouring more stars). 

*Gaia* may not be complete to magnitude $G=21$. To verify our results, we also compute the absolute magnitude for stars in gaia selected to have a likelihood ratio (CMD+PM) of at least 0.01 and $G < 20$. While the number of observed stars is fewer than before, the absolute magnitude is consistent with the $G=21$ sample. We conclude that incompleteness towards fainter magnitudes (which should make the galaxy appear fainter) appears to be slight. 



| Model          | observed # stars  | abs magnitude  | stellar mass          |
| -------------- | ----------------- | -------------- | --------------------- |
|                | number            | $M_V$ mag      | $\times10^3\,M_\odot$ |
| Plummer-bright | $62_{-19}^{+26}$  | $-5.0\pm0.5$   | $13_{-4}^{+7}$        |
| Plummer        | $112\pm15$        | $-5.1\pm0.2$   | $14\pm3$              |
| Exp            | $130 \pm 19$      | $-5.3 \pm 0.2$ | $17_{-3}^{+4}$        |
| Sérsic         | $148_{-17}^{+19}$ | $-5.4 \pm 0.2$ | $19_{-3}^{+4}$        |
| DELVE          | $185\pm40$        | $-4.8\pm0.4$   | $11\pm4$              |

Table: The results from our derivations of the absolute magnitude. 



## Kinematics

Given the range of literature values of velocity dispersions, we re-derive the velocity dispersion. 

We use the suggested sample of members from @geha+2026 with membership probabilities greater than 0.5. (We have also tried using sigma-clipping and a probabilistic mixture model extending J+24 with similar final membership lists). This list contains 19 likely members. 

We then fit the velocity dispersion using the following model following @walker++2008?:
$$
\ln {\cal L} = \sum_i -\frac{1}{2}\frac{(\text v_i - \mu)^2}{(\delta\text v_i)^2 + \sigma^2} + -\frac{1}{2} \ln2\pi \sigma^2
$$
where $\text v_i$ is the perspective-motion corrected velocity of the $i$-th star, for each star, $\delta \text v_i$ is the velocity uncertainty, $\mu$ is the systemic velocity, and $\sigma$ is the velocity dispersion. We adopt a prior of $\mu \sim N(0, 100)$ for km/s in the GSR frame, and $\sigma \sim U(0, 20)$. 

We use a NUTS sampler with 16 chains 10,000 steps each. The resulting distributions are shown in @fig:velocity_dispersion. Compared to @geha+2026, our velocity dispersion is slightly higher and we have

![image-20260407165057367](/Users/daniel/Library/Application Support/typora-user-images/image-20260407165057367.png)

Figure: the resulting velocity dispersion fits. 



![image-20260402110023952](/Users/daniel/Library/Application Support/typora-user-images/image-20260402110023952.png)

# N-body methods

## Potential and galactocentric frame

- Astropy v4.0 Galactocentric frame
- @EP2020 potential

## Point particle orbits

- Agama: 100,000 orbits by sampling observed properties from Table @tbl:obs_props
- Selected orbits 

  We use Agama to calculate 100,000 orbits in the EP2020 potential. 

  Orbits selected from quantile method, median properties of near $\pm2\sigma$ pericentre. 





## Initial conditions



| Orbit      | -2 (1.5 kpc) | -1 (4 kpc) | Mean (7kpc) | +1 (12kpc) | +2 (18kpc) | +3 (26kpc) |
| ---------- | ------------ | ---------- | ----------- | ---------- | ---------- | ---------- |
| ra         |              |            |             |            |            |            |
| dec        |              |            |             |            |            |            |
| dist       |              |            |             |            |            |            |
| pm ra      |              |            |             |            |            |            |
| pm_dec     |              |            |             |            |            |            |
| rv         |              |            |             |            |            |            |
| pericentre |              |            |             |            |            |            |
| apocentre  |              |            |             |            |            |            |





| halo    | $v_{\rm circ,\ max}$ | $r_{\rm max}$ | $\sigma$ conc. | $v_\text{circ, end, req}$ | $h  / {\rm kpc}$ | $z_\text{mean}$ |
| ------- | -------------------- | ------------- | -------------- | ------------------------- | ---------------- | --------------- |
| average | 22                   | 3.9           | 0              |                           |                  |                 |
| compact | 30                   | 3.0           | +2             | 15?                       |                  |                 |
| heavy ? | 30                   | 5.7           | 0              |                           |                  |                 |





| halo        | 26   | 18   | 12   | 7    | 4    | 1.5  |
| ----------- | ---- | ---- | ---- | ---- | ---- | ---- |
| average     | 3*   | 2    | 1*   | --   | --   | --   |
| heavy       | 3    | 2-3  | 1    | 1    | --   | --   |
| **compact** | --   | 5*   | 3-4  | 2    | 1*   | --   |

Table: the number of pericentres which a given halo must experience to match the present day velocity dispersion under a specified pericentre. 

1. Compact halo; 5 x 18kpc peris (default 10 Gyr)
   1. `bootes3/1e5_v30_r3.0/5_peri_18kpc`
2. Compact halo; 1 x 4 kpc pericentre (1.5 Gyr sim)
   1. `bootes3/1e5_v30_r3.0/5_peri_18kpc`
3. Average halo: 3 x 26 kpc pericentre (default 10 Gyr)
   1. `bootes3/1e5_v22_r3.9/3_peri_26kpc`
4. Average halo: 1 x 12 kpc peri (1.5 Gyr)
   1. `bootes3/1e5_v22_r3.9/1_peri_12kpc`



# Results

## Orbits

@fig:orbits illustrates the range of possible orbits for Boo III



@fig:pericentre_distance illustrates the tight correlation of pericentre with distance among our sampled orbits. In addition, the lower panel of @fig:pericentre_distance shows the tangental proper motion dependence on distance. WIth present uncertainties, the distance determines the pericentre. 



Appendix @sec:extra_orbits compares the effects of including a Large Magellanic Cloud and changing the Milky Way potential on the orbit of Boo III. Since Boo III is on the opposite side of the galaxy as the LMC, Boo III is very slightly affected by including the LMC. In addition, the uncertainties in the Milky Way potential are far less than the uncertainties due to the distance to Boo III.





## N-body results (guess)

- Small pericentric models disrupt quickly, can likely only survive *one* pericentre!
- Moderate pericentres undergoe heavy tidal stripping and cannot survive long either
- Requires large pericentre to be long-term satellite
- regardless of orbit, heavy mass loss
- total mass near the scale length of Boo III
- Stream density undetectable faint. 
- Consistent with LCDM galaxy formation



# Discussion & Conclusions









# Appendix



# Halo by halo

## NFW

### v30_r1.0

### v30_r2.2

### v30_r3.0

+2 sigma more concentrated than Ludlow+2016 at $z=0$. 

Relatively resilient, can survive

## Cored



