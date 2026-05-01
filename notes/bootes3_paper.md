*Abstract*. Boötes III (Boo3) is one of few suspected tidally-disrupting Milky Way satellites, and is among the largest and most diffuse faint dwarf galaxies known. While having a small Milky Way pericentre (<10 kpc), Boo3 remains spare in the literature and has no resolved stream stars reported. Here, we will characterize the tidal disruption of Boo3 with idealized N-body simulations. Specifically, we investigate the impact of Boo3's orbit, the dwarf's initial structure, and the inner MW potential. We then assess the observational characteristics of disruption in Boo3 under each scenario and predict the properties of a possible stellar stream. 





# Introduction

Only ${\sim}5$ dwarf galaxy satellites of the Milky Way show signs of tidal perturbations. While the Sagittarius dSph forms a stream wrapping across the entire sky (e.g., ) and the Large Magellanic Cloud (LMC)  has an extensive gas tail, the only other known tidal streams are from the Tucana III, Antlia II, and Crater II dwarf spheroidals (dSph). 

Boötes III is the next most promising candidate for tidal disruption. With a pericentres of ${\sim} 7\,$kpc, Boo III should be tidally disrupting (e.g., ). Indeed, @carlin+sand2018 claim that Boo III's orbit is consistent with the formation of the Styx stellar stream, detected alongside Boo III in @grillmair2009. In addition, the high initial velocity dispersions reported for Boo III of $\gtrsim 10$ km/s (@calin+2009, @carlin+sand2018) would potentially be consistent with active tidal disruption contaiminating these measurements (most similar ultra-faint dwarfs have velocity dispersions $\lesssim 5$). 

@carlin+sand2018 perform the first orbital modelling of Boo III. They derive a proper motion of (μα cos d, μδ)=(−1.14, −0.98)±(0.18, 0.20) mas yr−1 based on LOS velocity selected members. By excluding members from @carlin+2009 with inconsistent proper motions, they further re-derive the systematic velocity (197.1\pm3.6) and dispersion(10.7\pm3.5km/s). Based on their orbital analysis, the pericentre is 10-12\pm6 kpc and they assert that Boo III is likely disrupting and is consistent with the position of Styx.

@moskowitz+walker2020: Stellar density profiles and structural parameters for many dwarfs including Boo III. 

@vivas+martinez-vazquez+walker2020: Updated RRL census. derive new distance to Boo III ($18.34  \pm 0.19$) with 7 RRL stars. (One previously known from @seaser+2014). They find two RRL beyond the tidal radius, and mention that the disruption of the galaxy is a known fact. 



@tau+2024 find possible RRL stars nearby Bootes III out to very far distances, but some may have came from Sgr stream. They find 32 total RRL members, but not in line with Styx and with a gradient consistent with Sgr contaimination. 

@jensen+2024: finding Boo III extended structure

@Yang+2025, weird features/morphology in Boo III with new observations (photometry). 



# The observational status of Boötes III

We revisit the structural and line-of-sight velocity properties of the galaxy. @tbl:obs_props shows our recommended observational parameters for Boo III. 



| parameter                | Units              | value                       | Reference      |
| ------------------------ | ------------------ | --------------------------- | -------------- |
| $\alpha$                 | deg                | $209.47\pm0.11$             | @struct_params |
| $\delta$                 | deg                | $26.68\pm0.06$              | @struct_params |
| distance modulus         | mag                | $18.34 \pm 0.19$            | vivas+2020     |
| distance                 | mag                | $46.56 \pm 4$ kpc           | vivas+2020     |
| $\mu_\alpha \cos \delta$ | mas yr$^{-1}$      | $-1.16 \pm 0.02 \pm 0.017$  | battaglia+2022 |
| $\mu_\delta$             | mas yr$^{-1}$      | $-0.88  \pm 0.01 \pm 0.017$ | battaglia+2022 |
| $\mathrm{v}_\text{los}$  | km s$^{-1}$        | $188 \pm 2.2$               | @kin_params    |
| $\sigma_\mathrm{v}$      | km s$^{-1}$        | $7.7_{-1.5}^{+2.0}$ km/s    | @kin_params    |
| $R_h$*                   | arcmin             | $44_{-6}^{+7}$              | @struct_params |
| ellipticity              | --                 | $0.42_{-0.14}^{+0.11}$      | @struct_params |
| position angle           | deg                | $89\pm9$                    | @struct_params |
| $M_V$                    | mag                | $-5.1\pm0.24$               | @luminosity    |
| $M_\star$                | $\mathrm{M}_\odot$ | $(14\pm3) \times10^3$       | @luminosity    |
| [Fe/H]                   | dex                | $-2.5\pm0.1$                | @kin_params    |

Table: Recommended and derived properties of Boo III. Rows: right ascension ($\alpha$), declination $\delta$, distance modulus, distance, absolute proper motion in right ascension $\mu_\alpha \cos\delta$, proper motion in declination $\mu_\delta$, line of sight (los) velocity $\mathrm{v_{los}}$, LOS velocity dispersion $\sigma_\textrm{v}$, tentative half-light radius $R_h$, ellipticity, position angle, absolute V-band magnitude, total stellar mass, and metallicity. 

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



![image-20260429160516556](/Users/daniel/Library/Application Support/typora-user-images/image-20260429160516556.png)



| model       | $\Delta \xi/'$ | $\Delta \eta/'$ | $e$                    | $\theta / \deg$  | $R_h\ /\ '$      | $n_\text{S\'ersic}$ | $N_\text{sat}$ |
| ----------- | -------------- | --------------- | ---------------------- | ---------------- | ---------------- | ------------------- | -------------- |
| Plummer     | $9.0\pm5.7$    | $-7.2 \pm 3.4$  | $0.42_{-0.11}^{+0.14}$ | $88.5\pm9.5$     | $43.8_{-6}^{+7}$ | --                  | $112 \pm 15$   |
| Exponential | $9.4\pm6$      | $-6.1\pm5$      |                        |                  |                  |                     | $130 \pm20$    |
| Sérsic      | $12\pm4$       | $-9_{-3}^{+4}$  | $0.23\pm0.14$          | $93^{+20}_{-19}$ | $67_{-15}^{+23}$ | $2.9_{-0.7}^{+0.9}$ | $57_{-7}^{+8}$ |

## 



While the likelihoods above can provide a membership list of stars, using this to then derive a density profile of the satellite risks "double fitting" the surface density. Instead, we derive density profiles by taking a sample of stars selected based on the likelihood of the foreground to background ratios, excluding the spatial likelihood terms. The resulting density profiles, and comparison MCMC parameterized density profiles, are shown in Fig. X



![image-20260429155837759](/Users/daniel/Library/Application Support/typora-user-images/image-20260429155837759.png)

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

Figure: The resulting velocity dispersion fits. Black curves illustrate sampled Gaussian fits to the distribution from the MCMC model, and the orange points with errorbars are the observed velocity member stars.



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
| heavy ? | 40                   | 3.0           | +3             |                           |                  |                 |





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



![image-20260429143711521](/Users/daniel/Library/Application Support/typora-user-images/image-20260429143711521.png)

Figure: Orbits of Boo III with different pericentres. 



![image-20260429150209385](/Users/daniel/Library/Application Support/typora-user-images/image-20260429150209385.png)

Figure: The pericentre (top) and total angular momentum (bottom) as a function of heliocentric distance. Green points show the sampled MC orbits, and the large dots show the example orbits discussed in the text.



![image-20260429151837975](/Users/daniel/Library/Application Support/typora-user-images/image-20260429151837975.png)

Figure: The distribution of possible pericentres, with the examples considered here marked along the bottom.



## Tidal evolution

![image-20260430082213646](/Users/daniel/Library/Application Support/typora-user-images/image-20260430082213646.png)

![image-20260430082322360](/Users/daniel/Library/Application Support/typora-user-images/image-20260430082322360.png)

Figure: Two examples of the final dark matter and stellar distributions.



![image-20260430081434832](/Users/daniel/Library/Application Support/typora-user-images/image-20260430081434832.png)



![image-20260430081340817](/Users/daniel/Library/Application Support/typora-user-images/image-20260430081340817.png)

Figure: The final distributions of stars for each model considered above. 



Figure: The final stellar distribution of the cored model. Interestingly, a core of equal size does not strongly affect tidal evolution as the velocity dispersion constrains. 



![image-20260430080733367](/Users/daniel/Library/Application Support/typora-user-images/image-20260430080733367.png)

Figure: The along-stream density 





# Models to run

- Very extended initial stars
- More compact halo
- Low final velocity dispersion
- More cored models





# Discussion & Conclusions



- Caveat about orbital evolution
- More precise distance and velocity dispersion key to pinning down Boo III evolution—future/ongoing facilities.







# Appendix

## Structural parameters in DELVE

However, with a cut on mag < 23, the total number of main sequence stars is ~ 151. 



For DELVE (DR2), we select all stars within 6 arcminutes of BooIII's literature centre (209.3, 26.8).[^DELVE DR3 has incomplete coverage near Boo3. ] We require the sources to further satisfy:  

- g and r PSF magnitudes and their errors are defined
- g and r PSF magnitudes greater than the 5-$\sigma$ PSF depth (24.3, and 23.9)
  - In the field, these are 23.6 and 23.0 for the minimum around Boo III
  - As the turnover occurs at lower magnitudes, our reliable depth is limited to $r\sim 22$. 
- likely star source (extended_class_g <= 1).

In total, this returns 558,046 sources. We then correct the magnitudes for dust extinction using the @SFD dustmap with coefficients from @aboot+2018 with the dustmaps package. 



## Creating a matched filter

We use the matched-filter method to appropriately weight stars by their likelihood of belonging to the galaxy according to the matched filter

For the model of Boo III's CMD, we adopt a Padova isochrone of age 12 Gyr and metallicity -2.1, and weighted by a Kroupa (2001) IMF.  The isochrone is then convolved with a gaussian in magnitude of 0.19 mag (the distance modulus uncertainty) and a gaussian  in colour of width 0.05 mag plus the typical colour uncertainty added in quadrature [^colour_uncertainty]

[^colour_uncertainty]: We determine the typical colour uncertainty to a least-squares fit of the g-r colour uncertainty with g-magnitude assuming a exponential-quadratic formula, i.e. *fill in*] 

Fig. X shows the resulting distribution. Most of the stellar mass is in the MS as expected, with a fainter RGB and HB. 

The background density is derived from all stars between 5 and 6 degrees from the centre adopting a bandwidth of (0.05, 0.1) mag in colour and g-magnitude. 

The matched filter is calculated over a range from 23 to 15 mag and -0.8 to 1.7 in colour, with 200 colour bins and 300 magnitude bins. 

Finally, the likelihood ratio is calculated from the ratio of these two maps, and interpolated as appropriate. The right-hand panel of Fig. 1 illustrates this map. The HB and MS are the most critical regions for identifying the galaxy.



![image-20260327135427018](/Users/daniel/Library/Application Support/typora-user-images/image-20260327135427018.png)

Figure: The creation of a matched filter for Boo III in the *g-r* CMD.  In all cases, the colours are linear and arbitrarily scaled. **Left:** The expected density of stars belonging to Boo III based on our isochrone model. **Middle:** An empirical estimate of the background density of 



Figure: Matched filter of delve members



Figure: Density profile of DELVE stars





## The effects of the LMC and potential choice





## Additional Models

## 

![image-20260430075012123](/Users/daniel/Library/Application Support/typora-user-images/image-20260430075012123.png)

Figure: The final distribution of dark matter for the x.x simulation. 



![image-20260430075003149](/Users/daniel/Library/Application Support/typora-user-images/image-20260430075003149.png)

