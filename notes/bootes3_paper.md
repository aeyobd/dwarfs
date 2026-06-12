Boötes III (Boo3) is one of few suspected tidally-disrupting Milky Way satellites, and is among the largest and most diffuse faint dwarf galaxies known. While having a small Milky Way pericentre (<10 kpc), Boo3 remains spare in the literature and has no resolved stream stars reported. Here, we will characterize the tidal disruption of Boo3 with idealized N-body simulations. Specifically, we investigate the impact of Boo3's orbit, the dwarf's initial structure, and the inner MW potential. We then assess the observational characteristics of disruption in Boo3 under each scenario and predict the properties of a possible stellar stream. 

# Introduction

Only ${\sim}5$ dwarf galaxy satellites of the Milky Way show signs of tidal perturbations. The most dramatic examples are the Sagittarius dSph's extensive stream wrapping across the entire sky and the Large Magellanic Cloud (LMC) 30?? degree gas tail. Otherwise, the only other known tidal streams are from the Tucana III, Antlia II, and Crater II dwarf spheroidals (dSph). Indeed, most dwarf galaxies are expected to orbit too far from the Milky Way to be susceptible to tides (e.g., @pace+2020).

Boötes III is a unique candidate for tidal disruption. With a pericentres of ${\sim} 7\,$kpc, Boo III orbit implies the galaxy should be tidally disruption. In addition, Boo III appears to be situated within the *Styx* stream discovered in @grillmair2009. Boo III furthermore has had high previously measured velocity dispersions for an ultra-faint dwarf of $\gtrsim 10$ km/s (@calin+2009, @carlin+sand2018). Recently, extended structures around Boo III have been tentatively identified in @jensen+2024, @tau+2024, @yang+2025. The literature appears to coincide that Boo III is actively disrupting and the progenitor of the Styx stream. 

However, no resolved members of the *Styx* stellar stream have been identified. @jensen+2026 fail to confirm the Styx stellar stream with modern SDSS and DELVE data when repeating the methods in @grillmair2009. Furthermore, no one has modelled the tidal disruption of Boo III beyond test-particle orbits to test whether Styx has the appropriate morphology or luminosity to arise from Boo III. 

Idealized N-body simulations are an ideal framework to model the tidal disruption of gassless dwarf spheroidal galaxies. While cosmological simulations may provide a more realistic environment, they struggle to resolve even bright dwarf spheroidals (e.g., ). 

In this work, we revisit the observational properties of Boo III and then model its tidal disruption with N-body simulations to predict what a stellar stream around Boo III should appear as. 

# The observational status of Boötes III

We revisit the structural and line-of-sight velocity properties of the galaxy. @tbl:obs_props shows our recommended observational parameters for Boo III. 





| parameter                | Units              | value                       | Reference          |
| ------------------------ | ------------------ | --------------------------- | ------------------ |
| $\alpha$                 | deg                | $209.47\pm0.11$             | @struct_params     |
| $\delta$                 | deg                | $26.68\pm0.06$              | @struct_params     |
| distance modulus         | mag                | $18.34 \pm 0.19$            | @vivas+2020        |
| distance                 | mag                | $46.56 \pm 4$ kpc           | @vivas+2020        |
| $\mu_\alpha \cos \delta$ | mas yr$^{-1}$      | $-1.16 \pm 0.02 \pm 0.017$  | @battaglia+2022    |
| $\mu_\delta$             | mas yr$^{-1}$      | $-0.88  \pm 0.01 \pm 0.017$ | @battaglia+2022    |
| $\mathrm{v}_\text{los}$  | km s$^{-1}$        | $188 \pm 2.2*$              | Will likely change |
| $\sigma_\mathrm{v}$      | km s$^{-1}$        | $7.7_{-1.5}^{+2.0}*$ km/s   | Will likely change |
| $R_h$*                   | arcmin             | $44_{-6}^{+7}$              | @struct_params     |
| ellipticity              | --                 | $0.42_{-0.14}^{+0.11}$      | @struct_params     |
| position angle           | deg                | $89\pm9$                    | @struct_params     |
| $M_V$                    | mag                | $-5.1\pm0.24$               | @luminosity        |
| $M_\star$                | $\mathrm{M}_\odot$ | $(14\pm3) \times10^3$       | @luminosity        |
| [Fe/H]                   | dex                | $-2.5\pm0.1$                | @kin_params        |

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

Figure: The distribution of *Gaia* stars 



| model       | $\Delta \xi/'$ | $\Delta \eta/'$ | $e$                    | $\theta / \deg$  | $R_h\ /\ '$      | $n_\text{S\'ersic}$ | $N_\text{sat}$ |
| ----------- | -------------- | --------------- | ---------------------- | ---------------- | ---------------- | ------------------- | -------------- |
| Plummer     | $9.0\pm5.7$    | $-7.2 \pm 3.4$  | $0.42_{-0.11}^{+0.14}$ | $88.5\pm9.5$     | $43.8_{-6}^{+7}$ | --                  | $112 \pm 15$   |
| Exponential | $9.4\pm6$      | $-6.1\pm5$      |                        |                  |                  |                     | $130 \pm20$    |
| Sérsic      | $12\pm4$       | $-9_{-3}^{+4}$  | $0.23\pm0.14$          | $93^{+20}_{-19}$ | $67_{-15}^{+23}$ | $2.9_{-0.7}^{+0.9}$ | $57_{-7}^{+8}$ |



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

Our N-body methods are the same as in @boyea+2026, we recount the main details here. 

We adopt the Astropy v4.0 Galactocentric frame. 

For our MW-only model, we use the @EP2020 potential, an analytic approximation of @mcmillan2011. The potential consists of an Hernquist1990 buldge, a mw1995 thin and thick disk, and a NFW halo. 

For simplicity, we do not include the LMC, which does not affect Boo III's recent orbit (see Appendix X for a discussion.)

## Point particle orbits

- Agama: 100,000 orbits by sampling observed properties from Table @tbl:obs_props
- Selected orbits 

  We use Agama to calculate 100,000 orbits in the EP2020 potential. 

  Orbits selected from quantile method, median properties of near $\pm2\sigma$ pericentre. 



| name | Peri  | Apo    | Helio dist | pmra       | pmdec      | vlos     | period |      |
| ---- | ----- | ------ | ---------- | ---------- | ---------- | -------- | ------ | ---- |
|      | (kpc) | (kpc)  | (kpc)      | ($\masyr$) | ($\masyr$) | ($\kms$) | Gyr    |      |
| 1.5  | 1.52  | 74.36  | 39.1       | -1.15      | -0.88      | 188.1    | 0.89   |      |
| 4    | 3.74  | 84.08  | 42.68      | -1.15      | -0.88      | 188      | 1.04   |      |
| mean | 7.1   | 97.55  | 46.56      | -1.16      | -0.88      | 188      | 1.25   |      |
| 12   | 11.86 | 117.48 | 50.79      | -1.17      | -0.88      | 188      | 1.57   |      |
| 18   | 18.06 | 149.38 | 55.5       | -1.17      | -0.88      | 187.9    | 2.13   |      |
| 26   | 26.14 | 210.61 | 60.76      | -1.17      | -0.88      | 187.6    | 3.25   |      |

Table: The mean properties of point orbits for different pericentres. 

## Permissible initial conditions

In $\Lambda$-Cold Dark Matter (\LCDM), dark matter subhaloes (self-bound gravitating overdensities) all follow the NFW relation. In addition, the concentration of a dark matter halo (i.e. scale radius given the halo mass) is related to the halo's mass. We furthermore expect dark matter haloes to be related to stellar masses following the stellar mass--halo mass relation. Consequently, given an estimate of a dwarf's initial mass, we can infer its initial dark matter halo properties. 

From the models that follow, we find it unlikely that Boo III has lost more than about 90% of its stellar mass (unless the initial galaxy was unusually extended).

Based on a selected orbit, we can find which initial dark matter haloes then result in agreement with the present-day velocity dispersion using the tidal track machinery from @EN2021. @EN2021 Eq. 5 describes the relationship between rmax and vmax during tidal evolution, and their eqns. 10-16 predict the time evolution of a halo given the potential circular timescale at pericentre, initial halo circular timescale at rmax, and orbital eccentricity. 

Using the properties of select median orbits, we can thus solve for which halos produce a specified final velocity dispersion after 10 Gyr. We numerically solve for this root via a standard bisection algorithm. 

@fig:initial_halos shows which halos are predicted to produce the present day velocity dispersion $7km/s??$ after 10 Gyr on orbits with pericentres of 7, 12, and 18 kpc (the mean, 84th, and 95th percentile pericentres). The break in the curves occurs because of the slightly different fitting formulae in EN2021 between their modest and heavy mass loss regiemes. 





![image-20260612115749998](/Users/daniel/Library/Application Support/typora-user-images/image-20260612115749998.png)

Figure: The initial rmax and vmax es which EN2021 predicts will acchieve Boo III's velocity dispersion (curves), and the cosmologically expectations fo rBoo III's vmax given its stellar mass (green horizontal line), and Boo III's rmax given a vmax (grey slanted line). **Left:** the required halos matching the final velocity dispersion for different orbits. **Middle** The required halos for different final velocity dispersions on the fiducial 12kpc orbit. **Right:** The required halos for different lengths of orbit on the fiducial 12kpc orbit. In any case, Boo III requires either an unlikely orbit (e.g. 26kpc pericentre or a low number of pericentres) or an unexpectedly massive and concentrated initial halo. 



# Results



## Orbits

@fig:orbit_hist shows a histogram of randomly sampled point-particle orbits of Boo III. Boo III has a pericentre of 4--12 kpc (16-84th quantile) and anywhere from 1.5 -- 26 kpc (3-sigma range). 

@fig:orbits illustrates the range of possible orbits for Boo III. Despite the wide variation, all orbits are radial and planar (projecting to a line in the $x$--$y$ plane). Increasing the pericentre increases the apocentre and orbital period. The orbital pericentre is likely the most important quantity in determining the tidal evolution of Boo III. 

@fig:pericentre_distance illustrates the tight correlation of pericentre with distance among our sampled orbits. With present uncertainties, the distance determines the pericentre. The lower panel shows the relationship between total angular momentum today and distance. The total angular momentum nearly vanishes at heliocentric distances around 35kpc, resulting in a steep dependence of pericentre with assumed distance. This property results from the tangental velocity's dependence on heliocentric distance and since the solar motion nearly cancels the tangental motion at such a distance. Improving the distance measurement to Boo III will be key to improving our understanding of its orbital history. 

Appendix @sec:extra_orbits compares the effects of including a Large Magellanic Cloud and changing the Milky Way potential on the orbit of Boo III. Since Boo III is on the opposite side of the galaxy as the LMC, Boo III is very slightly affected by including the LMC. In addition, the uncertainties in the Milky Way potential are inconsequential next to the pericentre's sensitivity to assumed distance. 



![image-20260429151837975](/Users/daniel/Library/Application Support/typora-user-images/image-20260429151837975.png)

Figure: The distribution of possible pericentres, with the examples considered here marked along the bottom.



![image-20260429143711521](/Users/daniel/Library/Application Support/typora-user-images/image-20260429143711521.png)

Figure: Orbits of Boo III with different pericentres. 



![image-20260429150209385](/Users/daniel/Library/Application Support/typora-user-images/image-20260429150209385.png)

Figure: The pericentre (top) and total angular momentum (bottom) as a function of heliocentric distance. Green points show the sampled MC orbits, and the large dots show the example orbits discussed in the text.







## Tidal evolution

Given the wide range of pericentres, we must assume an orbital history then explore the tidal evolution of select dark matter and stellar components on that tidal history. 

Our fiducial model assumes Boo III follows mean orbit (12 kpc pericentre) for 9? Gyr. To acchieve a sufficiently close velocity dispersion, we require a high initial concentration, 3$\sigma$ from Ludlow+2016, and a high initial maximum circular velocity (30 km/s). Under this case, the final velocity dispersion is still below observations (7ish?? km/s) but is within the 2-sigma range which we find sufficient. We explore other orbital possibilities below.



![image-20260605104424706](/Users/daniel/Library/Application Support/typora-user-images/image-20260605104424706.png)

Figure: The final distribution of dark matter (purple) and stars (white) for the fiducial model in the $y$--$z$ Galactocentric plane. The density scale spans 5 orders of magnitude and is logarithmic??







![image-20260608154613744](/Users/daniel/Library/Application Support/typora-user-images/image-20260608154613744.png)

Figure: The tidal track of the fiducial model in terms of the stellarand dark matter component. The stars track is plotted as $r_h$, $\sqrt3\ \sigma_\text{v}$ (i.e. using the 3D velocity dispersion). As the stellar and DM track eventually coincide, the stellar component is *tidally limited* in this model, causing significant stellar mass loss and setting an upper limit on the extent of the stellar component. 



## Alternate orbits and halos

Tidal evolution is primarily determined by four properties: the pericentre, the orbital time (including eccentricity...), the halo mass, and the halo concentration. 

While this is a 4-dimensional problem, the present-day velocity dispersion constrains one variable, and assuming that the model survives for 10 Gyr constrains another. In this case, for a given halo scale radius (rmax), there is a unique solution to vmax which evolves to match the present day 

Generally, less concentrated halos and orbits with smaller pericentres should evolve more tidally. 

In any case, Boo III is unlikely to arise 



Another possibility is that Boo III was recently shifted onto an extreme orbit due to a perturber (e.g., s+++). We explore a model which evolves our fiducial halo on one pericentre of 4 kpc!? While such a pericentre exposes the galaxy to strong tidal forces, the short time evolved on such an orbit prevents the galaxy from developing significant stellar mass loss. *We have also explored other orbits not shown using yet more compact halos on tight orbits.* 



![image-20260430082213646](/Users/daniel/Library/Application Support/typora-user-images/image-20260430082213646.png)

Figure: The present-day dark matter and stellar distribution for the *shocked* model which experiences one extreme pericentre of 1.5 kpc. While dark matter is heavily disrupted, the stars remain relatively near the remnant but are stretched into an oblong/spiral shape. 







We may also explore halos which are more or less concentrated than our fiducial model on the same orbit. *TODO* a model with a smaller concentration and mass, e.g., 2.2 kpc, 22 km/s, is more resilient to tides. As we select models near the upper limit of cosmologically-reasonable initial conditions for a ultra-faint dwarf galaxy, we do not expect tides to be more severe than the models presented here. 



## Cored and other halo models. 





## Comparison of present-day properties

In all cases, we select models which roughly agree with the present-day properties of Boo III (half-light radius, velocity dispersion, sky position, distance, and kinematics).





![image-20260430081434832](/Users/daniel/Library/Application Support/typora-user-images/image-20260430081434832.png)



![image-20260430081340817](/Users/daniel/Library/Application Support/typora-user-images/image-20260430081340817.png)

Figure: The final distributions of stars for each model considered above. 



Figure: The final stellar distribution of the cored model. Interestingly, a core of equal size does not strongly affect tidal evolution as the velocity dispersion constrains. 





Finally, we calculate the surface density of stars along the stream axis. 



![image-20260430080733367](/Users/daniel/Library/Application Support/typora-user-images/image-20260430080733367.png)

Figure: The along-stream density 







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






