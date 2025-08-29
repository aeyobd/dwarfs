# Radial velocity modeling {#sec:rv_obs} 

## Data selection

For both Sculptor and Ursa Minor, we construct literature samples of radial velocity measurements. We combine these samples with J+24's members to produce RV consistent stars and to compute velocity dispersion, systematic velocities, and test for the appearance of velocity gradients. 

First, we crossmatch all catalogues to J+24 Gaia stars. If a study did not report GaiaDR3 source ID's, we match to the nearest star within 1-3 arcseconds (see REF @tbl:rv_measurements). We combine the mean RV measurement from each study using the inverse-variance weighted mean $\bar v$, standard uncertainty $\delta \bar v$, and (biased) variance $s^2$. We remove stars with significant velocity dispersions as measured between observations in a study or between studies. By using that $\chi^2=\frac{s^2}{\delta \bar v^2}$, we remove stars with a $\chi^2$ larger than the 99.9th percentile of the $\chi^2$ distribution with $N-1$ measurements. This cut typically removes stars with reduced chi-squared values $\tilde\chi^2  = \frac{s^2}{\nu\,\delta \bar v^2}\gtrsim 7$ (since the number of measurements is 1-3 typically). 

Next, we need to correct the coordinate frames for the solar motion and on-sky size of the galaxy. We transform the frame into the galactic standard of rest (GSR). The next step is to account for the slight differences in the direction of each radial velocity. Let the $\hat z$ be the direction from the sun to the dwarf galaxy. Then if $\phi$ is the angular distance between the centre of the galaxy and the individual star, the corrected radial velocity is then
$$
v_z = v_{\rm los, gsr}\cos\phi  - v_{\alpha}\cos\theta \sin\phi - v_\delta \sin\theta\sin\phi
$$
where $v_{\rm los, gsr}$ is the line of sight velocity in the GSR frame, $v_\alpha$ and $v_\delta$ are the tangental velocities in RA and Dec, and $\theta$ is the position angle of the star with respect to the centre of the dwarf. The correction from both effects induces an apparent gradient of about $1.3\,\kmsdeg$ for Sculptor and less for Ursa Minor [see also @WMO2008; @strigari2010]. We add the uncertainty in $v_z$ from the distance uncertainty and velocity dispersion in quadrature to the RV uncertainties for each star. We then use the $v_z$ values for the following modelling, however repeating with uncorrected, heliocentric velocities does not significantly affect the results. 

The combined likelihood, including RV information, becomes
$$
{\cal L} = {\cal L}_{\rm space} {\cal L}_{\rm CMD} {\cal L}_{\rm PM} {\cal L}_{\rm RV}
$$
where we assume that the satellite and background distributions are Gaussian. Specifically, 
$$
\begin{split}
{\cal L}_{\rm RV, sat} &= f\left( \frac{v_i -\mu_{v}}{\sqrt{\sigma_{v}^2 + (\delta v_i)^2}}\right) \\
{\cal L}_{\rm RV, bg} &= f\left( v_i /  \sigma_{\rm halo} \right)
\end{split}
$$
where $f$ is the probability density of a standard normal distribution, $\mu_v$ and $\sigma_v$ are the systemic velocity and dispersion of the satellite, and $\delta v_i$ is the individual measurement uncertainty. Typically, the velocity dispersion will dominate the uncertainty budget here. We assume a halo/background velocity dispersion of a constant $\sigma_{\rm halo} = 100\,\kms$  [e.g. @brown+2010].

Similar to above, we retain stars with the resulting membership probability of greater than 0.2. Because of the additional information from radial velocities, most stars have probabilities close to 1 or 0 so the probability cut is not too significant. 

We assume priors on systematic velocity and velocity dispersion of
$$
\begin{split}
\mu_{v} &= N(0\,\kms, \sigma_{\rm halo}^2) \\ 
\sigma_{v} &= U(0, 20\,\kms)
\end{split}
$$
where $\sigma_{\rm halo} = 100\,{\rm km\,s^{-1}}$ is the velocity dispersion of the MW halo adopted above, a reasonable assumption for dwarfs in orbit around the MW. 

## Results {#sec:rv_results}



![LOS velocity fit to Scl.](figures/scl_umi_rv_fits.pdf)

Figure: Velocity histogram of Scl and UMi in terms of $v_z$ (REF). Orange points are from our crossmatched RV membership sample.

For Sculptor, we combine radial velocity measurements from APOGEE, @sestito+2023a, @tolstoy+2023, and @WMO2009. @tolstoy+2023 and @WMO2009 provide the bulk of the measurements. We find that there is no significant velocity shift in crossmatched stars between catalogues. After crossmatching to high quality Gaia stars and excluding significant stellar velocity dispersions, we have a sample of 1918 members.

We derive a systemic velocity for Sculptor of $111.3\pm0.2\,\kms$with velocity dispersion $9.64\pm0.16\,\kms$. Our values are very consistent with previous work [e.g. @walker+2009, @arroyo-polonio+2024, @battaglia+2008]. See appendix REF for a more detailed comparison between different samples and additional tests.

We detect a moderately significant gradient of $4.3\pm1.3\,\kmsdeg$   at a position angle of $-149_{-13}^{+17}$ degrees (see appendix REF). Several past work has attempted to detect a gradient in Sculptor, but no consensus has been reached. @arroyo-polonio+2024 detect a velocity gradient of $4\pm1.5\,\kmsdeg$ in a similar direction using the @tolstoy+2023 sample, finding inconclusive statistical evidence. They additionally suggest a third chemodynamical component of the galaxy which may bias rotation measurements. @battaglia+2008 also detect a $-7.6_{-2.2}^{+3.0}\,\kmsdeg$ velocity gradient along the major axis, approximately the same direction.  Instead, @strigari2010; @martinez-garcia+2023 detect no significant gradient in Sculptor using @WMO2009 sample. Note that pre-*Gaia* work did not have as strong of a constraint on the proper motion of Scl, which limits conclusions of the intrinsic velocity gradient in Scl.

For UMi, we collect radial velocities from, APOGEE, @sestito+2023b, @pace+2020, and @spencer+2018. We shifted the velocities of @spencer+2018 ($-1.1\,\kms$) and @pace+2020 ($+1.1\,\kms$ ) to reach the same scale. We found 183 crossmatched common stars (passing 3$\sigma$ RV cut, velocity dispersion cut, and PSAT J+24 > 0.2 w/o velocities). Since the median difference in velocities in this crossmatch is about 2.2 km/s, we adopt 1 km/s as the approximate systematic error here. Our final sample includes 831 members.

We derive a mean $-245.8\pm0.3_{\rm stat}\,\kms$  and velocity dispersion of $8.8\pm0.2\,\kms$ for UMi.  This is consistent with @pace+2020 and to a lesser extent with @spencer+2018. We do not find evidence for a velocity gradient, consistent with past work [@pace+2020; @martinez-garcia+2023].

## Discussion and limitations

Our model here is relatively simple. Some things which we note as systematics or limitations:

- Inter-study systematics and biases. While basic crossmatches and a simple velocity shift, combining data from multiple instruments is challenging. This appears to be a minor issue (Sculptor) or is corrected for (Ursa Minor).
- Misrepresentative uncertainties. Inspection of the variances compared to the standard deviations within a study seems to imply that errors are accurately reported. APOGEE notes that their RV uncertainties are known to be underestimates but are a small proportion of our sample.
- Binarity. While not too large of a change for classical dwarfs, this could inflate velocity dispersions of about $9\,\kms$ by about $1\,\kms$[@spencer+2017]. Thus, our measurement is likely slightly inflated given the high binarity fractions measured in these systems [@arroyo-polonio+2023, @spencer+2018]. 
- Multiple populations. Both Sculptor and Ursa Minor likely contain multiple populations [@arroyo-polonio+2024, @pace+2020, @tolstoy+2004]. Since we only model a single population, and each population may have a different extent and velocity dispersion, this could result in biased velocity dispersions. However, it is unclear how to uniquely determine an overall velocity dispersion in a multi-population system.
- Selection effects. RV studies each have their own selection effects, which may affect the resulting dispersion, especially if different populations or regions of the galaxy have different velocities or velocity dispersions. We do not attempt to correct for this.

For both Ursa Minor and Sculptor, we also fit models to only data from individual surveys (see REF). Since the resulting parameters are very similar, we conclude that many of the systematic uncertainties are likely smaller than the present errors or that each large survey has similar biases. 



## Velocity modelling and comparisons

Here, we describe in additional detail, our methods and comparisons for RV modelling between studies.

Savage-Dickey calculated Bayes factor using Silverman-bandwidth KDE smoothed samples from posterior/prior.



|      | Study       | Instrument | Nspec | Nstar | Ngood | Nmemb | $\delta v_{\rm med}$ | $R_{\rm xmatch}$/arcmin |
| ---- | ----------- | ---------- | ----- | ----- | ----- | ----- | -------------------- | ----------------------- |
| Scl  | combined    |            | 8945  | 2280  | 2096  | 1981  | 0.9                  |                         |
|      | tolstoy+23  | FLAMES     | 3311  | 1701  | 1522  | 1482  | 0.65                 | --                      |
|      | sestito+23a | GMOS       | 2     | 2     | 2     | 2     | 13                   | --                      |
|      | walker+09   | MMFS       | 1818  | 1522  | 1417  | 1328  | 1.8                  | 3                       |
|      | APOGEE      | APOGEE     | 5082  | 253   | 170   | 164   | 0.5                  | --                      |
| UMi  | combined    |            | 4714  | 1225  | 1148  | 863   | 2.1                  |                         |
|      | sestito+23b | GRACES     | 5     | 5     | 5     | 5     | 1.8                  | --                      |
|      | pace+20     | DEIMOS     | 1716  | 1538  | 829   | 678   | 2.5                  | 1                       |
|      | spencer+18  | Hectoshell | 1407  | 970   | 596   | 406   | 0.9                  | 2                       |
|      | APOGEE      | APOGEE     | 9500  | 279   | 37    | 67    | 0.6                  | --                      |

Table: Summary of velocity measurements and derived properties. {#tbl:rv_measurements short="Spectroscopic LOS velocity measurements"}

measurement 

| study      | mean            | sigma           | $\partial \log\sigma / \partial \log R$ | $\partial v_z / \partial x$ (km/s/deg) | $\theta_{\rm grad} / ^{\circ}$ | $\log B_2/B_1$ |
| ---------- | --------------- | --------------- | --------------------------------------- | -------------------------------------- | ------------------------------ | -------------- |
| all        |                 |                 |                                         |                                        |                                |                |
|            | $111.3\pm0.2$   | $9.64\pm0.16$   | -                                       | -                                      | -                              | 0              |
|            | $111.3 \pm 0.2$ | $9.61 \pm0.16$  | -                                       | $4.3\pm1.3$                            | $-147_{-13}^{+17}$             | -1.9           |
|            | $111.2\pm0.2$   | $9.66\pm0.16$   | $0.07\pm0.02$                           | -                                      | -                              | -3.4           |
| tolstoy+23 |                 |                 |                                         |                                        |                                |                |
|            | $111.3 \pm 0.3$ | $9.79 \pm 0.18$ | -                                       | -                                      | -                              | 0              |
|            | $111.3\pm0.3$   | $9.77\pm0.19$   | --                                      | $4.3\pm1.4$                            | $-154_{-13}^{+19}$             | -1.3           |
|            | $111.2 \pm 0.3$ | $9.73\pm0.19$   | $0.085 \pm 0.023$                       | --                                     | --                             | -4.6           |
| walker+09  |                 |                 |                                         |                                        |                                |                |
|            | $111.1\pm0.3$   | $9.5\pm0.2$     | --                                      | --                                     | --                             | 0              |
|            | $111.1\pm0.3$   | $9.5\pm0.2$     | -                                       | $5.2_{-1.6}^{+1.7}$                    | $-134_{-16}^{+22}$             | -1.9           |
|            | $111.1\pm0.3$   | $9.6\pm0.2$     | $0.06\pm0.03$                           | --                                     | --                             | +0.3           |
| apogee     |                 |                 |                                         |                                        |                                |                |
|            | $111.2\pm0.7$   | $8.6\pm0.5$     | --                                      | --                                     | --                             | --             |
|            | $111.2\pm0.7$   | $8.5\pm0.5$     | --                                      | $6\pm3$                                | $-126_{-33}^{+45}$             | +0.1           |
|            | $111.1\pm0.7$   | $8.5\pm0.5$     | $0.07\pm0.06$                           | --                                     | --                             | +0.9           |

Table: MCMC fits for  different RV datasets for Sculptor among 3 different models. {#tbl:scl_rv_mcmc short="Sculptor RV fits"}





| study   | mean           | sigma                | $\log bf_{\rm sigma}$ | $\log bf_{\rm grad}$ |
| ------- | -------------- | -------------------- | --------------------- | -------------------- |
| all     | $-245.8\pm0.3$ | $8.8\pm0.2$          | +1.3                  | +0.9                 |
| pace    | $-244.6\pm0.4$ | $9.0\pm0.3$          | +0.3                  | +0.5                 |
| spencer | $-246.9\pm0.4$ | $8.8\pm0.3$          | +1.8                  | -0.06                |
| apogee  | $-245.6\pm1.2$ | $10.0_{-0.8}^{+1.0}$ | +1.0                  | +0.5                 |

Table: MCMC fits for UMi velocity dispersion. {#tbl:umi_rv_mcmc short="Ursa Minor RV fits"}



![Scl velocity sample](figures/scl_rv_scatter.pdf)

Figure: RV members of Sculptor plotted in the tangent plane coloured by corrected velocity difference from mean $v_z - \bar v_z$ . The black ellipse marks the half-light radius in @fig:scl_selection. The black and green arrows mark the proper motion (PM, GSR frame) and derived velocity gradient (rot) vectors (to scale).



![Scl velocity gradient](figures/scl_vel_gradient_scatter.pdf) 

Figure: The corrected LOS velocity along the best fit rotational axis. RV members are black points, the systematic $v_z$ is the horizontal grey line, blue lines represent the (projected) gradient from MCMC samples, and the orange line is a rolling median (with a window size of 50).