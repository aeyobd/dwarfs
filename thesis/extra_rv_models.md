# LOS velocity modeling {#sec:rv_obs} 

In this section, we analyze the observed line-of-sight (LOS) velocity distributions for Sculptor (Scl) and Ursa Minor (UMi). Our goal is to understand if the galaxies show evidence of tidal features---for example, a velocity gradient or rising outer velocity dispersion. Our derived systemic velocities and dispersions are consistent with literature. We find weak evidence for a velocity gradient in Scl and a rising outer velocity dispersion. However, Scl's velocity gradient is misaligned with its proper motion, likely inconsistent with a tidal interpretation.  Scl's velocity dispersion also mostly rises within the effective radius, more likely a feature of inner structure than tidal disruption. We find no evidence of a gradient in mean velocity or velocity dispersion in UMi. We conclude that Scl and UMi do not show clear features of tidal disruption given current velocity observations. 

## Data selection

For both Sculptor and Ursa Minor, we construct literature samples of LOS velocity measurements. We combine these samples with J+24's Bayesian likelihoods to produce LOS velocity-informed membership probabilities. We select stars with moderate to large membership probabilities for future analysis.



|      | Study       | Nspec | Nstar | Ngood | Nmemb | $\delta v_{\rm med}$ | $R_{\rm xmatch}$ |
| ---- | ----------- | ----- | ----- | ----- | ----- | -------------------- | ---------------- |
| Scl  | combined    | 8945  | 2280  | 2096  | 1981  | 0.9                  |                  |
|      | tolstoy+23  | 3311  | 1701  | 1522  | 1482  | 0.65                 | --               |
|      | sestito+23a | 2     | 2     | 2     | 2     | 13                   | --               |
|      | walker+09   | 1818  | 1522  | 1417  | 1328  | 1.8                  | 3"               |
|      | APOGEE      | 5082  | 253   | 170   | 164   | 0.5                  | --               |
| UMi  | combined    | 4714  | 1225  | 1148  | 863   | 2.1                  |                  |
|      | sestito+23b | 5     | 5     | 5     | 5     | 1.8                  | --               |
|      | pace+20     | 1716  | 1538  | 829   | 678   | 2.5                  | 1"               |
|      | spencer+18  | 1407  | 970   | 596   | 406   | 0.9                  | 2"               |
|      | APOGEE      | 9500  | 279   | 37    | 67    | 0.6                  | --               |

Table: Summary of velocity measurements and derived properties. {#tbl:rv_measurements short="Spectroscopic LOS velocity measurements"}



First, we crossmatch all catalogues to J+24 Gaia stars. If a study did not report GaiaDR3 source ID's, we match to the nearest star within 1--3 arcseconds (see REF @tbl:rv_measurements). We combine the mean RV measurement from each study using the inverse-variance weighted mean $\bar v$, standard uncertainty $\delta \bar v$, and (biased) variance $s^2$. 

We next remove stars with significant velocity dispersions as measured between observations in a study or between studies, to filter out unreliable measurements or possible spectroscopic binaries. By using that $\chi^2=\frac{s^2}{\delta \bar v^2}$, we remove stars with a $\chi^2$ larger than the 99.9th percentile of the $\chi^2$ distribution with $N-1$ measurements. This cut typically removes stars with reduced chi-squared values $\tilde\chi^2  = \frac{s^2}{\nu\,\delta \bar v^2}\gtrsim 7$ (as a star typically fewer than four velocity measurements).



We correct the velocities for the solar motion and on-sky size of the galaxy. First, we transform the velocities into the galactic standard of rest (GSR), a reference frame centred at the Sun's position but at rest with respect to the Galactic centre. Next, we transform velocities to correct for the apparent gradient induced by the dwarf's proper motion. Let the $\hat z$ be the direction from the sun to the dwarf galaxy. Then if $\phi$ is the angular distance between the centre of the galaxy and the individual star, the corrected radial velocity is then
$$
v_{\rm gsr}' = v_{\rm los, gsr}\cos\phi  - v_{\alpha}\cos\theta \sin\phi - v_\delta \sin\theta\sin\phi
$$ {eq:v_z}
where $v_{\rm los, gsr}$ is the LOS velocity in the GSR frame, $v_\alpha$ and $v_\delta$ are the tangental velocities in RA and Dec in the GSR frame, and $\theta$ is the position angle of the star with respect to the centre of the dwarf. The correction from both effects induces an apparent gradient of about $1.3\,\kmsdeg$ for Sculptor and less for Ursa Minor [see also @WMO2008; @strigari2010]. We add the uncertainty in $v_z$ from the distance uncertainty and velocity dispersion in quadrature to the RV uncertainties for each star. We use the $v_z$ values for the following modelling.

We build on J+24's likelihood framework to add velocity information (see @sec:extra_density). The total likelihood of a star belonging to the satellite or foreground/background population is
$$
{\cal L} = {\cal L}_{\rm space} {\cal L}_{\rm CMD} {\cal L}_{\rm PM} {\cal L}_{\rm los},
$$
where we have added ${\cal L}_{\rm los}$, LOS-velocity likelihood term. We assume that the satellite and background $v_{\rm los, gsr}$ distributions are Gaussian. Specifically, 
$$
\begin{split}
{\cal L}_{\rm los, sat} &= f\left( \frac{v_i -\mu_{v}}{\sqrt{\sigma_{v}^2 + (\delta v_i)^2}}\right) \\
{\cal L}_{\rm los, bg} &= f\left( v_i /  \sigma_{\rm halo} \right)
\end{split}
$$
where $f$ is the probability density of a standard normal distribution, $\mu_v$ and $\sigma_v$ are the systemic velocity and dispersion of the satellite, and $\delta v_i$ is the individual measurement uncertainty. Typically, the velocity dispersion dominates the uncertainty budget. We assume a halo/background velocity dispersion of a constant $\sigma_{\rm halo} = 100\,\kms$  [e.g. @brown+2010].

Similar to above, we retain stars with the resulting membership probability of greater than 0.2. Like in J+24, most stars have probabilities of nearly 1 or 0, so the precise choice of cut marginally affects the resulting sample.

## MCMC modeling

To model for possible gradients. We use Monte-Carlo Markov chains (MCMC) to calculate the posterior distributions of our models of the intrinsic velocity dispersion and possible velocity gradients in Scl and UMi.

The likelihood is based on the distribution of stars with respect to the model's predicted $\mu_\V$ and $\sigma_\V$ at each star's position in the tangent plane $\xi$ and $\eta$. Specifically, the log likelihood is 
$$
\log {\cal L} = \sum_i \frac{1}{\sigma_i} \log f\left(\frac{\V_i - \mu_i}{\sigma_i}\right),
$$
where $f$ is the probability distribution of a standard normal distribution, $\V_i$ is the velocity of the $i$th star, and $\mu_i$ and $\sigma_i$ are the model's predicted mean and velocity dispersion at the $i$th star's position $\xi_i$ and $\eta_i$. 

We consider a single model designed to derive systemic velocity and the velocity dispersion, detect velocity gradients, and detect a rising velocity dispersion. The mean velocity at a given point on the sky, $\mu$, is assumed to be,
$$
\mu(\xi, \eta) = \mu_0 + a\,\xi + b\,\eta
$$

for tangent-plane coordinates $\xi$ and $\eta$, systemtic velocity $\mu_0$, and velocity gradient slopes $a$ and $b$. The velocity dispersion at a given position, $\sigma$, is assumed to depend as a power-law on elliptical radius $R_{\rm ell}$ alone:
$$
\log \sigma = \log \sigma_0 + c\,\log(R_{\rm ell} / R_h)
$$
where $\sigma_0$ is the system's velocity dispersion at $R_h$, and $c$ is the velocity dispersion gradient slope.

Our model has parameters with priors
$$
\mu_{0} \sim N(0, \sigma^2_{\rm halo}) \,\kms\\
\sigma_{0} \sim U(0, 20)\ \kms \\
a \sim N(0, 6^2)\ \kms\,{\rm deg}^{-1}\\
b \sim N(0, 6^2)\ \kms\,{\rm deg}^{-1}\\
c\sim N(0, 0.3^2)\ {\rm dex} \\
$$
where $\sigma_{\rm halo} = 100\,{\rm km\,s^{-1}}$ is the velocity dispersion of the MW halo adopted above, $N(\mu, \sigma^2)$ is a standard normal distribution with mean $\mu$ and variance $\sigma^2$, and $U(l, h)$ is a uniform distribution between $l$ and $h$. 

Bayes factors quantify the relative evidence between Bayesian models. We can use Bayes factors to compare the above model to simpler models without gradients in $\sigma$ or $\mu$. Assuming that a gradient-free model represents the case $a=b=0$ (gradient free) or $c=0$ (fixed dispersion) models, we can quickly calculate the relative Bayes factors without running another model. We use the Savage-Dickey method [@dickey+lientz1970], where the Bayes factor is then the relative density of posterior versus prior samples when $a=b=0$ (gradient free) or $c=0$ (fixed dispersion). In order to calculate these densities, we use a kernel density estimator on the posterior samples with a Silvermann bandwidth in each dimension. 

While this is a more complex model, we find each parameter agrees with models of reduced complexity. 

## Results {#sec:rv_results}



![LOS velocity fit to Scl](figures/scl_umi_rv_fits.pdf)

Figure: Velocity histogram of Scl and UMi in terms of projected-correction GSR LOS velocity ([@eq:v_z]). Black points with error bars are from the crossmatched observed sample, and green lines represent MCMC fits to the velocity dispersion.

\input{rv_table.tex}



![Scl velocity gradient](figures/scl_rv_scatter_gradient.png){#fig:scl_velocity_gradient_scatter} 

Figure: **Top**  members of Sculptor plotted in the tangent plane coloured by corrected velocity difference from mean $v_z - \bar v_z$ . The black ellipse marks the half-light radius in @fig:scl_selection. The black and green arrows mark the proper motion (PM, GSR frame) and derived velocity gradient (rot) vectors (to scale). **Bottom**: The corrected LOS velocity along the best fit rotational axis. RV members are black points, the systematic $v_z$ is the horizontal grey line, blue lines represent the (projected) gradient from MCMC samples, and the orange line is a rolling median (with a window size of 50).







![Possible gradients in the velocity dispersion](figures/sigma_v_gradient.pdf)

Figure: The observed velocity dispersion gradient (black) in 10 equal number bins and the derived slopes from the model fitting. Scl shows moderate evidence for an increasing velocity dispersion with radius. UMi's evidence for a radial gradient is weaker, and Scl's velocity dispersion appears to be flat outside of $R_h$. 



For Sculptor, we combine radial velocity measurements from APOGEE, @sestito+2023a, @tolstoy+2023, and @WMO2009. @tolstoy+2023 and @WMO2009 provide the bulk of the measurements. We find that there is no significant velocity shift in crossmatched stars between catalogues. After crossmatching to high quality Gaia stars and excluding significant stellar velocity dispersions, we have a sample of 1918 members.

We derive a systemic velocity for Sculptor of $111.3\pm0.2\,\kms$with velocity dispersion $9.64\pm0.16\,\kms$. Our values are very consistent with previous work [e.g. @walker+2009, @arroyo-polonio+2024, @battaglia+2008]. 

We detect a moderately significant gradient of $4.3\pm1.3\,\kmsdeg$   at a position angle of $-149_{-13}^{+17}$ degrees. Several past work has attempted to detect a gradient in Sculptor, but no consensus has been reached. @arroyo-polonio+2024 detect a velocity gradient of $4\pm1.5\,\kmsdeg$ in a similar direction using the @tolstoy+2023 sample, finding inconclusive statistical evidence. They additionally suggest a third chemodynamical component of the galaxy which may bias rotation measurements. @battaglia+2008 also detect a $-7.6_{-2.2}^{+3.0}\,\kmsdeg$ velocity gradient along the major axis, approximately the same direction.  Instead, @strigari2010; @martinez-garcia+2023 detect no significant gradient in Sculptor using @WMO2009 sample. Note that pre-*Gaia* work did not have as strong of a constraint on the proper motion of Scl, which limits conclusions of the intrinsic velocity gradient in Scl.

@fig:scl_velocity_gradient_scatter plots the combined velocity sample and the MCMC samples for the velocity gradient. While some samples have a consistent gradient with 0, most samples have a positive velocity gradient, consistent with the rolling median trend. 

The velocity gradient in Scl is misaligned with the proper motion

For UMi, we collect radial velocities from, APOGEE, @sestito+2023b, @pace+2020, and @spencer+2018. We shifted the velocities of @spencer+2018 ($-1.1\,\kms$) and @pace+2020 ($+1.1\,\kms$ ) to reach the same scale. We found 183 crossmatched common stars (passing 3$\sigma$ RV cut, velocity dispersion cut, and PSAT J+24 > 0.2 w/o velocities). Since the median difference in velocities in this crossmatch is about 2.2 km/s, we adopt 1 km/s as the approximate systematic error here. Our final sample includes 831 members.

We derive a mean $-245.8\pm0.3_{\rm stat}\,\kms$  and velocity dispersion of $8.8\pm0.2\,\kms$ for UMi.  This is consistent with @pace+2020 and to a lesser extent with @spencer+2018. We do not find evidence for a velocity gradient, consistent with past work [@pace+2020; @martinez-garcia+2023].

For both Ursa Minor and Sculptor, we also fit models to only data from individual surveys. Since the resulting parameters are very similar, we conclude that many of the systematic uncertainties are likely smaller than the present errors or that each large survey has similar biases. 

## Discussion and limitations

Our model here is relatively simple. Some things which we note as systematics or limitations:

*Inter-study systematics and biases*. While basic crossmatches and a simple velocity shift, combining data from multiple instruments is challenging. This appears to be a minor issue (Sculptor) or is corrected for (Ursa Minor).

*Misrepresentative uncertainties*. Inspection of the variances compared to the standard deviations within a study seems to imply that errors are accurately reported. APOGEE notes that their RV uncertainties are known to be underestimates but are a small proportion of our sample. We find that most samples report internally-consistent uncertainties for multiple observations of the same star. 

*Binarity*. Stars in binary systems can inflate the inferred system's velocity dispersion. For example, a system with an apparent velocity dispersion of $9\,\kms$ may have a true velocity dispersion of $8\,\kms$ [@spencer+2017]. Thus, our measurement is likely slightly inflated given the high binary fractions measured in these systems [@spencer+2018; @arroyo-polonio+2023; @gration+2025]. 

*Multiple populations*. Both Sculptor and Ursa Minor likely contain multiple populations [@arroyo-polonio+2024, @pace+2020, @tolstoy+2004]. Since we only model a single population, and each population may have a different extent and velocity dispersion, this could result in biased velocity dispersions. However, it is unclear how to uniquely determine an overall velocity dispersion in a multi-population system.

*Selection effects*. RV studies each have their own selection effects, which may affect the resulting dispersion, especially if different populations or regions of the galaxy have different velocities or velocity dispersions. We do not attempt to correct for this. Given that samples extend through the half-light radius of each system, we do not expect current samples to miss large features except in the very outskirts. 

Could we observe velocity gradients instead in *Gaia*. Unfortunately, given Scl and UMi's distances ($\gtrsim 70\,\kpc$) *Gaia's* systemic proper motion uncertainty of $0.017\,\masyr$ corresponds to a velocity of $5\,\kms$, making the detection of gradients of amplitude less than $5\,\kms$ challenging. In addition, the typical uncertainties of member stars of Scl and UMi are even larger. While in future data releases, *Gaia* uncertainties may permit study in the internal motion of dwarfs, the current *Gaia* proper motions are not yet precise enough. 

## Summary

In this section, we analyze literature samples of LOS velocity measurements for both Scl and UMi. In each case, we find systemic LOS velocities and velocity dispersions consistent with past work. We detect weak evidence for a velocity gradient in Scl, however the gradient is mis-aligned with the proper motion of Scl so is unlikely to be of tidal origin. In both galaxies, the velocity dispersion profile is consistent with a constant velocity dispersion with radius. The lack of observational evidence for ongoing tidal disruption in the velocity distribution of stars in Scl and UMi further supports our interpretations in the main text.  
