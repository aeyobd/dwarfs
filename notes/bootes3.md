Among Milky Way satellites, Boötes III (Boo III) is a peculiar galaxy. Boo III is significantly larger than most dwarf galaxies (half-light radius of ~300 parsecs) while also being relatively faint. Boo III also has an orbit which may possibly cause disruption of the galaxy.



Boötes III (Boo3) is one of few suspected tidally-disrupting Milky Way satellites, and is among the largest and most diffuse faint dwarf galaxies known. While having a small Milky Way pericentre (<10 kpc), Boo3 remains spare in the literature and has no resolved stream stars reported. Here, we will characterize the tidal disruption of Boo3 with idealized N-body simulations. Specifically, we investigate the impact of Boo3's orbit, the dwarf's initial structure, and the inner MW potential. We then assess the observational characteristics of disruption in Boo3 under each scenario and predict the properties of a possible stellar stream. 



## History of the BooIII

Boo III was originally discovered in @grillmair2009 through matched filter analysis of SDSS data (along with four other streams). In particular, @grillmair2009 find Boo III within a stream dubbed "Styx" which they claim is the ongoing disruption of the dwarf galaxy. In their paper, the only derived properties are the approximate distance (46 kpc or $18.35\pm0.01$ for distance modulus), the centroid (ra, dec = 209.3˚, 26.8˚), a density profile ($\Sigma \sim r^{-1}$)

@correnti+bellazzini+ferraro2009 follow up the initial detection and detect and analyze a possible population of Blue Horizontal Branch (BHB) stars. From this population,they derive a similar density profile, a centroid of ($209.7\pm1.4, 26.8\pm0.6$ degrees), a distance modulus $18.58\pm0.05 \pm 0.14$. By extrapolating from the number of BHB stars detected, they infer $M_V = -5.8\pm0.5$ and an approximate ellipticity of $0.5$. 

@carlin+2009 present initial spectroscopic followup. They find 20 identified members (out of 300 targeted due to high foreground). Their main results are a high velocity dispersion $\sigma_\text{v} = 14. 0 \pm 3.2$km/s, and the derived systemic velocity $\text{v} = 197.5\pm3.8$km/s. They also derive [Fe/H] $\approx -2.1\pm0.2$, and a high metallicity dispersion. With a high dispersion, they interpret this system as likely undergoing active tidal disruption.

@massari+helmi2018 use *Gaia* DR 2 to derive proper motions for seven ultra faint dwarf galaxies, including Boo III. Their results from 34 possible members are:$\pmra = −1.21 \pm0.13$, $\pmdec =  −0.92 \pm0.17$, covariance = 0.23 (not includingg the 0.035mas/yr systematic).

@carlin+sand2018 perform the first orbital modelling of Boo III. They derive a proper motion of (μα cos d, μδ)=(−1.14, −0.98)±(0.18, 0.20) mas yr−1 based on LOS velocity selected members. By excluding members from @carlin+2009 with inconsistent proper motions, they further re-derive the systematic velocity (197.1\pm3.6) and dispersion(10.7\pm3.5km/s). Based on their orbital analysis, the pericentre is 10-12\pm6 kpc and they assert that Boo III is likely disrupting and is consistent with the position of Styx.

@moskowitz+walker2020: Stellar density profiles and structural parameters for many dwarfs including Boo III. 

@vivas+martinez-vazquez+walker2020: Updated RRL census. derive new distance to Boo III ($18.34  \pm 0.19$) with 7 RRL stars. (One previously known from @seaser+2014). They find two RRL beyond the tidal radius, and mention that the disruption of the galaxy is a known fact. 

@tau+2024 find possible RRL stars nearby Bootes III out to very far distances, but some may have came from Sgr stream. They find 32 total RRL members, but not in line with Styx and with a gradient consistent with Sgr contaimination. 

@jensen+2024: finding Boo III extended structure. However, their structural properties are dubious as the outer component may have a negative amplitude (within one-sigma) and the outer component likely exists in their *background-limited* regime where the algorithm begins to break down. 

@Yang+2025, weird features/morphology in Boo III with new observations (photometry). 



# Orbits of Boo III

## Effects of the MW potential

## Effects of the LMC





# Methodological challenges

The standard methodology to be published in my Scl / UMi paper is left nearly unchanged. 

Unfortunately, Boo III requires a more careful treatment of actions to use with the iterative action adjustment code (appendix of BNJE2026) and Agama fails to reverse map actions. This is likely because Boo III is a very low angular momentum orbit, so requires e.g. the Binney+2026 torus mapping code (Agamab). AGAMAb turns out to not be user friendly, easy to compile, backwards compatible, or well documented. I neglect this step for simplicity.



# Estimation of appropriate orbits

We can use the tidal track formalism in EN2021 to estimate which combinations of halos+orbits result in near-observed final velocity dispersions. 



## Procedure

EN2021 provides various fitting formula we can apply to our use-case. Everything is in terms of $\rmax$ and $\vmax$. The relative final $\rmax$ and $\vmax$ follows a universal pre-defined track. 
$$
\vmax/\vmax_0 = ...
$$
  where $x = \rmax/\rmax_0$. 

Next, the total evolution is split into two regiemes:
$$
\begin{cases}
regieme 1 & Tmax/Torb < 3/4 \\
regieme 2 & Tmax/Torb > 3/4
\end{cases}
$$


# Halo by halo

We evaluate models on a few different criteria:

1. Fit to Boo III (position/velocity, velocity dispersion)
2. Stellar and DM mass loss and tidal evolution
3. Stream generation and morphology
4. Anisotropic stellar components?



## NFW

### v22_r2.2

The mean halo. Only weakly forms tidal tails since can only lose a small amount of dark matter mass and remain consistent with the velocity dispersion. 

### 1e5_v35_r6.9

A much heavier but mean-concentration halo (at $z=0$). 

### v30_r3.0

+2 sigma more concentrated than Ludlow+2016 at $z=0$. 

Only somewhat resilient. Can survive one extreme pericentre of 1.5 kpc! (but with little evidence of tidal disruption), two mean pericentres of 7kpc, or until today at a pericentre of 18 kpc. 

### **v30_r2.2 (fiducial)**

This sits at $3\sigma$ more compact than ludlow+2016 at redshift $z=0$, or near the mean at redshift $z=2$ (about when Boo III infalls into the Milky Way)

The fiducial halo can survive until today 

### v30_r1.0

An extraordinarily compact halo, at the upper limits of what we might consider. At a redshift $z=2$, the $3-\sigma$ most compact halos have a scale radius of about $1\,\kpc$. 

At such a high compactness, this halo does successfully survive until today on the mean orbit with a pericentre of 7 kpc (and matching the present day observed velocity dispersion)). 



## Cored

Cored models interestingly do not fare worse than the cuspy ones, only requiring about a 1$\sigma$ increase in concentration. This is likely a result of the required compactness such that the stellar component survives

### v30_r2.2_c0.1



### v30_r2.2_c1



## Anisotropy  / oblate







# Plots

- [x] Final position/velocity/etc agreement
- [ ] Tidal and stellar tracks for each model (1/2)
- [ ] Final tangent plane distribution of all models 
- [ ] quantitative bound mass loss for each model (1/2)
- [ ] halo constraints plot





# Models to run





- Low final velocity dispersion!?
- More cored models

