# Sculptor

Sculptor (Scl) is one of the first discovered dwarf galaxies of the Milky Way (Shapley 1938; only preceded by the SMC and LMC!). As a classical dwarf spheroidal, Scl is relatively bright and compact. 

Since the initial discovery of Scl, most authors have noted that Scl has a slight ellipticity ($0.36$), often attributed to tidal effects. However, this does not align with the absolute proper motion or orbital path of the galaxy. 

While Scl has a relatively large pericentre (greater than 50kpc), [@Sestito+23] detect that the galaxy likely has an excess of stars in the outskirts (past about 60 arcminutes). This excess is perhaps one of the clearest indications that Scl may be affected by tides. Here, our goal is to determine if under a $\Lambda$ CMD paradigm with DM-only simulations, tidal effects of the Milky Way (or LMC) may indeed be consistent with these observations, or if these observations may reveal instead a extended stellar "halo" or second component of the galaxy—illustrating a complex formation for the galaxy.

# Observations

In the table below, we compile recent observations of Scl to provide the most modern and precise estimates of critical observational properties we then use to build our simulations. 

A brief comparison with other estimates of all properties reveals that most studies agree within 2$\sigma$ for the observed properties. 



| parameter                | value                            | Source    |
| ------------------------ | -------------------------------- | --------- |
| $\alpha$                 | $15.0183 \pm 0.0012$˚            | M+18      |
| $\delta$                 | $-33.7186 \pm 0.00072$˚          | M+18      |
| distance                 | $83.2 \pm 2$ kpc                 | Tran+22   |
| $\mu_\alpha \cos \delta$ | $0.099 \pm 0.002$ mas yr$^{-1}$  | MV20a     |
| $\mu_\delta$             | $-0.160 \pm 0.002$ mas yr$^{-1}$ | MV20a     |
| RV                       | $111.03 \pm 0.23$                | This work |
| $\sigma_v$               | $9.61\pm0.16$                    | This work |
| $r_h$                    | $12.33 \pm 0.05$ arcmin          | MV20*     |
| ell                      | $0.36 \pm 0.01$                  | M+18      |
| PA                       | $92\pm1$                         | M+18      |
| $M_V$                    | $-10.82\pm0.14$                  | M+18      |
| $\Upsilon_\star$         | $1.5 \pm 0.3$                    | assumed   |

In the GSR frame, the properties are

- $v_{\rm los} = 75.84 \pm 0.23$
- $\mu_{\alpha *} = -0.275 \pm 0.009$
- $\mu_\delta=0.332 \pm 0.012$

which means the orbit has a position angle of approximately -39.63 degrees on the sky (verified by simulations.)



![](figures/king_profile_fit.pdf)

## Tidal signatures in Gaia?

From the simulations (below), we expect a shift in the proper motions of order 0.05mas/yr, a shift in the radial velocities proportional to distance (xxx km/s/deg) as we move along the orbit. As a model-independent initial assumption, we can assume that the orbit traces along the system proper motion in GSR, so $\mu_\alpha, \mu_\delta = (c, c)$. If there are indeed tidal signatures in the Gaia observations, then we expect an overdensity of stars broadly consistent with Scl's CMD and proper motions and parallax. Therefore, we create a sample of Gaia stars with the cuts

- no NANs in PM, parallax, magnitudes
- RUWE < 1.4 (Gaia quality cut)
- parallax 3$\sigma$ consistency with sculptors distance (effectively same cut as consistent with zero)
- CMD within the polygon in BP - RP versus G described by  [(-0.21, 20.58), (-0.1, 20.1), (0.28, 19.77), (0.7, 19.64), (1.02, 18.5), (1.14, 17.59), (1.41, 16.73), (1.68, 15.94), (1.94, 15.94), (1.91, 16.8), (1.68, 17.01), (1.45, 17.87), (1.32, 19.16), (1.26, 20.1), (1.35, 21.0), (0.46, 21.0), (0.82, 20.27), (0.39, 20.33), (0.16, 20.76)]
- total proper motion less than 10 mas/yr (corresponding to a velocity of about 4000 km/s)



The matched filter of stars matching the above cuts. The red points match proper motions.

![matched filter search](figures/scl_matched_filter.pdf)

We also explore the CMD and PM distributions along the orbit and find little indication that there may be any over density of stars or change associated with tidal effects. Deeper surveys with proper motions like Vera Rubin observatory and Euclid will be essential to uncovering fainter features associated with dwarf galaxies.

## Radial velocity compilation

From now on, we consider members of sculptor to be stars in J+24's catalogue with a maximum PSAT > 0.2 (including 2component circular, elliptical and a 1 component spatial models.)

Note on purity: 97% of stars have consistent RV measurements with our adopted mean (accounting for intrinsic velocity dispersion and all.)

We present a combined sample crossmatching the studies [@sestito+23], apogee, [@walker+2009], and @tolstoy+23, containing a total of 2466 stars with measurements. Note that @battaglia+2009 is superceded by @tolstoy+23. If multiple samples measure the same star, then we take the weighted mean of the measurements. We flag stars with a measurement dispersion greater than 5$\sigma$ of the combined measurement error and stars outside of the range 60-150 km/s as these are likely too far off to be members (and the intrinsic velocity dispersion and systematic errors are the dominant contributors to dispersion.)



While the line of sight velocity dispersion is slightly different from the 1 dimensional velocity dispersion, the correction at 1 degree is $\sim 3\times 9\,\sin(2^\circ) = 1 $km/s.

Weak detection of RV gradient

## Derived properties

First, we calculate the stellar mass of Sculptor using our assumed M/L ratio.

Next, we use the emperical fit from @fattahi2018 which provides a stellar mass to maximum circular velocity relationship. Given the maximum circular velocity, we then calculate the radius of maximum circular velocity by solving for the radius using the mass-concentration relationship from @ludlow2016.

We assume 0.1dex uncertainties in each relationship, and draw random samples to visualize the approximate expected distribution of present-day observationally-consistent halos (figure below).

![](figures/r_max_v_max_in.pdf)



We assume that the anisotropy is zero ($\beta=0$). If the anisotropy is nonzero, than the real velocity dispersion may separate from 

| parameter           | value                                         |
| ------------------- | --------------------------------------------- |
| $L_\star$           | $1.82_{-0.22}^{+0.25}\times10^6\ L_\odot$     |
| $M_\star$           | $2.7_{-0.6}^{+0.7} \times10^6\ {\rm M}_\odot$ |
| $M_{200}$           | $0.48_{-0.25}^{+0.52}\ M_0$                   |
| $c_{\rm NFW}$       | $13.1_{-2.8}^{+3.6}$                          |
| $v_{\rm circ, max}$ | $31\pm8$                                      |
| $r_{\rm circ, max}$ | $5.9 \pm 2.9$                                 |

So to choose a halo, first we pick a v circ max. Then, the value for r circ max is given from the above relationship with a scatter of 0.1 dex. This is calculated by solve_rmax with deltac=0.1 * number of sigma we would like to deviate from the trend.

Methods

## Observational frames

We follow the Astropy v4.0 definition of the Galactocentric frame (i.e....)

We also work in a Galactic Standard of Rest (GSR) frame, which is identitical to the ICRS frame but with the assumed solar velocity subtracted from velocities.





## Orbital analysis

Given the uncertainties in the measurement errors, we can sample present-day measurments from these values, and then integrate the orbits to determine the pericentres. 



![](figures/peri_mc_orbits_corr.pdf)

# Methods

### Haloes

| halo    | $v_{\rm circ,\ max}$ | $r_{\rm max}$ | $\sigma_{vx, \rm best}$ | $h  / {\rm kpc}$ |
| ------- | -------------------- | ------------- | ----------------------- | ---------------- |
| average | 31.1                 | 5.96          | 9.25                    | 0.044            |
| compact | 31                   | 3.2           |                         | 0.024            |
| heavy?  | 47                   | 5.4           |                         |                  |
| both?   | 39                   | 4.5           |                         |                  |

Last column: present-day velocity dispersion for best visual fit.



## Setup

we use Agama to create initial n-body models. We then run models in isolation for ~5Gyr to ensure that the realization is in dynamical equilibrium. Our models are ran using Gadget-4. The softening is taken from @powers2003 and is
$$
h = \frac{4 R_{200}}{\sqrt{N_{200}}}
$$
which for our average halo equates to 0.14 kpc. 



Next, we integrate the orbit of Sculptor back 10 Gyr and place the halo in a static Milky Way potential described by @EP2021. 

The initial kinematic conditions for the galaxy are determined by integrating the orbit back in time 10 Gyr and determining the time of the next apocentre. We ignore dynamic friction for this calculation.



| Orbit      | 1        | 2       | 3     |
| ---------- | -------- | ------- | ----- |
| ra         | 15.0183  | ''      | ''    |
| dec        | -33.7186 | ''      | ''    |
| pm ra      | 0.0990   | 0.1007  |       |
| pm_dec     | -0.160   | -0.161  |       |
| dist       | 83.2     | 79.6    |       |
| rv         | 111.4    | 111.44  |       |
| pericentre | 53.0     | 51.5    | 54.5  |
| apocentre  | 102.0    | 97.8    | 106.0 |
| $x_0$      | 16.1574  | 15.4134 |       |
| $y_0$      | 92.5936  | 86.16   |       |
| $z_0$      | 39.5549  | 43.59   |       |
| $v_{x,0}$  | -2.31    | -3.12   |       |
| $v_{y,0}$  | -54.26   | -63.40  |       |
| $v_{z,0}$  | 129.015  | 126.7   |       |







# Results & Discussion







