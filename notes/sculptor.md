# Sculptor

Sculptor (Scl) is one of the first discovered dwarf galaxies of the Milky Way (Shapley 1938; only preceded by the SMC and LMC!). As a classical dwarf spheroidal, Scl is relatively bright and compact. 

Since the initial discovery of Scl, most authors have noted that Scl has a slight ellipticity ($0.36$), often attributed to tidal effects. However, this does not align with the absolute proper motion or orbital path of the galaxy. 

While Scl has a relatively large pericentre (greater than 50kpc), [@Sestito+23] detect that the galaxy likely has an excess of stars in the outskirts (past about 60 arcminutes). This excess is perhaps one of the clearest indications that Scl may be affected by tides. Here, our goal is to determine if under a $\Lambda$ CMD paradigm with DM-only simulations, tidal effects of the Milky Way (or LMC) may indeed be consistent with these observations, or if these observations may reveal instead a extended stellar "halo" or second component of the galaxy—illustrating a complex formation for the galaxy.

# Observations

In the table below, we compile recent observations of Scl to provide the most modern and precise estimates of critical observational properties we then use to build our simulations. 

A brief comparison with other estimates of all properties reveals that most studies agree within 2$\sigma$ for the observed properties. 

Note however, that the proper motions described in MV20a using Gaia DR3 also have a 0.017 mas/yr systematic error (angular scale of 1 deg using Lindegren+2021 eq. 25).

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
- parallax 3$\sigma$ consistency with zero (no offset applied)
- CMD within the polygon in BP - RP versus G described by ...
- total proper motion less than 10 mas/yr (corresponding to a velocity of about 4000 km/s)



The matched filter of stars matching the above cuts. The red points match proper motions.

![matched filter search](figures/scl_matched_filter.pdf)

We also explore the CMD and PM distributions along the orbit and find little indication that there may be any over density of stars or change associated with tidal effects. Deeper surveys with proper motions like Vera Rubin observatory and Euclid will be essential to uncovering fainter features associated with dwarf galaxies.

### DELVE



Select stars from DELVE DR2

```
SELECT TOP 1000
       *
FROM delve_dr2.objects
WHERE 11 < ra
and ra < 19
and -37.7 < dec
and dec < -29.7

```

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

### Milky Way Potential



Following @borukhovetskaya+2022, we define a multi-component MW potential as follows (in code units: 1 kpc, 10^10 Msun)

Next, we integrate the orbit of Sculptor back 10 Gyr and place the halo in a static Milky Way potential described by @EP2021. 



| Component  | Values                    |      |
| ---------- | ------------------------- | ---- |
| thin disk  | M=5.9, a=3.9, b=0.31      |      |
| thick disk | M=2, a=4.4, b=0.92        |      |
| bulge      | M = 2.1, a = 1.3          |      |
| halo       | Mvir=115, r=20.2, c=9.545 |      |



### Haloes

| halo    | $v_{\rm circ,\ max}$ | $r_{\rm max}$ | $\sigma_{vx, \rm best}$ | $h  / {\rm kpc}$ |
| ------- | -------------------- | ------------- | ----------------------- | ---------------- |
| average | 31.1                 | 5.96          | 8?                      | 0.014            |
| compact | 31                   | 3.2           | 9.2?                    | 0.017            |
| heavy   | 40                   | 5.9           |                         |                  |



## Setup

we use Agama to create initial n-body models. We then run models in isolation for ~5Gyr to ensure that the realization is in dynamical equilibrium. Our models are ran using Gadget-4. The softening is taken from @powers2003 and is
$$
h = \frac{4 R_{200}}{\sqrt{N_{200}}}
$$
which for our average halo equates to 0.14 kpc. 



The initial kinematic conditions for the galaxy are determined by integrating the orbit back in time 10 Gyr and determining the time of the next apocentre. We ignore dynamic friction for this calculation.



| Orbit                       | mean     | smallperi | largeperi |
| --------------------------- | -------- | --------- | --------- |
| ra                          | 15.0183  | ''        | ''        |
| dec                         | -33.7186 | ''        | ''        |
| pm ra                       | 0.0990   | 0.134     | 0.064     |
| pm_dec                      | -0.160   | -0.196    | -0.126    |
| heliocentric distance / kpc | 83.2     | 82.7      | 85.4      |
| rv                          | 111.4    | ''        | ''        |
| pericentre / kpc            | 53       | 43        | 63        |
| apocentre / kpc             | 102      | 96        | 114       |
| $x_0$ / kpc                 | 16.1574  | -2.4303   | 4.9579    |
| $y_0$                       | 92.5936  | -43.4768  | 57.7071   |
| $z_0$                       | 39.5549  | 86.0834   | -97.6210  |
| $v_{x,0}$ / km s            | -2.31    | -20.16    | 22.38     |
| $v_{y,0}$                   | -54.26   | -113.91   | 120.00    |
| $v_{z,0}$                   | 129.015  | -59.60    | 72.61     |



### LMC Orbits

Note that we set the model to begin at $T_0 = -4.88895$ Gyr and finish at $T_f=0.2444475$ Gyr (or in code units, 1036.892895 and 51.84464)

| Orbit                       | mean     | smallperi | largeperi |
| --------------------------- | -------- | --------- | --------- |
| ra                          | 15.0183  |           |           |
| dec                         | -33.7186 |           |           |
| pm ra                       | 0.0990   |           |           |
| pm_dec                      | -0.160   |           |           |
| heliocentric distance / kpc | 83.2     |           |           |
| rv                          | 111.4    |           |           |
| pericentre / kpc            | 53       |           |           |
| apocentre / kpc             | 102      |           |           |
| t last peri                 |          |           |           |
| perilmc                     |          |           |           |
| apolmc                      |          |           |           |
| t last apolmc               |          |           |           |
| $x_0$ / kpc                 | 20.0299  | 9.52      |           |
| $y_0$                       | 299.0933 | 236.706   |           |
| $z_0$                       | 100.9998 | 35.35     |           |
| $v_{x,0}$ / km s            | 8.63     | 6.48      |           |
| $v_{y,0}$                   | -50.58   | 11.93     |           |
| $v_{z,0}$                   | 52.99    | 73.08     |           |





# Results 



## Milky Way only



In both models which we run, the break radius (from @peñarrubia++) is between 90 and 110 arcminutes. This is not a satesfactory explanation for the aparent excess of stars. While each model does lose much of the extended dark matter halo, the interior component is rather unaffected, and the model remains in the slightly stripped regieme. 



![image-20241112091720413](/Users/daniel/Library/Application Support/typora-user-images/image-20241112091720413.png)

Density Profiles for the extreme orbit:

![image-20241112091539852](/Users/daniel/Library/Application Support/typora-user-images/image-20241112091539852.png)

![image-20241112091547626](/Users/daniel/Library/Application Support/typora-user-images/image-20241112091547626.png)

# Discussion

