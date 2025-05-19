## (if time) Density profile recovery test

To test if we can indeed recover the density profile we input, we can construct mock dwarf observations using the following framework

- Positions are drawn from a underlying density profile parameterized in terms of elliptical radius (for some centre, position angle, ellipticity)
- Member CMD are drawn from an isochrone with $Fe/H = [-2.5, -2.0, -1.5]$ and age uniform in [12, 8] Gyr ago. Use gaia uncertainties from >>>>>
- Proper motions drawn according to the assumed velocity dispersion and measurement error. Observed as norma distributions assuming PM error with magnitude ~ ...

We inject these fake stars into a region between coordinates ... and ... (nearby to Scl but with no obvious features in Gaia). 

In this section, we test several different density profiles

- Exponential (R_s = 2.5 arcmin) with 2,000 stars and 6,000 stars
- Sersic (n= 2) for an extended profile
- Faint Exponential (300 stars at a distance ..., i.e. Msun = ...)

Using this dataset, we "observe" the structural parameters with an uncertainty (...) and apply the J+24 algorithm to the 