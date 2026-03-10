Boötes III (Boo3) is an unusual system. Boo3 is one of few MW satellites wish



Boötes III (Boo3) is one of few suspected tidally-disrupting Milky Way satellites, and is among the largest and most diffuse faint dwarf galaxies known. While having a small Milky Way pericentre (<10 kpc), Boo3 remains spare in the literature and has no resolved stream stars reported. Here, we will characterize the tidal disruption of Boo3 with idealized N-body simulations. Specifically, we investigate the impact of Boo3's orbit, the dwarf's initial structure, and the inner MW potential. We then assess the observational characteristics of disruption in Boo3 under each scenario and predict the properties of a possible stellar stream. 





# Bootes III

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







| parameter                | value                               | Source          |
| ------------------------ | ----------------------------------- | --------------- |
| $\alpha$                 | 209.3                               | Grillmair2009   |
| $\delta$                 | 26.8                                | Grillmair2009   |
| distance                 | $46.5 \pm 0.4 \pm 2_\text{sys}$ kpc | Carlin+sand2018 |
| distance modulus         | 18.34 \pm 0.10                      | Carlin+sand2018 |
| $\mu_\alpha \cos \delta$ | $-1.160 \pm 0.037$ mas yr$^{-1}$    | battaglia+2022  |
| $\mu_\delta$             | $-0.88 \pm 0.037$ mas yr$^{-1}$     | battaglia+2022* |
| RV                       | $197.1\pm3.6$ km/s                  | Carlin+sand2018 |
| $\sigma_v$               | $10.7\pm3.5$ km/s                   | Carlin+sand2018 |
| $r_h$ (sphericalized)    | $27.04\pm$ arcminc                  | MW2020          |
| ell                      | $0.33\pm0.09$                       | MW2020          |
| PA                       | $-81\pm8$                           | MW2020          |
| $M_V$                    | $-5.8\pm0.5$                        | Correnti+2009   |





### Haloes

| halo    | $v_{\rm circ,\ max}$ | $r_{\rm max}$ | $\sigma_{vx, \rm best}$ | $h  / {\rm kpc}$ |
| ------- | -------------------- | ------------- | ----------------------- | ---------------- |
| average | 22                   | 3.9           |                         |                  |
| better  | 35                   | 3.0           | 15?                     |                  |





Orbits

| Orbit      | 1      | 2    | 3    |
| ---------- | ------ | ---- | ---- |
| ra         | 209.3  | ''   | ''   |
| dec        | 26.8   | ''   | ''   |
| dist       | 46.56  |      |      |
| pm ra      | -1.160 |      |      |
| pm_dec     | -0.88  |      |      |
| rv         | 197.5  |      |      |
| pericentre | 7.0    | 4.0  | 10.4 |
| apocentre  | 104.1  |      |      |
| $x_0$      | -10.33 |      |      |
| $y_0$      | -43.98 |      |      |
| $z_0$      | -93.79 |      |      |
| $v_{x,0}$  | 7.03   |      |      |
| $v_{y,0}$  | 29.78  |      |      |
| $v_{z,0}$  | -13.41 |      |      |











# Appendix

DELVE selection

select * 
FROM delve_dr3.coadd_objects
WHERE 
q3c_radial_query(ra, dec, 209.3, 26.8, 3)
