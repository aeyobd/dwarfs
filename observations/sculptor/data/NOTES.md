# Data Files

Many different data sources, here is a short description of what we will need.

## Radial velocity data
Documented in more depth in the `rv_xmatch.jl` notebook.

- `tolstoy+23.fits` This is a large sample of stars observed with VLT/FLAMES. The catalogue is publicly available on CDS as table E.1 (J/A+A/675/A49). Downloaded as a binary fits.
- `walker+09_summary.fits` and `walker+09_observations.fits`. From J/AJ/137/3100 (observation paper). The first table is the table which contains the data for each stars (e.g. ra, dec, galaxy, average velocity) but the second contains the individual observations (Tables 2-5 in the paper.) Downloaded from CDS to binary fits.

### Samples in sestito et al. (2023)
The next five entries are useful to exactly reproduce the plots in Sestito et al. (2023). I received each of these files from Fedrico. APOGEE should be publically available and the two stars in the paper are easy to look up, but DARTS does not seem to be public. However, DART is a subset of the `tolstoy+23.fits` file so I do not necessarily need these for my research project.

- `Sculptor_DARTS_BEST_psat40.csv` contains the Sculptor sample from DARTS with a PSAT > 0.40 as crossmatch from J+24.
- `Sculptor_DARTS_BEST_psat40_notAPOGEE.csv` same as above but removing apogee crossmatched stars
- `sculptor_apogeeDR17_xmatch.csv` contains the crossmatch between J+24 sculptor members and APOGEE DR17 RV measurements
- `sculptor_apogeeDR17_xmatch_notDART.csv` same as above but removing DARTS crossmatched stars
- `sestito+23_gmos.csv` contains the two stars that are followed up in the paper

## Gaia samples

- `battaglia+22_sculptor.fits` the data from Battaglia+22 for Sculptor, very similar to J+24 below.
- `gaia_4deg_cen.fits` is simply every stars within 4 degrees of sculptor on the sky downloaded from Gaia's ASQL tool.
- `gaia_6deg_ruwe.fits` is simply every stars within 6 degrees of the centre of Scl additionally imposing a ruwe cut so the size is not huge. The ASQL query is:
```sql
SELECT * FROM gaiadr3.gaia_source
where 1 = CONTAINS(
    POINT(15.0183, -33.7186),
    CIRCLE( ra, dec, 6)
    )
AND ruwe < 1.3
AND pmra IS NOT NULL
AND pmdec IS NOT NULL
```

### Jensen + 2024

- `jensen+24_1_comp_circ.fits`
- `jensen+24_2_comp_circ.fits`
- `jensen+24_2_comp_ell.fits`

These files are based on GaiaDR3, but applying J+24 membership selection to the stars. Ask J+24 for data access.


### Sgr stream
As a mini investigation, I looked into an excess of stars in proper motion space in the Gaia data near Scl. 
It turns out that this is likely just part of the Sgr stream. To double check this, I took a gaia cut centred on the stream closest to my Sculptor fields. This is stored in the file `sgr_stream_near_scl-result.fits` and the ASQL query is:

```
SELECT * FROM gaiadr3.gaia_source where 1 = CONTAINS( POINT(5.28, -15.62), CIRCLE( ra, dec, 4) ) AND ruwe < 1.3 AND pmra IS NOT NULL AND pmdec IS NOT NULL
```

