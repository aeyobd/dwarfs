| parameter                | value                                                        | Source  |
| ------------------------ | ------------------------------------------------------------ | ------- |
| $\alpha$                 | $ 227.2420 \pm 0.0045$˚                                      | M+18    |
| $\delta$                 | $67.2221 \pm 0.0158$˚                                        | M+18    |
| distance                 | $76 \pm 10$ kpc                                              |         |
| $\mu_\alpha \cos \delta$ | $-0.124 \pm 0.004 \pm 0.017$ mas yr$^{-1}$                   | MV20a   |
| $\mu_\delta$             | $0.078 \pm 0.004_{\rm stat} \pm 0.017_{\rm sys}$ mas yr$^{-1}$ | MV20a   |
| RV                       | $246.9 ± 0.1$ km s$^{-1}$                                    |         |
| $\sigma_v$               |                                                              |         |
| $r_h$                    | $18.2 \pm 0.1$ arcmin                                        | M+18    |
| ell                      | $0.55 \pm 0.1$                                               | M+18    |
| PA                       | $50 \pm 1^\circ$                                             | M+18    |
| $M_V$                    | $-19.03 \pm 0.05$                                            | M+18    |
| $\Upsilon_\star$         | $1.5 \pm 0.3$                                                | assumed |





### Gaia

```
SELECT * FROM gaiadr3.gaia_source 
WHERE 1 = CONTAINS( POINT(227.24, 67.22), CIRCLE(ra, dec, 6)) 
AND ruwe < 1.3 
AND pmra IS NOT NULL 
AND pmdec IS NOT NULL
```

