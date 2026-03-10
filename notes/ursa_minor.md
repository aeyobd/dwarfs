# Ursa Minor



| parameter                | value                                                        | Source    |
| ------------------------ | ------------------------------------------------------------ | --------- |
| $\alpha$                 | $ 227.2420 \pm 0.0045$˚                                      | M+18      |
| $\delta$                 | $67.2221 \pm 0.0158$˚                                        | M+18      |
| distance                 | $70.1 \pm 2.9$ kpc                                           |           |
| $\mu_\alpha \cos \delta$ | $-0.124 \pm 0.004 \pm 0.017$ mas yr$^{-1}$                   | MV20a     |
| $\mu_\delta$             | $0.078 \pm 0.004_{\rm stat} \pm 0.017_{\rm sys}$ mas yr$^{-1}$ | MV20a     |
| RV                       | $-245 ± 1$ km s$^{-1}$                                       | This work |
| $\sigma_v$               | $9 \pm 0.6$                                                  | This work |
| $r_h$                    | $18.2 \pm 0.1$ arcmin                                        | M+18      |
| ell                      | $0.55 \pm 0.1$                                               | M+18      |
| PA                       | $50 \pm 1^\circ$                                             | M+18      |
| $M_V$                    | $-19.03 \pm 0.05$                                            | M+18      |
| $\Upsilon_\star$         | $1.5 \pm 0.3$                                                | assumed   |





### Gaia

```
SELECT * FROM gaiadr3.gaia_source 
WHERE 1 = CONTAINS( POINT(227.24, 67.22), CIRCLE(ra, dec, 6)) 
AND ruwe < 1.3 
AND pmra IS NOT NULL 
AND pmdec IS NOT NULL
```



## SmallPeri V38R4.0 orbit

Orbit 1 (initial point orbit)

- x_i = -16.359333038330078
  y_i = 70.19319152832031
  z_i = 20.192514419555664
  v_x_i = 17.17299215197563
  v_y_i = 36.209068006277086
  v_z_i = -118.09322881698608

iteration 2 (after 1e5 particle simulation)

- x_i = -18.16843032836914
  y_i = 74.99334716796875
  z_i = 20.72937774658203
  v_x_i = 17.401021665334703
  v_y_i = 35.00183802843094
  v_z_i = -117.69138953685761
- dJ_R = 0.36
  dJ_z = 2.08
  dJ_phi = -0.69
  dtheta_R = 0.0
  dtheta_z = 0.01
  dtheta_phi = 0.0



Iteration 3 (1e5 as ell)

- dJ_R = 0.09
  dJ_z = 0.25
  dJ_phi = -0.1
  dtheta_R = 0.0
  dtheta_z = 0.0
  dtheta_phi = 0.0
- x_i = -18.361597061157227
  y_i = 75.73907470703125
  z_i = 20.892595291137695
  v_x_i = 17.448167255520822
  v_y_i = 34.86195577979088
  v_z_i = -117.41109347343445

Iteration 4 (1e6)

dJ_R = -0.07
dJ_z = -0.3
dJ_phi = 0.25
dtheta_R = 0.0
dtheta_z = 0.0
dtheta_phi = 0.0

x_i = -18.069520950317383
y_i = 74.80440521240234
z_i = 20.607799530029297
v_x_i = 17.12949180752039
v_y_i = 34.777782899141314
v_z_i = -117.65929777622223

Iteration 4

- dJ_R = 0.0
  dJ_z = 0.0
  dJ_phi = 0.0
  dtheta_R = -0.22
  dtheta_z = -0.1
  dtheta_phi = 0.09
- x_i = -17.401033401489258
  y_i = 74.5107650756836
  z_i = 21.337574005126953
  v_x_i = 14.276910091936589
  v_y_i = 48.61755864620209
  v_z_i = -114.075355219841
