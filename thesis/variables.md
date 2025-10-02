# Variables, terminology, and acronyms in this work

| Acronym |                                                            |                      |
| ------- | ---------------------------------------------------------- | -------------------- |
| PM      | Proper motion                                              |                      |
| CMD     | colour-magnitude diagram                                   |                      |
| \LCDM{} | Lambda cold dark matter: prevailing theory of cosmology    |                      |
| CDM     | cold dark matter                                           |                      |
| DM      | dark matter                                                |                      |
| Scl     | Sculptor dwarf spheroidal galaxy                           |                      |
| dSph    | a dwarf spheroidal galaxy                                  |                      |
| UMi     | Ursa Minor dwarf spheroidal galaxy                         |                      |
| MW      | The Milky Way                                              |                      |
| NFW     | The Navarro, Frenk, and White density profile.             | sec:NFW              |
| SMHM    | stellar-mass halo-mass relation                            | sec:galaxy_formation |
| LOS     | Line-of-sight velocity, (sometimes called radial velocity) |                      |
| RA      | Right ascension                                            |                      |
| Dec     | Declination                                                |                      |
| LMC     | Large Magellanic Cloud                                     |                      |
| RGB     | Red Giant Branch Stars                                     | (appendix RV)        |
| HB      | Horizontal branch star                                     | (appendix RV)        |
| MCMC    | Monte-carlo Markov Chain                                   | (appendix RV)        |
| RV      | Radial velocity                                            | (appendix RV)        |



| Terminology             |                                                              |      |
| ----------------------- | ------------------------------------------------------------ | ---- |
| Jacobi radius           |                                                              |      |
| pericentre              |                                                              |      |
| apocentre               |                                                              |      |
| baryon                  | Normal matter (i.e. electrons, protons, neutrons)            |      |
| mass-to-light ratio     | The ratio between the (dynamical) mass and luminosity of a stellar popultion, typically in units of solar masses per solar luminosity |      |
| gravitational softening |                                                              |      |
| N-body                  | Method for self-consistently evolving a gravitational system with a large number of particles (bodies). Often used to simulate, e.g., the gravitational structure and evolution of galaxies. |      |
| distribution function   | The phase-space density of finding a particle occupying a given configuration, i.e. what is the distribution of positions and velocities in a galaxy (or equivalently energy or actions and angles). |      |
| Action                  | Actions are conserved, momentum-like quantities useful for describing orbits. In a spherical potential, the two actions are the total energy and angular momentum. In cylendrical coordinates, the actions are the $z$-angular momentum, the radial action. |      |
| Action angles           |                                                              |      |
| (dark matter) halo      | A self-gravitating spheroidal structure composed of an overdensity of dark matter having collapsed in the early universe |      |
|                         |                                                              |      |
|                         |                                                              |      |
|                         |                                                              |      |
|                         |                                                              |      |
|                         |                                                              |      |
|                         |                                                              |      |
|                         |                                                              |      |



| Variable                  | Definition                                                   | Sections discussed |
| ------------------------- | ------------------------------------------------------------ | ------------------ |
| $\xi$, $\eta$             | Tangent plane coordinates in RA and Dec. Defined to be in arcminutes and such that $\xi, \eta$ represent changes in RA and Dec from the adopted dwarf galaxy centre. |                    |
| $\mu_{\alpha *}$          | Absolute proper motion in RA (corrected, so $\mu_{\alpha *} = \mu_\alpha \cos \delta$) |                    |
| $\mu_\delta$              | Proper motion in Dec                                         |                    |
| $\Sigma$                  | Surface density of stars                                     |                    |
| $\Sigma_h$                | Surface density at the half-light radius $R_h$               |                    |
| $\sigma_v$                | Dwarf galaxy line-of-sight velocity dispersion               |                    |
| $\rho_{\rm DM}$           | Dark matter 3D mass density                                  |                    |
| $\rho_{\rm crit}$         | Universe critical density. From @plank+2018, $\rho_{\rm crit} = ...$ |                    |
| $c$                       | NFW halo concentration, $c=r_{200} / r_s$                    |                    |
| $G$                       | Gaia G-band apparent magnitude                               |                    |
| $G_{\rm BP}  -G_{\rm RP}$ | Gaia BP-RP derived colour                                    |                    |
| ${\cal L}$                | Likelihood for J+24                                          |                    |
| ${\cal L}^{\rm PM}$       | Proper motion likelihood for satellite or Milky Way / background |                    |
| $\cal L^{\rm space}$      | J+24 spatial likelihood                                      |                    |
| ${\cal L}^{\rm CMD}$      | J+24 CMD likelihood                                          |                    |
| ${\cal L}_{\rm sat}$      | Satellite membership likelihood                              |                    |
| ${\cal L}_{\rm MW}$       | Milky Way / background membership likelihood                 |                    |
| $M_{200}$                 | Virial mass of NFW halo. The total mass contained within a region $r_{200}$ enclosing a mean density 200 times the universe density $\rho_{\rm crit}$. |                    |
| $P_{\rm sat}$             | Probability that a star belongs to the satellite from J+24   |                    |
| $R$                       | 2D (elliptical) radius.                                      |                    |
| $R_{\rm ell}$             | The sphericalized-elliptical radius                          |                    |
| $R_h$                     | The half light radius, radius which contains half of the galaxy's stars on the sky. |                    |
| $R_{\rm break}$           | The break radius                                             |                    |
| $r$                       | 3D spherical radius or distance                              |                    |
| $r_J$                     | the Jacobi radius                                            |                    |
| $r_s$                     | Scale radius of NFW halo                                     |                    |
| $r_{\rm max}$             | The radius at which $v_{\rm circ, max}$, the maximum circular velocity is reached |                    |
| $r_{200}$                 | "virial" radius of a spherical halo which encloses a mean density 200 times the universe critical density $\rho_{\rm crit}$ |                    |
| $\vmax$                   | The maximum circular velocity as determined by the spherical mass profile |                    |
| $v_{\rm gsr}'$            | GSR line-of-sight velocity corrected for perspective relative to dwarf galaxy centre. |                    |
| $v_{\rm gsr}$             | GSR line-of-sight velocity                                   |                    |
| $v_{\rm hel}$             | Heliocentric line-of-sight velocity                          |                    |
