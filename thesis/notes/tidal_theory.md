## Tidal evolution

A galaxy in equilibrium will remain in equilibrium, unless acted upon by an external force. As an example, @fig:lagrange_points illustrates the effective potential ($\Phi_{\rm eff} = \Phi + L^2 / 2mr^2$ for angular momentum $L$, host-satellite radius $r$ and satellite mass $m$). around a Sculptor-like NFW galaxy in a circular orbit around a Milky Way-like galaxy. The effective potential accounts for the centrifugal force in the rotating frame. There are two key saddle (Lagrange) points in @fig:lagrange_points, labeled $L_1$ and $L_2$. Since the trajectory of a bound particle is confined to be within the equipotential surface equal to the particle's total energy, the easiest (requiring the lowest energy) path for a particle to escape the satellites gravitational influence is through $L_1$ and $L_2$. Otherwise, the potential steeply increases. $L_1$ and $L_2$ are co-linear with the satellite and host-origin. In general, particles which have higher energy than the energy at $L_1$ and $L_2$ are most commonly unbound from the galaxy. 

The full evolution of a dwarf galaxy in a tidal follows an approximate sequence.

1. *Mass is lost*. In particular, particles and stars on weakly bound orbits are most likely to be removed by tides. Particles escape through $L_1$ and $L_2$
2. *Steams form*. Unbound mass becomes part of a stream or tidal tails. These particles follow similar orbits but are slightly higher or lower energy (depending on which side the particle escaped from).
3. *Mass redistributes*. Bound mass of the galaxy redistributes. This is visible as a wave of outward moving material, with the outermost material reaching equilibrium last. 
4. *A new equilibrium*. With mass loss, the gravitational potential decreases, resulting in a more compact dark matter halo and stars which adiabadically expand to a larger scale radius.

A tidally affected galaxy contains predictable observational clues to this process. Close to the galaxy, the density profile is flattened as a result of newly unbound material. Further away from the galaxy, this would be visible as a stream. Additionally, a satellite stream contains a velocity gradient along the path. Finally, due to mass loss, a satellite will have a reduced velocity dispersion, larger stellar size, and 

Dwarf galaxies evolve along *tidal tracks*. From N-body simulations across a numerous initial conditions and orbits, NFW halos evolve along a narrow track in terms of maximum circular velocity and radius $v_{\rm max}$ and $r_{\rm max}$ . @EN2021 derive an empirical fit to these tidal tracks, finding that 

$$
\frac{v_{\rm max}}{v_{\rm max, 0}} = 
2^\alpha 
\left(\frac{r_{\rm max}}{r_{\rm max, 0}}\right)^{\beta}\left[1 + \left(r_{\rm max} / r_{\rm max, 0}\right)^2\right]^{-\beta},
$$
where $\alpha=0.4$ and $\beta=0.65$. As illustrated in @fig:tidal_tracks, this formula works for both circular and elliptical orbits and is independent of the initial subhalo size or distance to the host.  

![Lagrange points](/Users/daniel/thesis/figures/lagrange_points.pdf){#fig:lagrange_points}

Figure: The potential contours for a 2-body galaxy system in the rotating frame. The saddle points, $L_1$ and $L_2$, define a characteristic radius inside which the effects of tides are relatively minor, and outside which tidal mass losses are expected to be substantial



![Tidal tracks of dwarf galaxies](/Users/daniel/Library/Application Support/typora-user-images/image-20250715095615423.png){#fig:tidal_tracks width=80%}

Figure: Tidal tracks of dwarf galaxies, the logarithm of maximum circular velocity and radius of relative to the initial conditions for satellites on a variety of orbits. Almost all satellites follow the tidal track suggested by @EN2021. Adapted from fig. 6 of @EN2021.