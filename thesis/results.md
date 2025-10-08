As discussed in Chapters -@sec:introduction and -@sec:observations, this thesis aims to test whether Galactic tides are responsible for the extended density profiles of Scl and UMi. In this Chapter, we analyze tailored N-body simulations, using the methods described in Chapter [-@sec:methods], to assess the tidal impact of the Galactic potential. To anticipate our main conclusion, we find that tides drive dark matter loss in both systems but leave their compact stellar components largely unaffected. The Large Magellanic Cloud (LMC) has the potential of substantially perturbing Scl's and even UMi's orbit, yet the resulting tidal effects are still too weak to account for the extended outer profiles. Our simulations demonstrate that recent tides are unlikely to have altered the stellar structure of Scl or UMi.

In this Chapter, we consider Scl first, describing tidal effects from the MW on its dark matter and stellar components. Next, we consider how accounting for the LMC may affect our conclusions. We then similarly analyze UMi, considering in turn the dark matter evolution, stellar evolution, and orbital effects of the LMC. 

# Tidal effects on Sculptor

## Evolution of Sculptor's dark matter halo

As a representation of an extreme tidal history, we initially investigate the \smallperi{} orbit, described in @sec:scl_smallperi, chosen to maximize possible effects of Galactic tides.

Scl experiences moderate tidal mass loss after 10 Gyrs of evolution. @fig:scl_sim_images shows the stripping of dark matter and the formation of diffuse streams trailing and leading Sculptor's orbit. Because tidal stripping can be described as a gradual removal of the least bound particles, most mass loss occurs in the outer halo. Instead, the inner regions of the galaxy may be relatively unaffected.

N-body models may deviate from a point-particle trajectory due to dynamical self-friction [e.g., @white1983; @miller+2020]. However, this effect is slight for Scl, which ends near the observed position, without adjusting the initial conditions (the green point in [@fig:scl_sim_images]).

The inner density cusp is tidally resilient. @fig:scl_tidal_track shows the initial and final circular velocity profiles, and the evolution of the maximum circular velocity. The maximum velocity drops from $31\,\kms$  to $22\,\kms$, evolving along the tidal track from @EN2021. The final circular velocity profile resembles the initial with an inner cusp, but has a sharper outer truncation. Quantitatively, the halo loses $>90\%$ of its initial mass (see @tbl:scl_sim_results). However, the inner structure is not expected to be affected, as the Jacobi radius is over 3kpc, outside of the initial and final $\rmax$ (see @tbl:scl_sim_results and @fig:scl_tidal_track). Thus, tides may remove significant amounts of mass, but mostly from the outer halo.

![Sculptor simulation snapshots](figures/scl_sim_images.png){#fig:scl_sim_images}

Figure: Images of the dark matter evolution over a selection of past apocentres and the present day. Limits range from -150 to 150 kpc in the $y$--$z$ ($\sim$orbital) plane, and the colourscale is logarithmic, spanning 5 orders of magnitude between the maximum and minimum values. The green dot represents the final expected position of the galaxy, and the solid and dotted grey curves represent the orbit over one previous or future radial oscillation, respectively.



![Sculptor tidal tracks](figures/scl_tidal_track.png){#fig:scl_tidal_track}

Figure: Dynamical evolution for the \smallperi{} model of Sculptor. Dotted and solid lines show the initial and final circular velocity profiles, and blue and orange lines show the dark matter and stellar (2D exponential) profiles. The points represent the evolution of the maximum circular velocity, and the dotted black line shows the tidal track from @EN2021. To calculate the velocity profiles, unbound particles are iteratively removed, recalculating the potential at each step assuming spherical symmetry.



| Property                                  | random samples  | \smallperi{} |
| ----------------------------------------- | --------------- | ------------ |
| pericentre                                | $53\pm3$        | 42           |
| apocentre                                 | $102\pm3$       | 94.4         |
| time of last pericentre / Gyr             | $-0.45 \pm 0.2$ | -0.47        |
| number of pericentres                     | 6               | 6            |
| Jacobi radius / kpc                       | $4.5 \pm 0.3$   | 3.5          |
| Jacobi radius / arcmin                    | $186\pm12$      | 101          |
| final heliocentric distance / kpc         | $83.2\pm2$      | 81.6         |
| $\V_\textrm{max, f} / \V_\textrm{max, i}$ |                 | 0.695        |
| $r_\textrm{max, f} / r_\textrm{max, i}$   |                 | 0.406        |
| fractional final bound mass               |                 | 0.0893       |

Table: The orbital and dark matter properties for the simulation of Sculptor. The random samples column shows the distributions from point orbits, and the \smallperi{} column contains the results from the N-body simulation. {#tbl:scl_sim_results short="Simulation results for Sculptor's dark matter"}



## Evolution of Sculptor's stars {#sec:scl_sim_stars}

Tides minimally affect the stellar component of Sculptor in the \smallperi{} orbit. In @fig:scl_smallperi_i_f, the projected stellar distribution displays no prominent distortions, and the radial density profile is nearly unchanged. Only at a surface density $\sim10^8$ times fainter than the centre do some faint tidal features emerge. The total stellar mass lost corresponds to $\sim 10$ stars in total (see [@tbl:scl_sim_results])---a formidable challenge with the best of observations.

This result implies that Scl's extended profile cannot be reproduced by Galactic tides operating on an initially exponential profile. The weak effect of tides suggests that the outer profile of Sculptor is innate, and not the result of tidal evolution. We check this assertion by choosing a different initial stellar profile which matches the observed profile and assessing how it evolves on the smallperi orbit. We show in @fig:scl_smallperi_plummer_i_f that a Plummer profile (instead of an exponential) provides an adequate fit to Scl's observed profile. The Plummer model loses more stellar mass and forms more luminous tidal tails. Observations reaching surface densities $\sim10$ times fainter than our data could reveal a stream in this case. Nevertheless, over the radial extent probed by our data, the stellar profile remains nearly unchanged by tidal evolution. 

The Jacobi and break radii further support that tidal effects should not be apparent in the observed stellar component. As calculated for this model (see [@tbl:scl_sim_results; @tbl:scl_sim_stars_results]), the break and Jacobi radii both fall outside of $\sim 100$ arcminutes for either stellar component. Indeed, the stellar component only begins to deviate from an exponential profile around this break radius ([@fig:scl_smallperi_i_f; @fig:scl_smallperi_plummer_i_f]). Since no orbits of Scl produce significantly smaller break or Jacobi radii,  it is unlikely any orbit would produce an observable density excess.

@tbl:scl_sim_stars_results quantifies the evolution of stellar properties. The stellar velocity dispersion decreases by only $\sim1\,\kms$ and the half-light radius expands by $\sim 10\%$. This is consistent with adiabatic expansion due to the reduction of the total mass [e.g., @stucker+2023]. In addition, the break and Jacobi radii are $\gtrsim 100$ arcminutes on the sky—tidal signatures would be beyond the reach of our data. Altogether, Galactic tides negligibly impact Scl's stellar component. 

![Sculptor initial and final density profiles](figures/scl_smallperi_i_f.pdf){#fig:scl_smallperi_i_f}

Figure: The tidal effects on Scl's stellar component, for the \smallperi{} orbit with the fiducial halo and exponential stars with $R_s=0.10\,\kpc$. **Top:** the initial (left) and final (right) 2D projected density of stars on the sky. The solid circle marks $6R_h$, the dotted circle the break radius, and the blue arrow the orbital direction. **Bottom:** The initial (dotted) and final (solid) stellar density profiles as compared to the observed stellar density profile. Arrows mark the half-light ($R_h$), break, and Jacobi radii ([@eq:r_break; @eq:r_jacobi]) . 





![Sculptor Plummer initial and final density profiles](figures/scl_plummer_i_f.pdf){#fig:scl_smallperi_plummer_i_f}

Figure: Similar to [@fig:scl_smallperi_i_f] except for Plummer initial stars with $R_h = 0.20\,\kpc$. While a faint stream may be visible with deeper observations, effects on the stellar profile are minimal. 



| Property                     | Exponential         | Plummer |
| ---------------------------- | ------------------- | ------- |
| $\sigma_{\V, i}\,/\,\kms$    | 9.8                 | 10.7    |
| $\sigma_{\V, f} \,/\,\kms$   | 8.8                 | 9.4     |
| fractional stellar mass loss | $2.1\times 10^{-6}$ | $0.024$ |
| $R_{h, i}\,/\,\kpc$          | 0.169               | 0.202   |
| $R_{h, f}\,/\,\kpc$          | 0.189               | 0.227   |
| break radius / arcmin        | $98$                | $105$   |
| break radius / kpc           | 2.3                 | 2.5     |

Table: The present-day stellar properties for the simulations of Sculptor. In each row, we have the initial stellar velocity dispersion (within 1kpc), the final velocity dispersion, the fraction of stellar mass unbound, the initial half-light radius, the final half-light radius, and the break radius in arcmin and kpc ([@eq:r_break]). {#tbl:scl_sim_stars_results short="Simulation results for Sculptor's stars"}



## Orbital effects of the LMC {#sec:scl_lmc}

The Milky Way isn't the only galaxy in town. Recently, work has shown that the infall of the LMC may substantially affect the Milky Way system [e.g., @erkal+2019; @cautun+2019; @garavito-camargo+2021; @vasiliev2023]. With a mass up to one fifth of the MW [e.g., @penarrubia+2015], the LMC infall affects conclusions about the MW properties and the orbits of satellites [see e.g., @patel+2020; @battaglia+2022; @correamagnus+vasiliev2022]. In this section, we examine how the LMC may affect the orbital history of Sculptor.

We use the `L3M11` model of the MW and LMC potential from @vasiliev2024. The `L3M11` potential is an evolving multipole approximation of an N-body simulation including a live MW and LMC dark matter halo. The potential includes a static MW bulge and disk, evolving MW and LMC halos, and the MW reflex motion. In their simulation, the MW was initially a NFW halo with $r_s=16.5\,$kpc and $M_{\rm 200}= 98.4\times10^{10}\Mo$, and the LMC a NFW halo with $r_s=11.7$ and $M_{200} = 24.6 \times 10^{10} \Mo$.  The total `L3M11` MW mass is lighter than our initial @EP2020 potential. 

The inclusion of the LMC reshapes Scl's orbital history, as shown in @fig:scl_lmc_orbits_effect. In the MW-only potential, Scl's orbit is typical of a long-term MW satellite. However, Scl's closest approach to the LMC $~\sim0.1\,\Gyr$ ago affects the long-term orbit---Scl is inferred to reach an apocentre of nearly $300\,\kpc$. Scl may even be on first infall, depending on the MW and LMC mass. Scl's is orbiting the Milky Way on a similar plane to the LMC, but in the opposite direction---Scl is unlikely to be an LMC satellite. 

Promisingly, the timing of the LMC encounter implies a break radius ($\sim 25'$, from @tbl:scl_lmc_sim_stars) consistent with the beginning of Scl's observed density excess (see @fig:classical_dwarfs_densities, and @sec:data_density_profiles). To probe this further, we select an orbit with the smallest LMC-Scl pericentre (20\,kpc) in the `L3M11` model, consistent with Scl's present-day positions and velocities. The orbit is selected following the procedure in [@sec:orbital_estimation] (with uncertainties doubled). This `LMC-flyby` orbit is integrated back in time $2\,\Gyr$ ago to isolate recent tidal effects. We modify Scl's initial halo to have $\rmax = 2.5\,\kpc$ and  $\vmax = 25\,\kms$, slightly reducing the initial stellar velocity dispersion. @fig:scl_lmc_orbits_effect shows this selected orbit in black and @tbl:orbit_ics records the initial conditions.



| Property                                  | random samples            | \texttt{LMC-flyby} |
| ----------------------------------------- | ------------------------- | ------------------ |
| pericentre                                | $44\pm 3$ ($29 \pm 2$)    | 39 (20)            |
| apocentre                                 | $218 \pm 8$               | --                 |
| time of last pericentre / Gyr             | $-0.38\pm0.01$ (-0.11)    | -0.33 (-0.10)      |
| number of pericentres                     | 1 (2)                     | 1 (1)              |
| Jacobi radius / kpc                       | $4.1\pm0.3$ ($4.5\pm0.2$) | 2.8 (2.6)          |
| Jacobi radius / arcmin                    | $168 \pm 11$ ($186\pm6$)  | 132 (121)          |
| final heliocentric distance / kpc         | $83.2\pm2$                | 72.9               |
| $\V_\textrm{max, f} / \V_\textrm{max, i}$ |                           | 0.928              |
| $r_\textrm{max, f} / r_\textrm{max, i}$   |                           | 0.763              |
| fractional final bound mass               |                           | 0.5402             |

Table: The orbital properties and dark matter evolution for the models including an LMC. Similar to @tbl:scl_sim_results except quantities with respect to the LMC are in parentheses. {#tbl:scl_lmc_sims short="Orbits and results for Sculptor in the MW+LMC potential."}



![Sculptor orbits with LMC](figures/scl_lmc_xyzr_orbits.png){#fig:scl_lmc_orbits_effect width=100%}

Figure: Similar to @fig:scl_orbits except for orbits with (orange) and without (green lines) the inclusion of an LMC (blue line) in the potential. The bottom row additionally shows the distance between Scl and the LMC over time.

## Tidal effects from the LMC

Perhaps surprisingly, the combined tidal effect of the MW and LMC is weaker for Scl than in the MW-only case. @fig:scl_lmc_sim_images shows the dark matter evolution of Scl and the passage of the LMC. With only one MW pericentre, Scl's dark matter is less disrupted than the previous MW-only model. The subsequent LMC passage modifies Scl's orbit but has otherwise little effect. The dark matter structure evolves mildly and $\sim 50\%$ of mass remains bound ([@tbl:scl_lmc_sims]).

Correspondingly, the stellar component is nearly unchanged by the combined MW and LMC tides. @fig:scl_lmc_i_f shows the projected stellar distributions and density profiles of this model. While the break radius is within the observed density profile, tidal effects are too weak to be detectable. Structural properties of the stars similarly evolve little ([@tbl:scl_lmc_sim_stars]). 

The reduced tides of the LMC-including model are likely a result of the altered orbit of Sculptor. Compared to the MW-only \smallperi{} orbit, the `LMC-flyby` model completes fewer orbits and therefore experiences a reduced net tidal effect. And while the instantaneous tidal force from the LMC is larger than the MW, Scl does not experience the LMC tidal field long enough to display disturbances. Furthermore, the Jacobi radius due to the LMC still falls outside the observed density profile (@fig:scl_lmc_i_f), and the MW Jacobi radius is even larger. As a result, tides in an MW and LMC potential are even weaker overall than for the MW-only orbit.



![Sculptor simulation snapshots with LMC](figures/scl_lmc_sim_images.png){#fig:scl_lmc_sim_images}

Figure: Similar to @fig:scl_sim_images for the case where the potential includes an LMC. The current position and path of the LMC are represented by the green dot and line, respectively. We also plot the full orbit (over the past 2Gyr) for both Scl and the LMC, as less than one radial period happens over this time frame. 



![Sculptor initial and final density with LMC](figures/scl_lmc_i_f.pdf){#fig:scl_lmc_i_f}

Figure: Similar to @fig:scl_smallperi_i_f for the \texttt{LMC-flyby} model. The Jacobi and break radii here are calculated with respect to the LMC; the corresponding radii with respect to the MW are larger. With only one MW pericentre and a recent, rapid LMC encounter, tidal forces do not appear to affect the stellar distribution.



| Property                     | Scl: LMC-exponential | LMC-Plummer     |
| ---------------------------- | -------------------- | --------------- |
| $\sigma_{\V, i}\,/\,\kms$    | 9.0                  | 9.4             |
| $\sigma_{\V, f} \,/\,\kms$   | 8.8                  | 9.2             |
| fractional stellar mass loss | $<10^{-12}$          | 0.0013          |
| $R_{h, i}$ / kpc             | 0.186                | 0.201           |
| $R_{h, f}$ / kpc             | 0.189                | 0.205           |
| break radius                 | $78'$, 1.6 kpc       | $81'$, 1.7 kpc  |
| LMC break radius             | $23'$, 0.49 kpc      | $24'$, 0.52 kpc |

Table: Similar to @tbl:scl_sim_stars_results, but for the properties of the stellar components of the \texttt{LMC-flyby} model of Sculptor.  {#tbl:scl_lmc_sim_stars short="Simulation results for Sculptor's stars in the MW+LMC potential"}



## Summary

We find, including only the MW potential, that tides only remove dark matter from the outskirts of Scl. The central cusp and compact stellar distribution are resilient to tides. Any tidal effects would also be well outside the reach of current observations. We have also found that the LMC strongly perturbs Scl's orbit---in this case, Scl may be on first infall. However, with only 1 pericentre each for the LMC and MW, the combined tides are weaker than for our initial model. In either case, we conclude that tides are unlikely to affect Sculptor's stellar component.





# Tidal effects on Ursa Minor

## Evolution of Ursa Minor's dark matter halo

![Ursa Minor simulation snapshots](figures/umi_sim_images.png){#fig:umi_sim_images}

Figure: Similar to @fig:scl_sim_images but for Ursa Minor on the \smallperi{} orbit. Dark matter evolution is more dramatic than for Scl. 



| Property                    | Random orbits    | \smallperi{} |
| --------------------------- | ---------------- | ------------ |
| pericentre                  | $37\pm3$         | 30           |
| apocentre                   | $83 \pm 4$       | 75           |
| time of last pericentre     | $-0.97 \pm 0.07$ | -0.80        |
| number of pericentres       | $\sim 8$         | 8            |
| Jacobi radius / kpc         | $3.7 \pm 0.3$    | 2.9          |
| Jacobi radius / arcmin      | $184 \pm 12$     | 156          |
| final heliocentric distance | $70.1 \pm 3.6$   | 64.7         |
| ${\vmax}_f / {\vmax}_i$     |                  | 0.511        |
| ${\rmax}_f / {\rmax}_i$     |                  | 0.249        |
| fractional dm final mass    |                  | 0.035        |

Table: The present-day properties for Ursa Minor's final dark matter halo. See @tbl:scl_sim_results for details. {#tbl:umi_sim_results short="Simulation results for Ursa Minor's dark matter"}



![Ursa Minor tidal tracks](figures/umi_tidal_track.png){#fig:umi_tidal_track}

Figure: Similar to @fig:scl_tidal_track except for Ursa Minor. Ursa Minor loses substantially more mass than Sculptor. 

The tidal evolution of Ursa Minor is similar to that of Sculptor in the MW-only potential. @Fig:umi_sim_images shows snapshots of the DM evolution. UMi loses significantly more DM mass than Scl, forming substantial dark matter streams encircling the MW several times. 

UMi only retains 3% of its total mass after 9 Gyr (@tbl:umi_sim_results). As a result, the final dark matter component is much smaller than the initial, but still evolves along the predicted tidal track (@fig:umi_tidal_track). Despite the more substantial tidal evolution, the Jacobi radius is still large, lying at $\sim 4\,\kpc$, well beyond the final $\rmax$. 

Because of UMi's mass loss, the orbit deviates substantially from a point orbit. Through our orbit-adjustment procedure in @sec:orbit_corrections, we recover nearly exactly the present-day position of Ursa Minor by changing the initial positions by $20\,\kpc$ and $\sim 9\,\kms$. These adjustments do not significantly affect the qualitative structure or pericentre of the orbit. 

## Evolution of Ursa Minor's stars



| Property                     | smallperi-exp       | smallperi-Plummer   |
| ---------------------------- | ------------------- | ------------------- |
| $\sigma_{\V, i}\,/\,\kms$    | 10.0                | 10.9                |
| $\sigma_{\V, f}\,/\,\kms$    | 8.2                 | 8.5                 |
| fractional stellar mass loss | $0.00015$           | 0.039               |
| $R_{h, i}\,/\,\kpc$          | 0.135               | 0.151               |
| $R_{h, f}\,/\,\kpc$          | 0.169               | 0.191               |
| break radius                 | 197 arcmin, 3.7 kpc | 204 arcmin, 3.8 kpc |

Table: Similar to @tbl:scl_sim_stars_results, the present-day stellar properties for the simulation of Ursa Minor for exponential and Plummer stars. {#tbl:umi_sim_stars_results short="Simulation results for Ursa Minor's stars"}





![Ursa Minor simulated density profiles](figures/umi_smallperi_i_f.pdf){#fig:umi_smallperi_i_f}

Figure: Similar to @fig:scl_smallperi_i_f: the tidal effects on the stellar surface density of Ursa Minor for exponential stars on the \smallperi{} orbit. 



![Ursa Minor Plummer model density](figures/umi_plummer_i_f.pdf){#fig:umi_plummer_i_f}

Figure: Similar to @fig:scl_smallperi_i_f: the tidal effects on the stellar surface density of Ursa Minor for Plummer stars on the \smallperi{} orbit. 



Tidal features in UMi's stellar component are still extremely faint, becoming apparent only outside 100 arcminutes in [@fig:umi_smallperi_i_f]. The observed size and velocity dispersion of Ursa Minor evolve little ([@tbl:umi_sim_stars_results]). For exponential initial conditions, tidal effects are unlikely to ever be observable in the near future. 

The break and Jacobi radii fall well outside the observed stellar profile. Tides would have to be far stronger to affect the observed stellar component. As a result, the minimal tidal evolution of this model is not unexpected. 

As for Scl (@sec:scl_sim_stars), we also consider a model where UMi's stars are initially a more extended Plummer profile, resembling the present-day density profile. The stellar evolution of this Plummer stellar component is similar (@fig:umi_plummer_i_f). However, because there are more loosely-bound stars, the Plummer model loses nearly 4% of its initial stellar mass to tides ([@tbl:umi_sim_stars_results]), and tidal features may be detectable if we measure densities 2 orders of magnitude fainter than our present data. We show the properties of a stream in the Appendix (@fig:umi_tidal_stream), but such a stream is unlikely to be observable in the near future. We conclude that tides do not strongly affect the stellar component of this model.



## Effects of the LMC 

![Ursa Minor orbits with LMC](figures/umi_lmc_xyzr_orbits.png){#fig:umi_orbits_lmc}

Figure: Orbits of Ursa Minor with (orange) and without (green) an LMC. The final positions of Ursa Minor and the LMC are plotted as scatter points, and the solid blue line represents the LMC trajectory. Note that the LMC mostly increases Ursa Minor's pericentres and apocentres. 

[@fig:umi_orbits_lmc] shows the effects of including an LMC on the orbit of Ursa Minor. Predominantly, the effect is to increase the orbital period, apocentre, and pericentre. Yet, the orbit remains in a similar plane and with similar morphology. As UMi is on the opposite side of the Galaxy of the LMC, and has a closest LMC approach of $\gtrsim 100\,\kpc$, this is not surprising. 

The deviation from the MW-only orbit is mostly due to the LMC-induced reflex motion of the Milky Way. Because the MW centre is accelerated towards the LMC and away from UMi, UMi's orbit increases in characteristic radius. 



## Summary

While tides affect UMi more strongly than Scl, the tidal effects are insufficient to reshape the observed stellar density profile. Faint tidal tails may be observable with deeper data. Finally, including the LMC in the potential further weakens the tides experienced by UMi.



# Modelling uncertainties

## Halo structure

As the above results show, tides only marginally affected the stellar components of Scl and UMi, even with orbits chosen to have the smallest observationally-consistent pericentres. While we have only presented select models, alternative initial conditions do not affect our qualitative conclusions on tidal effects expected for Scl and UMi in the Galactic (and LMC) potential.

Although our analysis neglects baryonic physics, Scl and UMi have predominantly stars older than $\sim 9$ Gyr [@carrera+2002; @deboer+2011; @weisz+2014; @delosreyes+2022; @sato+2025]. So, gas dynamics are unlikely to affect recent evolution. A collisionless dark-matter-only simulation should therefore be an excellent approximation.

Cored or less concentrated dark matter halos disrupt quicker [e.g., @stucker+2023]. Our fiducial UMi halo, in particular, is among the least concentrated halos consistent with UMi's velocity dispersion. Although Scl's fiducial halo is more concentrated, less concentrated and cored halos evolve similarly (see Appendix -@sec:extra_results). 

Galaxies are rarely perfect isotropic spheres. Sculptor and Ursa Minor are elliptical, and halos are expected to be radially anisotropic [e.g., @navarro+2010]. We test non-spherical and anisotropic models in Appendix -@sec:extra_results, finding that these assumptions likely do not alter our conclusions. 

While alternative initial conditions may influence the total mass evolution, they should produce a similar final stellar structure. A system's observed velocity dispersion directly constrains the mean density within $R_h$ [e.g., @wolf+2010]. Thus, the tidal force required to disrupt the stellar component does not strongly depend on the inner halo structure (via the Jacobi radius).

## Orbital uncertainties {#sec:scl_umi_orbit_uncert}

The long-term orbits of satellites are uncertain. Analytic Milky Way potentials neglect many unknown details, including triaxiality, halo twisting, mass evolution, and substructure. Due to these inadequacies, calculated orbits may diverge significantly from the true orbits of satellites [e.g., @dsouza+bell2022]. The mass-growth of the Milky Way and dynamical friction imply that orbits were typically less bound in the past (less affected by tides). Orbital energy and angular momentum of subhalos are not conserved in cosmological N-body simulations. Consequently, orbits in analytic potentials may overestimate the pericentre and underestimate the maximum tidal stress [although typically not by enough to change our conclusions, @santistevan+2023; @santistevan+2024]. 

As an example, @fig:scl_orbit_lmc_uncert illustrates how changes to the LMC potential modify the long-term orbital trajectories of Scl and UMi. More than 4 Gyr ago, the orbits of Scl diverge substantially. Some orbits are near apocentres of $\sim 300\,\kpc$ when others approach pericentres as small as $\sim 10\,\kpc$. Ursa Minor's orbit is more stable until the possible previous LMC pericentre. In some cases, Ursa Minor may have been bound to the LMC. 

Motivated by @fig:scl_orbit_lmc_uncert, we examine an extreme pericentre of $4\,\kpc$ of Scl with the MW in Appendix -@sec:extra_results, finding it still insufficient to produce the observed density profile. Regardless, we conclude our simulated orbits represent reasonable extremes for *recent* tidal effects. Past encounters with the LMC are revisited below as a form of "pre-processing" in @sec:stellar_halos.



![Long term orbital uncertainties](/Users/daniel/thesis/figures/scl_lmc_orbits_mass_effect.png){#fig:scl_orbit_lmc_uncert}

Figure: The long-term orbital history of Sculptor (**top**) and UMi (**bottom**) is uncertain. In both panels, light, transparent lines represent randomly sampled orbits of the satellites (following @sec:scl_smallperi) in three different LMC/MW mass models from @vasiliev2024. The LMC orbits are in solid, thick lines of the corresponding colour. The L2M11 has a lighter LMC mass, and the L3M10 model has a lighter MW mass than our fiducial L3M11 LMC model. 



## Summary

While the long-term tidal evolution is unconstrained, we conclude that our models are reasonable representations of recent tidal effects. As a result, recent tides are unlikely to affect the stellar distributions of Sculptor and Ursa Minor.

