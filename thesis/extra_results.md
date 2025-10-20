# Additional Simulation Results {#sec:extra_results}

In this section, we briefly explore additional simulations testing variations in the halo concentration, cored halos, alternative orbits, anisotropies, and ellipticities. While each of these variables influences total dark matter evolution, the resulting inner structure is similar across models. We conclude that consideration of these effects likely would not change our conclusions.

## Alternative initial conditions

The following models aim to reproduce the stellar velocity dispersion of Sculptor to within $\lesssim 1\,\kms$ with a similar present-day half-light radius.

### Halo concentration

Changing the initial concentration (or $\vmax$ and $\rmax$) primarily affects total dark matter evolution. @fig:tidal_tracks_concentration compares models on the \smallperi{} orbit with different initial halo parameters. The heavier halo has $\vmax=43$ and $\rmax=7$, whereas the lighter halo has $\vmax=25$ and $\rmax=2.5$. While the less concentrated halo loses more mass, the halo evolves to similar final structural parameters. 

The heavier halo diverges from the initial point orbit more substantially. We show the results of correcting the orbit (akin to Ursa Minor's correction, see @sec:orbit_corrections) in @fig:tidal_tracks_concentration as the "heavier, new orbit" halo. Since correcting the orbit has a relatively small effect on the tidal evolution, we neglect these corrections in further comparisons for simplicity. 



![Tidal dependence on halo concentration](figures/scl_mw_halo_boundmass.pdf){#fig:tidal_tracks_concentration}

Figure: A comparison of the evolution of different N-body models for different halo concentrations (heavier and lighter halo), and the heavy halo with the action-angle corrected orbit. **Left:** The maximum circular velocity $\vmax$ of the dark matter halo is plotted as a function of simulation time. **Right:** $\vmax$ is instead plotted as a function of $\rmax$. 

### Dark matter cores

Many dwarf galaxies appear to have dark matter cores. In Scl, the presence or absence of a core has been debated [e.g., @battaglia+2008; @walker+2009a;  @agnello+evans2012; @breddels+helmi2013; @amorisco+zavala+deboer2014; @richardson+fairbairn2014]. To simulate the evolution of a cored dark matter halo, we adopt a "cored-NFW" model,
$$
\rho/\rho	_s = \frac{1}{(1+r/r_s)^2 (r_c/r_s + r/r_s)},
$$
where $r_c$ is the core radius. We set $M_s=0.54\times10^{10}\,\Mo$ where $M_s = 4\pi/3\ r_s^3 \rho_s$, $r_s=1.08\,\kpc$, $r_c=r_s$, and use the same truncation as our fiducial halo (@eq:trunc_nfw). We plot the initial and final density,  as compared to a cuspy NFW, in @fig:cored_i_f.

In @fig:tidal_tracks_structure, we show the evolution of the cored model as compared to a heavier halo. The cored halo appears to evolve more mildly than the NFW halo. This is possibly because, to match the velocity dispersion, the cored halo has $\sim 50\%$ more mass within $1\,\kpc$ than a similar cuspy halo. 

![Tidal evolution of a cored density profile](figures/cored_density_i_f.pdf){#fig:cored_i_f}

Figure: The initial (dotted orange) and final (solid orange) 3D density profiles for the cored model of Scl on the \smallperi{} orbit. The blue thin line represents an NFW halo with the same $M_s$ and $r_s$. 



![Tidal tracks depending on halo substructure](figures/scl_mw_structure_boundmass.pdf){#fig:tidal_tracks_structure}

Figure: Similar to @fig:tidal_tracks_concentration, except testing the effects of including a core, velocity anisotropy, and evolving an oblate halo. While the normalization may differ, the tidal evolution is similar. 

### Velocity anisotropy

Radial velocity anisotropy may cause halos to disrupt faster [e.g. @chiang+bosch+schive2024]. To test the effects of moderate velocity anisotropy, we initialize a model with a velocity anisotropy with an Osipkov-Merritt profile rising from $\beta=0.2$ at the centre to $\beta=1$ and infinity, with scale length $4r_s$. 

@fig:anisotropy_i_f shows the initial and final anisotropy profiles after tidal evolution on Scl's \smallperi{} orbit. The dwarf galaxy becomes more isotropic with tidal evolution. Particles on more radially anisotropic orbits are more easily stripped as they have larger apocentres than more circular orbits of the equivalent energy. Regardless, the overall tidal evolution is very similar to our fiducial, isotropic case (see @fig:tidal_tracks_structure).

![Tidal evolution of anisotropy](figures/anisotropy_i_f.pdf){#fig:anisotropy_i_f}

Figure: The initial and final anisotropy profiles for the initially anisotropic model of Scl on the \smallperi{} orbit. The final profile is more isotropic (closer to $\beta=0$) than the initial.



### An ellipsoidal halo

Sculptor and Ursa Minor may be highly elliptical in 3D, possibly violating the assumption of spherical symmetry [e.g., @an+koposov2022]. Although, the shape of the underlying dark matter halo is unknown.

We create an oblate initial dark matter halo using \agama{}. @fig:oblate_i_f shows the initial equilibrium (after 5Gyr in isolation) isodensity contours of our model. The initial snapshot is sampled from the distribution function 
$$
\begin{split}
f = \frac{M_0}{(2\pi\,J_0)^3} \left[1 + \left(\frac{J_0}{h(J)}\right)^\eta\right]^{\Gamma / \eta} \ \left[1 + \left(\frac{g(J)}{J_0}\right)^\eta\right]^{-B/\eta} \\ \ \exp\left[-\left(\frac{g(J)}{J_{\rm cutoff}} \right)^\zeta\right] 
\end{split},
$$
(`doublePowerLaw` in \agama{}), with $g(J) = g_r J_r + g_z J_z + (3-g_r - g_z) |J_\phi|$ and $h(J) = h_rJ_r + h_zJ_z + (3-h_r-h_z) |J_\phi|$. The `example_doublepowerlaw.py` script in \agama{} solves for the best parameters matching a given density profile. We chose to create a model resembling an NFW but scaled by a factor of 0.5 in the $z$-axis. For a scale-free halo ($r_s=1$, $M_s=1$), the best-fit parameters are $J_0=0.890$, $\Gamma=1.46$, $\eta=0.568$, $B=2.97$, $h_r=0.845$, $h_z=1.66$, $g_r=0.753$, $g_z=1.69$, $M_0=0.965$, $J_{\rm cutoff}=2.40$, $\zeta=20$. The resulting initial conditions are stable in isolation and produce a halo with near the desired density profile and ellipticity. We scale the halo to have a major axis density profile equivalent to an NFW with $\rmax = 7\,\kpc$ and $\vmax=48\,\kpc$.

We find the evolution of the oblate halo to be nearly identical to the spherical halo (see @fig:tidal_tracks_structure), similar to @battaglia+sollima+nipoti2015. The oblate halo becomes spherical after tidal evolution (@fig:oblate_i_f). If Scl's dark matter halo is indeed elliptical today, this model may not be an adequate description. 

![Oblate halo projected density snapshots](figures/oblate_projected_2d.pdf){#fig:oblate_i_f}

Figure: The initial and final (after 9Gyr of tidal evolution) projected density profiles for the oblate halo on Scl's \smallperi{} orbit, projected on with $x'$ and $z'$ the major and minor axes. The contours are drawn assuming normal smoothing of 0.4 kpc and are log-spaced with intervals of 0.1 dex.





### Orbital variation

As discussed in @sec:scl_umi_orbit_uncert, the long-term orbit of Scl is uncertain. For a more massive LMC (e.g., the L3M11 model), Scl may have undergone an extreme pericentric passage with the MW $\sim 6\,\Gyr$ ago. We find that the 3$\sigma$ smallest pericentre is $4\,\kpc$, and simulate this model to test if such a pericentre may be sufficient.

Because of the strong tidal interaction with the MW, the trajectory is substantially perturbed from a point orbit. We adjust the initial conditions by comparing the change in actions (as calculated in the static MW-only potential) before and after the pericentre. Our final model, the `MW impact` model, can approximately reproduce the observed position of Scl (see @fig:scl_mw_impact_orbit). The model has a Galactocentric initial position of $[67.83, -352.2, 110.3]\,\kpc$ and velocity of $[-3.68, 30.79, -22.77]\,\kpc$.

@fig:tidal_tracks_umi compares the tidal evolution of Scl on the mean, `MW impact`, and \smallperi{} orbit. The mean orbit loses less mass than the \smallperi{} model. Instead, the `MW impact` orbit experiences most tidal evolution during its first MW pericentre. While evolving further along the tidal track, the stars of this model nevertheless remain exponential (@fig:scl_mw_impact_i_f). We suggest that the impulsive pericentric passage does not occur for long enough in this model to produce the expected extended density profile. A yet more extreme orbital history would be necessary to tidally transform Scl's stars. 

![Sculptor MW impact orbit](figures/scl_mw_impact_orbits.pdf){#fig:scl_mw_impact_orbit}

Figure: The orbit of Sculptor for the MW-impact model. The point (dotted) and n-body (solid) diverge by $\sim 50\,\kpc$ at early times. 



![Sculptor's tidal evolution for different orbits](figures/scl_orbits_boundmass.pdf){#fig:tidal_tracks_umi}

Figure: Similar to @fig:tidal_tracks_concentration, except testing the effects of orbits in the LMC. Most LMC models evolve more weakly than the MW models.

![Sculptor MW-impact density profiles](figures/scl_impact_i_f.pdf){#fig:scl_mw_impact_i_f}

Figure: Similar to @fig:scl_smallperi_i_f except for the orbit of Scl passing through the MW. 

## The formation of tidal tails

As discussed in @sec:results, any hints of a possible tidal stream around Scl and UMi are beyond the reach of current observational facilities. However, we can still predict the properties of such a stream. 

@fig:scl_tidal_stream and @fig:umi_tidal_stream show the resulting distributions of velocities and distance along the stream orbital axis $\xi'$ for Scl and UMi. The tidal tails may have detectible gradients in radial velocities ($\sim 10\,\kms$) and proper motions (mostly for UMi, of $\sim 0.1\,\masyr$). However, detecting such a gradient would require tracing stars across several degrees on the sky.  

![Sculptor predicted stream](figures/scl_sim_stream.pdf){#fig:scl_tidal_stream}

Figure: The predicted properties of a tidal tail in the Scl model. The panels are all as a function of $\xi'$, the distance along the stream as defined by the current GSR proper motion vector. The top panels show the GSR proper motions in RA and Dec, and the bottom two show the distance and GSR radial velocities. To sample the stream, we randomly draw 100,000 samples from the snapshot based on the stellar weights. A detectible gradient in $\mu_{\alpha*}$ and LOS velocity should be detectible if the stream is tracked across several degrees. 



![Ursa Minor predicted stream](/Users/daniel/thesis/figures/umi_sim_stream.pdf){#fig:umi_tidal_stream}

Figure: The properties of the stream around the UMi \smallperi{} orbit with Plummer stars. 



## Summary

While we have compared the evolution of orbits in the MW potential for simplicity, because the tidal evolution of Scl is very weak in the LMC model, the differences between models become even more slight. And while UMi has a lower pericentre than Scl, the differences between halo structure should apply similarly to UMi.
