# Additional results an tests {#sec:extra_results}

In this section, we briefly explore additional simulations testing assumptions not varied in the main text: variations in the halo concentration, cored halos, alternative orbits, anisotropy, and ellipsoidal models.

## Alternative initial conditions

The initial concentration (i.e. choice of $\vmax$ and $\rmax$) affects total dark matter evolution. @fig:tidal_tracks_concentration compares models on the \smallperi{} orbit with different initial halo parameters. The heavier halo has $\vmax=43$ and $\rmax=7$, whereas the lighter halo has $\vmax=25$ and $\rmax=2.5$. While the less concentrated halo loses more mass, the halo evolves to similar final structural parameters, as expected in @EN2021. We note these choices were constrained to match Scl's present-day velocity dispersion.



![Tidal dependence on halo concentration](figures/scl_mw_halo_boundmass.pdf){#fig:tidal_tracks_concentration}

Figure: A comparison of the evolution of different N-body models for different halo concentrations (heavier and lighter halo), and the mean orbit instead of the \smallperi{} orbit. In the left, the maximum circular velocity $\vmax$ of the dark matter halo is plotted as a function of simulation time. In the right, $\vmax$ is instead plotted as a function of $\rmax$. 



We have good observational evidence for some dwarf galaxies to be cored. Even Sculptor has debated claims for having a core [e.g., @amorisco+zavala+deboer2014; @breddels+helmi2013; @battaglia+2008; @walker+2009a; @richardson+fairbairn2014; @agnello+evans2012]. To test if a core matter, we adopt a "cored-NFW" model,
$$
\rho/\rho	_s = \frac{1}{(1+r/r_s)^2 (r_c/r_s + r/r_s)},
$$
where $r_c$ is the core radius. We adopt $r_c = 0.1r_s$ and adopt an extreme initial mass of $\vmax=49\,\kms$ and $\rmax=7.1\,\kpc$.



Sculptor and Ursa Minor may also be highly elliptical in 3D [e.g., @an+koposov2022], especially given their projected shapes. We consider 

We find the evolution of the oblate halo to be nearly identical to the spherical halo, similar to @battaglia+sollima+nipoti2015. 

While we have compared the evolution of orbits in the MW potential for simplicity, because the tidal evolution of Scl is very weak in the LMC model, the differences between models become even more slight.

![Tidal dependence on halo structure](figures/scl_mw_structure_boundmass.pdf){#fig:tidal_tracks_structure}

Figure: Similar to @fig:tidal_tracks_concentration, except testing the effects of including a core, velocity anisotropy, and evolving an oblate halo. While the normalization may differ, the tidal evolution is similar. 





![Ursa Minor tidal dependence on orbit](figures/umi_orbits_boundmass.pdf){#fig:tidal_tracks_umi}

Figure: Similar to @fig:tidal_tracks_concentration, except testing the effects of orbits in the LMC. Most LMC models evolve more weakly than the MW models, however these do not reach full agreement with the present-day position so misrepresent recent tidal evolution. 



![Scl MW impact stellar densities](figures/scl_impact_i_f.pdf)

Figure: The initial and final density profiles for the Sculptor model which passes through the MW previously. The impact does not substantially affect the inner density region. The hints of a tidal stream in this model are misaligned to the proper motion, likely a result of the LMC.





![Ursa Minor predicted stream](/Users/daniel/thesis/figures/umi_sim_stream.pdf){#fig:umi_tidal_stream}

Figure: The properties of the stream around the UMi \smallperi{} orbit with Plummer stars. The panels are all as a function of $\xi'$, the distance along the stream as defined by the current GSR proper motion vector. The top panels show the GSR proper motions in RA and Dec, and the bottom two show the distance and GSR radial velocities. To sample the stream, we randomly draw 100,000 samples from the snapshot based on the stellar weights. A detectible gradient in $\mu_{\alpha*}$ and LOS velocity should be detectible if the stream is tracked across several degrees. 

 

![Scl predicted stream](figures/scl_sim_stream.pdf){#fig:scl_tidal_stream}

Figure: The predicted properties of a tidal tail in the Scl model. 



![Scl predicted stream](figures/scl_mw_impact_stream.pdf){#fig:scl_mw_impact_tidal_stream}

Figure: the predicted properties of the stream of stars around Scl in the MW impact model. Note that, compared to @fig:scl_sim_stream, this model produces a much more scattered and misaligned stream, due to the long time since stream formation and nonlinear effects of the LMC encounter. **todo check 2d velocity gradients in this case**
