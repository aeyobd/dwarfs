This directory contains the code to generate all figures and tables in the Scl/UMi paper. 

Each figure is as follows (arXiv version):
1. classical_dwarf_densities.jl makes `classical_dwarf_profiles.pdf`
2. initial_vcirc.jl makes initial_velocity_nosigma.pdf
3. orbits.jl makes orbits.pdf
4. sims_pictures_zoom.jl makes umi_zoom_image.png, scl_zoom_image.png, scl_lmc_zoom_image.png
5. density_comparison.jl makes density_i_f_w_lmc.pdf
6. classical_rho_mean.jl makes rho_mean_pericentre.pdf
7. metallicity_gradient.jl makes scl_umi_fe_h_gradient.pdf

appendix figures
8. action_corrections.jl makes action_corrections.pdf
9. tidal_track_comparison.jl makes extra_tidal_tracks.pdf
10. multipop_literature_comparison.jl makes
11. made with step 4


Tables
1. read from /observations
2. table_orbit_ics.jl
3.  table_nbody_ics.jl, and table_stars_ics.jl
4. from actions_corrections.jl
5. table_2exp.jl
6. table_obs_props.jl
A. table_stellar_evo.jl (useful for discussion)


Utilities
0. paper_style.jl (basic formatting)
1. utils.jl
2. model_utils.jl
3. mw_isodensity.jl calculates the MW isocontour for the sims_pictures_zoom.jl notebook

