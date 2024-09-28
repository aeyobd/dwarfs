This folder contains several small items

- `make_nfw.jl` creates a NFW profile with scale radius and mass 1. This halo defaults to be placed in `halos` folder with a truncation radius of 100. These are the initial conditions for our models.
- `agama.patch` is the patch to apply to gadget4 to enable the use of Agama for the potential calculation. 
- `potentials` contains a variety of gravitational potentials to use for orbit integration and models.
