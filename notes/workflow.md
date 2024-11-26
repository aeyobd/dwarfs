This document introduces my workflow for running N-body simulations of dwarfs in the MW potential. Below, I discuss each tool used in more detail, but here is the top-level summary

1. Determine simulation orbits:
   1. Look up position, pm, rv, and distance to the object with uncertainties
   2. MC sample from the coordinates/parameters with uncertainties to get a range of possible final positions/velocities
   3. Integrate all the possible positions back in time with Gadget-4
   4. Analyze distribution of pericenters, pick the 16th, 84th percentiles as initial conditionals, along with mean value
   5. Plot these orbits to understand, e.g. # pericenter passages
2. Isolation Run
   1. Find literature estimated velocity dispersion
   2. Use that to estimate total halo mass/compactness
   3. Run a scaled NFW halo in isolation to dynamically relax
3. Orbit Run
   1. Set the output from the isolation run on an orbit from part 1. 
   2. Run under MW potential
4. Analyze results!
   1. See end for helper classes and scripts to make this easier
   2. Rapha has made scripts to assign stellar probabilities to each particle (https://rerrani.github.io/code.html)

## Detailed procedure

As a more detailed procedure, from the start

0. Install Gadget4 (with my agama patch), julia, my julia libraries (Arya, DensityEstimators, LilGuys). See description below
1. Find observational parameters to match for galaxy
   1. Need 6D kinematics
   1. Compile literature sources and adopt reasonably recent and consistent measurements. If scatter is large, add this to combined measurements
2. Chose a MW potential
   1. Here, EP2021

3. Find observationally consistent orbits
   1. Sample observations given uncertainties and integrate back in time to find initial conditions (time, position, and velocity of first apocentre/infall)
   2. One method: use the MCGadget binary (no selfgravity)
   3. Another method: Use agama or galpy to integrate orbit

4. Isolation Run
5. Orbit run
   1. Setup initial conditions
   2. Run Gadget
   3. Post: process (in output directory)
      1. `combine_outputs.py`
      2. `calc_centres.jl`
      3. `combine_outputs.py centres.hdf5` (to add centres to outputs)
      4. `calc_profiles.jl` to calculate 3D DM profiles
      
   4. Compare orbit to measurements and determine final snapshot time closest to today (`sim_analy/analyze_orbit.jl`)
   5. Science analysis: make plots of density profiles, mass evolution, scatter plots, phase space, etc.
   
6. Painting on stars
   1. Once the isolation run (or orbit run) are complete, we can add stars assuming that they are collisionless tracers of dark matter


# Fundamentals

Throughout this work, we use the following unit system:

| Quantity   | Units               |
| ---------- | ------------------- |
| G          | 1                   |
| Mass       | $10^{10}$ M$_\odot$ |
| radius     | 1 kpc               |
| velocities | 207.4 km/s          |
| time       | 4,715,000 yr        |



# Software

## Analysis tools



## Painting stars

As stars are sufficiently spread out in a dwarf galaxy, we can assume that they are simply attached to a dark matter particle, with probabilities governed by the ratio of the distribution functions. 

To find the distribution function, we use Eddington inversion (eq. 4-140b in BT87)
$$
f({\cal E}) = \frac{1}{\sqrt{8}\, \pi^2}\left( \int_0^{\cal E} \frac{d^2\rho}{d\Psi^2} \frac{1}{\sqrt{{\cal E} - \Psi}}\ d\Psi + \frac{1}{\sqrt{\cal E}} \left(\frac{d\rho}{d\Psi}\right)_{\Psi=0} \right)
$$


In practice the right term is zero as $\Psi \to 0$ as $r\to\infty$. For example, if $\rho \propto r^{-n}$ at large $r$ and $\Psi \sim r^{-1}$ then $d\rho / d\Psi \sim r^{-n+1}$ which goes to zero provided that $n > 1$. 

# 

## Agama

I use [Agama](http://agama.software/) for the calculation of the potential. Agama is nice because we also have an LMC N-body potential and because it is written in c++, Agama interfaces with Gadget4 well but also can do orbit integration and action calculations. 

Once the software is cloned, in the directory, 

```
python setup.py install --user
```

will install the library.



I have a custom fork of Agama which lets me calculate time-dependent quantities of potentials, which makes analyzing e.g. LMC potentials significantly easier. 

## Gadget4

See gadget4 documentation  https://wwwmpa.mpa-garching.mpg.de/gadget4/

The primary reason I am using Gadget4 now is that we can use Agama as the potential for the model, which allows for more complex

I find it easiest to store the small changes from gadget4 as a patch. To generate the patch, use

```
git diff --no-prefix > agama.patch
```

and to apply the patch, copy `agama/agama.patch` to the Gadget4 root directory and run in the Gadget directory. 

```
patch -p0 < agama.patch
```



### Building

The source code can be cloned from https://gitlab.mpcdf.mpg.de/vrs/gadget4. I am working from my own fork of the repository. To make, run the make command in the main folder of the repository, e.g.

```bash
make CONFIG=~/path/to/Config.sh DIR=~/path/to/out/
```

in my case

```
make CONFIG=~/dwarfs/gadget/Config.sh DIR=~/dwarfs/gadget
```

where Config.sh is the configfile

### The Configfile

Gadget contains many compile-time settings. For this run we don't need many

I do add in a `EXTERNALGRAVITY_MW` option to the makefile which allows the specification of our four component milky way potential (discussed elsewhere)

### Running

To run gadget with one processer

```bash
./Gadget4 param.txt
```

or for multithread

```
mpirun 
```

```bash
sbatch submit_isolation
```



### Parameterfile

Here is where all the other code options are stored. An overview of the important ones

### Numerical Setup

See @power+2003 for a discussion of this.

The gravitational softening is based on:
$$
h_{grav} = 4 \frac{R_{200}}{\sqrt{N_{200}}}
$$

I in fact find that dividing the softening by a factor of sqrt(10) gives slightly more precise results with a negligible increase in compute time.



Asya gave the similar relation
$$
h = 0.005R_{\rm max} (N/1e6)^{-1/2}
$$
which is $5/4$ of the value from @power+2003. In practice we have been using an additional factor 1/sqrt(10) as this appears to give slightly beter convergence.



# Isolation run

# Orbital runs

Gadget 4 will provide snapshots for each timestep, but now we have to compute useful observational and theoretical summaries of the galaxy's evolution.

## Initial conditions & parameters

Take a snapshot from the end of the isolation run, shift to the first apocentre on the point particle integrated orbit.



### Using Isolation results

I group all of the simulations with the same halo (including scale factors) into the same directory. Paths to the appropriate files for each directory are always set in `paths.sh`. In the case for my orbit simulations, the `iso_paths.sh` in the top level halo directory sets the path to the isolation run and the name of the snapshot to rescale. The parameters of the halo are set in `halo.toml`, which is then rescaled with `rerescale.sh` which calls `rescale_nfw.jl` from LilGuys.

**Note** that if we rescale the isolation run, the softening in the parameter file should be rescaled by the same amount

### Setting up an orbit

From the MCMC orbit analysis, the outputs produce orbit files named like `orbit.csv` which start from a pericentre and evolve to the present day. The orbit can also be manually set. Either way, in the main directory of an orbit, we can rescale the orbit.

To validate the initial conditions, there is a notebook called `notebooks/sim_analy/initial_conditions.jl` which calculates the centre and DM profile to check that these values are what we expect for the model.



## Post processing

We need to first calculate the centres for each snapshot and combine these so that the LilGuys.jl `Output` class can efficiently read and analyze the snapshot.

The first step is to paint the stars onto a snapshot of the isolation run. I have been painting them on around snapshot 20 just to make sure the stars are both in equilibrium but also have time to relax. I, however, calculate the distribution function from the analytic DM profile (the assumed truncated NFW profile), which seems to work a little better than the empirical distribution function. In the `stars` directory of the halo (not the orbit), the script `calc_dist.sh` calculates the distribution function, saving it to `distribution_function.hdf5`, and calculates all of the radii and energies of the snapshot we are painting stars onto, (stored in `energies.hdf5`). Next, we create a directory here for each different stellar profile, containing a `profile.toml` which will be read in from `lguys.load_profile`. The stellar probabilities (and some other diagnostic files) are then calculated when you run `bash paint.sh`. 

To make plotting easier, the 3D stellar profiles can be calculated periodically along the isolation run by running the `add_profiles.sh` script. To validate the accuracy & stability of the painted stellar system, the notebook `compare_stars_2d_isolation.ipynb` in `/analysis` will plot the profiles and velocity dispersion for each calculated profile, and if these seem to be constant in time, we are good to go!

Next, we can apply these stars to an orbit run. I also create a directory called `stars` in each orbit run analysis directory. The script I have called `project.sh` in these directories will project the stars from the first (and adopted final) snapshot onto the sky (in 2d) to be analyzed like Gaia data. Additionally, the 2d and 3d profiles for each snapshot can be calculated with `add_profiles.sh` and the 1D density profile can be calculated with `stellar_profiles.sh`. 

## Dark matter Analysis

### Orbit validation

All of the following discussion is carried out in `notebooks/sim_analy/analyze_orbit.jl`. 

This script will check if the orbit is approximantly the same as the point particle. In detail, they should differ a little as escaping particles carry away angular momentum from the particle causing the orbit to shrink. We also calculate the time of each pericentre and apocentre at this point.

The script also determines the best fit final time of the snapshot by calculating the $\chi^2$ of the position, proper motions, and distance/radial velocity of the snapshot in the ICRS frame. The properties are saved in the analysis directory as `orbital_properties.toml` and the projected orbit on the sky as `skyorbit.fits`.

### Mass & Velocity profiles.



## Stellar analysis





![image-20241109104358175](/Users/daniel/Library/Application Support/typora-user-images/image-20241109104358175.png)
