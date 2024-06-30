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

# Misc

Jaclyn's data 

https://www.canfar.net/storage/vault/list/jax/DSPHS_Gaia_eDR3/December_2022/1-Component_HB_Probabilities



## Numerical Setup

See @power+2003 for a discussion of this...
$$
h_{grav} = 4 \frac{R_{200}}{\sqrt{N}}
$$

Asya gave the similar relation
$$
h = 0.005R_{\rm max} (N/1e6)^{-1/2}
$$
which is $5/4$ of the value from @power+2003



## Units

| Quantity   | Units               |
| ---------- | ------------------- |
| G          | 1                   |
| Mass       | $10^{10}$ M$_\odot$ |
| radius     | 1 kpc               |
| velocities | 207.4 km/s          |
| time       | 4,715,000 yr        |



# COSMA

N-body simulations are computationally expensive. I have been using [COSMA](https://www.dur.ac.uk/icc/cosma/), but the Canadian supercomputers or even a local machine would also work.

Instructions to get an account: https://www.dur.ac.uk/resources/icc/cosma/gettingAnAccount.pdf. 

I typically ssh into cosma to edit any files or run scripts.

### Jupyter

See https://www.dur.ac.uk/icc/cosma/support/jupyterhub/. 

to launch a jupyter session, 

``` bash
ssh -i .ssh/id_ed25519_cosma -N -L 8443:login8b.cosma.dur.ac.uk:443 dc-boye1@login8b.cosma.dur.ac.uk
```

then open https://localhost:8443

### Modules/Dependencies

- `fftw` for fourier transforms
- `openmpi` required for gadget
- `gnu_comp` compiler
- `hdf5` for files
- `python`, using latest python version available
- `julia`: some of my newer scripts are written in julia 

### SLURM

https://www.dur.ac.uk/icc/cosma/support/slurm/

# Zeno

Zeno is a collection of scripts with documentation on [ifa](https://home.ifa.hawaii.edu/users/barnes/software.html). In particular, we want the following programs

- `halogsp` generates a 1d density profile given the halo parameters. First argument is the outfile. Then the other parameters specify halo shape
- `gspmodel` creates the actual model

### Installation

```
wget https://github.com/joshuabarnes/zeno/archive/refs/tags/v0.10.0-beta.tar.gz
```



add to bashrc

```
	export ZENOPATH="/Users/<yourname>/zeno"
	export ZCC="gcc"
	export ZCCFLAGS="-std=gnu99 -DLINUX -I$ZENOPATH/inc"
	export ZLDFLAGS="-L$ZENOPATH/lib"
	export ZENO_SAFE_SELECT="true"
	export ZENO_MSG_OPTION="all"
```

can comment out graphics part of makefile so does not depend on GLUT / openGL (much easier to build in some environments)

and run `make -f Zeno >& zenomake.log`.



### Converting to HDF5

Zeno produces a structured binary file, not useful for gadget.

There are two possible methods here: using Rapha's modified zeno code (with a utility) or my silly little python script.

- `tsf` writes out the details of a halo file as text, which can then be read in with python (a bit roundabout but works)
- `snapgadget` to convert to gadget format. 



in `/cosma/home/durham/dc-boru1/zeno/bin` using ICFormat=1.

https://github.com/kyleaoman/zeno

then run with gadget for short amount of time to convert to convert output

### Rescaling Halo

use `rescale.py`, put in mass/

uses some other functions/etc.



### Scaling

Zeno creates a profile with a scale radius $R_s=1$ and mass inside this radius $M_s=1$. The NFW profile is given by
$$
\rho(r)/\rho_{\rm crit} = \frac{\delta_c}{(r/r_s)(1+r/r_s)^2}
$$


## Orbits

## To Find Initial Orbit

use leapfrog to integrate backwards

give RA, DEC, distance, v_disp, proper motions, RV from GaiaDR2 and DR3, 

standard units mas/yr etc.?

Uses MC to run orbits observationally consistant, calculate peri/apocenters, and generate histogram. Gives median, 16th and 84th percentile pericenters. 

outputs in a text file

`setInOrbit.py`





# Agama



# Gadget4

documentation: https://wwwmpa.mpa-garching.mpg.de/gadget4/

## Building

The source code can be cloned from https://gitlab.mpcdf.mpg.de/vrs/gadget4. I am working from my own fork of the repository. To make, run the make command in the main folder of the repository, e.g.

```bash
make CONFIG=~/path/to/Config.sh DIR=~/path/to/out/
```

in my case

```
make CONFIG=~/dwarfs/gadget/Config.sh DIR=~/dwarfs/gadget
```

where Config.sh is the configfile

### Agama patch

run 

```
patch -r -u -N -d /path/to/gadget4/ -p 1 < example_nbody_simulation_gadget4.patch
```

from the `Agama/py` directory. This will apply the Agama 

### The Configfile

Gadget contains many compile-time settings. For this run we don't need many

I do add in a `EXTERNALGRAVITY_MW` option to the makefile which allows the specification of our four component milky way potential (discussed elsewhere)

### Agama potential

hard coded as "agama_potential.ini". 

## Running

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



- `InitCondFile` put name of hdf5 file (without extension) of initial particle positions/velocities



## Isolation Run

/cosma/home/durham/dc-boru1/Gadget2_Noah/Gadget2-ABLL

First, rescale profile with

```bash
python rescale.py snapshot_000.hdf5 iso_scaled.hdf5 2.11 2e8
```

for e.g. fornax with 2e8 Msun within 2.11 kpc scale radius



After running, need to centre the profile via centreSnapshot. 

change parameters 3, 4, 92, 99, 110, 111

ioptions

- UNEQUALSOFTENINGS
  - unequal softening lengths
- PEANOHILBERT
  - performance
- WALLCLOCK
- DOUBLEPRECISION
- SYNCHRONIZATION
  - hierarchy timesteps
- HAVE_HDF5
- H5_USE_16_PAI
- OUTPUTPOTENTIAL
- OUTPUTACCELLARATION

## For orbit run

we specify time in anchors.txt file, change third file in 

/cosma/home/durham/dc-boru1/Gadget-RAPHAMW/source/Gadget2



Use `setInOrbit.py` to put galaxy in galactocentric coordinates



Change lines 3, 4, 67, 129



## Units

| Quantity   | Units               |
| ---------- | ------------------- |
| G          | 1                   |
| Mass       | $10^{10}$ M$_\odot$ |
| radius     | 1 kpc               |
| velocities | 207.4 km/s          |
| time       | 4,715,000 yr        |



# Postprocessing

Gadget 4 will provide snapshots for each timestep, but now we have to compute useful observational and theoretical summaries of the galaxy's evolution.

## Painting stars

As stars are sufficiently spread out in a dwarf galaxy, we can assume that they are simply attached to a dark matter particle, with probabilities governed by the ratio of the distribution functions. 






# Analyzing observations

Observationally, we know the 6D kinimatic and have estimated total stellar mass 

## Orbits



## Initial profile parameters

to get the stellar mass, we need the absolute magnitude and stellar mass to light ratio. For example, Sculptor has an absolute magnitude of $M_V = -11.1$ [@mcconnachie2012], and the sun has an absolute magnitude of 4.83, so with a stellar mass to light ratio of 1.6, this gives a total stellar mass of $M_\star=3.8\times10^6\,$M$_\odot$ .

Next, we use that
$$
M_\star = m_0\,\nu^\alpha \exp({-\nu^\gamma})
$$
where $m_0 = 3\times10^8 M_\odot$, $\alpha=3.36$, and $\gamma=-2.4$. $\nu = V_{\rm max}/$50km/s, from an emperical fit to @fattahi2018. For the example of Sculptor, $\nu=0.643$ solves the above equation, so $V_{\rm max}=32.2$km/s. 

Finally, the characteristic radius is from @ludlow2016 (see their appendix C). The full equations are complex, but briefly, from a given M and z, you can calculate their fit to the expected concentration. We then need to solve for a value of $M_s$ which gives a value of $c$ which together predict the correct $V_{\rm circ, max}$. 

Using that the maximum circular velocity occurs at $\alpha r_s$ where $\alpha = 2.163$ (is a numerical solution), and the equation for circular velocity (way back above but repeated here:)
$$
\left(V_c/V_{200}\right)^2 = \frac{c}{r/r_s}\frac{A(r/r_s)}{A(c)}
$$
where $V_{200} = \sqrt{G\ M_{200} / r_{200}}$. The solution here is $M_{200} = 5.4e9$ where $c=13.0$, and $r_s = 2.86$ and $M_s = 6.12e8$.  





# Development

Here, I discuss my reorganizations of the projects

## snapshot

This class is simply a method to open and close hdf5 gadget files more easily. all new work should use this class as it makes things more consistent and easier

## coords

use these methods to perform coordinate transformations







# Tests

## Conservation laws

- Momentum: $\sum m\,v_i={\rm constant}$ 
- Angular momentum: $\sum m\,r\times v = {\rm constant}$ 
- Energy
  - $T = \sum \frac{1}{2} m v^2$
  - $ U = \frac{1}{2} \sum_i \sum_{j\neq i} \frac{G m_i m_j}{|x_i - x_j|} = \frac{1}{2} \sum_i m_i \Phi(x_i)$. The factor of 1/2 is so we only count each pair once.
  - $U_{\rm ext} = \sum_i m_i \Phi_{\rm ext}(x_i) $
  - $T + U + U_{\rm ext} = {\rm constant}$

## Two body

Here, I just test a keplarian 2-body system. The solution is the same for the one-body system, so we expect Kepler's laws to hold--the orbit should be elliptical with a period
$$
P^2 = \frac{4\pi^2}{\mu} a^3
$$
where $\mu = G(m+m)$ is the reduced gravitational mass of the system, and $a$ is twice (?) the semi-major axis of one of the orbits.



## Three body

To test the three body problem, I just use the special case (...)



## Hernquist  and NFW potential

For a spherical potential $\Phi$, we know that the roots of 
$$
0=\frac{1}{r^2} + 2\frac{\Phi(r) - E}{L^2}
$$
provide us with the maximum and minimum orbital radius, $r_{\rm min}$ and $r_{\rm max}$. We can also calculate the period in $r$ with 
$$
T_{r} = \int_{r_{\rm min}}^{r_{\rm max}} \frac{2}{\sqrt{2|E-\Phi(r)| - \frac{L^2}{r^2}}}\ dr
$$
I verify these relations for the hernquist and NFW spherical potentials. 

## Disk potential

Disk potentials are more complex but we can still check the overal relationships .

