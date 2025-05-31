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

Jaclyn's data 

https://www.canfar.net/storage/vault/list/jax/DSPHS_Gaia_eDR3/December_2022/1-Component_HB_Probabilities

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

Once zeno and my library is installed, the procedure to generate a halo is 

- calculate the profile by running `./make_profile.sh nfw`. This will write the zeno-file nfw.gsp to `profiles` and also outputs a .csv file with the same name for normal people to read in
- Create an nbody model by using `./make_halo.sh 1e4 profiles/nfw.gsp nfw_1e4` which in this example will generate a profile with 1e4 particles written to `halos/nfw_1e4`. 
- Use the NBody model for sciences!

### HaloGSP

Reverse engineered documentation.

HaloGSP generates the density profile ($\rho$, $M$, for each $r$). The parameters are:

- `m_a=1`, the mass within the scale radius `a`. Note that this is not $M_s$
- `a=1` is the scale radius
- `b=4` is the radius the taper begins
- `taper=exp` the taper formula (one of `exp`, `gauss`, or `sw`)
- `npoint=257` the number of points to calculate the density ate
- `rrange=1/4096:16` is the range the radii are calculated. The radii are evenly spaced in log r

All taper modes are defined as piecewise densities which are continuous at $\rho_b = \rho_{\rm NFW}(b)$; i.e.
$$
\rho(r) = \begin{cases}
\rho_{\rm NFW}(r) & r\le b \\
\rho_{\rm taper}(r) & r > b
\end{cases}
$$


Internally, zeno uses $\rho_{\rm NFW}(r) = m_a / (4\pi\,A(1)\,r\,(a + r)^2)$ which is as expected. 

The exponential taper is
$$
\rho_{e}(r) = \rho_b\,(b/r)^2\,\exp(-2\gamma\,(r/b-1))
$$
where $\gamma = b/(b+a) - 1/2$.

The gaussian taper is
$$
\rho_g(r) = \rho_b (b/r)^2 \exp(-\gamma ((r/b)^2 - 1))
$$
and finally, the `sw` taper (not sure the abbreviation here, but intended to be fast) is
$$
\rho_{\rm sw}(r)/\rho_b = (r/b)^{\gamma_s}\,\exp(-(r-b)/a)
$$
where $\gamma \ne \gamma_s  = b/a - (a+3b)/(a+b)$.



### GSPModel & GSPRealize

Zeno includes 2 (3?) programs to sample particles from the profile. `gspmodel` uses rejection sampling in velocities which does not actually create an equilibrium model. `gsprealize` instead uses the distribution functions to correctly sample velocities. There is also a `gsprealize_2.0` which I am entirely unsure of the differences.

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

There are two possible methods here: using Rapha's modified zeno code (with a utility) or I have a julia script

- `tsf` writes out the details of a halo file as text, which can then be read in with python (a bit roundabout but works)
- `snapgadget` to convert to gadget format. 

in `/cosma/home/durham/dc-boru1/zeno/bin` using ICFormat=1.

https://github.com/kyleaoman/zeno

then run with gadget for short amount of time to convert to convert output

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

### The Configfile

Gadget contains many compile-time settings. For this run we don't need many

I do add in a `EXTERNALGRAVITY_MW` option to the makefile which allows the specification of our four component milky way potential (discussed elsewhere)

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

## Numerical Setup

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
which is $5/4$ of the value from @power+2003



## Isolation Run

First, rescale profile with

```bash
python rescale.py snapshot_000.hdf5 iso_scaled.hdf5 2.11 2e8
```

for e.g. fornax with 2e8 Msun within 2.11 kpc scale radius

# Gadget tests

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

Disk potentials are more complex but we can still check the overall relationships & the behaviour of epicycles in the potential.

# Postprocessing

Gadget 4 will provide snapshots for each timestep, but now we have to compute useful observational and theoretical summaries of the galaxy's evolution.

## Painting stars

As stars are sufficiently spread out in a dwarf galaxy, we can assume that they are simply attached to a dark matter particle, with probabilities governed by the ratio of the distribution functions. 

To find the distribution function, we use Eddington inversion (eq. 4-140b in BT87)
$$
f({\cal E}) = \frac{1}{\sqrt{8}\, \pi^2}\left( \int_0^{\cal E} \frac{d^2\rho}{d\Psi^2} \frac{1}{\sqrt{{\cal E} - \Psi}}\ d\Psi + \frac{1}{\sqrt{\cal E}} \left(\frac{d\rho}{d\Psi}\right)_{\Psi=0} \right)
$$


In practice the right term is zero as $\Psi \to 0$ as $r\to\infty$, and if $\rho \propto r^{-n}$ at large $r$ and $\Psi \sim r^{-1}$ then $d\rho / d\Psi \sim r^{-n+1}$ which goes to zero provided that $n > 1$. 
