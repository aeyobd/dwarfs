*Note on templates*. I typically include `template` folders in all of the relevant directories. Ideally, these should have only a few files and a few lines of code total. When adding a new model or galaxy or orbit in any category, I recommend using the template as a starting point. These keep the parameters and files in the same format which makes analysis easier and comparing differences between models much easier as well.

In general, there should be a little `clean.sh` script which removes all output and non-essential data files in case the analysis needs to be repeated. 

This document, `workflow.md`, `Gadget` documentation, and documentation for any julia scripts are the goto references on the meaning of each file. To facilitate flexibility and concise template files, documentation in each file is minimal to non-existent.

Additionally, I keep 1 simulation per directory with all necessary initial conditions present in that directory. I also keep the analysis in the corresponding `analysis/...` directory (with same path) which also symbolically links to the simulation directory in `simulation`. The reason for the separation of simulations and analysis on the top level of the file tree is so that the simulation directories can be kept readonly for clarity and so simulation directories can be easily and efficiently rsync'ed or stored on a different filesystem.

**TODO**: Carefully construct final templates for all models & save logs (mosty done) 

**TODO**: All the analysis scripts should ideally delete the target file, but this is not very well implemented now.

**TODO**: Document each model with a short description (`README.md` file)

**TODO**: Version control of all simulation/analysis scripts.

# Observational properties

## Checklist

- [ ] Literature review galaxy
  - [ ] Find measurements of radial velocities, structural properties, density profiles, gas, etc.
  - [ ] Find simulations / dynamics papers
  - [ ] Look at chemistry / SFH for interesting details. 
- [ ] Data collection
  - [ ] Assemble sample of radial velocities
  - [ ] Assemble photometric / stellar membership catalogue (can fit either kinematic members or background + density profile)
- [ ] Membership selection
  - [ ] J+24



## J+24 analysis checklist

- [ ] Check how many PSAT members (should be >~ 10)
- [ ] Check correlation of L with space, pM, and CMD components for near-membership likelihoods
- [ ] Calculate density profiles
  - [ ] Changes to magnitudes?
  - [ ] Changes to PSAT?
  - [ ] Changes to centring/structural properties?
  - [ ] Changes to selection method?
  - [ ] BG selection?

# MonteCarlo Orbits

## Checklist

- [x] Setup directory `simulation/galaxyname/mc_orbits/simname/`
  - [ ] Set initial conditions (location and magnitude of uncertanties)
  - [ ] Setup potential `agama_potential.ini`
  - [ ] Double check gadget scripts (`param.txt`, `run.sh`)
- [ ] Submit run
- [ ] Create analysis directory `analysis/galaxyname/mc_orbits/simname/`
  - [ ] Add scripts (.....)
  - [ ] Run `analyze.sh`
- [ ] Analyze pericentre/apocentre distributions  & correlations with `analysis/galaxyname/mc_orbits/analyze.jl` 
- [ ] Use notebook to chose initial conditions for several interesting orbits (smaller pericentre, mean, etc)

Special cases

- [ ] Create `simulation/galaxyname/mc_orbits/simname_special_cases/`
  - [ ] Initial conditions `initial_conditions.toml` describes each orbit
  - [ ] 

A long term project is to transition to an in-house orbit integrator to allow for multi-component moving potentials and dynamical friction. Stay tuned :).

# Isolation Runs

## **checklist**

- [ ] Create model directory based on `template`. 
  - [ ] param.txt
    - [ ] Check **gravitational softening** is (scaled) from isolation run
    - [ ] Check simulation time 
    - [ ] Check output frequency, memory, and timelimit settings
    - [ ] Check accuracy, tree options...
  - [ ] make_init.sh
    - [ ] Set **halo** initial conditions path
  - [ ] halo.toml (rescaled to this value to be more physical.) 
  - [ ] run.sh (should contain ~1 bash line to run the simulation)
  - [ ] clean.sh (optional)
- [ ] Create initial conditions `make_init.sh`
- [ ] (rsync into submit directory) and submit with `submit_gadget.py`

## Convergence tests

Fortunately, there are only a few free parameters in the calculation of an N-body code. 

- How much does particle number affect the results
- Gravitational softening?
- Gravitational tree force accuracy? Method?
- Time integration accuracy?

# Orbital runs

Gadget 4 will provide snapshots for each timestep, but now we have to compute useful observational and theoretical summaries of the galaxy's evolution.



## checklist

(more specifics below, this is for file management & sanity)

### Simulation

- [x] Scale isolation to target halo size using rescale_nfw.jl (once per halo directory)
  - [ ] Check `iso_paths.sh` and `rerescale.sh`, running the latter
- [x] Create simulation directory. Ensure only following files
  - [x] agama_potential.ini (and any necessary dependencies, possibly make_pot.sh)
  - [ ] param.txt
    - [ ] Check consistency with isolation run
    - [ ] Check **gravitational softening** is (scaled) from isolation run
    - [ ] Check simulation time 
    - [ ] Check output frequency, memory, and timelimit settings
  - [ ] make_init.sh
    - [ ] Set **orbit.csv** file path 
  - [ ] run.sh (should contain ~1 bash line to run the simulation)
  - [ ] clean.sh (optional)
- [ ] Check initial conditions with `notebooks/sim_analy/initial_conditions.jl`, ensuring both vmax, rmax are reasonable, v circ & density profile matches expected halo, and centre is appropriate for orbit.
- [ ] (rsync into submit directory) and submit with `submit_gadget.py`

### Analysis

- [x] Check logs
- [x] Rsync model into analysis directory
- [x] make readonly
- [x] Create corresponding `analysis/galaxyname/haloname/orbitname/` directory with analysis scripts (`analyze.sh` typically)
- [x] Run analysis scripts
  - [ ] Combine outputs, calculate centres
  - [ ] Calculate DM profiles
- [x] Analyze orbit (`sim_analy/analyze_orbit.jl`)
  - [ ]  check degree of deviation from point-a orbit (in simulation directory, `orbit.csv`)
  - [ ] Determine pericentres / apocentres
  - [ ] Find best-fit final position of simulation

- [x] Analyze DM particles with `analyze.jl` notebook.

### Analysis (Stars)

(isolation run)

- [ ] Create `analysis/galaxyname/haloname/stars/` 
- [ ] Add template scripts (`calc_dist`, `paint`, `project`, `add profiles`)
- [ ] Calc distribution function & energies (once per halo)
- [ ] Paint stars
  - [ ] Check with `paint_stars_plots.jl` 
- [ ] Add  profiles
  - [ ] Check consistency with `iso_stars_plots.jl`?

(orbit run)

- [x] Create `analysis/galaxyname/haloname/orbitname/stars` directory
- [x] Add template scripts (`project`, `all_profiles`, `stellar_profiles`)
- [x] Project stars
- [x] Add final / initial stellar profiles
- [ ] Add snapshot stellar profiles
- [ ]  `analyze_stars.jl`
- [ ] `projected_stars.jl` 
- [ ]  `tidal_tails.jl` 