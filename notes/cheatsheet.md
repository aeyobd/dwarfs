# Observational properties



# Isolation Runs



# MonteCarlo Orbits

- [ ] Setup directory `simulation/galaxyname/mc_orbits/simname/`
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

# Orbital runs

Gadget 4 will provide snapshots for each timestep, but now we have to compute useful observational and theoretical summaries of the galaxy's evolution.



## checklist

(more specifics below, this is for file management & sanity)

### Simulation

- [ ] Scale isolation to target halo size using rescale_nfw.jl (once per halo directory)
  - [ ] Check `iso_paths.sh` and `rerescale.sh`, running the latter
- [ ] Create simulation directory. Ensure only following files
  - [ ] agama_potential.ini (and any necessary dependencies, possibly make_pot.sh)
  - [ ] param.txt
    - [ ] Check consistency with isolation run
    - [ ] Check **gravitational softening**
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
- [ ] make readonly
- [ ] Create corresponding `analysis/galaxyname/haloname/orbitname/` directory with analysis scripts (`analyze.sh` typically)
- [ ] Run analysis scripts
  - [ ] Combine outputs, calculate centres
  - [ ] Calculate DM profiles
- [ ] Analyze orbit (`sim_analy/analyze_orbit.jl`) and check degree of deviation from expected orbit (in simulation directory, `orbit.csv`)
- [ ] Analyze DM particles with `analyze.jl` notebook.

### Analysis (Stars)

(isolation)

- [ ] Create `analysis/galaxyname/haloname/stars/` 
- [ ] Add template scripts (`calc_dist`, `paint`, `project`, `add profiles`)
- [ ] Calc distribution function & energies (once per halo)
- [ ] Paint stars
  - [ ] Check with `paint_stars_plots.jl` 
- [ ] Add  profiles
  - [ ] Check consistency with `iso_stars_plots.jl`?



(simulation)

- [ ] Create `analysis/galaxyname/haloname/orbitname/stars` directory
- [ ] Add template scripts (`project`, `add_profiles`, `stellar_profiles`)
- [ ] Project stars
- [ ] Add final / initial stellar profiles
- [ ] Add snapshot stellar profiles
- [ ]  `analyze_stars.jl`
- [ ] `projected_stars.jl` 
- [ ]  `tidal_tails.jl` 