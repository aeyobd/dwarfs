## Introduction

- *opening*

  - dwarf galaxies are extreme DM environments
  - outline of this section

- Observational history

  - *Opening*

    - discovery and uncertain initial classification of dwarfs

    - high mass-to-light ratios in dwarfs, speculation and dark matter

    - what is a dwarf galaxy?

    - the local group system of galaxies

    - types of (local) dwarfs 

    - why do we focus on classicals?

- Cosmological context

  - *opening*
    - what is LCDM?
    - dark matter properties and general predictions
  - *structure formation*
    - growth of perturbations, formation of halos
    - hierarchical structure and accretion
    - role of dark matter properties in dwarfs
  - *structure of halos*
    - NFW definition
    - virial mass of NFW
    - vmax and rmax parameterization
    - why are mass and concentration related?

  - *galaxy formation inside LCDM*
    - abundance matching analysis
    - SMHM relation, steep in the dwarf regieme
    - extended dark matter halo dominates mass, especially for dwarfs
    - complexity in SMHM and galaxy formation

- Observing faint features in dwarfs

  - *The Gaia space telescope*
    - why is parallax / pm hard?
    - what is *Gaia* generally
    - how does *Gaia* determine astrometry, what does *Gaia measure
  - *Gaia and the local group*
    - *Gaia*'s impact on MW science
    - *Gaia* role in dwarf galaxy studies: motions & membership
  - *challenges in dwarf galaxies*
    - results of updated facilities
    - small-scale challenges for LCMD
    - **to add:** open questions
  - *dwarf galaxy light profiles*
    - why do we care about light profiles?
    - the four types of light profiles
    - Why is exponential the default?
    - Limitations to a pure exponential
    - moskowitz+walker 2020
    - summary

  - *Sculptor and Ursa Minor light profiles*
    - sestito++ reports of deviations from exponential
      - shown as a figure
      - thesis goal

- Tidal signatures

  - *Opening*
    - challenges in cosmological simulation, role of idealized experiments
    - early work on tides
    - **in progress** recent advancements
    - related, galaxy specific work
  - *Tidal and break radii*
    - define and intuit jacobi radius
    - and break radius

  - *A simple tidal simulation*
    - description of toy model
      - illustration of tidal effects, and break radius.

- Outline of thesis



## Gaia membership selection

- *opening*
  - exponential profiles, unusual scl umi, need critical analysis of density profiles
- Stellar membership with Gaia
  - Need to assess background density, *Gaia* provides more information
  - J+24 uses three planes of information, finding distant members
  - Likelihood of 3 components
  - How are membership probabilities determined?
  - How do we select our fiducial sample?
- The effects of membership criteria
  - rell and tangent plane
  - description of progressive membership samples
  - visual analysis of membership selection
  - differences between galaxies
  - the spatial component of the likelihood
  - ___ does not matter!
- Density profiles
  - How do we calculate these
  - comparison of profiles for each sample
  - differences between samples and the spatial component
  - does Muñoz 1 matter?
  - Comparison of classical dwarfs, 

## Numerical Methods

- *opening statement*
  - methods are to analyze tidal effects on Scl and UMi
- Orbital estimation
  - MC sampling of orbits
  - *Galactocentric frame*
  - *Milky Way potential specification*
  - *Orbits of Scl*
  - *Orbits of UMi*
- Initial conditions
  - agama method and density truncation
  - *Initial halo selection*
    - Using cosmological relations to derive (initial) halos
    - role of velocity dispersion as a constraint
    - does the choice of halo matter?
- Numerical methods
  - introduction to n-body simulations
  - *isolation run and simulation parameters*
    - Isolation setup
    - softening selection
  - *Numerical fidelity*
    - convergence of profiles, limitation of converged radius
  - *Orbital evolution*
    - Initializing simulation via shifting and rescaling
  - *Tidal mass losses*
    - How do we follow the centre?
    - limitations on the centring error?
  - *Stellar component*
    - Initial stellar profile
    - The energy tagging method



