## Background & Past Work

- What is dark matter? Why do we look at dwarfs?
- Forms of dark matter, lambda-CDM, and dwarf galaxies
- How does gravity affect dwarfs, theory of tidal perturbations
  - EN21, peñarrubia+09, etc.
- Instances of dwarfs undergoing weird processess
- Alternative processes and uncertainties in the evolution of dwarfs



## Gaia Membership Selection

![Sculptor selection criteria. ](figures/scl_selection.png){#fig:sculptor_selection}

We use J+24 data, updated PM from the MV 2020a catalogue.

J+24 select members using a multi-component Baysian algorithm. 

- First, stars with poor Gaia astrometry, color excess nearby parallaxes, and proper motions > 10 mas/yr are removed
- Stars are assigned a likelihood based on the location on the CMD (assuming)



Simple cuts, such as selecting PM, CMD based on purely geometric cuts reveal a similar density profile. 



One caveat is the original algorithm takes spatial position into account. When deriving a density profile, this assumption may influence the derived density profile, especially when the galaxy density is fainter than the background of similar appearing stars. To remidy this and estimate where the background begins to take over, we also explore a cut based on the likelihood ratio of only the CMD and PM components. 

## Searches for tidal tails

![image-20250313132102118](/Users/daniel/Library/Application Support/typora-user-images/image-20250313132102118.png)

## Density Profile Reliability and Uncertainties

Given a set of probable members, the density profile is given by EQUATION, where we use  ... 



- Using J+24 data, we validate
  - Check that PSAT, magnitude, no-space do not affect density profile shape too significantly
- Our "high quality" members all have > 50 member stars and do not depend too highly on the spatial component, mostly corresponding to the classical dwarfs
- We fit Sérsic profiles to each galaxy
  - The Sérsic index, $n$, is a measure of the deviation from an exponential. Exponentials have $n=1$, whereas more extended dwarf galaxies will have higher $n$
- To better estimate the uncertanties due to unknown galaxy properties and flexibility in the likelihood cut, we can 



![Density profiles](figures/density_profiles_medley.png){#fig:sculptor_observed_profiles}

- Above figure: add simple cuts & delve to convince more...

Dark Matter is one of the biggest open questions in modern astronomy.
Hello @fig:sculptor.





## Comparison of the Classical dwarfs

![image-20250313130050775](/Users/daniel/Library/Application Support/typora-user-images/image-20250313130050775.png)

![image-20250313130110550](/Users/daniel/Library/Application Support/typora-user-images/image-20250313130110550.png)

![image-20250313130043114](/Users/daniel/Library/Application Support/typora-user-images/image-20250313130043114.png)



## Summary

- Of the classical dwarfs, UMi & Scl stand out statistically, with high $n$ given their luminosity
- Including fainter dwarf galaxies, Boo 3 and Boo 1 appear to also have extended density profiles
  - Deeper data would be required to robustly measure this
