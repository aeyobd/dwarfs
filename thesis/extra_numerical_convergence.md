# Numerical convergence and parameters {#sec:extra_convergence}

Here, we describe the numerical specifics of our simulation and analysis. We present convergence tests supporting our choices of softening and additional parameters. We find alternate choices of numerical parameters are unlikely to improve our convergence or result in significantly different evolution.

## Softening

One challenge of N-body integration is close, collisional encounters, assuming point particles, cause divergences in the local force. However, this should not occur in a *collisionless* simulation. As a result, most collisionless N-body codes adopt a gravitational softening, a length scale below which the force of gravity begins to decrease between point particles.

@power+2003 empirically suggest that the ideal softening is
$$
h_{\rm grav} = 4 \frac{R_{200}}{\sqrt{N_{200}}},
$$

where $h_{\rm grav}$ is the softening length, and $N_{200}$ is the number of particles within $R_{200}$. This choice balances integration time and only compromises resolution in collisional regime. For our isolation halo ($\rmax=6\,\kpc$, $\vmax=31\,\kms$) with $10^7$ particles, this works out to be $0.044\,{\rm kpc}$, about $\sqrt{10}$ times larger than our adopted softening (@eq:softening_length).



![Softening convergence](figures/iso_converg_softening.png){#fig:softening_convergence}

Figure: Similar to @fig:numerical_convergence except for simulations with different softening lengths, $h_{\rm grav}$. $0.044\,\kpc$ is our fiducial softening length for this halo. All simulations here were ran in isolation of 5 Gyr with $10^{6}$ particles. 

@fig:softening_convergence illustrates the influence of softening length on isolation evolution. The simulation with  $h=0.14\,\kpc$ deviates the most from the expected profile. Instead, the final profiles for $h=0.044\,\kpc$ and $h=0.014\,\kpc$ are nearly identical. Smaller values of softening than $h=0.044\,\kpc$ likely will not affect the evolution in the inner regions. As computation time increases moderately with decreasing softening length, $h=0.044\,\kpc$ is an ideal choice balancing efficiency and accuracy. Our simulations therefore adopt this choice (as in @eq:softening_length).

## Time-stepping and force accuracy

In general, we use adaptive time-stepping and relative opening criteria for gravitational force computations. To verify that these choices and associated accuracy parameters minimally impact convergence or speed, we show a few more isolation runs (using only $10^5$ particles).

 We use in \gadget{} the relative tree opening criterion with the accuracy parameter set to $\alpha =0.005$ (so nodes are opened when $M\,l/r^3 < \alpha |a|$ for a distance from the node of mass $M$ of $r$, $l$ is the side length of the node, and $a$ is the particle's total acceleration), and adaptive time stepping with integration accuracy set to $\eta=0.01$ (particles must take time-steps smaller than $dt < \sqrt{2\,\eta\,h_{\rm grav} / a}$ for acceleration $a$). 

@fig:methods_convergence tests various additional simulation parameters in isolation. We compare against an equivalent model ran instead with `GADGET-2`, a model ran with constant timestep of 0.5Myr (`dt0.1` named after our dimensionless units), a simulation with the time-step accuracy set to $\eta=0.003$ (high accuracy), and a simulation using a geometric opening criterian instead of the force (geometric). In all cases, the final isolation velocity profile is nearly the same after 5Gyr. The numerical details likely do not limit our convergence (as marked by the arrow).



![Isolation method convergence](figures/iso_converg_methods.png){#fig:methods_convergence}

Figure: A comparison against different integration methods. See text for details on each model.
