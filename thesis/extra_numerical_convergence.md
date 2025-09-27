# Numerical convergence and parameters {#sec:extra_convergence}

Here, we describe some convergence tests to ensure our methods and results are minimally impacted by numerical limitations and assumptions. See @power+2003 for a detailed discussion of various assumptions and parameters used in N-body simulations. 

## Softening

One challenge of N-body integration is close, collisional encounters, assuming point particles, cause divergences in the local force. However, this should not occur in a *collisionless* simulation. As a result, most collisionless N-body codes adopt a gravitational softening, a length scale below which the force of gravity begins to decrease between point particles.

@power+2003 empirically suggest that the ideal softening is
$$
h_{\rm grav} = 4 \frac{R_{200}}{\sqrt{N_{200}}},
$$

where $h_{\rm grav}$ is the softening length, and $N_{200}$ is the number of particles within $R_{200}$. This choice balances integration time and only compromises resolution in collisional regime.

For our isolation halo ($M_s=2.7$, $r_s=2.76$) and with $10^7$ particles, this works out to be $0.044\,{\rm kpc}$.We adpoted the slightly smaller softening which was reduced by a factor of $\sqrt{10}$ which appears to improve agreement slightly in the innermost regions. 

![Softening convergence](figures/iso_converg_softening.png){#fig:softening_convergence}

## Time stepping and force accuracy

In general, we use adaptive timestepping and relative opening criteria for gravitational force computations. To verify that these choices and associated accuracy parameters minimally impact convergence or speed, we show a few more isolation runs (using only 1e5 particles)

- constant timestep (...), approximantly half of minimum timestep with adaptive timestepping
- geometric opening, with $\theta = 0.5$.
- strict integration accuracy, (facc = ....)

![Isolation method convergence](figures/iso_converg_methods.png){#fig:methods_convergence}
