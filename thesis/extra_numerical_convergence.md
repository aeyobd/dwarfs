# Numerical convergence and parameters {#sec:extra_convergence}

Here, we describe the numerical specifics of our simulation and analysis. We present convergence tests supporting our choices of softening and additional parameters. We find that alternate choices of numerical parameters neither improve convergence nor the resulting evolution. 

Due to the finite resolution of N-body simulations, collisional (close) encounters between particles are inevitable. Such collisional encounters are both unphysical (with possible arbitrary acceleration, violating the "collisionless" assumption) and computationally expensive (requiring many small time-steps). As a remedy, many N-body codes use a "softened" gravitational force law, where the force weakens when closer to a particle than a "softening length". The choice of softening length ideally balances resolution and computational speed.

Empirically, @power+2003 suggest that the ideal softening is
$$
h_{\rm grav} = 4 \frac{R_{200}}{\sqrt{N_{200}}},
$$ {#eq:power_softening}

where $h_{\rm grav}$ is the softening length, and $N_{200}$ is the number of particles within $R_{200}$. This choice balances integration time and only compromises resolution in collisional regime. For our isolation halo ($\rmax=6\,\kpc$, $\vmax=31\,\kms$) with $10^7$ particles, this works out to be $0.044\,{\rm kpc}$. Next we consider the effect of softening  as applied to our isolation halo.

The top panel of @fig:methods_convergence illustrates the influence of softening length for simulations of Scl's \smallperi{} orbit with $10^5$ particles. We show a model with the fiducial softening length ($h=0.14\,\kpc$), and larger and smaller softenings by a factor $\sqrt{10}$. The larger softening length is consistent with @eq:power_softening's prediction. With a larger softening length, the halo diverges from the expectation (using $100\times$ more particles) across most radii. On the other hand, the fiducial and smaller softening lengths predict similar final results. As computation time increases with decreasing softening length, the fiducial softening length balances efficiency and accuracy for this simulation. 



The two other key numerical parameters in \gadget{} are the tree-force and time-step accuracy parameters.  Specifically, \gadget{} opens a node if $M\,l/r^3 < \alpha |a|$, the node's mass is $M$, distance from the particle is $r$, side-length is $l$, and the particle's total acceleration is $a$. We adopt $\alpha =0.005$. We also elect to use adaptive time stepping with integration accuracy set to $\eta=0.01$. (particles must take time-steps smaller than $dt < \sqrt{2\,\eta\,h_{\rm grav} / a}$ for acceleration $a$).  The lower panel of @fig:methods_convergence tests changes to these parameters. We include a model with smaller integration accuracy ($\eta=0.003$ in "small timestep") and a stricter tree-force tolerance ($\alpha=0.001$, "high acc. force"). More precise tolerances on these parameters does not affect the evolution within the uncertainties of the final profile. 

As demonstrated in this subsection, stricter numerical accuracy and modified softening lengths do not affect our results. We thus conclude that our simulations are numerically well-converged (up to the "convergence" radius). 

![Isolation method convergence](figures/orbit_converg_methods.png){#fig:methods_convergence}

Figure: A comparison of the final profiles using different simulation methods for Sculptor's \smallperi{} model. The benchmark model is our fiducial $10^7$ particle run, and all other models use $10^5$ paritcles, with their "converged radius" marked by the black arrow and softening by the vertical bar(s). **Top:**  Models with $\sqrt{10}$ larger and smaller softening lengths than the fiducial. **Bottom:** Models with more precise timestep accuracy and gravitational force accuracy. **Summary:** Except for the model with a larger softening length, all simulations agree within uncertainties. 
