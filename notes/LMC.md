

| parameter                | value                                 | Source              |
| ------------------------ | ------------------------------------- | ------------------- |
| $\alpha$                 | 80.8                                  | McChonnachie2012    |
| $\delta$                 | -69.8                                 | McChonnachie2012    |
| distance                 | $49.59 \pm 0.54$ kpc                  | Pietrzyński+2019    |
| $\mu_\alpha \cos \delta$ | $1.910 ± 0.020$ mas yr$^{-1}$         | Kallivayalil+2013   |
| $\mu_\delta$             | $0.229 ± 0.047$                       | Kallivayalil+2013   |
| RV                       | $262\pm3$                             | in Mcconnachie 2012 |
| mass                     | $10-20 \times 10^{10}\,{\rm M}_\odot$ | Vasiliev2023        |

## Dynamical Friction



### Review of point-integration methods



The Chandrasekhar formula for an extended host
$$
\frac{d{\bf v}}{dt} = -\frac{4\pi\,G^2\,M\,\rho\,\ln\Lambda}{v^2} \left({\rm erf}(X) - \frac{2X}{\sqrt\pi} \exp(-X^2)\right) \frac{{\bf v}}{v}
$$
(from Binney & Tremaine 2008)

where $M$ is the satellite mass, $\bf v$ is the velocity of the satellite (relative to the host), $\rho$ is the host density at the satellite position, and $X = v/\sqrt2 \sigma$ where $\sigma$ is the local velocity dispersion of the host. The Coulomb factor is $\Lambda = r/\epsilon$ where $r$ is the distance between the host-satellite, and $\epsilon$ (from Jethwa et al. (2016)) is given by
$$
\epsilon = \begin{cases}
2.2 r_s - 14 kpc & r_s > 8kpc \\
0.45 r_s & r_s < 8kpc
\end{cases}
$$
(as described in Cautun+2019).

To calculate the velocity dispersion, we can use the (integrated) spherical Jeans equation assuming isotropy:
$$
\sigma_v^2(r) = \frac{1}{\rho} \int_r^\infty \rho(r') \frac{G\,M(r')}{r'^2} \, dr'
$$
In practice, we calculate all the values of $\sigma^2$ in a single integral and then interpolate over these to calculate the current dynamical friction of the system. V+21 actually simply set this value to zero.

## Non-intertial Acceleration

Since the LMC is 10-20% the mass of the Milky Way, the binary interaction between the pair should influence the resulting evolution. In this case, the acceleration of the Milky Way due to the LMC is 
$$
{\bf a}_{\rm MW} = -\nabla\Phi_{\rm LMC} ({\bf 0})
$$
This acceleration can then be subtracted from the current frame of reference uniformly so that the MW remains at the origin. I have also tried the other naive method of using Newton's third law, but this produces a worse match to the V+21 N-body simulation.



