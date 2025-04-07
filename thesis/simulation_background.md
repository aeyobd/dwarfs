# Introduction to Dark Matter simulations

In this section, we will cover

- How are dark matter simulations conducted
- Interpretations and uncertainties 
- Methods



![Break radius validation](figures/idealized_break_radius.png){#fig:idealized_break_radius}

Figure: The break radius of the simulations is set by 



From this argument, we note that the following properties must be approximately true for tides to occur:

- Close enough pericentre. The other break radius $r_J$ implies that if the host density is 3x the satellite, stars will be lost
- Corresponding time since last pericentre: If the time since last pericentre is not $\sim$ consistent with an observed break in the density profile, then tides 
- Halo evolution. As found in @EN2021, galaxies evolve along well defined tidal tracks (assuming spherical, isotropic, NFW halo, which may not be true, see ...). These tracks tend to "puff up" the stellar component while also removing dark matter mass, leaving a smaller, compacter DM halo with a more extended stellar component.
  - This information is mostly related to the statistical initial distribution of satellites from cosmology [ludlow+2016; @fattahi+2018]

# Potentials

From the @nfw1996 paper, eqns. 3 & 4
$$
\frac{\rho}{\rho_c} = \frac{\delta_c}{(r/r_s)(1+r/r_s)^2}
$$
where 
$$
\delta_c = \frac{200}{3}\frac{c^3}{[\ln(1+c)-c/(1+c)]}
$$
$c$ is concentration parameter, $\rho_c$ is the critical density of the universe, and $r_s$ is the characteristic scale length of the halo.

The NFW halo is sometimes described by $M_{200}$. $r_{200}$ is the radius at which the mean density of the halo interior to $r_{200}$ is 200 times the critical density of the universe, and $M_{200}$ is the mass contained inside $r_{200}$. As equations:
$$
\rho_{200} = 200\rho_{c} = \frac{M_{200}}{(4\pi/3) r_{200}^3}
$$

$$
r_{200} = \sqrt[3]{\frac{1}{200}\frac{3M_{200}}{4\pi \rho_{\rm c
}}}
$$


$$
M_{200} = \frac{4\pi}{3} r_{200}^3\ \rho_{200}
$$

$M_{200}$ is also sometimes called the virial mass of the halo. $r_{200}$ is directly related to $r_s$ by
$$
r_{200} = c\,r_s
$$


Another useful definition is
$$
A(x) \equiv \log (1+x) - \frac{x}{1+x}.
$$


We will also define a dimensionless radius
$$
x \equiv r/r_s.
$$

### 

A simple substitution to the definition gives


$$
\rho(x) =  \frac{\rho_s/3}{x\ (1+x)^2}
$$
where
$$
\rho_s = \frac{c^3}{A(c)} \rho_{200}
$$
and $A(c)$ is as above. The characteristic density can also be written in terms of scale mass,  $M_s = M_{200}/{A(c)}$  (see below), giving
$$
\rho_s = \frac{c^3}{A(c)} \frac{M_{200}}{(4\pi/3)\ r_{200}^3}  = \frac{3M_s}{4\pi\, {r_s}^3}
$$
Note that the NFW density profile is the same as an alpha-beta-gamma profile where $\alpha=\gamma=1$ and $\beta =3$.

### Circular velocity

The circular velocity in terms of $v_{200} = \sqrt{G M_{200} / R_{200}}$ is
$$
\left(v_{\rm circ}/v_{200}\right)^2 = \frac{A(x)/x}{A(c)/c},
$$

or in terms of $M_s$ and $r_s$, 
$$
v_{\rm circ}^2 = \frac{G M(r)}{r} = \frac{G M_s A(r/r_s)}{r}.
$$
Another parameterization of the NFW profile is in terms of the maximum circular velocity $v_{\rm circ}^{\rm max}$ and the radius at which it is reached $r_{\rm circ}^{\rm max}$. Given the scale radius, 
$$
r_{\rm circ}^{\rm max} = \alpha\ r_s
$$

where $\alpha\approx2.16258$, and $v_{\rm circ}^{\rm max}$ can be found from either of the equations above.

