## NFW (Halo)

### Source definition

From the @nfw1996 paper, eqns. 3 & 4
$$
\frac{\rho}{\rho_c} = \frac{\delta_c}{(r/r_s)(1+r/r_s)^2}
$$
where 
$$
\delta_c = \frac{200}{3}\frac{c^3}{[\ln(1+c)-c/(1+c)]}
$$
$c$ is concentration parameter, $\rho_c$ is the critical density of the universe, and $r_s$ is the characteristic scale length of the halo.



### Some useful definitions

The critical density of the universe is 
$$
\rho_{c} = 3H^2/8\pi G \approx127.35M_\odot/{\rm kpc}^3.
$$


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
### Density

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

### Contained mass

By integrating $\rho$, the contained mass interior to $r$ is
$$
M(x) = M_s\ A(x)
$$
where
$$
M_s \equiv (4\pi/3) \rho_s r_s^3 = \frac{M_{200}}{A(c)}
$$
 so $M(x) = M_s  A(x)$. Note that this differs from the mass contained within $r_s$ by $A(1) \approx 0.193147$.  

### Potential

From Poisson's equation, 
$$
\Phi(x) = -\frac{G M_s}{r_s} \ \frac{1}{x} \ \ln\left(1+x\right)
$$

 (The constant of integration are chosen such that $\Phi(0) = \Phi_0$ and $\Phi$ vanishes at $\infty$.)

### Acceleration 

By differentiating the potential,
$$
\vec{a}(x) = -\frac{G M_s}{r_s^2} \frac{A(x)}{x^2} \hat{r}
$$
Also is consistent with using contained mass ($F=G M(r)/r^2$)

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



## Truncated NFW



## Hernquist (Bulge)

The @hernquist1990 density profile for the galactic bulge is parameterized in terms of a characteristic mass, $M$, and radius $a$. 

Density Profile
$$
\rho(r) = \frac{M}{2\pi} \frac{a}{r} \frac{1}{(r+a)^3}
$$
Mass profile
$$
M(r) = M \frac{r^2}{(r+a)^2}
$$
Potential:
$$
\Phi(r) = - \frac{GM}{r + a}
$$
The acceleration is then
$$
\vec{a} = - \hat{r} \frac{G\,M}{\left(r + a\right)^2}
$$

Circular velocity

## Miyamoto-Nagai (Bulge)

$$
\Phi(r) = - \frac{M}{\sqrt{R^2 + b^2}}
$$

(special case of $a=0$ for Miyamoto Nagai Disk)

## Miyamoto-Nagai (Disk)

for the MW disk (thin and thick), the potential is derived from @miyamoto1975, 
$$
\Phi(R, z) = \frac{-GM}{\left(R^2 + \left[a + \sqrt{z^2 + b^2}\right]^{2}\right)^{1/2}}
$$
for convenience, we define the helper variables, $S = \sqrt{z^2 + b^2}$ and $D = \sqrt{R^2 + (a+S)^2}$, so the potential simplifies to 
$$
\Phi(R, z) = \frac{-GM}{D}
$$

$$
a_{\rm disk} = -\nabla \Phi = -\left(\hat{r}\frac{\partial}{\partial r} + \hat{z} \frac{\partial}{\partial{z}}\right) \Phi
$$

$$
a_{\rm disk} = -\frac{G\,M}{D^3}\,R\,\hat{R} - \frac{G\,M}{D^3} \left(\frac{a}{S}+1\right) \,z\,\hat{z}
$$



density
$$
\rho = \frac{b^2 M}{4\pi} \frac{aR^2 + (a+3r)(a + r)^2}{S^3 D^5}
$$




## Powerlaw Bulge (Vasily)


$$
\rho	\propto (1 + r/r_b)^{-\Gamma} \exp(-(r/u_b)^2)
$$



## Logarithmic

$$
\Phi = v^2/2 \ln\left(x^2 + y^2 + (z/q)^2 + d^2\right)
$$

acceleration
$$
\vec{a} = -\frac{v^2}{x^2 + y^2 + (z/q)^2 + d^2} 
\left(
	\vec{x} + \vec{y} + \vec{z} / q^2
\right)
$$




## Isothermal (disk)

$$
\rho	\sim \exp(-R/R_d)\ {\rm sech}^2(z/(2h))
$$



## $\alpha\beta\gamma$ profile (bulge)

 @zhau1996.


$$
\rho = r^{-\gamma} (1 + r^{\alpha})^{(\beta - \gamma)/\alpha}
$$


potential is compl

## $\alpha\beta\gamma$​​ profile (bulge?)


$$
\rho_h = (s / r_h)^{-\gamma} (1 + (s/r_h)^\alpha)^{(\gamma - \beta)/\alpha} \exp(-(s/u_h)^\eta)\\
s \equiv (pq)^{1/3} \sqrt{x^2 + (y/p)^2 + (z/q)^2}
$$

where $s$ is sphericalized (or hyperelliptical) radius 

with exponential drop off $u_h = 200$, $\eta=2$. Becomes as published with $u_h \to \infty$



# Milky Way Potentials





instead the @vasily+2021 potential is

| Component         | Values                                      |      |
| ----------------- | ------------------------------------------- | ---- |
| spherical bulge   | $\Gamma=1.8$, $r_b=0.2$, $u_b=1.8$, $M=1.2$ |      |
| disk (isothermal) | M=5, r_d=3, h_d=0.4                         |      |
| halo              |                                             |      |

## 



## Notes about derivations

For my sake, here are the equations to derive each of the following relations. 

Given a density profile, the contained mass is $M(r) = \int_0^r \rho dV$, which if assuming a spherical profile
$$
M(r) = 4\pi \int_0^r \rho\,r^2\,dr
$$
The potential is given by
$$
\Phi({\bf x}) = \iiint \frac{G\,\rho}{|{\bf x}- {\bf x}'|} d^3{\bf x'}
$$
For a spherica distribution (e.g. Galaxies book), the formula reduces to
$$
\Phi(r) = -4\pi G\left[\frac{1}{r}\int_0^r \rho(r')r'^2\,dr' + \int_r^\infty \rho(r') r' dr'\right]
$$
Note that the first term is simply $G M(r) / r$ so if we solve for $M(r)$ we only need one additional integral to evaluate $\Phi$. 



The acceleration is $a = -\nabla \Phi$. 





# 2D profiles

Given the three dimensional density, we can compute the 2D density with 
$$
\Sigma(R) = 2\int_R^{\infty} \rho(r)\ \frac{r}{\sqrt{r^2 - R^2}} \ dr
$$


## Exponential

$$
\Sigma(R) = \Sigma_0 \exp(-R/R_s)
$$


$$
M(R) =  \pi  R_s^2 \Sigma_0 \left(1 - x \exp(-x) - \exp(-x)\right)
$$


Solving $M(R) = 1/2M_{\rm tot}$,  the half mass radius is $1.6783\,R_s$​. 

for the 3D exponential, 
$$
M(r) = 4\pi r_s^3 \rho_0 \left(2 - 2e^{-x} - 2xe^{-x} - x^2\,e^{-x}\right)
$$
which gives $r_h = 2.674 r_s$. So, given $R_s$, $r_s = 0.62762\,R_s$. (Verify this!)
