# Milky Way Potentials

Following @borukhovetskaya2022, we define a multi-component MW potential as follows. 



| name      | component | reference          | potential                     |
| --------- | --------- | ------------------ | ----------------------------- |
| Hernquist | Bulge     | @hernquist1990     | $\Phi(r) = - \frac{G M}{r+a}$ |
| NFW       | Halo      | @nfw1996, @nfw1997 | $\Phi(r) \propto \ln(1+x)/x$  |
|           | Disk      | @miyamoto1975      |                               |
|           |           |                    |                               |



| Name      | Acceleration                                | Contained Mass         | Circular Velocity |
| --------- | ------------------------------------------- | ---------------------- | ----------------- |
| Hernquist | $\hat{a}(r) = -\hat{r} \frac{G M}{(r+a)^2}$ | $M(r) = M r^2/(r+a)^2$ |                   |
|           |                                             |                        |                   |
|           |                                             |                        |                   |





| Component  | Values                    |      |
| ---------- | ------------------------- | ---- |
| thin disk  | M=5.9, a=3.9, b=0.31      |      |
| thick disk | M=2, a=4.4, b=0.92        |      |
| bulge      | M = 2.1, a = 1.3          |      |
| halo       | Mvir=115, r=20.2, c=9.545 |      |

## Notes about derivations

For my sake, here are the equations to derive each of the following relations. 

Given a density profile, the contained mass is $M(r) = \int_0^r \rho dV$. The acceleration is $a = -\nabla \Phi$. 

## Bulge

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

## Disk

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



## NFW Halo

### Aside:Defined quantities

The NFW is parameterized in terms of 

with $c=r_{200}/r_s$

The NFW halo is written in terms of $M_{200}$, $R_{200}$. These are the radius and mass of a sphere with 200 times the average critical density, $\rho_{crit} = 3H^2/8\pi G=127.35M_\odot/{\rm kpc}^3 $, so
$$
\rho_{200} = 200\rho_{crit} = \frac{M_{200}}{4/3\ \pi R_{200}^3}
$$

$$
R_{200} = \sqrt[3]{\frac{1}{200}\frac{3M_{200}}{4\pi \rho_{\rm crit}}}
$$


$$
M_{200} = \frac{4\pi}{3} R_{200}^3\ \rho_{200}
$$


nother quantity which turns up often here is
$$
A(x) \equiv \log (1+x) - \frac{x}{1+x}
$$


We will also just define now
$$
x \equiv r/r_s
$$
as all quantities are simply scaled by the scale radius.

For example, Asya uses $M_{200} = 1.04e10$ and $c=12.5$ for Fornax, giving $M_s = 0.119\times10^{10}\,{\rm M}_\odot$ and $R_{s}=3.68\,$kpc.

### Density

From the @nfw1996 paper, 
$$
\rho(x) =  \frac{\rho_s/3}{x\ (1+x)^2}
$$
where
$$
\rho_s = \frac{c^3}{A(c)} \rho_{200}
$$
and $A(c)$ is as above. The characteristic density can also be written in terms of scale mass,  $M_s = M_{200}/{A(c)}$  (see below), giving
$$
\rho_s =  \frac{M_s}{4/3\ \pi\, {r_s}^3}
$$
### Contained mass

By integrating $\rho$, the contained mass interior to $r$ is
$$
M(x) = M_s\ A(x)
$$
where
$$
M_s \equiv 4\pi/3 \rho_s r_s^3 = \frac{M_{200}}{A(c)}
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

if $x = r/r_s$, then 
$$
\left(V_c/V_{200}\right)^2 = \frac{A(x)/x}{A(c)/c}
$$

from the radius and maximum circular velocity, 
$$
R_{\rm circ}^{\rm max} = \alpha\ r_s
$$

where $\alpha\approx2.16258$

# Misc



## Numerical Setup

$$
h_{grav} = \frac{0.005}{R_{\rm max}} \left(\frac{N}{10^6}\right)^{-0.5}
$$

10pc for $10^7$ particles

## Units

| Quantity   | Units               |
| ---------- | ------------------- |
| G          | 1                   |
| Mass       | $10^{10}$ M$_\odot$ |
| radius     | 1 kpc               |
| velocities | 207.4 km/s          |
| time       | 4,715,000 yr        |

# Tests

## Conservation laws

- Momentum: $\sum m\,v_i={\rm constant}$ 
- Angular momentum: $\sum m\,r\times v = {\rm constant}$ 
- Energy
  - $T = \sum \frac{1}{2} m v^2$
  - $ U = \frac{1}{2} \sum_i \sum_{j\neq i} \frac{G m_i m_j}{|x_i - x_j|} = \frac{1}{2} \sum_i m_i \Phi(x_i)$. The factor of 1/2 is so we only count each pair once.
  - $U_{\rm ext} = \sum_i m_i \Phi_{\rm ext}(x_i) $
  - $T + U + U_{\rm ext} = {\rm constant}$

## Two body

Here, I just test a keplarian 2-body system. The solution is the same for the one-body system, so we expect Kepler's laws to hold--the orbit should be elliptical with a period
$$
P^2 = \frac{4\pi^2}{\mu} a^3
$$
where $\mu = G(m+m)$ is the reduced gravitational mass of the system, and $a$ is twice (?) the semi-major axis of one of the orbits.



## Three body

To test the three body problem, I just use the special case (...)



## Hernquist  and NFW potential

For a spherical potential $\Phi$, we know that the roots of 
$$
0=\frac{1}{r^2} + 2\frac{\Phi(r) - E}{L^2}
$$
provide us with the maximum and minimum orbital radius, $r_{\rm min}$ and $r_{\rm max}$. We can also calculate the period in $r$ with 
$$
T_{r} = \int_{r_{\rm min}}^{r_{\rm max}} \frac{2}{\sqrt{2|E-\Phi(r)| - \frac{L^2}{r^2}}}\ dr
$$
I verify these relations for the hernquist and NFW spherical potentials. 

## Disk potential

Disk potentials are more complex but we can still check the overal relationships .




# Analyzing observations

Observationally, we know the 6D kinimatic and have estimated total stellar mass 

## Orbits



## Initial profile parameters

to get the stellar mass, we need the absolute magnitude and stellar mass to light ratio. For example, Sculptor has an absolute magnitude of $M_V = -11.1$ [@mcconnachie2012], and the sun has an absolute magnitude of 4.83, so with a stellar mass to light ratio of 1.6, this gives a total stellar mass of $M_\star=3.8\times10^6\,$M$_\odot$ .

Next, we use that
$$
M_\star = m_0\,\nu^\alpha \exp({-\nu^\gamma})
$$
where $m_0 = 3\times10^8 M_\odot$, $\alpha=3.36$, and $\gamma=-2.4$. $\nu = V_{\rm max}/$50km/s, from an emperical fit to @fattahi2018. For the example of Sculptor, $\nu=0.643$ solves the above equation, so $V_{\rm max}=32.2$km/s. 

Finally, the characteristic radius is from @ludlow2016 (see their appendix C). The full equations are complex, but briefly, from a given M and z, you can calculate their fit to the expected concentration. We then need to solve for a value of $M_s$ which gives a value of $c$ which together predict the correct $V_{\rm circ, max}$. 

Using that the maximum circular velocity occurs at $\alpha r_s$ where $\alpha = 2.163$ (is a numerical solution), and the equation for circular velocity (way back above but repeated here:)
$$
\left(V_c/V_{200}\right)^2 = \frac{c}{r/r_s}\frac{A(r/r_s)}{A(c)}
$$
where $V_{200} = \sqrt{G\ M_{200} / r_{200}}$. The solution here is $M_{200} = 5.4e9$ where $c=13.0$, and $r_s = 2.86$ and $M_s = 6.12e8$.  







