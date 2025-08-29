# MCMC Histograms

The J+24 model has some limitations, especially for background dominated systems. 

- J+24 assume that all dwarfs can be described by a single or double exponential
- This has the result that, especially for probability mass, the density profiles may be confidently reproducing the input profile.
- We would like to alleviate this assumption while retaining a bayesian framework



## Method

We extend the J+24 models using the following framework.



The total likelihood of membership of a star to either the background or the satellite is the product of CMD, PM, and spatial likelihood functions. In essence, the CMD/PM likelihoods are used as a prior on the likelihood of membership, and the spatial weighs into the final posterior likelihoods. Given the total likelihood of the satellite, the probability of membership is
$$
P_{\rm sat} = \frac{f_{\rm sat}{\cal L}_{\rm sat}}{f_{\rm sat}{\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}}
$$
Note that the prior likelihood is given by ${\cal L}_{\rm sat} = {\cal L}_{\rm cmd, sat}{\cal L}_{\rm pm, sat}$, the  and the denominator represents the marginal likelihood.
$$
{\cal L} = {\cal L}_{\rm sat} + (1-f_{\rm sat}){\cal L}_{\rm bg}
$$


The total (posterior) likelihood is
$$
{\cal L} = {\cal L}_{\rm pm} {\cal L}_{\rm cmd}  {\cal L}_{\rm space}
$$


However, note that the above formulae are framed in a per-star framework. For the density profile, the likelihood is instead given by
$$
{\cal L}_{\rm tot} = \prod_i {P}_{\rm s, sat} P_{\rm sat}^{\rm prior} + {P}_{\rm s, bg} (1-P_{\rm sat}^{\rm prior})
$$
By calculating the posterior distribution of the spatial likelihood as a function of radius, we can determine the 

## Histogram models

To parameterize the histogram model, we do
$$
{\cal L}_{\rm space} / {\cal L}_{\rm bg} = \Sigma
$$


where $\Sigma$ is a piecewise constant density profile 



## Including structural uncertainties 

For this model, we first chose bins for the data. We use a combination of equal width and equal number to allow for reduction of noise.

Bins are taken to be the larger of

- Bin containing next $n$ stars. We adpot the herusitic

$$
n =  2N^{2/5} / f
$$

where $f$ is the fraction of member stars within 1/3 half-light radii of the centre. 

- A bin of width $w$ in log R. We adopt $w=0.05$

Next, we construct the hierarchical bayesian model as
$$
\begin{cases}
P(\Sigma | \{S_i\}) \sim \prod_i\Sigma(R_i)\, P_{\rm sat}(S_i) + P_{\rm bg}(S_i)\\
\Sigma = {\rm piecewise}({\rm bins}, \Sigma_b) \\
\log \Sigma_b \sim U(-12, 6)\\
\end{cases}
$$
To account for structural uncertainties, we run a MCMC sampling the structural params 1024 times, creating new elliptical radii following 
$$
R_i = R_{\rm ell}(\xi_i + \delta\xi,  \eta_i+\delta\eta, e, \theta)\\
\theta \sim N(\theta_0, \delta\theta) \\
e \sim N(e_0, \delta e) \\
\delta\xi \sim N(0, \delta x)\\
\delta \eta \sim N(0, \delta x)
$$

## Results

Redo R_h using Sérsic fits...

Try generalized Sérsic fit.



| galaxy           | R_h (exp inner) | # candidates |
| ---------------- | --------------- | ------------ |
| Fornax           | $17.4\pm0.6$    | 23,154       |
| Leo I            | $3.6\pm0.2$     | 1,242        |
| Sculptor         | $9.3\pm0.2$     | 6,875        |
| Antlia 2         | 76.2pm7.2       | 649?         |
| Leo II           | $3.0\pm0.4$     | 347          |
| Carina           | $8.3\pm0.3$     | 2,389        |
| Draco            | $7.7\pm0.3$     | 1,781        |
| Ursa Minor       | $12.3 \pm 0.5$  | 2,122        |
| Canes Venatici I | 8.9?            | 156?         |
| Sextans I        | $20.0 \pm0.8$   | 1,830        |
| Crater 2         | 31.2?           | 507          |

Table: For faint dwarf galaxies

Note we skip the following galaxies mentioned in J+24



Faint dwarfs, onew which we could possibly rederive parameters for (>~ 30 stars)



| Galaxy         | # J+24 memb | R_h           | status |
| -------------- | ----------- | ------------- | ------ |
| Bootes 1       | **252**     | 11.26\pm0.27  | yy     |
| Bootes 3       | **123**     | 33\pm2.5      | yy     |
| Carina 2       | 69          | 8.68          | y      |
| Coma Berenices | 44          | 5.63          | n      |
| Draco 2        | 33          |               |        |
| Eridanus 4     | 47          | 4.9           | n      |
| Grus 2         | 53          | 6             | n      |
| Hercules       | 44          | 5.99\pm0.58   | n      |
| Hydrus 1       | **132**     | 7.42          | n      |
| Phoenix        | **207**     | 2.3\pm0.07    | n      |
| Reticulum 2    | 73          | 5.6\pm0.2     | n      |
| Sagittarius 2  | 72          | 1.7\pm0.05    | n      |
| Segue 1        | 37          | 3.95\pm0.48   | n      |
| Tucana 2       | 50          | 9.89          | n      |
| Tucana 3       | 83          | 6.0           | n      |
| Ursa Major 1   | 57          | 8.34 \pm 0.34 | n      |
| Ursa Major 2   | 69          | 13.95\pm0.46  | n      |





| Galaxy            | # J+24 memb |
| ----------------- | ----------- |
| Aquarius 2        | 15          |
| *Bootes 2*        | 24          |
| *Bootes 4* *      | 5           |
| *Bootes 5*        | 6           |
| Canes Venatici II | 15          |
| Carina 3          | 12          |
| Centaurus 1       | 29          |
| Cetus 2           | 8           |
| Cetus 3*          | 1           |
| Columba 1         | 9           |
| DELVE 2           | 9           |
| DESJ0225+0304*    | 2           |
| Eridanus 2        | 20          |
| Eridanus 3        | 3           |
| Grus 1            | 12          |
| Horologium 1      | 21          |
| Horologium 2*     | 5           |
| Hydra 2           | 21          |
| Indus 1*          | 3           |
| Leo 4             | 0           |
| Leo 5             | 9           |
| Leo Minor*        | 4           |
| Leo T             | 8           |
| Pegasus 3         | 2           |
| Pegasus 4         | 23          |
| Phoenix 2         | 11          |
| Pictor 2*         | 10          |
| Pictoris 1        | 12          |
| Pisces 2*         | 2           |
| Reticulum 3*      | 7           |
| Segue 2           | 27          |
| Triangulum 2      | 15          |
| Tucana 4          | 14          |
| Tucana 5*         | 6           |
| Virgo 1*          | 2           |
| Virgo 2*          | 2           |
| Willman 1         | 13          |

*: marked unreliable in J+24

Using the same methods above, we select members from J+24's stellar 





## Tests

- Does incorporating uncertainties with the MCMC matter?
- MCMC convergence & resolution?
- Recovery (long term project)



## Comparisons

In general, this method provides similar results to J+24. The primary differences are that the density profiles appear "smoother" and have larger, more representative uncertainties. The uncertainty quantification is critical for understanding limitations and the statistical significance of features in dwarf galaxies. However, the difference in uncertainties primary applies to either the innermost or regions where the density approaches the background density.

While for simplicity, we retain absolute probability cuts of the main sample in the main text, we note that these profiles likely better represent the observational knowledge of density profiles. 

In some cases, like Antlia 2, these methods produce substantially different density profiles. In a background dominated regime, the J+24 method can "hallucinate" a density profile. There is always some change that a foreground star will appear similar enough in Gaia to be associated with a galaxy. If the background density is much higher than the satellite density, then the resulting member density would be controlled primarly by the prior probability (given by the assumed density ratio at that region), however, it is likely that the stars being selected are just background stars which happened to stray into the satellite parameter space. This would result in a density profile exactly recovering the assumed spatial dependence.

