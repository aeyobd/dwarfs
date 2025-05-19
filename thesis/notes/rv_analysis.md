## Radial velocity modeling

For Scl and UMi, to derive robust RV measurements, we compile a variety of literature radial velocity data and compute the radial velocity and velocity dispersion. 

For Scl, we present a combined sample crossmatching the studies [@sestito+23], apogee (c.o. Federico), @walker+2009, and @tolstoy+23, containing a total of 2466 stars with measurements. Note that @battaglia+2009 is superceded by @tolstoy+23. If multiple samples measure the same star, then we take the weighted mean of the measurements. We flag stars with a measurement dispersion greater than 5$\sigma$ of the combined measurement error and stars outside of the range 60-150 km/s as these are likely too far off to be members (and the intrinsic velocity dispersion and systematic errors are the dominant contributors to dispersion.)

We estimate the combined standard deviation and uncertainty of a given star with
$$
\sigma_v^2 = \sum_i n_i\, \sigma_{v,i}^2 + (\mu_{v,i} - \mu_v)^2
$$
the probability that the star is a binary is then the probability that the standard deviation is significant. 



Note that there is, in detail, small projection effects to the radial velocities, especially as we move towards 1 degree from the centre of the galaxy. We define the perpendicular tangent plane velocity $v_z$ to be 
$$
v_z = v_{\rm los,\ gsr} \cos\theta + v_{\rm tan, rad} \sin\theta
$$
where $v_{\rm tan}$ is the tangental systematic proper motion of the system. Because most candidate members have large proper motion uncertainties, if the proper motion is within $3\sigma$ of the systemic proper motion (num for Scl and Num for UMi), we average the observed and systemic proper motions to reduce the noise in $v_z$. Note that this likely biases the velocity dispersion slightly, but given the uncertainty, the effect would be minimal. The tangental component only contributes ~1\% at most of $v_z$. 



To fit the systematic radial velocity, we use the Bayesian model
$$
v_z \sim N(v_0, \sigma_v + \delta v_z)\\
v_0 \sim U \\
\sigma_v \sim N?
$$
for Sculptor, and only changing the priors above for Ursa Minor
$$
v_0 = \\
\sigma_v = \\
$$




To test if there is a spatial gradient in velocity, we construct a new model
$$
v_z = N(v_0 + a \xi + b \eta, \sigma_v + \delta v_z)
$$




![Sculptor radial velocity observations](figures/scl_v_z_obs_scatter.png) 

![Sculptor radial velocity observations](figures/umi_v_z_obs_scatter.png)

