import DensityEstimators: histogram 
using DataFrames
import LilGuys as lguys
using Turing
using CairoMakie
using Arya
import Distributions: pdf, Normal, ccdf, Chisq

import KernelDensity
using Measurements: ±
import Measurements

P_CHI2_MIN = 0.001

@doc raw"add a and b in quadrature: $sqrt(a^2 + b^2)$"
⊕(a, b) = sqrt(a^2 + b^2)


"""
    xmatch(df1, df2, max_sep=2)

Cross matches 2 dataframes given angular seperation in arcmin.
Returns the filter of matched entries in df1 and the indicies of
the closest entry in df2 to each row in df1.
"""
function xmatch(df1::DataFrame, df2::DataFrame, max_sep=2)
	max_sep = max_sep / 3600
	dists = lguys.angular_distance.(df1.ra, df1.dec, df2.ra', df2.dec')

	idxs = [i[2] for i in dropdims(argmin(dists, dims=2), dims=2)]

	filt = dropdims(minimum(dists, dims=2), dims=2) .< max_sep
	return filt, idxs
end




"""
    sigma_clip(x, nσ=5)

Remove all values of x outsize of nσ standard deviations of the median,
progressively recalculating the standard deviation and median
at meach step
"""
function sigma_clip(x::AbstractVector{<:Real}, nσ::Real=5)
	dN = 1
	filt = .!isnan.(x)
	while dN > 0
		@info "iteration, removing $dN"
		σ = std(x[filt])
		μ = median(x[filt])
		N = sum(filt)
		filt .&= x .> μ - nσ * σ
		filt .&= x .< μ + nσ * σ
		dN = N - sum(filt)
	end
	
	return filt
end


"""
    plot_samples!(samples, x; thin=10, color=:black, alpha, kwargs...)

Plot the MCMC samples (dataframe with μ and σ)  as gaussians
across a range x. alpha defaults to 1/∛N
"""
function plot_samples!(samples, x;
		thin=10, color=:black, alpha=nothing, kwargs...)

	alpha = 1 / (size(samples, 1))^(1/3)
	for sample in eachrow(samples)[1:thin:end]
		y = lguys.gaussian.(x, sample.μ, sample.σ)
		lines!(x, y, color=color, alpha=alpha; kwargs...)
	end
end


"""
    plot_samples(measurements, samples; thin, bins=30)

Plot the samples and the RV measurements (as a histogram)..
See also `plot_samples!`
"""
function plot_samples(measurements, samples; thin=16, bins=30)
    x_obs = Float64.(measurements.RV)
    x = LinRange(extrema(x_obs)..., 100)

    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "RV", ylabel = "density")

    h = histogram(x_obs, bins, normalization=:pdf)
    plot_samples!(samples, x, thin=thin)
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

    fig

end


"""
    summarize(chain; p=0.16)

Summarizes the MCMC chain. 
Differing from turing in that we add the mean
and p-value quantile interval
"""
function summarize(chain; p=0.16, pp=0.025)
    df = Turing.summarize(chain) |> DataFrame
    Nvar = size(df, 1)
    lows = zeros(Nvar)
    lowlows = zeros(Nvar)
    mids = zeros(Nvar)
    highs = zeros(Nvar)
    highhighs = zeros(Nvar)

    for i in 1:Nvar
        x = chain[:, i, :].data
        llow, low, mid, high, hhigh = quantile(vec(x), [pp, p, 0.5, 1-p, 1-pp])
        lows[i] = low
        mids[i] = mid
        highs[i] = high
        lowlows[i] = llow
        highhighs[i] = hhigh
    end

    df[!, :median] = mids
    df[!, :error_lower] = mids .- lows
    df[!, :error_upper] = highs .- mids
    df[!, :error_lower_2] = mids .- lowlows
    df[!, :error_upper_2] = highhighs .- mids

    return df
end


"""
    filt_missing(col; low=-Inf, high=Inf)

Remove missing values from the array and exclude values outside of low and high.
"""
function filt_missing(col; low=-Inf, high=Inf, verbose=true)
	filt =  @. !ismissing(col) && !isnan(col)
	filt1 = high .> col .> low
    if verbose
        @info "excluding $(sum(.!(filt1)[filt])) outliers"
    end

	return filt .& filt1
end


@doc raw"""
    rv_correction(ra, dec, ra0, dec0, pmra, pmdec, distance)

Compute the correction to add to GSR radial velocities to calculate
$v_z$, the velocity in the direction parallel to the RV of the centre.
"""
function rv_correction(ra, dec, ra0, dec0, pmra, pmdec, distance)
    vra = lguys.pm2kms(pmra, distance)
    vdec = lguys.pm2kms(pmdec, distance)

    ϕ_pm = lguys.angular_distance.(ra, dec, ra0, dec0)

    xi, eta = lguys.to_tangent(ra, dec, ra0, dec0)

    θ_pm = @. atand(xi, eta) # PA relative to centre
    rv = 0 # want correction to vz

    sθ = @. sind(θ_pm)
    cθ = @. cosd(θ_pm)
    sϕ = @. sind(ϕ_pm)
    cϕ = @. cosd(ϕ_pm)

    vR = @. vra * sθ + vdec * cθ
    vθ = @. vra * cθ - vdec * sθ
    rv = 0 # just want shift to vz

    vx = @. vθ*cθ + vR*sθ*cϕ + rv*sθ*sϕ
    vy = @. -vθ*sθ + vR*cθ*cϕ + rv*cθ*sϕ
    vz = @. rv * cϕ - vR * sϕ

    vtot1 = @. sqrt(vx^2 + vy^2 + vz^2)
    vtot2 = @. sqrt(rv^2 + vra^2 + vdec^2)
    @assert all(isapprox.(vtot1, vtot2)) # make sure length is preserved

    return vz
end


@doc raw"""
    rv_correction(ra, dec, coord)

Compute the correction to add to GSR velocities (if coord is GSR)
or ICRS velocities (if coord is ICRS) to move to satellite $v_z$ frame.
"""
function rv_correction(ra, dec, coord::lguys.AbstractSkyCoord)
    ra0 = Measurements.value(coord.ra)
    dec0 = Measurements.value(coord.dec)
    return rv_correction(ra, dec, ra0, dec0, 
                         coord.pmra, coord.pmdec, coord.distance)
end


"""
    rv_gsr_shift(ra, dec)

Compute the apparent radial velocity of a stationary object at ra0, dec0
due to the solar motion.
"""
function rv_gsr_shift(ra0::Real, dec0::Real)
    rv_offset = lguys.transform(lguys.ICRS, lguys.GSR(ra=ra0, dec=dec0, radial_velocity=0)).radial_velocity
end


"""
    L_RV_SAT(RV, RV_err, μ, σ)

Compute the likelihood of a star belonging to the satellite 
based on its radial velocity, assuming a normal distribution.
""" 
function L_RV_SAT(RV::Real, RV_err::Real, μ_v::Real, σ_v::Real)
    σ = sqrt(σ_v^2 + RV_err^2)
    return pdf(Normal(μ_v, σ), RV)
end


"""
    L_RV_BKD(rv::Real, ra0, dec0)

Compute the likelihood of the velocity assuming 
a halo velocity dispersion of 100 km/s and mean velocity
stationary given the gsr radial velocity.
"""
function L_RV_BKD(rv::Real; σ=100)
    return pdf(Normal(0, σ), rv)
end


"""
    PSAT_RV(rv_meas, f_sat)

Compute the RV-informed satellite probability.
"""
function PSAT_RV(rv_meas::DataFrame, f_sat)
    L_SAT = @. rv_meas.L_CMD_SAT * rv_meas.L_PM_SAT * rv_meas.L_S_SAT * rv_meas.L_RV_SAT

    L_BKD = @. rv_meas.L_CMD_BKD * rv_meas.L_PM_BKD * rv_meas.L_S_BKD * rv_meas.L_RV_BKD

    return @. f_sat * L_SAT / (f_sat * L_SAT + (1-f_sat) * L_BKD)
end


"""
    add_PSAT_RV!(df; sigma_v, radial_velocity_gsr, f_sat, σ_halo=100)

Add PSAT column including radial velocity information and likelihoods
`L_RV_SAT` and `L_RV_BKD`
"""
function add_PSAT_RV!(df::DataFrame; sigma_v::Real, radial_velocity_gsr::Real, f_sat::Real, σ_halo=100)
    df[!, :L_RV_SAT] = L_RV_SAT.(df.vz, df.vz_err, radial_velocity_gsr, sigma_v)
    df[!, :L_RV_BKD] = L_RV_BKD.(df.vz, σ=σ_halo)
    df[!, :PSAT_RV] = PSAT_RV(df, f_sat)

    df
end


"""
    add_gsr!(df::DataFrame; distance, pmra, pmdec)

Add the GSR corrected proper motions and radial velocities to the dataframe.
Requires columns `ra`, `dec`, `RV`.
"""
function add_gsr!(df::DataFrame; distance, pmra, pmdec)
    obs = [lguys.ICRS(ra=row.ra, dec=row.dec, distance=distance,
                      pmra=row.pmra, pmdec=row.pmdec, radial_velocity=row.RV)
           for row in eachrow(df)
          ]

    obs_gsr = lguys.to_frame(lguys.transform.(lguys.GSR, obs))


    df[!, :pmra_gsr] = obs_gsr.pmra
    df[!, :pmdec_gsr] = obs_gsr.pmdec
    df[!, :radial_velocity_gsr] = obs_gsr.radial_velocity

    df
end


"""
    _new_pm(rv_meas, obs_props)

Use the observed properties to provide a PM prior for members. Not used currently
"""
function _new_pm(rv_meas, obs_props)
    σ_pm = LilGuys.kms2pm(obs_props["sigma_v"], obs_props["distance"])
    σ_ra = obs_props["pmra_err"] ⊕ σ_pm
    σ_dec = obs_props["pmdec_err"] ⊕ σ_pm

    pmra = map(eachrow(rv_meas)) do row
        xs = [row.pmra, obs_props["pmra"]]
        ws = [1/row.pmra_err^2, 1/σ_ra^2 * row.PSAT_RV]

        m = lguys.mean(xs, we)
        err = 1/sqrt(sum(ws))
        return m ± err
    end

    # ╔═╡ ad7ae1a8-c5a4-4d79-bd86-35431fa71856
    pmdec = map(eachrow(rv_meas)) do row
        xs = [row.pmdec, obs_props["pmdec"]]
        ws = [1/row.pmdec_err^2, 1/σ_dec^2 * row.PSAT_RV]

        m = lguys.mean(xs, we)
        err = 1/sqrt(sum(ws))
        return m ± err
    end

    return pmra, pmdec
end


"""
    get_error(df, key)

Retrieve the uncertainty of the key in the given dict (either max of key_em and key_ep or key_err). 
"""
function get_error(df, key)
    if key*"_em" ∈ keys(df)
        return max(df[key*"_em"], df[key*"_ep"])
    elseif key*"_err" ∈ keys(df)
        return df[key * "_err"]
    else
        @error "error for $key not found in dict"
    end
end


"""
    icrs(df::AbstractDict; add_sigma_pm_int)

Given the dictionary of observed properties, 
returns the ICRS coordinate with reported uncertanties, 
adding the velocity dispersion to the PM uncertainty. 

This formulation is most helpful for this analysis where we 
assume a constant PM for the satellite
"""
function icrs(df::AbstractDict; add_sigma_pm_int=true)
    σ_pm = lguys.kms2pm(df["sigma_v"], df["distance"])
    @info "σ_pm = $σ_pm"

    kwargs = Dict(
              Symbol(key) =>  df[key]  ± get_error(df, key) for key in ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity"]
             )

    if add_sigma_pm_int
        for k in ["pmra", "pmdec"]
            kwargs[Symbol(k)] = df[k] ± (get_error(df, k) ⊕ σ_pm)
        end
    end

    return lguys.ICRS(;kwargs...)
end


"""
    add_rv_means!(all_stars, all_studies)

Add the RV means, errors, std, and counts by averaging
over the columns RV_study, RV_err_study, RV_sigma_study, and RV_count_study.
"""
function add_rv_means!(all_stars, all_studies)

	for col in [:RV, :RV_err, :RV_sigma]
		all_stars[!, col] .= 0.
		allowmissing!(all_stars, col)
	end
	
	all_stars[!, :RV_count] .= 0
	all_stars[!, :RV_nstudy] .= 0
	
	for (i, row) in enumerate(eachrow(all_stars))
		xs = [row["RV_$study"] for study in all_studies]
		xs_err = [row["RV_err_$study"] for study in all_studies]
		xs_sigma = [row["RV_sigma_$study"] for study in all_studies]
		xs_count = [row["RV_count_$study"] for study in all_studies]

		m, m_err, m_sigma, m_count = safe_weighted_mean(xs, xs_err, xs_sigma, xs_count)
		n_study = sum(.!ismissing.(xs))


		all_stars[i, :RV] = m
		all_stars[i, :RV_err] = m_err
		all_stars[i, :RV_sigma] = m_sigma
		all_stars[i, :RV_count] = m_count
		all_stars[i, :RV_nstudy] = n_study
	end

	all_stars
end


@doc raw"""
    sem_inv_var_weights(weights)

Standard error of mean given inverse variance weights.

For observations $x_i$, the uncertanty on the mean
$\bar x$ is 
``
    \delta \bar x = \sqrt{\frac{1}{\sum_i w_i}}
``
given that $w_i = 1/\delta x_i^2$.
"""
function sem_inv_var_weights(weights)
    return 1 / sqrt(sum(weights))
end


"""
    save_weighted_mean(values, errors, sigmas, counts)

Compute the weighted mean, stderr, stdev, and counts 
given vectors of measurements with provided values, 
stderrors, sigmas (not used), and counts per measurement.
Skip over missing values, returning missings if no nonmissings.
"""
function safe_weighted_mean(values, errors, sigmas, counts)
	filt = filt_missing(values, verbose=false)
	filt .&= filt_missing(errors, verbose=false)
	if sum(filt) == 0
		return missing, missing, missing, 0
	end

	w = disallowmissing(1 ./ errors[filt] .^ 2)
	x = disallowmissing(values[filt])
	n = disallowmissing(counts[filt])

	x_mean = lguys.mean(x, w)
    std_err = sem_inv_var_weights(w)
	sigma = lguys.std(x, w)
	counts_tot = sum(n)
	
	return x_mean, std_err, sigma, counts_tot 
end


@doc raw"""
    model_vel_1c(x, xerr; <keyword arguments>)

Fits a 1c model to a set of velocities and uncertainties

# Arguments
- `μ_0_prior=0` prior mean on $\mu$
- `μ_s_prior=0` prior std on $\mu$
- `σ_max` maximum value of $\sigma$ (uniform prior)
"""
@model function model_vel_1c(x, xerr; μ_0_prior=0, μ_s_prior=100, σ_max=20)
    μ ~ Normal(μ_0_prior, μ_s_prior)
    σ ~ Uniform(0, σ_max)
	s = @. sqrt(σ^2 + xerr^2)

    m = fill(μ, length(x))
	x ~ MvNormal(m, s)
end



@doc raw"""
    model_vel_sigma_R(vz, vz_err, R_ell; <kwargs>)

Fits a 1c model to a set of velocities and uncertainties

# Args
- `μ_0_prior=0` prior mean on $\mu$
- `μ_s_prior=100` prior std on $\mu$
- `σ_max` maximum value of $\sigma$ (uniform prior)
- `R0` radius where $\sigma$ is not affected by `dlσ_dlR`.
"""
@model function model_vel_sigma_R(vz, vz_err, R_ell; μ_0_prior=0, μ_s_prior=100, σ_max=20, R0=10)
    log_R_ell = log10.(R_ell ./ R0)
    μ ~ Normal(μ_0_prior, μ_s_prior)
    σ ~ Uniform(0, σ_max)
    dlσ_dlR ~ Normal(0, 0.3)

    s_int = σ .* 10 .^ (dlσ_dlR * log_R_ell)
	s = @. sqrt(s_int^2 + vz_err^2)

    m = fill(μ, length(vz))
	vz ~ MvNormal(m, s)
end



@doc raw"""
    model_vel_gradient(vz, vz_err, ξ, η; <keyword arguments>)
Fits a normal (gaussian) distribution to 3d data with errors (to include spatial gradient)

# Arguments
- `μ_0_prior=0` prior mean on $\mu$
- `μ_s_prior=0` prior std on $\mu$
- `σ_max` maximum value of $\sigma$ (uniform prior)
- `rv_grad_s=0.1` Prior std on gradient in xi and eta (in units of km/s/arcmin).
"""
@model function model_vel_gradient(x, xerr, ξ, η; μ_0_prior=0, μ_s_prior=100, σ_max=20, rv_grad_s=0.1)
    μ ~ Normal(μ_0_prior, μ_s_prior)
    σ ~ Uniform(0, σ_max)

    A ~ Normal(0, rv_grad_s)
    B ~ Normal(0, rv_grad_s)

	m = @. μ + A*ξ + B*η
	s = @. sqrt(σ^2 + xerr^2)

	x ~ MvNormal(m, s)
end

@doc raw"""
    model_vel_gradient_both(vz, vz_err, ξ, η, R_ell; <keyword arguments>)
Fits a normal (gaussian) distribution to 3d data with errors (to include spatial gradient)

# Arguments
- `μ_0_prior=0` prior mean on $\mu$
- `μ_s_prior=0` prior std on $\mu$
- `σ_max` maximum value of $\sigma$ (uniform prior)
- `rv_grad_s=0.1` Prior std on gradient in xi and eta (in units of km/s/arcmin).
- `R_h` the half light radius.
"""
@model function model_vel_gradient_both(x, xerr, ξ, η, R_ell; μ_0_prior=0, μ_s_prior=100, σ_max=20, rv_grad_s=0.1, R_h)
    μ ~ Normal(μ_0_prior, μ_s_prior)
    σ ~ Uniform(0, σ_max)

    dlσ_dlR ~ Normal(0, 0.3)
    A ~ Normal(0, rv_grad_s)
    B ~ Normal(0, rv_grad_s)

	m = @. μ + A*ξ + B*η

    log_R_ell = log10.(R_ell ./ R_h)
    s_int = σ .* 10 .^ (dlσ_dlR * log_R_ell)
	s = @. sqrt(s_int^2 + xerr^2)

	x ~ MvNormal(m, s)
end

"""
    bayes_evidence(model, df_samples, arg)
    bayes_evidence(model, df_samples, args)

Compute the bayes evidence of a nested model given the Turing model 
and the samples. May give a single argument or a list of more than 
1 argument. Assumes that the base model is recovered when the parameters
are zeroed.
"""
function bayes_evidence(model, df_samples, arg::String)
	kde = KernelDensity.kde(df_samples[!, arg])
    Nsamples = size(df_samples, 1)
	df_prior = sample(model, Prior(), Nsamples) |> DataFrame
	kde_prior = KernelDensity.kde(df_prior[!, arg])

	return log(pdf(kde, 0.) / pdf(kde_prior, 0.))
end



function bayes_evidence(model, df_samples, args::AbstractVector{String})
	xs = tuple([df_samples[!, arg] for arg in args]...)
	kde = KernelDensity.kde(xs)

    Nsamples = size(df_samples, 1)
	df_prior = sample(model, Prior(), Nsamples) |> DataFrame
	
	xs = tuple([df_prior[!, arg] for arg in args]...)
	kde_prior = KernelDensity.kde(xs)

	x0 = zeros(length(args))
	return log(pdf(kde, x0...) / pdf(kde_prior, x0...))
end


"""
    prob_chi2(df_rv_meas)

Returns the p-values of the significance of the chi2 values for
each row in the DF, 
assuming
- `RV_err` weighted SEM on RV uncertainty
- `RV_sigma` weighted standard deviation of measurements
- `RV_count` number of observations.

In detail, this is a chi2 test where,
```

```
"""
function prob_chi2(df_rv_meas)
    return map(eachrow(df_rv_meas)) do row
        prob_chi2(row.RV_sigma, row.RV_err, row.RV_count)
    end
end


"""
    prob_chi2(sigma, sem, count)

Returns the p-values of the significance of the chi2 values for
each row in the DF, 
assuming
- `RV_err` weighted SEM on RV uncertainty
- `RV_sigma` weighted standard deviation of measurements
- `RV_count` number of observations.

In detail, this is a chi2 test where,
```

```
"""
function prob_chi2(sigma, sem, count)
    if count < 2
        return NaN
    end

    chi2 = sigma^2 / sem^2
    ν = count - 1
    ccdf(Chisq(ν), chi2)
end



"""
    filt_chi2(prob_chi2s)

Fiducial cut on chi2 (probability of 0.001).
"""
function filter_chi2(p_chi2)
    return @. isnan(p_chi2) || p_chi2 > 0.001
end



"""
Retrieve the value of F_BEST if source_id exists, otherwise return missing.
"""
function get_f_best(j24, source_id)
    if source_id |> ismissing
        return missing
    end
	if source_id ∉ j24.source_id
		return missing
	end
	return j24.F_BEST[j24.source_id .== source_id] |> only
end
