import DensityEstimators: histogram 
using DataFrames
import LilGuys as lguys
using Turing
using CairoMakie
using Arya
import Distributions: pdf, Normal

using Measurements: ±
import Measurements

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
    model_vel_1c(x, xerr)

Fits a 1c model to a set of velocities and uncertainties
"""
@model function model_vel_1c(x, xerr; μ_s_prior=100, μ_0_prior=0)
    μ ~ Normal(μ_0_prior, μ_s_prior)
    σ ~ Uniform(0, 20)
	s = @. sqrt(σ^2 + xerr^2)

	x ~ MvNormal(fill(μ, length(x)), s)
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
    plot_samples!(samples, x; thin, color, alpha kwargs...)

Plot the MCMC samples (dataframe with μ and σ)  as gaussians
across a range x.

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

    θ_pm = @. atand(xi, eta)
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
    df[!, :L_RV_SAT] = L_RV_SAT(df.v_z, df.v_z_err, radial_velocity_gsr, sigma_v)
    df[!, :L_RV_BKD] = L_RV_BKD(df.v_z, σ=σ_halo)
    df[!, :PSAT_RV] = PSAT_RV(df, f_sat)

    df
end


"""
    new_pm(rv_meas, obs_props)

Use the observed properties to provide a PM prior for members,
wei
"""
function new_pm(rv_meas, obs_props)
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



function get_error(df, key)
    if key*"_em" ∈ keys(df)
        return max(df[key*"_em", key*"_ep"])
    elseif key*"_err" ∈ keys(df)
        return df[key * "_err"]
    else
        @error "error for $key not found in dict"
    end
end


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

function sem_inv_var_weights(weights)
    return 1 / sqrt(sum(weights))
end

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
