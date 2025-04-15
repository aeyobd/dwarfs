import DensityEstimators: histogram 
using DataFrames
import LilGuys as lguys
using Turing
using CairoMakie
using Arya
import Distributions: pdf, Normal


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
    model_vel_1c(x, xerr; μ_min, μ_max)

Fits a 1c model to a set of velocities and uncertainties
"""
@model function model_vel_1c(x, xerr; μ_min=90, μ_max=120)
	μ ~ Uniform(μ_min, μ_max)
	σ ~ LogNormal(2.5, 1) # approx 1 - 100 km / s, very broad but should cover all
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
function summarize(chain; p=0.16)
    df = Turing.summarize(chain)
    Nvar = size(df, 1)
    lows = zeros(Nvar)
    mids = zeros(Nvar)
    highs = zeros(Nvar)

    for i in 1:Nvar
        x = chain[:, i, :].values
        low, mid, high = quantile(vec(x), [p, 0.5, 1-p])
        lows[i] = low
        mids[i] = mid
        highs[i] = high
    end

    df[!, median] = mids
    df[!, error_lower] = mids .- lows
    df[!, error_upper] = highs .- mids

    return df
end


"""
    filt_missing(col; low=-Inf, high=Inf)

Remove missing values from the array and exclude values outside of low and high.
"""
function filt_missing(col; low=-Inf, high=Inf)
	filt =  @. !ismissing(col) && !isnan(col)
	filt1 = high .> col .> low
    @info "excluding $(sum(.!(filt1)[filt])) outliers"

	return filt .& filt1
end


"""
    to_orthoganal_velocities(gsr, gsr0)

Convert a list of GSR velocities to be in the orthoganal cartesian frame relative
to gsr0. 
Return
- vx: velocity in direction of Ra at centre
- vy: velocity in direction of DEC at centre
- vz: velocity in direction of dwarf's centre
"""
function to_orthoganal_velocities(gsr::Vector{<lguys.GSR}, gsr0::lguys.GSR)

    vra = [lguys.pm2kms(o.pmra, gsr0.distance) for o in gsr]
    vdec = [lguys.pm2kms(o.pmdec, gsr0.distance) for o in gsr]

    ϕ_pm = lguys.angular_distance.(rv_meas.ra, rv_meas.dec, ra0, dec0)

    xi, eta = lguys.to_tangent([g.ra for g in gsr], [g.dec for g in gsr], gsr0.ra, gsr0.dec)

    θ_pm = @. atand(rv_meas.xi, rv_meas.eta)

    sθ = @. sind(θ_pm)
    cθ = @. cosd(θ_pm)
    sϕ = @. sind(ϕ_pm)
    cϕ = @. cosd(ϕ_pm)

    vR = @. vra * sθ + vdec * cθ
    vθ = @. vra * cθ + vdec * sθ

    vx = @. vθ*cθ + vR*sθ*cϕ + rv*sθ*sϕ
    vy = @. -vθ*sθ + vR*cθ*cϕ + rv*cθ*sϕ
    vz = @. rv * cϕ - vR * sϕ

    vtot1 = @. sqrt(vx^2 + vy^2 + vz^2)
    vtot2 = @. sqrt(rv^2 + vra^2 + vdec^2)
    @assert all(isapprox.(vtot1, vtot2)) # make sure length is preserved

    return vx, vy, vz
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
    rv_gsr_shift(ra, dec)

Compute the apparent radial velocity of a stationary object at ra0, dec0
due to the solar motion.
"""
function rv_gsr_shift(ra0::Real, dec0::Real)
    rv_offset = lguys.transform(lguys.ICRS, lguys.GSR(ra=ra0, dec=dec0, radial_velocity=0)).radial_velocity
end


"""
    L_RV_BKD(rv::Real, ra0, dec0)

Compute the likelihood of the velocity assuming 
a halo velocity dispersion of 100 km/s and mean velocity
stationary.
"""
function L_RV_BKD(rv::Real, ra0::Real, dec0::Real)
    μ = rv_gsr_shift(ra0, dec0)
    σ = 100
    return pdf(Normal(μ, σ), rv)
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
