using LilGuys
using Distributions
import DataFrames: DataFrame
import TOML
import StatsBase: sem

import DensityEstimators as DE


module GaiaFilters
    include("../../utils/gaia_filters.jl")
end

"""
Structural parameters to sample over
"""
Base.@kwdef struct StructuralParams
    position_err::Float64

    ellipticity::Float64
    ellipticity_err::Float64

    position_angle::Float64
    position_angle_err::Float64

    bins::Array{Float64}
end


Base.@kwdef struct GaiaData
    source_id::Vector{Float64}
    xi::Vector{Float64}
    eta::Vector{Float64}
    "Prior BG probability"
    P_bg::Vector{Float64}
    "Prior Satellite probability"
    P_sat::Vector{Float64}
end


function get_data_filename(galaxy)
	filenames = ["jensen+24_1c.fits", "jensen+24_2c.fits", "jensen+24_wide.fits", "j24_1c.fits", "j24_2c.fits"]
	filename = ""
    dir = "../$galaxy/data"
	for file in filenames
		if isfile(joinpath(dir, file))
			filename =  joinpath(dir, file)
		end
	end

	filename
end

function get_fits(galaxy, obs_props)
    filename = get_data_filename(galaxy)
    @info "reading fits from $filename"
    params = GaiaFilters.GaiaFilterParams(obs_props, filename=filename)
    all_stars =  GaiaFilters.read_gaia_stars(params)
    best = GaiaFilters.select_members(all_stars, params)

    best[!, :L_BKD_nospace] = best.L_CMD_BKD .* best.L_PM_BKD
    best[!, :L_SAT_nospace] = best.L_CMD_SAT .* best.L_PM_SAT
    return best
end


function get_obs_props(galaxy)
    return TOML.parsefile(joinpath("..", galaxy, "observed_properties.toml"))
end




function GaiaData(df::DataFrame)
    return GaiaData(
        df.source_id,
        df.xi,
        #df.xi_err,
        df.eta,
        #df.eta_err,
        df.L_BKD_nospace,
        df.L_SAT_nospace,
       )
end


function StructuralParams(data, obs_props; PSAT_min=0.99, kwargs...)
    cen_err = centring_error(data, PSAT_min=PSAT_min)
    bins = default_bins(data; kwargs...)
    return StructuralParams(cen_err, 
        obs_props["ellipticity"], read_error(obs_props, "ellipticity"),
        obs_props["position_angle"], read_error(obs_props, "position_angle"),
        bins
       )
end


"""
    read_error(obs_props, key)

Returns the maximum uncertainty of a key in the observed properties dataframe
"""
function read_error(obs_props, key)
	if key * "_em" ∈ keys(obs_props)
		return max(obs_props[key * "_em"], obs_props[key * "_ep"])
	else
		return obs_props[key * "_err"]
	end
end



"""
    centring_error(data; PSAT_min)

Given J+24 data, returns the centring error
"""
function centring_error(data; PSAT_min=0.99)
    stats = centring_stats(data, PSAT_min=PSAT_min)
    err_sys = max(abs(stats.eta_sys), abs(stats.xi_sys))
    err_stat = sqrt((stats.xi_stat^2 + stats.eta_stat^2)/2)
    return err_sys + err_stat
end


"""
    centring_stats(data; PSAT_min)

Given J+24 data, returns the centring statistics
"""
function centring_stats(data; PSAT_min=0.99)
    members = filter(r->r.PSAT>PSAT_min, data)


    return (;
        xi_stat = sem(members.xi),
        xi_sys = mean(members.xi),
        eta_stat = sem(members.eta),
        eta_sys = mean(members.eta),
        counts = size(members, 1),
       )
end




function get_R_max_cen(data, params::StructuralParams, n_sigma=3)
    R_max_circ = sqrt(maximum(@. data.xi^2 + data.eta^2))
    R_max_circ_min = sqrt(maximum(@. data.xi^2 + data.eta^2)) .- n_sigma*params.position_err
    ell_max = params.ellipticity + n_sigma * params.ellipticity_err

    return R_max_circ_min / (1 - ell_max)
end


function perturbed_radii(data, params::StructuralParams )

	d_xi = rand(Normal(0.0, params.position_err))
	d_eta = rand(Normal(0.0, params.position_err))
	pos_ang = rand(Normal(params.position_angle, params.position_angle_err))
	ell = rand(truncated(Normal(params.ellipticity, params.ellipticity_err), lower=0, upper=0.99))
	
    xi = data.xi .+ d_xi 
	eta = data.eta .+ d_eta

	radii = LilGuys.calc_R_ell(xi, eta, ell, pos_ang)

	return radii
end


function default_bins(stars; bin_width = 0.05, num_per_bin=nothing)
    if num_per_bin === nothing
        num_per_bin = max(round(Int, LilGuys.Interface.default_n_per_bin(stars.R_ell[stars.PSAT .> 0.2], nothing)), 2)
    end

    return 10 .^ LilGuys.bins_both(log10.(stars.R_ell), nothing,
        num_per_bin=num_per_bin, bin_width=bin_width)
end


function Σ_hist(radii::AbstractVector{<:Real}, bins::Vector{<:Real}, Σs::AbstractVector{<:Real})
    bin_idx = DE.bin_indices(radii, bins)
    Nbins = length(bins)

    Σs_safe = [0; Σs; 0] # 0 for values outsize bins
    Σ = Σs_safe[bin_idx .+ 1]
end
