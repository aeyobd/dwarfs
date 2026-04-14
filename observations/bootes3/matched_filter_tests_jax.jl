### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 4c482a6f-9cd3-4704-b8ce-b23e0f1a293f
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya, CairoMakie
end

# ╔═╡ 34c475e5-1eb9-482a-a29b-182caad50313
using DataFrames

# ╔═╡ e3cb612c-7c1f-4688-8aa0-3cc4de19ec9b
using PyFITS

# ╔═╡ 8fa62036-f19f-4f65-8f14-f0b1fb2ddc27
md"""
The goal of this notebook is to reproduce the satellite density maps in Jensen et al. 2024 and their thesis.

"""

# ╔═╡ 015eae5d-515c-4fa4-9a44-fc478692cd44
CairoMakie.activate!(type=:png)

# ╔═╡ 68cd2e80-a8d7-4fcd-87dd-a632818b77ac
import KernelDensity

# ╔═╡ e560b3a4-249d-11f1-a7d8-918587896934
module Utils
	include("delve_utils.jl")
end

# ╔═╡ d045eaad-a22c-4def-8000-0f767b72ffed
module GaiaUtils
	include("gaia_utils.jl")
end

# ╔═╡ 5246991c-805c-420a-9f38-0569955c8fb8
md"""
# utilitis
"""

# ╔═╡ 6dabd441-9f39-4e75-a651-7edba2da41e4
"""
    build_cmd_filter(
        iso_color, iso_mag, iso_mass,
        color_uncertainty,
        distance_uncertainty;
        color_range  = (-0.5, 2.5),
        mag_range    = (15.0, 28.0),
        n_color      = 200,
        n_mag        = 300,
        imf          = kroupa_imf,
    ) -> (density, color_edges, mag_edges)

Build a matched-filter CMD density map from an isochrone, weighting
each point by the IMF and convolving along both the colour axis (photometric
uncertainty) and the magnitude axis (distance modulus uncertainty).

# Arguments
- `iso_color`            : Vector of isochrone colours (e.g. g-r)
- `iso_mag`              : Vector of isochrone magnitudes (same length)
- `iso_mass`             : Vector of stellar masses in solar units (same length)
- `color_uncertainty`    : Callable σ_c(mag) -> Float64, 1-σ colour error
- `distance_uncertainty` : Callable σ_μ(mag) -> Float64, 1-σ distance modulus
                           uncertainty in magnitudes. For a fixed σ_μ, pass
                           `_ -> σ_μ_value`. For a depth-dependent spread
                           (e.g. HB stars at different distances along the
                           line of sight), make this mag-dependent.

# Returns
- `density`     : Matrix{Float64} (n_color × n_mag), normalised to sum=1
- `color_edges` : Vector of colour bin edges (length n_color + 1)
- `mag_edges`   : Vector of magnitude bin edges (length n_mag + 1)
"""
function build_cmd_filter(
    iso_color::AbstractVector{<:Real},
    iso_mag::AbstractVector{<:Real},
    iso_mass::AbstractVector{<:Real},
    color_uncertainty,           # callable: σ_c(mag)   -> Float64
    distance_uncertainty;        # callable: σ_μ(mag)   -> Float64
    color_range::Tuple{Real,Real} = color_range,
    mag_range::Tuple{Real,Real}   = (14.0, 24.0),
    n_color::Int = 200,
    n_mag::Int   = 200,
    imf = Utils.kroupa_imf,
	mass_weight = true
	
)
    @assert length(iso_color) == length(iso_mag) == length(iso_mass) > 0

    # --- 1. Build output grid ---
    color_edges  = range(color_range[1], color_range[2]; length = n_color + 1)
    mag_edges    = range(mag_range[1],   mag_range[2];   length = n_mag + 1)
    Δc = step(color_edges)
    Δm = step(mag_edges)

    color_centres = collect(color_edges[1:end-1] .+ Δc / 2)
    mag_centres   = collect(mag_edges[1:end-1]   .+ Δm / 2)

    density = zeros(Float64, n_color, n_mag)

    # --- 2. Loop over isochrone points ---
    for k in eachindex(iso_color)
        c0 = iso_color[k]
        m0 = iso_mag[k]
        ms = iso_mass[k]
		if (m0 > mag_range[2]) || (m0 < mag_range[1])
			continue
		end

        # IMF weight × mass interval
        Δmass = if k == 1
            abs(iso_mass[2]   - iso_mass[1])   / 2
        elseif k == lastindex(iso_mass)
            abs(iso_mass[end] - iso_mass[end-1]) / 2
        else
            abs(iso_mass[k+1] - iso_mass[k-1]) / 4
        end
		if !mass_weight
			Δmass = 1
		end
		
        w = imf(ms) * Δmass
		if k % 4000 == 0
			println(w)
		end
        iszero(w) && continue

        # Uncertainty kernels at this point
        σ_c = max(color_uncertainty(m0),    Δc)   # colour σ, floored to 1 pixel
        σ_μ = max(distance_uncertainty(m0), Δm)   # distance modulus σ, floored

        # Precompute 1D kernel values to avoid redundant exp() calls
        # — colour kernel (length n_color)
        inv2σc² = 1.0 / (2σ_c^2)
        norm_c  = 1.0 / (σ_c * sqrt(2π))
        kc = [norm_c * exp(-(color_centres[i] - c0)^2 * inv2σc²) for i in 1:n_color]

        # — magnitude kernel (length n_mag); only evaluate bins within ±4σ
        inv2σm² = 1.0 / (2σ_μ^2)
        norm_m  = 1.0 / (σ_μ * sqrt(2π))
        j_lo = max(1,    searchsortedlast(mag_edges, m0 - 4σ_μ))
        j_hi = min(n_mag, searchsortedlast(mag_edges, m0 + 4σ_μ))

        # --- 3. Outer product: colour kernel ⊗ magnitude kernel ---
        for j in j_lo:j_hi
            km = norm_m * exp(-(mag_centres[j] - m0)^2 * inv2σm²)
            wkm = w * km
            iszero(wkm) && continue
            for i in 1:n_color
                density[i, j] += wkm * kc[i]
            end
        end
    end

    # --- 4. Normalise ---
    s = sum(density)
    s > 0 && (density ./= s)

    return density, collect(color_edges), collect(mag_edges)
end

# ╔═╡ 7eab2c84-6335-475b-be25-e8c2fce07776
let
	m = LinRange(0.01, 4, 1000)
	lines(log10.(m), log10.(Utils.kroupa_imf.(m)),
		 axis=(; 
			  xlabel="mass",
			  ylabel="IMF",
			 ))
end

# ╔═╡ b2befe17-f705-44f2-bd2d-bddc3874b115
md"""
# Gaia: Sculptor
"""

# ╔═╡ 7f3a5c3f-fe9d-4a10-ae7d-b70bbdfccfcc
dm_scl = 19.67

# ╔═╡ 9d3d4d9f-72ef-48a8-8e95-41059f67e758
dm_err_scl = 0.14

# ╔═╡ 0a8c2527-6812-4293-b60e-8d91bed66712
iso_gaia = let
	df = GaiaUtils.get_isochrone(-1.68, 12)
	df[df.label .< 4, :]
end

# ╔═╡ c9b72c30-a6a8-45d0-8b89-86174370ccce
stars_gaia = let
	df = read_fits("data/j24_1c.fits")

	df[df.F_BEST .== 1, :]
end

# ╔═╡ df5e5629-8d6b-4cb7-92bf-bf4467b6c2f2
err_gaia = Utils.fit_color_err(stars_gaia.phot_g_mean_mag, stars_gaia.dBP)

# ╔═╡ 836339de-1267-4b14-aef9-56c37a62b0de
cmd_gaia = build_cmd_filter(
	iso_gaia.G_BPftmag .- iso_gaia.G_RPmag ,
	iso_gaia.Gmag .+ dm_scl,
	iso_gaia.Mini,
	x -> @.(sqrt(err_gaia.(x)^2 + 0.1^2)),
	x -> dm_err_scl,
	mag_range=(16, 28),
	color_range=(-0.5, 2.5)
)

# ╔═╡ a7d83f97-159b-404f-b2f2-b595a9fe9ee7
cmd_gaia2 = build_cmd_filter(
	iso_gaia.G_BPftmag .- iso_gaia.G_RPmag ,
	iso_gaia.Gmag .+ dm_scl,
	iso_gaia.Mini,
	x -> @.(sqrt(err_gaia.(x)^2 + 0.1^2)),
	x -> dm_err_scl,
	mag_range=(16, 28),
	color_range=(-0.5, 2.5),
	mass_weight=false,
)

# ╔═╡ a79f6e5d-f8bb-4710-b7a4-4d33326ef15e
md"""
# Delve
"""

# ╔═╡ e8c0e80d-336e-4f23-b826-278d90d732d1
import TOML

# ╔═╡ 8369d4f8-922b-4a64-8e56-d309e1ca491b
obs_props_scl = TOML.parsefile("../sculptor/observed_properties.toml")

# ╔═╡ be0654ee-4c84-49e5-a459-ce1f36a1f359
obs_props_scl["metallicity"]

# ╔═╡ 6ed9796f-b522-4e62-82b5-62f0297c1e35
obs_props = TOML.parsefile("observed_properties.toml")

# ╔═╡ 6c55e84e-39a9-444b-a067-95a3e3826cbd
distance_modulus = obs_props["distance_modulus"]

# ╔═╡ feceea26-7bc7-4789-9f4b-4df4c5c4c27a
iso_gaia

# ╔═╡ ef15705f-2390-40d5-b6d7-50064ca333ec
LilGuys.kpc2dm(LilGuys.Measurement(26.56, 0.25))

# ╔═╡ 53a1dd27-bcd8-481b-abd7-160da5dc643d
function resample_iso(iso; N=1000_000, cols=[:gmag, :rmag])
	massrange = extrema(iso.Mini)

	interps = Dict(col => LilGuys.lerp(iso.Mini, iso[!, col]) for col in cols)
	masses = LinRange(massrange..., N)
	df =  DataFrame(
		:Mini => masses,
	)

	for col in cols
		df[:, col] = interps[col].(masses)
	end
	df
end

# ╔═╡ 52792fb8-d9cf-4b2c-81f8-e3d3940c2250
iso_gaia_resampled = resample_iso(iso_gaia, cols=[:Gmag, :G_BPftmag, :G_RPmag]);

# ╔═╡ 777fb53b-a9d7-4218-9026-401c370bb415
cmd_gaia5 = build_cmd_filter(
	iso_gaia_resampled.G_BPftmag .- iso_gaia_resampled.G_RPmag ,
	iso_gaia_resampled.Gmag .+ dm_scl,
	iso_gaia_resampled.Mini,
	x -> @.(sqrt(err_gaia.(x)^2 + 0.1^2)),
	x -> dm_err_scl,
	mag_range=(16, 28),
	color_range=(-0.5, 2.5),

)

# ╔═╡ e3d35c58-3fad-4603-b52f-64d727c6b3ad
(diff(iso_gaia_resampled.Mini))

# ╔═╡ 53fc7954-4d68-4c78-bf8c-75ad7577b0e6
cmd_gaia3 = build_cmd_filter(
	iso_gaia_resampled.G_BPftmag .- iso_gaia_resampled.G_RPmag ,
	iso_gaia_resampled.Gmag .+ dm_scl,
	iso_gaia_resampled.Mini,
	x -> @.(sqrt(err_gaia.(x)^2 + 0.1^2)),
	x -> dm_err_scl,
	mag_range=(16, 28),
	color_range=(-0.5, 2.5),
	imf = x -> 1,
)

# ╔═╡ 6d62f006-1311-48f9-9fed-a8438fee222f
cmd_gaia4 = build_cmd_filter(
	iso_gaia_resampled.G_BPftmag .- iso_gaia_resampled.G_RPmag ,
	iso_gaia_resampled.Gmag .+ dm_scl,
	iso_gaia_resampled.Mini,
	x -> @.(sqrt(err_gaia.(x)^2 + 0.1^2)),
	x -> dm_err_scl,
	mag_range=(16, 28),
	color_range=(-0.5, 2.5),
	mass_weight=false,
	imf = x->1,

)

# ╔═╡ f63389bf-77fe-4495-9be5-b4f86f84dc57
stars = let
	df = read_fits("data/delve_dr2_good.fits")
	df[!, :bp_rp] = df.gmag .- df.rmag
	df[!, :G] = df.gmag

	df[!, :xi], df[!, :eta] = 60 .* LilGuys.to_tangent(df.ra, df.dec, obs_props["ra"], obs_props["dec"])
	df[!, :R] = @. sqrt(df.xi^2 + df.eta^2)
	df
end

# ╔═╡ baccd40f-4128-4031-973b-1389845f2972
stars_cen = stars[stars.R .< 30, :]

# ╔═╡ 4e7dacb3-366e-46f4-84f4-fd6206961c4e
let
	fig = Figure()
	ax = Axis(fig[1,1])
		
	hexbin!(stars.xi, stars.eta, bins=300)

	arc!((0,0), 30, 0, 2π)

	arc!((0,0), 5*obs_props["R_h"], 0, 2π)
	arc!((0,0), 7*obs_props["R_h"], 0, 2π)
	fig
end

# ╔═╡ 90a3ba08-28eb-4529-bfb7-c6d01b3dd979
R_h = obs_props["R_h"]

# ╔═╡ 22934ab9-2750-4826-b916-a51ecee71f68
stars_bkg = stars[stars.R .> 5R_h, :]

# ╔═╡ 88a5e8a3-00dd-47b3-ba2f-ab3ffc483c25
hexbin(stars_bkg.xi, stars_bkg.eta, bins=300)

# ╔═╡ c94244bd-31f0-4838-8495-05fc8be234e8
hist(stars.gmag)

# ╔═╡ bb0600b7-cf8a-484b-9d4a-7dab2c963da6
err_model(x, params) = @. (params[1] + x*params[2] + x^2 * params[3])

# ╔═╡ 764be90a-e602-44d8-9956-d407ee2a0d4e
gr_err = @. sqrt(stars.magerr_psf_g^2 + stars.magerr_psf_r^2)

# ╔═╡ 39d177ca-6402-4216-b503-c12f83f91d56
popt, covt = LilGuys.curve_fit(err_model, stars.gmag, log10.(gr_err), [0., 1., 1.])

# ╔═╡ 48993623-e10e-45c7-97df-93e3e0200007
color_err = Utils.fit_color_err(stars)

# ╔═╡ 3b521894-468e-4cad-bb0f-394567adcbdf
let
	fig = Figure()
	ax = Axis(fig[1,1])

	scatter!(stars.gmag, color_err.(stars.gmag))
	fig
end

# ╔═╡ bf0e16fe-94cd-41fd-a107-4a814dd448f4
color_err_gaia = Utils.fit_color_err(stars)

# ╔═╡ 080d5e15-dd13-423a-815c-e22edfc73dcc
iso = let 
	iso = Utils.get_isochrone(-2.1, 12)
	iso[iso.label .< 4, :]
end

# ╔═╡ 78b5414f-8961-4446-bfb7-9c0e852dd5b7
iso_new = resample_iso(iso)

# ╔═╡ eaa3ee0b-6b8f-4da5-a230-596613bbfa76
lines(iso_new.gmag .- iso_new.rmag, iso_new.gmag .+ distance_modulus)

# ╔═╡ 740afb89-f49d-4079-ac9e-f77a872bfc2c
lines(midpoints(iso.Mini), log10.(diff(iso.Mini)))

# ╔═╡ 71dfe223-539a-4e2b-84b1-04811b2533a0
color_range = (-0.5, 2.5)

# ╔═╡ 13896043-e020-409f-99fb-8d1254b994d2


# ╔═╡ 7cf639a4-cd5e-4f42-a0ac-a7cd768debdc
cmd_sat2 = build_cmd_filter(
	iso.gmag .- iso.rmag ,
	iso.gmag .+ distance_modulus,
	iso.Mini,
	x -> @.(sqrt(color_err.(x)^2 + 0.1^2)),
	x -> 0.02,
	n_color=80,
	n_mag=80,
	mass_weight=false,
	mag_range=(16, 24)
)

# ╔═╡ 109fdcda-b64e-477f-bfba-4f9588b7bf09
cmd_sat = build_cmd_filter(
	iso_new.gmag .- iso_new.rmag ,
	iso_new.gmag .+ distance_modulus,
	iso_new.Mini,
	x -> @.(sqrt(color_err.(x)^2 + 0.1^2)),
	x -> 0.05
)

# ╔═╡ a75c32b1-9011-4c88-8725-624f93fd0d98
let
	fig = Figure()

	ax = Axis(fig[1,1],
			 yreversed=true)
	y = log10.(cmd_sat[1] ./ maximum(cmd_sat[1]))
	
	p = heatmap!(cmd_sat[2], cmd_sat[3], y,colorrange=(-5, 0))


	scatter!(stars_cen.gmag .- stars_cen.rmag, stars_cen.gmag, markersize=2, color=COLORS[3])
	fig
end

# ╔═╡ 9ed01175-a59c-4b27-8b6c-5d286ea6faad
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 yreversed=true)
	y = (cmd_sat[1] ./ maximum(cmd_sat[1]))
	
	p = heatmap!(cmd_sat[2], cmd_sat[3], y,colorrange=(0, 1))

	Colorbar(fig[1,2], p, )
	fig
end

# ╔═╡ 16d32bcb-6f43-499b-9c88-94db8464e63e
cmd_sat_noimf = build_cmd_filter(
	iso_new.gmag .- iso_new.rmag ,
	iso_new.gmag .+ distance_modulus,
	iso_new.Mini,
	x -> @.(sqrt(color_err.(x)^2 + 0.1^2)),
	x -> 0.05,
	imf = x->1,
)

# ╔═╡ 79bb12ce-2a31-47c6-9265-2c5ec9dffa48
function cmd_axis(gs)
	return Axis(gs,
			   xlabel = "g - r",
			   ylabel = "g",
				yreversed = true,
			   )

end
	

# ╔═╡ e1eb418c-2da7-4a9a-a118-a04d10f1aa79
function plot_cmd(cmd; kwargs...)	
	fig = Figure()
	ax = cmd_axis(fig[1,1])
	
	y = (cmd[1] ./ maximum(cmd[1]))
	
	p = heatmap!(ax, cmd[2], cmd[3], y; kwargs...)

	Colorbar(fig[1,2], p, label="density")
	fig
end


# ╔═╡ 5fec3242-c0e2-480f-a941-d2540c7de79a
plot_cmd(cmd_gaia, colorrange=(1e-4, 1),  colorscale=log10)

# ╔═╡ a4159d32-8fd0-4577-ab14-8841ae328b52
plot_cmd(cmd_sat)

# ╔═╡ 8fc6ad6e-5cdd-48b6-94c4-97933d4d2d69
plot_cmd(cmd_sat2)

# ╔═╡ 2372feb6-7302-4297-bb7d-7c35935d0d43
plot_cmd(cmd_sat_noimf)

# ╔═╡ 75a23c31-7e10-4fe4-bcbd-1deb79a5211b
plot_cmd(cmd_sat_noimf)

# ╔═╡ 2d75bc5e-525a-47ad-8b2f-e6e3d4cf7024
plot_cmd(cmd_sat_noimf)

# ╔═╡ 74cefad0-c210-4db4-bd99-54308cedcf3a
function plot_log_cmd(cmd; colorrange=(-4, 0), kwargs...)	
	fig = Figure()
	ax = cmd_axis(fig[1,1])
	
	y = log10.(cmd[1] ./ maximum(cmd[1]))
	
	p = heatmap!(ax, cmd[2], cmd[3], y, colorrange=colorrange; kwargs...)

	Colorbar(fig[1,2], p, label="density")
	fig
end


# ╔═╡ 62a063a7-9e08-45cd-b187-b2c9cf587197
plot_log_cmd(cmd_gaia,)

# ╔═╡ a008c89c-3fae-4055-9d4f-1bfab2f79410
plot_log_cmd(cmd_gaia2)

# ╔═╡ ce27b5a3-372f-482e-b80e-b8b5fb4738bd
plot_log_cmd(cmd_gaia3,)

# ╔═╡ 7158c2e8-1af4-4bcf-bcd6-34222d3b9462
plot_log_cmd(cmd_gaia4,)

# ╔═╡ aeb7f058-02a7-494b-b24d-07dc646ed200
function plot_cmd_rel(cmd, cmd2; kwargs...)	
	fig = Figure()
	ax = cmd_axis(fig[1,1])
	
	y = (cmd[1] ./ cmd2[1])
	
	p = heatmap!(ax, cmd[2], cmd[3], y; kwargs...)

	Colorbar(fig[1,2], p, label="relative density")
	fig
end


# ╔═╡ 4b955d32-6e61-4d0d-a819-1eb307114015
plot_cmd_rel(cmd_gaia3, cmd_gaia, colorscale=log10, colorrange=(0.1, 10))

# ╔═╡ 65973035-10d5-4ece-8b21-c51de60daf90
let
	fig = Figure(size=(6, 2.5,) .* 72)
	ax = Axis(fig[1,1],
			 # limits=((15, 25), (0, 0.05)),
			  xlabel = "gmag",
			  ylabel = "density",
			  yscale = log10,
			  limits=(16, 23, 0.001, 1.0),
			  yticks = [0.01, 0.1, 1]
			 )
	x = midpoints(cmd_sat[3])
	y = sum(cmd_sat[1], dims=1) |> vec

	S = sum(diff(cmd_sat[3]) .* y)
	lines!(x, y ./ S, label="binned density map")

	m = LinRange(0.75, maximum(iso.Mini), 10000)
	w = Utils.kroupa_imf.(m)

	f = LilGuys.lerp(iso.Mini, iso.gmag .+ distance_modulus)

	stephist!(f.(m), bins=300, weights=w, normalization=:pdf, color=COLORS[2], label = "isochrone mag(kroupa IMF)")


	axislegend( ax, position=:rb)
	
	fig
end

# ╔═╡ 6288ddaf-517e-484a-ac6b-1838e828449e
md"""
# Testing against NGC 5466
"""

# ╔═╡ 44cab61f-1dcb-4ff4-8661-b8385ed0239e


# ╔═╡ 02224588-b6f9-4939-9caa-0da3920c59c5
obs_props_ngc5466 = Dict(
	"ra" => 211.36371,
	"dec" => 28.53444,
	"pmra" => -5.41,
	"pmdec" => -0.79,
	"distance_modulus" => LilGuys.kpc2dm(16.9),
	"R_h" => 10.72,
	"ellipticity" => 0,
	"position_angle" => 0,
	"metallicity" => -2.0,
)

# ╔═╡ 56f74781-7f4c-4eac-b6a6-729826a2a311
xi_ngc5466, eta_ngc5466 = 60 .* LilGuys.to_tangent(obs_props_ngc5466["ra"], obs_props_ngc5466["dec"], obs_props["ra"], obs_props["dec"])

# ╔═╡ a7b93b84-72e5-4f43-80d6-f5177cacc381
sample_ngc5466 = let
	r_ngc = @. sqrt(
		(stars.xi - xi_ngc5466)^2 
		+ (stars.eta - eta_ngc5466) ^ 2
	)

	df = stars[r_ngc .< obs_props_ngc5466["R_h"] * 1, :]

	df.xi .-= xi_ngc5466
	df.eta .-= eta_ngc5466

	df
end

# ╔═╡ 2a3f1432-0390-49b6-a539-e0b6c907da76
let
	fig = Figure()
	ax = cmd_axis(fig[1,1], )
	ax.title[] = "NGC 5466"
	
	scatter!(sample_ngc5466.gmag .- sample_ngc5466.rmag, sample_ngc5466.gmag, markersize=1, color=:black, alpha=0.1)

	fig

end

# ╔═╡ 8ded2ee2-e4c2-406d-9f88-d0bb1592e663
iso_ngc5466 = let 
	iso = Utils.get_isochrone(-1.9, 12)
	iso[iso.label .< 5, :] |> resample_iso
end

# ╔═╡ a60d65fd-da72-4dc8-bc7b-133d0bb7a51e


# ╔═╡ cf2317a0-0a82-4957-a6cb-d5aff70cf53a
cmd_ngc5466 = build_cmd_filter(
	iso_ngc5466.gmag .- iso_ngc5466.rmag ,
	iso_ngc5466.gmag .+ obs_props_ngc5466["distance_modulus"],
	iso_ngc5466.Mini,
	x -> @.(sqrt(color_err.(x)^2 + 0.03^2)),
	x -> 0.05,
	mag_range=(14, 21)
)

# ╔═╡ 8d52e45a-11bc-4334-bc51-2d9c7cc68e56


# ╔═╡ 05118762-a483-4fed-ad77-6789a8a1d494
function make_background_density(stars;    
	color_range::Tuple{Real,Real} = color_range,
    mag_range::Tuple{Real,Real}   = (14.0, 21.0),
    n_color::Int = 200,
	n_mag::Int = 200,
								 kwargs...
								 
)

	KernelDensity.kde((stars.gmag .- stars.rmag, stars.gmag), boundary=(color_range, mag_range), npoints=(n_color, n_mag),)
end

# ╔═╡ 4466fab1-8bc1-4b2f-b2e0-cb1cd6bb540f
cmd_emperical_ngc5466 = make_background_density(sample_ngc5466, bandwidth=(0.025, 0.05))

# ╔═╡ b89a98c7-8dd8-4ae1-82a1-a3501f423eca
let
	fig = Figure(
		size = (6, 3) .* 72
	)
	
	ax = cmd_axis(fig[1,1])
	ax.title[] = "synthetic NGC 5466"
	cmd =cmd_ngc5466
	
	y = (cmd[1] ./ maximum(cmd[1]))
	
	p = heatmap!(ax, cmd[2], cmd[3], y; )


	ax2 = cmd_axis(fig[1,2])
	ax2.title[] = "emperical NGC 5466"

	heatmap!(cmd_emperical_ngc5466, colorrange=(0, maximum(cmd_emperical_ngc5466.density)))

	Colorbar(fig[1,3], p, label="density (peak normalized)")


	linkaxes!(ax, ax2)
	xlims!(color_range...)
	ylims!(14, 21)
	ax2.yreversed[] = true
	fig

end

# ╔═╡ 32d10087-2acc-4abb-be53-7e7ce03f8de0
let
	fig = Figure(
		size = (6, 3) .* 72
	)
	
	ax = cmd_axis(fig[1,1])
	ax.title[] = "synthetic NGC 5466"
	cmd =cmd_ngc5466
	
	y = (cmd[1] ./ maximum(cmd[1]))
	
	p = heatmap!(ax, cmd[2], cmd[3], y; colorrange=(1e-5, 1), colorscale=log10)


	ax2 = cmd_axis(fig[1,2])
	ax2.title[] = "emperical NGC 5466"

	cmax = maximum(cmd_emperical_ngc5466.density)
	heatmap!(cmd_emperical_ngc5466, colorrange=(1e-5*cmax, cmax), colorscale=log10)

	Colorbar(fig[1,3], p, label="density (peak normalized)")


	linkaxes!(ax, ax2)
	xlims!(color_range...)
	ylims!(14, 21)
	ax2.yreversed[] = true
	fig

end

# ╔═╡ 937a539d-f664-488a-8d6f-9d6fde9df0eb
let
	fig = Figure(size=(6, 2.5,) .* 72)
	ax = Axis(fig[1,1],
			  xlabel = "gmag",
			  ylabel = "density",
			  yscale = log10,
			  limits=(12, 23, 0.001, 1.0),
			  yticks = [0.01, 0.1, 1]
			 )
	x = midpoints(cmd_ngc5466[3])
	y = sum(cmd_ngc5466[1], dims=1) |> vec

	S = sum(diff(cmd_ngc5466[3]) .* y)
	lines!(x, y ./ S, label="binned density map")



	x = cmd_emperical_ngc5466.y
	y = sum(cmd_emperical_ngc5466.density, dims=1) |> vec

	S = sum(LilGuys.gradient(cmd_emperical_ngc5466.x) .* y)
	lines!(x, y ./ S, label="binned density ")

	println(LilGuys.mean(y ./ S))
	println(LilGuys.mean(x))
	
	m = LinRange(0.75, maximum(iso.Mini), 10000)
	w = Utils.kroupa_imf.(m)

	f = LilGuys.lerp(iso.Mini, iso.gmag .+ obs_props_ngc5466["distance_modulus"])

	stephist!(f.(m), bins=300, weights=w, normalization=:pdf, color=COLORS[2], label = "isochrone mag(kroupa IMF)")


	axislegend( ax, position=:rb)
	
	fig
end

# ╔═╡ Cell order:
# ╟─8fa62036-f19f-4f65-8f14-f0b1fb2ddc27
# ╠═4c482a6f-9cd3-4704-b8ce-b23e0f1a293f
# ╠═015eae5d-515c-4fa4-9a44-fc478692cd44
# ╠═68cd2e80-a8d7-4fcd-87dd-a632818b77ac
# ╠═34c475e5-1eb9-482a-a29b-182caad50313
# ╠═e3cb612c-7c1f-4688-8aa0-3cc4de19ec9b
# ╠═e560b3a4-249d-11f1-a7d8-918587896934
# ╠═d045eaad-a22c-4def-8000-0f767b72ffed
# ╟─5246991c-805c-420a-9f38-0569955c8fb8
# ╠═6dabd441-9f39-4e75-a651-7edba2da41e4
# ╠═7eab2c84-6335-475b-be25-e8c2fce07776
# ╠═e1eb418c-2da7-4a9a-a118-a04d10f1aa79
# ╠═74cefad0-c210-4db4-bd99-54308cedcf3a
# ╠═aeb7f058-02a7-494b-b24d-07dc646ed200
# ╟─b2befe17-f705-44f2-bd2d-bddc3874b115
# ╠═8369d4f8-922b-4a64-8e56-d309e1ca491b
# ╠═7f3a5c3f-fe9d-4a10-ae7d-b70bbdfccfcc
# ╠═9d3d4d9f-72ef-48a8-8e95-41059f67e758
# ╠═be0654ee-4c84-49e5-a459-ce1f36a1f359
# ╠═0a8c2527-6812-4293-b60e-8d91bed66712
# ╠═52792fb8-d9cf-4b2c-81f8-e3d3940c2250
# ╠═c9b72c30-a6a8-45d0-8b89-86174370ccce
# ╠═df5e5629-8d6b-4cb7-92bf-bf4467b6c2f2
# ╠═836339de-1267-4b14-aef9-56c37a62b0de
# ╠═777fb53b-a9d7-4218-9026-401c370bb415
# ╠═e3d35c58-3fad-4603-b52f-64d727c6b3ad
# ╠═a7d83f97-159b-404f-b2f2-b595a9fe9ee7
# ╠═53fc7954-4d68-4c78-bf8c-75ad7577b0e6
# ╠═6d62f006-1311-48f9-9fed-a8438fee222f
# ╠═62a063a7-9e08-45cd-b187-b2c9cf587197
# ╠═a008c89c-3fae-4055-9d4f-1bfab2f79410
# ╠═ce27b5a3-372f-482e-b80e-b8b5fb4738bd
# ╠═7158c2e8-1af4-4bcf-bcd6-34222d3b9462
# ╠═4b955d32-6e61-4d0d-a819-1eb307114015
# ╠═5fec3242-c0e2-480f-a941-d2540c7de79a
# ╠═a79f6e5d-f8bb-4710-b7a4-4d33326ef15e
# ╠═e8c0e80d-336e-4f23-b826-278d90d732d1
# ╠═6ed9796f-b522-4e62-82b5-62f0297c1e35
# ╠═6c55e84e-39a9-444b-a067-95a3e3826cbd
# ╠═feceea26-7bc7-4789-9f4b-4df4c5c4c27a
# ╠═ef15705f-2390-40d5-b6d7-50064ca333ec
# ╠═3b521894-468e-4cad-bb0f-394567adcbdf
# ╠═eaa3ee0b-6b8f-4da5-a230-596613bbfa76
# ╠═78b5414f-8961-4446-bfb7-9c0e852dd5b7
# ╠═53a1dd27-bcd8-481b-abd7-160da5dc643d
# ╠═f63389bf-77fe-4495-9be5-b4f86f84dc57
# ╠═740afb89-f49d-4079-ac9e-f77a872bfc2c
# ╠═baccd40f-4128-4031-973b-1389845f2972
# ╠═4e7dacb3-366e-46f4-84f4-fd6206961c4e
# ╠═88a5e8a3-00dd-47b3-ba2f-ab3ffc483c25
# ╠═22934ab9-2750-4826-b916-a51ecee71f68
# ╠═90a3ba08-28eb-4529-bfb7-c6d01b3dd979
# ╠═c94244bd-31f0-4838-8495-05fc8be234e8
# ╠═a75c32b1-9011-4c88-8725-624f93fd0d98
# ╠═a4159d32-8fd0-4577-ab14-8841ae328b52
# ╠═8fc6ad6e-5cdd-48b6-94c4-97933d4d2d69
# ╠═2372feb6-7302-4297-bb7d-7c35935d0d43
# ╠═9ed01175-a59c-4b27-8b6c-5d286ea6faad
# ╠═75a23c31-7e10-4fe4-bcbd-1deb79a5211b
# ╠═2d75bc5e-525a-47ad-8b2f-e6e3d4cf7024
# ╠═bb0600b7-cf8a-484b-9d4a-7dab2c963da6
# ╠═764be90a-e602-44d8-9956-d407ee2a0d4e
# ╠═39d177ca-6402-4216-b503-c12f83f91d56
# ╠═48993623-e10e-45c7-97df-93e3e0200007
# ╠═bf0e16fe-94cd-41fd-a107-4a814dd448f4
# ╠═080d5e15-dd13-423a-815c-e22edfc73dcc
# ╠═71dfe223-539a-4e2b-84b1-04811b2533a0
# ╠═13896043-e020-409f-99fb-8d1254b994d2
# ╠═7cf639a4-cd5e-4f42-a0ac-a7cd768debdc
# ╠═109fdcda-b64e-477f-bfba-4f9588b7bf09
# ╠═16d32bcb-6f43-499b-9c88-94db8464e63e
# ╠═79bb12ce-2a31-47c6-9265-2c5ec9dffa48
# ╠═65973035-10d5-4ece-8b21-c51de60daf90
# ╠═6288ddaf-517e-484a-ac6b-1838e828449e
# ╠═44cab61f-1dcb-4ff4-8661-b8385ed0239e
# ╠═02224588-b6f9-4939-9caa-0da3920c59c5
# ╠═56f74781-7f4c-4eac-b6a6-729826a2a311
# ╠═a7b93b84-72e5-4f43-80d6-f5177cacc381
# ╠═2a3f1432-0390-49b6-a539-e0b6c907da76
# ╠═8ded2ee2-e4c2-406d-9f88-d0bb1592e663
# ╠═a60d65fd-da72-4dc8-bc7b-133d0bb7a51e
# ╠═cf2317a0-0a82-4957-a6cb-d5aff70cf53a
# ╠═8d52e45a-11bc-4334-bc51-2d9c7cc68e56
# ╠═05118762-a483-4fed-ad77-6789a8a1d494
# ╠═4466fab1-8bc1-4b2f-b2e0-cb1cd6bb540f
# ╠═b89a98c7-8dd8-4ae1-82a1-a3501f423eca
# ╠═32d10087-2acc-4abb-be53-7e7ce03f8de0
# ╠═937a539d-f664-488a-8d6f-9d6fde9df0eb
