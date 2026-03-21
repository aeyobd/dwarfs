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

# ╔═╡ 0a8c2527-6812-4293-b60e-8d91bed66712
iso_gaia = GaiaUtils.get_isochrone(-2, 12)

# ╔═╡ 7eab2c84-6335-475b-be25-e8c2fce07776
let
	m = LinRange(0.01, 4, 1000)
	lines(log10.(m), log10.(Utils.kroupa_imf.(m)),
		 axis=(; 
			  xlabel="mass",
			  ylabel="IMF",
			 ))
end

# ╔═╡ e8c0e80d-336e-4f23-b826-278d90d732d1
import TOML

# ╔═╡ 6ed9796f-b522-4e62-82b5-62f0297c1e35
obs_props = TOML.parsefile("observed_properties.toml")

# ╔═╡ 6c55e84e-39a9-444b-a067-95a3e3826cbd
distance_modulus = obs_props["distance_modulus"]

# ╔═╡ feceea26-7bc7-4789-9f4b-4df4c5c4c27a
iso_gaia

# ╔═╡ c9b72c30-a6a8-45d0-8b89-86174370ccce
stars_gaia = let
	df = read_fits("data/j24_1c.fits")

	df[df.F_BEST .== 1, :]
end

# ╔═╡ df5e5629-8d6b-4cb7-92bf-bf4467b6c2f2
err_gaia = Utils.fit_color_err(stars_gaia.phot_g_mean_mag, stars_gaia.dBP)

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
iso_gaia_resampled = resample_iso(iso_gaia, cols=[:Gmag, :G_BPftmag, :G_RPmag])

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
	# iso[iso.label .< 4, :]
end

# ╔═╡ 78b5414f-8961-4446-bfb7-9c0e852dd5b7
iso_new = resample_iso(iso)

# ╔═╡ eaa3ee0b-6b8f-4da5-a230-596613bbfa76
lines(iso_new.gmag .- iso_new.rmag, iso_new.gmag .+ distance_modulus)

# ╔═╡ 740afb89-f49d-4079-ac9e-f77a872bfc2c
lines(diff(iso.Mini))

# ╔═╡ 71dfe223-539a-4e2b-84b1-04811b2533a0
color_range = (-0.8, 1.8)

# ╔═╡ 13896043-e020-409f-99fb-8d1254b994d2


# ╔═╡ d222d357-8d82-494a-b85d-753314c044e4

function build_cmd_filter_noweight(
    iso_color::AbstractVector{<:Real},
    iso_mag::AbstractVector{<:Real},
    iso_mass::AbstractVector{<:Real},
    color_uncertainty,           # callable: σ_c(mag)   -> Float64
    distance_uncertainty;        # callable: σ_μ(mag)   -> Float64
    color_range::Tuple{Real,Real} = color_range,
    mag_range::Tuple{Real,Real}   = (15.0, 23.0),
    n_color::Int = 200,
    n_mag::Int   = 300,
    imf = Utils.kroupa_imf,
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
            abs(iso_mass[k+1] - iso_mass[k-1]) / 2
        end
        w = imf(ms)
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

# ╔═╡ a7d83f97-159b-404f-b2f2-b595a9fe9ee7
cmd_gaia2 = build_cmd_filter_noweight(
	iso_gaia_resampled.G_BPftmag .- iso_gaia_resampled.G_RPmag ,
	iso_gaia_resampled.Gmag .+ 19.60,
	iso_gaia_resampled.Mini,
	x -> @.(sqrt(err_gaia.(x)^2 + 0.1^2)),
	x -> 0.2,
	mag_range=(16, 21),
	color_range=(-0.5, 2.5)
)

# ╔═╡ 53fc7954-4d68-4c78-bf8c-75ad7577b0e6
cmd_gaia3 = build_cmd_filter_noweight(
	iso_gaia.G_BPftmag .- iso_gaia.G_RPmag ,
	iso_gaia.Gmag .+ 19.60,
	iso_gaia.Mini,
	x -> @.(sqrt(err_gaia.(x)^2 + 0.1^2)),
	x -> 0.2,
	mag_range=(16, 25),
	color_range=(-0.5, 2.5),
	n_color=200,
	n_mag=200
)

# ╔═╡ b8ae4ad0-afdf-4877-a02f-f9291d6411e9

function build_cmd_filter_neither(
    iso_color::AbstractVector{<:Real},
    iso_mag::AbstractVector{<:Real},
    iso_mass::AbstractVector{<:Real},
    color_uncertainty,           # callable: σ_c(mag)   -> Float64
    distance_uncertainty;        # callable: σ_μ(mag)   -> Float64
    color_range::Tuple{Real,Real} = color_range,
    mag_range::Tuple{Real,Real}   = (15.0, 23.0),
    n_color::Int = 200,
    n_mag::Int   = 300,
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
            abs(iso_mass[k+1] - iso_mass[k-1]) / 2
        end
        w = 1
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

# ╔═╡ 6dabd441-9f39-4e75-a651-7edba2da41e4

function build_cmd_filter_noiso(
    iso_color::AbstractVector{<:Real},
    iso_mag::AbstractVector{<:Real},
    iso_mass::AbstractVector{<:Real},
    color_uncertainty,           # callable: σ_c(mag)   -> Float64
    distance_uncertainty;        # callable: σ_μ(mag)   -> Float64
    color_range::Tuple{Real,Real} = color_range,
    mag_range::Tuple{Real,Real}   = (15.0, 23.0),
    n_color::Int = 200,
    n_mag::Int   = 300,
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
            abs(iso_mass[k+1] - iso_mass[k-1]) / 2
        end
        w = Δmass
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

# ╔═╡ cfc6be95-7c53-479f-81fc-ce4d9ca07b41
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
    mag_range::Tuple{Real,Real}   = (15.0, 23.0),
    n_color::Int = 200,
    n_mag::Int   = 300,
    imf = Utils.kroupa_imf,
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
            abs(iso_mass[k+1] - iso_mass[k-1]) / 2
        end
        w = imf(ms) * Δmass
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

# ╔═╡ 7cf639a4-cd5e-4f42-a0ac-a7cd768debdc
cmd_sat2 = build_cmd_filter(
	iso_new.gmag .- iso_new.rmag ,
	iso_new.gmag .+ distance_modulus,
	iso_new.Mini,
	x -> @.(sqrt(color_err.(x)^2 + 0.1^2)),
	x -> 0.02
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

# ╔═╡ 836339de-1267-4b14-aef9-56c37a62b0de
cmd_gaia = build_cmd_filter(
	iso_gaia_resampled.G_BPftmag .- iso_gaia_resampled.G_RPmag ,
	iso_gaia_resampled.Gmag .+ 19.60,
	iso_gaia_resampled.Mini,
	x -> @.(sqrt(err_gaia.(x)^2 + 0.1^2)),
	x -> 0.2,
	mag_range=(16, 21),
	color_range=(-0.5, 2.5)
)

# ╔═╡ a008c89c-3fae-4055-9d4f-1bfab2f79410
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 yreversed=true)
	y = (cmd_gaia[1] ./ maximum(cmd_gaia[1]))
	
	p = heatmap!(cmd_gaia[2], cmd_gaia[3], y,colorrange=(0, 1))

	Colorbar(fig[1,2], p, )
	fig
end

# ╔═╡ 75a1d70b-18d9-4b68-8646-6a21f7f36a43
function make_background_density(stars,     
	color_range::Tuple{Real,Real} = color_range,
    mag_range::Tuple{Real,Real}   = (15.0, 23.0),
    n_color::Int = 200,
	n_mag::Int = 300,
								 
)

	KernelDensity.kde((stars.gmag .- stars.rmag, stars.gmag), boundary=(color_range, mag_range), npoints=(n_color, n_mag),)
end

# ╔═╡ ff619ab8-3481-46bb-98b9-c706f0b1c507
dens_back = make_background_density(stars_bkg,)

# ╔═╡ 98829624-2eb3-4077-9995-cf45c52bfe81
heatmap(dens_back)

# ╔═╡ 79bb12ce-2a31-47c6-9265-2c5ec9dffa48
function cmd_axis(gs)
	return Axis(gs,
			   xlabel = "g - r",
			   ylabel = "g",
				yreversed = true,
			   )

end
	

# ╔═╡ a4159d32-8fd0-4577-ab14-8841ae328b52
let
	fig = Figure()
	ax = cmd_axis(fig[1,1])
	y = log10.(cmd_sat[1] ./ maximum(cmd_sat[1]))
	
	p = heatmap!(cmd_sat[2], cmd_sat[3], y,colorrange=(-3, 0))

	Colorbar(fig[1,2], p, label="log density")
	fig
end

# ╔═╡ 62a063a7-9e08-45cd-b187-b2c9cf587197
let
	fig = Figure()
	ax = cmd_axis(fig[1,1])
	y = log10.(cmd_gaia[1] ./ maximum(cmd_gaia[1]))
	
	p = heatmap!(cmd_gaia[2], cmd_gaia[3], y,colorrange=(-3, 0))

	Colorbar(fig[1,2], p, label="log density")
	fig
end

# ╔═╡ 07f0c3ba-964c-45d1-adce-b1c343c85d9f
let
	fig = Figure()
	
	ax = cmd_axis(fig[1,1])

	cmd = cmd_gaia3
	y = log10.(cmd[1] ./ maximum(cmd[1]))
	
	p = heatmap!(cmd[2], cmd[3], y, colorrange=(-3, 0))

	Colorbar(fig[1,2], p, label="density")
	fig
end

# ╔═╡ 4379bff6-5686-4709-b5b7-e8b4cced8664
let
	fig = Figure()
	ax = cmd_axis(fig[1,1])
		
	heatmap!(cmd_sat[2], cmd_sat[3], log10.(cmd_sat[1] ./ dens_back.density), colorrange=(-5, 1))

	
	# scatter!(stars_cen.gmag .- stars_cen.rmag, stars_cen.gmag, markersize=2, color=COLORS[3])

	fig

end

# ╔═╡ 11ce9995-58ac-4642-851b-092ce72c006e
import Interpolations

# ╔═╡ 767b40cd-5533-4fea-8a63-20e4c4d0d038
function make_ll_map(xbins, ybins, dens_sat, dens_back)
	dens_sat = dens_sat ./ sum(dens_sat)
	dens_back = dens_back ./ sum(dens_back)
	itp = Interpolations.interpolate((midpoints(xbins), midpoints(ybins)), 
							   dens_sat ./ dens_back, Interpolations.Gridded(Interpolations.Linear()))	

	return Interpolations.extrapolate(itp, 0)
end

# ╔═╡ 0870f9ab-56aa-4208-964b-32257048f38d


# ╔═╡ a4e94b36-dd41-4f05-9b16-007639650bd5
ll_map = make_ll_map(cmd_sat[2], cmd_sat[3], cmd_sat[1], dens_back.density)

# ╔═╡ 4b2ffbee-9820-4e48-8d1f-037207887093
heatmap(log10.(ll_map.itp.coefs), colorrange=(-5, 1))

# ╔═╡ eef36a62-34ce-4d46-9fce-6810cc430042
LLR = ll_map.(stars.gmag .- stars.rmag, stars.gmag)

# ╔═╡ b1857714-b3a4-4a38-a816-d288ba950dbd
let
	fig = Figure()
	ax = cmd_axis(fig[1,1])
		
	p = scatter!(stars.gmag .- stars.rmag, stars.gmag, color=LLR, markersize=1, 
				colorrange=(0, 100), colorscale=asinh)

	Colorbar(fig[1,2], p, label="LLR", minorticksvisible=false)

	fig

end
	

# ╔═╡ 663c10ad-e3b6-4f86-bdb4-97f12f163c53
let
	fig = Figure()
	ax = cmd_axis(fig[1,1])
		
	p = scatter!(stars.gmag .- stars.rmag, stars.gmag, color=LLR, markersize=1, 
				colorrange=(0, 100), colorscale=asinh)

	Colorbar(fig[1,2], p, label="LLR", minorticksvisible=false)

	fig

end
	

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

# ╔═╡ 0590d501-ff35-4a52-862b-e023143ea911
md"""
# Next stage
"""

# ╔═╡ 0942d5fd-f1b0-456a-883d-0d74494353af
sample = stars[LLR .> 1, :]

# ╔═╡ 54510ba4-4ecb-4c77-9437-2c200e3e8a60
scatter(sample.xi, sample.eta, markersize=1, color=:black)

# ╔═╡ 11f6e923-b608-43a8-bd7a-7b8b668549c0
hexbin(sample.xi, sample.eta, bins=80, colorscale=log10)

# ╔═╡ cf349fec-711c-4c58-86b2-74771e00f4bc
md"""
# Gaia
"""

# ╔═╡ Cell order:
# ╠═4c482a6f-9cd3-4704-b8ce-b23e0f1a293f
# ╠═015eae5d-515c-4fa4-9a44-fc478692cd44
# ╠═68cd2e80-a8d7-4fcd-87dd-a632818b77ac
# ╠═34c475e5-1eb9-482a-a29b-182caad50313
# ╠═e3cb612c-7c1f-4688-8aa0-3cc4de19ec9b
# ╠═e560b3a4-249d-11f1-a7d8-918587896934
# ╠═d045eaad-a22c-4def-8000-0f767b72ffed
# ╠═0a8c2527-6812-4293-b60e-8d91bed66712
# ╠═7eab2c84-6335-475b-be25-e8c2fce07776
# ╠═e8c0e80d-336e-4f23-b826-278d90d732d1
# ╠═6ed9796f-b522-4e62-82b5-62f0297c1e35
# ╠═6c55e84e-39a9-444b-a067-95a3e3826cbd
# ╠═7cf639a4-cd5e-4f42-a0ac-a7cd768debdc
# ╠═109fdcda-b64e-477f-bfba-4f9588b7bf09
# ╠═52792fb8-d9cf-4b2c-81f8-e3d3940c2250
# ╠═feceea26-7bc7-4789-9f4b-4df4c5c4c27a
# ╠═c9b72c30-a6a8-45d0-8b89-86174370ccce
# ╠═df5e5629-8d6b-4cb7-92bf-bf4467b6c2f2
# ╠═836339de-1267-4b14-aef9-56c37a62b0de
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
# ╠═62a063a7-9e08-45cd-b187-b2c9cf587197
# ╠═a7d83f97-159b-404f-b2f2-b595a9fe9ee7
# ╠═53fc7954-4d68-4c78-bf8c-75ad7577b0e6
# ╠═07f0c3ba-964c-45d1-adce-b1c343c85d9f
# ╠═9ed01175-a59c-4b27-8b6c-5d286ea6faad
# ╠═a008c89c-3fae-4055-9d4f-1bfab2f79410
# ╠═bb0600b7-cf8a-484b-9d4a-7dab2c963da6
# ╠═764be90a-e602-44d8-9956-d407ee2a0d4e
# ╠═39d177ca-6402-4216-b503-c12f83f91d56
# ╠═48993623-e10e-45c7-97df-93e3e0200007
# ╠═bf0e16fe-94cd-41fd-a107-4a814dd448f4
# ╠═080d5e15-dd13-423a-815c-e22edfc73dcc
# ╠═71dfe223-539a-4e2b-84b1-04811b2533a0
# ╠═13896043-e020-409f-99fb-8d1254b994d2
# ╠═d222d357-8d82-494a-b85d-753314c044e4
# ╠═b8ae4ad0-afdf-4877-a02f-f9291d6411e9
# ╠═6dabd441-9f39-4e75-a651-7edba2da41e4
# ╠═cfc6be95-7c53-479f-81fc-ce4d9ca07b41
# ╠═75a1d70b-18d9-4b68-8646-6a21f7f36a43
# ╠═ff619ab8-3481-46bb-98b9-c706f0b1c507
# ╠═98829624-2eb3-4077-9995-cf45c52bfe81
# ╠═79bb12ce-2a31-47c6-9265-2c5ec9dffa48
# ╠═4379bff6-5686-4709-b5b7-e8b4cced8664
# ╠═11ce9995-58ac-4642-851b-092ce72c006e
# ╠═767b40cd-5533-4fea-8a63-20e4c4d0d038
# ╠═0870f9ab-56aa-4208-964b-32257048f38d
# ╠═a4e94b36-dd41-4f05-9b16-007639650bd5
# ╠═4b2ffbee-9820-4e48-8d1f-037207887093
# ╠═eef36a62-34ce-4d46-9fce-6810cc430042
# ╠═b1857714-b3a4-4a38-a816-d288ba950dbd
# ╠═663c10ad-e3b6-4f86-bdb4-97f12f163c53
# ╠═65973035-10d5-4ece-8b21-c51de60daf90
# ╠═0590d501-ff35-4a52-862b-e023143ea911
# ╠═0942d5fd-f1b0-456a-883d-0d74494353af
# ╠═54510ba4-4ecb-4c77-9437-2c200e3e8a60
# ╠═11f6e923-b608-43a8-bd7a-7b8b668549c0
# ╠═cf349fec-711c-4c58-86b2-74771e00f4bc
