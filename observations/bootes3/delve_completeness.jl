### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 4c482a6f-9cd3-4704-b8ce-b23e0f1a293f
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya, CairoMakie
	FIGDIR = "./figures"
end

# ╔═╡ 34c475e5-1eb9-482a-a29b-182caad50313
using DataFrames

# ╔═╡ e3cb612c-7c1f-4688-8aa0-3cc4de19ec9b
using PyFITS

# ╔═╡ 015eae5d-515c-4fa4-9a44-fc478692cd44
CairoMakie.activate!(type=:png)

# ╔═╡ e8c0e80d-336e-4f23-b826-278d90d732d1
import TOML

# ╔═╡ e560b3a4-249d-11f1-a7d8-918587896934
module Utils
	include("delve_utils.jl")
end

# ╔═╡ 6ed9796f-b522-4e62-82b5-62f0297c1e35
obs_props = TOML.parsefile("observed_properties.toml")

# ╔═╡ cc053f1a-0151-4a18-91b2-8e403cb9d148
ra0, dec0 = obs_props["ra_original"], obs_props["dec_original"]

# ╔═╡ f63389bf-77fe-4495-9be5-b4f86f84dc57
stars = let
	df = read_fits("data/delve_dr2_good_6deg.fits")
	df[!, :bp_rp] = df.gmag .- df.rmag
	df[!, :G] = df.gmag

	df[!, :xi], df[!, :eta] = 60 .* LilGuys.to_tangent(df.ra, df.dec, ra0, dec0)
	df[!, :R] = @. sqrt(df.xi^2 + df.eta^2)
	df
end

# ╔═╡ 90a3ba08-28eb-4529-bfb7-c6d01b3dd979
R_h = obs_props["R_h"]

# ╔═╡ ed8b8c35-7afb-457d-a5c6-d9a16a497937
md"""
## Completeness
"""

# ╔═╡ c62a3752-30d6-425e-a8cd-13347bbf1983
import KernelDensity: kde

# ╔═╡ 72c9db08-0c9e-4188-a300-8f6034029c66
function kernel_mode(data::AbstractVector{<:Real}; bandwidth::Real)
    k = kde(data; bandwidth=bandwidth)
    mode_idx = argmax(k.density)
    return k.x[mode_idx]
end

# ╔═╡ 7a9a905a-8b3c-490c-ab9a-c8b3c581d1a7
kernel_mode(stars.gmag, bandwidth=0.05)

# ╔═╡ 89f8f6de-4569-4b96-9f3c-15841e98b49f
import StatsBase: median

# ╔═╡ 8c90c029-aedc-4745-84ed-911efe78b8d5
function five_sigma_depth(values)
	mags = first.(values)
	errs = last.(values)

	magerr_0 = 0.2
	magerr_w = 0.05
	filt = magerr_0 - magerr_w .<= errs .<= magerr_0 + magerr_w
	if sum(filt) == 0
		return NaN
	else
		return median(mags[filt])
	end
end

# ╔═╡ 7d1c5c26-0a44-42a6-b15b-c708dcd3a505
five_sigma_depth(zip(stars.gmag, stars.magerr_psf_g))

# ╔═╡ 0886ca9e-b31e-465a-aab9-fc46068ae91b
five_sigma_depth(zip(stars.rmag, stars.magerr_psf_r))

# ╔═╡ d28e2b85-79cd-4473-a735-6e46c77ec1eb
function binned_stats(x, y, z, stat_func; x_bins=45, y_bins=45)
    x_edges = range(minimum(x), maximum(x), length=x_bins + 1)
    y_edges = range(minimum(y), maximum(y), length=y_bins + 1)

    mode_grid = fill(NaN, x_bins, y_bins)

    for i in 1:x_bins
        for j in 1:y_bins
            mask = (x_edges[i] .<= x .< x_edges[i+1]) .&
                   (y_edges[j] .<= y .< y_edges[j+1])
            z_vals = z[mask]
			
            if length(z_vals) >= 2
                mode_grid[i, j] = stat_func(z_vals)
            end
        end
    end

    x_centers = collect(0.5 .* (x_edges[1:end-1] .+ x_edges[2:end]))
    y_centers = collect(0.5 .* (y_edges[1:end-1] .+ y_edges[2:end]))

    return Utils.DensityMap(x_centers, y_centers, mode_grid)
end

# ╔═╡ f68df67c-50fc-4bf3-8a2a-d5a810f869be
extrema(stars.eta)

# ╔═╡ 11d6dc5e-b762-4f99-9365-710d0015e580
function tangent_axis(gs)
	ax = Axis(gs,
			 xlabel = "xi / arcmin",
			 ylabel = "eta / arcmin",
			 xreversed = true, 
			 aspect=DataAspect())
end

# ╔═╡ 341d33bf-fc9c-4264-8513-05d7fa4ee45b
function plot_density_map(density_map; colorlabel="", kwargs...)
	fig = Figure()

	ax = Axis(fig[1,1],
			 xlabel = "xi / arcmin",
			 ylabel = "eta / arcmin",
			 xreversed = true, 
			 aspect=DataAspect())

	p = heatmap!(density_map.x, density_map.y, density_map.density; kwargs...)

	@info minimum(density_map.density[isfinite.(density_map.density)])
	@info diff(density_map.x)[1]
	Colorbar(fig[1,2], p, label=colorlabel)
	fig
end

# ╔═╡ 021bc8cb-0449-4cae-898a-b7258bf7d6ef
let 
	modes = binned_stats(stars.xi, stars.eta, stars.gmag, x -> kernel_mode(x, bandwidth=0.2))
		
	plot_density_map(modes, colorrange=(23.4, 24) .- 0.2)
end

# ╔═╡ 1d94e532-d99e-47e1-ab62-8b3c5692f4ee
let
	m = binned_stats(stars.xi, stars.eta, collect(zip(stars.gmag, stars.magerr_psf_g)), five_sigma_depth,)
	
	plot_density_map(m, colorlabel="5σ depth g",  colorrange=(23.5, 24))
end

# ╔═╡ 89f9f774-b1af-444c-bbba-059cdfe278e5
let
	m = binned_stats(stars.xi, stars.eta, collect(zip(stars.gmag, stars.magerr_psf_g)), five_sigma_depth)

	hist(vec(m.density[isfinite.(m.density)]), bins=60)
end

# ╔═╡ bb2ec2cc-4e82-42c6-a0a4-f2abedfaf2ce
let
	m = binned_stats(stars.xi, stars.eta, collect(zip(stars.rmag, stars.magerr_psf_r)), five_sigma_depth)

	hist(vec(m.density[isfinite.(m.density)]),  bins=60)
end

# ╔═╡ 0b1aab30-6042-4d59-bfef-e1ed7b17082e
let
	m = binned_stats(stars.xi, stars.eta, collect(zip(stars.rmag, stars.magerr_psf_r)), five_sigma_depth)

	plot_density_map(m, colorlabel="5σ depth r", colorrange=(22.9, 23.5))
end

# ╔═╡ c0459b64-2393-45dd-aac0-9d193866134a
let 
	modes = binned_stats(stars.xi, stars.eta, stars.rmag, x -> kernel_mode(x, bandwidth=0.1), x_bins=30, y_bins=30)
		
	plot_density_map(modes, colorrange=(22.2, 23.5) .- 0.3)
end

# ╔═╡ 642505d1-ef26-4a68-ba4a-2a39c3c06d6b
md"""
# Using dataframe properties
"""

# ╔═╡ 1bf9bdb4-418e-4333-936b-0b8ec18a735c
let
	fig = Figure()
	ax = tangent_axis(fig[1,1])
	p = scatter!(stars.xi, stars.eta, color=stars.nepochs_g, markersize=1)

	Colorbar(fig[1,2], p, label="# detections g")
	fig
end

# ╔═╡ c8484d13-4755-452e-8891-9419ff205b14
let
	fig = Figure()
	ax = tangent_axis(fig[1,1])
	p = scatter!(stars.xi, stars.eta, color=stars.nepochs_r, markersize=1)

	Colorbar(fig[1,2], p, label="# detections r")
	fig
end

# ╔═╡ d2a9b1a7-32b9-4fa0-9f7e-2838d2ef0f5c
let
	fig = Figure()
	ax = tangent_axis(fig[1,1])
	p = scatter!(stars.xi, stars.eta, color=stars.t_eff_g, markersize=1)

	Colorbar(fig[1,2], p, label="eff t g")
	fig
end

# ╔═╡ 8f780938-3a40-4239-ab8f-a1f233d7549c
let
	fig = Figure()
	ax = tangent_axis(fig[1,1])
	p = scatter!(stars.xi, stars.eta, color=stars.t_eff_r, markersize=1)

	Colorbar(fig[1,2], p, label="eff t r")
	fig
end

# ╔═╡ 25296d70-5981-41b9-898b-f7522741fb78


# ╔═╡ 198a4964-15bb-475c-a178-935746e72e7a


# ╔═╡ Cell order:
# ╠═4c482a6f-9cd3-4704-b8ce-b23e0f1a293f
# ╠═015eae5d-515c-4fa4-9a44-fc478692cd44
# ╠═34c475e5-1eb9-482a-a29b-182caad50313
# ╠═e8c0e80d-336e-4f23-b826-278d90d732d1
# ╠═e3cb612c-7c1f-4688-8aa0-3cc4de19ec9b
# ╠═e560b3a4-249d-11f1-a7d8-918587896934
# ╠═cc053f1a-0151-4a18-91b2-8e403cb9d148
# ╠═f63389bf-77fe-4495-9be5-b4f86f84dc57
# ╠═6ed9796f-b522-4e62-82b5-62f0297c1e35
# ╠═90a3ba08-28eb-4529-bfb7-c6d01b3dd979
# ╟─ed8b8c35-7afb-457d-a5c6-d9a16a497937
# ╠═c62a3752-30d6-425e-a8cd-13347bbf1983
# ╠═72c9db08-0c9e-4188-a300-8f6034029c66
# ╠═7a9a905a-8b3c-490c-ab9a-c8b3c581d1a7
# ╠═89f8f6de-4569-4b96-9f3c-15841e98b49f
# ╠═8c90c029-aedc-4745-84ed-911efe78b8d5
# ╠═7d1c5c26-0a44-42a6-b15b-c708dcd3a505
# ╠═0886ca9e-b31e-465a-aab9-fc46068ae91b
# ╠═d28e2b85-79cd-4473-a735-6e46c77ec1eb
# ╠═f68df67c-50fc-4bf3-8a2a-d5a810f869be
# ╠═11d6dc5e-b762-4f99-9365-710d0015e580
# ╠═341d33bf-fc9c-4264-8513-05d7fa4ee45b
# ╠═021bc8cb-0449-4cae-898a-b7258bf7d6ef
# ╠═1d94e532-d99e-47e1-ab62-8b3c5692f4ee
# ╠═89f9f774-b1af-444c-bbba-059cdfe278e5
# ╠═bb2ec2cc-4e82-42c6-a0a4-f2abedfaf2ce
# ╠═0b1aab30-6042-4d59-bfef-e1ed7b17082e
# ╠═c0459b64-2393-45dd-aac0-9d193866134a
# ╠═642505d1-ef26-4a68-ba4a-2a39c3c06d6b
# ╠═1bf9bdb4-418e-4333-936b-0b8ec18a735c
# ╠═c8484d13-4755-452e-8891-9419ff205b14
# ╠═d2a9b1a7-32b9-4fa0-9f7e-2838d2ef0f5c
# ╠═8f780938-3a40-4239-ab8f-a1f233d7549c
# ╠═25296d70-5981-41b9-898b-f7522741fb78
# ╠═198a4964-15bb-475c-a178-935746e72e7a
