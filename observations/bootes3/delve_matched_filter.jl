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

# ╔═╡ 6d3349e7-befc-40ca-981c-47ac3922c090
suffix = "_6deg"

# ╔═╡ c185afc3-f34b-4f75-b7e5-68c7d7e45a3a
iso_width = 0.05

# ╔═╡ 748de64f-58c3-4cd2-a86f-a9a32d7af339
mock = false

# ╔═╡ 7579731a-af68-4858-aa3d-08510fd4b5eb
mode = "g22"

# ╔═╡ 015eae5d-515c-4fa4-9a44-fc478692cd44
CairoMakie.activate!(type=:png)

# ╔═╡ e8c0e80d-336e-4f23-b826-278d90d732d1
import TOML

# ╔═╡ e560b3a4-249d-11f1-a7d8-918587896934
module Utils
	include("delve_utils.jl")
end

# ╔═╡ cfd3c3df-38b7-4f78-b798-a6cbc3bdceef
if mode == "all"
	mag_range = (15, 23.0)
	color_range = (-0.7, 1.5)
elseif mode == "g22"
	mag_range = (15, 22.0)
	color_range = (-0.7, 1.5)
elseif mode == "bhb"
	mag_range = (18, 19.5)
	color_range = (-0.7, 0.2)
elseif mode == "rgb"
	mag_range = (15, 21)
	color_range = (-0.7, 1.5)
elseif mode == "ms"
	mag_range = (21, 23.0)
	color_range = (-0.7, 1.5)
end

# ╔═╡ b1441903-47cc-471c-ae5f-04e0ba880aec
md"""
# Main drag
"""

# ╔═╡ 6ed9796f-b522-4e62-82b5-62f0297c1e35
obs_props = TOML.parsefile("observed_properties.toml")

# ╔═╡ 7cf639a4-cd5e-4f42-a0ac-a7cd768debdc
distance_modulus_err = obs_props["distance_modulus_em"]

# ╔═╡ 6c55e84e-39a9-444b-a067-95a3e3826cbd
distance_modulus =  obs_props["distance_modulus"]

# ╔═╡ cc053f1a-0151-4a18-91b2-8e403cb9d148
if mock
	ra0, dec0 = 205, 22
else
	ra0, dec0 = obs_props["ra_original"], obs_props["dec_original"]
end

# ╔═╡ 90a3ba08-28eb-4529-bfb7-c6d01b3dd979
R_h = obs_props["R_h"]

# ╔═╡ fdc022f6-f3f8-419c-a42e-4d695f00536d
@assert obs_props["distance_modulus_em"] == obs_props["distance_modulus_ep"]

# ╔═╡ f63389bf-77fe-4495-9be5-b4f86f84dc57
stars = let
	df = read_fits("data/delve_dr2_good$suffix.fits")
	df[!, :bp_rp] = df.gmag .- df.rmag
	df[!, :G] = df.gmag

	df[!, :xi], df[!, :eta] = 60 .* LilGuys.to_tangent(df.ra, df.dec, ra0, dec0)
	df[!, :R] = @. sqrt(df.xi^2 + df.eta^2)
	df
end

# ╔═╡ b346a53e-6162-4c5d-accf-4c59e327f58f
iso = Utils.get_isochrone(-2.1, 12, 
					 stage_max = mode=="rgb" ? 3 : 4
					)

# ╔═╡ 78b5414f-8961-4446-bfb7-9c0e852dd5b7
iso_new = Utils.resample_iso(iso, N=mode=="ms" ? 10_000 : 1_000_000)

# ╔═╡ 48993623-e10e-45c7-97df-93e3e0200007
color_err = Utils.delve_gr_err

# ╔═╡ 109fdcda-b64e-477f-bfba-4f9588b7bf09
cmd_sat = Utils.build_cmd_filter(
	iso_new.gmag .- iso_new.rmag ,
	iso_new.gmag .+ distance_modulus,
	iso_new.Mini,
	x -> @.(sqrt(color_err.(x)^2 + iso_width^2)),
	x -> distance_modulus_err,
	color_range=color_range,
	mag_range=mag_range,
)

# ╔═╡ 22934ab9-2750-4826-b916-a51ecee71f68
stars_bkg = stars[(stars.R .> 5*60), :]

# ╔═╡ ff619ab8-3481-46bb-98b9-c706f0b1c507
cmd_bkg = Utils.make_background_density(stars_bkg, cmd_sat.x, cmd_sat.y)

# ╔═╡ 8387e157-6963-422f-825b-bdebd5f795b5
@assert cmd_bkg.x ≈ cmd_sat.x

# ╔═╡ ed71897a-45a3-4b44-ab0f-fc8be5c584c2
@assert cmd_bkg.y ≈ cmd_sat.y

# ╔═╡ a576f974-d294-43de-a04a-6254b65343e5
l_sat = Utils.interpolate_density(cmd_sat)

# ╔═╡ 5353e8a6-abe0-431a-8dae-422b84430535
l_bkg = Utils.interpolate_density(cmd_bkg)

# ╔═╡ f61498bf-cc1f-48e1-981c-2e3b7ac88ede
df_out = let
	df = copy(stars)
	df[!, :L_sat] = l_sat.(stars.gmag .- stars.rmag, stars.gmag)
	df[!, :L_bg] = l_bkg.(stars.gmag .- stars.rmag, stars.gmag)
	df[!, :LLR_nospace] = @. log10(df.L_sat / df.L_bg)
	df[!, :F_BEST] = isfinite.(df.LLR_nospace)
	df
end

# ╔═╡ eef36a62-34ce-4d46-9fce-6810cc430042
LLR = log10.(df_out.L_sat ./ df_out.L_bg)

# ╔═╡ 263fd125-cb3c-4c33-80d9-8c8496f1ff41
let
	dx = diff(cmd_sat.x)[1]
	dy = diff(cmd_sat.y)[1]
	@assert sum(dx * dy * cmd_sat.density) ≈ 1
	@assert sum(dx * dy * cmd_bkg.density) ≈ 1
end

# ╔═╡ e9fa2fb1-e9f5-4f00-8ae7-f28dc40ff665
sum(df_out.F_BEST)

# ╔═╡ eafcd5f0-a300-4440-875e-ce618b449904
modesuffix = mode=="all" ? "" : ".$mode"

# ╔═╡ 49204205-cf52-4489-a2bd-972e1112b0a8
write_fits("samples/delve_matched_filter$suffix$modesuffix.fits", df_out, overwrite=true)

# ╔═╡ 708d7ad5-d284-41ac-8aff-2e74891e9f54
md"""
# Plots
"""

# ╔═╡ baccd40f-4128-4031-973b-1389845f2972
stars_cen = stars[stars.R .< R_h, :]

# ╔═╡ 79bb12ce-2a31-47c6-9265-2c5ec9dffa48
function cmd_axis(gs)
	return Axis(gs,
			   xlabel = L"$g - r$ [mag]",
			   ylabel = L"$g$ [mag]",
				yreversed = true,
				limits=(color_range, mag_range)
			   )

end
	

# ╔═╡ 98829624-2eb3-4077-9995-cf45c52bfe81
hist(stars.gmag, bins=100)

# ╔═╡ c2fc2020-35e7-406f-be56-6c9bbd994370
hist(stars.rmag, bins=100)

# ╔═╡ 88a5e8a3-00dd-47b3-ba2f-ab3ffc483c25
hexbin(stars_bkg.xi, stars_bkg.eta, bins=300)

# ╔═╡ c47505fc-01f6-4252-815a-2be038d6d3bc
hexbin(stars.xi, stars.eta, bins=300)

# ╔═╡ 0c964ad8-cc76-49cc-922c-12d20e5d259f
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"\xi / '", ylabel=L"\eta / '", aspect=DataAspect(), xreversed=true)
	scatter!(stars.xi, stars.eta, markersize=1, color=:black, alpha=0.05)
	fig
end

# ╔═╡ 4e7dacb3-366e-46f4-84f4-fd6206961c4e
let
	fig = Figure()
	ax = Axis(fig[1,1])
		
	hexbin!(stars.xi, stars.eta, bins=300)

	arc!((0,0), 30, 0, 2π)

	arc!((0,0), 5*R_h, 0, 2π)
	arc!((0,0), 10R_h, 0, 2π)
	fig
end

# ╔═╡ ba636601-87d4-4bde-b8a7-9e3102dbba6b


# ╔═╡ a75c32b1-9011-4c88-8725-624f93fd0d98
let
	fig = Figure()

	ax = cmd_axis(fig[1,1])
	y = log10.(cmd_sat.density )
	
	p = heatmap!(cmd_sat.x, cmd_sat.y, y,colorrange=(-3, 0))


	scatter!(stars_cen.gmag .- stars_cen.rmag, stars_cen.gmag, markersize=2, color=COLORS[3])
	fig
end

# ╔═╡ a4159d32-8fd0-4577-ab14-8841ae328b52
let
	fig = Figure()
	ax = cmd_axis(fig[1,1])
	y = log10.(cmd_sat.density )
	
	p = heatmap!(cmd_sat.x, cmd_sat.y, y,colorrange=(-3, 0))

	Colorbar(fig[1,2], p, label="log density")
	fig
end

# ╔═╡ 4379bff6-5686-4709-b5b7-e8b4cced8664
let
	fig = Figure(size=(6, 3) .* 72)
	ax = cmd_axis(fig[1,1])
	ax.title = "Boo III model"

	heatmap!(cmd_sat.x, cmd_sat.y, cmd_sat.density)
	
	ax = cmd_axis(fig[1,2])
	ax.title = "Background"

	heatmap!(cmd_bkg.x, cmd_bkg.y, cmd_bkg.density)
	hideydecorations!()


	ax = cmd_axis(fig[1,3])
	ax.title = L"L_\text{sat} / L_\text{bg}"

		
	heatmap!(cmd_sat.x, cmd_sat.y, (cmd_sat.density ./ cmd_bkg.density), colorrange=(0, 10))
	hideydecorations!()

	@savefig "matched_filter_model"
	fig

end

# ╔═╡ b1857714-b3a4-4a38-a816-d288ba950dbd
let
	fig = Figure()
	ax = cmd_axis(fig[1,1])
		
	p = scatter!(stars.gmag .- stars.rmag, stars.gmag, color=df_out.LLR_nospace, markersize=1, 
				colorrange=(-2, 2))

	Colorbar(fig[1,2], p, label="LLR", minorticksvisible=false)

	fig

end
	

# ╔═╡ 9e3ce97d-67f5-4a8d-b36d-d20e051715ad
let
	fig = Figure()
	ax = cmd_axis(fig[1,1])
		
	p = scatter!(stars.gmag .- stars.rmag, stars.gmag, color=log10.(df_out.L_bg), markersize=1, 
				colorrange=(-8, 2))

	Colorbar(fig[1,2], p, label="LLR", minorticksvisible=false)

	fig

end
	

# ╔═╡ 17522c59-622b-4fe3-88af-f48308f3137b
let
	fig = Figure()
	ax = cmd_axis(fig[1,1])
		
	p = scatter!(stars.gmag .- stars.rmag, stars.gmag, color=df_out.L_sat, markersize=1, 
				colorrange=(0, 1))

	Colorbar(fig[1,2], p, label="LR", minorticksvisible=false)

	fig

end
	

# ╔═╡ f884701a-80e7-458c-88d7-db626254608e
let
	fig = Figure()
	ax = cmd_axis(fig[1,1])
		
	p = scatter!(stars.gmag .- stars.rmag, stars.gmag, color=10 .^ LLR, markersize=1, 
				colorrange=(0, 10))

	Colorbar(fig[1,2], p, label="LR", minorticksvisible=false)

	fig

end
	

# ╔═╡ 1db09bab-cdc0-49e5-9b18-8b0eb82e3300
let
	fig = Figure()
	ax = cmd_axis(fig[1,1])
		
	p = scatter!(stars.gmag .- stars.rmag, stars.gmag, color=LLR .> 2, markersize=1, )

	Colorbar(fig[1,2], p, label="LR", minorticksvisible=false)

	fig

end
	

# ╔═╡ 2a76de21-3027-43ea-b947-56475a30e111
hist(LLR, bins=LinRange(-4, 2, 100))

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

# ╔═╡ Cell order:
# ╠═6d3349e7-befc-40ca-981c-47ac3922c090
# ╠═c185afc3-f34b-4f75-b7e5-68c7d7e45a3a
# ╠═748de64f-58c3-4cd2-a86f-a9a32d7af339
# ╠═7cf639a4-cd5e-4f42-a0ac-a7cd768debdc
# ╠═6c55e84e-39a9-444b-a067-95a3e3826cbd
# ╠═7579731a-af68-4858-aa3d-08510fd4b5eb
# ╠═4c482a6f-9cd3-4704-b8ce-b23e0f1a293f
# ╠═015eae5d-515c-4fa4-9a44-fc478692cd44
# ╠═34c475e5-1eb9-482a-a29b-182caad50313
# ╠═e8c0e80d-336e-4f23-b826-278d90d732d1
# ╠═e3cb612c-7c1f-4688-8aa0-3cc4de19ec9b
# ╠═e560b3a4-249d-11f1-a7d8-918587896934
# ╠═cfd3c3df-38b7-4f78-b798-a6cbc3bdceef
# ╠═cc053f1a-0151-4a18-91b2-8e403cb9d148
# ╟─b1441903-47cc-471c-ae5f-04e0ba880aec
# ╠═6ed9796f-b522-4e62-82b5-62f0297c1e35
# ╠═90a3ba08-28eb-4529-bfb7-c6d01b3dd979
# ╠═fdc022f6-f3f8-419c-a42e-4d695f00536d
# ╠═f63389bf-77fe-4495-9be5-b4f86f84dc57
# ╠═b346a53e-6162-4c5d-accf-4c59e327f58f
# ╠═78b5414f-8961-4446-bfb7-9c0e852dd5b7
# ╠═48993623-e10e-45c7-97df-93e3e0200007
# ╠═109fdcda-b64e-477f-bfba-4f9588b7bf09
# ╠═22934ab9-2750-4826-b916-a51ecee71f68
# ╠═ff619ab8-3481-46bb-98b9-c706f0b1c507
# ╠═8387e157-6963-422f-825b-bdebd5f795b5
# ╠═ed71897a-45a3-4b44-ab0f-fc8be5c584c2
# ╠═a576f974-d294-43de-a04a-6254b65343e5
# ╠═5353e8a6-abe0-431a-8dae-422b84430535
# ╠═eef36a62-34ce-4d46-9fce-6810cc430042
# ╠═f61498bf-cc1f-48e1-981c-2e3b7ac88ede
# ╠═263fd125-cb3c-4c33-80d9-8c8496f1ff41
# ╠═e9fa2fb1-e9f5-4f00-8ae7-f28dc40ff665
# ╠═eafcd5f0-a300-4440-875e-ce618b449904
# ╠═49204205-cf52-4489-a2bd-972e1112b0a8
# ╟─708d7ad5-d284-41ac-8aff-2e74891e9f54
# ╠═baccd40f-4128-4031-973b-1389845f2972
# ╠═79bb12ce-2a31-47c6-9265-2c5ec9dffa48
# ╠═98829624-2eb3-4077-9995-cf45c52bfe81
# ╠═c2fc2020-35e7-406f-be56-6c9bbd994370
# ╠═88a5e8a3-00dd-47b3-ba2f-ab3ffc483c25
# ╠═c47505fc-01f6-4252-815a-2be038d6d3bc
# ╠═0c964ad8-cc76-49cc-922c-12d20e5d259f
# ╠═4e7dacb3-366e-46f4-84f4-fd6206961c4e
# ╠═ba636601-87d4-4bde-b8a7-9e3102dbba6b
# ╠═a75c32b1-9011-4c88-8725-624f93fd0d98
# ╠═a4159d32-8fd0-4577-ab14-8841ae328b52
# ╠═4379bff6-5686-4709-b5b7-e8b4cced8664
# ╠═b1857714-b3a4-4a38-a816-d288ba950dbd
# ╠═9e3ce97d-67f5-4a8d-b36d-d20e051715ad
# ╠═17522c59-622b-4fe3-88af-f48308f3137b
# ╠═f884701a-80e7-458c-88d7-db626254608e
# ╠═1db09bab-cdc0-49e5-9b18-8b0eb82e3300
# ╠═2a76de21-3027-43ea-b947-56475a30e111
# ╠═ed8b8c35-7afb-457d-a5c6-d9a16a497937
# ╠═c62a3752-30d6-425e-a8cd-13347bbf1983
# ╠═72c9db08-0c9e-4188-a300-8f6034029c66
# ╠═7a9a905a-8b3c-490c-ab9a-c8b3c581d1a7
# ╠═89f8f6de-4569-4b96-9f3c-15841e98b49f
# ╠═8c90c029-aedc-4745-84ed-911efe78b8d5
# ╠═7d1c5c26-0a44-42a6-b15b-c708dcd3a505
# ╠═0886ca9e-b31e-465a-aab9-fc46068ae91b
# ╠═d28e2b85-79cd-4473-a735-6e46c77ec1eb
# ╠═f68df67c-50fc-4bf3-8a2a-d5a810f869be
# ╠═341d33bf-fc9c-4264-8513-05d7fa4ee45b
# ╠═021bc8cb-0449-4cae-898a-b7258bf7d6ef
# ╠═1d94e532-d99e-47e1-ab62-8b3c5692f4ee
# ╠═89f9f774-b1af-444c-bbba-059cdfe278e5
# ╠═bb2ec2cc-4e82-42c6-a0a4-f2abedfaf2ce
# ╠═0b1aab30-6042-4d59-bfef-e1ed7b17082e
# ╠═c0459b64-2393-45dd-aac0-9d193866134a
