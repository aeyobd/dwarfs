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

# ╔═╡ c185afc3-f34b-4f75-b7e5-68c7d7e45a3a
iso_width = 0.05

# ╔═╡ 11ce9995-58ac-4642-851b-092ce72c006e
import Interpolations

# ╔═╡ 015eae5d-515c-4fa4-9a44-fc478692cd44
CairoMakie.activate!(type=:png)

# ╔═╡ 68cd2e80-a8d7-4fcd-87dd-a632818b77ac
import KernelDensity

# ╔═╡ e8c0e80d-336e-4f23-b826-278d90d732d1
import TOML

# ╔═╡ e560b3a4-249d-11f1-a7d8-918587896934
module Utils
	include("delve_utils.jl")
end

# ╔═╡ 6ed9796f-b522-4e62-82b5-62f0297c1e35
obs_props = TOML.parsefile("observed_properties.toml")

# ╔═╡ 6c55e84e-39a9-444b-a067-95a3e3826cbd
distance_modulus = obs_props["distance_modulus"]

# ╔═╡ fdc022f6-f3f8-419c-a42e-4d695f00536d
@assert obs_props["distance_modulus_em"] == obs_props["distance_modulus_ep"]

# ╔═╡ 7cf639a4-cd5e-4f42-a0ac-a7cd768debdc
distance_modulus_err = obs_props["distance_modulus_em"]

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

# ╔═╡ c47505fc-01f6-4252-815a-2be038d6d3bc
hexbin(stars.xi, stars.eta, bins=300)

# ╔═╡ 0c964ad8-cc76-49cc-922c-12d20e5d259f
scatter(stars.xi, stars.eta, markersize=1, color=:black, alpha=0.05)

# ╔═╡ 90a3ba08-28eb-4529-bfb7-c6d01b3dd979
R_h = obs_props["R_h"]

# ╔═╡ 22934ab9-2750-4826-b916-a51ecee71f68
stars_bkg = stars[(stars.R .> 5R_h) .& (stars.eta .< 80), :]

# ╔═╡ 88a5e8a3-00dd-47b3-ba2f-ab3ffc483c25
hexbin(stars_bkg.xi, stars_bkg.eta, bins=300)

# ╔═╡ c94244bd-31f0-4838-8495-05fc8be234e8
hist(stars.gmag)

# ╔═╡ a4159d32-8fd0-4577-ab14-8841ae328b52


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
	scatter!(stars.gmag, sqrt.(stars.magerr_psf_g .^2 .+ stars.magerr_psf_r .^2), markersize=0.3, color=:black, alpha=0.03)

	scatter!(stars.gmag, color_err.(stars.gmag), markersize=1)
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

# ╔═╡ 109fdcda-b64e-477f-bfba-4f9588b7bf09
cmd_sat = Utils.build_cmd_filter(
	iso_new.gmag .- iso_new.rmag ,
	iso_new.gmag .+ distance_modulus,
	iso_new.Mini,
	x -> @.(sqrt(color_err.(x)^2 + iso_width^2)),
	x -> 0.19
)

# ╔═╡ a75c32b1-9011-4c88-8725-624f93fd0d98
let
	fig = Figure()

	ax = Axis(fig[1,1],
			 yreversed=true,
			 limits=(-0.5, 1.5, 16, 23))
	y = log10.(cmd_sat.density ./ maximum(cmd_sat.density))
	
	p = heatmap!(cmd_sat.x, cmd_sat.y, y,colorrange=(-3, 0))


	scatter!(stars_cen.gmag .- stars_cen.rmag, stars_cen.gmag, markersize=2, color=COLORS[3])
	fig
end

# ╔═╡ eaa3ee0b-6b8f-4da5-a230-596613bbfa76
lines(iso_new.gmag .- iso_new.rmag, iso_new.gmag .+ distance_modulus)

# ╔═╡ 71dfe223-539a-4e2b-84b1-04811b2533a0
color_range = (-0.8, 1.8)

# ╔═╡ ff619ab8-3481-46bb-98b9-c706f0b1c507
dens_back = Utils.make_background_density(stars_bkg, cmd_sat.x, cmd_sat.y)

# ╔═╡ 98829624-2eb3-4077-9995-cf45c52bfe81
heatmap(dens_back.x, dens_back.y, dens_back.density)

# ╔═╡ 79bb12ce-2a31-47c6-9265-2c5ec9dffa48
function cmd_axis(gs)
	return Axis(gs,
			   xlabel = "g - r",
			   ylabel = "g",
				yreversed = true,
			   )

end
	

# ╔═╡ 4379bff6-5686-4709-b5b7-e8b4cced8664
let
	fig = Figure()
	ax = cmd_axis(fig[1,1])
		
	heatmap!(cmd_sat.x, cmd_sat.y, log10.(cmd_sat.density ./ dens_back.density), colorrange=(-5, 0))


	fig

end

# ╔═╡ 767b40cd-5533-4fea-8a63-20e4c4d0d038
function make_ll_map(xbins, ybins, dens_sat, dens_back)
	dens_sat = dens_sat ./ sum(dens_sat)
	dens_back = dens_back ./ sum(dens_back)
	itp = Interpolations.interpolate((midpoints(xbins), midpoints(ybins)), 
							   dens_sat ./ dens_back, Interpolations.Gridded(Interpolations.Linear()))	

	return Interpolations.extrapolate(itp, 0)
end

# ╔═╡ 0870f9ab-56aa-4208-964b-32257048f38d
function make_l_sat(xbins, ybins, dens_sat)
	dens_sat = dens_sat ./ sum(dens_sat)
	itp = Interpolations.interpolate((midpoints(xbins), midpoints(ybins)), 
							   dens_sat, Interpolations.Gridded(Interpolations.Linear()))	

	return Interpolations.extrapolate(itp, 0)
end

# ╔═╡ 43feb35a-156f-4f84-ae33-47ed3f4ad67e
function make_l_bkg(xbins, ybins, dens_bkg)
	dens_bkg = dens_bkg ./ sum(dens_bkg)
	itp = Interpolations.interpolate((midpoints(xbins), midpoints(ybins)), 
							   dens_bkg, Interpolations.Gridded(Interpolations.Linear()))	

	return Interpolations.extrapolate(itp, 0)
end

# ╔═╡ a576f974-d294-43de-a04a-6254b65343e5
l_sat = Utils.interpolate_density(cmd_sat)

# ╔═╡ 5353e8a6-abe0-431a-8dae-422b84430535
l_bkg = Utils.interpolate_density(dens_back)

# ╔═╡ 0590d501-ff35-4a52-862b-e023143ea911
md"""
# Next stage
"""

# ╔═╡ 3ba9ff5f-9552-4a38-97a0-1cbf0b239012
import StatsBase: median

# ╔═╡ 66a21394-5edc-45ea-9b91-b2dcf618a880
function scatter_matched_filter(gs, LLR, thresh)
	sample = stars[LLR .> thresh, :]
	ax = Axis(gs, 
			 xreversed=true,
			 aspect=DataAspect())

	scatter!(sample.xi, sample.eta, markersize=0.5, color=:black, alpha=0.2)
	
	ax
end

# ╔═╡ 3ff9423f-b41b-4bdb-9f61-0c7a214029db
let
	fig = Figure(size=(4*72, 8*72))

	threshs = [-0.4, -0.2, 0, 0.2, 0.4, 0.6]
	cols = 2

	for (i, thresh) in enumerate(threshs)
		jj = (i-1)%cols + 1
		ii = (i+cols - 1) ÷ 2
		println(ii, ", ", jj)

		cmd_sat = Utils.build_cmd_filter(
			iso_new.gmag .- iso_new.rmag ,
			iso_new.gmag .+ distance_modulus .+ thresh,
			iso_new.Mini,
			x -> @.(sqrt(color_err.(x)^2 + 0.05^2)),
			x -> 0.05
		)
		l_sat = Utils.interpolate_density(cmd_sat)

		LLR = log10.(l_sat.(stars.gmag .- stars.rmag, stars.gmag) ./ l_bkg.(stars.gmag .- stars.rmag, stars.gmag))
		
		ax = scatter_matched_filter(fig[ii, jj], LLR, 0.0)
		ax.title = "dDM = $thresh"
		if ii < length(threshs) ÷ cols
			hidexdecorations!()
		end
		if jj > 1
			hideydecorations!()
		end
	end

	linkaxes!(fig.content...)
	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 5)

	fig

end

# ╔═╡ e91c226a-87e9-4481-b6dc-e3a1741ffe1f
let
	fig = Figure(size=(4*72, 5*72))

	threshs = [10, 12]
	cols = 2

	for (i, thresh) in enumerate(threshs)
		jj = (i-1)%cols + 1
		ii = (i+cols - 1) ÷ 2
		println(ii, ", ", jj)
		iso = Utils.get_isochrone(-2.1, thresh)

		iso_new = resample_iso(iso)
		cmd_sat = Utils.build_cmd_filter(
			iso_new.gmag .- iso_new.rmag ,
			iso_new.gmag .+ distance_modulus,
			iso_new.Mini,
			x -> @.(sqrt(color_err.(x)^2 + 0.05^2)),
			x -> 0.19
		)
		l_sat = Utils.interpolate_density(cmd_sat)

		LLR = log10.(l_sat.(stars.gmag .- stars.rmag, stars.gmag) ./ l_bkg.(stars.gmag .- stars.rmag, stars.gmag))
		
		ax = scatter_matched_filter(fig[ii, jj], LLR, 0.0)
		ax.title = "age = $thresh"
		if ii < length(threshs) ÷ cols
			hidexdecorations!()
		end
		if jj > 1
			hideydecorations!()
		end
	end

	linkaxes!(fig.content...)
	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 5)

	fig

end

# ╔═╡ 596f35c2-2953-4e41-92d2-0b8cbe01435f
let
	fig = Figure(size=(4*72, 8*72))

	threshs = [-2.5, -2, -1.5, -1.0, -0.5]
	cols = 2

	for (i, thresh) in enumerate(threshs)
		jj = (i-1)%cols + 1
		ii = (i+cols - 1) ÷ 2
		println(ii, ", ", jj)
		iso = Utils.get_isochrone(thresh, 12)

		iso_new = resample_iso(iso)
		cmd_sat = Utils.build_cmd_filter(
			iso_new.gmag .- iso_new.rmag ,
			iso_new.gmag .+ distance_modulus,
			iso_new.Mini,
			x -> @.(sqrt(color_err.(x)^2 + 0.05^2)),
			x -> 0.19
		)
		l_sat = Utils.interpolate_density(cmd_sat)

		LLR = log10.(l_sat.(stars.gmag .- stars.rmag, stars.gmag) ./ l_bkg.(stars.gmag .- stars.rmag, stars.gmag))
		
		ax = scatter_matched_filter(fig[ii, jj], LLR, 0.0)
		ax.title = "fe/h = $thresh"
		if ii < length(threshs) ÷ cols
			hidexdecorations!()
		end
		if jj > 1
			hideydecorations!()
		end
	end

	linkaxes!(fig.content...)
	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 5)

	fig

end

# ╔═╡ 0b46cd2f-54db-44e5-901a-e5f2e6e64fa4
let
	fig = Figure(size=(4*72, 5*72))

	threshs = [0, 0.025, 0.05, 0.075, 0.1, 0.2]
	cols = 2

	for (i, thresh) in enumerate(threshs)
		jj = (i-1)%cols + 1
		ii = (i+cols - 1) ÷ 2
		println(ii, ", ", jj)

		cmd_sat = Utils.build_cmd_filter(
			iso_new.gmag .- iso_new.rmag ,
			iso_new.gmag .+ distance_modulus,
			iso_new.Mini,
			x -> @.(sqrt(color_err.(x)^2 + thresh^2)),
			x -> 0.19
		)
		l_sat = Utils.interpolate_density(cmd_sat)

		LLR = log10.(l_sat.(stars.gmag .- stars.rmag, stars.gmag) ./ l_bkg.(stars.gmag .- stars.rmag, stars.gmag))
		
		ax = scatter_matched_filter(fig[ii, jj], LLR, 0.0)
		ax.title = "iso width = $thresh"
		if ii < length(threshs) ÷ cols
			hidexdecorations!()
		end
		if jj > 1
			hideydecorations!()
		end
	end

	linkaxes!(fig.content...)
	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 5)

	fig

end

# ╔═╡ 697b507f-14b9-47bb-9463-05f760c0f2f0
tangent_axis(gs; kwargs...) = Axis(gs,
	xlabel = L"$\xi$ / arcmin",
	ylabel = L"$\eta$ / arcmin",
	xreversed = true,
	aspect=DataAspect(),
								   
)

# ╔═╡ 6b608dcb-c853-4457-be72-60060786bb3c
R_max = maximum(stars.R)

# ╔═╡ f61498bf-cc1f-48e1-981c-2e3b7ac88ede
df_out = let
	df = copy(stars)
	df[!, :L_sat] = l_sat.(stars.gmag .- stars.rmag, stars.gmag)
	df[!, :L_bg] = l_bkg.(stars.gmag .- stars.rmag, stars.gmag)
	df
end

# ╔═╡ eef36a62-34ce-4d46-9fce-6810cc430042
LLR = log10.(df_out.L_sat ./ df_out.L_bg)

# ╔═╡ b1857714-b3a4-4a38-a816-d288ba950dbd
let
	fig = Figure()
	ax = cmd_axis(fig[1,1])
		
	p = scatter!(stars.gmag .- stars.rmag, stars.gmag, color=LLR, markersize=1, 
				colorrange=(-2, 2))

	Colorbar(fig[1,2], p, label="LLR", minorticksvisible=false)

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
	

# ╔═╡ 0942d5fd-f1b0-456a-883d-0d74494353af
sample = stars[LLR .> 0.5, :]

# ╔═╡ 54510ba4-4ecb-4c77-9437-2c200e3e8a60
scatter(sample.xi, sample.eta, markersize=0.3, color=:black)

# ╔═╡ 3169bac1-a7e9-4c35-9a82-0798316b262a
let
	fig = Figure(size=(4*72, 8*72))

	threshs = [-2, -1, 0, 0.5, 1, 1.5]
	cols = 2

	for (i, thresh) in enumerate(threshs)
		jj = (i-1)%cols + 1
		ii = (i+cols - 1) ÷ 2
		println(ii, ", ", jj)
		
		ax = scatter_matched_filter(fig[ii, jj], LLR, thresh)
		ax.title = "LLR > $thresh"
		if ii < length(threshs) ÷ cols
			hidexdecorations!()
		end
		if jj > 1
			hideydecorations!()
		end
	end

	linkaxes!(fig.content...)
	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 5)

	fig

end

# ╔═╡ 11f6e923-b608-43a8-bd7a-7b8b668549c0
function plot_matched_filter(thresh, q=0.99; weighted=false, bw=10)
	sample = stars[LLR .> thresh, :]
	if weighted 
		ws = 10 .^ (LLR[LLR .> thresh])
	else
		ws = KernelDensity.UniformWeights(1)
	end
	fig = Figure(size=(6, 5) .* 72)
	ax = tangent_axis(fig[1,1])
	xlims!(-R_max, R_max)
	ylims!(-R_max, R_max)

	k = KernelDensity.kde((sample.xi, sample.eta), weights = ws, bandwidth=(bw, bw))

	
	ax.xreversed[] = true

	c0 = median(k.density)
	c1 = LilGuys.quantile(vec(k.density), q)
	p = heatmap!(k, colorrange=(0, c1), colormap=:Greys)
	scatter!(0, 0, color=COLORS[1], marker=:+)

	Colorbar(fig[1,2], p)
	
	fig
end

# ╔═╡ cc8c9b6e-72c6-4aec-b650-1ff93b226cf1
plot_matched_filter(1.5, 0.995)

# ╔═╡ 57162c63-30c2-4a73-8b37-f962e83de060
plot_matched_filter(1, 0.995)

# ╔═╡ 87d03d79-38c9-40ff-83f3-28d4af397045
plot_matched_filter(-Inf, 0.98,  weighted=true)


# ╔═╡ 9f32320d-e1eb-4504-8efd-718ab01264eb
plot_matched_filter(1, weighted=true)

# ╔═╡ d4707ec8-e546-4e7b-b669-4793edbbb2ec
plot_matched_filter(1)

# ╔═╡ a1b4f511-2a73-4d5e-88eb-e516ba27c654
plot_matched_filter(0.5)

# ╔═╡ 2b42d10d-3a29-4c79-a44a-295f563617c5
plot_matched_filter(-0.5, 0.99)


# ╔═╡ a2b53ef1-d50d-42bf-bb3e-a1944d53aec1
plot_matched_filter(-1, 0.99)


# ╔═╡ ec22ea19-f66c-46f7-9973-d63a81fbda1a
plot_matched_filter(-2, 0.99)


# ╔═╡ abf77f38-3edf-4433-9676-7cdf46b001d8
plot_matched_filter(-Inf, 0.994)


# ╔═╡ c8fc5a1d-909a-4948-97db-eb8fc7722560
let
	fig = Figure()
	scatter_matched_filter(fig[1,1], LLR, 0.0 )

	fig

end

# ╔═╡ cf349fec-711c-4c58-86b2-74771e00f4bc
LLR

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
				colorrange=(0, 1e-4))

	Colorbar(fig[1,2], p, label="LR", minorticksvisible=false)

	fig

end
	

# ╔═╡ 49204205-cf52-4489-a2bd-972e1112b0a8
write_fits("samples/delve_matched_filter.fits", df_out, overwrite=true)

# ╔═╡ 674dd04b-b384-469e-80ce-13f3662387fb
md"""
# Density profile
"""

# ╔═╡ 688cb9b4-18b3-4d0e-b046-593b999b48f4
prof = LilGuys.SurfaceDensityProfile(
	stars.R, 
	weights = 10 .^ LLR,
	bins = -2:0.1:2.5
)

# ╔═╡ a13e0740-0ed6-4e83-abd8-756a5e2bfaa6
Σ_bg = median(prof.log_Sigma[prof.log_R .> 2])

# ╔═╡ f83a137b-a4c6-456b-8c0d-d856a78608ff
readdir()

# ╔═╡ 4e795b37-9ca9-421b-9b8f-8b8466921913
open("density_profiles/delve_mf_profile.toml", "w") do f
	print(f, prof)
end

# ╔═╡ 81a0004d-3764-4f8f-8845-7908b379fc28

let

	fig = Figure()
	ax = Axis(fig[1,1])
		
	LilGuys.plot_log_Σ!(ax, prof)

	hlines!(Σ_bg.middle)
	ylims!(-1, 2)
	fig
end

# ╔═╡ Cell order:
# ╠═c185afc3-f34b-4f75-b7e5-68c7d7e45a3a
# ╠═4c482a6f-9cd3-4704-b8ce-b23e0f1a293f
# ╠═11ce9995-58ac-4642-851b-092ce72c006e
# ╠═015eae5d-515c-4fa4-9a44-fc478692cd44
# ╠═68cd2e80-a8d7-4fcd-87dd-a632818b77ac
# ╠═34c475e5-1eb9-482a-a29b-182caad50313
# ╠═e8c0e80d-336e-4f23-b826-278d90d732d1
# ╠═e3cb612c-7c1f-4688-8aa0-3cc4de19ec9b
# ╠═e560b3a4-249d-11f1-a7d8-918587896934
# ╠═6ed9796f-b522-4e62-82b5-62f0297c1e35
# ╠═6c55e84e-39a9-444b-a067-95a3e3826cbd
# ╠═fdc022f6-f3f8-419c-a42e-4d695f00536d
# ╠═7cf639a4-cd5e-4f42-a0ac-a7cd768debdc
# ╠═109fdcda-b64e-477f-bfba-4f9588b7bf09
# ╠═ef15705f-2390-40d5-b6d7-50064ca333ec
# ╠═3b521894-468e-4cad-bb0f-394567adcbdf
# ╠═eaa3ee0b-6b8f-4da5-a230-596613bbfa76
# ╠═78b5414f-8961-4446-bfb7-9c0e852dd5b7
# ╠═53a1dd27-bcd8-481b-abd7-160da5dc643d
# ╠═f63389bf-77fe-4495-9be5-b4f86f84dc57
# ╠═baccd40f-4128-4031-973b-1389845f2972
# ╠═4e7dacb3-366e-46f4-84f4-fd6206961c4e
# ╠═88a5e8a3-00dd-47b3-ba2f-ab3ffc483c25
# ╠═c47505fc-01f6-4252-815a-2be038d6d3bc
# ╠═0c964ad8-cc76-49cc-922c-12d20e5d259f
# ╠═22934ab9-2750-4826-b916-a51ecee71f68
# ╠═90a3ba08-28eb-4529-bfb7-c6d01b3dd979
# ╠═c94244bd-31f0-4838-8495-05fc8be234e8
# ╠═a75c32b1-9011-4c88-8725-624f93fd0d98
# ╠═a4159d32-8fd0-4577-ab14-8841ae328b52
# ╠═bb0600b7-cf8a-484b-9d4a-7dab2c963da6
# ╠═764be90a-e602-44d8-9956-d407ee2a0d4e
# ╠═39d177ca-6402-4216-b503-c12f83f91d56
# ╠═48993623-e10e-45c7-97df-93e3e0200007
# ╠═bf0e16fe-94cd-41fd-a107-4a814dd448f4
# ╠═080d5e15-dd13-423a-815c-e22edfc73dcc
# ╠═71dfe223-539a-4e2b-84b1-04811b2533a0
# ╠═ff619ab8-3481-46bb-98b9-c706f0b1c507
# ╠═98829624-2eb3-4077-9995-cf45c52bfe81
# ╠═79bb12ce-2a31-47c6-9265-2c5ec9dffa48
# ╠═4379bff6-5686-4709-b5b7-e8b4cced8664
# ╠═767b40cd-5533-4fea-8a63-20e4c4d0d038
# ╠═0870f9ab-56aa-4208-964b-32257048f38d
# ╠═43feb35a-156f-4f84-ae33-47ed3f4ad67e
# ╠═a576f974-d294-43de-a04a-6254b65343e5
# ╠═5353e8a6-abe0-431a-8dae-422b84430535
# ╠═eef36a62-34ce-4d46-9fce-6810cc430042
# ╠═b1857714-b3a4-4a38-a816-d288ba950dbd
# ╠═9e3ce97d-67f5-4a8d-b36d-d20e051715ad
# ╠═17522c59-622b-4fe3-88af-f48308f3137b
# ╠═f884701a-80e7-458c-88d7-db626254608e
# ╟─0590d501-ff35-4a52-862b-e023143ea911
# ╠═0942d5fd-f1b0-456a-883d-0d74494353af
# ╠═54510ba4-4ecb-4c77-9437-2c200e3e8a60
# ╠═3ba9ff5f-9552-4a38-97a0-1cbf0b239012
# ╠═66a21394-5edc-45ea-9b91-b2dcf618a880
# ╠═3ff9423f-b41b-4bdb-9f61-0c7a214029db
# ╠═e91c226a-87e9-4481-b6dc-e3a1741ffe1f
# ╠═596f35c2-2953-4e41-92d2-0b8cbe01435f
# ╠═0b46cd2f-54db-44e5-901a-e5f2e6e64fa4
# ╠═3169bac1-a7e9-4c35-9a82-0798316b262a
# ╠═697b507f-14b9-47bb-9463-05f760c0f2f0
# ╠═6b608dcb-c853-4457-be72-60060786bb3c
# ╠═11f6e923-b608-43a8-bd7a-7b8b668549c0
# ╠═cc8c9b6e-72c6-4aec-b650-1ff93b226cf1
# ╠═57162c63-30c2-4a73-8b37-f962e83de060
# ╠═87d03d79-38c9-40ff-83f3-28d4af397045
# ╠═9f32320d-e1eb-4504-8efd-718ab01264eb
# ╠═d4707ec8-e546-4e7b-b669-4793edbbb2ec
# ╠═a1b4f511-2a73-4d5e-88eb-e516ba27c654
# ╠═c8fc5a1d-909a-4948-97db-eb8fc7722560
# ╠═2b42d10d-3a29-4c79-a44a-295f563617c5
# ╠═a2b53ef1-d50d-42bf-bb3e-a1944d53aec1
# ╠═ec22ea19-f66c-46f7-9973-d63a81fbda1a
# ╠═abf77f38-3edf-4433-9676-7cdf46b001d8
# ╠═cf349fec-711c-4c58-86b2-74771e00f4bc
# ╠═f61498bf-cc1f-48e1-981c-2e3b7ac88ede
# ╠═49204205-cf52-4489-a2bd-972e1112b0a8
# ╠═674dd04b-b384-469e-80ce-13f3662387fb
# ╠═688cb9b4-18b3-4d0e-b046-593b999b48f4
# ╠═a13e0740-0ed6-4e83-abd8-756a5e2bfaa6
# ╠═f83a137b-a4c6-456b-8c0d-d856a78608ff
# ╠═4e795b37-9ca9-421b-9b8f-8b8466921913
# ╠═81a0004d-3764-4f8f-8845-7908b379fc28
