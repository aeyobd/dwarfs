### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ bd5863e6-3db7-11f1-9b03-57d4b01f128c
begin
	import Pkg; Pkg.activate()


	using LilGuys
	using Arya, CairoMakie
	
	using PyFITS	
	import TOML
	using DensityEstimators
	using OrderedCollections
end

# ╔═╡ f868f5a0-91fa-49f2-9a71-033f2718cc4d
CairoMakie.activate!(type=:png)

# ╔═╡ e1305938-5d93-4776-85e2-a33ddde7a787
md"""
# Data Loading
"""

# ╔═╡ b2ee468f-e3e1-4973-a2ae-82829f39cb09
function load_final_stars(modelname, starsname)
	df = read_fits(joinpath(ENV["DWARFS_ROOT"], "analysis", "bootes3", modelname, "stars", starsname, "final.fits"))

	row = df[1, :]
	icrs_cen = ICRS(ra=row.ra, dec=row.dec, distance=row.distance, 
				   pmra = row.pmra, pmdec=row.pmdec, radial_velocity=row.radial_velocity)
	@assert df.index[1] == 0
	gsr_cen = LilGuys.transform(GSR, icrs_cen)
	xi_p, eta_p = LilGuys.to_orbit_coords(df.ra, df.dec, gsr_cen.ra, gsr_cen.dec,
										 atand(gsr_cen.pmra, gsr_cen.pmdec))

	df[!, :xi_p] = xi_p
	df[!, :eta_p] = eta_p
	df
end

# ╔═╡ a3641424-a6d3-485a-a61b-52dad7693130
modelnames = OrderedDict(
	"fiducial" => ["bootes3", "1e6_v30_r2.2/orbit_12kpc", "exp2d_rs0.20"],
	"fiducial_big" => ["bootes3", "1e6_v30_r2.2/orbit_12kpc", "exp2d_rs0.40"],
	"1x12kpc" => ["bootes3", "1e6_v22_r3.9/1_peri_12kpc", "exp2d_rs0.20"],
	"1x1.5kpc" => ["bootes3", "1e6_v30_r3.0/1_peri_1.5kpc", "exp2d_rs0.20"],
	"2x7kpc" => ["bootes3", "1e6_v30_r3.0/2_peri_7kpc", "exp2d_rs0.20"],
	"3x26kpc" => ["bootes3", "1e6_v22_r3.9/3_peri_26kpc", "exp2d_rs0.20"],
	# "5x18kpc" => ["bootes3", "1e6_v30_r3.0/5_peri_18kpc", "exp2d_rs0.20"],
	# "5x18kpc-P" => ["bootes3", "1e6_v30_r3.0/5_peri_18kpc", "plummer_rs0.30"],
)

# ╔═╡ 44e7711f-65b3-4cc5-9f66-98bb5e4bbb8c
allstars = Dict(
	label => load_final_stars(model[2], model[3]) for (label, model) in modelnames
)

# ╔═╡ a6a87544-5fdc-4118-8948-2ce0baa7724b
md"""
# Plots
"""

# ╔═╡ 9c31d05f-8c6e-4d5a-864e-83b640d13a5b
md"""
# Density Profiles...
"""

# ╔═╡ a5b0c806-c46d-4d89-b709-334138c28278
let 
	f = scatter(
		allstars["fiducial"].xi_p, allstars["fiducial"].eta_p, 
		markersize=1, alpha=0.02, color=:black,
	    axis=(; aspect=DataAspect())
	)
	
	hlines!(-1, color=COLORS[2])
	hlines!(1, color=COLORS[2])

	f

end

# ╔═╡ d2850fcc-746e-4bed-aa09-9edf872f807c
1500 / 60

# ╔═╡ 16b6cfdd-6851-475e-96e0-1da503d4ffa9
Ntot = 100

# ╔═╡ 01a9451b-c901-40ff-b25d-965b919138ef
function plot_stars_hist!(stars; 
						  limits=(-5, 5, -5, 5), bins=nothing, colorrange=nothing, kwargs...)
	w_scale = Ntot / sum(stars.weights)

	h = Arya.histogram2d(stars.xi * 60, stars.eta * 60, bins,  weights=stars.weights * w_scale,  normalization=:density, limits=limits)

	y = log10.(h.values)


	if colorrange isa Real
		colorrange = (maximum(y) - colorrange, maximum(y))
	end

	if isnothing(colorrange)
		colorrange = extrema(y)
	end
	if isnothing(colorrange[1])
		colorrange = (minimum(y), colorrange[2])
	end
	if isnothing(colorrange[2])
		colorrange = (colorrange[1]*maximum(y), maximum(y))
	end

	@info colorrange
	
	heatmap!(midpoints(h.xbins), midpoints(h.ybins), y; colorrange=colorrange, kwargs...)
				
end

# ╔═╡ 8466fec3-c787-4859-aaf6-8655fc5f8906
function plot_stars_hist(stars; limits=20 .* 60 .* (-1, 1, -1, 1), title="")
	fig = Figure(size=(4, 4) .* 72)
	ax = Axis(fig[1,1], 
			 aspect = DataAspect(), 
			 xlabel = "xi / deg", 
			 ylabel = "eta / deg", 
			 xreversed = true, 
			  title = title
			 )

	p = plot_stars_hist!(stars; 
						 limits=limits, colormap=:greys, bins=200, colorrange=5)

	Colorbar(fig[1,2], p, label="log stellar density")
	fig

end

# ╔═╡ 04cc08ae-154a-4f51-87d0-fcf59902502e
plot_stars_hist(allstars["fiducial"], title="fiducial")

# ╔═╡ b0f9f918-6e80-4d6a-9e7d-4e7500d26cdc
plot_stars_hist(allstars["fiducial_big"], title="fiducial")

# ╔═╡ 7436f81b-52b0-4848-9aec-7af5a29e4b23
plot_stars_hist(allstars["1x1.5kpc"], title="1x1.5kpc")

# ╔═╡ a74835c4-0128-4219-b72d-6637c9a3b5c0
plot_stars_hist(allstars["1x12kpc"], title="1x12kpc")

# ╔═╡ 23de9dfc-77d1-42c6-9438-6d584fc1d5b6
plot_stars_hist(allstars["2x7kpc"], title="2x7kpc")

# ╔═╡ bd6ff7e3-5c4f-427a-aae8-2a8661e824af
plot_stars_hist(allstars["3x26kpc"], title="3x26kpc")

# ╔═╡ 85f528f3-f79c-4c51-98f3-c90f55393be3
function along_axis_density_profile(stars; bins = -20:0.3:20, eta_max=1)
	w_scale = Ntot / sum(stars.weights)
	
	filt = abs.(stars.eta_p) .< eta_max
	_, h, _ = LilGuys.histogram(stars.xi_p[filt], bins, weights=stars.weights[filt] * w_scale, errors=:weighted)

	A = diff(bins) .* 2*eta_max .* 60^2
	return bins .* 60, h ./ A
end

# ╔═╡ 7963bb3a-53d8-424c-b013-ad8268f9dee5
function along_axis_density_profile_unweighted(stars; bins = -20:0.3:20, eta_max=1)
	
	filt = abs.(stars.eta_p) .< eta_max
	_, h, e = LilGuys.histogram(stars.xi_p[filt], bins)

	A = diff(bins) .* 2*eta_max .* 60^2
	return bins .* 60, h ./ A, e ./ A
end

# ╔═╡ 0a037d03-1f3d-4eaa-b569-fb384453d115
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/bootes3/observed_properties.toml"))

# ╔═╡ b7868c01-8901-4922-94ac-cc8400c43542
j24 = let
	df = read_fits(joinpath(ENV["DWARFS_ROOT"], "observations/bootes3/samples/jax_LLR_1_sample.fits"))

	gsr_cen = LilGuys.transform(GSR, ICRS(obs_props))
	xi_p, eta_p = LilGuys.to_orbit_coords(df.ra, df.dec, gsr_cen.ra, gsr_cen.dec,
										 atand(gsr_cen.pmra, gsr_cen.pmdec))

	df[!, :xi_p] = xi_p
	df[!, :eta_p] = eta_p
	df[!, :weights] .= 1
	df
end

# ╔═╡ 72213f55-c60f-4213-89e8-645f21c7ff89
let
	fig = Figure(size=(3.3, 2) .* 72)
	ax = Axis(
		fig[1,1], 
		xlabel = "along stream distance / deg", 
		ylabel = L"$\log\ \Sigma_\star$ / deg^2"
	)




	for (label, model) in modelnames
		stars = allstars[label]
	
		bins, Σ = along_axis_density_profile(stars)
	
		lines!(midpoints(bins) / 60, log10.(Σ ), label=label)
	end

	bins, Σ, Σ_e = along_axis_density_profile_unweighted(j24, bins=-4.5:1:4.5)

	Σ = LilGuys.Measurement.(Σ, Σ_e)
	log_Σ = log10.(Σ )
	errorscatter!(midpoints(bins) / 60, middle.(log_Σ), yerror=error_interval.(log_Σ), color=:black)

	ylims!(-5, 2)

	Legend(fig[1,2], ax)

	fig
end

# ╔═╡ Cell order:
# ╠═bd5863e6-3db7-11f1-9b03-57d4b01f128c
# ╠═f868f5a0-91fa-49f2-9a71-033f2718cc4d
# ╠═e1305938-5d93-4776-85e2-a33ddde7a787
# ╠═b2ee468f-e3e1-4973-a2ae-82829f39cb09
# ╠═a3641424-a6d3-485a-a61b-52dad7693130
# ╠═44e7711f-65b3-4cc5-9f66-98bb5e4bbb8c
# ╟─a6a87544-5fdc-4118-8948-2ce0baa7724b
# ╠═8466fec3-c787-4859-aaf6-8655fc5f8906
# ╠═01a9451b-c901-40ff-b25d-965b919138ef
# ╠═04cc08ae-154a-4f51-87d0-fcf59902502e
# ╠═b0f9f918-6e80-4d6a-9e7d-4e7500d26cdc
# ╠═7436f81b-52b0-4848-9aec-7af5a29e4b23
# ╠═a74835c4-0128-4219-b72d-6637c9a3b5c0
# ╠═23de9dfc-77d1-42c6-9438-6d584fc1d5b6
# ╠═bd6ff7e3-5c4f-427a-aae8-2a8661e824af
# ╟─9c31d05f-8c6e-4d5a-864e-83b640d13a5b
# ╠═a5b0c806-c46d-4d89-b709-334138c28278
# ╠═d2850fcc-746e-4bed-aa09-9edf872f807c
# ╠═16b6cfdd-6851-475e-96e0-1da503d4ffa9
# ╠═85f528f3-f79c-4c51-98f3-c90f55393be3
# ╠═7963bb3a-53d8-424c-b013-ad8268f9dee5
# ╠═b7868c01-8901-4922-94ac-cc8400c43542
# ╠═0a037d03-1f3d-4eaa-b569-fb384453d115
# ╠═72213f55-c60f-4213-89e8-645f21c7ff89
