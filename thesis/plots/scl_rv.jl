### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ 5df6928e-59aa-4f38-ba3a-6428855934ec
using PyFITS

# ╔═╡ 6562dd6e-397b-4632-bc44-f227201d9012
using KernelDensity

# ╔═╡ 880b78c0-574f-454e-ad81-4525aa8d9713
include("./paper_style.jl")

# ╔═╡ c6631ad6-cc41-4a67-82c9-9d2b20115c9c
include("utils.jl")

# ╔═╡ 88025a55-f541-4aa1-8426-14e83bd790b7
CairoMakie.activate!(type=:png)

# ╔═╡ 63d5e3f6-6cec-4e50-8fa1-8907508aa47f
import TOML

# ╔═╡ 4c551542-92dc-4a61-b4aa-e272fd29d565
stars = read_fits(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/rv_members_all.fits")

# ╔═╡ b2f9ed4b-52be-4eeb-a419-d6406148e6c3
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", "sculptor", "observed_properties.toml"))

# ╔═╡ 1f07b964-a443-4f70-abf5-35512f9d6b69
gsr = LilGuys.transform(GSR, ICRS(obs_props))

# ╔═╡ 23989a87-6aad-4ca1-bd49-f6f13977234b
v0 = gsr.radial_velocity

# ╔═╡ 87c47d8f-4639-49ca-a4e7-467f583c2e3e
σv =  obs_props["sigma_v"]

# ╔═╡ 0e7ed733-f30b-4e5d-88b3-211b454125b3
extrema(stars.vz)

# ╔═╡ 88fbbe64-e5a1-489e-9cf3-8882bc650417
hist(stars.vz, bins=60)

# ╔═╡ 6343f079-f595-4969-889d-e734c7ad4878
import StatsBase: weights

# ╔═╡ bcff6bcf-f59c-47e9-9dd2-40f4dcaa4d1d
r_max = 60

# ╔═╡ 1d6d4397-2f2e-43d2-9761-3100b2a31018
outside_bins = .!( (-r_max .< stars.xi .< r_max) .& (-r_max .< stars.eta .< r_max))

# ╔═╡ 1523b300-076f-48e3-8c65-1a493a17433b
skipmissing(stars.RV_gmos) |> collect

# ╔═╡ 5e44fbcc-aa87-4ccf-a149-8d75575cb50d
df_gmos = stars[.!ismissing.(stars.RV_gmos), :]

# ╔═╡ 21c9e84e-80d2-4122-b541-3c5f283ae6a3
sort(abs.(stars.xi), order=Base.Order.Reverse)

# ╔═╡ d0f91075-c9db-40f3-94f1-275273556879
sort(abs.(stars.eta), order=Base.Order.Reverse)

# ╔═╡ ace4143d-a8ff-40f4-ae16-76cce672d860
df_gmos.xi

# ╔═╡ 762ab1b8-7cbe-4b37-8e81-cb3024bc5f40
tangent_bins = -r_max:6:r_max

# ╔═╡ 4a17d821-ae04-45e9-ab37-0a4962f4a153
pm_gradient = LilGuys.pm2kms.([gsr.pmra, gsr.pmdec], gsr.distance) ./ (180/π)

# ╔═╡ 528f8b9a-67b9-4010-ab04-0c6f8d7202f2
derived_gradient = 5 * [sind(-149), cosd(-149)]

# ╔═╡ a4d317d9-7816-4d2a-ae74-748dc7c1c23b
scatter(stars.xi, stars.eta, axis=(;
		xreversed=true,
		aspect=DataAspect(),
		))

# ╔═╡ 5fce21b3-a4bb-4e47-a4e8-7a02c7590a9f
let
	fig, ax = FigAxis(
		xlabel = L"\xi\ /\ \textrm{arcmin}",
		ylabel = L"\eta\ /\ \textrm{arcmin}",
		aspect=DataAspect(),
		xreversed=true,
		limits=(-70, 70, -60, 60)
	)

	n_min = 3
	memb_stars = stars


	dv = memb_stars.vz .- v0
	x = memb_stars.xi
	y = memb_stars.eta
	w = 1 ./ memb_stars.vz_err .^ 2

	colorrange=(-3σv, 3σv)

	k1 = Arya.histogram2d(x, y, tangent_bins, weights= w .* dv)
	k2 = Arya.histogram2d(x, y, tangent_bins, weights=w)
	k3 = Arya.histogram2d(x, y, tangent_bins)

	k1.values ./= k2.values
	k1.values[k3.values .< n_min] .= NaN


	scatter!(x, y,
		color = dv,
		colormap=:bluesreds,
		colorrange=colorrange,
		markersize=5.
	)


	p = heatmap!(k1.xbins, k1.ybins, k1.values, 
		colormap=:bluesreds,
		colorrange=colorrange,
		rasterize=true,
		)

	arrows!([0], [0], [3pm_gradient[1]], 3pm_gradient[2], linewidth=2, arrowsize=8)
	# arrows!([0], [0], [5derived_gradient[1]], 5derived_gradient[2])
	
	Colorbar(fig[1, 2], p, label=L"$(v_{z} - \bar{v}_z)$ / km/s",)
	ellipse!(1obs_props["r_h"], obs_props["ellipticity"], obs_props["position_angle"], color=:black)
	text!(-15, 0, text=L"R_h")

	@savefig "scl_rv_2dhist"
	fig
end

# ╔═╡ 0b19e9e2-e7f4-427b-a784-52401f5b94ff
let
	fig, ax = FigAxis(
		xlabel = L"\xi\ /\ \textrm{arcmin}",
		ylabel = L"\eta\ /\ \textrm{arcmin}",
		aspect=DataAspect(),
		xreversed=true,
		limits=(-70, 70, -60, 60)
	)

	n_min = 3
	memb_stars = stars


	dv = memb_stars.vz .- v0
	x = memb_stars.xi
	y = memb_stars.eta

	colorrange=(-3σv, 3σv)

	k1 = Arya.histogram2d(x, y, tangent_bins, weights= dv)
	k3 = Arya.histogram2d(x, y, tangent_bins)

	k1.values ./= k3.values
	k1.values[k3.values .< n_min] .= NaN

	scatter!(x, y,
		color = dv,
		colormap=:bluesreds,
		colorrange=colorrange,
		markersize=3.
	)


	p = heatmap!(k1.xbins, k1.ybins, k1.values, 
		colormap=:bluesreds,
		colorrange=colorrange,
		rasterize=true,
		)

	Colorbar(fig[1, 2], p, label=L"$v_{z} - \bar{v}_z$ / km/s",
)
	fig
end

# ╔═╡ 0d32dc0e-ab0e-46d9-9f51-0684c41d97d1
import DensityEstimators as DE

# ╔═╡ 4a4b77de-f1dd-408f-9707-06f62420db8f
let
	fig, ax = FigAxis(
		xlabel = L"\xi\ /\ \textrm{arcmin}",
		ylabel = L"\eta\ /\ \textrm{arcmin}",
		aspect=DataAspect(),
		xreversed=true,
		limits=(-70, 70, -60, 60)
	)

	density_min = 3e-5
	memb_stars = stars


	dv = memb_stars.vz .- v0
	x = memb_stars.xi
	y = memb_stars.eta
	w = 1 ./ memb_stars.vz_err .^ 0

	colorrange=(-3σv, 3σv)
	bandwidth=3

	bins = 100
	xlims = (-60, 60)
	ylims = (-60, 60)
	k1 = DE.kde2d(x, y, bandwidth, weights= w .* dv, bins=bins, limits=(xlims, ylims))
	k2 = DE.kde2d(x, y,bandwidth,  weights=w,  bins=bins, limits=(xlims, ylims))
	k3 = DE.kde2d(x, y,bandwidth,  bins=bins, limits=(xlims, ylims) )

	density = k1.values ./ k2.values
	density[k3.values .< density_min] .= NaN

	scatter!(x, y,
		color = dv,
		colormap=:bluesreds,
		colorrange=colorrange,
		markersize=5.
	)


	p = heatmap!(k1.x, k1.y, density, 
		colormap=:bluesreds,
		colorrange=colorrange,
		rasterize=true,
		)

	Colorbar(fig[1, 2], p, label=L"$v_{z} - \bar{v}_z$ / km/s",
)
	fig
end

# ╔═╡ 16ea0587-7f70-405c-a556-0f1b2991ba5b
let
	fig, ax = FigAxis(
		xlabel = L"\xi\ /\ \textrm{arcmin}",
		ylabel = L"\eta\ /\ \textrm{arcmin}",
		aspect=DataAspect(),
		xreversed=true,
		limits=(-70, 70, -60, 60)
	)

	density_min = 1e-5
	memb_stars = stars


	dv = memb_stars.vz .- v0
	x = memb_stars.xi
	y = memb_stars.eta
	w = 1 ./ memb_stars.vz_err .^ 0

	colorrange=(-3σv, 3σv)
	bandwidth=DE.bandwidth_knn(x, y; η=3, k=3)

	bins = 100
	xlims = (-70, 70)
	ylims = (-60, 60)
	k1 = DE.kde2d(x, y, bandwidth, weights= w .* dv, bins=bins, limits=(xlims, ylims))
	k2 = DE.kde2d(x, y,bandwidth,  weights=w, bins=bins, limits=(xlims, ylims))
	k3 = DE.kde2d(x, y,bandwidth, bins=bins,limits=(xlims, ylims) )

	density = k1.values ./ k2.values
	density[k3.values .< density_min] .= NaN


	p = heatmap!(k1.x, k1.y, density, 
		colormap=:bluesreds,
		colorrange=colorrange,
		rasterize=true,
		)

	# scatter!(x, y, color=dv, strokecolor=:grey, strokewidth=0.2, markersize=1,
	# 		colorrange=colorrange, colormap=:bluesreds
	# 		)

	Colorbar(fig[1, 2], p, label=L"$v_{z} - \bar{v}_z$ / km/s",
)
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═88025a55-f541-4aa1-8426-14e83bd790b7
# ╠═5df6928e-59aa-4f38-ba3a-6428855934ec
# ╠═63d5e3f6-6cec-4e50-8fa1-8907508aa47f
# ╠═880b78c0-574f-454e-ad81-4525aa8d9713
# ╠═c6631ad6-cc41-4a67-82c9-9d2b20115c9c
# ╠═4c551542-92dc-4a61-b4aa-e272fd29d565
# ╠═b2f9ed4b-52be-4eeb-a419-d6406148e6c3
# ╠═1f07b964-a443-4f70-abf5-35512f9d6b69
# ╠═23989a87-6aad-4ca1-bd49-f6f13977234b
# ╠═87c47d8f-4639-49ca-a4e7-467f583c2e3e
# ╠═0e7ed733-f30b-4e5d-88b3-211b454125b3
# ╠═88fbbe64-e5a1-489e-9cf3-8882bc650417
# ╠═6343f079-f595-4969-889d-e734c7ad4878
# ╠═bcff6bcf-f59c-47e9-9dd2-40f4dcaa4d1d
# ╠═1d6d4397-2f2e-43d2-9761-3100b2a31018
# ╠═1523b300-076f-48e3-8c65-1a493a17433b
# ╠═5e44fbcc-aa87-4ccf-a149-8d75575cb50d
# ╠═21c9e84e-80d2-4122-b541-3c5f283ae6a3
# ╠═d0f91075-c9db-40f3-94f1-275273556879
# ╠═ace4143d-a8ff-40f4-ae16-76cce672d860
# ╠═762ab1b8-7cbe-4b37-8e81-cb3024bc5f40
# ╠═4a17d821-ae04-45e9-ab37-0a4962f4a153
# ╠═528f8b9a-67b9-4010-ab04-0c6f8d7202f2
# ╠═a4d317d9-7816-4d2a-ae74-748dc7c1c23b
# ╠═5fce21b3-a4bb-4e47-a4e8-7a02c7590a9f
# ╠═0b19e9e2-e7f4-427b-a784-52401f5b94ff
# ╠═6562dd6e-397b-4632-bc44-f227201d9012
# ╠═4a4b77de-f1dd-408f-9707-06f62420db8f
# ╠═16ea0587-7f70-405c-a556-0f1b2991ba5b
# ╠═0d32dc0e-ab0e-46d9-9f51-0684c41d97d1
