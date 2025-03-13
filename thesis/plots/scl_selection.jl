### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ bff50014-bfa9-11ee-33f0-0f67e543c2d4
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames 
	using CSV
	using CairoMakie
end

# ╔═╡ ae29bed0-6700-47f1-8952-35e867ce126b
using OrderedCollections

# ╔═╡ 1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
begin 
	import LilGuys as lguys
	using LilGuys
	FIGDIR = "figures"
end

# ╔═╡ 69c98029-165c-407b-9a63-a27e06e30e45
include("paper_style.jl")

# ╔═╡ 47b8b3b0-0228-4f50-9da4-37d388ef9e9f
md"""
# Jensen et al. 2024 sammple

Some plots to understand the (unmodified) J+24 data sample.
The goals here are to investigate the main density profile, likelihoods, and what the cuts appear as in CMD space. 

"""

# ╔═╡ eca9c1c5-e984-42d4-8854-b227fdec0a8a
galaxyname = "sculptor"

# ╔═╡ 0004f638-c57a-4dab-8b97-c77840cafbbf
import TOML

# ╔═╡ cf7aeb0a-d453-462a-b43d-a832567440fd
diverging_cmap = Reverse(:bluesreds)

# ╔═╡ 2eb4aa78-0fea-460b-a18e-06a129c41504
md"""
# Inputs & Data loading
"""

# ╔═╡ 5db5adc4-9b98-4025-9b3c-65b40e5d4c59
obs_dir = ENV["DWARFS_ROOT"] * "/observations/" 

# ╔═╡ 9ecf79a8-2ed3-40c6-b555-a102250ecbd4
observed_properties = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/" * galaxyname * "/observed_properties.toml")

# ╔═╡ 26cf1867-02be-4d36-8c35-6c58a1feca27
datafile = obs_dir * "/$galaxyname/data/jensen+24_wide.fits"

# ╔═╡ ec227641-86e6-46b7-8019-9b02072ed9f7
begin 
	all_stars = LilGuys.read_fits(datafile)

	all_stars[:, :r_ell_old] = all_stars[:, :r_ell]
	all_stars.xi, all_stars.eta = LilGuys.to_tangent(all_stars.ra, all_stars.dec, observed_properties["ra"], observed_properties["dec"])

	all_stars.r_ell = 60LilGuys.calc_r_ell(all_stars.xi, all_stars.eta, observed_properties["ellipticity"], observed_properties["position_angle"])
	
	all_stars
end

# ╔═╡ 082a06dd-eeb5-4761-a233-1ee89e8cb819
best_stars = all_stars[all_stars.F_BEST .== 1.0, :]

# ╔═╡ 88fbdd09-30be-4fc3-95ae-acce6e0018e1
members = best_stars[best_stars.PSAT .> 0.2, :]

# ╔═╡ 731ea468-5003-44e9-95b8-7fa7ef4b475b
Nmemb = size(members, 1)

# ╔═╡ 60d0e593-88fd-4b4c-9009-cc24a597c6d5
members_nospace = best_stars[best_stars.PSAT_NOSPACE .> 0.2, :]

# ╔═╡ a9d94121-ea6e-416a-bae8-aa93c16bde72
md"""
# Utilities
"""

# ╔═╡ 965217b8-b2a5-485b-83de-cac887065b19
plot_labels = OrderedDict(
	:xi => L"$\xi$\,/\,degrees",
	:eta => L"$\eta$\,/\,degrees",
	:xi_am => L"$\xi$\,/\,arcmin",
	:eta_am => L"$\eta$\,/\,arcmin",
	:G => "G (mag)",
	:bp_rp => "BP – RP (mag)",
	:pmra => L"$\mu_{\alpha*}$ / mas yr$^{-1}$",
	:pmdec => L"$\mu_{\delta}$ / mas yr$^{-1}$",
)

# ╔═╡ 2d4da56b-5d0a-49d3-83ae-a90f85192101
θmax = maximum(sqrt.(all_stars.xi.^2 .+ all_stars.eta .^ 2))

# ╔═╡ 5a9b5529-50f1-4cb5-b716-42180ea77d5f
md"""
# Simple selection plots
"""

# ╔═╡ 77f69d97-f71a-48e9-a048-1bb520222855
md"""
The plots below show various subsamples of the dataset in different planes to get a handle on how generally members are selected using the J+24 algorithm.
"""

# ╔═╡ 12aa5708-dfee-4b48-8835-69055ad82680
member_markersize = 3

# ╔═╡ 57f558c1-ca31-44bb-a681-a2ed083a0b70
bg_markersize = 2

# ╔═╡ e4d42cbd-0166-43d8-838c-27a14956623b
datasets = OrderedDict(
	:best => best_stars,
	:members => members,
)

# ╔═╡ 6d343a21-5bf9-4329-a65d-f6ad30c00c5e
scatter_kwargs = Dict(
	:best => (;	alpha=0.1, markersize=bg_markersize, color=:black, 
		label="all" => (alpha=1, markersize=3),
		rasterize=true,
	),
	:members => (;
		alpha=1, markersize=member_markersize,
		label = "members",
		color=COLORS[2]
	),
)

# ╔═╡ 8a4f4283-5699-44ca-956a-869a41177f05
Arya.update_fontsize!(12)

# ╔═╡ 77b7678f-158f-46cb-82b1-2e719ec6895a
function compare_samples(datasets, scatter_kwargs) 
	fig = Figure(
		size = (5.39*72, 5.39*72),
	)

	# tangent
	dθ = 60θmax
	ax = Axis(fig[1, 1], 
		xlabel=plot_labels[:xi_am], 
		ylabel=plot_labels[:eta_am], 
		aspect=1, 
		limits=(-dθ, dθ, -dθ, dθ), 
		xgridvisible=false, ygridvisible=false
	)


	for (label, df) in datasets
		scatter!(60df.xi, 60df.eta; scatter_kwargs[label]...)
	end

	# cmd
	ax =  Axis(fig[1,2], 
		yreversed=true,
		xlabel=plot_labels[:bp_rp],
		ylabel=plot_labels[:G],
		limits=(-0.5, 2.5, 15, 22),
		xgridvisible=false,
		ygridvisible=false,
		aspect = 1,
	)

	for (label, df) in datasets
		scatter!(df.bp_rp, df.phot_g_mean_mag; scatter_kwargs[label]...)
	end
	axislegend(position=:lt)

	# proper motions
	ax = Axis(fig[2,1],
		xlabel = plot_labels[:pmra],
		ylabel = plot_labels[:pmdec],
		aspect=DataAspect(),
		limits=(-10, 10, -10, 10),
		xgridvisible=false,
		ygridvisible=false,
		)
	
	for (label, df) in datasets
		scatter!(df.pmra, df.pmdec; scatter_kwargs[label]...)
	end



	# parallax
	ax = Axis(fig[2,2],
		xlabel=L"\varpi\ /\ \textrm{mas}", 
		ylabel=L"\delta\varpi\ /\ \textrm{mas}",
		limits=(-5, 5, 0, 3),
		aspect = 1,
		xgridvisible=false, ygridvisible=false,
	)

	for (label, df) in datasets
		scatter!(df.parallax, df.parallax_error; scatter_kwargs[label]...)
	end
	
	
	fig
end

# ╔═╡ b94901b0-ecc7-480b-b24a-fc526c9491c8
@savefig "scl_selection" compare_samples(
		(
		:best => best_stars,
		:members => members,
	),
	Dict(
		:best => (;	alpha=0.1, markersize=1, color=:black, 
			label="all" => (alpha=1, markersize=3),
			rasterize=true,
		),
		:members => (;
			alpha=1, markersize=1,
			label = "members" =>(alpha=1, markersize=3),
			color=COLORS[6]
		),
	)
)

# ╔═╡ Cell order:
# ╟─47b8b3b0-0228-4f50-9da4-37d388ef9e9f
# ╠═eca9c1c5-e984-42d4-8854-b227fdec0a8a
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═69c98029-165c-407b-9a63-a27e06e30e45
# ╠═0004f638-c57a-4dab-8b97-c77840cafbbf
# ╠═ae29bed0-6700-47f1-8952-35e867ce126b
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╠═cf7aeb0a-d453-462a-b43d-a832567440fd
# ╟─2eb4aa78-0fea-460b-a18e-06a129c41504
# ╠═5db5adc4-9b98-4025-9b3c-65b40e5d4c59
# ╠═9ecf79a8-2ed3-40c6-b555-a102250ecbd4
# ╠═26cf1867-02be-4d36-8c35-6c58a1feca27
# ╠═ec227641-86e6-46b7-8019-9b02072ed9f7
# ╠═731ea468-5003-44e9-95b8-7fa7ef4b475b
# ╠═88fbdd09-30be-4fc3-95ae-acce6e0018e1
# ╠═60d0e593-88fd-4b4c-9009-cc24a597c6d5
# ╠═082a06dd-eeb5-4761-a233-1ee89e8cb819
# ╟─a9d94121-ea6e-416a-bae8-aa93c16bde72
# ╠═965217b8-b2a5-485b-83de-cac887065b19
# ╠═2d4da56b-5d0a-49d3-83ae-a90f85192101
# ╟─5a9b5529-50f1-4cb5-b716-42180ea77d5f
# ╟─77f69d97-f71a-48e9-a048-1bb520222855
# ╠═12aa5708-dfee-4b48-8835-69055ad82680
# ╠═57f558c1-ca31-44bb-a681-a2ed083a0b70
# ╠═e4d42cbd-0166-43d8-838c-27a14956623b
# ╠═6d343a21-5bf9-4329-a65d-f6ad30c00c5e
# ╠═8a4f4283-5699-44ca-956a-869a41177f05
# ╠═77b7678f-158f-46cb-82b1-2e719ec6895a
# ╠═b94901b0-ecc7-480b-b24a-fc526c9491c8
