### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ bff50014-bfa9-11ee-33f0-0f67e543c2d4
begin 
	import Pkg; Pkg.activate()

	import PythonCall
	using DataFrames 
	using CSV
	using CairoMakie
end

# ╔═╡ c3af8de5-611a-457e-8ba0-636803c5a76c
using PyFITS

# ╔═╡ 2d5297cd-6a01-4b26-ac77-995b878d765d
using Arya

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

# ╔═╡ d3bd7158-ea70-47a0-9800-9bfc08f3557c
include(ENV["DWARFS_ROOT"] * "/utils/gaia_filters.jl")

# ╔═╡ 47b8b3b0-0228-4f50-9da4-37d388ef9e9f
md"""
# Jensen et al. 2024 sammple

Some plots to understand the (unmodified) J+24 data sample.
The goals here are to investigate the main density profile, likelihoods, and what the cuts appear as in CMD space. 

"""

# ╔═╡ eca9c1c5-e984-42d4-8854-b227fdec0a8a
galaxyname = "fornax"

# ╔═╡ a4848423-9be3-48d7-98c2-8d1096eb2560
module Utils
	include("utils.jl")
end

# ╔═╡ 0004f638-c57a-4dab-8b97-c77840cafbbf
import TOML

# ╔═╡ 2eb4aa78-0fea-460b-a18e-06a129c41504
md"""
# Inputs & Data loading
"""

# ╔═╡ 5db5adc4-9b98-4025-9b3c-65b40e5d4c59
obs_dir = ENV["DWARFS_ROOT"] * "/observations/" 

# ╔═╡ 9ecf79a8-2ed3-40c6-b555-a102250ecbd4
observed_properties = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/" * galaxyname * "/observed_properties.toml")

# ╔═╡ 26cf1867-02be-4d36-8c35-6c58a1feca27
datafile = obs_dir * "/$galaxyname/data/jensen+24_1c.fits"

# ╔═╡ 90cec348-1947-4091-a5dd-ae67cf80fddb
filt_params = GaiaFilterParams(observed_properties, filename=datafile)

# ╔═╡ ec227641-86e6-46b7-8019-9b02072ed9f7
all_stars = read_gaia_stars(filt_params)

# ╔═╡ 082a06dd-eeb5-4761-a233-1ee89e8cb819
best_stars = all_stars[all_stars.F_BEST .== 1.0, :]

# ╔═╡ 88fbdd09-30be-4fc3-95ae-acce6e0018e1
members = best_stars[best_stars.PSAT .> 0.2, :]

# ╔═╡ 731ea468-5003-44e9-95b8-7fa7ef4b475b
Nmemb = size(members, 1)

# ╔═╡ 60d0e593-88fd-4b4c-9009-cc24a597c6d5
members_nospace = best_stars[best_stars.LLR_nospace .> 0.0, :]

# ╔═╡ cc2de56d-9fbc-4bee-a52c-c9034a246b92
rv_members = read_fits(ENV["DWARFS_ROOT"] * "/observations/fornax/velocities/processed/rv_combined_psat_0.2.fits")

# ╔═╡ 5a9b5529-50f1-4cb5-b716-42180ea77d5f
md"""
# Simple selection plots
"""

# ╔═╡ 77f69d97-f71a-48e9-a048-1bb520222855
md"""
The plots below show various subsamples of the dataset in different planes to get a handle on how generally members are selected using the J+24 algorithm.
"""

# ╔═╡ 8a4f4283-5699-44ca-956a-869a41177f05
Arya.update_fontsize!(12)

# ╔═╡ 45559d19-2e95-4ada-bba6-32c2da2b8953
my_colors = Utils.SELECTION_COLORS

# ╔═╡ 4084c8fa-bbf7-4ce7-8207-d185cc090f34
import Random

# ╔═╡ b94901b0-ecc7-480b-b24a-fc526c9491c8
@savefig "fornax_selection" Utils.compare_j24_samples(
		OrderedDict(
		:best => best_stars,
		:members_nospace => members_nospace,
		:members => members,
		:rv => rv_members[Random.randperm(size(rv_members, 1)), :],
		),
	Dict(
		:best => (;	alpha=1, markersize=0.5, color=my_colors[1], 
			label="all" => (;markersize=2),
			rasterize=10,
		),
		:members_nospace => (;
			alpha=1, markersize=1,
			label = "CMD + PM" =>(alpha=1, markersize=2),
			color=my_colors[2],
			strokecolor=my_colors[2],
			strokewidth=0.3
		),
		:members => (;
			markersize=3,
			label = L"P_\textrm{sat} > 0.2" =>(alpha=1, markersize=2*2),
			marker=:x,
			#color=:transparent,
			#strokewidth=0.3,
			color = my_colors[3],
			alpha=1,
			strokecolor = my_colors[2],
			strokewidth=0.0
		),
		:rv => (;
			markersize=4,
			marker=:diamond,
			label = "RV members" =>(alpha=1, markersize=2.5*2),
			color = my_colors[4],
			strokecolor = :black,
			strokewidth=0.0
		),
		:rv_distant => (;
			markersize=8,
			marker=:star5,
			color = my_colors[4],
			strokewidth=1,
			strokecolor=COLORS[4]
		),
	),
	observed_properties,
	legend_position=:rt,
	age = 2,
)

# ╔═╡ 15bf78a4-f246-4849-af09-db7a9f249f75
md"""
# Extra
"""

# ╔═╡ 7ab9871b-1384-409b-b2fc-bc267d20ce3c
Ag, Ab, Ar = Utils.get_extinction(best_stars.ra, best_stars.dec, best_stars.bp_rp)

# ╔═╡ a2e13309-88b5-48da-b65d-1e711bdb316c
rv_nonmembers = read_fits(ENV["DWARFS_ROOT"] * "/observations/fornax/velocities/processed/rv_combined_psat_0.2_nonmemb.fits")

# ╔═╡ 59ca6804-5d00-4ae9-9f23-bde3c911b4c1
scatter(best_stars.xi, best_stars.eta, color=Ag, markersize=5)

# ╔═╡ 63881146-5d8a-4285-b7a6-60bbde98e9da
Utils.compare_j24_samples(
		OrderedDict(
			:members => rv_members,
		  :nonmemb => rv_nonmembers,	
		),
	Dict(
		:nonmemb => (;	alpha=1, markersize=2, color=:red, 
			label="nonmemb" => (alpha=1, markersize=2),
		),
		:members => (;	alpha=1, markersize=2, color=:black, 
			label="memb" => (alpha=1, markersize=2),
		),

	),
	observed_properties,
	age= 2
)

# ╔═╡ Cell order:
# ╟─47b8b3b0-0228-4f50-9da4-37d388ef9e9f
# ╠═eca9c1c5-e984-42d4-8854-b227fdec0a8a
# ╠═c3af8de5-611a-457e-8ba0-636803c5a76c
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═2d5297cd-6a01-4b26-ac77-995b878d765d
# ╠═69c98029-165c-407b-9a63-a27e06e30e45
# ╠═a4848423-9be3-48d7-98c2-8d1096eb2560
# ╠═0004f638-c57a-4dab-8b97-c77840cafbbf
# ╠═ae29bed0-6700-47f1-8952-35e867ce126b
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╠═d3bd7158-ea70-47a0-9800-9bfc08f3557c
# ╟─2eb4aa78-0fea-460b-a18e-06a129c41504
# ╠═5db5adc4-9b98-4025-9b3c-65b40e5d4c59
# ╠═9ecf79a8-2ed3-40c6-b555-a102250ecbd4
# ╠═26cf1867-02be-4d36-8c35-6c58a1feca27
# ╠═90cec348-1947-4091-a5dd-ae67cf80fddb
# ╠═ec227641-86e6-46b7-8019-9b02072ed9f7
# ╠═731ea468-5003-44e9-95b8-7fa7ef4b475b
# ╠═88fbdd09-30be-4fc3-95ae-acce6e0018e1
# ╠═60d0e593-88fd-4b4c-9009-cc24a597c6d5
# ╠═082a06dd-eeb5-4761-a233-1ee89e8cb819
# ╠═cc2de56d-9fbc-4bee-a52c-c9034a246b92
# ╟─5a9b5529-50f1-4cb5-b716-42180ea77d5f
# ╟─77f69d97-f71a-48e9-a048-1bb520222855
# ╠═8a4f4283-5699-44ca-956a-869a41177f05
# ╠═45559d19-2e95-4ada-bba6-32c2da2b8953
# ╠═4084c8fa-bbf7-4ce7-8207-d185cc090f34
# ╠═b94901b0-ecc7-480b-b24a-fc526c9491c8
# ╟─15bf78a4-f246-4849-af09-db7a9f249f75
# ╠═7ab9871b-1384-409b-b2fc-bc267d20ce3c
# ╠═a2e13309-88b5-48da-b65d-1e711bdb316c
# ╠═59ca6804-5d00-4ae9-9f23-bde3c911b4c1
# ╠═63881146-5d8a-4285-b7a6-60bbde98e9da
