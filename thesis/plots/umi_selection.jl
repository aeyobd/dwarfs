### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ bff50014-bfa9-11ee-33f0-0f67e543c2d4
begin 
	import Pkg; Pkg.activate()

	using PyFITS
	using DataFrames 
	using CSV
	using CairoMakie
end

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

# ╔═╡ 9c7a3d3b-d4a4-4e15-a183-d3b4a8cf39bc
include("utils.jl")

# ╔═╡ d3bd7158-ea70-47a0-9800-9bfc08f3557c
include(ENV["DWARFS_ROOT"] * "/utils/gaia_filters.jl")

# ╔═╡ 47b8b3b0-0228-4f50-9da4-37d388ef9e9f
md"""
# Jensen et al. 2024 sample

Some plots to understand the (unmodified) J+24 data sample.
The goals here are to investigate the main density profile, likelihoods, and what the cuts appear as in CMD space. 

"""

# ╔═╡ eca9c1c5-e984-42d4-8854-b227fdec0a8a
galaxyname = "ursa_minor"

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
datafile = obs_dir * "/$galaxyname/data/jensen+24_2c.fits"

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

# ╔═╡ bc87bc28-167d-493d-9553-e90afeaee2ee
rv_members = read_fits(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/velocities/processed/rv_members_all.fits")

# ╔═╡ a9d94121-ea6e-416a-bae8-aa93c16bde72
md"""
# Utilities
"""

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

# ╔═╡ e7fd8470-9225-4dfc-b4ff-08a8e40069fa
R_h = observed_properties["r_h"]

# ╔═╡ 96360af3-e877-4c83-88ef-464ea9eb7b14
rv_distant = rv_members[rv_members.R_ell .> 3R_h, :]

# ╔═╡ 9eaaf2be-936b-4d26-9e2e-c04f551c515c
rv_distant[:, [:xi, :eta, :R_ell]]

# ╔═╡ ec974caa-59e6-4555-b9b4-48f33d31d00b
rv_distant.R_ell ./ R_h

# ╔═╡ 43e9b9bd-b7ba-4d2c-b97e-2d28c1ba270d
import Random

# ╔═╡ a47cbe4c-15ec-4bb0-8992-e7650e1b072f
my_colors = SELECTION_COLORS

# ╔═╡ b94901b0-ecc7-480b-b24a-fc526c9491c8
@savefig "umi_selection" compare_j24_samples(
		OrderedDict(
		:best => best_stars,
		:members_nospace => members_nospace,
		:members => members,
		:rv => rv_members[Random.randperm(size(rv_members, 1)), :],
		:rv_distant => rv_distant,
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
			marker=:diamond,
			#color=:transparent,
			#strokewidth=0.3,
			color = my_colors[3],
			alpha=1,
			strokecolor = my_colors[2],
			strokewidth=0.0
		),
		:rv => (;
			markersize=4,
			marker=:xcross,
			label = "RV members" =>(alpha=1, markersize=2.5*2),
			color = my_colors[4],
			strokecolor = :black,
			strokewidth=0.0
		),
		:rv_distant => (;
			markersize=6,
			marker=:star5,
			color = my_colors[4],
			strokewidth=0.3,
			strokecolor=:black,
		),
	),
	observed_properties
)

# ╔═╡ 90436a7d-e357-4365-922a-42f61e3dfb22
md"""
# Extra
"""

# ╔═╡ 5d0a5768-85ba-4913-9108-9302cc63cf4b
rv_all = read_fits(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/velocities/processed/rv_combined.fits")

# ╔═╡ a776c46b-4a7e-4857-aac2-4120c6d1d4d2
id_nonmemb = setdiff(rv_all.source_id, rv_members.source_id)

# ╔═╡ 55ab1364-bcf5-4c29-8081-cc3808b4e53c
rv_nonmemb = rv_all[rv_all.source_id .∈ [id_nonmemb], :]

# ╔═╡ 0f5143a6-f892-4a3e-967d-5e6b04cb12f0
compare_j24_samples(
		OrderedDict(
			:members => rv_members,
			:nonmemb => rv_nonmemb,	
		),
	Dict(
		:nonmemb => (;	alpha=1, markersize=2, color=:red, 
			label="nonmemb" => (alpha=1, markersize=2),
		),
		:members => (;	alpha=1, markersize=2, color=:black, 
			label="memb" => (alpha=1, markersize=2),
		),
	),
	observed_properties
)

# ╔═╡ Cell order:
# ╠═47b8b3b0-0228-4f50-9da4-37d388ef9e9f
# ╠═eca9c1c5-e984-42d4-8854-b227fdec0a8a
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═2d5297cd-6a01-4b26-ac77-995b878d765d
# ╠═69c98029-165c-407b-9a63-a27e06e30e45
# ╠═9c7a3d3b-d4a4-4e15-a183-d3b4a8cf39bc
# ╠═0004f638-c57a-4dab-8b97-c77840cafbbf
# ╠═ae29bed0-6700-47f1-8952-35e867ce126b
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╠═cf7aeb0a-d453-462a-b43d-a832567440fd
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
# ╠═bc87bc28-167d-493d-9553-e90afeaee2ee
# ╟─a9d94121-ea6e-416a-bae8-aa93c16bde72
# ╟─5a9b5529-50f1-4cb5-b716-42180ea77d5f
# ╟─77f69d97-f71a-48e9-a048-1bb520222855
# ╠═12aa5708-dfee-4b48-8835-69055ad82680
# ╠═57f558c1-ca31-44bb-a681-a2ed083a0b70
# ╠═e7fd8470-9225-4dfc-b4ff-08a8e40069fa
# ╠═96360af3-e877-4c83-88ef-464ea9eb7b14
# ╠═9eaaf2be-936b-4d26-9e2e-c04f551c515c
# ╠═ec974caa-59e6-4555-b9b4-48f33d31d00b
# ╠═43e9b9bd-b7ba-4d2c-b97e-2d28c1ba270d
# ╠═a47cbe4c-15ec-4bb0-8992-e7650e1b072f
# ╠═b94901b0-ecc7-480b-b24a-fc526c9491c8
# ╠═90436a7d-e357-4365-922a-42f61e3dfb22
# ╠═5d0a5768-85ba-4913-9108-9302cc63cf4b
# ╠═a776c46b-4a7e-4857-aac2-4120c6d1d4d2
# ╠═55ab1364-bcf5-4c29-8081-cc3808b4e53c
# ╠═0f5143a6-f892-4a3e-967d-5e6b04cb12f0
