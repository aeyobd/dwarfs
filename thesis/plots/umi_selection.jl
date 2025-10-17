### A Pluto.jl notebook ###
# v0.20.19

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

# ╔═╡ 3083e55a-05f9-45aa-96ba-8f817d97c19c
CairoMakie.activate!(type=:png, px_per_unit=2)

# ╔═╡ 9c7a3d3b-d4a4-4e15-a183-d3b4a8cf39bc
module Utils 
	include("gaia_utils.jl")
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
rv_members = read_fits(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/velocities/processed/rv_combined_x_2c_psat_0.2.fits")

# ╔═╡ 5a9b5529-50f1-4cb5-b716-42180ea77d5f
md"""
# Simple selection plots
"""

# ╔═╡ 77f69d97-f71a-48e9-a048-1bb520222855
md"""
The plots below show various subsamples of the dataset in different planes to get a handle on how generally members are selected using the J+24 algorithm.
"""

# ╔═╡ e7fd8470-9225-4dfc-b4ff-08a8e40069fa
R_h = observed_properties["R_h"]

# ╔═╡ 96360af3-e877-4c83-88ef-464ea9eb7b14
rv_distant = rv_members[rv_members.R_ell .> 6R_h, :]

# ╔═╡ 1b31378d-480d-4a7a-9007-9e0b5fd512e1
md"""
Compare the above table against sestito+23b
"""

# ╔═╡ 9eaaf2be-936b-4d26-9e2e-c04f551c515c
rv_distant[:, [:xi, :eta, :R_ell]]

# ╔═╡ 52506697-0d8b-4242-97c0-1e1f56da18a1
rv_distant.R_ell ./ R_h

# ╔═╡ 6fa9d3a7-3762-409c-8dc8-9219add216f9
xi_m1, eta_m1 = -42.086475336200905, -14.6280146094025

# ╔═╡ da9ce4d1-0af1-417b-bb57-197bd101d189
samples = OrderedDict(
	:best => best_stars,
	:members_nospace => members_nospace,
	:members => members,
	:rv => rv_members,
	:rv_distant => rv_distant,
)

# ╔═╡ b94901b0-ecc7-480b-b24a-fc526c9491c8
@savefig "umi_selection" let

	fig = Utils.compare_j24_samples(samples, Utils.default_styles,
		observed_properties,
		legend_position=:rt,
		title = "Ursa Minor"
	)


	ax = fig.content[2]

	Makie.current_axis!(ax)
	Utils.ellipse!(3*0.5, 0, 0, x0=xi_m1, y0=eta_m1, color=COLORS[5])
	text!(ax, xi_m1, eta_m1, color=COLORS[5], text="Muñoz 1", fontsize = 0.8 * theme(:fontsize)[]*0.8, align=(:left, :center), offset=(12, 0))

	fig
end

# ╔═╡ Cell order:
# ╟─47b8b3b0-0228-4f50-9da4-37d388ef9e9f
# ╠═eca9c1c5-e984-42d4-8854-b227fdec0a8a
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═2d5297cd-6a01-4b26-ac77-995b878d765d
# ╠═69c98029-165c-407b-9a63-a27e06e30e45
# ╠═3083e55a-05f9-45aa-96ba-8f817d97c19c
# ╠═9c7a3d3b-d4a4-4e15-a183-d3b4a8cf39bc
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
# ╠═bc87bc28-167d-493d-9553-e90afeaee2ee
# ╟─5a9b5529-50f1-4cb5-b716-42180ea77d5f
# ╟─77f69d97-f71a-48e9-a048-1bb520222855
# ╠═e7fd8470-9225-4dfc-b4ff-08a8e40069fa
# ╠═96360af3-e877-4c83-88ef-464ea9eb7b14
# ╟─1b31378d-480d-4a7a-9007-9e0b5fd512e1
# ╠═9eaaf2be-936b-4d26-9e2e-c04f551c515c
# ╠═52506697-0d8b-4242-97c0-1e1f56da18a1
# ╠═6fa9d3a7-3762-409c-8dc8-9219add216f9
# ╠═da9ce4d1-0af1-417b-bb57-197bd101d189
# ╠═b94901b0-ecc7-480b-b24a-fc526c9491c8
