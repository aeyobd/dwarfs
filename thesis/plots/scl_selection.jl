### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ bff50014-bfa9-11ee-33f0-0f67e543c2d4
begin 
	import Pkg; Pkg.activate()


	using DataFrames 
	using CSV
	using CairoMakie
	using PyFITS
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
# Jensen et al. 2024 sammple

Some plots to understand the (unmodified) J+24 data sample.
The goals here are to investigate the main density profile, likelihoods, and what the cuts appear as in CMD space. 

"""

# ╔═╡ eca9c1c5-e984-42d4-8854-b227fdec0a8a
galaxyname = "sculptor"

# ╔═╡ 60985df9-2423-4350-b2f9-4ce7c8e6589a
CairoMakie.activate!(type=:png, px_per_unit=2)

# ╔═╡ 0004f638-c57a-4dab-8b97-c77840cafbbf
import TOML

# ╔═╡ a4848423-9be3-48d7-98c2-8d1096eb2560
module Utils 
	include("utils.jl")
end 

# ╔═╡ ea8e6489-27c8-421b-9512-b0f67a79a294
import Random

# ╔═╡ 2eb4aa78-0fea-460b-a18e-06a129c41504
md"""
# Inputs & Data loading
"""

# ╔═╡ 5db5adc4-9b98-4025-9b3c-65b40e5d4c59
obs_dir = ENV["DWARFS_ROOT"] * "/observations/" 

# ╔═╡ 9ecf79a8-2ed3-40c6-b555-a102250ecbd4
observed_properties = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/" * galaxyname * "/observed_properties.toml")

# ╔═╡ 7865b1b1-3a55-41cf-b566-a6770c471077
observed_properties

# ╔═╡ 26cf1867-02be-4d36-8c35-6c58a1feca27
datafile = obs_dir * "/$galaxyname/data/jensen+24_wide_2c.fits"

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

# ╔═╡ 7b635d02-617f-41df-9c7c-221c04c10e87
size(members, 1) * 0.0005

# ╔═╡ 60d0e593-88fd-4b4c-9009-cc24a597c6d5
members_nospace = best_stars[best_stars.LLR_nospace .> 0.0, :]

# ╔═╡ bc87bc28-167d-493d-9553-e90afeaee2ee
rv_members = read_fits(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/rv_combined_x_wide_2c_psat_0.2.fits")

# ╔═╡ 5a9b5529-50f1-4cb5-b716-42180ea77d5f
md"""
# Simple selection plots
"""

# ╔═╡ 77f69d97-f71a-48e9-a048-1bb520222855
md"""
The plots below show various subsamples of the dataset in different planes to get a handle on how generally members are selected using the J+24 algorithm.
"""

# ╔═╡ 3ab5488e-4c7c-4b8d-9c4b-cc0a04620298
 my_colors = Utils.SELECTION_COLORS

# ╔═╡ 9c66468e-c357-4268-875b-ec83510fd982
md"""
# Extra
"""

# ╔═╡ c1fe9907-cdd8-4b69-a9b3-f2553b25cdf6
R_h = observed_properties["R_h"]

# ╔═╡ bbae75e1-b317-4e4f-a62e-817715a5a095
observed_properties["R_h_inner"]

# ╔═╡ b83ebaf6-cae2-40c0-bd71-98bbbc10eb92
Ag, Ab, Ar = Utils.get_extinction(best_stars.ra, best_stars.dec, best_stars.bp_rp)

# ╔═╡ f37a09e1-3c0f-44a7-b7ad-b48c5553871e
scatter(best_stars.xi, best_stars.eta, color=Ag, markersize=5)

# ╔═╡ 430fae9d-c708-4447-80ce-aabf19b161d2
rv_distant = rv_members[rv_members.R_ell .> 7R_h, :]

# ╔═╡ aeb5e17f-6479-4504-ae4f-c9218f374d48
rv_distant.PSAT

# ╔═╡ b94901b0-ecc7-480b-b24a-fc526c9491c8
@savefig "scl_selection" Utils.compare_j24_samples(
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
			rasterize=2,
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
			label = L"fiducial ($P_\textrm{sat} > 0.2$)" =>(alpha=1, markersize=2*2),
			marker=:rect,
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
	observed_properties
)

# ╔═╡ bc4ad5db-3e90-46e8-ad54-674b02f124c0
rv_members[.!ismissing.(rv_members.RV_gmos), [:xi, :eta, :R_ell]]

# ╔═╡ a373edbf-d441-46ce-bc0f-08ac0d693c44
rv_distant.R_ell ./ observed_properties["R_h"]

# ╔═╡ 9376c4d1-4237-408f-a845-1decf183cf10
 12.33 .* sqrt(1-0.33)

# ╔═╡ 13f558a3-a42e-4384-ac6e-2a036f6e634f
LilGuys.mean(skipmissing(rv_members.RV_t23))

# ╔═╡ 2d474904-ec96-41e7-bd17-8969ea5e3c40
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log R ell",
		ylabel = "[Fe/H]"
			 )

	for (i, col) in enumerate([ :fe_h_t23, :fe_h_apogee, :fe_h_gmos,])
		errcol = "fe_h_err_" * split(string(col), "_")[end]
		filt = .!ismissing.(rv_members[!, col])

		errorscatter!(disallowmissing(log10.(rv_members.R_ell[filt])), disallowmissing(rv_members[filt, col]), yerror=disallowmissing(rv_members[filt, errcol]), color=COLORS[i])

	end
	fig
end

# ╔═╡ 9380d93d-34bb-4e67-a423-a624182d6a56
rv_nonmemb = read_fits(ENV["DWARFS_ROOT"] * "/observations/sculptor/velocities/processed/rv_combined_x_wide_2c_psat_0.2_nonmemb.fits")

# ╔═╡ 4e06c085-3024-4c8a-ba1b-833e3f8630b1
Utils.compare_j24_samples(
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
	observed_properties,
	age= 2
)

# ╔═╡ 065175ac-e754-4629-b56c-a4895c8aec70
alpha_R = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 24609348-e4c3-437e-b36a-fb45d7a7352c
(1 - mass(LilGuys.Exp2D(), 9alpha_R)) * sum(best_stars.PSAT)

# ╔═╡ 28fbacb4-8078-4af2-bf8f-23adab46669b


# ╔═╡ f3a2ae1c-92fc-42cf-89da-c46d34d0c5aa


# ╔═╡ 6b4f3e96-d525-4885-97fd-a3d35171ad76


# ╔═╡ 7f6ab558-627e-46c5-8f47-b1ec03fd2b63


# ╔═╡ c3eef548-6678-46a3-96a1-01246f5f70e0
Arya.get_arya_cmap()

# ╔═╡ f6fe794a-e925-4544-b5ac-fda9eb96629b
import StatsBase: median

# ╔═╡ 508f8d47-f731-4a0a-a761-c1ef2403ec8b
median(members.R_ell)

# ╔═╡ fefd2f20-7116-4fff-b9c5-8bc3628b508f
median(members.dRP)

# ╔═╡ 45fb74e8-de19-44b6-8abd-6e2e9d3af76f
median(members.dBP .+ members.dRP)

# ╔═╡ a8263f41-7b79-4f45-98d7-2d96324bc8fd
median(members.dG)

# ╔═╡ 8b1539ed-a52a-4d12-91e7-50c3690b29b8
isochrone = Utils.get_isochrone(-0.5, 5)

# ╔═╡ 4a9f3b39-ee38-4624-a87b-d15489642f7d
let
	fig = Figure()
		
	ax = Axis(fig[1,1], yreversed=true)
		
	
	p = lines!(isochrone.G_BPbrmag .- isochrone.G_RPmag, isochrone.Gmag .+ observed_properties["distance_modulus"], color = isochrone.Mini, )


	Colorbar(fig[1,2], p)

	fig

end

# ╔═╡ 50d90ce5-99df-46af-b309-a8ccd0ccf921


# ╔═╡ Cell order:
# ╟─47b8b3b0-0228-4f50-9da4-37d388ef9e9f
# ╠═eca9c1c5-e984-42d4-8854-b227fdec0a8a
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═60985df9-2423-4350-b2f9-4ce7c8e6589a
# ╠═2d5297cd-6a01-4b26-ac77-995b878d765d
# ╠═0004f638-c57a-4dab-8b97-c77840cafbbf
# ╠═ae29bed0-6700-47f1-8952-35e867ce126b
# ╠═69c98029-165c-407b-9a63-a27e06e30e45
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╠═a4848423-9be3-48d7-98c2-8d1096eb2560
# ╠═7865b1b1-3a55-41cf-b566-a6770c471077
# ╠═d3bd7158-ea70-47a0-9800-9bfc08f3557c
# ╠═ea8e6489-27c8-421b-9512-b0f67a79a294
# ╟─2eb4aa78-0fea-460b-a18e-06a129c41504
# ╠═5db5adc4-9b98-4025-9b3c-65b40e5d4c59
# ╠═9ecf79a8-2ed3-40c6-b555-a102250ecbd4
# ╠═26cf1867-02be-4d36-8c35-6c58a1feca27
# ╠═90cec348-1947-4091-a5dd-ae67cf80fddb
# ╠═ec227641-86e6-46b7-8019-9b02072ed9f7
# ╠═731ea468-5003-44e9-95b8-7fa7ef4b475b
# ╠═88fbdd09-30be-4fc3-95ae-acce6e0018e1
# ╠═7b635d02-617f-41df-9c7c-221c04c10e87
# ╠═60d0e593-88fd-4b4c-9009-cc24a597c6d5
# ╠═082a06dd-eeb5-4761-a233-1ee89e8cb819
# ╠═bc87bc28-167d-493d-9553-e90afeaee2ee
# ╟─5a9b5529-50f1-4cb5-b716-42180ea77d5f
# ╟─77f69d97-f71a-48e9-a048-1bb520222855
# ╠═3ab5488e-4c7c-4b8d-9c4b-cc0a04620298
# ╠═aeb5e17f-6479-4504-ae4f-c9218f374d48
# ╠═b94901b0-ecc7-480b-b24a-fc526c9491c8
# ╟─9c66468e-c357-4268-875b-ec83510fd982
# ╠═c1fe9907-cdd8-4b69-a9b3-f2553b25cdf6
# ╠═bbae75e1-b317-4e4f-a62e-817715a5a095
# ╠═b83ebaf6-cae2-40c0-bd71-98bbbc10eb92
# ╠═f37a09e1-3c0f-44a7-b7ad-b48c5553871e
# ╠═430fae9d-c708-4447-80ce-aabf19b161d2
# ╠═bc4ad5db-3e90-46e8-ad54-674b02f124c0
# ╠═a373edbf-d441-46ce-bc0f-08ac0d693c44
# ╠═9376c4d1-4237-408f-a845-1decf183cf10
# ╠═508f8d47-f731-4a0a-a761-c1ef2403ec8b
# ╠═13f558a3-a42e-4384-ac6e-2a036f6e634f
# ╠═2d474904-ec96-41e7-bd17-8969ea5e3c40
# ╠═4e06c085-3024-4c8a-ba1b-833e3f8630b1
# ╠═9380d93d-34bb-4e67-a423-a624182d6a56
# ╠═065175ac-e754-4629-b56c-a4895c8aec70
# ╠═24609348-e4c3-437e-b36a-fb45d7a7352c
# ╠═28fbacb4-8078-4af2-bf8f-23adab46669b
# ╠═f3a2ae1c-92fc-42cf-89da-c46d34d0c5aa
# ╠═6b4f3e96-d525-4885-97fd-a3d35171ad76
# ╠═7f6ab558-627e-46c5-8f47-b1ec03fd2b63
# ╠═c3eef548-6678-46a3-96a1-01246f5f70e0
# ╠═f6fe794a-e925-4544-b5ac-fda9eb96629b
# ╠═fefd2f20-7116-4fff-b9c5-8bc3628b508f
# ╠═45fb74e8-de19-44b6-8abd-6e2e9d3af76f
# ╠═a8263f41-7b79-4f45-98d7-2d96324bc8fd
# ╠═8b1539ed-a52a-4d12-91e7-50c3690b29b8
# ╠═4a9f3b39-ee38-4624-a87b-d15489642f7d
# ╠═50d90ce5-99df-46af-b309-a8ccd0ccf921
