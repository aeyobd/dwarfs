### A Pluto.jl notebook ###
# v0.20.8

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
include("style.jl")

# ╔═╡ d3bd7158-ea70-47a0-9800-9bfc08f3557c
include(ENV["DWARFS_ROOT"] * "/utils/gaia_filters.jl")

# ╔═╡ 47b8b3b0-0228-4f50-9da4-37d388ef9e9f
md"""
# Jensen et al. 2024 sammple

This notebook creates plots of Jaclyn's sample for sculptor or ursa minor (controled with `galaxyname`).
The idea is to show the background compared to selected stars. 
Additionally, versions of the plots are made with very distant member stars highlighted.

"""

# ╔═╡ eca9c1c5-e984-42d4-8854-b227fdec0a8a
galaxyname = "sculptor"

# ╔═╡ 0004f638-c57a-4dab-8b97-c77840cafbbf
import TOML

# ╔═╡ a4848423-9be3-48d7-98c2-8d1096eb2560
module Utils 
	include("utils.jl")
end 

# ╔═╡ 8a77b19e-f760-4899-88e5-1c53b1406e87
CairoMakie.activate!(type=:png)

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
datafile = if galaxyname == "sculptor"
	obs_dir * "/$galaxyname/data/jensen+24_wide_2c.fits"
elseif galaxyname == "ursa_minor"
	obs_dir * "/$galaxyname/data/jensen+24_2c.fits"
end

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

# ╔═╡ 29d16d4e-a6b8-4945-860e-cc8690426d49
rv_file = Dict(
	"sculptor" => "rv_combined_x_wide_2c_psat_0.2.fits",
	"ursa_minor" => "rv_combined_x_2c_psat_0.2.fits",
)[galaxyname]

# ╔═╡ bc87bc28-167d-493d-9553-e90afeaee2ee
rv_members = read_fits(ENV["DWARFS_ROOT"] * "/observations/$galaxyname/velocities/processed/$rv_file")

# ╔═╡ 766ea9d3-d1e1-4a4d-b6ef-5a3528cb1579
R_h = observed_properties["R_h"]

# ╔═╡ a04a394e-4d4b-4d95-baab-9886e151ec44
if galaxyname == "sculptor"
	rv_distant = rv_members[rv_members.R_ell .> 6R_h, :]
else
	rv_distant = rv_members[rv_members.R_ell .> 6R_h, :]
end

# ╔═╡ b7244cfb-7326-4e81-b751-0121f202f2b3
rv_distant.R_ell ./ R_h

# ╔═╡ 5a9b5529-50f1-4cb5-b716-42180ea77d5f
md"""
# Simple selection plots
"""

# ╔═╡ 77f69d97-f71a-48e9-a048-1bb520222855
md"""
The plots below show various subsamples of the dataset in different planes to get a handle on how generally members are selected using the J+24 algorithm.
"""

# ╔═╡ dc0fa286-0dbb-4da5-bfae-ebe3655f7d8a
grey = colorant"#808080"

# ╔═╡ 74a3920c-f93d-4a7c-932a-2938cd6bb020
samples = OrderedDict(
	:best => best_stars,
	:members => members,
	:distant => rv_distant
)

# ╔═╡ 5af9ddf6-7d14-4120-9d20-bcd20e756a09
samples_norv = OrderedDict(
	:best => best_stars,
	:members => members,
)

# ╔═╡ 31bac4e9-4f19-4391-8a2a-4408c746a753
styles = Dict(
	:best => (;	markersize=3, color=grey, 
		label="all" => (;markersize=12, alpha=0.3),
		rasterize=10,
		alpha=0.1
	),
	:members => (;
		markersize=12,
		label = "members" =>(; alpha=1, markersize=24),
		marker=:diamond,
		color = COLORS[1],
		alpha=1,
	),
	:distant => (;
		 markersize=36,
		 label = "distant RV members", 
		 marker = :star5,
		 color = COLORS[4],
		 alpha=1, 
		 strokecolor=:black, 
		 strokewidth = 0.0
	)
)

# ╔═╡ b94901b0-ecc7-480b-b24a-fc526c9491c8
# @savefig "$(galaxyname)_selection_w_rv" let
# 	fig = Utils.compare_j24_samples(
# 		samples, styles,
# 		observed_properties
# 	)

# 	Makie.resize_to_layout!(fig)
# 	fig

# end

# ╔═╡ 650f211d-ad9f-4e83-8920-caf0a3bb7dce
extrema(log10.(members.L_CMD_SAT))

# ╔═╡ 0ef75901-9182-4bed-9aa2-c9c1bc574c69
@savefig "$(galaxyname)_selection_gradient" let
	fig = Figure(size=(1720, 800))
	colorrange=(-5, 0)
	age = 12
	kwargs_bg = styles[:best]

	dθ = maximum(sqrt.(all_stars.xi.^2 .+ all_stars.eta .^ 2))
	ax = Axis(fig[1,1],
		xlabel=Utils.plot_labels[:xi_am], 
		ylabel=Utils.plot_labels[:eta_am], 
		#aspect=1, 
		limits=(-dθ, dθ, -dθ, dθ), 
		xgridvisible=false, ygridvisible=false,
		xreversed = true
	)
	scatter!(best_stars.xi, best_stars.eta; kwargs_bg... )
	scatter!(members.xi, members.eta, color=log10.(members.L_S_SAT), label="members", colorrange=colorrange)
	
	Utils.ellipse!(6observed_properties["R_h"], observed_properties["ellipticity"], observed_properties["position_angle"], color=:black)
	
	b = 6observed_properties["R_h"] * sqrt(1 - observed_properties["ellipticity"])
    θ = observed_properties["position_angle"]
	
    text!(b*cosd(θ), -b*sind(θ), text=L"6R_h", rotation=deg2rad(θ-90), 
          align = (:center, :bottom), color=:black, fontsize=0.8*theme(:fontsize)[])

	# CMD

	isochrone = Utils.get_isochrone(observed_properties["metallicity"], age)
    G_iso = isochrone.Gmag .+ observed_properties["distance_modulus"]
    BP_RP_iso = isochrone.G_BPftmag .- isochrone.G_RPmag

	Gmin = minimum(best_stars.G) - 0.2
    Gmax = 21
	ax_cmd =  Axis(fig[1,2], 
		yreversed=true,
		xlabel=Utils.plot_labels[:bp_rp],
		ylabel=Utils.plot_labels[:G],
		limits=(-0.5, 2.5, Gmin, Gmax),
		xticks = 0:2,
		xgridvisible=false,
		ygridvisible=false,
		#aspect = 1,
	)

	Ag, Ab, Ar = Utils.get_extinction(best_stars.ra, best_stars.dec, best_stars.bp_rp)
	scatter!(best_stars.bp_rp .- Ab .+ Ar, best_stars.phot_g_mean_mag .- Ag; kwargs_bg...)

	Ag, Ab, Ar = Utils.get_extinction(members.ra, members.dec, members.bp_rp)

	scatter!(members.bp_rp .- Ab .+ Ar, members.G .- Ag, color=log10.(members.L_CMD_SAT), colorrange=colorrange)

	lines!(BP_RP_iso, G_iso, color=:black)


	# PM
	 ax_pm = Axis(fig[1,3],
		xlabel = Utils.plot_labels[:pmra],
		ylabel = Utils.plot_labels[:pmdec],
		#aspect=DataAspect(),
		limits=(-10, 10, -10, 10),
		xgridvisible=false,
		ygridvisible=false,
		)
	
	scatter!(best_stars.pmra, best_stars.pmdec; kwargs_bg... )
	p = scatter!(members.pmra, members.pmdec, color=log10.(members.L_PM_SAT), colorrange=colorrange)
	
	scatter!([observed_properties["pmra"]], [observed_properties["pmdec"]],
		 color=:black
   )


	Colorbar(fig[1,4], p, width=theme(:fontsize), label="log likelihood (component)")

	
	rowsize!(fig.layout, 1, Aspect(1, 1))
	Legend(fig[2, 2], ax, tellwidth=false, tellheight=true, nbanks=3)

	Makie.resize_to_layout!(fig)

	fig
end

# ╔═╡ b4cef660-5ec1-429d-972b-fb50dd0b155c
@savefig "$(galaxyname)_selection" let
	fig = Utils.compare_j24_samples(
		samples_norv, styles,
		observed_properties
	)

	Makie.resize_to_layout!(fig)
	fig

end

# ╔═╡ 04073727-7f96-4d6b-a0d4-0b3975955461
@savefig "$(galaxyname)_rv_stars" let
	fig = Figure()
	
	Utils.plot_tangent!(fig[1,1],
		samples, styles,
		observed_properties
	)

	fig.content[1].aspect = DataAspect()
	
	fig

end

# ╔═╡ 889cc1ea-f453-42a3-8036-8ee12191be7f
md"""
# Extra
"""

# ╔═╡ 26574dbe-b201-4ea4-aa37-239c4d38b271
R_h

# ╔═╡ 0363e1b8-c9be-4059-a2a5-166f01d7fe56
prof = LilGuys.Exp2D()

# ╔═╡ e52ae552-1221-4871-8877-65cf33eb5a4c
α = LilGuys.R_h(prof)

# ╔═╡ d7c254be-6030-451b-a947-d74a17ee62f7
(1 - LilGuys.mass_2D(prof, 6α)) *  sum(best_stars.PSAT)

# ╔═╡ Cell order:
# ╟─47b8b3b0-0228-4f50-9da4-37d388ef9e9f
# ╠═eca9c1c5-e984-42d4-8854-b227fdec0a8a
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═2d5297cd-6a01-4b26-ac77-995b878d765d
# ╠═0004f638-c57a-4dab-8b97-c77840cafbbf
# ╠═ae29bed0-6700-47f1-8952-35e867ce126b
# ╠═69c98029-165c-407b-9a63-a27e06e30e45
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╠═a4848423-9be3-48d7-98c2-8d1096eb2560
# ╠═8a77b19e-f760-4899-88e5-1c53b1406e87
# ╠═7865b1b1-3a55-41cf-b566-a6770c471077
# ╠═d3bd7158-ea70-47a0-9800-9bfc08f3557c
# ╟─2eb4aa78-0fea-460b-a18e-06a129c41504
# ╠═5db5adc4-9b98-4025-9b3c-65b40e5d4c59
# ╠═9ecf79a8-2ed3-40c6-b555-a102250ecbd4
# ╠═26cf1867-02be-4d36-8c35-6c58a1feca27
# ╠═90cec348-1947-4091-a5dd-ae67cf80fddb
# ╠═ec227641-86e6-46b7-8019-9b02072ed9f7
# ╠═731ea468-5003-44e9-95b8-7fa7ef4b475b
# ╠═88fbdd09-30be-4fc3-95ae-acce6e0018e1
# ╠═082a06dd-eeb5-4761-a233-1ee89e8cb819
# ╠═29d16d4e-a6b8-4945-860e-cc8690426d49
# ╠═bc87bc28-167d-493d-9553-e90afeaee2ee
# ╠═766ea9d3-d1e1-4a4d-b6ef-5a3528cb1579
# ╠═a04a394e-4d4b-4d95-baab-9886e151ec44
# ╠═b7244cfb-7326-4e81-b751-0121f202f2b3
# ╟─5a9b5529-50f1-4cb5-b716-42180ea77d5f
# ╟─77f69d97-f71a-48e9-a048-1bb520222855
# ╠═dc0fa286-0dbb-4da5-bfae-ebe3655f7d8a
# ╠═74a3920c-f93d-4a7c-932a-2938cd6bb020
# ╠═5af9ddf6-7d14-4120-9d20-bcd20e756a09
# ╠═31bac4e9-4f19-4391-8a2a-4408c746a753
# ╠═b94901b0-ecc7-480b-b24a-fc526c9491c8
# ╠═650f211d-ad9f-4e83-8920-caf0a3bb7dce
# ╠═0ef75901-9182-4bed-9aa2-c9c1bc574c69
# ╠═b4cef660-5ec1-429d-972b-fb50dd0b155c
# ╠═04073727-7f96-4d6b-a0d4-0b3975955461
# ╟─889cc1ea-f453-42a3-8036-8ee12191be7f
# ╠═26574dbe-b201-4ea4-aa37-239c4d38b271
# ╠═0363e1b8-c9be-4059-a2a5-166f01d7fe56
# ╠═e52ae552-1221-4871-8877-65cf33eb5a4c
# ╠═d7c254be-6030-451b-a947-d74a17ee62f7
