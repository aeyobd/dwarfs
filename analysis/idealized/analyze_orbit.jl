### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 061b1886-1878-11ef-3806-b91643300982
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	using LilGuys

	using Arya
end

# ╔═╡ 4018f9cb-b54f-408a-8047-cc995a167034
using OrderedCollections

# ╔═╡ cd43d649-1d52-46a7-a621-68c8122fce6d
using PyFITS

# ╔═╡ 2c702eb7-ebb6-44c9-8e01-ca52d011c014
using HDF5

# ╔═╡ 6b3df4f3-de5b-40b1-8aec-1cfce9aa843a
using PlutoUI

# ╔═╡ 8b41af50-9ae0-475b-bacc-3799e2949b30
md"""
Analyzes the orbit of a n-body halo in a gravitational potential.
Requires the centres to be calculated prior.

This notebook only uses the centres of the orbit, but then calculates useful quantities such as the pericentre, the orbital period, the time of each peri and apocentre and outputs this into the model analysis directory.

"""

# ╔═╡ 27577252-53dc-415c-b9a1-82155ef9e4ca
md"""
# Setup
"""

# ╔═╡ 882d4fc5-07ae-4b06-8da5-67f0894595db
import LinearAlgebra: dot

# ╔═╡ ab57edae-2292-4aef-9c1f-53802dbc0600
import TOML

# ╔═╡ 56d15948-4fd8-4541-9925-99837d9584f5
CairoMakie.activate!(type=:png)

# ╔═╡ 643cd0bf-77b3-4201-9ff7-09dd5aee277c
md"""
# inputs
"""

# ╔═╡ 0e2e7f93-09d1-4902-ad24-223df50d37cb
function notebook_inputs(; kwargs...)
	return PlutoUI.combine() do Child
		
		user_inputs = [
			md""" $(string(name)): $(
				Child(name, obj)
			)"""
			
			for (name, obj) in kwargs
		]
		
		md"""
		#### Inputs
		$(user_inputs)
		"""
	end
end

# ╔═╡ d94663e8-b30e-4712-8a3e-6ef7f954f141
@bind inputs confirm(notebook_inputs(;
	haloname = TextField(default="1e6"),
	orbitname = TextField(default="orbit_"),
))

# ╔═╡ 3953b76f-a726-4211-a6b8-5cf38149dcdf
galaxyname = "idealized"

# ╔═╡ 24ed3bc1-37c6-4601-b941-0780c53a9630
haloname = inputs.haloname

# ╔═╡ 079afa70-c7dc-4347-9ae3-8459ba2fa941
orbitname =  inputs.orbitname

# ╔═╡ 69d83e00-7eb6-4271-838f-80e4d1654dac
modelname = "$galaxyname/$haloname/$orbitname"

# ╔═╡ ac2c7484-9acd-4fda-9699-fdf17da507c2
parentdir = ENV["DWARFS_ROOT"]

# ╔═╡ d142b7bd-3002-4331-a725-577873c42f28
properties_file = joinpath(parentdir, "analysis", modelname, "simulation/orbit.toml")

# ╔═╡ 0dd476fd-be53-4e9b-a686-a4462485c64c
orbit_file = if isfile(joinpath(parentdir, "analysis", modelname, "orbit.csv"))
	joinpath(parentdir, "analysis", modelname, "orbit.csv")
else
	joinpath(parentdir, "analysis", modelname, "simulation/orbit.csv")
end

# ╔═╡ 2bc762ad-e590-443e-b3c2-91dc42a8a4d9
outfile = joinpath(parentdir, "analysis", modelname, "orbital_properties.toml")

# ╔═╡ bb9ef388-bb5a-45a3-836e-4c06dbe0ab65
centresfile = joinpath(parentdir, "analysis", modelname, "centres.hdf5")

# ╔═╡ 30969f77-667e-4ae4-9897-82c1c1182652
md"""
# File loading
"""

# ╔═╡ b250bf10-c228-4b14-938a-35561ae871d7
h5open(centresfile, "r") do  f
	global x_cen, v_cen, t
	x_cen = read(f["positions"])
	v_cen = read(f["velocities"])
	t = read(f["times"])
end

# ╔═╡ bb0cb8c2-2cbd-4205-a39e-4c4c0ff98b8a
begin 
	orbit_expected = CSV.read(orbit_file, DataFrame)
	x_cen_exp = transpose(hcat(orbit_expected.x, orbit_expected.y, orbit_expected.z))
	v_cen_exp = transpose(hcat(orbit_expected.v_x, orbit_expected.v_y, orbit_expected.v_z))
end

# ╔═╡ 08c3df42-738b-47c4-aa6b-fc39a9cfc02f
md"""
# plots
"""

# ╔═╡ a1c992c6-ad12-4968-b105-adfa1f327e76
let
	fig = LilGuys.plot_xyz(x_cen, x_cen_exp, labels=["n body", "point particle"])
	@savefig "centre_xyz"
	fig
end

# ╔═╡ ad1782db-17a2-4385-ab7b-0ce038600a0d
LilGuys.plot_xyz(x_cen, x_cen_exp, labels=["n body", "point particle"], limits=((-200., 200.), (-200., 200.), (-200., 200.)))

# ╔═╡ 5255c605-56ea-4eb3-bd20-5134e3a96705
LilGuys.plot_xyz(v_cen, v_cen_exp, units=" / km/ s")

# ╔═╡ 15293cb8-61d3-478d-a2ae-5a5b2006db44
T2GYR = LilGuys.T2GYR

# ╔═╡ f88b909f-c3dc-41e0-bdb1-25e229964d27
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "radius / kpc")
	r = radii(x_cen)
	lines!(t * T2GYR, r, label="n-body")
	lines!(T2GYR*(orbit_expected.t .- orbit_expected.t[1] .+ t[1]), radii(x_cen_exp),
		label="point particle"
	)

	axislegend(ax)

	@savefig "centre_r_t"
	fig
end

# ╔═╡ f134b3ce-53f0-47e9-84e9-1e73064d5191
snap_cen = Snapshot(x_cen, v_cen, zeros(size(x_cen, 2)))

# ╔═╡ 319b905c-2d08-4a95-9d95-9cd26e2f5b1f
times = t * T2GYR

# ╔═╡ 9530e936-1225-4cfc-aa9a-bf7644d612f5
r = radii(x_cen)

# ╔═╡ 7a30bd90-946e-418c-8339-be64c37cda76
vr = [dot(x_cen[:, i], v_cen[:, i]) / r[i] for i in 1:size(x_cen, 2)]

# ╔═╡ 88a9bb2d-9f86-4a96-a3b8-7a057480c7c3
"""

Given the velocities, finds when the velocity changes sign, 
i.e. a local extrema in r
"""
function find_all_peris(r)
	is_local_max(i) = (r[i] >= r[i-1]) && (r[i] >= r[i+1])
	is_local_min(i) = (r[i] <= r[i-1]) && (r[i] <= r[i+1])
	
	peris = []
	apos = []
	
	for i in 2:length(r)-1
		if is_local_max(i)
			push!(apos, i)
		elseif is_local_min(i)
			push!(peris, i)
		end
	end
	return peris, apos
end

# ╔═╡ 8d1508af-1715-4ef6-aab9-e95a02265913
idx_peris, idx_apos = find_all_peris(r)

# ╔═╡ ddd74bbb-df27-4150-b497-b1df5243f518
r[idx_apos]

# ╔═╡ 04d29fcb-70a0-414b-a487-7a18c44b9d58
let
	fig = Figure(size=(400, 150))
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "radius / kpc")
	r = radii(x_cen)
	lines!(t * T2GYR, r)
	
	scatter!(t[idx_peris] * T2GYR, r[idx_peris], 
		label="pericentrs"
	)

	scatter!(t[idx_apos] * T2GYR, r[idx_apos], 
		label="apocentres"
	)

	Legend(fig[1, 2], ax)

	@savefig "r_t_orbit"

	fig
end

# ╔═╡ 5704daca-a8c4-4292-a6c0-ea294f4373fd
md"""
## Saving
"""

# ╔═╡ 76d5bed6-f1ba-4a2d-8425-ceb40d18abdc
let
	
	orbital_properties = OrderedDict(
		"pericentre" => minimum(r),
		"apocentre" => maximum(r),
		"pericentres" => r[idx_peris],
		"apocentres" => r[idx_apos],
		"idx_peris" => idx_peris, 
		"idx_apos" => idx_apos,
		"t_peris" => t[idx_peris], 
		"t_apos" => t[idx_apos]
	)


	open(outfile, "w") do f
		TOML.print(f, orbital_properties)
	end

	println("saved properties to $outfile")
	orbital_properties
end

# ╔═╡ Cell order:
# ╟─8b41af50-9ae0-475b-bacc-3799e2949b30
# ╠═d94663e8-b30e-4712-8a3e-6ef7f954f141
# ╟─27577252-53dc-415c-b9a1-82155ef9e4ca
# ╠═061b1886-1878-11ef-3806-b91643300982
# ╠═4018f9cb-b54f-408a-8047-cc995a167034
# ╠═cd43d649-1d52-46a7-a621-68c8122fce6d
# ╠═882d4fc5-07ae-4b06-8da5-67f0894595db
# ╠═ab57edae-2292-4aef-9c1f-53802dbc0600
# ╠═2c702eb7-ebb6-44c9-8e01-ca52d011c014
# ╠═6b3df4f3-de5b-40b1-8aec-1cfce9aa843a
# ╠═56d15948-4fd8-4541-9925-99837d9584f5
# ╟─643cd0bf-77b3-4201-9ff7-09dd5aee277c
# ╠═0e2e7f93-09d1-4902-ad24-223df50d37cb
# ╠═3953b76f-a726-4211-a6b8-5cf38149dcdf
# ╠═24ed3bc1-37c6-4601-b941-0780c53a9630
# ╠═079afa70-c7dc-4347-9ae3-8459ba2fa941
# ╠═69d83e00-7eb6-4271-838f-80e4d1654dac
# ╠═d142b7bd-3002-4331-a725-577873c42f28
# ╠═ac2c7484-9acd-4fda-9699-fdf17da507c2
# ╠═0dd476fd-be53-4e9b-a686-a4462485c64c
# ╠═2bc762ad-e590-443e-b3c2-91dc42a8a4d9
# ╠═bb9ef388-bb5a-45a3-836e-4c06dbe0ab65
# ╟─30969f77-667e-4ae4-9897-82c1c1182652
# ╠═b250bf10-c228-4b14-938a-35561ae871d7
# ╠═bb0cb8c2-2cbd-4205-a39e-4c4c0ff98b8a
# ╟─08c3df42-738b-47c4-aa6b-fc39a9cfc02f
# ╠═a1c992c6-ad12-4968-b105-adfa1f327e76
# ╠═ad1782db-17a2-4385-ab7b-0ce038600a0d
# ╠═5255c605-56ea-4eb3-bd20-5134e3a96705
# ╠═15293cb8-61d3-478d-a2ae-5a5b2006db44
# ╠═ddd74bbb-df27-4150-b497-b1df5243f518
# ╠═f88b909f-c3dc-41e0-bdb1-25e229964d27
# ╠═f134b3ce-53f0-47e9-84e9-1e73064d5191
# ╠═319b905c-2d08-4a95-9d95-9cd26e2f5b1f
# ╠═9530e936-1225-4cfc-aa9a-bf7644d612f5
# ╠═7a30bd90-946e-418c-8339-be64c37cda76
# ╠═88a9bb2d-9f86-4a96-a3b8-7a057480c7c3
# ╠═8d1508af-1715-4ef6-aab9-e95a02265913
# ╠═04d29fcb-70a0-414b-a487-7a18c44b9d58
# ╟─5704daca-a8c4-4292-a6c0-ea294f4373fd
# ╠═76d5bed6-f1ba-4a2d-8425-ceb40d18abdc
