### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
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

# ╔═╡ ab57edae-2292-4aef-9c1f-53802dbc0600
import TOML

# ╔═╡ 7834415a-b479-495b-998a-1d12f42f0dc6
Slider = PlutoUI.Slider

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
	haloname = TextField(default="1e6_v31_r3.2"),
	orbitname = TextField(default="orbit_"),
	t_min = NumberField(-10:0.1:10),
))

# ╔═╡ 3953b76f-a726-4211-a6b8-5cf38149dcdf
galaxyname = "idealized"

# ╔═╡ 24ed3bc1-37c6-4601-b941-0780c53a9630
haloname = inputs.haloname

# ╔═╡ 079afa70-c7dc-4347-9ae3-8459ba2fa941
orbitname =  inputs.orbitname

# ╔═╡ 94344455-d1d2-4ef9-af11-2d79ee4729ee
t_min = inputs.t_min

# ╔═╡ dd56b7ec-be11-447f-acc1-12750d82879b
md"""
##### the below should hopefully be always the same
"""

# ╔═╡ bc6deac8-b70a-483b-9fd7-1413c6f17aa7
mc_name = ""

# ╔═╡ 69d83e00-7eb6-4271-838f-80e4d1654dac
modelname = "$galaxyname/$haloname/$orbitname"

# ╔═╡ ac2c7484-9acd-4fda-9699-fdf17da507c2
parentdir = ENV["DWARFS_ROOT"]

# ╔═╡ d142b7bd-3002-4331-a725-577873c42f28
properties_file = joinpath(parentdir, "analysis", modelname, "simulation/orbit.toml")

# ╔═╡ 61d788ec-3518-4e3a-8eef-59c86ae5fc1a
obs_file =  "$parentdir/observations/$galaxyname/observed_properties.toml"

# ╔═╡ 0dd476fd-be53-4e9b-a686-a4462485c64c
orbit_file = joinpath(parentdir, "analysis", modelname, "simulation/orbit.csv")

# ╔═╡ 2bc762ad-e590-443e-b3c2-91dc42a8a4d9
outfile = joinpath(parentdir, "analysis", modelname, "orbital_properties.toml")

# ╔═╡ bb9ef388-bb5a-45a3-836e-4c06dbe0ab65
centresfile = joinpath(parentdir, "analysis", modelname, "centres.hdf5")

# ╔═╡ 9c427388-c657-4bb7-bc0a-b4de3597c645
skyorbit_outfile = joinpath(parentdir, "analysis", modelname, "skyorbit.fits")

# ╔═╡ 4ceac504-5ad2-4cbb-ac15-e094f80ffdbc
FIGDIR = joinpath(parentdir, "analysis", modelname, "figures")

# ╔═╡ b8dec95d-452b-47af-b3f2-2893f970b2f5
if !isdir(FIGDIR)
	mkdir(FIGDIR)
end

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
	orbit_expected.t .-= orbit_expected.t[1]
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

# ╔═╡ 5255c605-56ea-4eb3-bd20-5134e3a96705
LilGuys.plot_xyz(v_cen, v_cen_exp, units=" / km/ s")

# ╔═╡ 15293cb8-61d3-478d-a2ae-5a5b2006db44
T2GYR = LilGuys.T2GYR

# ╔═╡ f88b909f-c3dc-41e0-bdb1-25e229964d27
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "radius / kpc")
	r = calc_r(x_cen)
	lines!(t * T2GYR, r, label="n-body")
	lines!(T2GYR*(orbit_expected.t), calc_r(x_cen_exp),
		label="point particle"
	)

	axislegend(ax)

	@savefig "centre_r_t"
	fig
end

# ╔═╡ 7d29a3bd-dc83-4eb3-ae65-fce5270ed8d5
md"""
# Sky Properties
"""

# ╔═╡ f134b3ce-53f0-47e9-84e9-1e73064d5191
snap_cen = Snapshot(x_cen, v_cen, zeros(size(x_cen, 2)))

# ╔═╡ 64e558da-2928-4815-ad5a-7528516311f9


# ╔═╡ 319b905c-2d08-4a95-9d95-9cd26e2f5b1f
times = t * T2GYR

# ╔═╡ d9df3376-6ca1-4701-afb5-2df994bb3442
idx_f = length(times)

# ╔═╡ aa2c3a93-19a3-43d8-82de-ae6ed8c4b9f7
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="R / kpc", ylabel="z / kpc",
		aspect=DataAspect(),
		xgridvisible=false,
		ygridvisible=false,
	)
	x = x_cen[1, :]
	y = x_cen[2, :]
	z = x_cen[3, :]
	R = @. sqrt(x^2 + y^2)
	lines!(R, z, color=times)

	scatter!(R[idx_f], z[idx_f], color=COLORS[2])
	scatter!(R[1], z[1], color = COLORS[3], marker=:rtriangle, )

	@savefig "R_z_centre"
	fig
end

# ╔═╡ 9530e936-1225-4cfc-aa9a-bf7644d612f5
r = calc_r(x_cen)

# ╔═╡ 6d0612d6-8609-4832-80ae-4e8e78c557cc
minimum(r)

# ╔═╡ 882d4fc5-07ae-4b06-8da5-67f0894595db
import LinearAlgebra: dot

# ╔═╡ 7a30bd90-946e-418c-8339-be64c37cda76
vr = [dot(x_cen[:, i], v_cen[:, i]) / r[i] for i in 1:size(x_cen, 2)]

# ╔═╡ 0c69519c-9650-46b9-89f9-cc37227f5b1a
v_cen

# ╔═╡ d95c457b-c9be-4570-bc90-b4bbb7de56e2
plot(t*T2GYR, vr * V2KMS, 
	axis = (;
	xlabel = "time / Gyr",
	ylabel = "galactocentric radial velocity / km/s"
	)
)

# ╔═╡ 88a9bb2d-9f86-4a96-a3b8-7a057480c7c3
"""

Given the velocities, finds when the velocity changes sign, 
i.e. a local extrema in r
"""
function find_all_peris(r)
	is_local_max(i) = (r[i] >= r[i-1]) && (r[i] >= r[i+1])
	is_local_min(i) = (r[i] <= r[i-1]) && (r[i] <= r[i+1])
	
	peris = Int[]
	apos = Int[]
	
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

# ╔═╡ efcbae60-cf7c-4e74-aae4-39d19b74b6fa
idx_peri = maximum(idx_peris)

# ╔═╡ d81c1455-728a-4023-ad65-e3cce37a69f9
r[idx_peris[end]], minimum(r)

# ╔═╡ a5ce5442-73ca-4aaf-915a-72fe9936e791
d_idx = 20

# ╔═╡ 7e14d1c7-b37f-4c0c-8e2a-1ac7cce27aaa
begin
	peri_filt = idx_f-d_idx:idx_f
	t_last_peri_arg = argmin(r[peri_filt])
	t_last_peri = t[peri_filt[t_last_peri_arg]] * T2GYR
	delta_t_peri = t[idx_f] * T2GYR - t_last_peri
end

# ╔═╡ f64dcd49-b0ca-4319-b615-5520b23d7818
lcm

# ╔═╡ 04d29fcb-70a0-414b-a487-7a18c44b9d58
let
	fig = Figure(size=(700, 300))
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "radius / kpc")
	r = calc_r(x_cen)
	lines!(t * T2GYR, r)
	scatter!(t[idx_f] * T2GYR, r[idx_f], 
		label="adpoted end", marker=:rect
	)
	
	scatter!(t[idx_peris] * T2GYR, r[idx_peris], 
		label="pericentrs"
	)

	scatter!(t[idx_apos] * T2GYR, r[idx_apos], 
		label="apocentres"
	)
	
	scatter!(t[idx_peri] * T2GYR, r[idx_peri], 
		label=" last pericentre"
	)
	
	Legend(fig[1, 2], ax)

	@savefig "r_t_orbit"

	fig
end

# ╔═╡ af8a50bd-e761-4439-9fc9-80048c264d5b
begin 
	if idx_peri > 0
		t_peri = T2GYR * t[idx_peri]
		r_peri = r[idx_peri]

	else 
		t_peri = NaN
		r_peri = NaN
	end

end

# ╔═╡ 73bb2d61-37f3-4782-ae89-d36d1ff8f9ff
t_f = t[idx_f] * T2GYR

# ╔═╡ 5676b72f-a981-4efb-8328-c25d3c5d6fb0
periods = [diff(times[idx_peris]); diff(times[idx_apos])]

# ╔═╡ 592a18e9-ee9a-4638-9454-f0bda3a0a3f2
period = LilGuys.mean(periods)

# ╔═╡ 179c3c32-1368-4a58-b4b8-26d9d3f19f8c
md"""
## Identification of interesting locations to search....
- for this project, want to find stars between 2 and 10 degrees away along stream path.
"""

# ╔═╡ 5704daca-a8c4-4292-a6c0-ea294f4373fd
md"""
## Saving
"""

# ╔═╡ 76d5bed6-f1ba-4a2d-8425-ceb40d18abdc
let
	
	orbital_properties = Dict(
		"pericentre" => r_peri,
		"apocentre" => r[idx_apos[idx_apos .< idx_f]],
		"period" => period,
		"idx_f" => idx_f,
		"r_f" => r[idx_f],
		"idx_peri" => idx_peri,
		"t_last_peri" => t_f - t_peri,
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
# ╠═ab57edae-2292-4aef-9c1f-53802dbc0600
# ╠═2c702eb7-ebb6-44c9-8e01-ca52d011c014
# ╠═6b3df4f3-de5b-40b1-8aec-1cfce9aa843a
# ╠═7834415a-b479-495b-998a-1d12f42f0dc6
# ╟─643cd0bf-77b3-4201-9ff7-09dd5aee277c
# ╠═0e2e7f93-09d1-4902-ad24-223df50d37cb
# ╠═3953b76f-a726-4211-a6b8-5cf38149dcdf
# ╠═24ed3bc1-37c6-4601-b941-0780c53a9630
# ╠═079afa70-c7dc-4347-9ae3-8459ba2fa941
# ╠═94344455-d1d2-4ef9-af11-2d79ee4729ee
# ╟─dd56b7ec-be11-447f-acc1-12750d82879b
# ╠═bc6deac8-b70a-483b-9fd7-1413c6f17aa7
# ╠═69d83e00-7eb6-4271-838f-80e4d1654dac
# ╠═d142b7bd-3002-4331-a725-577873c42f28
# ╠═61d788ec-3518-4e3a-8eef-59c86ae5fc1a
# ╠═ac2c7484-9acd-4fda-9699-fdf17da507c2
# ╠═0dd476fd-be53-4e9b-a686-a4462485c64c
# ╠═2bc762ad-e590-443e-b3c2-91dc42a8a4d9
# ╠═bb9ef388-bb5a-45a3-836e-4c06dbe0ab65
# ╠═9c427388-c657-4bb7-bc0a-b4de3597c645
# ╠═4ceac504-5ad2-4cbb-ac15-e094f80ffdbc
# ╠═b8dec95d-452b-47af-b3f2-2893f970b2f5
# ╟─30969f77-667e-4ae4-9897-82c1c1182652
# ╠═b250bf10-c228-4b14-938a-35561ae871d7
# ╠═bb0cb8c2-2cbd-4205-a39e-4c4c0ff98b8a
# ╟─08c3df42-738b-47c4-aa6b-fc39a9cfc02f
# ╠═a1c992c6-ad12-4968-b105-adfa1f327e76
# ╠═5255c605-56ea-4eb3-bd20-5134e3a96705
# ╠═aa2c3a93-19a3-43d8-82de-ae6ed8c4b9f7
# ╠═15293cb8-61d3-478d-a2ae-5a5b2006db44
# ╠═ddd74bbb-df27-4150-b497-b1df5243f518
# ╠═f88b909f-c3dc-41e0-bdb1-25e229964d27
# ╠═6d0612d6-8609-4832-80ae-4e8e78c557cc
# ╟─7d29a3bd-dc83-4eb3-ae65-fce5270ed8d5
# ╠═f134b3ce-53f0-47e9-84e9-1e73064d5191
# ╠═64e558da-2928-4815-ad5a-7528516311f9
# ╠═319b905c-2d08-4a95-9d95-9cd26e2f5b1f
# ╠═d9df3376-6ca1-4701-afb5-2df994bb3442
# ╠═9530e936-1225-4cfc-aa9a-bf7644d612f5
# ╠═882d4fc5-07ae-4b06-8da5-67f0894595db
# ╠═7a30bd90-946e-418c-8339-be64c37cda76
# ╠═0c69519c-9650-46b9-89f9-cc37227f5b1a
# ╟─d95c457b-c9be-4570-bc90-b4bbb7de56e2
# ╠═88a9bb2d-9f86-4a96-a3b8-7a057480c7c3
# ╠═8d1508af-1715-4ef6-aab9-e95a02265913
# ╠═efcbae60-cf7c-4e74-aae4-39d19b74b6fa
# ╠═d81c1455-728a-4023-ad65-e3cce37a69f9
# ╠═a5ce5442-73ca-4aaf-915a-72fe9936e791
# ╠═7e14d1c7-b37f-4c0c-8e2a-1ac7cce27aaa
# ╠═f64dcd49-b0ca-4319-b615-5520b23d7818
# ╠═04d29fcb-70a0-414b-a487-7a18c44b9d58
# ╠═af8a50bd-e761-4439-9fc9-80048c264d5b
# ╠═73bb2d61-37f3-4782-ae89-d36d1ff8f9ff
# ╠═5676b72f-a981-4efb-8328-c25d3c5d6fb0
# ╠═592a18e9-ee9a-4638-9454-f0bda3a0a3f2
# ╟─179c3c32-1368-4a58-b4b8-26d9d3f19f8c
# ╠═5704daca-a8c4-4292-a6c0-ea294f4373fd
# ╠═76d5bed6-f1ba-4a2d-8425-ceb40d18abdc
