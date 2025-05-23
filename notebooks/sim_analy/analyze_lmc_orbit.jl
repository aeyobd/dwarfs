### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

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

# ╔═╡ 3f1890df-eea5-4903-b991-b68dd6b85cbf
using OrderedCollections

# ╔═╡ 8b41af50-9ae0-475b-bacc-3799e2949b30
md"""
Analyzes the orbit of a n-body halo in a gravitational potential.
Requires the centres to be calculated prior.

This notebook only uses the centres of the orbit, but then calculates useful quantities such as the pericentre, the orbital period, the time of each peri and apocentre and outputs this into the model analysis directory.

"""

# ╔═╡ ab57edae-2292-4aef-9c1f-53802dbc0600
import TOML

# ╔═╡ 643cd0bf-77b3-4201-9ff7-09dd5aee277c
md"""
# inputs
"""

# ╔═╡ 597cd4fe-8e50-4ab9-b933-c7d0a2c0a156
galaxyname = "sculptor"

# ╔═╡ 69d83e00-7eb6-4271-838f-80e4d1654dac
modelname = "$galaxyname/1e7_V31_r4.2/vasiliev24_L3M11_2x_smallperilmc"

# ╔═╡ dd56b7ec-be11-447f-acc1-12750d82879b
md"""
the below should hopefully be always the same
"""

# ╔═╡ ac2c7484-9acd-4fda-9699-fdf17da507c2
parentdir = ENV["DWARFS_ROOT"]

# ╔═╡ d142b7bd-3002-4331-a725-577873c42f28
properties_file =  "$parentdir/observations/$galaxyname/observed_properties.toml"

# ╔═╡ 0dd476fd-be53-4e9b-a686-a4462485c64c
orbit_file = joinpath(parentdir, "analysis", modelname, "simulation/orbit.csv")

# ╔═╡ 460dfdff-7f5f-46bc-bd5c-4a43d916d157
lmc_file = joinpath(parentdir, "analysis", modelname, "lmc_traj.csv")

# ╔═╡ 3711412d-ae30-449b-8671-88fe0d128d20
lmc_exp_file = joinpath(parentdir, "analysis", modelname, "lmc_traj_exp.csv")

# ╔═╡ 48a8bc72-fd81-44f5-abab-9afded7df15f
modeldir = joinpath(parentdir, "analysis", modelname)

# ╔═╡ 2bc762ad-e590-443e-b3c2-91dc42a8a4d9
outfile = joinpath(parentdir, "analysis", modelname, "orbital_properties_lmc.toml")

# ╔═╡ bb9ef388-bb5a-45a3-836e-4c06dbe0ab65
centresfile = joinpath(parentdir, "analysis", modelname, "centres.hdf5")

# ╔═╡ 9c427388-c657-4bb7-bc0a-b4de3597c645
skyorbit_outfile = joinpath(parentdir, "analysis", modelname, "skyorbit.fits")

# ╔═╡ 63d28179-452f-42e3-9ded-938efa7bb144


# ╔═╡ 30969f77-667e-4ae4-9897-82c1c1182652
md"""
# File loading
"""

# ╔═╡ 96a57df5-a7b7-447a-a4a6-2b05e391a5c6
obs_today = TOML.parsefile(properties_file)

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

# ╔═╡ fdf491f5-a540-4fd7-8092-3040695d4f7e
begin 
	orbit_lmc = CSV.read(lmc_file, DataFrame)
	x_lmc = transpose(hcat(orbit_lmc.x, orbit_lmc.y, orbit_lmc.z))
	v_lmc = transpose(hcat(orbit_lmc.v_x, orbit_lmc.v_y, orbit_lmc.v_z))
end

# ╔═╡ 1b16b7ab-fbad-40d9-8dfe-7980410a8363
begin 
	orbit_lmc_exp = CSV.read(lmc_exp_file, DataFrame)
	x_lmc_exp = transpose(hcat(orbit_lmc_exp.x, orbit_lmc_exp.y, orbit_lmc_exp.z))
	v_lmc_exp = transpose(hcat(orbit_lmc_exp.v_x, orbit_lmc_exp.v_y, orbit_lmc_exp.v_z))
end

# ╔═╡ ae315001-c5a5-4517-8fd3-97c9cd7c36ec
@assert orbit_lmc.time ≈ t

# ╔═╡ 0deba830-429b-4097-9bc9-b5195dd21e30
x_scl_lmc = x_cen .- x_lmc

# ╔═╡ 1e1b3641-6454-4f9b-b96a-336406dab56a
v_scl_lmc = v_cen .- v_lmc

# ╔═╡ 06e931d8-5518-4594-818c-4dccba83b4c1
x_scl_lmc_exp = x_cen_exp .- x_lmc_exp

# ╔═╡ 08c3df42-738b-47c4-aa6b-fc39a9cfc02f
md"""
# plots
"""

# ╔═╡ 78b1e507-14c3-4279-b502-d92ed5bfb27b
idx_interest = [1, argmin(radii(x_scl_lmc))]

# ╔═╡ 8535eb24-a43b-4680-a54a-a7bddb535999
x_cen[:, idx_interest]

# ╔═╡ d3738bfb-0fdb-4a02-a186-8fc4180bbe35
x_lmc[:, idx_interest]

# ╔═╡ d6086a06-0d65-42b4-a125-7b64dc167197
x_lmc[:, end]

# ╔═╡ 3410ae62-b91f-42c1-a8ee-9445565923d3
LilGuys.radii((v_cen[:, idx_interest] .- v_lmc[:, idx_interest]) * V2KMS)

# ╔═╡ 79147e9e-f09b-4e05-882f-367b31cbcef2


# ╔═╡ a1c992c6-ad12-4968-b105-adfa1f327e76
LilGuys.plot_xyz(x_cen, x_cen_exp, x_lmc, labels=["n body", "point particle", "LMC"])

# ╔═╡ b4c12b61-6827-4c63-bbbc-e708fdb8ba4b
LilGuys.plot_xyz(x_cen, x_lmc, labels=["Scl", "LMC"], idx_scatter=fill(idx_interest, 2))

# ╔═╡ c96035d2-8860-4609-96d3-116f0928fa14
LilGuys.struct_to_dict((; a=1))

# ╔═╡ e500cd04-30be-4de2-890b-cc7f460cc176
LilGuys.plot_xyz(x_scl_lmc, x_scl_lmc_exp, labels=["Scl - LMC", "point particle"])

# ╔═╡ 5255c605-56ea-4eb3-bd20-5134e3a96705
LilGuys.plot_xyz(v_cen, v_cen_exp, v_lmc, units=" / km/ s")

# ╔═╡ 15293cb8-61d3-478d-a2ae-5a5b2006db44
T2GYR = LilGuys.T2GYR

# ╔═╡ f88b909f-c3dc-41e0-bdb1-25e229964d27
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "radius / kpc")
	r = radii(x_cen)
	lines!(t * T2GYR, r, label="model")
	lines!(T2GYR*(orbit_expected.t), radii(x_cen_exp),
		label="expected"
	)

	axislegend(ax)
	fig
end

# ╔═╡ 01b57b06-c658-49e1-b61d-cf1c7e223723
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "Scl - LMC distance / kpc",
		limits=(nothing, nothing, 0, nothing)
	)
	r = radii(x_scl_lmc)
	lines!(t * T2GYR, r, label="nbody")
	@info minimum(r)

	r = radii(x_scl_lmc_exp)
	lines!(orbit_expected.t * T2GYR, r, label="point particle")
	axislegend(ax, position=:lb)
	fig
end

# ╔═╡ 319b905c-2d08-4a95-9d95-9cd26e2f5b1f
times = t * T2GYR

# ╔═╡ 7bbbc936-7626-423b-bcaa-9870f320c432
times

# ╔═╡ 88ae0d01-2a5d-4434-8452-5f89ba9d55c3
idx_gyr = [argmin(abs.(times .- i )) for i in -5:0]

# ╔═╡ a82343a4-e41b-44a0-9f3f-b0fccc0bc2b2
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "Scl - LMC distance / kpc",
		limits=(nothing, nothing, 0, nothing)
	)
	r = radii(x_scl_lmc)
	lines!(t * T2GYR, r, label="nbody")
	scatter!(t[idx_gyr] * T2GYR, r[idx_gyr], label="nbody")

	fig
end

# ╔═╡ 8a17fad5-c41b-48d6-be30-ef52973880d1
times[idx_interest]

# ╔═╡ eb7a269a-2a20-4ea7-b52c-b74277df5617
times[idx_interest]

# ╔═╡ 761750fd-4b18-4c5f-8b79-0f87307e08a1
LilGuys.plot_xyz(x_scl_lmc, times=times / T2GYR, idx_scatter=[idx_gyr])

# ╔═╡ b7c3051f-7aec-472f-906d-c29b8e078475
LilGuys.plot_xyz(v_scl_lmc, label="scl - lmc", units = " / km / s", times=times / T2GYR, idx_scatter=[idx_gyr])

# ╔═╡ a8e051d1-39dd-444b-85de-f235339efa26
LilGuys.plot_xyz(v_cen, label="scl - lmc", units = " / km / s", times=times / T2GYR, idx_scatter=[idx_gyr])

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

	scatter!(R[1], z[1], color = COLORS[3], marker=:rtriangle, )
	fig
end

# ╔═╡ cb4241fe-30d0-4a27-bc90-e3bdf15aa32e
length(times)

# ╔═╡ 9530e936-1225-4cfc-aa9a-bf7644d612f5
r = radii(x_cen)

# ╔═╡ 76a51188-2168-4cc6-b9e2-ede9f6ec34e4
times[argmin(r)]

# ╔═╡ 4d696d88-216c-40e9-9159-2e23441b61d3
r_lmc = radii(x_scl_lmc)

# ╔═╡ cae63625-df19-4efc-95f2-58887e8fb759
md"""
For a specific mode
"""

# ╔═╡ 53d482ae-1e74-4b33-9392-90bc110c113d
idx_f_2 = argmin(abs.(times .- -0.25))

# ╔═╡ 335763ff-bca3-406c-b2a9-7769260cabad
times[idx_f_2]

# ╔═╡ 9d901528-b95c-4d91-899e-033ae68ae73d
times[argmin(r)] - times[idx_f_2]

# ╔═╡ d574d14c-c48b-4439-b533-7352c2858b8a
radii(x_cen[:, idx_f_2] .- [8.122, 0, 0])

# ╔═╡ b45cced1-51ef-432b-a800-24ba5d543f2a
radii(x_scl_lmc[:, idx_f_2])

# ╔═╡ 619a9aae-dcdc-4264-8086-c8149387ed6b
md"""
# Saving
"""

# ╔═╡ 35a185f0-ef91-4db0-abba-607005b6bc1a
idx_peri = argmin(r_lmc)

# ╔═╡ 4ac8a97d-9142-42e1-9826-7b81a5b9a28a
t_peri = times[idx_peri]

# ╔═╡ c289ee03-5926-4ce0-80e4-02f297f0636c
r_peri = r_lmc[idx_peri]

# ╔═╡ 6ef4980c-8987-4263-bedf-273bb065724b
delta_v_peri = radii(v_scl_lmc)[idx_peri] * V2KMS

# ╔═╡ 8f8a9115-d7bf-4ce5-9132-bb6dbd71c8c7
orbital_props = TOML.parsefile(modeldir * "/orbital_properties.toml")

# ╔═╡ e901b8c7-9989-4dd9-9298-420ab96d78f1
t_f = t[orbital_props["idx_f"]]

# ╔═╡ 38e23f2e-962d-4edb-90c8-4a9e16b1d793
let
	
	orbital_properties = OrderedDict(
		"pericentre" => r_peri,
		"idx_peri" => idx_peri,
		"t_last_peri" => t_f - t_peri,
		"delta_v_peri" => delta_v_peri
	)


	open(outfile, "w") do f
		TOML.print(f, orbital_properties)
	end

	println("saved properties to $outfile")
	orbital_properties
end

# ╔═╡ Cell order:
# ╟─8b41af50-9ae0-475b-bacc-3799e2949b30
# ╠═061b1886-1878-11ef-3806-b91643300982
# ╠═ab57edae-2292-4aef-9c1f-53802dbc0600
# ╟─643cd0bf-77b3-4201-9ff7-09dd5aee277c
# ╠═597cd4fe-8e50-4ab9-b933-c7d0a2c0a156
# ╠═69d83e00-7eb6-4271-838f-80e4d1654dac
# ╟─dd56b7ec-be11-447f-acc1-12750d82879b
# ╠═ac2c7484-9acd-4fda-9699-fdf17da507c2
# ╠═d142b7bd-3002-4331-a725-577873c42f28
# ╠═0dd476fd-be53-4e9b-a686-a4462485c64c
# ╠═460dfdff-7f5f-46bc-bd5c-4a43d916d157
# ╠═3711412d-ae30-449b-8671-88fe0d128d20
# ╠═48a8bc72-fd81-44f5-abab-9afded7df15f
# ╠═2bc762ad-e590-443e-b3c2-91dc42a8a4d9
# ╠═bb9ef388-bb5a-45a3-836e-4c06dbe0ab65
# ╠═9c427388-c657-4bb7-bc0a-b4de3597c645
# ╠═63d28179-452f-42e3-9ded-938efa7bb144
# ╟─30969f77-667e-4ae4-9897-82c1c1182652
# ╠═96a57df5-a7b7-447a-a4a6-2b05e391a5c6
# ╠═2c702eb7-ebb6-44c9-8e01-ca52d011c014
# ╠═b250bf10-c228-4b14-938a-35561ae871d7
# ╠═bb0cb8c2-2cbd-4205-a39e-4c4c0ff98b8a
# ╠═fdf491f5-a540-4fd7-8092-3040695d4f7e
# ╠═1b16b7ab-fbad-40d9-8dfe-7980410a8363
# ╠═ae315001-c5a5-4517-8fd3-97c9cd7c36ec
# ╠═0deba830-429b-4097-9bc9-b5195dd21e30
# ╠═1e1b3641-6454-4f9b-b96a-336406dab56a
# ╠═06e931d8-5518-4594-818c-4dccba83b4c1
# ╟─08c3df42-738b-47c4-aa6b-fc39a9cfc02f
# ╠═7bbbc936-7626-423b-bcaa-9870f320c432
# ╠═88ae0d01-2a5d-4434-8452-5f89ba9d55c3
# ╠═78b1e507-14c3-4279-b502-d92ed5bfb27b
# ╠═8a17fad5-c41b-48d6-be30-ef52973880d1
# ╠═8535eb24-a43b-4680-a54a-a7bddb535999
# ╠═d3738bfb-0fdb-4a02-a186-8fc4180bbe35
# ╠═d6086a06-0d65-42b4-a125-7b64dc167197
# ╠═3410ae62-b91f-42c1-a8ee-9445565923d3
# ╠═79147e9e-f09b-4e05-882f-367b31cbcef2
# ╠═eb7a269a-2a20-4ea7-b52c-b74277df5617
# ╠═a1c992c6-ad12-4968-b105-adfa1f327e76
# ╠═b4c12b61-6827-4c63-bbbc-e708fdb8ba4b
# ╠═761750fd-4b18-4c5f-8b79-0f87307e08a1
# ╠═c96035d2-8860-4609-96d3-116f0928fa14
# ╠═e500cd04-30be-4de2-890b-cc7f460cc176
# ╠═5255c605-56ea-4eb3-bd20-5134e3a96705
# ╠═b7c3051f-7aec-472f-906d-c29b8e078475
# ╠═a8e051d1-39dd-444b-85de-f235339efa26
# ╠═aa2c3a93-19a3-43d8-82de-ae6ed8c4b9f7
# ╠═15293cb8-61d3-478d-a2ae-5a5b2006db44
# ╠═f88b909f-c3dc-41e0-bdb1-25e229964d27
# ╠═01b57b06-c658-49e1-b61d-cf1c7e223723
# ╠═a82343a4-e41b-44a0-9f3f-b0fccc0bc2b2
# ╠═319b905c-2d08-4a95-9d95-9cd26e2f5b1f
# ╠═76a51188-2168-4cc6-b9e2-ede9f6ec34e4
# ╠═cb4241fe-30d0-4a27-bc90-e3bdf15aa32e
# ╠═9530e936-1225-4cfc-aa9a-bf7644d612f5
# ╠═4d696d88-216c-40e9-9159-2e23441b61d3
# ╠═cae63625-df19-4efc-95f2-58887e8fb759
# ╠═53d482ae-1e74-4b33-9392-90bc110c113d
# ╠═335763ff-bca3-406c-b2a9-7769260cabad
# ╠═9d901528-b95c-4d91-899e-033ae68ae73d
# ╠═d574d14c-c48b-4439-b533-7352c2858b8a
# ╠═b45cced1-51ef-432b-a800-24ba5d543f2a
# ╠═619a9aae-dcdc-4264-8086-c8149387ed6b
# ╠═3f1890df-eea5-4903-b991-b68dd6b85cbf
# ╠═35a185f0-ef91-4db0-abba-607005b6bc1a
# ╠═4ac8a97d-9142-42e1-9826-7b81a5b9a28a
# ╠═c289ee03-5926-4ce0-80e4-02f297f0636c
# ╠═6ef4980c-8987-4263-bedf-273bb065724b
# ╠═8f8a9115-d7bf-4ce5-9132-bb6dbd71c8c7
# ╠═e901b8c7-9989-4dd9-9298-420ab96d78f1
# ╠═38e23f2e-962d-4edb-90c8-4a9e16b1d793
