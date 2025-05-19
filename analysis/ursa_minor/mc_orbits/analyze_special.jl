### A Pluto.jl notebook ###
# v0.20.8

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

# ╔═╡ e9e2c787-4e0e-4169-a4a3-401fea21baba
begin 
	import Pkg
	Pkg.activate()
	using CairoMakie
	using DataFrames, CSV

	using Printf
	import LilGuys as lguys

	using Arya
	using PlutoUI
end


# ╔═╡ 91cac7f1-e3a7-4f4e-a136-d99cbfef75ed
using PyFITS

# ╔═╡ d975d00c-fd69-4dd0-90d4-c4cbe73d9754
using Statistics, Distributions

# ╔═╡ 1be546d7-8b11-4dc5-9295-393d39f65123
using OrderedCollections

# ╔═╡ 83fca0a6-5693-4b06-89ff-2d097231c402
include("./mc_analysis_utils.jl")

# ╔═╡ 7450144e-5464-4036-a215-b6e2cd270405
md"""
This notebook analyzes the result of the MC samples of orbits in the same potential to determine the plausable range of pericentres and apocentres
"""

# ╔═╡ baf5a0e0-3ee1-4056-b3a5-15058b4a0e53
@bind modelname TextField(default="EP2020")

# ╔═╡ abac7f99-fc5d-4536-a845-52cc269deee7
FIGDIR = "./$(modelname)_special_cases/figures"

# ╔═╡ b3e118f4-0503-4a6c-8f0d-d096bd0236a3
using LilGuys; FIGDIR

# ╔═╡ 2559bc4c-3bad-4d70-b02f-8fe48b772fb6
modelname_special = "$(modelname)_special_cases"

# ╔═╡ e5f728b8-8412-4f57-ad38-a0a35bb08a48
orbit_labels = ["mean", "smallperi", "largeperi"]

# ╔═╡ 3be492fa-0cb9-40f1-a39e-e25bf485fd3e
md"""
# Setup
"""

# ╔═╡ 3b83205d-91c1-481e-9305-0d59bc692135
coord_labels = Dict(
	:ra => "ra / degrees",
	:dec => "dec / degrees",
	:pmra => L"$\mu_{\alpha *}$ / mas\,yr$^{-1}$",
	:pmdec => L"$\mu_{\delta}$ / mas\,yr$^{-1}$",
	:radial_velocity => L"$v_\textrm{los}$ / km\,s$^{-1}$",
	:distance => "distance / kpc",
)

# ╔═╡ 88536e86-cf2a-4dff-ae64-514821957d40
md"""
warning, the simulation output and peris_apos.fits must analyze the same simulation. Checking the individual orbits below should confirm this and the loading code should fail if there are any issues.
"""

# ╔═╡ 26d616da-95ec-4fb9-b9a8-2f095d74c722
"""
	sort_snap(snap)

returns the snapshot sorted by the index
"""
function sort_snap(snap)
	return snap[sortperm(snap.index)]
end

# ╔═╡ 38fd9158-7807-4a26-b534-b0bcb5dbb0e5
md"""
# Data Loading
"""

# ╔═╡ 9a22d47b-8474-4596-b418-de33eb07c627
begin 
	out = lguys.Output("$modelname");

	df_peris_apos = read_fits("$modelname/peris_apos.fits")
	snap = out[1] |> sort_snap

	@assert all(snap.index  .== df_peris_apos.index) "snapshot and peri apo index must match"
end

# ╔═╡ dcdca4a7-08a6-4d1b-a972-c386493207d0
begin 
	out_special = lguys.Output(modelname_special);

	df_peris_apos_special = read_fits("$modelname_special/peris_apos.fits")
	snap_special = out_special[1] |> sort_snap

	@assert all(snap_special.index  .== df_peris_apos_special.index) "snapshot and peri apo index must match"
end

# ╔═╡ 56631a4d-323e-42ac-ad71-3e35face895d
obs_props_filename = ENV["DWARFS_ROOT"] * 
	"/observations/ursa_minor/observed_properties.toml"


# ╔═╡ c60b4d12-642f-4247-a022-f2c0bb16a092
delta_t = df_peris_apos.t_last_peri

# ╔═╡ bbf3f229-fc3c-46ae-af28-0f8bd81e7d32
df_peris_apos_special.pericentre

# ╔═╡ 4481a5e1-d635-4fea-a5c5-c85f0d6df62f
peris = df_peris_apos.pericentre

# ╔═╡ 384be6a6-f9d9-47e0-9792-aef6689dcbdb
apos = df_peris_apos.apocentre

# ╔═╡ 1acef60e-60d6-47ba-85fd-f9780934788b
md"""
# plots
"""

# ╔═╡ ca1c236e-795a-408b-845b-9c13bc838619
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "pericentre / kpc",
		ylabel = "count"
	)

	bins, counts, err = lguys.histogram(peris)
	scatter!(lguys.midpoints(bins), counts)
	vlines!(df_peris_apos_special.pericentre, color=COLORS[1:3])

	fig
end

# ╔═╡ 46b4242b-8af7-4233-8ecf-d86740b4c884
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "apocentre / kpc",
		ylabel = "count"
	)

	bins, counts, err = lguys.histogram(apos)
	scatter!(lguys.midpoints(bins), counts)
	vlines!(df_peris_apos_special.apocentre, color=COLORS[1:3])

	fig
end

# ╔═╡ 377628e5-67a7-45fa-9e98-992acf7171ea
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "time of last peri / Gyr",
		ylabel = "count"
	)

	bins, counts, err = lguys.histogram(df_peris_apos.t_last_peri * lguys.T2GYR)
	scatter!(lguys.midpoints(bins), counts)
	vlines!(df_peris_apos_special.t_last_peri * lguys.T2GYR, color=COLORS[1:3])

	fig
end

# ╔═╡ 2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
ϵ = -lguys.specific_kinetic_energy(out[1]) .- out[1].potential_ext

# ╔═╡ 8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
extrema(ϵ)

# ╔═╡ 5f11f6ab-c9ab-4600-acca-a0bb84d81a12
begin
	points = [lguys.Galactocentric(
		snap.positions[:, i]*lguys.R2KPC, 
		-snap.velocities[:, i]*lguys.V2KMS,
	)
		for i in 1:length(snap)]
	
	observations = lguys.transform.(lguys.ICRS, points)
end


# ╔═╡ 92aac8e8-d599-4a1e-965e-e653bc54c509
dists = getproperty.(observations, :distance)

# ╔═╡ 3e18be79-9403-4104-9d4f-9860163ad8b9
begin
	points_special = [lguys.Galactocentric(
		snap_special.positions[:, i]*lguys.R2KPC, 
		-snap_special.velocities[:, i]*lguys.V2KMS,
	)
		for i in 1:length(snap_special)]
	
	observations_special = lguys.transform.(lguys.ICRS, points_special)
end


# ╔═╡ ede3836c-740d-4ac7-bbc7-3165981a1878
normal_dist(x, μ, σ) = 1/√(2π) * 1/σ * exp(-(x-μ)^2/2σ^2)

# ╔═╡ 44660b2f-6220-473b-bb2f-07e23b176491
columns = [:ra, :dec, :distance, :pmra, :pmdec, :radial_velocity]

# ╔═╡ 67ca9626-8416-4254-87f8-d1c5c88bf2de
orbit_points_kwargs = [Dict(
	:color => COLORS[i],
	:alpha => 1,
	:markersize => 10,
	:label => orbit_labels[i]
) for i in 1:length(snap_special)
]

# ╔═╡ 8a863d67-2535-4117-ae19-437b3e4c6c71
observations_special

# ╔═╡ c48b4e73-480e-4a50-b5fc-db5f6c5b040e
let
	fig = Figure(size=(600, 600))
	plot_kwargs = Dict(
		:color => :black,
		:alpha => 0.1,
		:markersize => 1,
	)

	ax_kwargs = Dict(
		:xgridvisible => false,
		:ygridvisible => false,
		:ylabel => "pericentre / kpc",
	)

	ax_idx = Dict(
		:pmra => [1, 1],
		:pmdec => [1, 2],
		:distance => [2, 1],
		:radial_velocity => [2, 2],
	)

	for sym in [:pmra, :pmdec, :distance, :radial_velocity]
		ax = Axis(fig[ax_idx[sym]...],
			xlabel=coord_labels[sym];
			ax_kwargs...
		)
		x = [getproperty(o, sym) for o in observations]
		scatter!(x, peris; rasterize=true, plot_kwargs...)
		
		for i in eachindex(orbit_labels)
				x = getproperty(observations_special[i], sym)
				y = df_peris_apos_special[i, :pericentre]
				scatter!(x, y; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)


	@savefig "peri_mc_orbits_corr"
	fig
end

# ╔═╡ 43d43f63-4c13-4b23-950e-ada59aa86bc9
let
	fig = Figure(size=(600, 600))
	plot_kwargs = Dict(
		:color => :black,
		:alpha => 0.1,
		:markersize => 1,
	)

	ax_kwargs = Dict(
		:xgridvisible => false,
		:ygridvisible => false,
		:ylabel => "apocenter / kpc",
	)

	ax_idx = Dict(
		:pmra => [1, 1],
		:pmdec => [1, 2],
		:distance => [2, 1],
		:radial_velocity => [2, 2],
	)

	for sym in [:pmra, :pmdec, :distance, :radial_velocity]
		ax = Axis(fig[ax_idx[sym]...],
			xlabel=coord_labels[sym];
			ax_kwargs...
		)
		x = [getproperty(o, sym) for o in observations]
		scatter!(x, apos; plot_kwargs...)
		
		for i in eachindex(orbit_labels)
				x = getproperty(observations_special[i], sym)
				y = df_peris_apos_special[i, :apocentre]
				scatter!(x, y; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)


	@savefig "apo_mc_orbits_corr"
	fig
end

# ╔═╡ f14fb67d-37c1-46f8-9d63-3f0d4fd93b3c
let
	fig = Figure(size=(600, 600))
	plot_kwargs = Dict(
		:color => :black,
		:alpha => 0.1,
		:markersize => 1,
	)

	ax_kwargs = Dict(
		:xgridvisible => false,
		:ygridvisible => false,
		:ylabel => L"$\Delta\,t_\textrm{peri}$ / Gyr",
	)

	ax_idx = Dict(
		:pmra => [1, 1],
		:pmdec => [1, 2],
		:distance => [2, 1],
		:radial_velocity => [2, 2],
	)

	for sym in [:pmra, :pmdec, :distance, :radial_velocity]
		ax = Axis(fig[ax_idx[sym]...],
			xlabel=coord_labels[sym];
			ax_kwargs...
		)
		x = [getproperty(o, sym) for o in observations]
		scatter!(x, df_peris_apos.t_last_peri * lguys.T2GYR .+ 0.001 .* randn(length(snap)); plot_kwargs...)
				for i in eachindex(orbit_labels)
				x = getproperty(observations_special[i], sym)
				y = lguys.T2GYR * df_peris_apos_special[i, :t_last_peri]
				scatter!(x, y; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)


	@savefig "t_last_peri"
	fig
end

# ╔═╡ ac81acd8-4a78-4230-bc70-3b78a861b618
let

	for sym in columns
		
		fig = Figure()
		ax = Axis(fig[1,1], 
			xlabel=String(sym),
			ylabel="density",
			#limits=((μ - 5σ, μ + 5σ), nothing),
		)
		
	    x = getproperty.(observations, sym)
		
	    stephist!(x, bins=50, normalization=:pdf, label="MC samples", color=:black)

	    
	    x_mod = LinRange(minimum(x), maximum(x), 1000)

		vlines!(getproperty.(observations_special, sym), color=COLORS[1:length(orbit_labels)])
		axislegend()

		@info fig
	end

end

# ╔═╡ 16f4ac20-d8cf-4218-8c01-c15e04e567fb
md"""
# The selected orbits
"""

# ╔═╡ d31f91e8-6db6-4771-9544-8e54a816ecc1
begin
	
	positions = [lguys.extract_vector(out_special, :positions, i) for i in eachindex(orbit_labels)]
	velocities = [lguys.extract_vector(out_special, :velocities, i) for i in eachindex(orbit_labels)]
	Φs_ext = [lguys.extract(out_special, :potential_ext, i) for i in eachindex(orbit_labels)]

end

# ╔═╡ 5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
begin
	rs = lguys.radii.(positions)
	vs = lguys.radii.(velocities)
end

# ╔═╡ e5d40e2f-ac47-4827-853d-2f94bc39a624
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="radius / kpc",
	)

	for i in eachindex(orbit_labels)
		lines!(-out_special.times * lguys.T2GYR, rs[i], label=orbit_labels[i])
	
		hlines!([df_peris_apos_special.pericentre[i], df_peris_apos_special.apocentre[i]], linestyle=:dot)

		# scatter!(df_peris_apos_special.t_last_peri[i] * lguys.T2GYR, df_peris_apos_special.pericentre[i])
		
		# scatter!(df_peris_apos_special.t_last_apo[i] * lguys.T2GYR, df_peris_apos_special.apocentre[i])
	end

	Legend(fig[1, 2], ax)
	lguys.hide_grid!(ax)

	@savefig "r_time_orbits"
	fig
end

# ╔═╡ ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
positions

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
let

	fig = lguys.plot_xyz(positions..., labels=orbit_labels)

	resize_to_layout!(fig)

	@savefig "xyz_orbit"

	fig
end

# ╔═╡ 57a8d1c8-3940-4430-8b46-375fb2bf1695
let
	fig, ax = FigAxis(
		aspect=DataAspect(),
		xlabel="R / kpc",
		ylabel = "z / kpc",
	)

	for i in eachindex(positions)
		x = positions[i][1, :]
		y = positions[i][2, :]
		z = positions[i][3, :]
		R = @. sqrt(x^2 + y^2)

		lines!(R, z, label=orbit_labels[i])

	end
	@savefig  "R_z_orbit"

	axislegend()
	fig
end

# ╔═╡ 34efe422-e672-4d17-a558-ce32fb704a8e
lguys.plot_xyz(velocities..., units=" / km / s")

# ╔═╡ ad078920-225d-436e-835b-d87a9db53c49
let
	fig = Figure()
	ax = Axis(fig[1,1])

	for i in eachindex(orbit_labels)
		scatter!(rs[i], vs[i], color=out.times)
	end

	fig
end

# ╔═╡ 09bbae0d-ca3e-426d-b77c-69dd68ca42cc
begin
	M_b = 115
	r_b = 20
	
	a_exp(r) =  M_b / (r+r_b)^2
	phi_exp(r) = -M_b / (r+r_b)
end

# ╔═╡ 35f0ea14-a945-4745-910c-365b730676c5
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "time / Gyr",
		ylabel = "difference in potential from NFW",
	)
	
	for i in 1:length(orbit_labels)
		scatter!(out_special.times * lguys.T2GYR, Φs_ext[i] .- phi_exp.(rs[i]))
	end
	fig
end

# ╔═╡ 2e7c1798-4066-4c46-b5ed-732263728ac0
out[14]

# ╔═╡ f6b27164-ee7c-428b-aefb-75e89d178f3e
let
	fig = Figure()
	ax = Axis3(fig[1, 1], 
		xlabel = L"$x$ / kpc",
		ylabel = L"$y$ / kpc",
		zlabel = L"$z$ / kpc",
	)
	
	scatter!(snap.positions,
)
	fig
end

# ╔═╡ 8b818798-69fb-481d-ade1-9fd436b1f281
kms_label = L" / km\,s$^{-1}$"

# ╔═╡ 5fdd8307-d528-4cd7-a5e4-1f15aba75cd5
let
	fig = Figure()
	ax = Axis3(fig[1, 1], 
		xlabel = L"$v_x$ %$kms_label",
		ylabel = L"$v_y$ %$kms_label",
		zlabel = L"$v_z$ %$kms_label",
	)
	
	scatter!(-snap.velocities .* lguys.V2KMS,
)
	fig
end

# ╔═╡ 2dfe9a85-6553-4632-81e0-33c148fd1102
reverse(out.times)

# ╔═╡ b8c9823f-ca6b-48bf-9140-40440562dac0
import TOML

# ╔═╡ 3639400c-e5eb-4ce6-9ec5-7a4eafa9f458
obs_props = TOML.parsefile(obs_props_filename)

# ╔═╡ d51da98c-d3df-491d-a14c-bc6575e6c22a
err = get_coord_err(obs_props)

# ╔═╡ 106ebb4c-300f-4b7e-b9a2-85c0da493c8c
σv = obs_props["sigma_v"] .+ randn(length(peris)) * obs_props["sigma_v_err"]

# ╔═╡ 54e8d50a-2e01-465c-b423-afced2272a51
r_breaks = lguys.break_radius.(σv/lguys.V2KMS, delta_t)

# ╔═╡ a4fe3b8b-ca53-40fc-b729-e70e95bd9fd9
r_breaks_arcmin = lguys.kpc2arcmin.(r_breaks, [o.distance for o in observations])

# ╔═╡ 98b33c5a-b6cf-4472-a87b-3b8a278e0e36
mean(r_breaks_arcmin), std(r_breaks_arcmin)

# ╔═╡ cd0ea8cc-882b-479c-9f2a-3a433c43e0a0
hist(r_breaks)

# ╔═╡ 515fdab4-012a-4a16-a4e6-392fb203c7b9
mean(r_breaks), std(r_breaks)

# ╔═╡ 97ac4c60-ddb2-4165-a739-0e726d3bc690
obs_props["radial_velocity"]

# ╔═╡ 519a88f0-8e2d-4c09-83e0-3cc2ee147e35
function get_initial_t(rs)
	N = length(rs)
	t_ini = -1
	for t in N-2:-1:1
		if diff(rs)[t] >= 0 && diff(rs)[t+1] <= 0
			t_ini = t
			break
		end
	end
	return t_ini
end
	

# ╔═╡ 5f45e7c7-e447-48bf-ade4-38f516df2dad
for i in 1:length(orbit_labels)
	label = orbit_labels[i]
	fname = "$modelname_special/orbit_$label.csv"
	vel = -reverse(velocities[i], dims=2)
	pos = reverse(positions[i], dims=2)
	t0 = 1 + size(pos, 2) - get_initial_t(rs[i])
	println(t0)
	t = out_special.times[end] .- reverse(out_special.times)

	orbit_df = DataFrame(
		t = t[t0:end],
		x = pos[1, t0:end],
		y = pos[2, t0:end],
		z = pos[3, t0:end],
		v_x = vel[1, t0:end],
		v_y = vel[2, t0:end],
		v_z = vel[3, t0:end],
	)

	println("saving to $fname")
	CSV.write(fname, orbit_df)
end

# ╔═╡ 1152cd63-baab-426a-b464-b10857eed4ec
for i in 1:length(orbit_labels)
	label = orbit_labels[i]
	fname = "$modelname_special/orbit_$label.toml"
	o = observations_special[i]

	t = get_initial_t(rs[i])

	properties = OrderedDict(
		"name" => orbit_labels[i],
		"ra" => o.ra,
		"dec" => o.dec,
		"distance" => o.distance,
		"pmra" => o.pmra,
		"pmdec" => o.pmdec,
		"radial_velocity" => o.radial_velocity,
		"distance_err" => err["distance"],
		"pmra_err" => err["pmra"],
		"pmdec_err" => err["pmdec"],
		"radial_velocity_err" => err["radial_velocity"],
		"pericentre" => df_peris_apos_special.pericentre[i],
		"apocentre" => df_peris_apos_special.apocentre[i],
		"t_apo_i" => lguys.T2GYR*(out_special.times[end] - out.times[t]),
		"r_i" => rs[i][t],
		"x_i" => positions[i][1, t],
		"y_i" => positions[i][2, t],
		"z_i" => positions[i][3, t],
		"v_x_i" => -1*lguys.V2KMS* velocities[i][1, t],
		"v_y_i" => -1*lguys.V2KMS* velocities[i][2, t],
		"v_z_i" => -1*lguys.V2KMS* velocities[i][3, t],
	)

	println("saving to $fname")
	println(properties)
	open(fname, "w") do f
		TOML.print(f, properties)
	end
end

# ╔═╡ 5316884b-3971-4ca7-9106-f638241d3388
get_initial_t(rs[1])

# ╔═╡ 9bb4eba3-e6f1-43a7-a202-7dc9375661b9
rs[1][1853:1856]

# ╔═╡ a88df07c-d6d0-49be-8bc9-d9adc5647838
rs[1][1854]

# ╔═╡ a23c4b5c-8752-4188-b917-ebfbd82e4422
rs[1][end:-1:1][length(rs[1]) + 1 - 1854]

# ╔═╡ de1e5245-0946-47cd-8e2c-ba1914cfeb74
begin 
	# orbit info
	for i in 1:length(orbit_labels)
		t = get_initial_t(rs[i])
		@printf "orbit: \t\t %i\n" i
		
		@printf "pericentre:\t %0.1f\n" df_peris_apos_special.pericentre[i]
		@printf "apocentre: \t %0.1f\n" df_peris_apos_special.apocentre[i]

		@printf "time of first apocentre: %f \n" out_special.times[end] - out.times[t]
		@printf "radius of first apocentre: %f\n" rs[i][t]
		@printf "intial position: [%f, %f, %f]\n" positions[i][:, t]...
		@printf "intial velocity: [%f, %f, %f]\n" -1* velocities[i][:, t]...
		@printf "final position: [%f, %f, %f]\n" positions[i][:, 1]...
		@printf "final velocity: [%f, %f, %f]\n" -lguys.V2KMS * velocities[i][:, 1]...

		o = observations_special[i]
		@printf "ra: \t %f\n" o.ra
		@printf "dec: \t %f\n" o.dec
		@printf "pmra: \t %f\n" o.pmra
		@printf "pmdec: \t %f\n" o.pmdec
		@printf "dist: \t %f\n" o.distance
		@printf "rv: \t %f\n" o.radial_velocity

		println()
	end
end

# ╔═╡ Cell order:
# ╟─7450144e-5464-4036-a215-b6e2cd270405
# ╠═baf5a0e0-3ee1-4056-b3a5-15058b4a0e53
# ╠═abac7f99-fc5d-4536-a845-52cc269deee7
# ╠═2559bc4c-3bad-4d70-b02f-8fe48b772fb6
# ╠═e5f728b8-8412-4f57-ad38-a0a35bb08a48
# ╠═3be492fa-0cb9-40f1-a39e-e25bf485fd3e
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═91cac7f1-e3a7-4f4e-a136-d99cbfef75ed
# ╠═b3e118f4-0503-4a6c-8f0d-d096bd0236a3
# ╠═d975d00c-fd69-4dd0-90d4-c4cbe73d9754
# ╠═3b83205d-91c1-481e-9305-0d59bc692135
# ╟─88536e86-cf2a-4dff-ae64-514821957d40
# ╟─26d616da-95ec-4fb9-b9a8-2f095d74c722
# ╟─38fd9158-7807-4a26-b534-b0bcb5dbb0e5
# ╠═9a22d47b-8474-4596-b418-de33eb07c627
# ╠═dcdca4a7-08a6-4d1b-a972-c386493207d0
# ╠═56631a4d-323e-42ac-ad71-3e35face895d
# ╠═3639400c-e5eb-4ce6-9ec5-7a4eafa9f458
# ╠═83fca0a6-5693-4b06-89ff-2d097231c402
# ╠═d51da98c-d3df-491d-a14c-bc6575e6c22a
# ╠═106ebb4c-300f-4b7e-b9a2-85c0da493c8c
# ╠═c60b4d12-642f-4247-a022-f2c0bb16a092
# ╠═54e8d50a-2e01-465c-b423-afced2272a51
# ╠═a4fe3b8b-ca53-40fc-b729-e70e95bd9fd9
# ╠═cd0ea8cc-882b-479c-9f2a-3a433c43e0a0
# ╠═515fdab4-012a-4a16-a4e6-392fb203c7b9
# ╠═98b33c5a-b6cf-4472-a87b-3b8a278e0e36
# ╠═bbf3f229-fc3c-46ae-af28-0f8bd81e7d32
# ╠═4481a5e1-d635-4fea-a5c5-c85f0d6df62f
# ╠═384be6a6-f9d9-47e0-9792-aef6689dcbdb
# ╟─1acef60e-60d6-47ba-85fd-f9780934788b
# ╠═ca1c236e-795a-408b-845b-9c13bc838619
# ╠═46b4242b-8af7-4233-8ecf-d86740b4c884
# ╠═377628e5-67a7-45fa-9e98-992acf7171ea
# ╠═92aac8e8-d599-4a1e-965e-e653bc54c509
# ╠═2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
# ╠═8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
# ╠═5f11f6ab-c9ab-4600-acca-a0bb84d81a12
# ╠═3e18be79-9403-4104-9d4f-9860163ad8b9
# ╠═ede3836c-740d-4ac7-bbc7-3165981a1878
# ╠═44660b2f-6220-473b-bb2f-07e23b176491
# ╠═67ca9626-8416-4254-87f8-d1c5c88bf2de
# ╠═8a863d67-2535-4117-ae19-437b3e4c6c71
# ╠═97ac4c60-ddb2-4165-a739-0e726d3bc690
# ╠═c48b4e73-480e-4a50-b5fc-db5f6c5b040e
# ╟─43d43f63-4c13-4b23-950e-ada59aa86bc9
# ╠═f14fb67d-37c1-46f8-9d63-3f0d4fd93b3c
# ╠═ac81acd8-4a78-4230-bc70-3b78a861b618
# ╟─16f4ac20-d8cf-4218-8c01-c15e04e567fb
# ╠═e5d40e2f-ac47-4827-853d-2f94bc39a624
# ╠═d31f91e8-6db6-4771-9544-8e54a816ecc1
# ╠═5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
# ╠═ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
# ╠═57a8d1c8-3940-4430-8b46-375fb2bf1695
# ╠═34efe422-e672-4d17-a558-ce32fb704a8e
# ╠═ad078920-225d-436e-835b-d87a9db53c49
# ╠═09bbae0d-ca3e-426d-b77c-69dd68ca42cc
# ╠═35f0ea14-a945-4745-910c-365b730676c5
# ╟─2e7c1798-4066-4c46-b5ed-732263728ac0
# ╟─f6b27164-ee7c-428b-aefb-75e89d178f3e
# ╠═8b818798-69fb-481d-ade1-9fd436b1f281
# ╟─5fdd8307-d528-4cd7-a5e4-1f15aba75cd5
# ╠═2dfe9a85-6553-4632-81e0-33c148fd1102
# ╠═5f45e7c7-e447-48bf-ade4-38f516df2dad
# ╠═b8c9823f-ca6b-48bf-9140-40440562dac0
# ╠═1be546d7-8b11-4dc5-9295-393d39f65123
# ╠═1152cd63-baab-426a-b464-b10857eed4ec
# ╠═519a88f0-8e2d-4c09-83e0-3cc2ee147e35
# ╠═5316884b-3971-4ca7-9106-f638241d3388
# ╠═9bb4eba3-e6f1-43a7-a202-7dc9375661b9
# ╠═a88df07c-d6d0-49be-8bc9-d9adc5647838
# ╠═a23c4b5c-8752-4188-b917-ebfbd82e4422
# ╠═de1e5245-0946-47cd-8e2c-ba1914cfeb74
