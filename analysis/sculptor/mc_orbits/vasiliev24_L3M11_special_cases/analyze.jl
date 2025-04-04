### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ e9e2c787-4e0e-4169-a4a3-401fea21baba
begin 
	import Pkg
	Pkg.activate()
	using CairoMakie
	using DataFrames, CSV

	using Printf
	import LilGuys as lguys

	using Arya
end


# ╔═╡ acce474f-5b88-42dc-aa1b-840f3d3967b7
using LilGuys

# ╔═╡ 1be546d7-8b11-4dc5-9295-393d39f65123
using OrderedCollections

# ╔═╡ d975d00c-fd69-4dd0-90d4-c4cbe73d9754
using Statistics, Distributions

# ╔═╡ 7450144e-5464-4036-a215-b6e2cd270405
md"""
This notebook analyzes the result of the MC samples of orbits in the same potential to determine the plausable range of pericentres and apocentres
"""

# ╔═╡ 3be492fa-0cb9-40f1-a39e-e25bf485fd3e
md"""
# Setup
"""

# ╔═╡ 0f1eadfc-bff3-42b8-beee-2765589671ac
save = Makie.save

# ╔═╡ b8c9823f-ca6b-48bf-9140-40440562dac0
import TOML

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

# ╔═╡ fdbe2478-32eb-4264-b484-a7b072a4028d
sim_dir = ENV["DWARFS_ROOT"] * "/simulations/sculptor/mc_orbits/vasiliev24_L3M11_special_cases/"

# ╔═╡ cd64bd83-fee6-4fed-94ab-131ed49ebef0
if !isdefined(Main, :obs)
	include(joinpath(sim_dir, "sample.jl"))
end

# ╔═╡ 380b8e74-2837-4620-b600-ae6674b35d16
p_value = 0.001349898031630093 # 3sigma

# ╔═╡ 46348ecb-ee07-4b6a-af03-fc4f2635f57b
fig_dir = "./figures"

# ╔═╡ 9a22d47b-8474-4596-b418-de33eb07c627
begin 
	out = lguys.Output("../vasiliev24_L3M11/");

	df_peris_apos = lguys.read_fits("../vasiliev24_L3M11/peris_apos.fits")
	snap = out[1] |> sort_snap

	@assert all(snap.index  .== df_peris_apos.index) "snapshot and peri apo index must match"
end

# ╔═╡ 2559bc4c-3bad-4d70-b02f-8fe48b772fb6
special_orbits_dir = "."

# ╔═╡ e5f728b8-8412-4f57-ad38-a0a35bb08a48
orbit_labels = ["mean", "smallperi", "largeperi"]

# ╔═╡ dcdca4a7-08a6-4d1b-a972-c386493207d0
begin 
	out_special = lguys.Output(special_orbits_dir);

	df_peris_apos_special = lguys.read_fits("$special_orbits_dir/peris_apos.fits")
	snap_special = out_special[1] |> sort_snap

	@assert all(snap_special.index  .== df_peris_apos_special.index) "snapshot and peri apo index must match"
end

# ╔═╡ bbf3f229-fc3c-46ae-af28-0f8bd81e7d32
df_peris_apos_special

# ╔═╡ 4481a5e1-d635-4fea-a5c5-c85f0d6df62f
peris = df_peris_apos.pericentre

# ╔═╡ 413d4e5d-c9cd-4aca-be1e-d132b2bd616d
peri_qs = lguys.quantile(peris, [p_value, 1-p_value, 0.5])

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
	vlines!(df_peris_apos_special.pericentre)

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
	vlines!(df_peris_apos_special.apocentre)

	fig
end

# ╔═╡ 73f827a2-a529-49ce-b9be-0b256d4ac768
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "dt lmc / Gyr",
		ylabel = "count"
	)

	bins, counts, err = lguys.histogram(df_peris_apos.t_last_peri_lmc * T2GYR)
	scatter!(lguys.midpoints(bins), counts)
	vlines!(df_peris_apos_special.t_last_peri_lmc * T2GYR)

	fig
end

# ╔═╡ a57d212d-e5b4-4fc4-aace-de3b1ab67c57
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "perilmc",
		ylabel = "count"
	)

	bins, counts, err = lguys.histogram(df_peris_apos.peri_lmc)
	scatter!(lguys.midpoints(bins), counts)
	vlines!(df_peris_apos_special.peri_lmc )

	fig
end

# ╔═╡ 2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
ϵ = -lguys.calc_K_spec(out[1]) .- out[1].Φs_ext

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
		scatter!(x, peris; plot_kwargs...)
		
		for i in eachindex(orbit_labels)
				x = getproperty(observations_special[i], sym)
				y = df_peris_apos_special[i, :pericentre]
				scatter!(x, y; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)


	save(joinpath(fig_dir, "peri_mc_orbits_corr.pdf"), fig)
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


	save(joinpath(fig_dir, "apo_mc_orbits_corr.pdf"), fig)
	fig
end

# ╔═╡ 16ce3cab-d92f-41ef-8d6a-a449fe6a8e58
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
		:ylabel => "perilmc / kpc",
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
		scatter!(x, df_peris_apos.peri_lmc; plot_kwargs...)
		
		for i in eachindex(orbit_labels)
				x = getproperty(observations_special[i], sym)
				y = df_peris_apos_special[i, :peri_lmc]
				scatter!(x, y; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)


	save(joinpath(fig_dir, "peri_lmc_mc_corr.pdf"), fig)
	fig
end

# ╔═╡ 4dbbf202-bc06-4f23-bddc-8aeb1ae1bd7a
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
		:ylabel => "t last peri lmc / Gyr",
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
		scatter!(x, df_peris_apos.t_last_peri * T2GYR; plot_kwargs...)
		
		for i in eachindex(orbit_labels)
				x = getproperty(observations_special[i], sym)
				y = df_peris_apos_special[i, :t_last_peri] * T2GYR
				scatter!(x, y; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)


	save(joinpath(fig_dir, "t_peri_mc_corr.pdf"), fig)
	fig
end

# ╔═╡ d1801170-9efa-4234-8666-fa7b38a35c9b
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
		:ylabel => "t last peri lmc / Gyr",
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
		scatter!(x, df_peris_apos.t_last_peri_lmc * T2GYR; plot_kwargs...)
		
		for i in eachindex(orbit_labels)
				x = getproperty(observations_special[i], sym)
				y = df_peris_apos_special[i, :t_last_peri_lmc] * T2GYR
				scatter!(x, y; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)


	save(joinpath(fig_dir, "t_perilmc_mc_corr.pdf"), fig)
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

# ╔═╡ 4b29ef51-5047-4c8a-9692-a388a6e51c2b
traj_lmc = CSV.read("lmc_traj.csv", DataFrame)

# ╔═╡ 40cd8af5-70cd-45eb-9841-a9fce83e9f9d
traj_lmc.time == out.times

# ╔═╡ 921b889a-376a-43dc-bf86-302e0e2521d0
pos_lmc = [traj_lmc.x traj_lmc.y traj_lmc.z]'

# ╔═╡ 6cfc3e05-17e0-4054-863f-9aedf215af4f


# ╔═╡ 43b8d382-6226-46bb-9cbe-6463eb0415d0
vel_lmc = [traj_lmc.v_x traj_lmc.v_y traj_lmc.v_z]'

# ╔═╡ d31f91e8-6db6-4771-9544-8e54a816ecc1
begin
	
	positions = [lguys.extract_vector(out_special, :positions, i) for i in eachindex(orbit_labels)]
	velocities = [lguys.extract_vector(out_special, :velocities, i) for i in eachindex(orbit_labels)]
	accelerations = [lguys.extract_vector(out_special, :accelerations, i) for i in eachindex(orbit_labels)]
	Φs_ext = [lguys.extract(out_special, :Φs_ext, i) for i in eachindex(orbit_labels)]

end

# ╔═╡ 3831432c-488f-4cfd-90d9-d137930c6bb2
positions_scl_lmc = [pos .- pos_lmc for pos in positions]

# ╔═╡ a7053872-c788-4d9f-a6a5-e53f83e95aff
velocities_scl_lmc = [vel .- vel_lmc for vel in velocities]

# ╔═╡ 5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
begin
	rs = lguys.calc_r.(positions)
	vs = lguys.calc_r.(velocities)
end

# ╔═╡ e5d40e2f-ac47-4827-853d-2f94bc39a624
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="radius / kpc",
	)

	for i in eachindex(orbit_labels)
		lines!(out_special.times * lguys.T2GYR, rs[i], label=orbit_labels[i])
	
		hlines!([df_peris_apos_special.pericentre[i], df_peris_apos_special.apocentre[i]], linestyle=:dot)
	end

	Legend(fig[1, 2], ax)
	lguys.MakieExtension.hide_grid!(ax)

	save(joinpath(fig_dir, "r_time_orbits.pdf"), fig)
	fig
end

# ╔═╡ 27754eca-1c6f-4cc4-8153-a4bacf4d473b
accs = lguys.calc_r.(accelerations)

# ╔═╡ cf4e92a9-80a3-4ef8-9404-0ee6707a3401
rs_scl_lmc = lguys.calc_r.(positions_scl_lmc)

# ╔═╡ 7afbee36-7c8b-4d73-b5c0-079d2346f020
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="radius / kpc",
	)

	for i in eachindex(orbit_labels)
		lines!(out_special.times * lguys.T2GYR, rs_scl_lmc[i], label=orbit_labels[i])
	
		hlines!([df_peris_apos_special.peri_lmc[i], df_peris_apos_special.apo_lmc[i]], linestyle=:dot)
	end

	Legend(fig[1, 2], ax)
	lguys.MakieExtension.hide_grid!(ax)

	fig
end

# ╔═╡ ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
positions

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
let

	fig = lguys.MakieExtension.plot_xyz(positions..., labels=orbit_labels)

	resize_to_layout!(fig)

	save(joinpath(fig_dir, "xyz_orbit.pdf"), fig)

	fig
end

# ╔═╡ e2e31f30-bbf3-41b5-ba87-5b3a6db152bb
let

	fig = lguys.MakieExtension.plot_xyz(positions_scl_lmc..., labels=orbit_labels)

	resize_to_layout!(fig)

	fig
end

# ╔═╡ c6746c48-4547-483e-807d-07865f0576c6
 lguys.MakieExtension.plot_xyz(velocities_scl_lmc..., labels=orbit_labels)

# ╔═╡ 83aa4296-047a-455d-a5ad-d9f9d111165b
 lguys.MakieExtension.plot_xyz(velocities..., labels=orbit_labels)

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
	save(joinpath(fig_dir, "R_z_orbit.pdf"), fig)

	axislegend()
	fig
end

# ╔═╡ ad078920-225d-436e-835b-d87a9db53c49
let
	fig = Figure()
	ax = Axis(fig[1,1])

	for i in eachindex(orbit_labels)
		scatter!(rs[i], vs[i], color=out.times)
	end

	fig
end

# ╔═╡ 12eb9758-74dc-44ba-be9f-7bb9eff7ab43
let
	fig = Figure()
	ax = Axis(fig[1,1])

	for i in eachindex(orbit_labels)
		scatter!(rs[i], accs[i], color=out.times)
	end

	fig
end

# ╔═╡ 2665bd17-fa3d-4e00-bbdf-b5b7cb0f67c7
let
	fig = Figure()
	ax = Axis(fig[1,1])

	for i in eachindex(orbit_labels)
		scatter!(rs_scl_lmc[i], accs[i], color=out.times)
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

# ╔═╡ 5f45e7c7-e447-48bf-ade4-38f516df2dad
for i in 1:length(orbit_labels)
	label = orbit_labels[i]
	fname = "orbit_$label.csv"
	vel = -reverse(velocities[i], dims=2)
	pos = reverse(positions[i], dims=2)
	t0 = 1
	
	t = -reverse(out_special.times)

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

# ╔═╡ aeb3c2ca-5caf-4b75-a6a9-4699e29240f6


# ╔═╡ 1152cd63-baab-426a-b464-b10857eed4ec
for i in 1:length(orbit_labels)
	label = orbit_labels[i]
	fname = "orbit_$label.toml"
	o = observations_special[i]

	t = 1

	properties = OrderedDict(
		"name" => orbit_labels[i],
		"ra" => o.ra,
		"dec" => o.dec,
		"distance" => o.distance,
		"pmra" => o.pmra,
		"pmdec" => o.pmdec,
		"radial_velocity" => o.radial_velocity,
		"distance_err" => err.distance,
		"pmra_err" => err.pmra,
		"pmdec_err" => err.pmdec,
		"radial_velocity_err" => err.radial_velocity,
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

# ╔═╡ de1e5245-0946-47cd-8e2c-ba1914cfeb74
begin 
	# orbit info
	for i in 1:length(orbit_labels)
		t = length(rs[i])
		@printf "orbit: \t\t %i\n" i
		
		@printf "pericentre:\t %0.1f\n" df_peris_apos_special.pericentre[i]
		@printf "apocentre: \t %0.1f\n" df_peris_apos_special.apocentre[i]
		@printf "t last peri:\t %0.2f\n" df_peris_apos_special.t_last_peri[i] * T2GYR
		
		@printf "perilmc:\t %0.1f\n" df_peris_apos_special.peri_lmc[i]
		@printf "apolmc: \t %0.1f\n" df_peris_apos_special.apo_lmc[i]
		@printf "t last perilmc:\t %0.2f\n" df_peris_apos_special.t_last_peri_lmc[i] * T2GYR

	
		@printf "intial position: [%0.4f, %0.4f, %0.4f]\n" positions[i][:, t]...
		@printf "intial velocity: [%0.2f, %0.2f, %0.2f]\n" -lguys.V2KMS* velocities[i][:, t]...
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

# ╔═╡ 7796da0f-1efb-4757-8610-53a9f9e9854f
out.times * lguys.T2GYR

# ╔═╡ Cell order:
# ╟─7450144e-5464-4036-a215-b6e2cd270405
# ╠═3be492fa-0cb9-40f1-a39e-e25bf485fd3e
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═0f1eadfc-bff3-42b8-beee-2765589671ac
# ╠═acce474f-5b88-42dc-aa1b-840f3d3967b7
# ╠═b8c9823f-ca6b-48bf-9140-40440562dac0
# ╠═1be546d7-8b11-4dc5-9295-393d39f65123
# ╠═d975d00c-fd69-4dd0-90d4-c4cbe73d9754
# ╠═cd64bd83-fee6-4fed-94ab-131ed49ebef0
# ╠═3b83205d-91c1-481e-9305-0d59bc692135
# ╟─88536e86-cf2a-4dff-ae64-514821957d40
# ╠═26d616da-95ec-4fb9-b9a8-2f095d74c722
# ╠═38fd9158-7807-4a26-b534-b0bcb5dbb0e5
# ╠═fdbe2478-32eb-4264-b484-a7b072a4028d
# ╠═380b8e74-2837-4620-b600-ae6674b35d16
# ╠═46348ecb-ee07-4b6a-af03-fc4f2635f57b
# ╠═9a22d47b-8474-4596-b418-de33eb07c627
# ╠═2559bc4c-3bad-4d70-b02f-8fe48b772fb6
# ╠═e5f728b8-8412-4f57-ad38-a0a35bb08a48
# ╠═dcdca4a7-08a6-4d1b-a972-c386493207d0
# ╠═413d4e5d-c9cd-4aca-be1e-d132b2bd616d
# ╠═bbf3f229-fc3c-46ae-af28-0f8bd81e7d32
# ╠═4481a5e1-d635-4fea-a5c5-c85f0d6df62f
# ╠═384be6a6-f9d9-47e0-9792-aef6689dcbdb
# ╟─1acef60e-60d6-47ba-85fd-f9780934788b
# ╠═ca1c236e-795a-408b-845b-9c13bc838619
# ╠═46b4242b-8af7-4233-8ecf-d86740b4c884
# ╠═73f827a2-a529-49ce-b9be-0b256d4ac768
# ╠═a57d212d-e5b4-4fc4-aace-de3b1ab67c57
# ╠═92aac8e8-d599-4a1e-965e-e653bc54c509
# ╠═2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
# ╠═8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
# ╠═5f11f6ab-c9ab-4600-acca-a0bb84d81a12
# ╠═3e18be79-9403-4104-9d4f-9860163ad8b9
# ╠═ede3836c-740d-4ac7-bbc7-3165981a1878
# ╠═44660b2f-6220-473b-bb2f-07e23b176491
# ╠═67ca9626-8416-4254-87f8-d1c5c88bf2de
# ╠═c48b4e73-480e-4a50-b5fc-db5f6c5b040e
# ╟─43d43f63-4c13-4b23-950e-ada59aa86bc9
# ╠═16ce3cab-d92f-41ef-8d6a-a449fe6a8e58
# ╠═4dbbf202-bc06-4f23-bddc-8aeb1ae1bd7a
# ╠═d1801170-9efa-4234-8666-fa7b38a35c9b
# ╠═ac81acd8-4a78-4230-bc70-3b78a861b618
# ╟─16f4ac20-d8cf-4218-8c01-c15e04e567fb
# ╠═4b29ef51-5047-4c8a-9692-a388a6e51c2b
# ╠═40cd8af5-70cd-45eb-9841-a9fce83e9f9d
# ╠═921b889a-376a-43dc-bf86-302e0e2521d0
# ╠═6cfc3e05-17e0-4054-863f-9aedf215af4f
# ╠═43b8d382-6226-46bb-9cbe-6463eb0415d0
# ╠═e5d40e2f-ac47-4827-853d-2f94bc39a624
# ╠═7afbee36-7c8b-4d73-b5c0-079d2346f020
# ╠═d31f91e8-6db6-4771-9544-8e54a816ecc1
# ╠═3831432c-488f-4cfd-90d9-d137930c6bb2
# ╠═a7053872-c788-4d9f-a6a5-e53f83e95aff
# ╠═5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
# ╠═27754eca-1c6f-4cc4-8153-a4bacf4d473b
# ╠═cf4e92a9-80a3-4ef8-9404-0ee6707a3401
# ╠═ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
# ╠═e2e31f30-bbf3-41b5-ba87-5b3a6db152bb
# ╠═c6746c48-4547-483e-807d-07865f0576c6
# ╠═83aa4296-047a-455d-a5ad-d9f9d111165b
# ╠═57a8d1c8-3940-4430-8b46-375fb2bf1695
# ╠═ad078920-225d-436e-835b-d87a9db53c49
# ╠═12eb9758-74dc-44ba-be9f-7bb9eff7ab43
# ╠═2665bd17-fa3d-4e00-bbdf-b5b7cb0f67c7
# ╠═09bbae0d-ca3e-426d-b77c-69dd68ca42cc
# ╠═35f0ea14-a945-4745-910c-365b730676c5
# ╟─2e7c1798-4066-4c46-b5ed-732263728ac0
# ╟─f6b27164-ee7c-428b-aefb-75e89d178f3e
# ╠═8b818798-69fb-481d-ade1-9fd436b1f281
# ╟─5fdd8307-d528-4cd7-a5e4-1f15aba75cd5
# ╠═2dfe9a85-6553-4632-81e0-33c148fd1102
# ╠═5f45e7c7-e447-48bf-ade4-38f516df2dad
# ╠═aeb3c2ca-5caf-4b75-a6a9-4699e29240f6
# ╠═1152cd63-baab-426a-b464-b10857eed4ec
# ╠═de1e5245-0946-47cd-8e2c-ba1914cfeb74
# ╠═7796da0f-1efb-4757-8610-53a9f9e9854f
