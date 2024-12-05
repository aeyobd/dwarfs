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


# ╔═╡ d975d00c-fd69-4dd0-90d4-c4cbe73d9754
using Statistics, Distributions

# ╔═╡ 7450144e-5464-4036-a215-b6e2cd270405
md"""
This notebook analyzes the result of the MC samples of orbits in the same potential to determine the plausable range of pericentres and apocentres
"""

# ╔═╡ a7ce5b0c-84a6-4d63-94f1-68e7a0d9e758
obs_prop_filename = ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml"

# ╔═╡ 1e50dbfe-09e5-4f42-83c3-a8291b8e1b1a
obs = lguys.coord_from_file(obs_prop_filename)

# ╔═╡ 0146ee17-de5f-4877-aaa6-83a898e01416
err = lguys.coord_err_from_file(obs_prop_filename)

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

# ╔═╡ 9a22d47b-8474-4596-b418-de33eb07c627
begin 
	out = lguys.Output(".");

	df_peris_apos = lguys.read_fits("peris_apos.fits")
	snap = out[1] |> sort_snap

	@assert all(snap.index  .== df_peris_apos.index) "snapshot and peri apo index must match"
end

# ╔═╡ cd2c1a82-8b0b-4abc-8831-f621e5ec483a
df_peris_apos

# ╔═╡ b9f469ed-6e4e-41ee-ac75-1b5bfa0a114a
p_value = 0.0014

# ╔═╡ fb6debf2-0161-477f-b29b-5a0f1f70f340
[snap.positions[:, 1]; snap.velocities[:, 1] * lguys.V2KMS]

# ╔═╡ 4481a5e1-d635-4fea-a5c5-c85f0d6df62f
peris = df_peris_apos.pericentre

# ╔═╡ e3084451-1492-49c5-b5dc-8a085e513640
length(peris)

# ╔═╡ 384be6a6-f9d9-47e0-9792-aef6689dcbdb
apos = df_peris_apos.apocentre

# ╔═╡ a1b2edd5-20e4-456f-a8e5-4634b0a433c2
perilmc = df_peris_apos.peri_lmc

# ╔═╡ 413d4e5d-c9cd-4aca-be1e-d132b2bd616d
peri_qs = lguys.quantile(perilmc, [p_value, 1-p_value])

# ╔═╡ 17a63cc8-84f4-4248-a7b0-c8378454b1f7
idx = [1, 
	argmin(abs.(peri_qs[1] .- perilmc)),
	argmin(abs.(peri_qs[2] .- perilmc)),
]

# ╔═╡ bbf3f229-fc3c-46ae-af28-0f8bd81e7d32
peris[idx]

# ╔═╡ b0209748-88d7-4b60-b6b9-7a367a5d1e58
quantile(perilmc, [p_value, 0.5, 1-p_value])

# ╔═╡ e5f728b8-8412-4f57-ad38-a0a35bb08a48
orbit_labels = ["mean", "smallperi", "largeperi"]

# ╔═╡ 46348ecb-ee07-4b6a-af03-fc4f2635f57b
fig_dir = "./figures"

# ╔═╡ 6651d141-f6ca-4e8f-a785-5b14275b367c
T2GYR = lguys.T2GYR

# ╔═╡ b2fffcba-c474-430c-8100-037a0d352d69
Δt_lmc = df_peris_apos.t_last_peri_lmc * T2GYR

# ╔═╡ 850c29d0-f395-448f-87be-d788210594c2
Δt_mw = df_peris_apos.t_last_peri * T2GYR

# ╔═╡ 73af5376-3b38-4d5f-b6e2-1cafba27dabd
begin
	# loads in trajectory of lmc in Vasiliev 2021
	lmc_file = ENV["DWARFS_ROOT"] * "/agama/potentials/vasiliev+21/trajlmc.txt"
	lmc_traj = CSV.read(lmc_file, DataFrame, delim=" ", header = [:time, :x, :y, :z, :v_x, :v_y, :v_z])
	
	lmc_x = lguys.lerp(lmc_traj.time, lmc_traj.x)
	lmc_y = lguys.lerp(lmc_traj.time, lmc_traj.y)
	lmc_z = lguys.lerp(lmc_traj.time, lmc_traj.z)
	lmc_v_x = lguys.lerp(lmc_traj.time, lmc_traj.v_x)
	lmc_v_y = lguys.lerp(lmc_traj.time, lmc_traj.v_y)
	lmc_v_z = lguys.lerp(lmc_traj.time, lmc_traj.v_z)

	V_T2GYR = 1 / lguys.kpc_per_Gyr_per_kms
	
	times_v = -out.times * T2GYR / V_T2GYR
	
	lmc_pos = reshape([lmc_x.(times_v) lmc_y.(times_v) lmc_z.(times_v)]', (3, :))
	lmc_vel = reshape([lmc_v_x.(times_v) lmc_v_y.(times_v) lmc_v_z.(times_v)]', (3, :))

end

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

	fig
end

# ╔═╡ 7b7b8357-ebf1-4257-8a38-272b245842cc
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "perilmc / kpc",
		ylabel = "count"
	)

	bins, counts, err = lguys.histogram(perilmc)
	scatter!(lguys.midpoints(bins), counts)

	fig
end

# ╔═╡ 6e658f28-0100-41b8-af8a-bb832baa1d12
unique(Δt_lmc)

# ╔═╡ c43412d1-3a9b-4cd0-901b-05048c0b07ec
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "delta t MW / Gyr",
		ylabel = "count"
	)

	hist!(Δt_mw)

	fig
end

# ╔═╡ 88cbd5b4-29a9-4c2f-ad9e-657460d6314c
0.1 / T2GYR, 0.12/ T2GYR

# ╔═╡ 7d7bc811-ced1-44a1-8378-d6c20b49d0e9


# ╔═╡ 7901c3ae-91b4-4fc1-8dc3-9b617c268533
open("outputlist.txt", "w") do f
	t = 0
	t_max = 1036.892895
	i = 0
	while t < t_max
		println(f, t)
		if 0.1  < t*T2GYR< 0.12
			t += 0.1
		else
			t += 1
		end

		i += 1
	end

	println(f, t_max)

	println("i = ", i)
end


# ╔═╡ 4cbbebec-f231-4af4-a10f-e01d090ea2ca
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "delta t LMC / Gyr",
		ylabel = "count"
	)

	hist!(Δt_lmc)

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

	fig
end

# ╔═╡ 471a5501-bc7b-4114-be1a-9088b4e674cd
hist(lguys.calc_r(snap.positions),
	axis=(; xlabel="initial galactocentric radius / kpc")
)

# ╔═╡ 68b0383a-3d5a-4b94-967c-f0e31e8a0ce1
hist(lguys.calc_r(snap.velocities) * lguys.V2KMS,
	axis=(; xlabel="initial galactocentric velocity / km/s")
)

# ╔═╡ 2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
ϵ = -lguys.calc_K_spec(out[1]) .- out[1].Φs_ext

# ╔═╡ 8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
extrema(ϵ)

# ╔═╡ 5f11f6ab-c9ab-4600-acca-a0bb84d81a12
begin
	points = [lguys.Galactocentric(
		snap.positions[:, i]*lguys.R2KPC, 
		-snap.velocities[:, i]*lguys.V2KMS)
		for i in 1:length(snap)]
	
	observations = lguys.transform.(lguys.ICRS, points)
end


# ╔═╡ 92aac8e8-d599-4a1e-965e-e653bc54c509
dists = getproperty.(observations, :distance)

# ╔═╡ ede3836c-740d-4ac7-bbc7-3165981a1878
normal_dist(x, μ, σ) = 1/√(2π) * 1/σ * exp(-(x-μ)^2/2σ^2)

# ╔═╡ 44660b2f-6220-473b-bb2f-07e23b176491
columns = [:ra, :dec, :distance, :pmra, :pmdec, :radial_velocity]

# ╔═╡ d3063e30-2cb3-4f1b-8546-0d5e81d90d9f
let
	
	for sym in columns[1:2]
	    x = [getproperty(o, sym) for o in observations]
	    y = peris

		
	    p = scatter(x, y, alpha=0.1,
			axis=(; xlabel=coord_labels[sym], ylabel="pericentre / kpc")
		)
	    @info p 
	end

end

# ╔═╡ 643676b7-ae9b-4d5b-b99d-f8bd338cf8a6
mkdir("figures")

# ╔═╡ c48b4e73-480e-4a50-b5fc-db5f6c5b040e
let
	fig = Figure(size=(600, 600))
	plot_kwargs = Dict(
		:color => :black,
		:alpha => 0.1,
		:markersize => 1,
		:rasterize => true,
	)

	
	orbit_points_kwargs = [Dict(
		:color => COLORS[i],
		:alpha => 1,
		:markersize => 10,
		:label => orbit_labels[i]
	) for i in eachindex(idx)
		]
	
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
		for i in eachindex(idx)
			scatter!(x[idx[i]], peris[idx[i]]; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)


	save(joinpath(fig_dir, "peri_mc_orbits_corr.pdf"), fig)
	fig
end

# ╔═╡ 99cbc08e-5b0e-4748-9727-667e0eca2f74
let
	fig = Figure(size=(600, 600))
	plot_kwargs = Dict(
		:color => :black,
		:alpha => 0.1,
		:markersize => 1,
		:rasterize => true,
	)

	
	orbit_points_kwargs = [Dict(
		:color => COLORS[i],
		:alpha => 1,
		:markersize => 10,
		:label => orbit_labels[i]
	) for i in eachindex(idx)
		]
	
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
		scatter!(x, perilmc; plot_kwargs...)
		for i in eachindex(idx)
			scatter!(x[idx[i]], perilmc[idx[i]]; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)

	fig
end

# ╔═╡ 456acc60-e7e1-428e-bd5d-bc72a18a146f
let
	fig = Figure(size=(600, 600))
	plot_kwargs = Dict(
		:color => :black,
		:alpha => 0.1,
		:markersize => 1,
	)

	
	orbit_points_kwargs = [Dict(
		:color => COLORS[i],
		:alpha => 1,
		:markersize => 10,
		:label => orbit_labels[i]
	) for i in eachindex(idx)
		]
	
	ax_kwargs = Dict(
		:xgridvisible => false,
		:ygridvisible => false,
		:ylabel => "t last peri",
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
		scatter!(x, Δt_lmc; plot_kwargs...)
		for i in eachindex(idx)
			scatter!(x[idx[i]], Δt_lmc[idx[i]]; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)

	fig
end

# ╔═╡ e9311003-b6af-4e09-907f-c3388d59287d
let
	fig = Figure(size=(600, 600))
	plot_kwargs = Dict(
		:color => :black,
		:alpha => 0.1,
		:markersize => 1,
	)

	
	orbit_points_kwargs = [Dict(
		:color => COLORS[i],
		:alpha => 1,
		:markersize => 10,
		:label => orbit_labels[i]
	) for i in eachindex(idx)
		]
	
	ax_kwargs = Dict(
		:xgridvisible => false,
		:ygridvisible => false,
		:ylabel => "t last peri",
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
		scatter!(x, Δt_mw; plot_kwargs...)
		for i in eachindex(idx)
			scatter!(x[idx[i]], Δt_mw[idx[i]]; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)

	fig
end

# ╔═╡ 43d43f63-4c13-4b23-950e-ada59aa86bc9
let
	
	for sym in columns
	    x = [getproperty(o, sym) for o in observations]
	    y = apos

		
	    p = scatter(x, y, alpha=0.1, 
			axis=(; xlabel=String(sym), ylabel="apocentre / kpc")
		)
	    @info p 
	end

end

# ╔═╡ 69e77193-29cc-4304-98a1-44828eaedf9f
md"""
# Validation
"""

# ╔═╡ 6b95d3b2-38db-4376-83b5-8c6e6f1fdfa2
let
	fig = Figure(size=(600, 600))

	ax_kwargs = Dict(
		:xgridvisible => false,
		:ygridvisible => false,
		:ylabel => "density",
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
	    x = getproperty.(observations, sym)
		
		    stephist!(x, bins=50, normalization=:pdf, label="MC samples", color=:black)
	
		
		x_mod = LinRange(minimum(x), maximum(x), 1000)

		mu_exp = getproperty(obs, sym)
		err_exp = getproperty(err, sym)
		y_mod = normal_dist.(x_mod, mu_exp, err_exp)
		lines!(x_mod, y_mod, label="expected gaussian")
			
		axislegend(labelsize=10, padding=(6, 6, 6, 6), patchlabelgap=1, patchsize=(6, 6))
	end

	

	save(joinpath(fig_dir, "peri_mc_orbits_corr.pdf"), fig)
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


		μ = getproperty(obs, sym)
		σ = getproperty(err, sym)
		x_mod = LinRange(μ - 3σ, μ + 3σ, 10_000)
		y_mod = normal_dist.(x_mod, μ, σ) 
		lines!(x_mod, y_mod, label="expected")
			
		axislegend()

		@info fig
	end

end

# ╔═╡ de2f3380-90df-48f5-ba60-8417e91f4818
function median_residual(observations)
	for sym in [:distance, :pmra, :pmdec, :radial_velocity, :ra, :dec]
		md = median(getproperty.(observations, sym))
		res = (md - getproperty(obs, sym) ) / getproperty(err, sym)
		@printf "Δ ln %-15s  = %6.2f \t \n"  sym res
		@printf "     %-15s  = %6.2f \t \n"  sym md
	end
end

# ╔═╡ 0bb7cc30-69f2-496f-965f-0c331a928de0
median(peris)

# ╔═╡ 13f03c68-0edb-418b-b1f7-7a7f24659be2
let 
	peri1 = lguys.quantile(perilmc, 0.5p_value)
	peri2 = lguys.quantile(perilmc, 1.5p_value)

	idx_s = @. peri1 < perilmc < peri2

	println(median(perilmc[idx_s]))
	median_residual(observations[idx_s])
end

# ╔═╡ 4ee33ce2-c00a-4fcf-b7fc-b78c1677c9e4
let 
	peri1 = lguys.quantile(peris, 0.5p_value)
	peri2 = lguys.quantile(peris, 1.5p_value)

	idx_s = @. peri1 < peris < peri2

	println(median(peris[idx_s]))
	median_residual(observations[idx_s])
end

# ╔═╡ 34104429-05e0-40a6-83e5-078dbe346504
let
	peri2 = lguys.quantile(perilmc, 1 - 0.5p_value)
	peri1 = lguys.quantile(perilmc, 1 - 1.5p_value)

	idx_s = @. peri1 < perilmc < peri2

	median_residual(observations[idx_s])
end

# ╔═╡ e7fc6024-374b-422d-837b-26448e06e1db
observations[1].radial_velocity

# ╔═╡ 853b4153-7e5d-498b-bccd-26729d5371d9
obs.radial_velocity

# ╔═╡ e5825c4a-b446-44a3-8fd5-d94664965aca
median_residual(observations[[1]])

# ╔═╡ 16f4ac20-d8cf-4218-8c01-c15e04e567fb
md"""
# The selected orbits
"""

# ╔═╡ d31f91e8-6db6-4771-9544-8e54a816ecc1
begin

	positions = [lguys.extract_vector(out, :positions, i) for i in idx]
	velocities = [lguys.extract_vector(out, :velocities, i) for i in idx]
	accelerations = [lguys.extract_vector(out, :accelerations, i) for i in idx]
end

# ╔═╡ d3de3a33-8f9a-44a5-9851-1992c88f04ae
begin
	Φs_ext = [lguys.extract(out, :Φs_ext, i) for i in idx]
	Φs = [lguys.extract(out, :Φs, i) for i in idx]
end

# ╔═╡ 5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
begin
	rs = lguys.calc_r.(positions)
	vs = lguys.calc_r.(velocities)
	accs = lguys.calc_r.(accelerations)
end

# ╔═╡ e5d40e2f-ac47-4827-853d-2f94bc39a624
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="radius / kpc",
		limits=(nothing, nothing, 0, nothing)
	)

	for i in eachindex(idx)
		lines!(-out.times * lguys.T2GYR, rs[i], label=orbit_labels[i])
	
		#hlines!([peris[idx[i]], apos[idx[i]]], linestyle=:dot)
	end

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 2c073544-d24c-4583-99e9-30a57d82b154
rs_lmc = lguys.calc_r.(positions, [lmc_pos])

# ╔═╡ f3a53066-2edc-4c22-9132-ef0ec4c66db0
times_v

# ╔═╡ ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
positions

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
lguys.Plots.plot_xyz(positions..., labels=orbit_labels)

# ╔═╡ aec104ef-a4ca-4d69-aa72-1657c83c037d
lguys.Plots.plot_xyz([p .- lmc_pos for p in positions]..., labels=orbit_labels)

# ╔═╡ 34efe422-e672-4d17-a558-ce32fb704a8e
lguys.Plots.plot_xyz(velocities..., units=" / km / s")

# ╔═╡ 756f6d16-bc47-47d8-a090-32b4d7b3f3a0
times = out.times * lguys.T2GYR

# ╔═╡ 52e8c357-1269-4165-a6d7-d7e3c129f123
idx_perilmc = [argmin(r) for r in rs_lmc]

# ╔═╡ c82d023f-8123-4b20-ace1-cfb48dd94a5e
idx_peri = [argmin(r) for r in rs]

# ╔═╡ 57a8d1c8-3940-4430-8b46-375fb2bf1695
let
	x = positions[1][1, :]
	y = positions[1][2, :]
	z = positions[1][3, :]
	R = @. sqrt(x^2 + y^2)

	fig = Figure()
	ax = Axis(fig[1, 1])
	 p = plot!(R, z, color=times)

	scatter!(R[idx_peri[1]], z[idx_peri[1]])
	scatter!(R[idx_perilmc[1]], z[idx_perilmc[1]])
	
	Colorbar(fig[1, 2], p)

	fig
end

# ╔═╡ a2ea5f2e-16b0-4b20-8831-83a61f813750
let
	x = positions[1][1, :] .- lmc_pos[1, :]
	y = positions[1][2, :] .- lmc_pos[2, :]
	z = positions[1][3, :] .- lmc_pos[3, :]
	R = @. sqrt(x^2 + y^2)

	fig = Figure()
	ax = Axis(fig[1, 1])
	p = lines!(R, z, color=times)

	scatter!(R[idx_peri[1]], z[idx_peri[1]])
	scatter!(R[idx_perilmc[1]], z[idx_perilmc[1]])
	Colorbar(fig[1, 2], p)

	fig
end

# ╔═╡ 385eb53c-011b-431e-993b-ca12fece56e2
let
	fig = Figure()
	ax_acc = Axis(fig[1, 1], 
		yticksmirrored = false,
		xlabel="time / Gyr", 
		ylabel="Scl acceleration"
	)
	ax_dist = Axis(fig[1, 1], yaxisposition=:right,
		yticksmirrored = false,
		ylabel="distance",
		
	)
	
	hidexdecorations!(ax_dist)
	linkxaxes!(ax_acc, ax_dist)

	i = 1
	lines!(ax_acc, out.times * lguys.T2GYR, accs[1],
		color="black"
	)

	lines!(ax_dist, out.times * lguys.T2GYR, rs[1],
	)
	vlines!(ax_dist, times[idx_peri[i]], alpha=0.3)

	lines!(ax_dist, out.times * lguys.T2GYR, rs_lmc[1],
	)

	vlines!(ax_dist, times[idx_perilmc[i]], alpha=0.3)

	
	lguys.Plots.hide_grid!(ax_acc)
	lguys.Plots.hide_grid!(ax_dist)
	
	fig
end

# ╔═╡ ad078920-225d-436e-835b-d87a9db53c49
let
	fig = Figure()
	ax = Axis(fig[1,1])

	for i in eachindex(idx)
		scatter!(rs[i], vs[i], color=out.times)
	end

	fig
end

# ╔═╡ 8f8ea990-6ae8-45c2-a4f5-b7fc7a425b1d
let
	scatter(rs[1], accs[1])
end

# ╔═╡ 5d5f719c-e5b7-4e05-94b0-71523da46b66
let 
	fig = Figure()
	ax = Axis(fig[1,1])
	for i in 1:length(idx)
		scatter!(positions[i][1, :], accs[i])
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

# ╔═╡ fa4e7992-1ac6-4d71-a923-8b3cf81d0030
let 
	fig = Figure()
	ax = Axis(fig[1,1])
	
	for i in 1:length(idx)
		scatter!(out.times, (accs[i] .- a_exp.(rs[i])) ./ (a_exp.(rs[i])) )
	end
	fig
end 

# ╔═╡ 35f0ea14-a945-4745-910c-365b730676c5
let 
	fig = Figure()
	ax = Axis(fig[1,1])
	
	for i in 1:length(idx)
		scatter!(out.times, Φs_ext[i, :] .- phi_exp.(rs[i]))
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
for i in 1:length(idx)
	fname = "orbit$i.csv"
	vel = -reverse(velocities[i], dims=2)
	pos = reverse(positions[i], dims=2)
	t = -reverse(out.times)

	t0 = 1
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

# ╔═╡ b8c9823f-ca6b-48bf-9140-40440562dac0
import TOML

# ╔═╡ 94405d26-1eca-4dcc-b331-b47e9903b07d
times

# ╔═╡ 1152cd63-baab-426a-b464-b10857eed4ec
for i in 1:length(idx)
	fname = "orbit$i.toml"
	o = observations[idx[i]]
	properties = Dict(
		"ra" => o.ra,
		"dec" => o.dec,
		"distance" => o.distance,
		"pmra" => o.pmra,
		"pmdec" => o.pmdec,
		"radial_velocity" => o.radial_velocity,
		"distance_err" => err.distance,
		"pmra_err" => err.pmra,
		"pmdec_err" => err.pmdec,
		"radial_velocity_err" => err.radial_velocity
	)

	println("saving to $fname")
	open(fname, "w") do f
		TOML.print(f, properties)
	end
end

# ╔═╡ 519a88f0-8e2d-4c09-83e0-3cc2ee147e35
function get_initial_t(j)
	i = idx[j]
	N = length(rs[j])
	t_ini = -1
	for t in N-2:-1:1
		if diff(rs[j])[t] >= 0 && diff(rs[j])[t+1] <= 0
			t_ini = t
			break
		end
	end
	return t_ini
end
	

# ╔═╡ 5316884b-3971-4ca7-9106-f638241d3388
get_initial_t(1)

# ╔═╡ de1e5245-0946-47cd-8e2c-ba1914cfeb74
begin 
	# orbit info
	for i in 1:length(idx)
		t =1
		@printf "orbit: \t\t %i\n" i
		@printf "index: \t\t %i\n" idx[i]
		
		@printf "pericentre:\t %0.1f\n" peris[idx[i]]
		@printf "apocentre: \t %0.1f\n" apos[idx[i]]

		@printf "time of first apocentre: %f \n" out.times[end] - out.times[t]
		@printf "radius of first apocentre: %f\n" rs[i][t]
		@printf "intial position: [%f, %f, %f]\n" positions[i][:, t]...
		@printf "intial velocity: [%f, %f, %f]\n" -lguys.V2KMS* velocities[i][:, t]...
		@printf "final position: [%f, %f, %f]\n" positions[i][:, 1]...
		@printf "final velocity: [%f, %f, %f]\n" -lguys.V2KMS * velocities[i][:, 1]...

		o = observations[idx[i]]
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
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═d975d00c-fd69-4dd0-90d4-c4cbe73d9754
# ╠═a7ce5b0c-84a6-4d63-94f1-68e7a0d9e758
# ╠═1e50dbfe-09e5-4f42-83c3-a8291b8e1b1a
# ╠═0146ee17-de5f-4877-aaa6-83a898e01416
# ╠═3b83205d-91c1-481e-9305-0d59bc692135
# ╟─88536e86-cf2a-4dff-ae64-514821957d40
# ╠═26d616da-95ec-4fb9-b9a8-2f095d74c722
# ╠═9a22d47b-8474-4596-b418-de33eb07c627
# ╠═cd2c1a82-8b0b-4abc-8831-f621e5ec483a
# ╠═b9f469ed-6e4e-41ee-ac75-1b5bfa0a114a
# ╠═17a63cc8-84f4-4248-a7b0-c8378454b1f7
# ╠═e3084451-1492-49c5-b5dc-8a085e513640
# ╠═413d4e5d-c9cd-4aca-be1e-d132b2bd616d
# ╠═bbf3f229-fc3c-46ae-af28-0f8bd81e7d32
# ╠═fb6debf2-0161-477f-b29b-5a0f1f70f340
# ╠═4481a5e1-d635-4fea-a5c5-c85f0d6df62f
# ╠═384be6a6-f9d9-47e0-9792-aef6689dcbdb
# ╠═a1b2edd5-20e4-456f-a8e5-4634b0a433c2
# ╠═b0209748-88d7-4b60-b6b9-7a367a5d1e58
# ╠═b2fffcba-c474-430c-8100-037a0d352d69
# ╠═850c29d0-f395-448f-87be-d788210594c2
# ╠═e5f728b8-8412-4f57-ad38-a0a35bb08a48
# ╠═46348ecb-ee07-4b6a-af03-fc4f2635f57b
# ╠═6651d141-f6ca-4e8f-a785-5b14275b367c
# ╠═73af5376-3b38-4d5f-b6e2-1cafba27dabd
# ╟─1acef60e-60d6-47ba-85fd-f9780934788b
# ╠═ca1c236e-795a-408b-845b-9c13bc838619
# ╠═7b7b8357-ebf1-4257-8a38-272b245842cc
# ╠═6e658f28-0100-41b8-af8a-bb832baa1d12
# ╠═c43412d1-3a9b-4cd0-901b-05048c0b07ec
# ╠═88cbd5b4-29a9-4c2f-ad9e-657460d6314c
# ╠═7d7bc811-ced1-44a1-8378-d6c20b49d0e9
# ╠═7901c3ae-91b4-4fc1-8dc3-9b617c268533
# ╠═4cbbebec-f231-4af4-a10f-e01d090ea2ca
# ╠═46b4242b-8af7-4233-8ecf-d86740b4c884
# ╠═92aac8e8-d599-4a1e-965e-e653bc54c509
# ╠═471a5501-bc7b-4114-be1a-9088b4e674cd
# ╠═68b0383a-3d5a-4b94-967c-f0e31e8a0ce1
# ╠═2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
# ╠═8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
# ╠═5f11f6ab-c9ab-4600-acca-a0bb84d81a12
# ╠═ede3836c-740d-4ac7-bbc7-3165981a1878
# ╠═44660b2f-6220-473b-bb2f-07e23b176491
# ╠═d3063e30-2cb3-4f1b-8546-0d5e81d90d9f
# ╠═643676b7-ae9b-4d5b-b99d-f8bd338cf8a6
# ╠═c48b4e73-480e-4a50-b5fc-db5f6c5b040e
# ╠═99cbc08e-5b0e-4748-9727-667e0eca2f74
# ╠═456acc60-e7e1-428e-bd5d-bc72a18a146f
# ╠═e9311003-b6af-4e09-907f-c3388d59287d
# ╠═43d43f63-4c13-4b23-950e-ada59aa86bc9
# ╟─69e77193-29cc-4304-98a1-44828eaedf9f
# ╠═6b95d3b2-38db-4376-83b5-8c6e6f1fdfa2
# ╠═ac81acd8-4a78-4230-bc70-3b78a861b618
# ╠═de2f3380-90df-48f5-ba60-8417e91f4818
# ╠═0bb7cc30-69f2-496f-965f-0c331a928de0
# ╠═13f03c68-0edb-418b-b1f7-7a7f24659be2
# ╠═4ee33ce2-c00a-4fcf-b7fc-b78c1677c9e4
# ╠═34104429-05e0-40a6-83e5-078dbe346504
# ╠═e7fc6024-374b-422d-837b-26448e06e1db
# ╠═853b4153-7e5d-498b-bccd-26729d5371d9
# ╠═e5825c4a-b446-44a3-8fd5-d94664965aca
# ╟─16f4ac20-d8cf-4218-8c01-c15e04e567fb
# ╠═e5d40e2f-ac47-4827-853d-2f94bc39a624
# ╠═d31f91e8-6db6-4771-9544-8e54a816ecc1
# ╠═d3de3a33-8f9a-44a5-9851-1992c88f04ae
# ╠═5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
# ╠═2c073544-d24c-4583-99e9-30a57d82b154
# ╠═f3a53066-2edc-4c22-9132-ef0ec4c66db0
# ╠═ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
# ╠═aec104ef-a4ca-4d69-aa72-1657c83c037d
# ╠═57a8d1c8-3940-4430-8b46-375fb2bf1695
# ╠═a2ea5f2e-16b0-4b20-8831-83a61f813750
# ╠═34efe422-e672-4d17-a558-ce32fb704a8e
# ╠═756f6d16-bc47-47d8-a090-32b4d7b3f3a0
# ╠═385eb53c-011b-431e-993b-ca12fece56e2
# ╠═52e8c357-1269-4165-a6d7-d7e3c129f123
# ╠═c82d023f-8123-4b20-ace1-cfb48dd94a5e
# ╠═ad078920-225d-436e-835b-d87a9db53c49
# ╠═8f8ea990-6ae8-45c2-a4f5-b7fc7a425b1d
# ╠═5d5f719c-e5b7-4e05-94b0-71523da46b66
# ╠═09bbae0d-ca3e-426d-b77c-69dd68ca42cc
# ╠═fa4e7992-1ac6-4d71-a923-8b3cf81d0030
# ╠═35f0ea14-a945-4745-910c-365b730676c5
# ╟─2e7c1798-4066-4c46-b5ed-732263728ac0
# ╟─f6b27164-ee7c-428b-aefb-75e89d178f3e
# ╠═8b818798-69fb-481d-ade1-9fd436b1f281
# ╟─5fdd8307-d528-4cd7-a5e4-1f15aba75cd5
# ╠═2dfe9a85-6553-4632-81e0-33c148fd1102
# ╠═5f45e7c7-e447-48bf-ade4-38f516df2dad
# ╠═b8c9823f-ca6b-48bf-9140-40440562dac0
# ╠═94405d26-1eca-4dcc-b331-b47e9903b07d
# ╠═1152cd63-baab-426a-b464-b10857eed4ec
# ╠═519a88f0-8e2d-4c09-83e0-3cc2ee147e35
# ╠═5316884b-3971-4ca7-9106-f638241d3388
# ╠═de1e5245-0946-47cd-8e2c-ba1914cfeb74
