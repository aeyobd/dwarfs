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

# ╔═╡ fdbe2478-32eb-4264-b484-a7b072a4028d
sim_dir = ENV["DWARFS_ROOT"] * "/simulations/sculptor/mc_orbits/iorio/"

# ╔═╡ cd64bd83-fee6-4fed-94ab-131ed49ebef0
if !isdefined(Main, :obs)
	include(joinpath(sim_dir, "sample.jl"))
end

# ╔═╡ 380b8e74-2837-4620-b600-ae6674b35d16
p_value = 0.001349898031630093 # 3sigma

# ╔═╡ 17a63cc8-84f4-4248-a7b0-c8378454b1f7
idx = [1, 2
]

# ╔═╡ e5f728b8-8412-4f57-ad38-a0a35bb08a48
orbit_labels = ["FH18", "EH18"]

# ╔═╡ 46348ecb-ee07-4b6a-af03-fc4f2635f57b
fig_dir = "./figures"

# ╔═╡ 9a22d47b-8474-4596-b418-de33eb07c627
begin 
	out = lguys.Output(".");

	df_peris_apos = lguys.read_fits("peris_apos.fits")
	snap = out[1] |> sort_snap

	@assert all(snap.index  .== df_peris_apos.index) "snapshot and peri apo index must match"
end

# ╔═╡ 2559bc4c-3bad-4d70-b02f-8fe48b772fb6
frames_df = lguys.read_fits(ENV["DWARFS_ROOT"] * "/simulations/sculptor/mc_orbits/iorio/gc_frames.fits")

# ╔═╡ bd89c450-084c-4a28-9c37-b4f03e5b2ca9
frames = [lguys.GalactocentricFrame(
	d=row.d, z_sun=row.z_sun, 
	v_sun=[row.v_sun_x, row.v_sun_y, row.v_sun_z],
	ra=row.ra,
	dec=row.dec
) for row in eachrow(frames_df)]

# ╔═╡ aa6809d1-45ee-4973-a80e-61066389db35
μ_frame = Dict(k => dist.μ for (k, dist) in frame_distributions)

# ╔═╡ ab397f1a-51c6-45eb-a8cd-d47759361f9e
σ_frame = Dict(k => dist.σ for (k, dist) in frame_distributions)

# ╔═╡ fb6debf2-0161-477f-b29b-5a0f1f70f340
[snap.positions[:, 1]; snap.velocities[:, 1] * lguys.V2KMS]

# ╔═╡ 4481a5e1-d635-4fea-a5c5-c85f0d6df62f
peris = df_peris_apos.pericenter

# ╔═╡ 413d4e5d-c9cd-4aca-be1e-d132b2bd616d
peri_qs = lguys.quantile(peris, [p_value, 1-p_value])

# ╔═╡ bbf3f229-fc3c-46ae-af28-0f8bd81e7d32
peris[idx]

# ╔═╡ 384be6a6-f9d9-47e0-9792-aef6689dcbdb
apos = df_peris_apos.apocenter

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
	axis=(; xlabel="initial galactocentric velocity / km/s", ylabel="count")
)

# ╔═╡ 2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
ϵ = -lguys.calc_K_spec(out[1]) .- out[1].Φs_ext

# ╔═╡ 8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
extrema(ϵ)

# ╔═╡ 75967d46-614b-422e-a0fd-2660c9f070fc
frames

# ╔═╡ 5f11f6ab-c9ab-4600-acca-a0bb84d81a12
begin
	points = [lguys.Galactocentric(
		snap.positions[:, i]*lguys.R2KPC, 
		-snap.velocities[:, i]*lguys.V2KMS,
		frame=frames[i]
	)
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

# ╔═╡ 8932a1ea-fedc-47c6-82e2-7d9b0580eb6f
let
	orbit_points_kwargs = [Dict(
		:color => COLORS[i],
		:alpha => 1,
		:markersize => 10,
		:label => orbit_labels[i]
	) for i in eachindex(idx)
	]
	for sym in names(frames_df)
	    x = frames_df[:, sym]
	    y = peris

		fig = Figure(size=(300, 300),)
			
		ax = Axis(fig[1, 1],
			aspect=1,
			xlabel=String(sym), ylabel="pericentre / kpc",
		)
		
	    p = scatter!(x, y, alpha=0.1, markersize=1, color=:black,)
		lguys.Plots.hide_grid!(ax)

		for i in eachindex(idx)
			scatter!(x[idx[i]], y[idx[i]]; orbit_points_kwargs[i]...)
		end
	    @info fig
	end

end

# ╔═╡ c48b4e73-480e-4a50-b5fc-db5f6c5b040e
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

# ╔═╡ 43d43f63-4c13-4b23-950e-ada59aa86bc9
let
	
	for sym in columns
	    x = [getproperty(o, sym) for o in observations]
	    y = peris

		
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

	    
	    x_mod = LinRange(minimum(x), maximum(x), 1000)

			
		axislegend()

		@info fig
	end

end

# ╔═╡ de2f3380-90df-48f5-ba60-8417e91f4818
function median_residual(observations)
	for sym in [:distance, :pmra, :pmdec, :radial_velocity]
		md = median(getproperty.(observations, sym))
		res = (md - getproperty(obs, sym) ) / getproperty(err, sym)
		@printf "Δ ln %-15s  = %6.2f \t \n"  sym res
	end

	for sym in [:d, :v_sun_x, :v_sun_y, :v_sun_z]
		md = median(frames_df[:, sym])
		res = (md - μ_frame[sym]) / σ_frame[sym]
		@printf "Δ ln %-15s  = %6.2f \t \n"  sym res
	end
end

# ╔═╡ 4ee33ce2-c00a-4fcf-b7fc-b78c1677c9e4
let 
	peri1 = lguys.quantile(peris, p_value * 0.8)
	peri2 = lguys.quantile(peris, p_value * 1.2)

	idx_s = @. peri1 < peris < peri2

	median_residual(observations[idx_s])
end

# ╔═╡ 34104429-05e0-40a6-83e5-078dbe346504
let
	peri2 = lguys.quantile(peris, 1 - p_value * 0.8)
	peri1 = lguys.quantile(peris, 1 - p_value * 1.2)

	idx_s = @. peri1 < peris < peri2

	median_residual(observations[idx_s])
end

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
	Φs_ext = [lguys.extract(out, :Φs_ext, i) for i in idx]

end

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
		title = "Iorio orbits"
	)

	for i in eachindex(idx)
		lines!(out.times * lguys.T2GYR, rs[i], label=orbit_labels[i])
	
		hlines!([peris[idx[i]], apos[idx[i]]], linestyle=:dot)
	end

	Legend(fig[1, 2], ax)
	lguys.Plots.hide_grid!(ax)
	fig
end

# ╔═╡ ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
positions

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
let

	fig = lguys.Plots.plot_xyz(positions..., labels=orbit_labels)

	resize_to_layout!(fig)
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

	axislegend()
	fig
end

# ╔═╡ 34efe422-e672-4d17-a558-ce32fb704a8e
lguys.Plots.plot_xyz(velocities..., units=" / km / s")

# ╔═╡ ad078920-225d-436e-835b-d87a9db53c49
let
	fig = Figure()
	ax = Axis(fig[1,1])

	for i in eachindex(idx)
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
	ax = Axis(fig[1,1])
	
	for i in 1:length(idx)
		scatter!(out.times, Φs_ext[i] .- phi_exp.(rs[i]))
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

# ╔═╡ 1152cd63-baab-426a-b464-b10857eed4ec
for i in 1:length(idx)
	fname = "orbit$i.toml"
	o = observations[idx[i]]
	frame = frames_df[i, :]
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

	for name in names(frames_df)
		properties["gc_$name"] = frame[name]
	end
	println("saving to $fname")
	println(properties)
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
	

# ╔═╡ 5f45e7c7-e447-48bf-ade4-38f516df2dad
for i in 1:length(idx)
	fname = "orbit$i.csv"
	vel = -reverse(velocities[i], dims=2)
	pos = reverse(positions[i], dims=2)
	t0 = size(pos, 2) - get_initial_t(i)
	println(t0)
	t = out.times[end] .- reverse(out.times)

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

# ╔═╡ 5316884b-3971-4ca7-9106-f638241d3388
get_initial_t(1)

# ╔═╡ de1e5245-0946-47cd-8e2c-ba1914cfeb74
begin 
	# orbit info
	for i in 1:length(idx)
		t = get_initial_t(i)
		@printf "orbit: \t\t %i\n" i
		@printf "index: \t\t %i\n" idx[i]
		
		@printf "pericentre:\t %0.1f\n" peris[idx[i]]
		@printf "apocentre: \t %0.1f\n" apos[idx[i]]

		@printf "time of first apocentre: %f \n" out.times[end] - out.times[t]
		@printf "radius of first apocentre: %f\n" rs[i][t]
		@printf "intial position: [%f, %f, %f]\n" positions[i][:, t]...
		@printf "intial velocity: [%f, %f, %f]\n" -1* velocities[i][:, t]...
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
# ╠═3be492fa-0cb9-40f1-a39e-e25bf485fd3e
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═d975d00c-fd69-4dd0-90d4-c4cbe73d9754
# ╠═cd64bd83-fee6-4fed-94ab-131ed49ebef0
# ╠═3b83205d-91c1-481e-9305-0d59bc692135
# ╟─88536e86-cf2a-4dff-ae64-514821957d40
# ╠═26d616da-95ec-4fb9-b9a8-2f095d74c722
# ╠═38fd9158-7807-4a26-b534-b0bcb5dbb0e5
# ╠═fdbe2478-32eb-4264-b484-a7b072a4028d
# ╠═380b8e74-2837-4620-b600-ae6674b35d16
# ╠═17a63cc8-84f4-4248-a7b0-c8378454b1f7
# ╠═e5f728b8-8412-4f57-ad38-a0a35bb08a48
# ╠═46348ecb-ee07-4b6a-af03-fc4f2635f57b
# ╠═9a22d47b-8474-4596-b418-de33eb07c627
# ╠═2559bc4c-3bad-4d70-b02f-8fe48b772fb6
# ╠═bd89c450-084c-4a28-9c37-b4f03e5b2ca9
# ╠═aa6809d1-45ee-4973-a80e-61066389db35
# ╠═ab397f1a-51c6-45eb-a8cd-d47759361f9e
# ╠═413d4e5d-c9cd-4aca-be1e-d132b2bd616d
# ╠═bbf3f229-fc3c-46ae-af28-0f8bd81e7d32
# ╠═fb6debf2-0161-477f-b29b-5a0f1f70f340
# ╠═4481a5e1-d635-4fea-a5c5-c85f0d6df62f
# ╠═384be6a6-f9d9-47e0-9792-aef6689dcbdb
# ╟─1acef60e-60d6-47ba-85fd-f9780934788b
# ╠═ca1c236e-795a-408b-845b-9c13bc838619
# ╠═46b4242b-8af7-4233-8ecf-d86740b4c884
# ╠═92aac8e8-d599-4a1e-965e-e653bc54c509
# ╠═471a5501-bc7b-4114-be1a-9088b4e674cd
# ╠═68b0383a-3d5a-4b94-967c-f0e31e8a0ce1
# ╠═2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
# ╠═8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
# ╠═75967d46-614b-422e-a0fd-2660c9f070fc
# ╠═5f11f6ab-c9ab-4600-acca-a0bb84d81a12
# ╠═ede3836c-740d-4ac7-bbc7-3165981a1878
# ╠═44660b2f-6220-473b-bb2f-07e23b176491
# ╠═d3063e30-2cb3-4f1b-8546-0d5e81d90d9f
# ╠═8932a1ea-fedc-47c6-82e2-7d9b0580eb6f
# ╠═c48b4e73-480e-4a50-b5fc-db5f6c5b040e
# ╠═43d43f63-4c13-4b23-950e-ada59aa86bc9
# ╟─69e77193-29cc-4304-98a1-44828eaedf9f
# ╠═6b95d3b2-38db-4376-83b5-8c6e6f1fdfa2
# ╠═ac81acd8-4a78-4230-bc70-3b78a861b618
# ╠═de2f3380-90df-48f5-ba60-8417e91f4818
# ╠═4ee33ce2-c00a-4fcf-b7fc-b78c1677c9e4
# ╠═34104429-05e0-40a6-83e5-078dbe346504
# ╠═e5825c4a-b446-44a3-8fd5-d94664965aca
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
# ╠═1152cd63-baab-426a-b464-b10857eed4ec
# ╠═519a88f0-8e2d-4c09-83e0-3cc2ee147e35
# ╠═5316884b-3971-4ca7-9106-f638241d3388
# ╠═de1e5245-0946-47cd-8e2c-ba1914cfeb74