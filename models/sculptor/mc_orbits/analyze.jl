### A Pluto.jl notebook ###
# v0.19.43

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
using Statistics

# ╔═╡ a7ce5b0c-84a6-4d63-94f1-68e7a0d9e758
if !@isdefined obs # read in sample (but only once)
	include("sample.jl")
	obs
end

# ╔═╡ cd9dea36-93f8-4461-85b1-87030d9367bb
pwd()

# ╔═╡ 9c7e5bb4-8db0-4527-a8ec-331e2aed958b
begin
	filename = "out/combined.hdf5"
	out = lguys.Output(filename);
end

# ╔═╡ 4481a5e1-d635-4fea-a5c5-c85f0d6df62f
begin
	df = CSV.read("peris_apos.csv", DataFrame)
	df_idx = sortperm(df[!, :index])
	peris = df[df_idx, :pericenter]
	apos = df[df_idx, :apocenter];
end

# ╔═╡ 9c259aa6-bc9e-4351-ac68-eccd05ff4c3d
115/lguys.A_NFW(9.545)

# ╔═╡ 4d4a18fc-8990-4dcc-97c5-c9e01708ea2e
begin 
	snap = out[1]
	snap = snap[sortperm(snap.index)]
end

# ╔═╡ fb6debf2-0161-477f-b29b-5a0f1f70f340
[snap.positions[:, 1]; snap.velocities[:, 1] * lguys.V2KMS]

# ╔═╡ fcadcc96-1da1-4e6f-9d3b-2c56e55488b7
let
	hist(peris)
end

# ╔═╡ 46b4242b-8af7-4233-8ecf-d86740b4c884
	hist(apos)


# ╔═╡ e10615ab-0ece-4ebb-81d0-47ffdd6ea6c5
begin
	percies = [16, 84]
	quantiles = lguys.percentile(peris, [16, 84])
	idx_qs = [findfirst(isequal(p), peris) for p in quantiles]
	idx = vcat([1], idx_qs)
end

# ╔═╡ 9fad353f-751e-4577-b704-aa2eadeb0969
quantiles

# ╔═╡ 92b28dd8-c353-4959-b181-a843367b3223
hist(lguys.calc_r(snap.velocities))

# ╔═╡ 471a5501-bc7b-4114-be1a-9088b4e674cd
hist(lguys.calc_r(snap.positions),
	axis=(; xlabel="radius")
)

# ╔═╡ 2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
ϵ = lguys.calc_E_spec_kin(out[1]) + out[1].Φs_ext

# ╔═╡ 8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
maximum(ϵ)

# ╔═╡ 5f11f6ab-c9ab-4600-acca-a0bb84d81a12
begin
	points = [lguys.Galactocentric(snap.positions[:, i]*lguys.R2KPC, -snap.velocities[:, i]*lguys.V2KMS) for i in 1:length(snap)]
	observations = lguys.transform.(lguys.ICRS, points)
end


# ╔═╡ 92aac8e8-d599-4a1e-965e-e653bc54c509
dists = getproperty.(observations, :distance)

# ╔═╡ ede3836c-740d-4ac7-bbc7-3165981a1878
normal_dist(x, μ, σ) = 1/√(2π) * 1/σ * exp(-(x-μ)^2/2σ^2)

# ╔═╡ 44660b2f-6220-473b-bb2f-07e23b176491
columns = [:ra, :dec, :distance, :pmra, :pmdec, :radial_velocity]

# ╔═╡ ac81acd8-4a78-4230-bc70-3b78a861b618
let

	for sym in columns
	
	    μ = getproperty(obs, sym)
	    σ = getproperty(err, sym)
		
		fig = Figure()
		ax = Axis(fig[1,1], 
			xlabel=String(sym),
			ylabel="density",
			#limits=((μ - 5σ, μ + 5σ), nothing),
		)
		
	    x = getproperty.(observations, sym)
		
	    stephist!(x, normalization=:pdf)

	    
	    x_mod = LinRange(μ - 3σ, μ + 3σ, 1000)
	    y_mod = normal_dist.(x_mod, μ, σ)
		
	    lines!(x_mod, y_mod)
		@info fig
	end

end

# ╔═╡ d3063e30-2cb3-4f1b-8546-0d5e81d90d9f
let
	
	for sym in columns
	    x = [getproperty(o, sym) for o in observations]
	    y = peris

		
	    p = scatter(x, y, alpha=0.1,
			axis=(; xlabel=String(sym), ylabel="pericentre / kpc")
		)
	    @info p 
	end

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

# ╔═╡ 4ee33ce2-c00a-4fcf-b7fc-b78c1677c9e4
begin 
	peri1 = lguys.percentile(peris, 10)
	peri2 = lguys.percentile(peris, 20)

	idx_s = @. peri1 < peris < peri2

	for sym in [:distance, :pmra, :pmdec, :radial_velocity, :ra, :dec]
		md = median(getproperty.(observations[idx_s], sym))
		println("d $sym = $(md - getproperty(obs, sym))")
	end
end

# ╔═╡ e5825c4a-b446-44a3-8fd5-d94664965aca
for f in fieldnames(lguys.ICRS)
	a = getproperty(observations[1], f)
	b = getproperty(obs, f)
	println("error of $f:\t", (a-b)/b)
end

# ╔═╡ a644a33e-57a4-4af4-98be-1e0e84948069


# ╔═╡ d31f91e8-6db6-4771-9544-8e54a816ecc1
begin
	
	positions = lguys.extract_vector(out, :positions, idx)
	velocities = lguys.extract_vector(out, :velocities, idx)
	accelerations = lguys.extract_vector(out, :accelerations, idx)
	positions = [positions[:, i, :] for i in 1:length(idx)]
	velocities = [velocities[:, i, :] for i in 1:length(idx)]
	accelerations = [accelerations[:, i, :] for i in 1:length(idx)]
	Φs_ext = lguys.extract(out, :Φs_ext, idx)
	Φs = lguys.extract(out, :Φs, idx)


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
		ylabel="radius / kpc"
	)

	for i in eachindex(idx)
		lines!(out.times * lguys.R0, rs[i])
	
		hlines!([peris[idx[i]], apos[idx[i]]])
	end

	fig
end

# ╔═╡ ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
positions

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
lguys.plot_xyz(positions...)

# ╔═╡ 57a8d1c8-3940-4430-8b46-375fb2bf1695
let
	x = positions[1][1, :]
	y = positions[1][2, :]
	z = positions[1][3, :]
	R = @. sqrt(x^2 + y^2)

	plot(R, z)
end

# ╔═╡ 34efe422-e672-4d17-a558-ce32fb704a8e
lguys.plot_xyz(velocities...)

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
		scatter!(out.times, accs[i] .- a_exp.(rs[i]))
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
scatter(snap.positions)

# ╔═╡ 5fdd8307-d528-4cd7-a5e4-1f15aba75cd5
scatter(-snap.velocities .* lguys.V0)

# ╔═╡ 2dfe9a85-6553-4632-81e0-33c148fd1102
reverse(out.times)

# ╔═╡ b8c9823f-ca6b-48bf-9140-40440562dac0
import TOML

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

# ╔═╡ ca334fc0-3840-4182-8bbf-b78375bb7ed9
positions

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

# ╔═╡ 8e6959cf-b93e-48e8-adf4-5578fe43ab53
lguys.T2GYR

# ╔═╡ Cell order:
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═a7ce5b0c-84a6-4d63-94f1-68e7a0d9e758
# ╠═cd9dea36-93f8-4461-85b1-87030d9367bb
# ╠═9c7e5bb4-8db0-4527-a8ec-331e2aed958b
# ╠═fb6debf2-0161-477f-b29b-5a0f1f70f340
# ╠═4481a5e1-d635-4fea-a5c5-c85f0d6df62f
# ╠═9c259aa6-bc9e-4351-ac68-eccd05ff4c3d
# ╠═4d4a18fc-8990-4dcc-97c5-c9e01708ea2e
# ╠═fcadcc96-1da1-4e6f-9d3b-2c56e55488b7
# ╠═46b4242b-8af7-4233-8ecf-d86740b4c884
# ╠═e10615ab-0ece-4ebb-81d0-47ffdd6ea6c5
# ╠═9fad353f-751e-4577-b704-aa2eadeb0969
# ╠═92aac8e8-d599-4a1e-965e-e653bc54c509
# ╠═e5d40e2f-ac47-4827-853d-2f94bc39a624
# ╠═92b28dd8-c353-4959-b181-a843367b3223
# ╠═471a5501-bc7b-4114-be1a-9088b4e674cd
# ╠═2c094ad9-23b4-40f1-a1ec-3b61bf96bffe
# ╠═8501b0a7-a71f-41b4-b6f6-5f34b37f24d5
# ╠═5f11f6ab-c9ab-4600-acca-a0bb84d81a12
# ╠═ede3836c-740d-4ac7-bbc7-3165981a1878
# ╠═44660b2f-6220-473b-bb2f-07e23b176491
# ╠═ac81acd8-4a78-4230-bc70-3b78a861b618
# ╠═d3063e30-2cb3-4f1b-8546-0d5e81d90d9f
# ╠═43d43f63-4c13-4b23-950e-ada59aa86bc9
# ╠═d975d00c-fd69-4dd0-90d4-c4cbe73d9754
# ╠═4ee33ce2-c00a-4fcf-b7fc-b78c1677c9e4
# ╠═e5825c4a-b446-44a3-8fd5-d94664965aca
# ╟─a644a33e-57a4-4af4-98be-1e0e84948069
# ╠═d31f91e8-6db6-4771-9544-8e54a816ecc1
# ╠═5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
# ╠═ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
# ╠═57a8d1c8-3940-4430-8b46-375fb2bf1695
# ╠═34efe422-e672-4d17-a558-ce32fb704a8e
# ╠═ad078920-225d-436e-835b-d87a9db53c49
# ╠═8f8ea990-6ae8-45c2-a4f5-b7fc7a425b1d
# ╠═5d5f719c-e5b7-4e05-94b0-71523da46b66
# ╠═09bbae0d-ca3e-426d-b77c-69dd68ca42cc
# ╠═fa4e7992-1ac6-4d71-a923-8b3cf81d0030
# ╠═35f0ea14-a945-4745-910c-365b730676c5
# ╟─2e7c1798-4066-4c46-b5ed-732263728ac0
# ╠═f6b27164-ee7c-428b-aefb-75e89d178f3e
# ╠═5fdd8307-d528-4cd7-a5e4-1f15aba75cd5
# ╠═2dfe9a85-6553-4632-81e0-33c148fd1102
# ╠═5f45e7c7-e447-48bf-ade4-38f516df2dad
# ╠═b8c9823f-ca6b-48bf-9140-40440562dac0
# ╠═1152cd63-baab-426a-b464-b10857eed4ec
# ╠═ca334fc0-3840-4182-8bbf-b78375bb7ed9
# ╠═519a88f0-8e2d-4c09-83e0-3cc2ee147e35
# ╠═5316884b-3971-4ca7-9106-f638241d3388
# ╠═de1e5245-0946-47cd-8e2c-ba1914cfeb74
# ╠═8e6959cf-b93e-48e8-adf4-5578fe43ab53
