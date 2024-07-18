### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ e502e544-4409-11ef-0edc-ff1c38cb579d
begin
	using Pkg; Pkg.activate()


	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ f6064906-de30-4fc4-9e6b-2e965c635aa4
models_dir = "/astro/dboyea/sculptor/orbits"

# ╔═╡ f6031981-ac78-4276-99a3-231720ab5931
load_profile(name) = LilGuys.Profiles3D(joinpath(models_dir, "$name/out/profiles.hdf5"))

# ╔═╡ 56b7f1f4-a44a-4c68-acbc-a56279e6de47
labels = ["V32, r=04", "V=50, r=10.", "V50, r=5"]

# ╔═╡ e6427b32-64d9-4b9b-a250-2a0ea6ae5c7e
model_names = ["orbit1", "orbit1_V50", "V50_r0.5"]

# ╔═╡ c1dbfc5c-4d3f-401a-93cd-784653b7f791
profiles = load_profile.(model_names)

# ╔═╡ ce879308-d7f4-42c0-9257-bcbfdc163e55
0.15319 * V2KMS

# ╔═╡ 9a85c819-0e42-45f9-85bb-025830612d75
function t_max_of(prof)
	return 2π * prof.r_circ_max / prof.v_circ_max
end

# ╔═╡ d014b999-4ab2-40a7-90c6-c95ae289b60f
function ρ_max_of(prof)
	return 3π / t_max_of(prof) ^ 2
end

# ╔═╡ a5f6ada4-ca00-4bfe-a5a0-bc739a5f4f0e
let
	fig, ax = FigAxis(
		xlabel="time / Gyr",
		ylabel = L"$v_\textrm{circ, max}$ / km\,s$^{-1}$"
	)


	for i in eachindex(profiles)
		
		vs = [p.v_circ_max for p in profiles[i].profiles]
		lines!(profiles[i].snapshot_index * T2GYR, (vs * V2KMS), label=labels[i])
	end

	axislegend()
	fig
end

# ╔═╡ f2ad6d7f-47bd-480e-9aa9-9d1e37ce597f
let
	fig, ax = FigAxis(
		xlabel="time / Gyr",
		ylabel = L"$t_\textrm{circ, max}$ / km\,s$^{-1}$"
	)


	for i in eachindex(profiles)
		
		ts = t_max_of.(profiles[i].profiles)
		lines!(profiles[i].snapshot_index * T2GYR, ts * T2GYR, label=labels[i])
	end

	axislegend()
	fig
end

# ╔═╡ 77494685-ed93-4216-a69b-8f89e5781378
let
	fig, ax = FigAxis(
		xlabel="time / Gyr",
		ylabel = L"$r_\textrm{circ, max}$ / kpc"
	)

	for profile in profiles
		vs = [p.r_circ_max for p in profile.profiles]
		lines!(profile.snapshot_index* T2GYR, vs)
	end
	
	fig
end

# ╔═╡ 79e1f031-a3a0-4c70-bb8c-32ab97f77d92
t_max_host = 300

# ╔═╡ 11192afb-1200-4995-a3a8-f4905dd23473
let
	fig, ax = FigAxis(
		xlabel="time / Gyr",
		ylabel = L"$t_\textrm{circ, max}$ / Gyr"
	)


	for i in eachindex(profiles)
		
		ts = t_max_of.(profiles[i].profiles)
		lines!(profiles[i].snapshot_index * T2GYR, ts * T2GYR, label=labels[i])
	end

	hlines!(t_max_host * T2GYR / 2)

	axislegend()
	fig
end

# ╔═╡ b8da09a6-bbf3-4489-b368-d0362eeede4d
R2KPC / T2GYR / V2KMS

# ╔═╡ 627ebe42-adc4-4939-9d7e-40375e93e7a6
let
	fig, ax = FigAxis(
		xlabel=L"$r_\textrm{circ max}$ / kpc",
		ylabel = L"$v_\textrm{circ, max}$ / kpc"
	)

	for profile in profiles
		vs = [p.v_circ_max for p in profile.profiles]

		rs = [p.r_circ_max for p in profile.profiles]
		lines!(log10.(rs), log10.(vs * V2KMS))
	end

	r = LinRange(1, 10, 100)
	v = 2π * r ./ t_max_host * 3/2

	lines!(log10.(r), log10.(v * V2KMS))
	
	fig
end

# ╔═╡ 1420b193-cc05-4b27-a7bd-b7aa331e18b3
let
	fig, ax = FigAxis(
		xlabel="time / Gyr",
		ylabel = "bound number of particles",
	)

	for i in eachindex(profiles)
		
		vs = [p.N_bound for p in profiles[i].profiles]
		lines!(profiles[i].snapshot_index * T2GYR, log10.(vs), label=labels[i])
	end

	axislegend()
	
	fig
end

# ╔═╡ Cell order:
# ╠═e502e544-4409-11ef-0edc-ff1c38cb579d
# ╠═f6064906-de30-4fc4-9e6b-2e965c635aa4
# ╠═f6031981-ac78-4276-99a3-231720ab5931
# ╠═56b7f1f4-a44a-4c68-acbc-a56279e6de47
# ╠═e6427b32-64d9-4b9b-a250-2a0ea6ae5c7e
# ╠═c1dbfc5c-4d3f-401a-93cd-784653b7f791
# ╠═ce879308-d7f4-42c0-9257-bcbfdc163e55
# ╠═9a85c819-0e42-45f9-85bb-025830612d75
# ╠═d014b999-4ab2-40a7-90c6-c95ae289b60f
# ╠═a5f6ada4-ca00-4bfe-a5a0-bc739a5f4f0e
# ╠═f2ad6d7f-47bd-480e-9aa9-9d1e37ce597f
# ╠═11192afb-1200-4995-a3a8-f4905dd23473
# ╠═77494685-ed93-4216-a69b-8f89e5781378
# ╠═79e1f031-a3a0-4c70-bb8c-32ab97f77d92
# ╠═b8da09a6-bbf3-4489-b368-d0362eeede4d
# ╠═627ebe42-adc4-4939-9d7e-40375e93e7a6
# ╠═1420b193-cc05-4b27-a7bd-b7aa331e18b3
