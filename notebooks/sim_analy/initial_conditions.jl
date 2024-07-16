### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 6e08e538-bc82-11ee-1a75-d97f506d18c5
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya

	import LilGuys as lguys
end

# ╔═╡ 01165464-35b2-43e9-9d52-d06fa9edf3cf
using LilGuys

# ╔═╡ f979f2a8-3420-4ede-a739-7d727dfdf818
md"""
# Initial conditions
Check the initial DM halo. Loads a snapshot and checks if the halo has the expected density & velocity profile
"""

# ╔═╡ 920546bd-4838-413c-b687-f891a7f5e985
import TOML

# ╔═╡ 4833d7f7-2c0a-4e15-a47b-cb2474351283
12 * 3600

# ╔═╡ b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
md"""
To make this both a cml utility and interactive, we take inputs in the following cells
"""

# ╔═╡ 7eb3e35f-c2a5-499e-b884-85fb59060ec5
md"""
# Input
"""

# ╔═╡ 405c2a84-cfaf-469f-8eaa-0765f30a21de
model_dir = "/arc7/home/dboyea/sculptor/isolation/1e5/"

# ╔═╡ d7a04cc7-369e-4687-b423-deda779f1c57
name = "out/snapshot_101"

# ╔═╡ 3dd35dd4-8e3c-458b-a6ce-b1c957266ce4
begin 
	#params = TOML.parsefile(joinpath(model_dir, "$name.toml"))
	params = TOML.parsefile(joinpath(model_dir, "halo.toml"))
	
	halo = lguys.NFW(; lguys.dict_to_tuple(params["profile"])...)

	# halo = lguys.NFW(; M_s=1/lguys.A_NFW(1), r_s=1)
	# halo = lguys.NFW(M_s=0.289934, r_s=2.76)
end

# ╔═╡ d3313f08-7e4e-43b8-b55d-ea099d031bfe
begin 
	using DataFrames, CSV

	zeno_prof = CSV.read("/astro/dboyea/dwarfs/zeno/profiles/nfw.csv", DataFrame)

	zeno_prof.Radius *= halo.r_s
	zeno_prof.Mass *= halo.M_s * LilGuys.A_NFW(1)
	zeno_prof.Density *=  halo.M_s / halo.r_s^3 * LilGuys.A_NFW(1)

	zeno_prof[!, :V_circ] = sqrt.(zeno_prof.Mass ./ zeno_prof.Radius)

	zeno_prof
end

# ╔═╡ 54a2a708-d8ba-4c5c-9e67-ac656dd8e9f4
lguys.calc_M(halo, 1)

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
begin 
	println("loading model from $model_dir")
	snap = lguys.Snapshot(joinpath(model_dir, "$name.hdf5"))

	snap.x_cen = lguys.centroid(snap.positions)
	snap.v_cen = lguys.centroid(snap.velocities)

	println(snap.x_cen)
	println(snap.v_cen)

	snap.positions .-= snap.x_cen
	snap.velocities .-= snap.v_cen

	snap.x_cen = zeros(3)
	snap.v_cen = zeros(3)
	
	lguys.calc_r(snap)
end

# ╔═╡ 0ccb9018-d88c-4cec-a8da-625be1289bfe
snap.x_cen

# ╔═╡ 5ebe92b8-602e-42be-8751-58898b7323b0
snap.v_cen

# ╔═╡ 14e3b593-17b9-4acd-a9cb-d5923662a02c
prof = lguys.calc_profile(snap, filt_bound=false)

# ╔═╡ 0e89851e-763f-495b-b677-b664501a17ef
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"\log \; r / \textrm{kpc}", ylabel=L"$V_\textrm{circ}$ / km s$^{-1}$")

	rc, Vc = lguys.calc_v_circ(snap)
	lines!(log10.(rc), Vc * lguys.V2KMS, label="initial")

	log_r = LinRange(-2, 2.5, 1000)

	lines!(log_r, lguys.V2KMS * lguys.calc_v_circ.(halo, 10 .^ log_r))

	lines!(prof.log_r, prof.v_circ .* V2KMS)
	lines!(log10.(zeno_prof.Radius), zeno_prof.V_circ .* V2KMS)
	scatter!(log10.(lguys.calc_r_circ_max(halo)), lguys.calc_v_circ_max(halo) * lguys.V2KMS)
	fig
end

# ╔═╡ a49d1735-203b-47dd-81e1-500ef42b054e
md"""
phase space distribution of star particles initial and final snapshot
"""

# ╔═╡ 72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log radius / kpc", ylabel="velocity (km/s)" )
	Arya.hist2d!(ax, log10.(lguys.calc_r(snap.positions)), lguys.calc_r(snap.velocities) * lguys.V2KMS, bins=100)


	fig
end

# ╔═╡ a35b5f3d-ed9e-48f9-b96f-0a3c00ff2410
md"""
the dark matter distribution of  the snapshot (initial and final)
"""

# ╔═╡ b9746093-0f2f-4478-82ba-00911c8fcceb
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc", title="initial")

	bins = LinRange(-10, 10, 30)
	Arya.hist2d!(ax, snap.positions[1, :], snap.positions[2, :], bins = bins, colorscale=log10)

	fig
end

# ╔═╡ 24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
md"""
stellar distribution of initial and final snapshot
"""

# ╔═╡ e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
function plot_ρ_dm!(snap; kwargs...)
	pos = lguys.extract_vector(snap, :positions)
	mass = lguys.extract(snap, :masses)
	rs = lguys.calc_r(pos)
	r, ρ = lguys.calc_ρ_hist(rs, 100, weights=mass)
	lines!(log10.(lguys.midpoint(r)), log10.(ρ); kwargs...)
end

# ╔═╡ 60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log\, r / \textrm{kpc}", ylabel = L"\log\, \rho_\textrm{DM}\quad [10^{10} M_\odot / \textrm{kpc}^3]",
	#limits=(-2, 5, -12, 2)
	)

	log_r = LinRange(-2, 3, 1000)
	y = log10.(lguys.calc_ρ.(halo, 10 .^ log_r))
	lines!(log_r, y, label="NFW", color="black", linestyle=:dot)

	lines!(log10.(zeno_prof.Radius), 
		log10.(zeno_prof.Density ), 
		label="zeno", linestyle=:dash)
	
	lines!(prof.log_r, log10.(prof.rho))

	axislegend(ax)
	fig

end

# ╔═╡ aca95a0a-98e0-4b7a-bca2-e3c30f9df6e9


# ╔═╡ 84b3759e-a598-4afc-a2b4-ce841e80ff96
let
	skip = 10
	snap_i = snap
	idx = 1:skip:length(snap_i)
	r = sort(lguys.calc_r(snap_i))[idx]
	
	M = cumsum(snap_i.masses)[idx]

	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"\log\;r\;/\;\textrm{kpc}",
		ylabel="log M in",
		#limits=((-3, 3), nothing)
	)

	lines!(log10.(r), log10.(M))


	lines!(log10.(zeno_prof.Radius), 
		log10.(zeno_prof.Mass), 
		label="zeno", linestyle=:dash)

	
	log_r = LinRange(-2, 3, 1000)
	y = log10.(lguys.calc_M.(halo, 10 .^ log_r) )
	lines!(log_r, y, label="expected", color="black", linestyle=:dot)

	fig
end

# ╔═╡ 4e45e756-8a9c-43b4-aac7-2016347f5afb
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"\log\;r\;/\;\textrm{kpc}",
		ylabel="log number of particles",
		#limits=((-3, 3), nothing)
	)

	lines!(prof.log_r_bins[2:end], log10.(cumsum(prof.counts)))

	fig
end

# ╔═╡ f01efc63-c4ac-45ae-8dac-209819a6249e
lguys.calc_M(halo, 2)

# ╔═╡ 34d9fdea-8961-44ca-a92f-2f48a281f2cd
let
	fig, ax = FigAxis( ylabel="log counts / bin", xlabel="log r / kpc",
		limits=(nothing, (-0.5, 5))
	)
	
	scatter!(prof.log_r, log10.(prof.counts))

	fig
end

# ╔═╡ Cell order:
# ╟─f979f2a8-3420-4ede-a739-7d727dfdf818
# ╠═6e08e538-bc82-11ee-1a75-d97f506d18c5
# ╠═01165464-35b2-43e9-9d52-d06fa9edf3cf
# ╠═920546bd-4838-413c-b687-f891a7f5e985
# ╠═4833d7f7-2c0a-4e15-a47b-cb2474351283
# ╟─b7ef1dbd-1865-4ac3-a4d7-26fc9b443c45
# ╟─7eb3e35f-c2a5-499e-b884-85fb59060ec5
# ╠═405c2a84-cfaf-469f-8eaa-0765f30a21de
# ╠═d7a04cc7-369e-4687-b423-deda779f1c57
# ╠═3dd35dd4-8e3c-458b-a6ce-b1c957266ce4
# ╠═d3313f08-7e4e-43b8-b55d-ea099d031bfe
# ╠═0ccb9018-d88c-4cec-a8da-625be1289bfe
# ╠═5ebe92b8-602e-42be-8751-58898b7323b0
# ╠═54a2a708-d8ba-4c5c-9e67-ac656dd8e9f4
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╠═14e3b593-17b9-4acd-a9cb-d5923662a02c
# ╠═0e89851e-763f-495b-b677-b664501a17ef
# ╟─a49d1735-203b-47dd-81e1-500ef42b054e
# ╟─72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
# ╟─a35b5f3d-ed9e-48f9-b96f-0a3c00ff2410
# ╠═b9746093-0f2f-4478-82ba-00911c8fcceb
# ╟─24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
# ╠═e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
# ╠═60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
# ╠═aca95a0a-98e0-4b7a-bca2-e3c30f9df6e9
# ╠═84b3759e-a598-4afc-a2b4-ce841e80ff96
# ╠═4e45e756-8a9c-43b4-aac7-2016347f5afb
# ╠═f01efc63-c4ac-45ae-8dac-209819a6249e
# ╠═34d9fdea-8961-44ca-a92f-2f48a281f2cd
