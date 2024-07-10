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

# ╔═╡ 8b21a4a9-ab04-4f24-b260-9b66d4cfc58e
using LilGuys

# ╔═╡ cf70db72-c259-4edf-9ebf-e3f3696c0f3d
using StatsBase

# ╔═╡ 170116fb-9f01-4a10-b989-93ed9a200d48
using HDF5

# ╔═╡ cddb5789-8738-4a31-bc24-e382d1742c52
md"""
# Initial stars
A quick notebook to check the initial stellar profile from a simulation
"""

# ╔═╡ 7eb3e35f-c2a5-499e-b884-85fb59060ec5
md"""
# Input
"""

# ╔═╡ 405c2a84-cfaf-469f-8eaa-0765f30a21de
model_dir = "/arc7/home/dboyea/sculptor/isolation/1e6/halos"

# ╔═╡ d7a04cc7-369e-4687-b423-deda779f1c57
name = "V50"

# ╔═╡ dcd21d08-cd8b-4d8c-a42d-fc0979093157
stars_dir = "/arc7/home/dboyea/sculptor/isolation/1e6/stars/"

# ╔═╡ 5a808597-89d7-4ba3-b330-69d654d0ef7b
stars_name = "exp2d_rs0.1"

# ╔═╡ 0ada1949-5322-471c-9a03-8b83921bac92
md"""
# File loading
"""

# ╔═╡ 920546bd-4838-413c-b687-f891a7f5e985
import TOML

# ╔═╡ 900fda29-c87e-41e4-a46c-b7faea21c0ce
params = TOML.parsefile(joinpath(model_dir, "$name.used.toml"))

# ╔═╡ 3dd35dd4-8e3c-458b-a6ce-b1c957266ce4
halo = lguys.NFW(; lguys.dict_to_tuple(params["profile"])...)

# ╔═╡ 954a44b3-afdf-4be4-bb74-ff50f33c9468
lguys.calc_M200(halo)

# ╔═╡ 204de061-5d08-4c52-94cd-7b40c710886e
stars_params = TOML.parsefile(joinpath(stars_dir, "$stars_name.toml"))

# ╔═╡ d3490582-d9d9-41df-accb-de687ba6ba8a
Mss = stars_params["profile_kwargs"]["M"] * params["M_scale"]

# ╔═╡ bda4d698-7dc5-40ae-bddc-e29d11a8daa8
R_s_s = stars_params["profile_kwargs"]["R_s"] * params["r_scale"]

# ╔═╡ 9b14ad87-426e-4486-a605-677436fd0d64
prof = lguys.Exp2D(; M=Mss, R_s=R_s_s)

# ╔═╡ 80605831-e74b-41f0-928c-02067ac9f2fc
h5open(joinpath(stars_dir, "$(stars_name)_stars.hdf5")) do f
	global probabilities = f["probabilities"][:]
	global p_idx = f["index"][:]
end

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
begin 
	println("loading model from $model_dir")
	snap = lguys.Snapshot(joinpath(model_dir, "$name.hdf5"))

	snap.weights = probabilities[snap.index]
end

# ╔═╡ 0948c690-d9ee-4ac6-b5d7-6c434dda0d40
md"""
# Plots
"""

# ╔═╡ a49d1735-203b-47dd-81e1-500ef42b054e
md"""
phase space distribution of star particles initial and final snapshot
"""

# ╔═╡ 72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log radius / kpc", ylabel="velocity (km/s)" )
	Arya.hist2d!(ax, log10.(lguys.calc_r(snap.positions)), lguys.calc_r(snap.velocities) * V2KMS, bins=100, weights=probabilities[snap.index])


	fig
end

# ╔═╡ a35b5f3d-ed9e-48f9-b96f-0a3c00ff2410
md"""
the dark matter distribution of  the snapshot (initial and final)
"""

# ╔═╡ d57fdbe3-d610-41d0-a613-d8ac3ed42ede
lguys.Plots.ProjectedDensity(snap, limits=((-1., 1.), (-1, 1)), bins=50)

# ╔═╡ b2597e2f-eebf-4a80-9e61-552a46f3b771
lguys.calc_E_tot(snap)

# ╔═╡ 4769a811-c89a-4b37-9e1f-599e2d1e8afa
virial_ratio = -2 * sum(lguys.calc_E_spec_kin(snap)) ./ sum(0.5 * lguys.calc_radial_discrete_Φ(snap))

# ╔═╡ 6550afcc-40fa-4386-a414-b87825ab6a12
# ╠═╡ disabled = true
#=╠═╡
using StatsBase
  ╠═╡ =#

# ╔═╡ 24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
md"""
stellar distribution of initial and final snapshot
"""

# ╔═╡ e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
function plot_ρ_dm!(snap; kwargs...)
	pos = lguys.extract_vector(snap, :positions)
	mass = lguys.extract(snap, :masses)
	rs = lguys.calc_r(pos)
	r, ρ = lguys.calc_ρ_hist(rs, 300, weights=mass)
	lines!(log10.(lguys.midpoint(r)), log10.(ρ); kwargs...)
end

# ╔═╡ 27f8deff-96ae-4d9a-a110-d146ac34965a
begin 
	function plot_ρ_dm(snap)
		plot()
		plot_ρ_dm!(snap)
	end

end

# ╔═╡ 60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log\, r / \textrm{kpc}", ylabel = L"\log\, \rho_\textrm{DM}\quad [10^{10} M_\odot / \textrm{kpc}^3]",
	limits=(-3, 3, -12, 0))
	plot_ρ_dm!(snap, label="initial")

	# ρ_0 = M_s / (4 * π * R_s^3)

	# ρ_nfw(r) = ρ_0  / (r/R_s) * 1/(r/R_s + 1)^2

	log_r = LinRange(-2, 3, 1000)
	y = log10.(lguys.calc_ρ.(halo, 10 .^ log_r))
	lines!(log_r, y, label="expected", color="black", linestyle=:dot)

	axislegend(ax)
	fig

end

# ╔═╡ 4e45e756-8a9c-43b4-aac7-2016347f5afb
let
	skip = 10
	snap_i = snap
	idx = 1:skip:length(snap_i)
	r = sort(lguys.calc_r(snap_i))[idx]
	
	M = idx

	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"\log\;r\;/\;\textrm{kpc}",
		ylabel="log number of particles",
		limits=((-3, 3), nothing)
	)

	ax.yticks = 0:2:7
	lines!(log10.(r), log10.(M))

	log_r = LinRange(-2, 3, 1000)
	y = log10.(lguys.calc_M.(halo, 10 .^ log_r) ./ snap.masses[1])
	lines!(log_r, y, label="expected", color="black", linestyle=:dot)

	fig
end

# ╔═╡ f01efc63-c4ac-45ae-8dac-209819a6249e
lguys.calc_M(halo, 2)

# ╔═╡ b9572089-54bd-41ea-b45d-38d3a0e9d13f
let
	fig, ax = FigAxis(
		xlabel = "vx / km s",
		ylabel = "density"
	)
	
	lguys.Plots.vx_hist_fit!(snap, bins=50, direction=2)

	fig
end

# ╔═╡ e93cb11b-3655-42aa-819e-97967000b1df
let
	fig = Figure()
	ax = lguys.Plots.Axis_rho(fig[1, 1], limits=((-1, 2), (-25, 5)))

	lguys.Plots.plot_ρ_s!(snap)
	x = LinRange(-1, 1, 1000)

	y = log10.(lguys.calc_ρ.(prof, 10 .^ x))

	lines!(x, y)
	fig
end

# ╔═╡ 6c857d5b-41c0-4028-aad7-15a5d8d8244a
TOML.print

# ╔═╡ Cell order:
# ╟─cddb5789-8738-4a31-bc24-e382d1742c52
# ╠═6e08e538-bc82-11ee-1a75-d97f506d18c5
# ╠═8b21a4a9-ab04-4f24-b260-9b66d4cfc58e
# ╠═cf70db72-c259-4edf-9ebf-e3f3696c0f3d
# ╟─7eb3e35f-c2a5-499e-b884-85fb59060ec5
# ╠═405c2a84-cfaf-469f-8eaa-0765f30a21de
# ╠═d7a04cc7-369e-4687-b423-deda779f1c57
# ╠═dcd21d08-cd8b-4d8c-a42d-fc0979093157
# ╠═5a808597-89d7-4ba3-b330-69d654d0ef7b
# ╟─0ada1949-5322-471c-9a03-8b83921bac92
# ╠═920546bd-4838-413c-b687-f891a7f5e985
# ╠═900fda29-c87e-41e4-a46c-b7faea21c0ce
# ╠═3dd35dd4-8e3c-458b-a6ce-b1c957266ce4
# ╠═954a44b3-afdf-4be4-bb74-ff50f33c9468
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╠═204de061-5d08-4c52-94cd-7b40c710886e
# ╠═d3490582-d9d9-41df-accb-de687ba6ba8a
# ╠═bda4d698-7dc5-40ae-bddc-e29d11a8daa8
# ╠═9b14ad87-426e-4486-a605-677436fd0d64
# ╠═170116fb-9f01-4a10-b989-93ed9a200d48
# ╠═80605831-e74b-41f0-928c-02067ac9f2fc
# ╟─0948c690-d9ee-4ac6-b5d7-6c434dda0d40
# ╟─a49d1735-203b-47dd-81e1-500ef42b054e
# ╠═72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
# ╟─a35b5f3d-ed9e-48f9-b96f-0a3c00ff2410
# ╠═d57fdbe3-d610-41d0-a613-d8ac3ed42ede
# ╠═b2597e2f-eebf-4a80-9e61-552a46f3b771
# ╠═4769a811-c89a-4b37-9e1f-599e2d1e8afa
# ╠═6550afcc-40fa-4386-a414-b87825ab6a12
# ╟─24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
# ╠═e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
# ╠═27f8deff-96ae-4d9a-a110-d146ac34965a
# ╠═60f8d0cd-ca8e-457c-98b9-1ee23645f9dd
# ╠═4e45e756-8a9c-43b4-aac7-2016347f5afb
# ╠═f01efc63-c4ac-45ae-8dac-209819a6249e
# ╠═b9572089-54bd-41ea-b45d-38d3a0e9d13f
# ╠═e93cb11b-3655-42aa-819e-97967000b1df
# ╠═6c857d5b-41c0-4028-aad7-15a5d8d8244a
