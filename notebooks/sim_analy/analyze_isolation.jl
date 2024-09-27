### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 6e08e538-bc82-11ee-1a75-d97f506d18c5
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya

	import LilGuys as lguys
end

# ╔═╡ 07ce4cef-df87-494b-8f16-02489e88bca6
using LilGuys

# ╔═╡ 374489bc-627f-4fc9-9734-7c49456710ac
begin 
	import DataFrames, CSV
	using HDF5
	import NaNMath as nm
end

# ╔═╡ 7e548032-56a4-4b5c-bc47-ea86bb4cf917
md"""
This notebook creates plots of the dark matter properties for a halo in isolation.

In particular, we plot the density profiles, velocity profiles, and evolution of the calculated centre of the snapshot.
"""

# ╔═╡ 96c91860-f3cc-4531-a8cf-39c85887b394
import TOML

# ╔═╡ 6dc6b811-5993-42d1-a6ac-07920aa4f564
save = CairoMakie.save

# ╔═╡ 7eb3e35f-c2a5-499e-b884-85fb59060ec5
md"""
# Inputs
"""

# ╔═╡ 405c2a84-cfaf-469f-8eaa-0765f30a21de
name = "/arc7/home/dboyea/sculptor/isolation/1e6/fiducial"

# ╔═╡ a29c993a-c7eb-4b57-a474-50bdbd0ce1ec
halo = lguys.load_profile(joinpath(name, "halo.toml"))

# ╔═╡ 6435d239-4a83-4990-a25e-cd03bb0e2022
softening = 0.044 # kpc

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
out = lguys.Output(joinpath(name, "out/combined.hdf5"))

# ╔═╡ dd31d3ee-7fdf-46f5-b213-45faae93ae5e
idxs = [1, 3, 6, length(out)]

# ╔═╡ 38ab4838-960e-4321-b70c-6f14584e9e27
snaps = [out[i] for i in idxs]

# ╔═╡ 97f89831-00e6-49a2-a712-ac47fd2dee47
out.times[idxs] * lguys.T2GYR

# ╔═╡ c4fc0a87-c75f-4528-aa20-75cc7aff856c
figure_dir = joinpath(name, "figures/")

# ╔═╡ 327f790d-e652-48b2-92e5-e2dffd5b15e2
mkpath(figure_dir)

# ╔═╡ 5de2aa65-86ed-46fc-99c6-2cb53ca6f5c5
profs = LilGuys.read_structs_from_hdf5(joinpath(name, "out/profiles.hdf5"), LilGuys.MassProfile3D)

# ╔═╡ 97e98ab8-b60b-4b48-b465-a34a16858f88
md"""
# Initial/Final
"""

# ╔═╡ 9a9f4dc1-3573-41ee-be1a-eee39d3371b0
fit = lguys.fit_v_r_circ_max(snaps[end])

# ╔═╡ 0e89851e-763f-495b-b677-b664501a17ef
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"\log \; r / \textrm{kpc}", ylabel=L"V_\textrm{circ}")


	for i in eachindex(snaps)
		lines!(profs[i].log_r, profs[i].v_circ * V2KMS, label="snapshot $(idxs[i])" )
	end


	scatter!(log10.(fit.r_circ_max), fit.v_circ_max * V2KMS)

	V_nfw(x) = lguys.calc_v_circ(halo, x)
	log_r = LinRange(-2, 2.5, 1000)
	y = V_nfw.(10 .^ log_r)
	lines!(log_r, y * V2KMS, label="NFW", linestyle=:dot, color=:black)

	axislegend(position=:lt)

	save(figure_dir * "v_circ.pdf", fig)
	fig
end

# ╔═╡ 7299dbaf-b332-4e53-85de-7acbfa0c3853
import StatsBase: percentile

# ╔═╡ a49d1735-203b-47dd-81e1-500ef42b054e
md"""
phase space distribution of star particles initial and final snapshot
"""

# ╔═╡ 06cafe6c-98ed-42ce-8b6c-b3e12afab896
import DensityEstimators as DE

# ╔═╡ e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
function plot_ρ_dm!(snap; N=100, kwargs...)
	pos = lguys.extract_vector(snap, :positions)
	mass = lguys.extract(snap, :masses)
	rs = lguys.calc_r(pos)
	bins = 10 .^ LinRange(log10(minimum(rs)), log10(maximum(rs)), N)
	r, ρ = lguys.calc_ρ_hist(rs, bins, weights=mass)
	lines!(log10.(lguys.midpoints(r)), log10.(ρ); kwargs...)
end

# ╔═╡ 2edbc22e-1432-441a-9fba-ef256321fe25


# ╔═╡ 264aedc1-b624-475e-af28-1d31b533839d
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log\, r / \textrm{kpc}", ylabel = L"\log\, \rho_\textrm{DM}\quad [10^{10} M_\odot / \textrm{kpc}^3]",
	limits=(-2.5, 4, -20, 1))

	for i in idxs
		lines!(profs[i].log_r, log10.(profs[i].rho), label="snapshot $(profs.snapshot_index[i])")
	end

	log_r = LinRange(-2, 4, 1000)
	y = log10.(lguys.calc_ρ.(halo, 10 .^ log_r))
	lines!(log_r, y, label="trunc-NFW", color="black", linestyle=:dot)

	vlines!(log10(softening), color="grey")
	axislegend(ax)
	
	save(figure_dir * "dm_density.pdf", fig)

	fig

end

# ╔═╡ 4e45e756-8a9c-43b4-aac7-2016347f5afb
let
	skip = 10
	snap = snaps[1]
	idx = 1:skip:length(snap)
	r = sort(lguys.calc_r(snap))[idx]
	
	M = idx

	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"\log\;r\;/\;\textrm{kpc}",
		ylabel="log number of particles",
		limits=((-3, 3), nothing)
	)

	ax.yticks = 0:2:7
	lines!(log10.(r), log10.(M))

	fig
end

# ╔═╡ e33d56a7-7a0e-4fa9-8f0d-041b43584d59
sum(lguys.calc_r(snaps[end]) .< 0.1)

# ╔═╡ 34d9fdea-8961-44ca-a92f-2f48a281f2cd
let
	fig, ax = FigAxis( ylabel="count", xlabel="log r / kpc", yscale=log10)
	
	hist!(log10.(lguys.calc_r(snaps[1])), label="", bins=100)

	fig
end

# ╔═╡ 5f1c61f9-50d4-43cb-aa78-fa85314f26b7
let
	fig, ax = FigAxis( ylabel="count", xlabel="e spec", yscale=log10, 
		limits=(-4, -1, 1, 1e6)
	)

	for i in eachindex(snaps)
		e = -lguys.calc_E_spec(snaps[i])
		e = e[e .> 0]
		stephist!(log10.(e), bins=100, label="$i")
	end
	
	fig
end

# ╔═╡ 1ce0aa3c-953e-4f7b-bf7e-7c815c505b5e
println(sum(lguys.calc_r(snaps[1]) .< softening))

# ╔═╡ fb0dec74-aaab-43a4-9b37-d13634c5dcad
md"""
# Time variation
"""

# ╔═╡ 8e87001c-4bc0-4580-a0b4-43fe24d92c99
skip = 30

# ╔═╡ ebdd5430-c7c1-4fc7-82f5-8acd8ca99070
# ╠═╡ disabled = true
#=╠═╡
let 
	skip = 30
	idx = 1:skip:length(out)

	Ls = hcat([lguys.calc_L_tot(snap) for snap in out[idx]]...)

	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="total angular momentum component"
	)

	t = out.times[idx] * T2GYR
	
	for i in 1:3
		scatter!(t, Ls[i, :], label=["x", "y", "z"][i])
		lines!(t, Ls[i, :], label=["x", "y", "z"][i])
	end

	axislegend(ax)
	fig
end
  ╠═╡ =#

# ╔═╡ d6bc6ccb-15aa-415a-aa48-8a3cff89749e
md"""
TODO: Why are these methods slightly different?
"""

# ╔═╡ 36b741d4-5a47-4d2e-8c93-0b29c53e2a0e
out[3].x_cen

# ╔═╡ 17435494-a855-4295-9ce7-e60d937c2aa8
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="total angular momentum component"
	)

	idx = profs.snapshot_index
	t = out.times[idx] * T2GYR
	Ls = [profs[i].L for i in eachindex(profs.snapshot_index)]
	Ls = hcat(Ls...)

	ii = 1:3:length(idx)
	t = t[ii]
	Ls = Ls[:, ii]
	
	for i in 1:3
		scatter!(t, Ls[i, :], label=["x", "y", "z"][i])
		lines!(t, Ls[i, :])
	end

	axislegend(ax)
	fig
end

# ╔═╡ 9de516f4-968f-44e9-b35a-74882a802dc0
profs.snapshot_index[1]

# ╔═╡ bfa7593c-4915-4e03-83e6-8f790de4c1a5
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="snapshot", ylabel="relative change in energy")
	idx = 1:skip:length(out)
	Es = [profs[i].E for i in eachindex(profs.snapshot_index)]
	scatter!((Es ./ Es[1]))

	save(figure_dir * "energy.pdf", fig)

	fig
end

# ╔═╡ e3a45a8e-cc52-4e9d-9db3-97109b59fc77
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="snapshot", ylabel="-W / 2T (virial ratio)")

	idx = 1:1skip:length(out)
	
	E_kin = [profs[i].K for i in eachindex(profs.profiles)]
	E_pot = -[profs[i].W for i in eachindex(profs.profiles)]

	scatter!(-E_pot ./ 2E_kin)
	hlines!(1, color=:black)

	save(figure_dir * "virial.pdf", fig)

	fig
end

# ╔═╡ 449493f7-6d80-407d-ab85-d734203245af
profs[1].K

# ╔═╡ 91a44ed4-8466-4a58-b3ff-1e7630b8ac8c
let
	fig = lguys.Plots.plot_xyz(out.x_cen)
	fig.content[1].title = "centre"

	save(figure_dir * "centre.pdf", fig)
	
	fig
end

# ╔═╡ e61c095e-a763-466b-b419-755fd0aadd0d
lguys.Plots.plot_xyz(out.v_cen * V2KMS, units=" / km s⁻¹")

# ╔═╡ b5c71290-d2de-424d-b026-f1ae15d7d86e
percentile(lguys.calc_r(out[1]), [5, 10, 50, 90, 95])

# ╔═╡ dc221349-eb61-4ace-8de3-a6c50249aca0
function find_radii_fracs(out, x_cen; skip=10) 
	rs = Vector[]
	Ms = Vector[]
	rs_s = Vector[]

	percens = [0.0001, 0.001, 0.01, 0.1, 1, 3, 10, 50, 90, 97]
	
	for i in 1:skip:length(out)
		r = lguys.calc_r(out[i])
		
		push!(rs, percentile(r, percens))
	end

	rs = hcat(rs...)
	rs_s = hcat(rs_s...)

	return percens, rs, rs_s

end

# ╔═╡ 34244a2e-9501-451c-bd77-bebfebde2a78
percens, rs, rs_s = find_radii_fracs(out, out.x_cen, skip=1)

# ╔═╡ ed206b2b-7ee6-4b77-a5bf-bd3dcc3f976f
# TODO: Add these calculations to the profiles3D

# ╔═╡ 967136d3-8d58-4fdc-9537-aa3a85a92528
times = out.times * T2GYR

# ╔═╡ f21cfe22-95f3-485d-902b-b022a41548c2
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel="log r containing fraction of DM")
	
	for i in eachindex(percens)
		label = "$(percens[i]) "
		y = log10.(rs[i, :])
		scatterlines!(times[1:1:end], y, 				
			label=label, color=i, colorrange=(1, length(percens))
		)
	end

	Legend(fig[1,2], ax, "f = ")

	save(figure_dir * "mass_fraction_w_time.pdf", fig)

	fig
end

# ╔═╡ 3b2bb553-0130-4c8a-80ad-6e1f7071a293
#lguys.Plots.plot_xyz(lguys.extract_vector(out, :positions, 100_000))

# ╔═╡ 5153654a-567f-4663-9b4a-6f08f9e49d1a
md"""
# 2D Density plots
"""

# ╔═╡ a35b5f3d-ed9e-48f9-b96f-0a3c00ff2410
md"""
the dark matter distribution of  the snapshot (initial and final)
"""

# ╔═╡ b9746093-0f2f-4478-82ba-00911c8fcceb
function plot_xy_dens(snap_i; r_max=10)
	fig = Figure()
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc")

	bins = LinRange(-r_max, r_max, 100)
	hist_kwargs = (; bins=bins, colorscale=log10, colorrange=(1, nothing))
	Arya.hist2d!(ax, snap_i.positions[1, :], snap_i.positions[2, :]; hist_kwargs...)
	
	fig
end

# ╔═╡ 3cc0df1f-dca3-4355-a969-3874aa414d1a
for snap in snaps
	@info plot_xy_dens(snap)
end

# ╔═╡ 2dedc9e1-9b68-4c76-97da-4aba767de32d
for snap in snaps
	@info plot_xy_dens(snap, r_max=300)
end

# ╔═╡ 24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
md"""
stellar distribution of initial and final snapshot
"""

# ╔═╡ e5ca8db2-2c3d-4b97-9242-ab1d2ebf51b3
function calc_phase_dens(snap; bins=100)
	rs = lguys.calc_r(snap)
	vs = lguys.calc_v(snap)
	xbins = 10 .^ DE.make_bins(log10.(rs), DE.calc_limits(log10.(rs)), bins)
	ybins = DE.make_bins(vs, DE.calc_limits(vs), bins)
	h = DE.histogram2d(rs, vs, (xbins, ybins), weights=snap.masses, limits=(xbins[1], xbins[end], ybins[1], ybins[end]))

	# since f is normalized to /dx/dy/dz, we do need to correct h by r^2 and v^2

	h.values ./= diff(4π/3 * xbins .^ 3)
	h.values ./= diff(4π/3 * ybins .^ 3)

	h
end

# ╔═╡ 8a097a37-a903-4627-ba25-0a1f0289955f
function plot_phase_dens!(snap)
	h = calc_phase_dens(snap)

	heatmap!(h, colorscale=log10, colorrange=(1e-15 * maximum(h.values), maximum(h.values)))
end

# ╔═╡ 9996a264-743f-4d39-a5fe-1cda5b99930b
function plot_phase_dens(snap)
	fig, ax = FigAxis(
		xscale=log10
	)

	h = plot_phase_dens!(snap)

	Colorbar(fig[1, 2], h)
	fig
end

# ╔═╡ 97b87dc9-1fcc-4b8c-a6bb-80cd573ed45a
for snap in snaps
	@info plot_phase_dens(snap)
end

# ╔═╡ 72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
function plot_rv_dens(snap_i)
	fig = Figure(size=(700, 300))
	ax = Axis(fig[1,1], xlabel="log radius / kpc", ylabel="velocity (km/s)" )

	hist_kwargs = (; bins=100, colorscale=log10, colorrange=(1, nothing))
	
	h = Arya.hist2d!(ax, log10.(lguys.calc_r(snap_i.positions)), lguys.calc_r(snap_i.velocities) * V2KMS; hist_kwargs...)

	Colorbar(fig[1, 2], h)

	fig
end

# ╔═╡ 3f625247-0d7d-4f6d-bf84-a40f19f1849a
plot_rv_dens(snaps[end])

# ╔═╡ 9b88e37a-93eb-4c03-b654-0cd9284fc811
function plot_phase_rvr!(snap; bins=100)
	x = lguys.calc_r(snap)
	v = lguys.calc_v_rad(snap)

	h = Arya.histogram2d(log10.(x), v *V2KMS, bins)

	heatmap!(h, colorscale=log10, colorrange=(1, maximum(h.values)))
end

# ╔═╡ 564c6bcb-04be-4a1c-8f01-f7d76be74eb8
function plot_phase_xvx!(snap; bins=100)
	x = asinh.(snap.positions[1, :])
	v = snap.velocities[1, :]

	h = Arya.histogram2d(x, v * V2KMS, bins)

	heatmap!(h, colorscale=log10, colorrange=(1, maximum(h.values)))
end

# ╔═╡ c9ffe8ad-97d2-40e7-8ba7-e27d3708d723
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="asinh x / kpc", ylabel="x velocity (km/s)" )

	h = plot_phase_xvx!(snaps[1])

	Colorbar(fig[1, 2], h)

	fig
end

# ╔═╡ 1a320f74-5cb7-44c5-8a59-c6fced771f52
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log radius / kpc", ylabel="radial velocity (km/s)" )

	h = plot_phase_rvr!(snaps[end])

	Colorbar(fig[1, 2], h)

	fig
end

# ╔═╡ Cell order:
# ╠═7e548032-56a4-4b5c-bc47-ea86bb4cf917
# ╠═6e08e538-bc82-11ee-1a75-d97f506d18c5
# ╠═07ce4cef-df87-494b-8f16-02489e88bca6
# ╠═374489bc-627f-4fc9-9734-7c49456710ac
# ╠═96c91860-f3cc-4531-a8cf-39c85887b394
# ╠═6dc6b811-5993-42d1-a6ac-07920aa4f564
# ╟─7eb3e35f-c2a5-499e-b884-85fb59060ec5
# ╠═405c2a84-cfaf-469f-8eaa-0765f30a21de
# ╠═a29c993a-c7eb-4b57-a474-50bdbd0ce1ec
# ╠═dd31d3ee-7fdf-46f5-b213-45faae93ae5e
# ╠═38ab4838-960e-4321-b70c-6f14584e9e27
# ╠═6435d239-4a83-4990-a25e-cd03bb0e2022
# ╠═97f89831-00e6-49a2-a712-ac47fd2dee47
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╠═c4fc0a87-c75f-4528-aa20-75cc7aff856c
# ╠═327f790d-e652-48b2-92e5-e2dffd5b15e2
# ╠═5de2aa65-86ed-46fc-99c6-2cb53ca6f5c5
# ╟─97e98ab8-b60b-4b48-b465-a34a16858f88
# ╠═9a9f4dc1-3573-41ee-be1a-eee39d3371b0
# ╠═0e89851e-763f-495b-b677-b664501a17ef
# ╠═7299dbaf-b332-4e53-85de-7acbfa0c3853
# ╟─a49d1735-203b-47dd-81e1-500ef42b054e
# ╠═06cafe6c-98ed-42ce-8b6c-b3e12afab896
# ╠═e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
# ╠═2edbc22e-1432-441a-9fba-ef256321fe25
# ╠═264aedc1-b624-475e-af28-1d31b533839d
# ╠═4e45e756-8a9c-43b4-aac7-2016347f5afb
# ╠═e33d56a7-7a0e-4fa9-8f0d-041b43584d59
# ╠═34d9fdea-8961-44ca-a92f-2f48a281f2cd
# ╠═5f1c61f9-50d4-43cb-aa78-fa85314f26b7
# ╠═1ce0aa3c-953e-4f7b-bf7e-7c815c505b5e
# ╟─fb0dec74-aaab-43a4-9b37-d13634c5dcad
# ╠═8e87001c-4bc0-4580-a0b4-43fe24d92c99
# ╠═ebdd5430-c7c1-4fc7-82f5-8acd8ca99070
# ╠═d6bc6ccb-15aa-415a-aa48-8a3cff89749e
# ╠═36b741d4-5a47-4d2e-8c93-0b29c53e2a0e
# ╠═17435494-a855-4295-9ce7-e60d937c2aa8
# ╠═9de516f4-968f-44e9-b35a-74882a802dc0
# ╠═bfa7593c-4915-4e03-83e6-8f790de4c1a5
# ╠═e3a45a8e-cc52-4e9d-9db3-97109b59fc77
# ╠═449493f7-6d80-407d-ab85-d734203245af
# ╠═91a44ed4-8466-4a58-b3ff-1e7630b8ac8c
# ╠═e61c095e-a763-466b-b419-755fd0aadd0d
# ╠═b5c71290-d2de-424d-b026-f1ae15d7d86e
# ╠═dc221349-eb61-4ace-8de3-a6c50249aca0
# ╠═34244a2e-9501-451c-bd77-bebfebde2a78
# ╠═ed206b2b-7ee6-4b77-a5bf-bd3dcc3f976f
# ╠═f21cfe22-95f3-485d-902b-b022a41548c2
# ╠═967136d3-8d58-4fdc-9537-aa3a85a92528
# ╠═3b2bb553-0130-4c8a-80ad-6e1f7071a293
# ╟─5153654a-567f-4663-9b4a-6f08f9e49d1a
# ╟─a35b5f3d-ed9e-48f9-b96f-0a3c00ff2410
# ╠═b9746093-0f2f-4478-82ba-00911c8fcceb
# ╠═3cc0df1f-dca3-4355-a969-3874aa414d1a
# ╠═2dedc9e1-9b68-4c76-97da-4aba767de32d
# ╟─24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
# ╠═e5ca8db2-2c3d-4b97-9242-ab1d2ebf51b3
# ╠═8a097a37-a903-4627-ba25-0a1f0289955f
# ╠═9996a264-743f-4d39-a5fe-1cda5b99930b
# ╠═97b87dc9-1fcc-4b8c-a6bb-80cd573ed45a
# ╠═72dfab8a-c6c8-4dcc-b399-a0cf6cb0dea0
# ╠═3f625247-0d7d-4f6d-bf84-a40f19f1849a
# ╠═9b88e37a-93eb-4c03-b654-0cd9284fc811
# ╠═564c6bcb-04be-4a1c-8f01-f7d76be74eb8
# ╠═c9ffe8ad-97d2-40e7-8ba7-e27d3708d723
# ╠═1a320f74-5cb7-44c5-8a59-c6fced771f52
