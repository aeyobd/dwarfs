### A Pluto.jl notebook ###
# v0.20.18

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

# ╔═╡ 6e08e538-bc82-11ee-1a75-d97f506d18c5
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya

	import LilGuys as lguys
end

# ╔═╡ 06d06ab0-4585-4cc6-b875-5cdd2e47657a
using PyFITS

# ╔═╡ 374489bc-627f-4fc9-9734-7c49456710ac
begin 
	import DataFrames, CSV
	using HDF5
	import NaNMath as nm
end

# ╔═╡ 7e232e3c-34b9-4b08-a297-7aa19aca644b
using PlutoUI

# ╔═╡ 4f1e87b3-f776-4550-a0f9-45eec1b79340
function notebook_inputs(; kwargs...)
	return PlutoUI.combine() do Child
		
		user_inputs = [
			md""" $(string(name)): $(
				Child(name, obj)
			)"""
			
			for (name, obj) in kwargs
		]
		
		md"""
		#### Inputs
		$(user_inputs)
		"""
	end
end

# ╔═╡ 03a1bd03-dc40-42c0-bbcb-4aa562f47bbf
@bind inputs confirm(notebook_inputs(;
	haloname = TextField(default="1e5/example"),
))

# ╔═╡ 7e548032-56a4-4b5c-bc47-ea86bb4cf917
md"""
This notebook creates plots of the dark matter properties for a halo in isolation.

In particular, we plot the density profiles, velocity profiles, and evolution of the calculated centre of the snapshot.
"""

# ╔═╡ 96c91860-f3cc-4531-a8cf-39c85887b394
import TOML

# ╔═╡ 7eb3e35f-c2a5-499e-b884-85fb59060ec5
md"""
# Inputs
"""

# ╔═╡ 03164025-b3e0-4535-abba-7cf7da5360d7
outpath = joinpath(ENV["DWARFS_ROOT"], "analysis/isolation", inputs.haloname)

# ╔═╡ c4fc0a87-c75f-4528-aa20-75cc7aff856c
FIGDIR = joinpath(outpath, "figures/")

# ╔═╡ 07ce4cef-df87-494b-8f16-02489e88bca6
using LilGuys; FIGDIR

# ╔═╡ a29c993a-c7eb-4b57-a474-50bdbd0ce1ec
halo = lguys.load_profile(joinpath(outpath, "halo.toml"))

# ╔═╡ 6435d239-4a83-4990-a25e-cd03bb0e2022
softening = 0.044 # kpc

# ╔═╡ 9104ed25-9bc8-4582-995b-37595b539281
out = lguys.Output(joinpath(outpath, "combined.hdf5"))

# ╔═╡ dd31d3ee-7fdf-46f5-b213-45faae93ae5e
idxs = [1, 3, 6, length(out)]

# ╔═╡ 38ab4838-960e-4321-b70c-6f14584e9e27
snaps = [out[i] for i in idxs]

# ╔═╡ 97f89831-00e6-49a2-a712-ac47fd2dee47
out.times[idxs] * lguys.T2GYR

# ╔═╡ 382d9256-5a83-4a70-be12-4818fe6c19b4
x_cen_err, v_cen_err = h5open(joinpath(outpath, "centres.hdf5"), "r") do f
	xe = f["position_errs"][:]
	ve = f["velocity_errs"][:]
	return xe, ve
end

# ╔═╡ 13fddaff-5379-47ba-a0cf-aba1c87124a0
scalars = read_fits(joinpath(outpath, "profiles_scalars.fits"))

# ╔═╡ 669967a7-0f29-4537-adc1-e96712793db7
profs = LilGuys.read_ordered_structs(joinpath(outpath, "profiles.hdf5"), LilGuys.MassProfile)

# ╔═╡ 537906ed-1c2d-4be6-9e5b-4696dd3aa854
density_profs = LilGuys.DensityProfile.(snaps)

# ╔═╡ 97e98ab8-b60b-4b48-b465-a34a16858f88
md"""
# Initial/Final
"""

# ╔═╡ 9a9f4dc1-3573-41ee-be1a-eee39d3371b0
fit = lguys.fit_v_r_circ_max(profs[end].second.radii, LilGuys.circular_velocity(profs[end].second))

# ╔═╡ 9dc4f57d-2bb2-4718-9f9d-86445e23bd75
 v_circ_max(halo) * V2KMS

# ╔═╡ 72e84e62-b6c6-4afe-9708-7be49064d081
r_circ_max(halo)

# ╔═╡ 0e89851e-763f-495b-b677-b664501a17ef
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"\log \; r / \textrm{kpc}", ylabel=L"V_\textrm{circ}")


	for i in eachindex(profs)
		lines!(log10.(profs[i].second.radii), LilGuys.circular_velocity(profs[i].second) * V2KMS, color=i, colorrange=(1, length(profs)))
	end


	scatter!(log10.(fit.r_circ_max), float.(fit.v_circ_max)* V2KMS)

	V_nfw(x) = lguys.v_circ(halo, x)
	log_r = LinRange(-2, 2.5, 1000)
	y = V_nfw.(10 .^ log_r)
	lines!(log_r, y * V2KMS, label="expected", linestyle=:dot, color=:black)

	axislegend(position=:lt)

	fig
end

# ╔═╡ 42a2814f-5a2b-4255-995a-b50472736e4f
mass(halo, 30)

# ╔═╡ 184c3f9b-d039-469c-b318-507c85b0a988
sum(snaps[1].masses)

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

# ╔═╡ f1717efb-2749-41dc-87f9-f265918a5984
density_profs

# ╔═╡ 264aedc1-b624-475e-af28-1d31b533839d
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log\, r / \textrm{kpc}", ylabel = L"\log\, \rho_\textrm{DM}\quad [10^{10} M_\odot / \textrm{kpc}^3]",
	limits=(-2.5, 4, -20, 1))

	for i in eachindex(density_profs)
		lines!(density_profs[i].log_r, log10.(density_profs[i].rho))
	end

	log_r = LinRange(-2, 4, 1000)
	y = log10.(lguys.density.(halo, 10 .^ log_r))
	lines!(log_r, y, label="analytic", color="black", linestyle=:dot)

	vlines!(log10(softening), color="grey")
	axislegend(ax)
	
	@savefig "dm_density"

	fig

end

# ╔═╡ 4e45e756-8a9c-43b4-aac7-2016347f5afb
let
	skip = 10
	snap = snaps[1]
	idx = 1:skip:length(snap)
	r = sort(lguys.radii(snap))[idx]
	
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
sum(lguys.radii(snaps[end]) .< 0.1)

# ╔═╡ 34d9fdea-8961-44ca-a92f-2f48a281f2cd
let
	fig, ax = FigAxis( ylabel="count", xlabel="log r / kpc", yscale=log10, yticks=Makie.automatic)
	
	hist!(log10.(lguys.radii(snaps[1])), label="", bins=100)

	fig
end

# ╔═╡ 5f1c61f9-50d4-43cb-aa78-fa85314f26b7
let
	fig, ax = FigAxis( ylabel="count", xlabel="e spec", yscale=log10, 
		limits=(-4, -1, 1, 1e6),
			yticks = Makie.automatic,
	)

	for i in eachindex(snaps)
		e = lguys.specific_energy(snaps[i])
		e = e[e .> 0]
		stephist!(log10.(e), bins=100, label="$i")
	end
	
	fig
end

# ╔═╡ 1ce0aa3c-953e-4f7b-bf7e-7c815c505b5e
println(sum(lguys.radii(snaps[1]) .< softening))

# ╔═╡ fb0dec74-aaab-43a4-9b37-d13634c5dcad
md"""
# Time variation
"""

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

	t = scalars.time
	Ls = hcat(scalars.L_1, scalars.L_2, scalars.L_3)'

	for i in 1:3
		scatter!(t, Ls[i, :], label=["x", "y", "z"][i])
		lines!(t, Ls[i, :])
	end

	axislegend(ax)
	fig
end

# ╔═╡ bfa7593c-4915-4e03-83e6-8f790de4c1a5
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="snapshot", ylabel="relative change in energy")
	Es = scalars.E
	scatter!((Es ./ Es[1]))

	fig
end

# ╔═╡ e3a45a8e-cc52-4e9d-9db3-97109b59fc77
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="snapshot", ylabel="-W / 2T (virial ratio)")
	
	E_kin = scalars.K
	E_pot = -scalars.W

	scatter!(E_pot ./ 2E_kin)
	hlines!(1, color=:black)


	@savefig "virial_ratio"
	fig
end

# ╔═╡ 91a44ed4-8466-4a58-b3ff-1e7630b8ac8c
let
	fig = lguys.plot_xyz(out.x_cen)
	fig.content[1].title = "centre"

	
	fig
end

# ╔═╡ 4d2da507-7250-4685-9602-6bd65247592a
cen_vel_mean = (out.x_cen[:, end]) ./ (out.times[end] - out.times[begin])

# ╔═╡ 7588c578-499b-4dae-817e-ecff230ab6be
let
	fig = lguys.plot_xyz((out.x_cen .- cen_vel_mean .* out.times') ./ x_cen_err')
	fig.content[1].title = "centre sigma"

	
	fig
end

# ╔═╡ faa5a990-71da-4971-afeb-b56ac04ad191
LilGuys.std((out.x_cen .- cen_vel_mean .* out.times') ./ x_cen_err')

# ╔═╡ bfe804ef-2080-431d-9544-ab151a3c609b
LilGuys.std((out.x_cen ) ./ x_cen_err')

# ╔═╡ e61c095e-a763-466b-b419-755fd0aadd0d
lguys.plot_xyz(out.v_cen * V2KMS, units=" / km s⁻¹")

# ╔═╡ fe25dc38-621b-47f8-aac5-f81c69ad7d8f
lguys.plot_xyz(out.v_cen ./ v_cen_err', units=" vel sigma")

# ╔═╡ 46bd928f-cbe3-4cb6-b09c-532164a64ce2
LilGuys.std((out.v_cen ) ./ v_cen_err')

# ╔═╡ 84f0bac9-5655-4acd-88db-8ba3114f712f
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel="log r containing fraction of DM")

	qs = profs[1].second.quantiles
	t = parse.(Int, first.(profs))

	for i in eachindex(qs)
		label = "$(qs[i]) "
		y = log10.([prof.second.r_quantile[i] for prof in profs])
		
		scatter!(t ./ 20, y, 				
			label=label, color=i, colorrange=(1, length(qs))
		)
	end



	Legend(fig[1,2], ax, "f = ")

	save(figure_dir * "mass_fraction_w_time.pdf", fig)

	fig
end

# ╔═╡ 967136d3-8d58-4fdc-9537-aa3a85a92528
times = out.times * T2GYR

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

# ╔═╡ ba823177-97e7-43df-8b0a-2ab1f2bafbc0
import LinearAlgebra: ⋅, ×

# ╔═╡ fc34af5e-069c-47c5-8599-a97d520c1426
NFW(v_circ_max = 0.1419, r_circ_max=4.2).r_s / 2.76 * 0.14

# ╔═╡ bb20e906-636b-4807-98fe-64453f746697
scalars.r_circ_max

# ╔═╡ 03da8ffc-247a-4580-8715-4d8e294b02a8
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = "log r", 
		ylabel = L"\beta",
		
	)

	idx = eachindex(out)[1:10:end]
	for i in idx
		snap = out[i]
		bins = lguys.quantile(radii(snap), LinRange(0, 1, 20+1))
		σ, β = LilGuys.β_prof(snap, r_bins=(bins))
		x = midpoints(log10.(bins))
		#vlines!(log10(5 * halo.r_s * 0.75))
		
		scatterlines!((x), β, color=i, colorrange=(1, length(out)))
	end
	
	fig
end

# ╔═╡ e5ca8db2-2c3d-4b97-9242-ab1d2ebf51b3
function calc_phase_dens(snap; bins=100)
	rs = lguys.radii(snap)
	vs = lguys.speeds(snap)
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
	
	h = Arya.hist2d!(ax, log10.(lguys.radii(snap_i.positions)), lguys.radii(snap_i.velocities) * V2KMS; hist_kwargs...)

	Colorbar(fig[1, 2], h)

	fig
end

# ╔═╡ 3f625247-0d7d-4f6d-bf84-a40f19f1849a
plot_rv_dens(snaps[end])

# ╔═╡ 9b88e37a-93eb-4c03-b654-0cd9284fc811
function plot_phase_rvr!(snap; bins=100)
	x = lguys.radii(snap)
	v = lguys.radial_velocities(snap)

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
# ╠═03a1bd03-dc40-42c0-bbcb-4aa562f47bbf
# ╟─4f1e87b3-f776-4550-a0f9-45eec1b79340
# ╠═7e548032-56a4-4b5c-bc47-ea86bb4cf917
# ╠═6e08e538-bc82-11ee-1a75-d97f506d18c5
# ╠═06d06ab0-4585-4cc6-b875-5cdd2e47657a
# ╠═07ce4cef-df87-494b-8f16-02489e88bca6
# ╠═374489bc-627f-4fc9-9734-7c49456710ac
# ╠═96c91860-f3cc-4531-a8cf-39c85887b394
# ╠═7e232e3c-34b9-4b08-a297-7aa19aca644b
# ╠═c4fc0a87-c75f-4528-aa20-75cc7aff856c
# ╟─7eb3e35f-c2a5-499e-b884-85fb59060ec5
# ╠═03164025-b3e0-4535-abba-7cf7da5360d7
# ╠═a29c993a-c7eb-4b57-a474-50bdbd0ce1ec
# ╠═dd31d3ee-7fdf-46f5-b213-45faae93ae5e
# ╠═38ab4838-960e-4321-b70c-6f14584e9e27
# ╠═6435d239-4a83-4990-a25e-cd03bb0e2022
# ╠═97f89831-00e6-49a2-a712-ac47fd2dee47
# ╠═9104ed25-9bc8-4582-995b-37595b539281
# ╠═382d9256-5a83-4a70-be12-4818fe6c19b4
# ╠═13fddaff-5379-47ba-a0cf-aba1c87124a0
# ╠═669967a7-0f29-4537-adc1-e96712793db7
# ╠═537906ed-1c2d-4be6-9e5b-4696dd3aa854
# ╟─97e98ab8-b60b-4b48-b465-a34a16858f88
# ╠═9a9f4dc1-3573-41ee-be1a-eee39d3371b0
# ╠═9dc4f57d-2bb2-4718-9f9d-86445e23bd75
# ╠═72e84e62-b6c6-4afe-9708-7be49064d081
# ╠═0e89851e-763f-495b-b677-b664501a17ef
# ╠═42a2814f-5a2b-4255-995a-b50472736e4f
# ╠═184c3f9b-d039-469c-b318-507c85b0a988
# ╠═7299dbaf-b332-4e53-85de-7acbfa0c3853
# ╟─a49d1735-203b-47dd-81e1-500ef42b054e
# ╠═06cafe6c-98ed-42ce-8b6c-b3e12afab896
# ╠═e5b9ce74-4d2d-4c5d-ad45-b6e233a4ec50
# ╠═f1717efb-2749-41dc-87f9-f265918a5984
# ╠═264aedc1-b624-475e-af28-1d31b533839d
# ╠═4e45e756-8a9c-43b4-aac7-2016347f5afb
# ╠═e33d56a7-7a0e-4fa9-8f0d-041b43584d59
# ╠═34d9fdea-8961-44ca-a92f-2f48a281f2cd
# ╠═5f1c61f9-50d4-43cb-aa78-fa85314f26b7
# ╠═1ce0aa3c-953e-4f7b-bf7e-7c815c505b5e
# ╟─fb0dec74-aaab-43a4-9b37-d13634c5dcad
# ╠═d6bc6ccb-15aa-415a-aa48-8a3cff89749e
# ╠═36b741d4-5a47-4d2e-8c93-0b29c53e2a0e
# ╠═17435494-a855-4295-9ce7-e60d937c2aa8
# ╠═bfa7593c-4915-4e03-83e6-8f790de4c1a5
# ╠═e3a45a8e-cc52-4e9d-9db3-97109b59fc77
# ╠═91a44ed4-8466-4a58-b3ff-1e7630b8ac8c
# ╠═4d2da507-7250-4685-9602-6bd65247592a
# ╠═7588c578-499b-4dae-817e-ecff230ab6be
# ╠═faa5a990-71da-4971-afeb-b56ac04ad191
# ╠═bfe804ef-2080-431d-9544-ab151a3c609b
# ╠═e61c095e-a763-466b-b419-755fd0aadd0d
# ╠═fe25dc38-621b-47f8-aac5-f81c69ad7d8f
# ╠═46bd928f-cbe3-4cb6-b09c-532164a64ce2
# ╠═84f0bac9-5655-4acd-88db-8ba3114f712f
# ╠═967136d3-8d58-4fdc-9537-aa3a85a92528
# ╟─5153654a-567f-4663-9b4a-6f08f9e49d1a
# ╟─a35b5f3d-ed9e-48f9-b96f-0a3c00ff2410
# ╠═b9746093-0f2f-4478-82ba-00911c8fcceb
# ╠═3cc0df1f-dca3-4355-a969-3874aa414d1a
# ╠═2dedc9e1-9b68-4c76-97da-4aba767de32d
# ╟─24c1b4c5-4be3-4ea0-8b0e-a0b6fb8647e9
# ╠═ba823177-97e7-43df-8b0a-2ab1f2bafbc0
# ╠═fc34af5e-069c-47c5-8599-a97d520c1426
# ╠═bb20e906-636b-4807-98fe-64453f746697
# ╠═03da8ffc-247a-4580-8715-4d8e294b02a8
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
