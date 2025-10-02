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

# ╔═╡ 340ffbbe-17bd-11ef-35c6-63505bb128b7
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	import LilGuys as lguys

	using Arya
end

# ╔═╡ a8691ea6-a672-4aa9-b058-9c6abc35bf31
using OrderedCollections

# ╔═╡ d401ec4b-048e-4aae-85a8-f7f0d8e44a79
using LilGuys

# ╔═╡ 0e405727-b322-4cf5-98f9-bb3e695f616d
using PyFITS

# ╔═╡ 283d335c-2941-4dc6-9a8e-c47b0864087c
using PlutoUI

# ╔═╡ 377284f2-dcee-44d3-9a04-728605cea92a
md"""
Given the stellar probabilty file, makes plots based on the 3D properties of the stars in the sample.
"""

# ╔═╡ 30ba4458-f21b-4777-987c-e65ecfd34258
md"""
# Setup
"""

# ╔═╡ faeaf38d-8c06-4646-8179-57ffb05f720e
import DensityEstimators as DE

# ╔═╡ f0d2b68a-fae2-4486-a434-a8816e400e84
import TOML

# ╔═╡ 5c7e710a-c7ed-452a-891b-eae8bf9e6d49
CairoMakie.activate!(type=:png)

# ╔═╡ da3a9778-0f7c-482b-831b-d93a0d80a48c
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

# ╔═╡ ab06c999-3ff6-4580-a979-f0ddeb466569
@bind inputs confirm(notebook_inputs(;
	galaxyname = TextField(default="sculptor"),
	modelname = TextField(60, default="1e6_v38_r7_oblate_0.5/orbit_"),
	starsname = TextField(default="exp2d_rs0.13"),
))

# ╔═╡ 2106bfe1-a53d-4ef8-9111-e191a8056351
starsname = inputs.starsname

# ╔═╡ f0d74eaa-81e9-4b04-9765-24a0935b1430
model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", inputs.galaxyname, inputs.modelname)

# ╔═╡ 9c42eb0a-029d-46f7-afb0-de03f82c5889
obs_today_filename =  ENV["DWARFS_ROOT"] * "/observations/$(inputs.galaxyname)/observed_properties.toml"

# ╔═╡ 08aa0f76-3d74-45b5-b9e9-6abbf6350910
stars_dir_in = joinpath(model_dir, "../stars/$starsname")

# ╔═╡ 7c64a7c7-232b-47b6-973a-62e1ac21df3a
stars_dir_out = joinpath(model_dir, "stars/$starsname")

# ╔═╡ 64350409-6bae-4e1f-be11-b2ec7d48d1f1
FIGDIR = joinpath(stars_dir_out,  "figures")

# ╔═╡ 84755c6a-8a91-4d9a-a6ec-895debcf6608
readdir(stars_dir_out)

# ╔═╡ 1b5c00d2-9df6-4a9c-ae32-05abcbf0e41a
paramsfile = joinpath(stars_dir_in, "profile.toml")

# ╔═╡ 3acc8caf-5a7e-4a6d-84b2-b26680ffe377
md"""
initial & final snapshots
"""

# ╔═╡ 396cd0a8-1d73-44dd-89db-3243fb9e8ac4
md"""
# File loading
"""

# ╔═╡ 436a5be3-f597-4fc4-80a8-dc5af302ad66
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ d86a8f96-a6ee-4fdc-a8d7-322085c06e78
orbit_props["idx_peris"]

# ╔═╡ d374b131-6333-4bea-bcbf-f51e26e2438b
orbit_props["idx_apos"]

# ╔═╡ 84dc77f7-14b3-4a2e-a556-c025d7df0095
params = TOML.parsefile(paramsfile)

# ╔═╡ 2612e3a2-6e6e-494e-b140-720dd2db6ec2
obs_expected = TOML.parsefile(obs_today_filename)

# ╔═╡ f9fe37ef-de81-4d69-9308-cda968851ed2
df_probs = lguys.read_hdf5_table(stars_dir_in * "/probabilities_stars.hdf5")

# ╔═╡ e76a6fa2-c020-48bf-b065-bc7ca51ecd98
p_idx = df_probs.index; probabilities = df_probs.probability

# ╔═╡ 07defcc4-8dca-4caf-a34a-a341818d80ad
length(probabilities)

# ╔═╡ 6fb4b7e1-a22c-4ff8-bbe9-dbf5de5acd37
out =  lguys.Output(model_dir, weights=probabilities)

# ╔═╡ 93694113-9b3d-4a11-a34b-2a73350d2d87
@bind idxs_in confirm(RangeSlider(1:length(out), default=1:orbit_props["idx_f"]))

# ╔═╡ 973955ad-3214-42cf-831f-a782f1d2434a
idx_i = idxs_in[1]

# ╔═╡ 8f9ee2a7-de43-43e0-8257-93ecc630044f
idx_f = idxs_in[end]

# ╔═╡ 396a53a3-de0f-4d97-9693-40f3757d66f9
snap_i = out[idx_i]

# ╔═╡ 44cbb2ce-5f43-4bc9-a4c4-c0b9df692cd2
length(snap_i.masses)

# ╔═╡ 6feeaae2-cb01-46ad-ad1d-daaca1caf7ec
snap_f = out[idx_f]

# ╔═╡ 6f7c27ce-1328-45e5-b42c-768b76c47607
t_f = snap_f.time

# ╔═╡ 1d96e0ea-17b4-4964-9ca3-9b0f8195ee3f
@assert t_f == out.times[idx_f]

# ╔═╡ 5ee4f95d-0587-44ab-b543-9b7039d545e6
md"""
# Plots
"""

# ╔═╡ 77479cd4-513c-4603-9aa0-1acd964c403a
let
	fig = Figure(size=(700, 700))
	r_max = 100

	x = lguys.x_position(snap_f)
	y = lguys.y_position(snap_f)
	z = lguys.z_position(snap_f)
	ps = snap_f.weights
	bins = LinRange(-r_max, r_max, 300)

	h = Arya.histogram2d(y, z, bins, weights=ps)
	cmax = maximum(h.values)
	kwargs = (colorscale=log10, colorrange=(1e-6*cmax, cmax), weights=ps, bins=bins)

	
	ax_yz = Axis(fig[2,2], aspect=1,
		xlabel = "y / kpc", ylabel="z / kpc",
	)
	hm = Arya.hist2d!(ax_yz, y, z; kwargs...)
	hideydecorations!(ax_yz)
		

	ax_xz = Axis(fig[2,1], aspect=1,
		xlabel = "x / kpc", ylabel="z / kpc",
	)
	Arya.hist2d!(ax_xz, x, z; kwargs...)
	


	ax_xy = Axis(fig[1,1], aspect=1,
		xlabel = "x / kpc", ylabel="y / kpc", 
	)
	hidexdecorations!(ax_xy)

	Arya.hist2d!(ax_xy, x, y; kwargs...)


	linkyaxes!(ax_xz, ax_yz)
	linkxaxes!(ax_xz, ax_xz)

	Colorbar(fig[1, 2], hm, tellwidth=false, label="stellar density")

    rowsize!(fig.layout, 1, ax_yz.scene.viewport[].widths[1])
	rowgap!(fig.layout, 1, -50.0)
	colgap!(fig.layout, 1, 50.)

	resize_to_layout!(fig)
	
	@savefig "density_xy_zoomout"

	
	fig
end

# ╔═╡ 9e36808a-43b3-4a73-aa5e-36d05ae8bfcc
import DensityEstimators: histogram2d

# ╔═╡ b1817495-c364-4a75-aa12-150c5afe4323
function plot_isocontours(fig, snap; filter_bound=false, bins = -2:0.1:2)
	
	for idx in 1:3
		i, j = [(1,2), (2,3), (1,3)][idx]
		if filter_bound
			filt = LilGuys.bound_particles(snap)
		else
			filt = trues(length(snap))
		end
		x = snap.positions[i, filt] .- snap.x_cen[i]
		y = snap.positions[j, filt] .- snap.x_cen[j]

		h = histogram2d(x, y, (bins, bins), weights=snap.weights)
		h_min = 1e-10

		h_max = log10(maximum(h.values))
		
		ax = Axis(fig[1,idx])
		contour!(midpoints(h.xbins), midpoints(h.ybins), log10.(h.values .+ h_min), levels=LinRange(h_max-5, h_max, 10), linewidth=theme(:linewidth)[]/2, linestyle=:solid)

		if idx == 1
			text!(0.1, 0.9, space=:relative, text="$(round(snap.time * T2GYR, digits=2)) Gyr", align=(:left, :top))
		end
		hidedecorations!(ax)
	end

	#rowsize!(fig.layout, 1, Aspect(1, 1.0))
	fig
end

# ╔═╡ f743edac-e31a-452b-9243-31abbc64a506
let 
	fig = Figure(figsize=(2*72, 10*72))
	
	for (i, idx) in enumerate([1; orbit_props["idx_apos"]])
		gs = GridLayout(fig[i, 1])
		
		plot_isocontours(gs, out[idx])
		rowsize!(fig.layout, i, Aspect(1,1/3))
		colgap!(gs, 0)
	end

	resize_to_layout!(fig)
	fig
end

# ╔═╡ 865b194d-9efa-41c2-8184-63bcfcd90f6f
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc", title="initial")

	bin_range = LinRange(-10, 10, 100)
	bins = (out.x_cen[1, idx_i]  .+ bin_range,  out.x_cen[2, idx_i]  .+ bin_range)

	h = Arya.histogram2d( snap_i.positions[1, :], snap_i.positions[2, :], bins, weights=snap_i.weights)
	
	colorrange =(1e-7, 1) .* maximum(h.values)

		

	Arya.hist2d!(ax, snap_i.positions[1, :], snap_i.positions[2, :], weights=snap_i.weights, bins = bins, colorscale=log10, colorrange=colorrange)

	ax2 = Axis(fig[1,2], aspect=1,
	xlabel = "x / kpc", 
	title="final")


	bins = (out.x_cen[1, idx_f]  .+ bin_range,  out.x_cen[2, idx_f]  .+ bin_range)

	h = Arya.hist2d!(ax2, snap_f.positions[1, :], snap_f.positions[2, :], weights=snap_f.weights, bins = bins, colorscale=log10, colorrange=colorrange)

	Colorbar(fig[1, 3], h, label="stellar density")

    rowsize!(fig.layout, 1, ax.scene.viewport[].widths[1])

	fig
end

# ╔═╡ 87c9a6dc-cb71-4476-b7d7-2fa9b5cd2d10
filt_inner = radii(snap_f) .< 1

# ╔═╡ 9fc833c7-da43-4f97-a106-a8e6aa3c383f
[LilGuys.std(snap_f.positions[i, filt_inner], snap_f.weights[filt_inner]) for i in 1:3]

# ╔═╡ 702c40e9-f82f-4ffd-8a1e-35b129eb64f3
[LilGuys.std(snap_f.velocities[i, filt_inner], snap_f.weights[filt_inner]) for i in 1:3] .* V2KMS

# ╔═╡ Cell order:
# ╟─377284f2-dcee-44d3-9a04-728605cea92a
# ╠═ab06c999-3ff6-4580-a979-f0ddeb466569
# ╟─30ba4458-f21b-4777-987c-e65ecfd34258
# ╠═a8691ea6-a672-4aa9-b058-9c6abc35bf31
# ╠═340ffbbe-17bd-11ef-35c6-63505bb128b7
# ╠═faeaf38d-8c06-4646-8179-57ffb05f720e
# ╠═d401ec4b-048e-4aae-85a8-f7f0d8e44a79
# ╠═0e405727-b322-4cf5-98f9-bb3e695f616d
# ╠═f0d2b68a-fae2-4486-a434-a8816e400e84
# ╠═283d335c-2941-4dc6-9a8e-c47b0864087c
# ╠═5c7e710a-c7ed-452a-891b-eae8bf9e6d49
# ╠═da3a9778-0f7c-482b-831b-d93a0d80a48c
# ╠═64350409-6bae-4e1f-be11-b2ec7d48d1f1
# ╠═84755c6a-8a91-4d9a-a6ec-895debcf6608
# ╠═2106bfe1-a53d-4ef8-9111-e191a8056351
# ╠═f0d74eaa-81e9-4b04-9765-24a0935b1430
# ╠═9c42eb0a-029d-46f7-afb0-de03f82c5889
# ╠═08aa0f76-3d74-45b5-b9e9-6abbf6350910
# ╠═7c64a7c7-232b-47b6-973a-62e1ac21df3a
# ╠═1b5c00d2-9df6-4a9c-ae32-05abcbf0e41a
# ╟─3acc8caf-5a7e-4a6d-84b2-b26680ffe377
# ╠═93694113-9b3d-4a11-a34b-2a73350d2d87
# ╠═d86a8f96-a6ee-4fdc-a8d7-322085c06e78
# ╠═d374b131-6333-4bea-bcbf-f51e26e2438b
# ╠═973955ad-3214-42cf-831f-a782f1d2434a
# ╠═8f9ee2a7-de43-43e0-8257-93ecc630044f
# ╟─396cd0a8-1d73-44dd-89db-3243fb9e8ac4
# ╠═436a5be3-f597-4fc4-80a8-dc5af302ad66
# ╠═84dc77f7-14b3-4a2e-a556-c025d7df0095
# ╠═2612e3a2-6e6e-494e-b140-720dd2db6ec2
# ╠═f9fe37ef-de81-4d69-9308-cda968851ed2
# ╠═e76a6fa2-c020-48bf-b065-bc7ca51ecd98
# ╠═07defcc4-8dca-4caf-a34a-a341818d80ad
# ╠═44cbb2ce-5f43-4bc9-a4c4-c0b9df692cd2
# ╠═6fb4b7e1-a22c-4ff8-bbe9-dbf5de5acd37
# ╠═396a53a3-de0f-4d97-9693-40f3757d66f9
# ╠═6feeaae2-cb01-46ad-ad1d-daaca1caf7ec
# ╠═6f7c27ce-1328-45e5-b42c-768b76c47607
# ╠═1d96e0ea-17b4-4964-9ca3-9b0f8195ee3f
# ╟─5ee4f95d-0587-44ab-b543-9b7039d545e6
# ╠═77479cd4-513c-4603-9aa0-1acd964c403a
# ╠═b1817495-c364-4a75-aa12-150c5afe4323
# ╠═9e36808a-43b3-4a73-aa5e-36d05ae8bfcc
# ╠═f743edac-e31a-452b-9243-31abbc64a506
# ╠═865b194d-9efa-41c2-8184-63bcfcd90f6f
# ╠═87c9a6dc-cb71-4476-b7d7-2fa9b5cd2d10
# ╠═9fc833c7-da43-4f97-a106-a8e6aa3c383f
# ╠═702c40e9-f82f-4ffd-8a1e-35b129eb64f3
