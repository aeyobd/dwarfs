### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 2a0d82b2-3fe5-11ef-224e-e352922a2e5a
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames, CSV
	
	using CairoMakie

	
	import LilGuys as lguys
	using Arya
	
end

# ╔═╡ ae6cdb42-4193-436e-ab6a-008f9c6b7b32
using KernelDensity

# ╔═╡ d1667518-70c0-4235-bbbb-6e42d9ec7854
include("filter_utils.jl")

# ╔═╡ 218c549f-a8a9-4795-ac07-6e7aa3a5e758
filenames = ["centre", "leading", "trailing", "background"]

# ╔═╡ eb97c46b-5b9d-4ae6-82c1-2f7c409b69ab
import TOML

# ╔═╡ 13834fc5-04c5-44ed-ad34-a87408bfbe9b
scl_props = TOML.parsefile("/astro/dboyea/dwarfs/sculptor_obs_properties.toml")

# ╔═╡ 1aec8806-29bd-43e3-9b39-5719d9012c7c
function load_sample(filename)
	return lguys.load_fits("sculptor/$(filename)_sample.fits")
end

# ╔═╡ 72bd86f1-41ce-4727-979f-9262658f1873
load_sample("simple")

# ╔═╡ 5ec3651f-653d-4a43-bd09-da523c49f4cb
samples = load_sample.(filenames)

# ╔═╡ 67ed76ef-0da3-4c70-bf46-593b819f981e
ra_all = vcat([sample.ra for sample in samples]...)

# ╔═╡ 01f7fec7-ef67-45d7-a651-525f83eb0579
dec_all = vcat([sample.dec for sample in samples]...)

# ╔═╡ a8627570-50c3-48bd-a343-4ad6218ed064
allstars = unique(
	vcat(samples...), [:source_id], keep=:first)

# ╔═╡ b5593f3c-bc6d-4e8a-acc4-22eef851859f
unique(allstars.source_id)

# ╔═╡ 4683ce6c-4fce-4e05-afde-b50f2cbb3c9c
xi, eta = lguys.to_tangent(allstars.ra, allstars.dec, scl_props["ra"], scl_props["dec"])

# ╔═╡ ec2fea08-c8b6-4bbd-88d6-972a079d5cd8
obs_props = TOML.parsefile("/astro/dboyea/dwarfs/sculptor_obs_properties.toml")

# ╔═╡ 5eef7ffe-b473-4cc1-95d5-eeb887b97d8a
begin 
	orbit = lguys.load_fits("/astro/dboyea/sculptor/orbits/orbit1/skyorbit.fits")

	orbit[!, :xi], orbit[!, :eta] = lguys.to_tangent(orbit.ra, orbit.dec, obs_props["ra"], obs_props["dec"])
end

# ╔═╡ f328d46f-8ede-4c91-97f5-724dac51064d
bandwidth = (0.05, 0.05)

# ╔═╡ 7e9fd578-2789-428f-acb9-aa700b8569b4
kd = kde([xi eta], bandwidth=bandwidth)

# ╔═╡ 680c83c4-ebae-446c-8171-224b180370cc
idx_orbit=1775:1800

# ╔═╡ 146894f9-0da8-4be5-90be-677566abdfcf
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xreversed = true,
		xlabel="ra / degrees",
		ylabel="dec / degrees"
	)


	scatter!(ra_all, dec_all, alpha=0.05)

	lines!(orbit.ra[idx_orbit], orbit.dec[idx_orbit], color=COLORS[2])
	fig
end

# ╔═╡ 836befe3-9f91-49d0-b607-41506586bbe4
function xieta_axis()
	fig = Figure()
	ax = Axis(fig[1, 1], aspect=DataAspect(),
		xlabel = L"$\xi$ / degrees",
		ylabel = L"$\eta$ / degrees",
		xreversed = true,
		#limits=(-r_ell_max, r_ell_max, -r_ell_max, r_ell_max) ./ 60
	)

	return fig, ax
end

# ╔═╡ d7edb06b-e54a-46cc-9b74-0914931ac507
let
	fig = Figure()
	ax = Axis(fig[1, 1], aspect=DataAspect(),
		xlabel = L"$\xi$ / degrees",
		ylabel = L"$\eta$ / degrees",
		xreversed = true,
		#limits=(-r_ell_max, r_ell_max, -r_ell_max, r_ell_max) ./ 60
	)
	scale = 0.008


	h = heatmap!(kd.x, kd.y, asinh.(kd.density ./ scale))

	lines!(orbit.xi[idx_orbit], orbit.eta[idx_orbit])
	Colorbar(fig[1, 2], h, label="asinh density / $scale")

	resize_to_layout!(fig)
	fig
end

# ╔═╡ 8fdb7e29-9f99-40ca-a0fb-13c51da49977
md"""
# CMD along the orbital trajectory
"""

# ╔═╡ 1148c29c-b4ed-4bad-9077-0eaaf6ff9e1f
md"""
As a second test to if there are any signs of tidal tails, we can construct the CMD along the orbital trajectories for stars matching Sculptors proper motions.
"""

# ╔═╡ b546fc5c-3ca6-4c2f-9946-1777aecaacf5
allstars_cen = samples[1]

# ╔═╡ e6528acc-a86b-412a-a668-7886b324a9ea
let
	fig, ax = xieta_axis()
	
	scatter!(allstars_cen.xi, allstars_cen.eta, markersize=2)

	fig
end

# ╔═╡ e5721c60-9037-4203-9dee-6d54d27113fc
let
	r_ell_max = 240
	fig = Figure()
	ax = Axis(fig[1, 1], aspect=DataAspect(),
		xlabel = L"$\xi$ / degrees",
		ylabel = L"$\eta$ / degrees",
		xreversed = true,
		limits=(-r_ell_max, r_ell_max, -r_ell_max, r_ell_max) ./ 60
	)
	scale = 0.008

	kd = kde([allstars_cen.xi allstars_cen.eta], bandwidth=bandwidth)
	
	h = heatmap!(kd.x, kd.y, log10.(kd.density), colorrange=(-2.5, 1))

	lines!(orbit.xi[idx_orbit], orbit.eta[idx_orbit])
	
	Colorbar(fig[1, 2], h, label="log density")

	resize_to_layout!(fig)
	fig
end

# ╔═╡ b28bd5ab-c023-40b9-9c96-4f00bfeeaabb
begin 
	pm_only_params = read_file("sculptor/simple.toml")
	#pop!(pm_only_params["cmd_cut"])
	pm_only_params = DensityParams(pm_only_params)
end

# ╔═╡ 3cf5efd9-1c41-416e-b133-ee488ee02904
readdir("sculptor")

# ╔═╡ 5a51d8ff-ff1c-4995-add6-a7ad4fad5db4
allstars_pm = load_stars("sculptor/" * pm_only_params.filename, pm_only_params)

# ╔═╡ 067d17d3-4b0b-4868-bf08-b3a139487189


# ╔═╡ d1925d0c-030a-4d4f-b932-74185c870aec
import LinearAlgebra: normalize, ×

# ╔═╡ 1afeb2ce-3da0-44ec-9840-9351e61e7758
orbit_vector = normalize([
	orbit.xi[idx_orbit][end] - orbit.xi[idx_orbit][1],
	orbit.eta[idx_orbit][end] - orbit.eta[idx_orbit][1]
])

# ╔═╡ 767fc400-6d9d-49a4-97ed-f492036e9a94
function select_region(allstars, centre; radius=0.5)
	x_cen, y_cen = centre

	r = @. sqrt((allstars.xi - x_cen)^2 + (allstars.eta - y_cen)^2)

	filt = r .< radius

	return allstars[filt, :]
end

# ╔═╡ 249b8e2d-0401-4102-9b97-4701403f6569
function plot_region!(df; kwargs...)
	scatter!(df.xi, df.eta; kwargs...)
end

# ╔═╡ 37a056f7-0294-46a7-992a-d7ad4a210c52
function cmd_axis(gs)
	return Axis(gs,
		xlabel = "Bp-Rp",
		ylabel = "G",
		yreversed=true,
		limits = (-0.5, 2, 16, 22)
	)
end

# ╔═╡ c63495d0-a748-44b4-a9a8-6d2e5036a776
function plot_cmd!(df; kwargs...)
	scatter!(df.bp_rp, df.phot_g_mean_mag, alpha=0.3; kwargs...)
end

# ╔═╡ 54093519-551f-45af-bc32-9ae374f41da4
function plot_cmd(df; kwargs...)

	fig = Figure()
	ax = cmd_axis(fig[1, 1])

	plot_cmd!(df; kwargs...)
	return fig, ax
end

# ╔═╡ 0d93e04f-41e3-4934-84af-aaceddfc0ac3
bg_vector = ([orbit_vector; 0] × [0, 0, 1])[1:2]

# ╔═╡ 881194b9-394d-4f71-8557-b8b962ad8261
rs_test = [0.5, 1, 1.25, 1.5, 3.5]

# ╔═╡ 50b5fead-5387-4ee2-a0a5-03c305ba14d6
max_r_test = 0.25

# ╔═╡ 79a7a2a7-6d4d-457e-b857-4684b9381f12
function plot_cmd_members!(centre, radius)
	df = select_region(allstars_pm, centre, radius=radius)

	plot_cmd!(df, color=:grey)

	df = select_region(allstars, centre, radius=radius)
	
	plot_cmd!(df, color=COLORS[2])

	N = size(df, 1)
	text!(0.05, 0.9, text="$N stars in CMD", color=:black, space=:relative)
end

# ╔═╡ 6033a274-b0be-45c6-8fdd-50852c770252
let
	fig = Figure()
	ax = cmd_axis(fig[1, 1])
	
	plot_cmd_members!([0, 0], max_r_test)[1]
	fig
end

# ╔═╡ 11109ee7-bfaa-4359-8204-3dc9e329c76e
let
	rs = rs_test

	N = length(rs)

	fig = Figure(size=(600, 1200))
	
	for i in 1:N
		r = rs[i]
		ax = cmd_axis(fig[i, 1])
		plot_cmd_members!(orbit_vector * r, max_r_test)
		
		text!(0.05, 0.5, text="$r degrees along orbit", space=:relative)
		if i < N
			hidexdecorations!(ax, grid=false)
		end
	end

	linkxaxes!(fig.content...)
	fig
end

# ╔═╡ 40eab3d3-df4c-41bb-994c-ebc7ba210fdd
let
	rs = -rs_test

	N = length(rs)

	fig = Figure(size=(600, 1200))
	
	for i in 1:N
		r = rs[i]
		ax = cmd_axis(fig[i, 1])
		plot_cmd_members!(orbit_vector * r, max_r_test)
		
		text!(0.05, 0.5, text="$r degrees along orbit", space=:relative)
		if i < N
			hidexdecorations!(ax, grid=false)
		end
	end

	linkxaxes!(fig.content...)
	fig
end

# ╔═╡ 636ac8df-0ec6-4a6e-b484-b68e13fb837b
let
	rs = rs_test
	
	N = length(rs)

	fig = Figure(size=(600, 1200))
	
	for i in 1:N
		r = rs[i]
		ax = cmd_axis(fig[i, 1])
		plot_cmd_members!(bg_vector * r,max_r_test)
		
		text!(0.05, 0.5, text="$r degrees ⟂ to orbit", space=:relative)
		if i < N
			hidexdecorations!(ax, grid=false)
		end
	end

	linkxaxes!(fig.content...)
	fig
end

# ╔═╡ 3630d70d-4e61-4ef5-9fb2-3f513201f07a
let
	rs = -rs_test

	N = length(rs)

	fig = Figure(size=(600, 1200))
	
	for i in 1:N
		r = rs[i]
		ax = cmd_axis(fig[i, 1])
		plot_cmd_members!(bg_vector * r, max_r_test)
		
		text!(0.05, 0.5, text="$r degrees ⟂ to orbit", space=:relative)
		if i < N
			hidexdecorations!(ax, grid=false)
		end
	end

	linkxaxes!(fig.content...)
	fig
end

# ╔═╡ dd277c41-cb0a-4d4a-866d-92af40e1f775
let
	fig, ax = xieta_axis()

	rs = rs_test

	println("background, +")
	for i in eachindex(rs)
		r = rs[i]
		df = select_region(allstars, r * bg_vector, radius=max_r_test)
		plot_region!(df, color=r, colorrange=extrema(rs), colormap=:blues, label="⟂, +; r=$r")

		println(size(df, 1))
	end

	println("background, -")
	for i in eachindex(rs)
		r = -rs[i]
		df = select_region(allstars, r * bg_vector, radius=max_r_test)
		plot_region!(df, color=r, 
			colorrange=extrema(-rs), colormap=Reverse(:greens), label="⟂, -; r=$r")

		println(size(df, 1))
	end

	println("along orbit, +")
	for i in eachindex(rs)
		r = rs[i]
		df = select_region(allstars, r * orbit_vector, radius=max_r_test)
		plot_region!(df, color=r, 
			colorrange=extrema(rs), colormap=:reds, label="along orbit, +; r=$r")
		println(size(df, 1))
	end

	println("along orbit, -")
	for i in eachindex(rs)
		r = -rs[i]
		df = select_region(allstars, r * orbit_vector, radius=max_r_test)
		plot_region!(df, color=r, colorrange=(-rs[end], 0), colormap=(:greys),
		label = "along orbit, -; r=$r")
		println(size(df, 1))
	end
	
	lines!(orbit.xi[idx_orbit], orbit.eta[idx_orbit])

	#Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 428b5220-3fb6-423e-ab88-a9337e52a9c1


# ╔═╡ Cell order:
# ╠═2a0d82b2-3fe5-11ef-224e-e352922a2e5a
# ╠═218c549f-a8a9-4795-ac07-6e7aa3a5e758
# ╠═eb97c46b-5b9d-4ae6-82c1-2f7c409b69ab
# ╠═13834fc5-04c5-44ed-ad34-a87408bfbe9b
# ╠═1aec8806-29bd-43e3-9b39-5719d9012c7c
# ╠═72bd86f1-41ce-4727-979f-9262658f1873
# ╠═5ec3651f-653d-4a43-bd09-da523c49f4cb
# ╠═67ed76ef-0da3-4c70-bf46-593b819f981e
# ╠═01f7fec7-ef67-45d7-a651-525f83eb0579
# ╠═a8627570-50c3-48bd-a343-4ad6218ed064
# ╠═ae6cdb42-4193-436e-ab6a-008f9c6b7b32
# ╠═b5593f3c-bc6d-4e8a-acc4-22eef851859f
# ╠═4683ce6c-4fce-4e05-afde-b50f2cbb3c9c
# ╠═7e9fd578-2789-428f-acb9-aa700b8569b4
# ╠═ec2fea08-c8b6-4bbd-88d6-972a079d5cd8
# ╠═5eef7ffe-b473-4cc1-95d5-eeb887b97d8a
# ╠═f328d46f-8ede-4c91-97f5-724dac51064d
# ╠═680c83c4-ebae-446c-8171-224b180370cc
# ╠═146894f9-0da8-4be5-90be-677566abdfcf
# ╠═836befe3-9f91-49d0-b607-41506586bbe4
# ╠═d7edb06b-e54a-46cc-9b74-0914931ac507
# ╟─8fdb7e29-9f99-40ca-a0fb-13c51da49977
# ╟─1148c29c-b4ed-4bad-9077-0eaaf6ff9e1f
# ╠═b546fc5c-3ca6-4c2f-9946-1777aecaacf5
# ╠═e6528acc-a86b-412a-a668-7886b324a9ea
# ╠═e5721c60-9037-4203-9dee-6d54d27113fc
# ╠═d1667518-70c0-4235-bbbb-6e42d9ec7854
# ╠═b28bd5ab-c023-40b9-9c96-4f00bfeeaabb
# ╠═3cf5efd9-1c41-416e-b133-ee488ee02904
# ╠═5a51d8ff-ff1c-4995-add6-a7ad4fad5db4
# ╠═067d17d3-4b0b-4868-bf08-b3a139487189
# ╠═d1925d0c-030a-4d4f-b932-74185c870aec
# ╠═1afeb2ce-3da0-44ec-9840-9351e61e7758
# ╠═767fc400-6d9d-49a4-97ed-f492036e9a94
# ╠═249b8e2d-0401-4102-9b97-4701403f6569
# ╠═37a056f7-0294-46a7-992a-d7ad4a210c52
# ╠═54093519-551f-45af-bc32-9ae374f41da4
# ╠═c63495d0-a748-44b4-a9a8-6d2e5036a776
# ╠═0d93e04f-41e3-4934-84af-aaceddfc0ac3
# ╠═881194b9-394d-4f71-8557-b8b962ad8261
# ╠═50b5fead-5387-4ee2-a0a5-03c305ba14d6
# ╠═79a7a2a7-6d4d-457e-b857-4684b9381f12
# ╠═6033a274-b0be-45c6-8fdd-50852c770252
# ╠═11109ee7-bfaa-4359-8204-3dc9e329c76e
# ╠═40eab3d3-df4c-41bb-994c-ebc7ba210fdd
# ╠═636ac8df-0ec6-4a6e-b484-b68e13fb837b
# ╠═3630d70d-4e61-4ef5-9fb2-3f513201f07a
# ╠═dd277c41-cb0a-4d4a-866d-92af40e1f775
# ╠═428b5220-3fb6-423e-ab88-a9337e52a9c1
