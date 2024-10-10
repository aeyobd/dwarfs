### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ e9d818c0-1897-11ef-1913-1b38c9d4b4d6
begin 
	import Pkg; Pkg.activate()
	
	using CairoMakie
	
	import LilGuys as lguys
	using Arya

	Arya.__init__()
end

# ╔═╡ 91564a7b-71f9-4619-9d51-69caf49af642
using DataFrames

# ╔═╡ d9ed11d2-0f91-4e27-8ff7-a546eb469719
using LilGuys

# ╔═╡ 9abf6c6e-20c4-48bc-b394-9ec58dc76c75
using KernelDensity

# ╔═╡ 0f7680b9-4529-4438-8a5e-a9d6e6530eaf
import TOML

# ╔═╡ b967b388-bb0a-4664-bef0-0236e37798f0
import StatsBase: sample, Weights

# ╔═╡ 576b1ea1-05bd-48ce-b977-276b9271fc1f
font = "/astro/dboyea/fonts/Arev.ttf"

# ╔═╡ b2f25a23-71bc-4217-8bef-c49d0cbc6676
begin 
	Arya.__init__()

	fg = "#000000ff"
	bg = "#ffffffff"
	
	update_theme!(Theme(
		fontsize=24,
		backgroundcolor=bg,
		pt_per_unit=0.352777, # units are mm
		px_per_unit=0.3527777, # px are 1 pt large
		textcolor=fg,
		Axis = (;
			xgridvisible=false,
			ygridvisible=false,
			backgroundcolor=:transparent,
			xmajortickwidth=20,
			splinescolor=fg,
			xtickcolor=fg,
			xminortickcolor=fg,
			ytickcolor=fg,
			yminortickcolor=fg,
			leftspinecolor=fg,
			rightspinecolor=fg,
			topspinecolor=fg,
			bottomspinecolor=fg
		),
		Legend = (;
			backgroundcolor=:transparent,
			framecolor=fg,
		),
		fonts = (;
			regular = font
		)
	)
	)
		
end

# ╔═╡ d794992b-328b-4327-96fa-67564dc8db98
Makie.load_font(font)

# ╔═╡ c64b9fe9-b3ff-498f-85dd-55e227443b63
name = "exp2d_rs0.10"

# ╔═╡ 96eafa68-e517-44ea-bf10-0649e82c8c9a
model_dir = ENV["DWARFS_ROOT"] * "/analysis/sculptor/1e6_V31_r3.2/orbit1/"

# ╔═╡ 9b6b3d13-d390-457c-b7fa-45f8fba3e4d0
readdir(model_dir)

# ╔═╡ be631a25-dc1a-4e28-9dfd-ed89c1de605d


# ╔═╡ 6d207751-695f-4247-94da-ced5a146092f
prof_expected = lguys.StellarProfile(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/processed/fiducial_sample_profile.toml"))

# ╔═╡ 88f31bfc-a67d-4654-af8b-46dc91500558
r_b = 58.69

# ╔═╡ de26b439-d9dd-449e-9289-9d3daed87cc7
fig_dir = "/astro/dboyea/figures"; mkpath(fig_dir)

# ╔═╡ 171b307c-cc13-4928-9526-a632031554e0
# ╠═╡ disabled = true
#=╠═╡
model_dir = "/astro/dboyea/sculptor/orbits/orbit1"
  ╠═╡ =#

# ╔═╡ c47ce6da-46b1-46e0-bd7a-acc5f6bdb92b
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ 570e2840-977b-46c6-b974-52b0ed7decf3
starsfile = "$model_dir/../stars/exp2d_rs0.10/probabilities_stars.hdf5"

# ╔═╡ 22323266-c04e-4865-bd1a-a699dfa1a7d3
begin 
	using HDF5

	f = h5open(starsfile)
	p_idx = f["index"][:]
	probabilities = f["probability"][:][sortperm(p_idx)]
	p_idx = sort(p_idx)
	close(f)
	
end

# ╔═╡ 85f4981d-82e0-45f7-ad9f-8bb17b5827b2
md"""
This notebook creates the plots I presented at the Durham Small Galaxies Cosmic Questions conference in 2024.
"""

# ╔═╡ 3df5d132-1a1a-40d6-b2d8-dcce2fcc4c34
out = lguys.Output("$model_dir", weights=probabilities)

# ╔═╡ a25f3edf-8d5b-401d-b727-8d8cd2d14910
snap_i = out[1]

# ╔═╡ 5c2534ca-ed11-47b4-9635-a132fd2a35e6
snap_f = out[orbit_props["idx_f"]]

# ╔═╡ a9e83ce1-86d2-4426-a12d-51dc85f2b152
params = TOML.parsefile(joinpath(model_dir, "../stars/$name/profile.toml"))

# ╔═╡ 41d07a5b-32a7-43f2-9a34-eab54c8ab4e0
R_s = params["Exp2D"]["R_s"]

# ╔═╡ 3b3a4fc6-5bf3-4b81-8621-8f425c7dd517
R_s_arcmin = lguys.kpc_to_arcmin(R_s, 82.3)

# ╔═╡ 937b4fd2-648d-4154-ae53-c06976d415b5
obs_df_all = lguys.read_fits("$model_dir/stars/$(name)/final.fits")

# ╔═╡ c9f0021e-ce16-4dcd-8480-c67c3c242430


# ╔═╡ f700b0bb-64b5-4be0-bf8a-88e72a7500e1
prof_f = lguys.StellarProfile("$model_dir/stars/$(name)/final_profile.toml")

# ╔═╡ 2e81adf7-47f4-4f61-857e-8997c15dc943
prof_i = lguys.StellarProfile("$model_dir/stars/$(name)/initial_profile.toml")

# ╔═╡ 2e69a41e-1b32-43e4-a23c-4550539275e6
prof_a_sky = lguys.Exp2D(; R_s=R_s_arcmin)

# ╔═╡ 10792def-f5cc-4963-9822-88ff8eabed95
function plot_model_and_exp!(prof, R_s_kpc; 
		y_offset=0, distance=distance, kwargs...)
	R_s_arcmin = lguys.kpc_to_arcmin(R_s_kpc, distance)
	println(R_s_arcmin)

	profile2 = lguys.Exp2D(R_s = R_s_arcmin)

	log_Σ(r) = log10(lguys.calc_Σ(profile2, r))

	log_R = LinRange(-2, maximum(prof.log_r), 1000)
	y = log_Σ.(10 .^ log_R)
	y .-= y[1]
	
	lines!(log_R, y .+ y_offset; kwargs...)

	errscatter!(prof.log_r, prof.log_Sigma .+ y_offset,
		yerr=prof.log_Sigma_err,
	)
end

# ╔═╡ 148c2f18-ddbb-4324-9fff-6619da29fe73
save = CairoMakie.save

# ╔═╡ de69d265-cda1-4543-ad02-2ee3091964d6
log_r_label = "log R / arcmin"

# ╔═╡ bfab4ae8-a94b-4a81-aad3-9b706f2474bb
function sigma_axis(; kwargs...) 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=log_r_label,
		ylabel = "log Σ / (fraction/arcmin^2)"
		;kwargs...
	)

	return fig, ax
end

# ╔═╡ 2d0a6d7d-c047-4688-ab4f-fbb489c7974f
fig_dir

# ╔═╡ 932c4fef-992b-4518-80d0-59c8e126ccb5
let 
	fig = Figure(
	)
	ax = Axis(fig[1,1], 
		xlabel=log_r_label,
		ylabel = "log Σ / Σ₀",
		limits=((0, 2.0), (-4.5, 0.5)),
	)

	errscatter!(prof_expected.log_r, prof_expected.log_Sigma .- 1.4,
		yerr=prof_expected.log_Sigma_err,
		label="Jensen+24",
		color=fg
	)


	ds = 2.4

	lines!(prof_i.log_r, prof_i.log_Sigma .+ ds, 
			label="initial", color=COLORS[3])
	
	lines!(prof_f.log_r, prof_f.log_Sigma .+ ds, 
			label="final", color=COLORS[2])

	
	vlines!(log10(r_b), color=:grey, label="break radius")
	# axislegend(position=:lb,
	# 	backgroundcolor=:transparent
	# )
	
	colsize!(fig.layout, 1, Aspect(1, 1.3))

	#resize_to_layout!(fig)

	save(joinpath(fig_dir, "density_i_f.svg"), fig)

	fig
end

# ╔═╡ 693cb85b-8f19-414f-8a26-da2035232df0
abspath(fig_dir)

# ╔═╡ 15cc735b-4da6-496e-acf6-3b090758935d
idx_f = orbit_props["idx_f"]

# ╔═╡ 165468ca-36e0-4024-97d5-14406022a4d8
let
	fig = Figure(
		figure_padding=0
	)

	obs_df = obs_df_all
	dec = obs_df.dec[1]
	ra = obs_df.ra[1]
	bw = 0.05

	dra = 10
	ddec = dra * cosd(dec)
	limits = ((ra - dra, ra + dra), 
	(dec - ddec, dec + ddec))
	
	ax = Axis(fig[1, 1], limits=limits,
		aspect=1,
		xreversed=true,
		xlabel="RA / degrees",
		ylabel="Dec / degrees",
	)
	
	x = obs_df.ra
	y = obs_df.dec
	hi = kde([x y], bandwidth=(bw, bw*cosd(dec)), boundary=limits, 
		weights=obs_df.weights, npoints=(8_000, 8_000))
	# areas = diff(hi.xbins) .* (diff(hi.ybins)')
	# hi.values ./= areas
	
		
	h = heatmap!(hi.x, hi.y, hi.density, colorscale=log10, colorrange=(1e-10*maximum(hi.density), maximum(hi.density)), colormap=(:greys))
	println(maximum(hi.density))
	
	
	# hidedecorations!(ax)
	# hidespines!(ax)

	fig
end

# ╔═╡ 6c94b8e5-3979-42e8-ae8f-61906c72c3db
function xi_eta_axis(dx=10, dy=5; kwargs...)
	fig = Figure(;kwargs...)
	

	limits = ((-dx, dx), (-dy, dy))

	ax = Axis(fig[1,1],
		xlabel="ξ / degrees",
		ylabel="η / degrees",
		limits=limits,
		aspect = DataAspect(),
		xgridvisible=false,
		ygridvisible=false
		
	)

	return fig, ax
end

# ╔═╡ 4a5adacf-db5e-4dc0-90c9-619e9b571ba8
function mean_2d(obs_df, values; bins=100, centre=false, val_min = 0, limits=nothing)
	if centre
		val_mean = lguys.mean(values, lguys.weights(obs_df.weights))
		println("centre: ", val_mean)
		val = values .- val_mean
	else
		val = values
	end
	
	x = obs_df.xi
	y = obs_df.eta
	weights = obs_df.weights

	h_vel = Arya.histogram2d(x, y, bins, weights=weights .* val, limits=limits)

	h_mass = Arya.histogram2d(x, y, bins, weights=weights, limits=limits)
	h_vel.values ./= h_mass.values
	h_vel.values[h_mass.values .< val_min] .= NaN

	return h_vel

end

# ╔═╡ ce3b2336-fce2-473f-bfe5-c0f07cfab6e0
df_j24 = lguys.read_fits("/astro/dboyea/dwarfs/observations/sculptor/processed/fiducial_sample.fits")

# ╔═╡ 01b0f4e8-9f0d-457b-900a-db4c85ec602d
let
	dr = 5
	fig, ax = xi_eta_axis(dr, dr)
	obs_df = obs_df_all

	bins = LinRange(-dr, dr, 100)
	dv = 15

	v0 = obs_df.radial_velocity[1]
	x = obs_df.xi
	y = obs_df.eta
	weights = obs_df.weights

	limits=((-dr, dr), (-dr, dr))
	val_min = 0
	val = obs_df.radial_velocity
	h_vel = Arya.histogram2d(x, y, bins, weights=weights .* val, limits=limits)

	h_mass = Arya.histogram2d(x, y, bins, weights=weights, limits=limits)
	h_vel.values ./= h_mass.values
	h_vel.values[h_mass.values .< val_min] .= NaN


	h = h_mass
	
	p = heatmap!(h.xbins, h.ybins, h.values, 
		colorscale=log10, colorrange=(1e-10 * maximum(h.values), maximum(h.values))
		#colormap=:redsblues,  colorrange=(v0-dv, v0+dv)
	)

	Colorbar(fig[1, 2], p, label="stellar density")

	#scatter!(df_j24.xi, df_j24.eta, color=:black, markersize=1)

	arc!(Point2f(0, 0), r_b/60, 0, 2π, color=:black)

	save("$fig_dir/skydensity.svg", fig)
	fig
end

# ╔═╡ 5e49981d-b456-4ae2-825d-bd30503de4eb
let
	ratio = (3600, 3600)
	fig = Figure(
		padding=0,
		backgroundcolor=:transparent

	)
	
	r_max = 160
	bandwidth=0.5

	dpcm = 5

	r_scale = 1/5
	limits = ((-ratio[1] * r_scale, ratio[1] * r_scale),
		(-ratio[2] * r_scale, ratio[2] * r_scale)
	)
	
	ax = Axis(fig[1,1], aspect=DataAspect(),
		xlabel = "y / kpc", ylabel="z / kpc", 
		#padding=0,
		backgroundcolor=:transparent,
		#title="dark matter", 
		)

	x, y = snap_f.positions[2, :], snap_f.positions[3, :]
	
	hi = kde([x y], bandwidth=(bandwidth, bandwidth), boundary=limits,
		npoints= dpcm .* ratio
	)

	z = hi.density
	z = log10.(1e-5  .+ hi.density ./ maximum(hi.density))
	
	image!(z,
		colormap=Reverse(:greys)
	)

	hidedecorations!(ax)
	hidespines!(ax)

	save(joinpath(fig_dir, "xy_projection.png"), fig)

	#colsize!(fig.layout, 1, Aspect(1, 1/sqrt(2)))

	resize_to_layout!(fig)
	fig
end

# ╔═╡ 5529178c-0e4e-4d10-ae8c-613b2fa2af3e
let
	ratio = (100, 100)
	fig = Figure(
		padding=0,
		backgroundcolor=:transparent
	)
	
	r_max = 160
	bandwidth=0.5

	dpcm = 5

	r_scale = 1/5
	limits = ((-ratio[1] * r_scale, ratio[1] * r_scale),
		(-ratio[2] * r_scale, ratio[2] * r_scale)
	)
	
	ax = Axis(fig[1,1], aspect=DataAspect(),
		xlabel = "y / kpc", ylabel="z / kpc", 
		backgroundcolor=:transparent,
		)

	x, y = out.x_cen[2, :], out.x_cen[3, :]
	
	lines!(x, y, color=:grey)

	save(joinpath(fig_dir, "xy_projection.png"), fig)

	#colsize!(fig.layout, 1, Aspect(1, 1/sqrt(2)))

	resize_to_layout!(fig)
	fig
end

# ╔═╡ 4cbf370e-b630-49ea-bc29-f253c1e42311
idx_samples = sample(1:size(obs_df_all, 1), Weights(obs_df_all.weights), 50_000)

# ╔═╡ ca16a975-35e8-46c0-acfb-5fe2ed40c16a
let

	s = 3
	fig = Figure(
		padding=0,
		backgroundcolor=:white,
		size=(300, 300) .* s
		
	)
	
	
	ax = Axis(fig[1,1], aspect=DataAspect(),
		backgroundcolor=:white,
		limits=(-15, 15, -15, 15) ./ 30,
		)


	df = obs_df_all[idx_samples, :]
	dx = 0.01
	N = size(df, 1)
	scatter!(df.xi .+ dx*randn(N), df.eta.+ dx*randn(N), 
		color=:black, alpha=0.05, markersize=3 * s, strokewidth=0
	)

	hidedecorations!(ax)

	save(joinpath(fig_dir, "stars_samples.png"), fig)
	fig
end

# ╔═╡ Cell order:
# ╟─85f4981d-82e0-45f7-ad9f-8bb17b5827b2
# ╠═e9d818c0-1897-11ef-1913-1b38c9d4b4d6
# ╠═0f7680b9-4529-4438-8a5e-a9d6e6530eaf
# ╠═91564a7b-71f9-4619-9d51-69caf49af642
# ╠═b967b388-bb0a-4664-bef0-0236e37798f0
# ╠═b2f25a23-71bc-4217-8bef-c49d0cbc6676
# ╠═576b1ea1-05bd-48ce-b977-276b9271fc1f
# ╠═d794992b-328b-4327-96fa-67564dc8db98
# ╠═c64b9fe9-b3ff-498f-85dd-55e227443b63
# ╠═96eafa68-e517-44ea-bf10-0649e82c8c9a
# ╠═9b6b3d13-d390-457c-b7fa-45f8fba3e4d0
# ╠═a25f3edf-8d5b-401d-b727-8d8cd2d14910
# ╠═be631a25-dc1a-4e28-9dfd-ed89c1de605d
# ╠═6d207751-695f-4247-94da-ced5a146092f
# ╠═88f31bfc-a67d-4654-af8b-46dc91500558
# ╠═de26b439-d9dd-449e-9289-9d3daed87cc7
# ╠═171b307c-cc13-4928-9526-a632031554e0
# ╠═3df5d132-1a1a-40d6-b2d8-dcce2fcc4c34
# ╠═c47ce6da-46b1-46e0-bd7a-acc5f6bdb92b
# ╠═5c2534ca-ed11-47b4-9635-a132fd2a35e6
# ╠═570e2840-977b-46c6-b974-52b0ed7decf3
# ╠═22323266-c04e-4865-bd1a-a699dfa1a7d3
# ╠═a9e83ce1-86d2-4426-a12d-51dc85f2b152
# ╠═41d07a5b-32a7-43f2-9a34-eab54c8ab4e0
# ╠═3b3a4fc6-5bf3-4b81-8621-8f425c7dd517
# ╠═937b4fd2-648d-4154-ae53-c06976d415b5
# ╠═c9f0021e-ce16-4dcd-8480-c67c3c242430
# ╠═f700b0bb-64b5-4be0-bf8a-88e72a7500e1
# ╠═2e81adf7-47f4-4f61-857e-8997c15dc943
# ╠═2e69a41e-1b32-43e4-a23c-4550539275e6
# ╠═10792def-f5cc-4963-9822-88ff8eabed95
# ╠═d9ed11d2-0f91-4e27-8ff7-a546eb469719
# ╠═148c2f18-ddbb-4324-9fff-6619da29fe73
# ╠═bfab4ae8-a94b-4a81-aad3-9b706f2474bb
# ╠═de69d265-cda1-4543-ad02-2ee3091964d6
# ╠═2d0a6d7d-c047-4688-ab4f-fbb489c7974f
# ╠═932c4fef-992b-4518-80d0-59c8e126ccb5
# ╠═693cb85b-8f19-414f-8a26-da2035232df0
# ╠═9abf6c6e-20c4-48bc-b394-9ec58dc76c75
# ╠═15cc735b-4da6-496e-acf6-3b090758935d
# ╠═165468ca-36e0-4024-97d5-14406022a4d8
# ╠═6c94b8e5-3979-42e8-ae8f-61906c72c3db
# ╠═4a5adacf-db5e-4dc0-90c9-619e9b571ba8
# ╠═ce3b2336-fce2-473f-bfe5-c0f07cfab6e0
# ╠═01b0f4e8-9f0d-457b-900a-db4c85ec602d
# ╠═5e49981d-b456-4ae2-825d-bd30503de4eb
# ╠═5529178c-0e4e-4d10-ae8c-613b2fa2af3e
# ╠═4cbf370e-b630-49ea-bc29-f253c1e42311
# ╠═ca16a975-35e8-46c0-acfb-5fe2ed40c16a
