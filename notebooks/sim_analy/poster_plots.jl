### A Pluto.jl notebook ###
# v0.19.45

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

# ╔═╡ 576b1ea1-05bd-48ce-b977-276b9271fc1f
font = "/astro/dboyea/fonts/Arev.ttf"

# ╔═╡ b2f25a23-71bc-4217-8bef-c49d0cbc6676
begin 
	Arya.__init__()

	fg = "#90cacfff"
	bg = "#041729ff"
	
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
name = "exp2d_rs0.1"

# ╔═╡ 6d207751-695f-4247-94da-ced5a146092f
prof_expected = lguys.ObsProfile("/arc7/home/dboyea/dwarfs/notebooks/density_fits/sculptor/fiducial_profile.toml")

# ╔═╡ 88f31bfc-a67d-4654-af8b-46dc91500558
r_b = 58.69

# ╔═╡ de26b439-d9dd-449e-9289-9d3daed87cc7
fig_dir = "/astro/dboyea/figures"; mkpath(fig_dir)

# ╔═╡ e032f8e2-28b5-46c4-896e-27da9dfff22f
params = TOML.parsefile("/astro/dboyea/sculptor/isolation/1e6/stars/$name.toml")

# ╔═╡ 171b307c-cc13-4928-9526-a632031554e0
model_dir = "/astro/dboyea/sculptor/orbits/orbit1"

# ╔═╡ c47ce6da-46b1-46e0-bd7a-acc5f6bdb92b
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ 570e2840-977b-46c6-b974-52b0ed7decf3
starsfile = "/arc7/home/dboyea/sculptor/isolation/1e6/halos/V50_r0.5/stars/$(name)_stars.hdf5"

# ╔═╡ 22323266-c04e-4865-bd1a-a699dfa1a7d3
begin 
	using HDF5

	f = h5open(starsfile)
	p_idx = f["index"][:]
	probabilities = f["probabilities"][:][sortperm(p_idx)]
	p_idx = sort(p_idx)
	close(f)
	
end

# ╔═╡ 3df5d132-1a1a-40d6-b2d8-dcce2fcc4c34
out = lguys.Output("$model_dir/out", weights=probabilities)

# ╔═╡ a25f3edf-8d5b-401d-b727-8d8cd2d14910
snap_i = out[1]

# ╔═╡ fb302f97-6068-4fd8-b458-925aa952a8c2
x_cen = out.x_cen

# ╔═╡ 5c2534ca-ed11-47b4-9635-a132fd2a35e6
snap_f = out[orbit_props["idx_f"]]

# ╔═╡ 7042352d-fff5-4985-97e9-26b3f1b2fd6a
function make_sample(snap; 
	cen=nothing, rel_p_cut=0, r_max=Inf,
	frame=lguys.GSR
)
	filt = snap.weights .>= rel_p_cut * maximum(snap.weights)
	
	obs_df = lguys.to_sky(snap[filt], SkyFrame=frame, add_centre=true)

	obs_df[!, :xi], obs_df[!, :eta] = lguys.to_tangent(obs_df.ra, obs_df.dec, obs_df.ra[1], obs_df.dec[1])

	obs_df[!, :r_ell] = 60lguys.calc_r_ell(obs_df.xi, obs_df.eta, 0, 0)

	return obs_df[:, Not(:frame)]
end

# ╔═╡ bc729be3-3af9-447e-ab4f-061297fef3b0
obs_df_all = make_sample(snap_f, rel_p_cut=0, r_max=Inf)

# ╔═╡ 41d07a5b-32a7-43f2-9a34-eab54c8ab4e0
R_s = params["profile_kwargs"]["R_s"]

# ╔═╡ 3b3a4fc6-5bf3-4b81-8621-8f425c7dd517
R_s_arcmin = lguys.kpc_to_arcmin(R_s, 82.3)

# ╔═╡ f700b0bb-64b5-4be0-bf8a-88e72a7500e1
prof_f = lguys.ObsProfile("$model_dir/stars/$(name)_today_profile.toml")

# ╔═╡ 2e81adf7-47f4-4f61-857e-8997c15dc943
prof_i = lguys.ObsProfile("$model_dir/stars/$(name)_i_today_profile.toml")

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

	errscatter!(prof_expected.log_r, prof_expected.log_Sigma,
		yerr=prof_expected.log_Sigma_err,
		label="Jensen+24",
		color=fg
	)

	

	lines!(prof_i.log_r, prof_i.log_Sigma, 
			label="initial", color=COLORS[3])
	
	lines!(prof_f.log_r, prof_f.log_Sigma, 
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
df_j24 = lguys.load_fits("/astro/dboyea/dwarfs/notebooks/density_fits/sculptor/fiducial_sample.fits")

# ╔═╡ 03fba963-206a-4025-aef3-afe4ee50d5a6
let
	dr = 5
	fig, ax = xi_eta_axis(dr, dr)
	obs_df = obs_df_all

	bins = LinRange(-dr, dr, 100)
	dv = 15

	v0 = obs_df.radial_velocity[1]
	h = mean_2d(obs_df,  obs_df.radial_velocity, bins=bins, val_min=1e-20)

	p = heatmap!(h.xbins, h.ybins, h.values, 
	colormap=:redsblues,  colorrange=(v0-dv, v0+dv))

	Colorbar(fig[1, 2], p, label="GSR radial velocity")

	scatter!(df_j24.xi, df_j24.eta, color=:black, markersize=1)

	arc!(Point2f(0, 0), r_b/60, 0, 2π, color=:black)
	fig
end

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

	scatter!(df_j24.xi, df_j24.eta, color=:black, markersize=1)

	arc!(Point2f(0, 0), r_b/60, 0, 2π, color=:black)

	save("$fig_dir/skydensity.svg", fig)
	fig
end

# ╔═╡ 5e49981d-b456-4ae2-825d-bd30503de4eb
let
	ratio = (594, 841)
	fig = Figure(
		padding=0,
		backgroundcolor=:transparent

	)
	
	r_max = 160
	bandwidth=0.5

	dpcm = 5

	r_scale = 1/2
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
		colormap=Arya.get_arya_cmap()
	)

	hidedecorations!(ax)
	hidespines!(ax)

	save(joinpath(fig_dir, "xy_projection.png"), fig)

	#colsize!(fig.layout, 1, Aspect(1, 1/sqrt(2)))

	resize_to_layout!(fig)
	fig
end

# ╔═╡ Cell order:
# ╠═e9d818c0-1897-11ef-1913-1b38c9d4b4d6
# ╠═0f7680b9-4529-4438-8a5e-a9d6e6530eaf
# ╠═91564a7b-71f9-4619-9d51-69caf49af642
# ╠═b2f25a23-71bc-4217-8bef-c49d0cbc6676
# ╠═576b1ea1-05bd-48ce-b977-276b9271fc1f
# ╠═d794992b-328b-4327-96fa-67564dc8db98
# ╠═c64b9fe9-b3ff-498f-85dd-55e227443b63
# ╠═a25f3edf-8d5b-401d-b727-8d8cd2d14910
# ╠═6d207751-695f-4247-94da-ced5a146092f
# ╠═88f31bfc-a67d-4654-af8b-46dc91500558
# ╠═de26b439-d9dd-449e-9289-9d3daed87cc7
# ╠═fb302f97-6068-4fd8-b458-925aa952a8c2
# ╠═e032f8e2-28b5-46c4-896e-27da9dfff22f
# ╠═171b307c-cc13-4928-9526-a632031554e0
# ╠═3df5d132-1a1a-40d6-b2d8-dcce2fcc4c34
# ╠═c47ce6da-46b1-46e0-bd7a-acc5f6bdb92b
# ╠═5c2534ca-ed11-47b4-9635-a132fd2a35e6
# ╠═570e2840-977b-46c6-b974-52b0ed7decf3
# ╠═22323266-c04e-4865-bd1a-a699dfa1a7d3
# ╠═7042352d-fff5-4985-97e9-26b3f1b2fd6a
# ╠═bc729be3-3af9-447e-ab4f-061297fef3b0
# ╠═41d07a5b-32a7-43f2-9a34-eab54c8ab4e0
# ╠═3b3a4fc6-5bf3-4b81-8621-8f425c7dd517
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
# ╠═03fba963-206a-4025-aef3-afe4ee50d5a6
# ╠═01b0f4e8-9f0d-457b-900a-db4c85ec602d
# ╠═5e49981d-b456-4ae2-825d-bd30503de4eb
