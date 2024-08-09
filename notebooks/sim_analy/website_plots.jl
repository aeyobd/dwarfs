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

# ╔═╡ 605fdd1f-f402-416b-bd33-e2c611b502c1
fg = :black

# ╔═╡ 3c20733c-0369-4a45-b984-6167879718e7
bg = :transparent

# ╔═╡ d26647a6-8d24-48ec-86e4-de582874eaec
web_theme = Theme(
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
		# fonts = (;
		# 	regular = font
		# )
	)

# ╔═╡ b2f25a23-71bc-4217-8bef-c49d0cbc6676
begin 
	Arya.__init__()

	
	update_theme!(web_theme
	)
		
end

# ╔═╡ 4a0bafed-85d0-45d8-b00f-21a41d572955


# ╔═╡ 576b1ea1-05bd-48ce-b977-276b9271fc1f
font = "/astro/dboyea/fonts/Arev.ttf"

# ╔═╡ d794992b-328b-4327-96fa-67564dc8db98
Makie.load_font(font)

# ╔═╡ c64b9fe9-b3ff-498f-85dd-55e227443b63
name = "exp2d_rs0.1"

# ╔═╡ 6d207751-695f-4247-94da-ced5a146092f
prof_expected = lguys.ObsProfile("/arc7/home/dboyea/dwarfs/notebooks/density_fits/sculptor/fiducial_profile.toml")

# ╔═╡ 88f31bfc-a67d-4654-af8b-46dc91500558
r_b = 58.69

# ╔═╡ de26b439-d9dd-449e-9289-9d3daed87cc7
fig_dir = "figures"; mkpath(fig_dir)

# ╔═╡ 116795da-77f9-417f-8897-e771d4c2a94d
abspath(fig_dir)

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

# ╔═╡ 436045c2-e85e-40e8-b40b-955fc3c18848
extrema(lguys.calc_r(out.x_cen))

# ╔═╡ a25f3edf-8d5b-401d-b727-8d8cd2d14910
snap_i = out[1]

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

# ╔═╡ e0ca4c4c-8b2c-4049-a7d3-3ac4fa26ae3f
let 
	x_cen = out.x_cen
	times = out.times * T2GYR
	
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr", ylabel="r / kpc",
		xgridvisible=false,
		ygridvisible=false,
	)
	x = x_cen[1, :]
	y = x_cen[2, :]
	z = x_cen[3, :]
	r = @. sqrt(x^2 + y^2 + z^2)
	lines!(times, r, color=times)

	save("$fig_dir/orbit_r_t.svg", fig)
	fig
end

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
	axislegend(position=:lb,
	)
	
	colsize!(fig.layout, 1, Aspect(1, 1.3))

	#resize_to_layout!(fig)

	save(joinpath(fig_dir, "density_i_f.svg"), fig)

	fig
end

# ╔═╡ 693cb85b-8f19-414f-8a26-da2035232df0
abspath(fig_dir)

# ╔═╡ 15cc735b-4da6-496e-acf6-3b090758935d
idx_f = orbit_props["idx_f"]

# ╔═╡ aee39b8e-86ef-48e4-90e6-4432c67d453c
let 
	x_cen = out.x_cen
	times = out.times * T2GYR
	
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="R / kpc", ylabel="z / kpc",
		aspect=DataAspect(),
		xgridvisible=false,
		ygridvisible=false,
	)
	x = x_cen[1, :]
	y = x_cen[2, :]
	z = x_cen[3, :]
	R = @. sqrt(x^2 + y^2)
	lines!(R, z, color=times)

	scatter!(R[idx_f], z[idx_f], color=COLORS[2])
	scatter!(R[1], z[1], color = COLORS[3], marker=:rtriangle, )

	save("$fig_dir/orbit.svg", fig)
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

# ╔═╡ 96dc1860-3df3-4803-a8bf-b9228161583a
sum(df_j24.source_id .∈ [[5006419626331394048, 5026130884816022016]])

# ╔═╡ 03fba963-206a-4025-aef3-afe4ee50d5a6
let
	dr = 5
	fig, ax = xi_eta_axis(dr, dr)
	obs_df = obs_df_all

	bins = LinRange(-dr, dr, 100)
	dv = 15

	v0 = obs_df.radial_velocity[1]
	h = mean_2d(obs_df,  obs_df.radial_velocity, bins=bins, val_min=1e-20)

	h.values .-= v0
	v0 = 0
	p = heatmap!(h.xbins, h.ybins, h.values, 
	colormap=Reverse(:redsblues),  colorrange=(v0-dv, v0+dv))

	Colorbar(fig[1, 2], p, label=L"$\delta v_\textrm{rad, GSR}$")

	#scatter!(df_j24.xi, df_j24.eta, color=:black, markersize=1)
	save("$fig_dir/skydensity_velocity.svg", fig)

	arc!(Point2f(0, 0), r_b/60, 0, 2π, color=:black)
	fig
end

# ╔═╡ 232bc6db-e6d7-49a7-9894-ee8f782a90da


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
	fig = Figure(
		padding=0,
		backgroundcolor=:transparent

	)
	
	r_max = 100
	bandwidth=0.5

	dpcm = 5

	r_scale = 1/2
	limits = ((-r_max, r_max), (-r_max, r_max)
	)
	
	ax = Axis(fig[1,1], aspect=DataAspect(),
		xlabel = "y / kpc", ylabel="z / kpc", 
		#padding=0,
		backgroundcolor=:transparent,
		#title="dark matter", 
		)

	x, y = snap_f.positions[2, :], snap_f.positions[3, :]
	
	hi = kde([x y], bandwidth=(bandwidth, bandwidth), boundary=limits,
		npoints= (dpcm * r_max, dpcm * r_max) 
	)

	z = hi.density
	#z = log10.(1e-5  .+ hi.density ./ maximum(hi.density))
	
	h = heatmap!(hi.x, hi.y, z,
		colormap=Arya.get_arya_cmap(),
		colorscale=log10,
		colorrange=(1e-5*maximum(z), maximum(z)),
	)

	Colorbar(fig[1, 2], h, label=L"log $\rho_\textrm{DM}$ / $10^{10}$\,M$_\odot$ kpc$^{-3}$")

	save(joinpath(fig_dir, "xy_projection.png"), fig)

	#colsize!(fig.layout, 1, Aspect(1, 1/sqrt(2)))

	resize_to_layout!(fig)
	fig
end

# ╔═╡ 3083a583-e16a-44fc-af01-f68e7b54948c


# ╔═╡ 1c668565-0996-418d-983f-53aec9324e7e
let
	fg = "#6AA5D5"

	
	with_theme(Theme(
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
			regular = "Avenir"
		)
	)) do
	
		fig = Figure()
	
	r_max = 100
	bandwidth=0.5

	dpcm = 5

	r_scale = 1/2
	limits = ((-r_max, r_max), (-r_max, r_max)
	)
	
	ax = Axis(fig[1,1], aspect=DataAspect(),
		xlabel = "y / kpc", ylabel="z / kpc", 
		#padding=0,
		backgroundcolor=:transparent,
		#title="dark matter", 
		)

	x, y = snap_f.positions[2, :], snap_f.positions[3, :]
	
	hi = kde([x y], bandwidth=(bandwidth, bandwidth), boundary=limits,
		npoints= (dpcm * r_max, dpcm * r_max) 
	)

	z = hi.density
	#z = log10.(1e-5  .+ hi.density ./ maximum(hi.density))
	
	h = heatmap!(hi.x, hi.y, z,
		colormap=Arya.get_arya_cmap(),
		colorscale=log10,
		colorrange=(1e-5*maximum(z), maximum(z)),
	)

	Colorbar(fig[1, 2], h, label="dark matter density")

	#colsize!(fig.layout, 1, Aspect(1, 1/sqrt(2)))

	resize_to_layout!(fig)
	fig
	end
end

# ╔═╡ 3a569321-de5b-4e8e-b3ab-0a0d9bd9c964
prof_dm = lguys.Profiles3D("$model_dir/out/profiles.hdf5")

# ╔═╡ 2c697fc9-bc0c-4d7f-bed8-077be11b94d1
let 
	prof_i = prof_dm.profiles[1]

	prof_f = prof_dm.profiles[length(prof_dm.profiles)]
	r_max = prof_i.r_circ_max
	v_max = prof_i.v_circ_max
	
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=L"\log \; r / \textrm{kpc}", 
		ylabel=L"$v_\textrm{circ}$ / km s$^{-1}$",
		yscale=log10,
		yticks=[1, 10, 20, 30, 40, 50, 60],
		yminorticks=[1:9; 10:2:60],
		limits=((-1.2, 3), (5, 40)),
		xgridvisible=false,
		ygridvisible=false
	)


	lines!(prof_i.log_r, prof_i.v_circ * V2KMS, label="initial")

	lines!(prof_f.log_r, prof_f.v_circ * V2KMS, label="final")


	α = 0.4
	β = 0.65
	x = LinRange(1, 0.1, 100)

	y = @. 2^α * x^β * (1 + x^2)^(-α)
	lines!(log10.(x .* r_max[1]), y .* v_max[1] * V2KMS,  label="EN21",
	color=:black, linestyle=:dash)

	lines!(log10.(r_max), v_max * V2KMS, color=Arya.COLORS[4], label=L"v_\textrm{circ,\ max}")



	vs = [prof.v_circ_max for prof in prof_dm.profiles]
	rs = [prof.r_circ_max for prof in prof_dm.profiles]

	scatter!(log10.(rs), vs * V2KMS, 
		color=COLORS[3],
		label="V circ max"
	)

	
	axislegend(ax, position=:rt)
	save(joinpath(fig_dir, "v_circ_profiles.svg"), fig)
	fig
end

# ╔═╡ Cell order:
# ╠═e9d818c0-1897-11ef-1913-1b38c9d4b4d6
# ╠═0f7680b9-4529-4438-8a5e-a9d6e6530eaf
# ╠═91564a7b-71f9-4619-9d51-69caf49af642
# ╠═d26647a6-8d24-48ec-86e4-de582874eaec
# ╠═605fdd1f-f402-416b-bd33-e2c611b502c1
# ╠═3c20733c-0369-4a45-b984-6167879718e7
# ╠═b2f25a23-71bc-4217-8bef-c49d0cbc6676
# ╠═e0ca4c4c-8b2c-4049-a7d3-3ac4fa26ae3f
# ╠═aee39b8e-86ef-48e4-90e6-4432c67d453c
# ╠═436045c2-e85e-40e8-b40b-955fc3c18848
# ╠═4a0bafed-85d0-45d8-b00f-21a41d572955
# ╠═116795da-77f9-417f-8897-e771d4c2a94d
# ╠═576b1ea1-05bd-48ce-b977-276b9271fc1f
# ╠═d794992b-328b-4327-96fa-67564dc8db98
# ╠═c64b9fe9-b3ff-498f-85dd-55e227443b63
# ╠═a25f3edf-8d5b-401d-b727-8d8cd2d14910
# ╠═6d207751-695f-4247-94da-ced5a146092f
# ╠═88f31bfc-a67d-4654-af8b-46dc91500558
# ╠═de26b439-d9dd-449e-9289-9d3daed87cc7
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
# ╠═6c94b8e5-3979-42e8-ae8f-61906c72c3db
# ╠═4a5adacf-db5e-4dc0-90c9-619e9b571ba8
# ╠═ce3b2336-fce2-473f-bfe5-c0f07cfab6e0
# ╠═96dc1860-3df3-4803-a8bf-b9228161583a
# ╠═03fba963-206a-4025-aef3-afe4ee50d5a6
# ╠═232bc6db-e6d7-49a7-9894-ee8f782a90da
# ╠═01b0f4e8-9f0d-457b-900a-db4c85ec602d
# ╠═5e49981d-b456-4ae2-825d-bd30503de4eb
# ╠═3083a583-e16a-44fc-af01-f68e7b54948c
# ╟─1c668565-0996-418d-983f-53aec9324e7e
# ╠═3a569321-de5b-4e8e-b3ab-0a0d9bd9c964
# ╠═2c697fc9-bc0c-4d7f-bed8-077be11b94d1
