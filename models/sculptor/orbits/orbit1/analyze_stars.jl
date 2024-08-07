### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 340ffbbe-17bd-11ef-35c6-63505bb128b7
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	import LilGuys as lguys

	using Arya
end

# ╔═╡ cf6a7cbb-1034-4026-a3f3-1e854d2929e2
using FITSIO, Tables

# ╔═╡ f0d2b68a-fae2-4486-a434-a8816e400e84
import TOML

# ╔═╡ 8ed2e715-1095-4c99-9ecb-54fecbc27d7f
pwd()

# ╔═╡ 436a5be3-f597-4fc4-80a8-dc5af302ad66
orbit_props = TOML.parsefile("orbital_properties.toml")

# ╔═╡ 29988108-b02c-418c-a720-5766f47c39ff
starsname = "exp2d_rs0.1"

# ╔═╡ f0d74eaa-81e9-4b04-9765-24a0935b1430
starsfile = "../../isolation/1e6/stars/$(starsname)_stars.hdf5"

# ╔═╡ f9fe37ef-de81-4d69-9308-cda968851ed2
begin 
	using HDF5

	f = h5open(starsfile)
	p_idx = f["index"][:]
	probabilities = f["probabilities"][:][sortperm(p_idx)]
	p_idx = sort(p_idx)
	close(f)
	
end

# ╔═╡ 377284f2-dcee-44d3-9a04-728605cea92a
md"""
Given a stellar probability file, calculates initial-final density profiles, 
and projects stars onto the sky
"""

# ╔═╡ 1b5c00d2-9df6-4a9c-ae32-05abcbf0e41a
paramsfile = "../../isolation/1e6/stars/$starsname.toml"

# ╔═╡ 84dc77f7-14b3-4a2e-a556-c025d7df0095
params = TOML.parsefile(paramsfile)

# ╔═╡ 5e12d306-a430-4b15-b3a7-d4806a5856cd
name = splitext(basename(starsfile))[1]

# ╔═╡ ef3481f8-2505-4c04-a04a-29bdf34e9abb
outfile = "$(starsname)_today.fits"

# ╔═╡ a80d9e83-6b11-4c55-94ec-294d4247af42
outfile_i = "$(starsname)_i_today.fits"

# ╔═╡ 172588cc-ae22-440e-8488-f508aaf7ce96
rel_p_cut = 1e-15

# ╔═╡ 973955ad-3214-42cf-831f-a782f1d2434a
idx_i = 1 

# ╔═╡ 8f9ee2a7-de43-43e0-8257-93ecc630044f
idx_f = orbit_props["idx_f"]

# ╔═╡ 6fb4b7e1-a22c-4ff8-bbe9-dbf5de5acd37
begin 
	out =  lguys.Output("out/combined.hdf5")
	
	cens = CSV.read("out/centres.csv", DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
	out.x_cen .= x_cen
	out.v_cen .= v_cen
	cens[!, "v_x"] .= cens.vx
	cens[!, "v_y"] .= cens.vy
	cens[!, "v_z"] .= cens.vz
	out
end

# ╔═╡ 6feeaae2-cb01-46ad-ad1d-daaca1caf7ec
snap_f = out[idx_f]

# ╔═╡ 396a53a3-de0f-4d97-9693-40f3757d66f9
snap_i = out[idx_i]

# ╔═╡ 77479cd4-513c-4603-9aa0-1acd964c403a
let
	fig = Figure()
	r_max = 100
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "y / kpc", ylabel="z / kpc", title="stars",
	limits=(-r_max, r_max, -r_max, r_max))
	bins = LinRange(-r_max, r_max, 300)

	hm = Arya.hist2d!(ax, snap_f.positions[2, :], snap_f.positions[3, :], bins = bins, colorscale=log10, colorrange=(1e-10, 0.06), weights=probabilities[snap_f.index])

	Colorbar(fig[:, end+1], hm)
	fig
end

# ╔═╡ e33245ad-7562-47af-94b9-9708af45b2a4


# ╔═╡ 7bc2c15c-33f7-43f3-a47f-ca39ffc22071
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc", title="initial")

	bin_range = LinRange(-2, 2, 100)
	colorrange =(1e-5, nothing)

	bins = (x_cen[1, idx_i]  .+ bin_range,  x_cen[2, idx_i]  .+ bin_range)
		
	probs = probabilities[snap_i.index]

	Arya.hist2d!(ax, snap_i.positions[1, :], snap_i.positions[2, :], weights=probs, bins = bins, colorscale=log10, colorrange=colorrange)

	ax2 = Axis(fig[1,2], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc",
	title="final")

	probs = probabilities[snap_f.index]

	bins = (x_cen[1, idx_f]  .+ bin_range,  x_cen[2, idx_f]  .+ bin_range)

	Arya.hist2d!(ax2, snap_f.positions[1, :], snap_f.positions[2, :], weights=probs, bins = bins, colorscale=log10, colorrange=colorrange)
	
	fig
end

# ╔═╡ 8df71dd5-ba0b-4918-9dc0-791c6c3bad5f
let 
	fig = Figure()
	ax = Axis(fig[1,1], yscale=log10, limits=(nothing, (1e-8, 1)),
		xlabel=L"\epsilon", ylabel="count")
	stephist!(lguys.calc_ϵ(snap_i), weights=probabilities[snap_i.index])
	es = lguys.calc_ϵ(snap_f)
	filt = es .> 0
	es = es[filt]
	probs = probabilities[snap_f.index][filt]
	stephist!(es, weights = probs)

	fig

end

# ╔═╡ c6ce37f4-09f6-43c2-9ee3-286d9e6baf7f
md"""
# 3D density profiles
"""

# ╔═╡ a1138e64-5fa5-4a0e-aeef-487ee78a7adc
function plot_ρ_s!(snap; bins=40, kwargs...)
	pos = lguys.extract_vector(snap, :positions, p_idx)
	rs = lguys.calc_r(pos, snap.x_cen)
	r, ρ = lguys.calc_ρ_hist(rs, bins, weights=probabilities)
	lines!(log10.(lguys.midpoint(r)), log10.(ρ); kwargs...)
end

# ╔═╡ 91daae57-94dc-4a8e-b981-75f1406e0768
begin
	profile_class = getproperty(lguys, Symbol(params["profile"]))
	profile = profile_class(;lguys.dict_to_tuple(params["profile_kwargs"])...)
	ρ_s(r) = lguys.calc_ρ(profile, r)
end

# ╔═╡ f5a3ea2f-4f19-4b0e-af55-74902f2c6485
md"""
# Sky projection
"""

# ╔═╡ 51dea031-015b-4506-879d-9245c122d608
snap_cen = lguys.Snapshot(x_cen, v_cen, ones(size(x_cen, 2)))

# ╔═╡ f87e4610-34ac-49f9-9858-0b3ef72eef15
cen = (cens[idx_f, :])

# ╔═╡ c351295d-1fb9-4088-8a27-5c732924959e
cen.vx

# ╔═╡ 24636618-7310-4557-baac-01d3d5d076dd
cen

# ╔═╡ 46f1792a-53a4-4540-b8f3-93747f071ed0
obs_c_galcen = lguys.Galactocentric(x=cen.x, y=cen.y, z=cen.z, 
	v_x=cen.vx*lguys.V0, v_y=cen.vy*lguys.V0, v_z=cen.vz*lguys.V0)

# ╔═╡ 7e588ae3-89f3-4b91-8963-f6bf4391a859
function save_obs(obs_df, outfile)
	FITS(outfile, "w") do f
		df = Dict(String(name) => obs_df[:, name] for name in names(obs_df))
		write(f, df)
	
		println("written to $outfile")

		df
	end
end

# ╔═╡ 2612e3a2-6e6e-494e-b140-720dd2db6ec2
obs_today_file = TOML.parsefile("../../mc_orbits/orbit1.toml")

# ╔═╡ 7c4c5136-17d5-4dc5-9e5c-25e2348d2a84
obs_today_icrs = lguys.ICRS(;
	ra=obs_today_file["ra"], dec=obs_today_file["dec"],
	distance=obs_today_file["distance"],
	pm_ra=obs_today_file["pm_ra"],
	pm_dec=obs_today_file["pm_dec"],
	radial_velocity=obs_today_file["radial_velocity"],
)

# ╔═╡ 494a0dda-0a31-4049-a089-d78437c209cc
obs_today_err = lguys.ICRS(;
	ra=0, dec=0,
	distance=obs_today_file["distance_err"],
	pm_ra=obs_today_file["pm_ra_err"],
	pm_dec=obs_today_file["pm_dec_err"],
	radial_velocity=obs_today_file["radial_velocity_err"],
)

# ╔═╡ d2de33be-b037-4255-b9f4-da29abb23754
frame = lguys.HelioRest

# ╔═╡ 01b5ad85-9f37-4d8b-a29d-e47526f112ec
function make_sample(snap; 
	cen=nothing, rel_p_cut=rel_p_cut, r_max=Inf,
	Frame=frame
)

	ps = probabilities[snap.index]
	
	p_min = rel_p_cut * maximum(ps)
	println("adopting p_min = $p_min")

	filt = ps .> p_min
	snap_stars = snap[filt]
	ps = ps[filt]
	println("sanity check: ", ps == probabilities[snap_stars.index])
	println("number of final stars: ", length(snap_stars))
	
	obs_pred = lguys.to_sky(snap_stars, SkyFrame=Frame)
	
	obs_df = DataFrame(; 
	collect(key => [getproperty(o, key) for o in obs_pred]
			for key in [:ra, :dec, :pm_ra, :pm_dec, :distance, :radial_velocity])...
		
	)

	obs_df[!, "index"] = snap_stars.index
	obs_df[!, "probability"] = ps

	if cen !== nothing
		obs_c_galcen = lguys.Galactocentric(x=cen.x, y=cen.y, z=cen.z, 
			v_x=cen.v_x, v_y=cen.v_y, v_z=cen.v_z)
		obs_c = lguys.transform(lguys.ICRS, obs_c_galcen)
	
		cen_df = DataFrame(ra=obs_c.ra, dec=obs_c.dec, pm_ra=obs_c.pm_ra, pm_dec=obs_c.pm_dec, radial_velocity=obs_c.radial_velocity, distance=obs_c.distance, index=-1, probability=0.0)
	
		obs_df = append!(cen_df, obs_df)
		
		obs_df[!, "xi"], obs_df[!, "eta"] = lguys.to_tangent(obs_df.ra, obs_df.dec, obs_c.ra, obs_c.dec)
		obs_df[!, "r_ell"] = @. 60 * sqrt(obs_df.xi^2 + obs_df.eta^2)

	end
	
	rename!(obs_df, "pm_ra"=>"pmra")
	rename!(obs_df, "pm_dec"=>"pmdec")

	return obs_df
end

# ╔═╡ 1f722acb-f7b9-4d6c-900e-11eae85e0708
obs_df = make_sample(snap_f, cen = cens[idx_f, :], Frame=frame)

# ╔═╡ 0215fe76-f9c5-4274-ae43-89960a0caaef
obs_c = obs_df[1, :]

# ╔═╡ 5e17bf16-2e3d-4372-8283-c43b8ad2550d
save_obs(obs_df, outfile)

# ╔═╡ 249351c0-6fb6-49ee-ab44-e464f34a1bbe
sky_orbit = make_sample(snap_cen, Frame=frame)

# ╔═╡ b1b218ff-694b-48bc-b999-806330ad4308
obs_today = lguys.transform(frame, obs_today_icrs)

# ╔═╡ 8e32cd34-e97e-4cb7-95ba-67f4ed0c9aed
obs_gc = lguys.transform(lguys.Galactocentric, obs_today
)
	

# ╔═╡ 1ef7e0c1-56f3-4be7-a138-c088ff973ec6
obs_gc.x

# ╔═╡ 1b528c3f-877b-40f5-bc05-7b1c96c8d5c9
let
	global obs_df_i

	snap_i_shifted = copy(snap_i)
	x_c = [obs_gc.x, obs_gc.y, obs_gc.z]
	v_c = [obs_gc.v_x, obs_gc.v_y, obs_gc.v_z] / lguys.V0

	snap_i_shifted.positions .+= x_c .- snap_i.x_cen

	snap_i_shifted.velocities .+= v_c .- snap_i.v_cen

	obs_df_i = make_sample(snap_i_shifted, cen=obs_gc)
end

# ╔═╡ ff759ae0-35c3-460e-97df-f865e26a0f48
save_obs(obs_df_i, outfile_i)

# ╔═╡ 32852a2a-1742-4c37-9430-d1425e4f2521
m_star_f = obs_df.probability

# ╔═╡ 8ad01781-8b5d-4d57-a0b5-7a445fb09b5b
let 
	fig = Figure()

	dy = 5
	bins = 100

	
	dx = dy * 1/cosd(obs_today.dec)
	limits = (obs_today.ra .+ (-dx, dx), obs_today.dec .+ (-dy, dy))
	
	ax = Axis(fig[1,1],
		xlabel="RA / degrees",
		ylabel="dec / degrees",
		limits=limits,
		aspect = 1,
		
	)

	x = obs_df.ra
	y = obs_df.dec


	hi = Arya.histogram2d(x, y, bins, weights=m_star_f, limits=limits)
	areas = diff(hi.xbins) .* (diff(hi.ybins)')
	hi.values ./= areas
	
		
	h = heatmap!(hi, colorscale=log10, colorrange=(1e-10, maximum(hi.values)))
	errscatter!([obs_today.ra], [obs_today.dec], color=COLORS[3], size=10)

	idx = idx_f - 20: idx_f + 20
	lines!(sky_orbit.ra[idx], sky_orbit.dec[idx])
	
	Colorbar(fig[1, 2], h,
		label="stellar density"
	)
	fig
end

# ╔═╡ 0487a830-3152-4d43-89ed-7562dc2dea52
let 
	fig = Figure()

	bins = 100
	dy = 5
	dx = dy * 1/cosd(obs_today.dec)
	dv = 10
	limits = (obs_today.ra .+ (-dx, dx), obs_today.dec .+ (-dy, dy))

	ax = Axis(fig[1,1],
		xlabel="RA / degrees",
		ylabel="dec / degrees",
		limits=limits,
		aspect = 1,
		xgridvisible=false,
		ygridvisible=false
		
	)

	x = obs_df.ra
	y = obs_df.dec

	v_mean = lguys.mean(obs_df.radial_velocity, lguys.weights(obs_df.probability))

	delta_v = obs_df.radial_velocity .- v_mean
	
	h_vel = Arya.histogram2d(x, y, bins, weights=m_star_f .* delta_v, limits=limits)
	h_mass = Arya.histogram2d(x, y, bins, weights=m_star_f, limits=limits)

	h_vel.values ./= h_mass.values
	
	h = heatmap!(h_vel, colormap=Reverse(:bluesreds), colorrange=(-dv, dv))

	idx = idx_f - 20: idx_f + 20
	lines!(sky_orbit.ra[idx], sky_orbit.dec[idx], color=:black)

	pm_vec = 5 .* (sky_orbit.pmra[idx_f], sky_orbit.pmdec[idx_f])
	arrows!([sky_orbit.ra[idx_f]], 
		[sky_orbit.dec[idx_f]], 
		[pm_vec[1]], [pm_vec[2]], 
		color=COLORS[3],
		linewidth=2
	)
	
	Colorbar(fig[1, 2], h, 
		label = L"$\Delta \tilde{v}_\textrm{rad}$ / km s$^{-1}$",
	)
	fig
end

# ╔═╡ d89d2f21-e6c2-4911-8dd3-480b8e51020b
out.times * lguys.T0

# ╔═╡ 2ef43371-bfae-4a65-8fa9-d1ab5ade32f1
let
	println("testing pm if works")

	i1 = idx_f
	i2 = idx_f - 1
	dt = (out.times[i1] - out.times[i2]) * 1e9 * lguys.T0 # years
	
 	ddec = sky_orbit.dec[i1] - sky_orbit.dec[i2] 
	dra = sky_orbit.ra[i1] - sky_orbit.ra[i2] 

	deg_to_mas = 60^2 * 1e3
	
	pmdec = ddec * deg_to_mas / dt 
	pmra = dra * deg_to_mas / dt  * cos(deg2rad(sky_orbit.dec[i1]))

	println(pmra, ", ", pmdec)
	println(sky_orbit.pmra[i1], ", ", sky_orbit.pmdec[i1])
end

# ╔═╡ 4537d2a8-dc90-4106-81c5-4a230734b182
let 
	fig = Figure()

	dy = 5
	dx = dy * 1/cosd(obs_today.dec)
	ax = Axis(fig[1,1],
		xlabel="RA / degrees",
		ylabel="dec / degrees",
		limits=(obs_today.ra .+ (-dx, dx), obs_today.dec .+ (-dy, dy)),
		aspect = 1,
		
	)

	x = obs_df_i.ra
	y = obs_df_i.dec

	
	h = Arya.hist2d!(ax, x, y, weights=obs_df_i.probability, bins=100, colorscale=log10, colorrange=(1e-10, nothing))

	Colorbar(fig[1, 2], h)
	fig
end

# ╔═╡ 1c242116-e66d-453b-ad62-b6a65cdbe284
sky_orbit[idx_f, :]

# ╔═╡ b3e68e32-c058-467d-b214-aab6a4cd1e19
r_cut = 2 # degrees

# ╔═╡ e5c84f04-6df0-4f17-977c-7d31cd0b37dc
obs_df[!, :r] = sqrt.(obs_df.xi .^ 2 .+ obs_df.eta .^ 2)

# ╔═╡ 52abcdc2-4b9e-408f-9c1f-9dc0678c60ca
obs_df.xi

# ╔═╡ 15e303d8-6020-4cfd-ac01-419a253a4a0b
sum(obs_df.probability[obs_df.r .< r_cut])

# ╔═╡ 6858dc1b-6f60-4ca1-b2a2-80b4ca345b8c
sum(obs_df.probability[obs_df.r .> r_cut])

# ╔═╡ 9fed63d6-c139-4d28-b00c-37dc1b8dc004
scatter(obs_df.ra[obs_df.r .< r_cut], obs_df.dec[obs_df.r .< r_cut])

# ╔═╡ 9d74ffbb-4c31-4062-9434-7755f53e4da0
let
	dr = 0.1
	r_max = 2
	
	limits = (obs_today.pm_ra .+ (-dr, dr), obs_today.pm_dec .+ (-dr, dr))
	
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
	limits=limits,
	aspect=DataAspect())

	filt = obs_df.r .< r_max

	x = obs_df.pmra[filt]
	y = obs_df.pmdec[filt]
	
	h = Arya.hist2d!(ax, x, y, weights=m_star_f[filt], bins=100, 
		colorscale=log10, colorrange=(1e-8, nothing))
	
	errscatter!([obs_today.pm_ra], [obs_today.pm_dec], xerr=[obs_today_err.pm_ra], yerr=[obs_today_err.pm_ra], color=COLORS[3])
	
	scatter!([sky_orbit.pmra[idx_f]], 
		[sky_orbit.pmdec[idx_f]], markersize=10)

	idx = idx_f - 20: idx_f + 20
	lines!(sky_orbit.pmra[idx], sky_orbit.pmdec[idx])

	Colorbar(fig[1, 2], h, label="stellar density")
	
	fig
end

# ╔═╡ 975a2008-cf02-4442-9ee9-0b1bbb20889d
let
	dr = 0.1
	r_max = 2
	bins = 80
	
	limits = (obs_today.pm_ra .+ (-dr, dr), nothing)
	
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel="density",
	limits=limits)

	filt = obs_df.r .< r_max

	x = obs_df.pmra[filt]
	y = obs_df.pmdec[filt]
	
	h = Arya.histogram(x, bins, weights=m_star_f[filt], normalization=:pdf)

	
	barplot!(h)
	
	fig
end

# ╔═╡ 12a30334-899e-4061-b7d2-af8c2346721d
let 
	x = log10.(probabilities[probabilities .> 0])
	h = Arya.histogram(x)

	barplot(h)
end

# ╔═╡ 6fe6deb4-ae44-4ca0-9617-95841fdaf791
let
	
	fig = Figure()
	dx = 5
	dy = 40
	r_max = 2
	
	limits = (sky_orbit.distance[idx_f] .+ (-dx, dx), sky_orbit.radial_velocity[idx_f] .+ (-dy, dy))
	
	ax = Axis(fig[1,1],
		xlabel="distance / kpc",
		ylabel = L"$\tilde{v}_\textrm{rad}$ / km s$^{-1}$",
		limits=limits
	)



	filt = obs_df.xi .^ 2 .+ obs_df.eta .^ 2 .< r_max ^ 2
	x = obs_df.distance[filt]
	y = obs_df.radial_velocity[filt]
	
	h = Arya.hist2d!(ax, x, y, weights=m_star_f[filt], bins=100,
		colorscale=log10, colorrange=(1e-10, nothing)
	)

	errscatter!([obs_today.distance], [obs_today.radial_velocity], xerr=[obs_today_err.distance], yerr=[obs_today_err.radial_velocity], color=:red)

	scatter!([sky_orbit.distance[idx_f]], 
		[sky_orbit.radial_velocity[idx_f]], markersize=10)
	
	idx = idx_f - 20: idx_f + 20
	lines!(sky_orbit.distance[idx], sky_orbit.radial_velocity[idx])

	
	Colorbar(fig[1, 2], h)

	fig
end

# ╔═╡ 20e742e9-a909-4be0-822b-f4ca6015b8aa
let
	fig, ax = FigAxis(aspect=1,
		limits=(-2, 2, -2, 2),
		xlabel=L"\xi",
		ylabel=L"\eta",
		xgridvisible=false,
		ygridvisible=false
	)
	
	scatter!(obs_df.xi, obs_df.eta, 
		alpha=0.1, color=:black, markersize=3)
	fig
end

# ╔═╡ 25ccc096-1d27-43ce-8fe1-5a3a9d7cd54e
let
	fig = Figure()
	dx = 4
	dy = 20
	# limits = (obs["distance"] .+ (-dx, dx), obs["radial_velocity"] .+ (-dy, dy))
	
	ax = Axis(fig[1,1],
		xlabel="distance / kpc",
		ylabel = "radial velocity / km/s",
		# limits=limits
	)


	
	x = obs_df_i.distance
	y = obs_df_i.radial_velocity
	
	h = Arya.hist2d!(ax, x, y, weights=obs_df_i.probability, bins=100,
		colorscale=log10, colorrange=(1e-10, nothing)
	)

	Colorbar(fig[1, 2], h)
	
	fig
end

# ╔═╡ 55211531-6eb7-4689-9ecd-74854b8ad884
hist(obs_df_i.radial_velocity)

# ╔═╡ 487c6fda-5cb3-4ef4-a988-da2400c96cd4
obs_df_i.radial_velocity

# ╔═╡ e23e85f9-7667-4f6e-8af6-2516fa292e2b
md"""
# Velocity profile
"""

# ╔═╡ 0dd09c1e-67c0-4f23-bd12-41cbef62e4de
import StatsBase: weights, std

# ╔═╡ f348c1b6-5196-4546-ba22-72677611a420
begin 
	ra0 = obs_c.ra
	dec0 = obs_c.dec

	
	r_ell = obs_df.r_ell
	rv = obs_df.radial_velocity
end

# ╔═╡ d32757d2-dc08-488c-9ddd-d3eefefa2db7
function calc_σv(r_ell, rv, mass;  r_max = 30)
	filt = r_ell .< r_max
	vel = rv[filt]

	return std(vel, weights(mass)[filt])
end

# ╔═╡ 194ff30a-31a7-44bc-ac54-722a629457fc
σv = calc_σv(r_ell, rv, m_star_f, r_max=240)

# ╔═╡ fbd46bd2-79d7-460e-b0ab-0e34a68a1f0a
gaussian(x, μ, σ) = 1/sqrt(2π)* 1/σ * exp(-(x-μ)^2/(2σ^2))

# ╔═╡ ed097e29-d2dc-4fa7-94e5-483b380600cc
let
	r_max = 240
	
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="RV (km/s)",
		ylabel="density",
		title="within $r_max arcmin"
	)
	
	filt = r_ell .< r_max
	
	rv = rv[filt]
	mass = m_star_f[filt]

	μ = Arya.mean(rv, Arya.sb.weights(mass))
	
	h = Arya.histogram(rv, weights=mass, normalization=:pdf)

	scatter!(lguys.midpoint(h.bins), h.values)

	x_model = LinRange(80, 130, 100)
	y_model = gaussian.(x_model, μ, σv)
	lines!(x_model, y_model)
	fig
end

# ╔═╡ b9a4b2b7-be95-4ccb-ad74-9b761abfae8a
@doc raw"""
	calc_rb(\sigma, delta_t)

Given a radial velocity dispersion (in km/s), and the time since the last pericentre (in Gyr), returns the break radius (in kpc).

See peñarrubia et al. (200X) for the original equation.
```math
r_{\rm break} = C\,\sigma_{v}\,\Delta t
```
where $C=0.55$ is an empirical fit
"""
function calc_rb(σ, delta_t)
	kpc_per_kms_Gyr = 0.9777
	return 0.55 * σ * (delta_t)
end

# ╔═╡ 750185c9-a317-4978-aa04-486e5bfb7a63
r_b_kpc = calc_rb(σv, orbit_props["t_last_peri"]) # kpc

# ╔═╡ 1866c280-89c3-4a71-9dbf-50b583360145
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log r / \textrm{kpc}", ylabel = L"\log \rho_\star", 
		limits=((-1.9, 1), (-8, 2)))

	plot_ρ_s!(snap_i, bins=500, label="initial")
	plot_ρ_s!(snap_f, bins=500, label="final")

	x = LinRange(-2, 1, 1000)
	r = 10 .^ x
	y = log10.(ρ_s.(r))
	
	lines!(x, y)
	
	vlines!(log10.(r_b_kpc))
	
	axislegend(ax)
	fig
end

# ╔═╡ 02183d68-bf7a-4774-9b48-f50712eeb552
r_b_arcmin = 206265 / 60 * obs_c.distance / 1e3

# ╔═╡ 9b75409d-55f5-47c3-ab63-8168d31d3d54
md"""
# properties
"""

# ╔═╡ 6408828d-d570-4585-8fa8-24661857fb35
function get_M_h(output::lguys.Output, radius; idxs=(1:10:length(output)))
	N = length(idxs)
	M = Vector{Float64}(undef, N)
	
	for i in eachindex(idxs)
		snap = output[idxs[i]]
		ϵ = lguys.calc_ϵ(snap)
		filt = ϵ .> 0
		
		rs = lguys.calc_r(snap[filt])
		filt2 = rs .< radius

		ps = probabilities[snap.index[filt][filt2]]
		M[i] = sum(ps)
	end

	return M
end

# ╔═╡ df2fc50f-3a03-424e-ae88-02ae207a04ee
probabilities[out[1].index]

# ╔═╡ b36b594c-23b4-4683-8a12-3fa1b4c8b0d9
r_h = 0.1

# ╔═╡ 03cdff7d-6e86-45ab-9837-9029293639fa
probabilities

# ╔═╡ cdbed68e-0da2-4648-a5d2-61b5b07fb4a2
# ╠═╡ disabled = true
#=╠═╡
let
	fig = Figure()

	n_rh = 10
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="normalized stellar mass within $n_rh rh",
		# yscale=log10,
		#yticks=[1, 0.1],
	)
	

	M_dm_h = get_M_h(out, n_rh * r_h)
	scatter!(out.times[1:10:end] * lguys.T0, M_dm_h ./ M_dm_h[1])
	
	
	fig
end
  ╠═╡ =#

# ╔═╡ f3db0a1d-78bc-4ac4-84b6-d53f27a28795
out[10].x_cen

# ╔═╡ d42795d0-bd69-4c2c-be5b-e27e85199ee3
# ╠═╡ disabled = true
#=╠═╡
let
	fig = Figure()

	n_rh = 10
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel=L"\langle v_{x,\,\star} \rangle_{r \leq 1\,\textrm{kpc}} \ / \ \textrm{km\,s^{-1}}",
		# yscale=log10,
		#yticks=[1, 0.1],
	)

	idx = 1:10:length(out)
	vs = [calc_σ_v(out[i]) for i in idx]

	scatter!(out.times[idx] * lguys.T0, vs)
	
	fig
end
  ╠═╡ =#

# ╔═╡ a0391689-66a2-473f-9704-e12a3d033d13
import LinearAlgebra: dot

# ╔═╡ 44ab0a25-ab4c-4e90-8619-2a068a285755
function calc_v_rad(snap)
	x_vec = snap.positions .- snap.x_cen
	v_vec = snap.velocities .- snap.v_cen

	x_hat = x_vec / lguys.calc_r(x_vec)'

	v_rad = sum(x_hat .* v_vec, dims=1)
	return vec(v_rad)
end

# ╔═╡ c2ccf9de-e3cd-4950-9a4d-6f425d261ccb
function calc_v_tan(snap)
	v = lguys.calc_v(snap)
	v_rad = calc_v_rad(snap)
	v_tan = @. sqrt(v^2 - v_rad^2)
	return vec(v_tan)
end

# ╔═╡ 89ea28dc-2677-4e6f-9066-d5ffff3b474c
function calc_v_los(snap)
	(snap.velocities[1, :] .- snap.v_cen[1])
end

# ╔═╡ b4f7cbb0-fb5b-4dc5-b66e-679a8e5b630d
function calc_σ_v(snap; r_max=1)
	rs = lguys.calc_r(snap)
	filt = rs .< r_max
	vs = calc_v_los(snap[filt])
	masses = probabilities[snap[filt].index]
	return lguys.mean(vs .^ 2, weights(masses))^0.5 * lguys.V0
end

# ╔═╡ 76438e61-d227-41cf-b9ea-e658bc389772
let
	snap = snap_f
	
	mass = probabilities[snap.index]
	v_rad = calc_v_los(snap)  * lguys.V0
	logr = log10.(lguys.calc_r(snap))
	filt = logr .< 1

	h = Arya.histogram(v_rad[filt], 
		weights=mass[filt], normalization=:pdf, limits=(-20, 20))


	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"v_{x} \quad [\textrm{km\,s^{-1}}]",
		ylabel="stellar density"
	)
	
	scatter!(lguys.midpoint(h.bins), h.values)

	σ = calc_σ_v(snap)
	x_model = LinRange(-20, 20, 100)
	y_model = gaussian.(x_model, 0, σ)
	lines!(x_model, y_model, 
		color=COLORS[2], label="gaussian (σ = $(round(σ, digits=2)) km / s)"
	)

	axislegend()
	fig
	
end

# ╔═╡ 04629bcf-5c19-41eb-8903-72947c209cbf
let
	snap = snap_f
	
	mass = probabilities[snap.index]
	v_rad = calc_v_rad(snap) * lguys.V0
	logr = log10.(lguys.calc_r(snap))
	filt = logr .< 1

	h = Arya.histogram(v_rad[filt], 
		weights=mass[filt], normalization=:pdf, limits=(-20, 20))


	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"v_\textrm{rad} = \vec{v} \cdot \hat{r} \  / \ \textrm{km\,s^{-1}}",
		ylabel="stellar density"
	)
	
	scatter!(lguys.midpoint(h.bins), h.values)

	σ = calc_σ_v(snap)
	x_model = LinRange(-20, 20, 100)
	y_model = gaussian.(x_model, 0, σ*√3/2)
	lines!(x_model, y_model)
	fig
	
end

# ╔═╡ d664ab12-a2c1-4531-a4ff-250ffa3ce9eb
calc_σ_v(snap_i) / sqrt(3)

# ╔═╡ 6e34b91c-c336-4538-a961-60833d37f070
function v_rad_hist(snap, bins=40)

	mass = probabilities[snap.index]
	v_rad = calc_v_rad(snap)
	logr = log10.(lguys.calc_r(snap))
	h1 = Arya.histogram(logr, bins, weights=v_rad .* mass)
	h2 = Arya.histogram(logr, bins, weights=mass)

	x_bins = h1.bins
	v_bins = h1.values
	counts = h2.values

	return x_bins, v_bins ./ counts
end

# ╔═╡ 17c880ba-d6de-4d9c-b44b-d5ecdb4da226
calc_v_rad(snap_i)

# ╔═╡ 57d5decd-8259-4ec2-87ac-44d28625cd7b
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="log r",
		ylabel="mean radial velocity",
		limits=((nothing, 1), (-60, 60))
	)
	
	x, y= v_rad_hist(snap_f, 100)
	scatter!(lguys.midpoint(x), y * lguys.V0)

	vlines!(log10.(r_b_kpc))
	fig
end

# ╔═╡ 98bcfb2a-8eb5-4c7c-8dea-46d000a6d245
r_b_arcmin

# ╔═╡ 0c26f965-0381-4f79-a6ce-0772ae922b3f
let
	snap = snap_f
	
	mass = probabilities[snap.index]
	v_rad = calc_v_rad(snap) * lguys.V0
	logr = log10.(lguys.calc_r(snap))

	limits = (-2, 1, -100, 100)
	fig = Figure()
	ax = Axis(fig[1,1],
		limits=limits,
		xlabel=L"\log r / \textrm{kpc}",
		ylabel=L"v_\textrm{rad} / \textrm{km\,s^{-1}}"
	)
	
	h = Arya.hist2d!(ax, logr, v_rad, weights=mass, bins=200, colorrange=(1e-10, 1), colorscale=log10,
	)
	vlines!(log10.(r_b_kpc))

	Colorbar(fig[1, 2], h, label="stellar mass density")
	fig
end

# ╔═╡ Cell order:
# ╠═377284f2-dcee-44d3-9a04-728605cea92a
# ╠═340ffbbe-17bd-11ef-35c6-63505bb128b7
# ╠═f0d2b68a-fae2-4486-a434-a8816e400e84
# ╠═8ed2e715-1095-4c99-9ecb-54fecbc27d7f
# ╠═436a5be3-f597-4fc4-80a8-dc5af302ad66
# ╠═29988108-b02c-418c-a720-5766f47c39ff
# ╠═f0d74eaa-81e9-4b04-9765-24a0935b1430
# ╠═1b5c00d2-9df6-4a9c-ae32-05abcbf0e41a
# ╠═84dc77f7-14b3-4a2e-a556-c025d7df0095
# ╠═5e12d306-a430-4b15-b3a7-d4806a5856cd
# ╠═ef3481f8-2505-4c04-a04a-29bdf34e9abb
# ╠═a80d9e83-6b11-4c55-94ec-294d4247af42
# ╠═172588cc-ae22-440e-8488-f508aaf7ce96
# ╠═973955ad-3214-42cf-831f-a782f1d2434a
# ╠═6feeaae2-cb01-46ad-ad1d-daaca1caf7ec
# ╠═8f9ee2a7-de43-43e0-8257-93ecc630044f
# ╠═396a53a3-de0f-4d97-9693-40f3757d66f9
# ╠═f9fe37ef-de81-4d69-9308-cda968851ed2
# ╠═6fb4b7e1-a22c-4ff8-bbe9-dbf5de5acd37
# ╠═77479cd4-513c-4603-9aa0-1acd964c403a
# ╠═e33245ad-7562-47af-94b9-9708af45b2a4
# ╠═7bc2c15c-33f7-43f3-a47f-ca39ffc22071
# ╠═8df71dd5-ba0b-4918-9dc0-791c6c3bad5f
# ╟─c6ce37f4-09f6-43c2-9ee3-286d9e6baf7f
# ╠═a1138e64-5fa5-4a0e-aeef-487ee78a7adc
# ╠═91daae57-94dc-4a8e-b981-75f1406e0768
# ╠═1866c280-89c3-4a71-9dbf-50b583360145
# ╟─f5a3ea2f-4f19-4b0e-af55-74902f2c6485
# ╠═01b5ad85-9f37-4d8b-a29d-e47526f112ec
# ╠═1f722acb-f7b9-4d6c-900e-11eae85e0708
# ╠═51dea031-015b-4506-879d-9245c122d608
# ╠═249351c0-6fb6-49ee-ab44-e464f34a1bbe
# ╠═f87e4610-34ac-49f9-9858-0b3ef72eef15
# ╠═c351295d-1fb9-4088-8a27-5c732924959e
# ╠═8e32cd34-e97e-4cb7-95ba-67f4ed0c9aed
# ╠═1ef7e0c1-56f3-4be7-a138-c088ff973ec6
# ╠═24636618-7310-4557-baac-01d3d5d076dd
# ╠═46f1792a-53a4-4540-b8f3-93747f071ed0
# ╠═1b528c3f-877b-40f5-bc05-7b1c96c8d5c9
# ╠═0215fe76-f9c5-4274-ae43-89960a0caaef
# ╠═7e588ae3-89f3-4b91-8963-f6bf4391a859
# ╠═5e17bf16-2e3d-4372-8283-c43b8ad2550d
# ╠═ff759ae0-35c3-460e-97df-f865e26a0f48
# ╠═cf6a7cbb-1034-4026-a3f3-1e854d2929e2
# ╠═2612e3a2-6e6e-494e-b140-720dd2db6ec2
# ╠═7c4c5136-17d5-4dc5-9e5c-25e2348d2a84
# ╠═494a0dda-0a31-4049-a089-d78437c209cc
# ╠═d2de33be-b037-4255-b9f4-da29abb23754
# ╠═b1b218ff-694b-48bc-b999-806330ad4308
# ╠═32852a2a-1742-4c37-9430-d1425e4f2521
# ╠═8ad01781-8b5d-4d57-a0b5-7a445fb09b5b
# ╠═0487a830-3152-4d43-89ed-7562dc2dea52
# ╠═d89d2f21-e6c2-4911-8dd3-480b8e51020b
# ╠═2ef43371-bfae-4a65-8fa9-d1ab5ade32f1
# ╠═4537d2a8-dc90-4106-81c5-4a230734b182
# ╠═1c242116-e66d-453b-ad62-b6a65cdbe284
# ╠═b3e68e32-c058-467d-b214-aab6a4cd1e19
# ╠═e5c84f04-6df0-4f17-977c-7d31cd0b37dc
# ╠═52abcdc2-4b9e-408f-9c1f-9dc0678c60ca
# ╠═15e303d8-6020-4cfd-ac01-419a253a4a0b
# ╠═6858dc1b-6f60-4ca1-b2a2-80b4ca345b8c
# ╠═9fed63d6-c139-4d28-b00c-37dc1b8dc004
# ╠═9d74ffbb-4c31-4062-9434-7755f53e4da0
# ╠═975a2008-cf02-4442-9ee9-0b1bbb20889d
# ╠═12a30334-899e-4061-b7d2-af8c2346721d
# ╠═6fe6deb4-ae44-4ca0-9617-95841fdaf791
# ╠═20e742e9-a909-4be0-822b-f4ca6015b8aa
# ╠═25ccc096-1d27-43ce-8fe1-5a3a9d7cd54e
# ╠═55211531-6eb7-4689-9ecd-74854b8ad884
# ╠═487c6fda-5cb3-4ef4-a988-da2400c96cd4
# ╟─e23e85f9-7667-4f6e-8af6-2516fa292e2b
# ╠═0dd09c1e-67c0-4f23-bd12-41cbef62e4de
# ╠═f348c1b6-5196-4546-ba22-72677611a420
# ╠═d32757d2-dc08-488c-9ddd-d3eefefa2db7
# ╠═194ff30a-31a7-44bc-ac54-722a629457fc
# ╠═fbd46bd2-79d7-460e-b0ab-0e34a68a1f0a
# ╟─ed097e29-d2dc-4fa7-94e5-483b380600cc
# ╟─76438e61-d227-41cf-b9ea-e658bc389772
# ╟─04629bcf-5c19-41eb-8903-72947c209cbf
# ╟─b9a4b2b7-be95-4ccb-ad74-9b761abfae8a
# ╠═750185c9-a317-4978-aa04-486e5bfb7a63
# ╠═02183d68-bf7a-4774-9b48-f50712eeb552
# ╠═9b75409d-55f5-47c3-ab63-8168d31d3d54
# ╠═6408828d-d570-4585-8fa8-24661857fb35
# ╠═df2fc50f-3a03-424e-ae88-02ae207a04ee
# ╠═b36b594c-23b4-4683-8a12-3fa1b4c8b0d9
# ╠═03cdff7d-6e86-45ab-9837-9029293639fa
# ╠═cdbed68e-0da2-4648-a5d2-61b5b07fb4a2
# ╠═f3db0a1d-78bc-4ac4-84b6-d53f27a28795
# ╠═b4f7cbb0-fb5b-4dc5-b66e-679a8e5b630d
# ╠═d664ab12-a2c1-4531-a4ff-250ffa3ce9eb
# ╠═d42795d0-bd69-4c2c-be5b-e27e85199ee3
# ╠═a0391689-66a2-473f-9704-e12a3d033d13
# ╠═44ab0a25-ab4c-4e90-8619-2a068a285755
# ╠═c2ccf9de-e3cd-4950-9a4d-6f425d261ccb
# ╠═89ea28dc-2677-4e6f-9066-d5ffff3b474c
# ╠═6e34b91c-c336-4538-a961-60833d37f070
# ╠═17c880ba-d6de-4d9c-b44b-d5ecdb4da226
# ╠═57d5decd-8259-4ec2-87ac-44d28625cd7b
# ╠═98bcfb2a-8eb5-4c7c-8dea-46d000a6d245
# ╠═0c26f965-0381-4f79-a6ce-0772ae922b3f
