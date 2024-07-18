### A Pluto.jl notebook ###
# v0.19.43

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

# ╔═╡ d401ec4b-048e-4aae-85a8-f7f0d8e44a79
using LilGuys

# ╔═╡ cf6a7cbb-1034-4026-a3f3-1e854d2929e2
using FITSIO

# ╔═╡ f0d2b68a-fae2-4486-a434-a8816e400e84
import TOML

# ╔═╡ cb6a58a6-9ba9-44b5-95a6-062965c13259
models_dir = "/arc7/home/dboyea/sculptor"

# ╔═╡ 0a73bf88-3f46-4864-97f5-41705ea6913d
model_dir = "/arc7/home/dboyea/sculptor/orbits/V50_r0.5"

# ╔═╡ 29988108-b02c-418c-a720-5766f47c39ff
starsname = "exp2d_rs0.07"

# ╔═╡ d7f5a3ed-ae4a-4ea3-b776-00dba6506a88
r_scale = 1

# ╔═╡ f0d74eaa-81e9-4b04-9765-24a0935b1430
starsfile = "/arc7/home/dboyea/sculptor/isolation/1e6/halos/V50_r0.5/stars/$(starsname)_stars.hdf5"

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

# ╔═╡ b3a16249-b8d9-4a6b-9294-cd654a17dc17
md"""
# Inputs
"""

# ╔═╡ 1b5c00d2-9df6-4a9c-ae32-05abcbf0e41a
paramsfile = "/astro/dboyea/sculptor/isolation/1e6/halos/V50_r0.5/stars/$starsname.toml"

# ╔═╡ 172588cc-ae22-440e-8488-f508aaf7ce96
rel_p_cut = 1e-20

# ╔═╡ ef3481f8-2505-4c04-a04a-29bdf34e9abb
outfile = "stars/$(starsname)_today.fits"

# ╔═╡ 973955ad-3214-42cf-831f-a782f1d2434a
idx_i = 1 

# ╔═╡ 5e12d306-a430-4b15-b3a7-d4806a5856cd
name = splitext(basename(starsfile))[1]

# ╔═╡ a80d9e83-6b11-4c55-94ec-294d4247af42
outfile_i = "stars/$(starsname)_i_today.fits"

# ╔═╡ 9c42eb0a-029d-46f7-afb0-de03f82c5889
obs_today_filename = "../../models/sculptor/mc_orbits/orbit1.toml"

# ╔═╡ 396cd0a8-1d73-44dd-89db-3243fb9e8ac4
md"""
# File loading
"""

# ╔═╡ 436a5be3-f597-4fc4-80a8-dc5af302ad66
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ 8f9ee2a7-de43-43e0-8257-93ecc630044f
idx_f = orbit_props["idx_f"]

# ╔═╡ 84dc77f7-14b3-4a2e-a556-c025d7df0095
params = TOML.parsefile(paramsfile)

# ╔═╡ 2612e3a2-6e6e-494e-b140-720dd2db6ec2
obs_today_file = TOML.parsefile(obs_today_filename)

# ╔═╡ 07defcc4-8dca-4caf-a34a-a341818d80ad
length(probabilities)

# ╔═╡ 6fb4b7e1-a22c-4ff8-bbe9-dbf5de5acd37
begin 
	out =  lguys.Output(model_dir * "/out", weights=probabilities)
	out
end

# ╔═╡ 396a53a3-de0f-4d97-9693-40f3757d66f9
snap_i = out[idx_i]

# ╔═╡ 44cbb2ce-5f43-4bc9-a4c4-c0b9df692cd2
length(snap_i.masses)

# ╔═╡ 6feeaae2-cb01-46ad-ad1d-daaca1caf7ec
snap_f = out[idx_f]

# ╔═╡ 5ee4f95d-0587-44ab-b543-9b7039d545e6
md"""
# Plots
"""

# ╔═╡ 1fd7b586-75ad-42a8-b752-ea82b290cb47
out.weights

# ╔═╡ 77479cd4-513c-4603-9aa0-1acd964c403a
let
	fig = Figure(size=(700, 700))
	r_max = 100

	x = lguys.get_x(snap_f)
	y = lguys.get_y(snap_f)
	z = lguys.get_z(snap_f)
	ps = snap_f.weights
	bins = LinRange(-r_max, r_max, 300)
	
	kwargs = (colorscale=log10, colorrange=(1e-10, nothing), weights=ps, bins=bins)

	
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
	fig
end

# ╔═╡ 7bc2c15c-33f7-43f3-a47f-ca39ffc22071
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc", title="initial")

	bin_range = LinRange(-2, 2, 100)
	colorrange =(1e-5, nothing)

	bins = (out.x_cen[1, idx_i]  .+ bin_range,  out.x_cen[2, idx_i]  .+ bin_range)
		

	Arya.hist2d!(ax, snap_i.positions[1, :], snap_i.positions[2, :], weights=snap_i.weights, bins = bins, colorscale=log10, colorrange=colorrange)

	ax2 = Axis(fig[1,2], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc",
	title="final")


	bins = (out.x_cen[1, idx_f]  .+ bin_range,  out.x_cen[2, idx_f]  .+ bin_range)

	Arya.hist2d!(ax2, snap_f.positions[1, :], snap_f.positions[2, :], weights=snap_f.weights, bins = bins, colorscale=log10, colorrange=colorrange)
	
	fig
end

# ╔═╡ 8df71dd5-ba0b-4918-9dc0-791c6c3bad5f
let 
	fig = Figure()
	ax = Axis(fig[1,1], yscale=log10, limits=(nothing, (1e-8, 1)),
		xlabel=L"\epsilon", ylabel="count")
	stephist!(lguys.calc_ϵ(snap_i), weights=snap_i.weights)
	es = lguys.calc_ϵ(snap_f)
	filt = es .> 0
	es = es[filt]
	stephist!(es, weights = snap_f.weights[filt])

	fig

end

# ╔═╡ c6ce37f4-09f6-43c2-9ee3-286d9e6baf7f
md"""
# 3D density profiles
"""

# ╔═╡ a1138e64-5fa5-4a0e-aeef-487ee78a7adc
function plot_ρ_s!(snap; bins=100, kwargs...)
	rs = lguys.calc_r(snap.positions, snap.x_cen)
	r, ρ = lguys.calc_ρ_hist(rs, bins, weights=snap.weights)
	lines!(log10.(lguys.midpoints(r)), log10.(ρ); kwargs...)
end

# ╔═╡ 91daae57-94dc-4a8e-b981-75f1406e0768
begin
	profile_class = getproperty(lguys, Symbol(params["profile"]))
	params["profile_kwargs"]["R_s"] *= r_scale
	profile = profile_class(;lguys.dict_to_tuple(params["profile_kwargs"])...)
	ρ_s(r) = lguys.calc_ρ(profile, r)
end

# ╔═╡ a0391689-66a2-473f-9704-e12a3d033d13
import LinearAlgebra: dot

# ╔═╡ 44ab0a25-ab4c-4e90-8619-2a068a285755
""" 
	calc_v_rad(snap)

returns the radial velocities relative to the snapshot centre in code units
"""
function calc_v_rad(snap)
	x_vec = snap.positions .- snap.x_cen
	v_vec = snap.velocities .- snap.v_cen

	# normalize
	x_hat = x_vec ./ lguys.calc_r(x_vec)'

	# dot product
	v_rad = sum(x_hat .* v_vec, dims=1)

	# matrix -> vector
	v_rad = dropdims(v_rad, dims=1)
	
	return v_rad 
end

# ╔═╡ 6e34b91c-c336-4538-a961-60833d37f070
function v_rad_hist(snap, bins=40)

	mass = snap.weights
	v_rad = calc_v_rad(snap)
	logr = log10.(lguys.calc_r(snap))
	h1 = Arya.histogram(logr, bins, weights=v_rad .* mass, normalization=:none)
	h2 = Arya.histogram(logr, bins, weights=mass, normalization=:none)

	x_bins = h1.bins
	v_bins = h1.values
	counts = h2.values

	return x_bins, v_bins ./ counts
end

# ╔═╡ 227a4b71-afbd-4121-930b-696d06ccc9ba
md"""
double checking the velocity radial 3D calculation. Blue arrows should all point inward and red outward. Tiny slice in the x-y vx-vy plane...
"""

# ╔═╡ 253e43df-58fc-4dee-b1c4-35e273499ab7
let
	vs = calc_v_rad(snap_f)

	filt = snap_f.weights .> 0.01 * maximum(snap_f.weights)
	x = lguys.get_x(snap_f)[filt]
	y = lguys.get_y(snap_f)[filt]
	w = snap_f.weights[filt]
	scatter(x, y, color=vs[filt])
end

# ╔═╡ 7a94107d-2a45-4458-8d26-5cf836501a1e
let
	snap = snap_f
	xc = snap.x_cen
	vc = snap.v_cen

	dx = 0.1

	filt = lguys.get_z(snap) .< 1e-2
	filt .&= lguys.get_v_z(snap) .< 1e-2
	filt .&= snap.weights .> 0.1*maximum(snap.weights)
	println(sum(filt))

	snap = snap[filt]
	vs = calc_v_rad(snap)

	
	fig, ax = FigAxis(
		aspect=DataAspect(),
		limits=dx .* (-1, 1, -1, 1),
		xlabel="x/kpc",
		ylabel="y/kpc"
	)

	vx = lguys.get_v_x(snap) .- vc[1]
	vy = lguys.get_v_y(snap) .- vc[2]

	x = lguys.get_x(snap) .- xc[1]
	y = lguys.get_y(snap) .- xc[2]
	h = arrows!(x, y, vx, vy, color=vs * V2KMS, colorrange=(-1, 1), colormap=:bluesreds,
	label="3D radial velocity km/s")

	Colorbar(fig[1, 2], h)

	fig
end

# ╔═╡ 17b8f17b-0801-45ca-a86b-bba1d78f9ecd
lguys.calc_r(snap_f)

# ╔═╡ 51ba2a8e-6f6f-43bd-ac5f-5d238bd41165
lguys.calc_r(snap_f, snap_f.x_cen)

# ╔═╡ 56be7b5a-3e46-4162-94cb-3f5783efd183
snap_f.x_cen

# ╔═╡ f5a3ea2f-4f19-4b0e-af55-74902f2c6485
md"""
# Sky projection
"""

# ╔═╡ 12a30334-899e-4061-b7d2-af8c2346721d
let 
	x = log10.(out.weights[out.weights .> 0])
	h = Arya.histogram(x, normalization=:none)

	barplot(h)
end

# ╔═╡ 51dea031-015b-4506-879d-9245c122d608
snap_cen = lguys.Snapshot(out.x_cen, out.v_cen, ones(size(out.x_cen, 2)))

# ╔═╡ f87e4610-34ac-49f9-9858-0b3ef72eef15
cen = [out.x_cen[:, idx_f]; out.v_cen[:, idx_f]]

# ╔═╡ 24636618-7310-4557-baac-01d3d5d076dd
cen

# ╔═╡ 46f1792a-53a4-4540-b8f3-93747f071ed0
obs_c_galcen = lguys.Galactocentric(x=cen[1], y=cen[2], z=cen[3], 
	v_x=cen[4]*V2KMS, v_y=cen[5]*V2KMS, v_z=cen[6]*V2KMS)

# ╔═╡ 7e588ae3-89f3-4b91-8963-f6bf4391a859
function save_obs(obs_df, outfile)
	rm(outfile, force=true)
	FITS(outfile, "w") do f
		df = Dict(String(name) => obs_df[:, name] for name in names(obs_df))
		write(f, df)
	
		println("written to $outfile")

		df
	end
end

# ╔═╡ 77918cc8-bfd7-4a84-a1d2-b7d5c721f5ba
mkpath(joinpath(model_dir,"stars") )

# ╔═╡ 7c4c5136-17d5-4dc5-9e5c-25e2348d2a84
obs_today_icrs = lguys.ICRS(;
	ra=obs_today_file["ra"], dec=obs_today_file["dec"],
	distance=obs_today_file["distance"],
	pmra=obs_today_file["pmra"],
	pmdec=obs_today_file["pmdec"],
	radial_velocity=obs_today_file["radial_velocity"],
)

# ╔═╡ 494a0dda-0a31-4049-a089-d78437c209cc
obs_today_err = lguys.ICRS(;
	ra=0, dec=0,
	distance=obs_today_file["distance_err"],
	pmra=obs_today_file["pmra_err"],
	pmdec=obs_today_file["pmdec_err"],
	radial_velocity=obs_today_file["radial_velocity_err"],
)

# ╔═╡ d2de33be-b037-4255-b9f4-da29abb23754
frame = lguys.GSR

# ╔═╡ 01b5ad85-9f37-4d8b-a29d-e47526f112ec
function make_sample(snap; 
	cen=nothing, rel_p_cut=rel_p_cut, r_max=Inf,
	frame=frame
)
	filt = snap.weights .> rel_p_cut * maximum(snap.weights)
	
	obs_df = lguys.to_sky(snap[filt], SkyFrame=frame, add_centre=true)

	obs_df[!, :xi], obs_df[!, :eta] = lguys.to_tangent(obs_df.ra, obs_df.dec, obs_df.ra[1], obs_df.dec[1])

	obs_df[!, :r_ell] = 60lguys.calc_r_ell(obs_df.xi, obs_df.eta, 0, 0)

	return obs_df[:, Not(:frame)]
end

# ╔═╡ 1f722acb-f7b9-4d6c-900e-11eae85e0708
obs_df = make_sample(snap_f)

# ╔═╡ 0215fe76-f9c5-4274-ae43-89960a0caaef
obs_c = obs_df[1, :]

# ╔═╡ 5e17bf16-2e3d-4372-8283-c43b8ad2550d
lguys.write_fits(model_dir * "/" * outfile, obs_df)

# ╔═╡ 249351c0-6fb6-49ee-ab44-e464f34a1bbe
sky_orbit = lguys.to_sky(snap_cen, SkyFrame=frame)

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

	snap_i_shifted = deepcopy(snap_i)
	x_c = [obs_gc.x, obs_gc.y, obs_gc.z]
	v_c = [obs_gc.v_x, obs_gc.v_y, obs_gc.v_z] / V2KMS

	snap_i_shifted.positions .+= x_c .- snap_i.x_cen

	snap_i_shifted.velocities .+= v_c .- snap_i.v_cen
	snap_i_shifted.x_cen = x_c
	snap_i_shifted.v_cen = v_c

	obs_df_i = make_sample(snap_i_shifted)
end

# ╔═╡ ff759ae0-35c3-460e-97df-f865e26a0f48
lguys.write_fits(model_dir * "/" * outfile_i, obs_df_i)

# ╔═╡ 9f9b525d-f6f6-4fc0-b5b9-036662fe8ba8
md"""
# Plots of simulated sample
"""

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

	
	h = Arya.hist2d!(ax, x, y, weights=obs_df_i.weights, bins=100, colorscale=log10, colorrange=(1e-10, nothing))

	Colorbar(fig[1, 2], h)
	fig
end

# ╔═╡ d7fece88-3327-4435-ab61-b45ff62b3b2e
function mean_2d(obs_df, values; bins=100, centre=false, limits=nothing)
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

	return h_vel

end

# ╔═╡ 7f9d1703-8bef-4eb4-9864-50724ff168af
obs_c

# ╔═╡ e904f104-2d01-45f0-a6f1-2040131d8780
function ra_dec_axis(ddeg=5; kwargs...)
	fig = Figure(;kwargs...)
	
	dy = ddeg
	dx = dy * 1/cosd(obs_c.dec)
	limits = (obs_c.ra .+ (-dx, dx), obs_c.dec .+ (-dy, dy))

	ax = Axis(fig[1,1],
		xlabel="RA / degrees",
		ylabel="dec / degrees",
		limits=limits,
		aspect = 1,
		xgridvisible=false,
		ygridvisible=false
		
	)

	return fig, ax
end

# ╔═╡ 6ebfff07-c43f-4d4d-8604-9fd4f1de5d25
function xi_eta_axis(dx=10, dy=5; kwargs...)
	fig = Figure(;kwargs...)
	

	limits = ((-dx, dx), (-dy, dy))

	ax = Axis(fig[1,1],
		xlabel=L"$\xi'$ / degrees",
		ylabel=L"$\eta'$ / degrees",
		limits=limits,
		aspect = DataAspect(),
		xgridvisible=false,
		ygridvisible=false
		
	)

	return fig, ax
end

# ╔═╡ 8ad01781-8b5d-4d57-a0b5-7a445fb09b5b
let 
	fig, ax = ra_dec_axis()

	bins = 100
	limits = ax.limits.val
	x = obs_df.ra
	y = obs_df.dec


	hi = Arya.histogram2d(x, y, bins, weights=obs_df.weights, limits=limits)
	areas = diff(hi.xbins) .* (diff(hi.ybins)')
	hi.values ./= areas
	
		
	h = heatmap!(hi, colorscale=log10, colorrange=(1e-12, maximum(hi.values)))
	errscatter!([obs_today.ra], [obs_today.dec], color=COLORS[3], size=10)

	idx = idx_f - 20: idx_f + 20
	lines!(sky_orbit.ra[idx], sky_orbit.dec[idx])
	
	Colorbar(fig[1, 2], h,
		label="stellar density"
	)
	fig
end

# ╔═╡ edf68b42-4fe9-4e14-b7ed-739e89a1541a
let
	dr = 10
	fig, ax = xi_eta_axis(dr, dr)

	bins = LinRange(-dr, dr, 100)

	h = Arya.histogram2d(obs_df.xi, obs_df.eta, bins, weights=obs_df.weights)

	p = heatmap!(h.xbins, h.ybins, h.values, colorscale=log10, colorrange=(1e-20, maximum(h.values)))

	Colorbar(fig[1, 2], p)

	fig
end

# ╔═╡ 7fa19be2-4014-4df0-821e-4c509eca4f28
let
	dr = 10
	fig, ax = xi_eta_axis(dr, dr)

	bins = LinRange(-dr, dr, 100)

	h = mean_2d(obs_df,  obs_df.radial_velocity, bins=bins)

	p = heatmap!(h.xbins, h.ybins, h.values, 
	colormap=:redsblues,  colorrange=(50, 80))

	Colorbar(fig[1, 2], p)

	fig
end

# ╔═╡ 96307998-07a0-45bf-bf10-cd14bfcfe20a
let
	dr = 1
	fig, ax = xi_eta_axis(dr, dr)

	bins = LinRange(-dr, dr, 100)

	h = mean_2d(obs_df,  obs_df.radial_velocity, bins=bins)

	p = heatmap!(h.xbins, h.ybins, h.values, 
	colormap=:redsblues,  colorrange=(50, 80))

	Colorbar(fig[1, 2], p)

	fig
end

# ╔═╡ 35c52910-089c-4507-956a-2b0649507495
filt_cen = obs_df.r_ell .< 2

# ╔═╡ 2ef43371-bfae-4a65-8fa9-d1ab5ade32f1
let
	println("testing pm if works")

	i1 = idx_f
	i2 = idx_f - 1
	dt = (out.times[i1] - out.times[i2]) * 1e9 * T2GYR # years
	
 	ddec = sky_orbit.dec[i1] - sky_orbit.dec[i2] 
	dra = sky_orbit.ra[i1] - sky_orbit.ra[i2] 

	deg_to_mas = 60^2 * 1e3
	
	pmdec = ddec * deg_to_mas / dt 
	pmra = dra * deg_to_mas / dt  * cos(deg2rad(sky_orbit.dec[i1]))

	println(pmra, ", ", pmdec)
	println(sky_orbit.pmra[i1], ", ", sky_orbit.pmdec[i1])
end

# ╔═╡ 1c242116-e66d-453b-ad62-b6a65cdbe284
sky_orbit[idx_f, :]

# ╔═╡ b3e68e32-c058-467d-b214-aab6a4cd1e19
r_cut = 120 # arcmin

# ╔═╡ 9fed63d6-c139-4d28-b00c-37dc1b8dc004
let 
	fig, ax = FigAxis(aspect=1,
		xlabel="RA / degrees",
		ylabel="Dec / degrees",
		xgridvisible=false,
		ygridvisible=false
	)
		
	scatter!(obs_df.ra[obs_df.r_ell .< r_cut], obs_df.dec[obs_df.r_ell .< r_cut], 
				alpha=0.1, color=:black, markersize=3)

	fig
end

# ╔═╡ e93f2ee2-5fb3-43cb-ae34-58a0941534f1
hist(obs_df_i.r_ell)

# ╔═╡ b7cd6c89-4534-4222-bb27-1ff511692ee1
let
	dr = 0.2
	r_max = r_cut
	
	limits = (obs_df_i.pmra[1] .+ (-dr, dr), obs_df_i.pmdec[1] .+ (-dr, dr))
	
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		title="initial",
	limits=limits,
	aspect=DataAspect())

	filt = obs_df_i.r_ell .< r_max

	x = obs_df_i.pmra[filt]
	y = obs_df_i.pmdec[filt]
	
	h = Arya.hist2d!(ax, x, y, weights=obs_df_i.weights[filt], bins=100, 
		colorscale=log10, colorrange=(1e-15, nothing))

	Colorbar(fig[1, 2], h, label="stellar density")
	
	fig
end

# ╔═╡ 9d74ffbb-4c31-4062-9434-7755f53e4da0
let
	dr = 0.1
	r_max = r_cut
	
	limits = (obs_df.pmra[1] .+ (-dr, dr), obs_df.pmdec[1] .+ (-dr, dr))
	
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		limits=limits,
		aspect=DataAspect(),
		title="today",
	)

	filt = obs_df.r_ell .< r_max

	x = obs_df.pmra[filt]
	y = obs_df.pmdec[filt]
	
	h = Arya.hist2d!(ax, x, y, weights=obs_df.weights[filt], bins=100, 
		colorscale=log10, colorrange=(1e-10, nothing))
	
	errscatter!([obs_today.pmra], [obs_today.pmdec], xerr=[obs_today_err.pmra], yerr=[obs_today_err.pmra], color=COLORS[3])
	
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
	bins = 80
	
	limits = (obs_c.pmra .+ (-dr, dr), nothing)
	
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel="density",
	limits=limits)

	filt = obs_df.r_ell .< r_cut

	x = obs_df.pmra[filt]
	y = obs_df.pmdec[filt]
	
	h = Arya.histogram(x, bins, weights=obs_df.weights[filt], normalization=:pdf)

	
	barplot!(h)
	
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
		# limits=limits,
		title="initial",
	)


	
	x = obs_df_i.distance
	y = obs_df_i.radial_velocity
	
	h = Arya.hist2d!(ax, x, y, weights=obs_df_i.weights, bins=100,
		colorscale=log10, colorrange=(1e-10, nothing)
	)

	
	Colorbar(fig[1, 2], h)
	
	fig
end

# ╔═╡ 4eae6b35-e71d-47da-bde9-d67b30ce143a
hist(obs_df.radial_velocity, weights=obs_df.weights, bins=100)

# ╔═╡ 6fe6deb4-ae44-4ca0-9617-95841fdaf791
let
	
	fig = Figure()
	dx = 10
	dy = 40
	r_max = 2
	
	limits = (sky_orbit.distance[idx_f] .+ (-dx, dx), sky_orbit.radial_velocity[idx_f] .+ (-dy, dy))
	
	ax = Axis(fig[1,1],
		xlabel="distance / kpc",
		ylabel = L"$\tilde{v}_\textrm{rad}$ / km s$^{-1}$",
		limits=limits,
		title="today",
	)



	filt = obs_df.xi .^ 2 .+ obs_df.eta .^ 2 .< r_max ^ 2
	x = obs_df.distance[filt]
	y = obs_df.radial_velocity[filt]
	
	h = Arya.hist2d!(ax, x, y, weights=obs_df.weights[filt], bins=100,
		colorscale=log10, colorrange=(1e-10, nothing)
	)
	
	errscatter!([obs_today.distance], [obs_today.radial_velocity], xerr=[obs_today_err.distance], yerr=[obs_today_err.radial_velocity], color=COLORS[3])

	scatter!([sky_orbit.distance[idx_f]], 
		[sky_orbit.radial_velocity[idx_f]], markersize=10)
	
	idx = idx_f - 20: idx_f + 20
	lines!(sky_orbit.distance[idx], sky_orbit.radial_velocity[idx])

	
	Colorbar(fig[1, 2], h)

	fig
end

# ╔═╡ 55211531-6eb7-4689-9ecd-74854b8ad884
hist(obs_df_i.radial_velocity, weights=obs_df_i.weights)

# ╔═╡ e23e85f9-7667-4f6e-8af6-2516fa292e2b
md"""
# Velocity profile & Break Radius Calculation
"""

# ╔═╡ 0dd09c1e-67c0-4f23-bd12-41cbef62e4de
import StatsBase: weights, std

# ╔═╡ 5355bcc3-2494-473d-9da7-fc38930b8ee7
let
	fig, ax = FigAxis(
		limits = (-1, 2.5, 30, 100),
		xlabel = "log r",
		ylabel="radial velocity"
		
	)

	hist2d!(log10.(obs_df.r_ell), obs_df.radial_velocity, weights=obs_df.weights,
		colorscale=log10,
		colorrange=(1e-8, 0.3),
		bins=100
	)
	
	fig
end

# ╔═╡ 82a491cd-e351-4715-90be-5a2d17e960ce
obs_df.r_ell

# ╔═╡ d32757d2-dc08-488c-9ddd-d3eefefa2db7
"""
Given a set of radii and velocity
"""
function calc_σv(r_ell, rv, mass;  r_max = 30)
	filt = r_ell .< r_max
	vel = rv[filt]

	return lguys.std(vel, weights(mass)[filt])
end

# ╔═╡ b4f7cbb0-fb5b-4dc5-b66e-679a8e5b630d
function calc_σv(snap::lguys.Snapshot; r_max=1)
	rs = lguys.calc_r(snap)
	filt = rs .< r_max
	vs = lguys.calc_v(snap[filt])
	masses = snap.weights[filt]
	σ = lguys.std(vs, weights(masses), mean=0) # zero mean

	return σ * V2KMS / sqrt(3)
end

# ╔═╡ ca228ca7-ff7c-4f73-8fa7-28707c61d8e5
let

	fig = Figure()
	ax = Axis(fig[2, 1],
		xlabel = "log r / arcmin",
		ylabel=L"$\sigma_\textrm{v, los}$ / km s$^{-1}$",
		limits=((0.8, 2.5), (0, 20))
	)

	x = log10.(obs_df.r_ell)
	prob = obs_df.weights
	v = obs_df.radial_velocity

	filt = @. !isnan(x)
	filt .&= x .> -0.5
	filt .&= x .< 3

	x = x[filt]
	prob = prob[filt]
	v = v[filt]

	bins = 50
	
	r_bins = Arya.make_bins(x, (-0.5, 2.3), Arya.bins_equal_number, n=bins)

	println(r_bins)
	
	σs = Vector{Float64}(undef, bins)
	err = Vector{Float64}(undef, bins)
	Ns = Vector{Float64}(undef, bins)
	
	for i in 1:bins
		filt = r_bins[i] .<= x .< r_bins[i+1]
		σs[i] = calc_σv(x[filt], v[filt], prob[filt], r_max=Inf)
		N = sum(filt)
		if N > 1
			err[i] = σs[i] * sqrt(2/(N - 1))
		else
			err[i] = NaN
		end
		Ns[i] = N
	end

	errscatter!(lguys.midpoints(r_bins), σs, yerr=err)

	# ax2 = Axis(fig[1, 1], ylabel="count / bin")
	# scatter!(ax2, lguys.midpoints(r_bins), Ns)
	# hidexdecorations!(ax2, grid=false)

	# linkxaxes!(ax, ax2)

	fig
end

# ╔═╡ 422839f0-6da4-46b9-8689-2dd13b03188b
function calc_σvx(snap::lguys.Snapshot; r_max=1)
	rs = lguys.calc_r(snap)
	filt = rs .< r_max
	vs = snap[filt].velocities[1, :] .- snap.v_cen[1]
	masses = snap.weights[filt]
	
	σ = lguys.std(vs, weights(masses))
	return σ * V2KMS
end

# ╔═╡ 194ff30a-31a7-44bc-ac54-722a629457fc
σv = calc_σv(obs_df.r_ell, obs_df.radial_velocity, obs_df.weights, r_max=r_cut)

# ╔═╡ dcb868c4-2619-4437-af8d-37842bfb6e63
σv_i = calc_σv(obs_df_i.r_ell, obs_df_i.radial_velocity, obs_df_i.weights, r_max=r_cut * 60)

# ╔═╡ d664ab12-a2c1-4531-a4ff-250ffa3ce9eb
calc_σv(snap_i, r_max=6) 

# ╔═╡ 08b66f99-f81a-4495-a933-9291e986373a
calc_σvx(snap_i, r_max=6) 

# ╔═╡ cfccf1a1-22ac-4fb8-8b91-c7af98ad3c4d
lguys.arcmin_to_kpc(240, 86)

# ╔═╡ fbd46bd2-79d7-460e-b0ab-0e34a68a1f0a
gaussian(x, μ, σ) = 1/sqrt(2π)* 1/σ * exp(-(x-μ)^2/(2σ^2))

# ╔═╡ 04629bcf-5c19-41eb-8903-72947c209cbf
let
	snap = snap_f
	
	mass = snap_f.weights
	v_rad = calc_v_rad(snap)
	logr = log10.(lguys.calc_r(snap))
	filt = logr .< 1

	h = Arya.histogram(v_rad[filt] * V2KMS, 
		weights=mass[filt], normalization=:pdf, limits=(-20, 20))


	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"v_\textrm{rad} = \vec{v} \cdot \hat{r} \  / \ \textrm{km\,s^{-1}}",
		ylabel="stellar density"
	)
	
	scatter!(lguys.midpoints(h.bins), h.values)

	σ = calc_σv(snap)
	x_model = LinRange(-20, 20, 100)
	y_model = gaussian.(x_model, 0, σ*√3/2)
	lines!(x_model, y_model)
	fig
	
end

# ╔═╡ ed097e29-d2dc-4fa7-94e5-483b380600cc
let
	
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="RV (km/s)",
		ylabel="density",
		title="within $r_cut arcmin"
	)
	
	filt = obs_df.r_ell .< r_cut
	
	rv = obs_df.radial_velocity[filt]
	mass = obs_df.weights[filt]

	μ = Arya.mean(rv, Arya.sb.weights(mass))
	
	h = Arya.histogram(rv, weights=mass, normalization=:pdf)

	scatter!(lguys.midpoints(h.bins), h.values)

	x_model = LinRange(μ - 3σv, μ+3σv, 100)
	y_model = gaussian.(x_model, μ, σv)
	lines!(x_model, y_model, 
		color=COLORS[2], 
		label="gaussian (σ = $(round(σv, digits=2)) km / s)"
	)

	axislegend()
	fig
end

# ╔═╡ c2ccf9de-e3cd-4950-9a4d-6f425d261ccb
function calc_v_tan(snap)
	v = lguys.calc_v(snap)
	v_rad = calc_v_rad(snap)
	v_tan = @. sqrt(v^2 - v_rad^2)
	return vec(v_tan)
end

# ╔═╡ 76438e61-d227-41cf-b9ea-e658bc389772
let
	snap = snap_i
	
	mass = snap.weights
	v_rad = (snap.velocities[1, :] .- snap.v_cen[1]) * V2KMS
	logr = log10.(lguys.calc_r(snap))
	filt = logr .< 1

	h = Arya.histogram(v_rad[filt], 
		weights=mass[filt], normalization=:pdf, limits=(-20, 20))


	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"v_{x} \quad [\textrm{km\,s^{-1}}]",
		ylabel="stellar density"
	)
	
	scatter!(lguys.midpoints(h.bins), h.values)

	x_model = LinRange(-20, 20, 100)
	y_model = gaussian.(x_model, 0, σv)
	lines!(x_model, y_model, 
		color=COLORS[2], label="gaussian (σ = $(round(σv, digits=2)) km / s)"
	)

	axislegend()
	fig
	
end

# ╔═╡ b9a4b2b7-be95-4ccb-ad74-9b761abfae8a
@doc raw"""
	calc_rb(σ, delta_t)

Given a radial velocity dispersion σ (in km/s), and the time since the last pericentre (in Gyr), returns the break radius (in kpc).

See peñarrubia et al. (200X) for the original equation.
```math
r_{\rm break} = C\,\sigma_{v}\,\Delta t
```
where $C=0.55$ is an empirical fit
"""
function calc_rb(σ, delta_t)
	kpc_per_Gyr_per_kms = 1.0227
	return 0.55 * σ * (delta_t) * kpc_per_Gyr_per_kms
end

# ╔═╡ 54d0ee8e-52d6-4b8b-84a9-1ddc66659137
orbit_props["t_last_peri"]

# ╔═╡ 750185c9-a317-4978-aa04-486e5bfb7a63
r_b_kpc = calc_rb(σv, orbit_props["t_last_peri"]) # kpc

# ╔═╡ 1866c280-89c3-4a71-9dbf-50b583360145
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log r / \textrm{kpc}", ylabel = L"\log \rho_\star", 
		limits=((-1.9, 1), (-12, 2)))

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

# ╔═╡ 0c26f965-0381-4f79-a6ce-0772ae922b3f
let
	snap = snap_f
	
	mass = snap.weights
	v_rad = calc_v_rad(snap) 
	logr = log10.(lguys.calc_r(snap))

	limits = (-2, 3, -100, 100)
	fig = Figure()
	ax = Axis(fig[1,1],
		limits=limits,
		xlabel=L"\log r / \textrm{kpc}",
		ylabel=L"v_\textrm{rad} / \textrm{km\,s^{-1}}"
	)
	
	h = Arya.hist2d!(ax, logr, v_rad * V2KMS, weights=mass, bins=200, colorrange=(1e-10, 1), colorscale=log10,
	)
	vlines!(log10.(r_b_kpc))

	Colorbar(fig[1, 2], h, label="stellar mass density")
	fig
end

# ╔═╡ 57d5decd-8259-4ec2-87ac-44d28625cd7b
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="log r",
		ylabel="mean 3D radial velocity (km / s)",
		limits=((nothing, 2), (-40, 50))
	)
	
	x, y= v_rad_hist(snap_f, 50)
	scatter!(lguys.midpoints(x), y * V2KMS)

	vlines!(log10.(r_b_kpc), linestyle=:dash, color=:black)
	fig
end

# ╔═╡ 02183d68-bf7a-4774-9b48-f50712eeb552
r_b_arcmin = lguys.kpc_to_arcmin(r_b_kpc, obs_today.distance)

# ╔═╡ 9b75409d-55f5-47c3-ab63-8168d31d3d54
md"""
# Evolutionary Properties
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
		filt[filt] .= rs .< radius

		ps = snap.weights[filt]
		M[i] = sum(ps)
	end

	return M
end

# ╔═╡ b36b594c-23b4-4683-8a12-3fa1b4c8b0d9
r_h = 0.1

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

# ╔═╡ d42795d0-bd69-4c2c-be5b-e27e85199ee3
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
	vs = [calc_σv(out[i]) for i in idx]

	scatter!(out.times[idx] * T2GYR, vs)
	
	fig
end

# ╔═╡ Cell order:
# ╠═377284f2-dcee-44d3-9a04-728605cea92a
# ╠═340ffbbe-17bd-11ef-35c6-63505bb128b7
# ╠═d401ec4b-048e-4aae-85a8-f7f0d8e44a79
# ╠═f0d2b68a-fae2-4486-a434-a8816e400e84
# ╟─b3a16249-b8d9-4a6b-9294-cd654a17dc17
# ╠═cb6a58a6-9ba9-44b5-95a6-062965c13259
# ╠═0a73bf88-3f46-4864-97f5-41705ea6913d
# ╠═29988108-b02c-418c-a720-5766f47c39ff
# ╠═d7f5a3ed-ae4a-4ea3-b776-00dba6506a88
# ╠═f0d74eaa-81e9-4b04-9765-24a0935b1430
# ╠═1b5c00d2-9df6-4a9c-ae32-05abcbf0e41a
# ╠═172588cc-ae22-440e-8488-f508aaf7ce96
# ╠═ef3481f8-2505-4c04-a04a-29bdf34e9abb
# ╠═973955ad-3214-42cf-831f-a782f1d2434a
# ╠═8f9ee2a7-de43-43e0-8257-93ecc630044f
# ╠═5e12d306-a430-4b15-b3a7-d4806a5856cd
# ╠═a80d9e83-6b11-4c55-94ec-294d4247af42
# ╠═9c42eb0a-029d-46f7-afb0-de03f82c5889
# ╟─396cd0a8-1d73-44dd-89db-3243fb9e8ac4
# ╠═436a5be3-f597-4fc4-80a8-dc5af302ad66
# ╠═84dc77f7-14b3-4a2e-a556-c025d7df0095
# ╠═2612e3a2-6e6e-494e-b140-720dd2db6ec2
# ╠═f9fe37ef-de81-4d69-9308-cda968851ed2
# ╠═07defcc4-8dca-4caf-a34a-a341818d80ad
# ╠═44cbb2ce-5f43-4bc9-a4c4-c0b9df692cd2
# ╠═6fb4b7e1-a22c-4ff8-bbe9-dbf5de5acd37
# ╠═396a53a3-de0f-4d97-9693-40f3757d66f9
# ╠═6feeaae2-cb01-46ad-ad1d-daaca1caf7ec
# ╟─5ee4f95d-0587-44ab-b543-9b7039d545e6
# ╠═1fd7b586-75ad-42a8-b752-ea82b290cb47
# ╟─77479cd4-513c-4603-9aa0-1acd964c403a
# ╠═7bc2c15c-33f7-43f3-a47f-ca39ffc22071
# ╠═8df71dd5-ba0b-4918-9dc0-791c6c3bad5f
# ╟─c6ce37f4-09f6-43c2-9ee3-286d9e6baf7f
# ╠═a1138e64-5fa5-4a0e-aeef-487ee78a7adc
# ╠═91daae57-94dc-4a8e-b981-75f1406e0768
# ╠═1866c280-89c3-4a71-9dbf-50b583360145
# ╠═0c26f965-0381-4f79-a6ce-0772ae922b3f
# ╠═6e34b91c-c336-4538-a961-60833d37f070
# ╠═57d5decd-8259-4ec2-87ac-44d28625cd7b
# ╠═a0391689-66a2-473f-9704-e12a3d033d13
# ╠═44ab0a25-ab4c-4e90-8619-2a068a285755
# ╟─227a4b71-afbd-4121-930b-696d06ccc9ba
# ╠═253e43df-58fc-4dee-b1c4-35e273499ab7
# ╠═7a94107d-2a45-4458-8d26-5cf836501a1e
# ╠═17b8f17b-0801-45ca-a86b-bba1d78f9ecd
# ╠═51ba2a8e-6f6f-43bd-ac5f-5d238bd41165
# ╠═56be7b5a-3e46-4162-94cb-3f5783efd183
# ╠═04629bcf-5c19-41eb-8903-72947c209cbf
# ╟─f5a3ea2f-4f19-4b0e-af55-74902f2c6485
# ╠═01b5ad85-9f37-4d8b-a29d-e47526f112ec
# ╠═12a30334-899e-4061-b7d2-af8c2346721d
# ╠═1f722acb-f7b9-4d6c-900e-11eae85e0708
# ╠═51dea031-015b-4506-879d-9245c122d608
# ╠═249351c0-6fb6-49ee-ab44-e464f34a1bbe
# ╠═f87e4610-34ac-49f9-9858-0b3ef72eef15
# ╠═8e32cd34-e97e-4cb7-95ba-67f4ed0c9aed
# ╠═1ef7e0c1-56f3-4be7-a138-c088ff973ec6
# ╠═24636618-7310-4557-baac-01d3d5d076dd
# ╠═46f1792a-53a4-4540-b8f3-93747f071ed0
# ╠═1b528c3f-877b-40f5-bc05-7b1c96c8d5c9
# ╠═0215fe76-f9c5-4274-ae43-89960a0caaef
# ╠═7e588ae3-89f3-4b91-8963-f6bf4391a859
# ╠═5e17bf16-2e3d-4372-8283-c43b8ad2550d
# ╠═77918cc8-bfd7-4a84-a1d2-b7d5c721f5ba
# ╠═ff759ae0-35c3-460e-97df-f865e26a0f48
# ╠═cf6a7cbb-1034-4026-a3f3-1e854d2929e2
# ╠═7c4c5136-17d5-4dc5-9e5c-25e2348d2a84
# ╠═494a0dda-0a31-4049-a089-d78437c209cc
# ╠═d2de33be-b037-4255-b9f4-da29abb23754
# ╠═b1b218ff-694b-48bc-b999-806330ad4308
# ╟─9f9b525d-f6f6-4fc0-b5b9-036662fe8ba8
# ╠═9fed63d6-c139-4d28-b00c-37dc1b8dc004
# ╠═20e742e9-a909-4be0-822b-f4ca6015b8aa
# ╠═4537d2a8-dc90-4106-81c5-4a230734b182
# ╠═d7fece88-3327-4435-ab61-b45ff62b3b2e
# ╠═7f9d1703-8bef-4eb4-9864-50724ff168af
# ╠═e904f104-2d01-45f0-a6f1-2040131d8780
# ╠═6ebfff07-c43f-4d4d-8604-9fd4f1de5d25
# ╠═8ad01781-8b5d-4d57-a0b5-7a445fb09b5b
# ╠═edf68b42-4fe9-4e14-b7ed-739e89a1541a
# ╠═7fa19be2-4014-4df0-821e-4c509eca4f28
# ╠═96307998-07a0-45bf-bf10-cd14bfcfe20a
# ╠═35c52910-089c-4507-956a-2b0649507495
# ╠═2ef43371-bfae-4a65-8fa9-d1ab5ade32f1
# ╠═1c242116-e66d-453b-ad62-b6a65cdbe284
# ╠═b3e68e32-c058-467d-b214-aab6a4cd1e19
# ╠═e93f2ee2-5fb3-43cb-ae34-58a0941534f1
# ╠═b7cd6c89-4534-4222-bb27-1ff511692ee1
# ╠═9d74ffbb-4c31-4062-9434-7755f53e4da0
# ╠═975a2008-cf02-4442-9ee9-0b1bbb20889d
# ╠═25ccc096-1d27-43ce-8fe1-5a3a9d7cd54e
# ╠═4eae6b35-e71d-47da-bde9-d67b30ce143a
# ╠═6fe6deb4-ae44-4ca0-9617-95841fdaf791
# ╠═55211531-6eb7-4689-9ecd-74854b8ad884
# ╟─e23e85f9-7667-4f6e-8af6-2516fa292e2b
# ╠═0dd09c1e-67c0-4f23-bd12-41cbef62e4de
# ╠═5355bcc3-2494-473d-9da7-fc38930b8ee7
# ╠═82a491cd-e351-4715-90be-5a2d17e960ce
# ╠═ca228ca7-ff7c-4f73-8fa7-28707c61d8e5
# ╠═ed097e29-d2dc-4fa7-94e5-483b380600cc
# ╠═d32757d2-dc08-488c-9ddd-d3eefefa2db7
# ╠═b4f7cbb0-fb5b-4dc5-b66e-679a8e5b630d
# ╠═422839f0-6da4-46b9-8689-2dd13b03188b
# ╠═194ff30a-31a7-44bc-ac54-722a629457fc
# ╠═dcb868c4-2619-4437-af8d-37842bfb6e63
# ╠═d664ab12-a2c1-4531-a4ff-250ffa3ce9eb
# ╠═08b66f99-f81a-4495-a933-9291e986373a
# ╠═cfccf1a1-22ac-4fb8-8b91-c7af98ad3c4d
# ╠═fbd46bd2-79d7-460e-b0ab-0e34a68a1f0a
# ╠═c2ccf9de-e3cd-4950-9a4d-6f425d261ccb
# ╠═76438e61-d227-41cf-b9ea-e658bc389772
# ╟─b9a4b2b7-be95-4ccb-ad74-9b761abfae8a
# ╠═54d0ee8e-52d6-4b8b-84a9-1ddc66659137
# ╠═750185c9-a317-4978-aa04-486e5bfb7a63
# ╠═02183d68-bf7a-4774-9b48-f50712eeb552
# ╟─9b75409d-55f5-47c3-ab63-8168d31d3d54
# ╠═6408828d-d570-4585-8fa8-24661857fb35
# ╠═b36b594c-23b4-4683-8a12-3fa1b4c8b0d9
# ╠═cdbed68e-0da2-4648-a5d2-61b5b07fb4a2
# ╠═d42795d0-bd69-4c2c-be5b-e27e85199ee3
