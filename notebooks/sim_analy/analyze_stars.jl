### A Pluto.jl notebook ###
# v0.19.46

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

# ╔═╡ 377284f2-dcee-44d3-9a04-728605cea92a
md"""
Given the stellar probabilty file, makes plots based on the 3D properties of the stars in the sample.
"""

# ╔═╡ faeaf38d-8c06-4646-8179-57ffb05f720e
import DensityEstimators

# ╔═╡ f0d2b68a-fae2-4486-a434-a8816e400e84
import TOML

# ╔═╡ b3a16249-b8d9-4a6b-9294-cd654a17dc17
md"""
# Inputs
"""

# ╔═╡ cb6a58a6-9ba9-44b5-95a6-062965c13259
models_dir = "/arc7/home/dboyea/sculptor"

# ╔═╡ 0a73bf88-3f46-4864-97f5-41705ea6913d
model_dir = "/arc7/home/dboyea/sculptor/orbits/orbit1/1e6_V31_r3.2"

# ╔═╡ 29988108-b02c-418c-a720-5766f47c39ff
starsname = "fiducial/stars/exp2d_rs0.10"

# ╔═╡ d7f5a3ed-ae4a-4ea3-b776-00dba6506a88
r_scale = 1

# ╔═╡ f0d74eaa-81e9-4b04-9765-24a0935b1430
starsfile = "/arc7/home/dboyea/sculptor/isolation/1e6/$(starsname)_stars.hdf5"

# ╔═╡ 1b5c00d2-9df6-4a9c-ae32-05abcbf0e41a
paramsfile = "/astro/dboyea/sculptor/isolation/1e6/$starsname.toml"

# ╔═╡ 64350409-6bae-4e1f-be11-b2ec7d48d1f1
fig_dir = joinpath(dirname(model_dir),  "figures"); mkpath(fig_dir)

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

# ╔═╡ f9fe37ef-de81-4d69-9308-cda968851ed2
df_probs = lguys.read_hdf5_table(starsfile)

# ╔═╡ e76a6fa2-c020-48bf-b065-bc7ca51ecd98
p_idx = df_probs.index; probabilities = df_probs.probability

# ╔═╡ 07defcc4-8dca-4caf-a34a-a341818d80ad
length(probabilities)

# ╔═╡ 6fb4b7e1-a22c-4ff8-bbe9-dbf5de5acd37
out =  lguys.Output(model_dir * "/out", weights=probabilities)

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

# ╔═╡ 04e315d5-8684-4206-bc8f-a9483d721b2a
profile = lguys.load_profile(params)

# ╔═╡ 91daae57-94dc-4a8e-b981-75f1406e0768
ρ_s(r) = lguys.calc_ρ(profile, r)

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
	h1 = DensityEstimators.histogram(logr, bins, weights=v_rad .* mass, normalization=:none)
	h2 = DensityEstimators.histogram(logr, bins, weights=mass, normalization=:none)

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

# ╔═╡ 17b8f17b-0801-45ca-a86b-bba1d78f9ecd
lguys.calc_r(snap_f)

# ╔═╡ 51ba2a8e-6f6f-43bd-ac5f-5d238bd41165
lguys.calc_r(snap_f, snap_f.x_cen)

# ╔═╡ 56be7b5a-3e46-4162-94cb-3f5783efd183
snap_f.x_cen

# ╔═╡ 9f9b525d-f6f6-4fc0-b5b9-036662fe8ba8
md"""
# Plots of simulated sample
"""

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

# ╔═╡ e23e85f9-7667-4f6e-8af6-2516fa292e2b
md"""
# Velocity profile & Break Radius Calculation
"""

# ╔═╡ 0dd09c1e-67c0-4f23-bd12-41cbef62e4de
import StatsBase: weights, std

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

# ╔═╡ 422839f0-6da4-46b9-8689-2dd13b03188b
function calc_σvx(snap::lguys.Snapshot; r_max=1)
	rs = lguys.calc_r(snap)
	filt = rs .< r_max
	vs = snap[filt].velocities[1, :] .- snap.v_cen[1]
	masses = snap.weights[filt]
	
	σ = lguys.std(vs, weights(masses))
	return σ * V2KMS
end

# ╔═╡ d664ab12-a2c1-4531-a4ff-250ffa3ce9eb
σv = calc_σv(snap_i, r_max=6) 

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

	h = DensityEstimators.histogram(v_rad[filt] * V2KMS, 20, 
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

	h = DensityEstimators.histogram(v_rad[filt], 15,
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

# ╔═╡ 13a87549-1318-494c-9147-3f71095bf2ef
r_b_kpc = calc_rb(σv, orbit_props["t_last_peri"])

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

# ╔═╡ 54d0ee8e-52d6-4b8b-84a9-1ddc66659137
orbit_props["t_last_peri"]

# ╔═╡ d53669fc-84a1-4138-8445-5f31c3ec44a5
@info lguys.kpc_to_arcmin(r_b_kpc, stars)

# ╔═╡ 9b75409d-55f5-47c3-ab63-8168d31d3d54
md"""
# Evolutionary Properties
"""

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

	idx = 1:30:length(out)
	vs = [calc_σv(out[i]) for i in idx]

	scatter!(out.times[idx] * T2GYR, vs)
	
	fig
end

# ╔═╡ Cell order:
# ╠═377284f2-dcee-44d3-9a04-728605cea92a
# ╠═340ffbbe-17bd-11ef-35c6-63505bb128b7
# ╠═faeaf38d-8c06-4646-8179-57ffb05f720e
# ╠═d401ec4b-048e-4aae-85a8-f7f0d8e44a79
# ╠═f0d2b68a-fae2-4486-a434-a8816e400e84
# ╟─b3a16249-b8d9-4a6b-9294-cd654a17dc17
# ╠═cb6a58a6-9ba9-44b5-95a6-062965c13259
# ╠═0a73bf88-3f46-4864-97f5-41705ea6913d
# ╠═29988108-b02c-418c-a720-5766f47c39ff
# ╠═d7f5a3ed-ae4a-4ea3-b776-00dba6506a88
# ╠═f0d74eaa-81e9-4b04-9765-24a0935b1430
# ╠═1b5c00d2-9df6-4a9c-ae32-05abcbf0e41a
# ╠═64350409-6bae-4e1f-be11-b2ec7d48d1f1
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
# ╠═e76a6fa2-c020-48bf-b065-bc7ca51ecd98
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
# ╠═04e315d5-8684-4206-bc8f-a9483d721b2a
# ╠═91daae57-94dc-4a8e-b981-75f1406e0768
# ╠═1866c280-89c3-4a71-9dbf-50b583360145
# ╠═0c26f965-0381-4f79-a6ce-0772ae922b3f
# ╠═6e34b91c-c336-4538-a961-60833d37f070
# ╠═57d5decd-8259-4ec2-87ac-44d28625cd7b
# ╠═a0391689-66a2-473f-9704-e12a3d033d13
# ╠═44ab0a25-ab4c-4e90-8619-2a068a285755
# ╟─227a4b71-afbd-4121-930b-696d06ccc9ba
# ╠═253e43df-58fc-4dee-b1c4-35e273499ab7
# ╠═17b8f17b-0801-45ca-a86b-bba1d78f9ecd
# ╠═51ba2a8e-6f6f-43bd-ac5f-5d238bd41165
# ╠═56be7b5a-3e46-4162-94cb-3f5783efd183
# ╠═04629bcf-5c19-41eb-8903-72947c209cbf
# ╟─9f9b525d-f6f6-4fc0-b5b9-036662fe8ba8
# ╠═d7fece88-3327-4435-ab61-b45ff62b3b2e
# ╠═e904f104-2d01-45f0-a6f1-2040131d8780
# ╠═6ebfff07-c43f-4d4d-8604-9fd4f1de5d25
# ╟─e23e85f9-7667-4f6e-8af6-2516fa292e2b
# ╠═0dd09c1e-67c0-4f23-bd12-41cbef62e4de
# ╠═d32757d2-dc08-488c-9ddd-d3eefefa2db7
# ╠═b4f7cbb0-fb5b-4dc5-b66e-679a8e5b630d
# ╠═422839f0-6da4-46b9-8689-2dd13b03188b
# ╠═d664ab12-a2c1-4531-a4ff-250ffa3ce9eb
# ╠═08b66f99-f81a-4495-a933-9291e986373a
# ╠═cfccf1a1-22ac-4fb8-8b91-c7af98ad3c4d
# ╠═fbd46bd2-79d7-460e-b0ab-0e34a68a1f0a
# ╠═c2ccf9de-e3cd-4950-9a4d-6f425d261ccb
# ╠═76438e61-d227-41cf-b9ea-e658bc389772
# ╟─b9a4b2b7-be95-4ccb-ad74-9b761abfae8a
# ╠═13a87549-1318-494c-9147-3f71095bf2ef
# ╠═54d0ee8e-52d6-4b8b-84a9-1ddc66659137
# ╠═d53669fc-84a1-4138-8445-5f31c3ec44a5
# ╟─9b75409d-55f5-47c3-ab63-8168d31d3d54
# ╠═cdbed68e-0da2-4648-a5d2-61b5b07fb4a2
# ╠═d42795d0-bd69-4c2c-be5b-e27e85199ee3
