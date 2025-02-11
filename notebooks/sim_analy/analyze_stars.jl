### A Pluto.jl notebook ###
# v0.20.4

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
import DensityEstimators as DE

# ╔═╡ f0d2b68a-fae2-4486-a434-a8816e400e84
import TOML

# ╔═╡ b3a16249-b8d9-4a6b-9294-cd654a17dc17
md"""
# Inputs
"""

# ╔═╡ 2106bfe1-a53d-4ef8-9111-e191a8056351
starsname = "exp2d_rs0.10"

# ╔═╡ f0d74eaa-81e9-4b04-9765-24a0935b1430
#model_dir = ENV["DWARFS_ROOT"] * "/analysis/sculptor/1e7_V31_r4.2/vasiliev24_L3M11_2x_smallperilmc"
model_dir = ENV["DWARFS_ROOT"] * "/analysis/ursa_minor/1e6_v37_r5.0/orbit_mean.2"

# ╔═╡ 9c42eb0a-029d-46f7-afb0-de03f82c5889
obs_today_filename =  ENV["DWARFS_ROOT"] * "/observations/ursa_minor/observed_properties.toml"

# ╔═╡ 08aa0f76-3d74-45b5-b9e9-6abbf6350910
stars_dir_in = joinpath(model_dir, "../stars/$starsname")

# ╔═╡ 7c64a7c7-232b-47b6-973a-62e1ac21df3a
stars_dir_out = joinpath(model_dir, "stars/$starsname")

# ╔═╡ 64350409-6bae-4e1f-be11-b2ec7d48d1f1
figdir = joinpath(stars_dir_out,  "figures"); mkpath(figdir)

# ╔═╡ 1b5c00d2-9df6-4a9c-ae32-05abcbf0e41a
paramsfile = joinpath(stars_dir_in, "profile.toml")

# ╔═╡ 973955ad-3214-42cf-831f-a782f1d2434a
idx_i = 1 

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
df_probs = lguys.read_hdf5_table(stars_dir_in * "/probabilities_stars.hdf5")

# ╔═╡ e76a6fa2-c020-48bf-b065-bc7ca51ecd98
p_idx = df_probs.index; probabilities = df_probs.probability

# ╔═╡ 07defcc4-8dca-4caf-a34a-a341818d80ad
length(probabilities)

# ╔═╡ 6fb4b7e1-a22c-4ff8-bbe9-dbf5de5acd37
out =  lguys.Output(model_dir, weights=probabilities)

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

# ╔═╡ 77479cd4-513c-4603-9aa0-1acd964c403a
let
	fig = Figure(size=(700, 700))
	r_max = 100

	x = lguys.get_x(snap_f)
	y = lguys.get_y(snap_f)
	z = lguys.get_z(snap_f)
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

# ╔═╡ 7bc2c15c-33f7-43f3-a47f-ca39ffc22071
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc", title="initial")

	bin_range = LinRange(-2, 2, 100)
	bins = (out.x_cen[1, idx_i]  .+ bin_range,  out.x_cen[2, idx_i]  .+ bin_range)

	h = Arya.histogram2d( snap_i.positions[1, :], snap_i.positions[2, :], bins, weights=snap_i.weights)
	
	colorrange =(1e-5, 1) .* maximum(h.values)

		

	Arya.hist2d!(ax, snap_i.positions[1, :], snap_i.positions[2, :], weights=snap_i.weights, bins = bins, colorscale=log10, colorrange=colorrange)

	ax2 = Axis(fig[1,2], aspect=1,
	xlabel = "x / kpc", 
	title="final")


	bins = (out.x_cen[1, idx_f]  .+ bin_range,  out.x_cen[2, idx_f]  .+ bin_range)

	h = Arya.hist2d!(ax2, snap_f.positions[1, :], snap_f.positions[2, :], weights=snap_f.weights, bins = bins, colorscale=log10, colorrange=colorrange)

	Colorbar(fig[1, 3], h, label="stellar density")

    rowsize!(fig.layout, 1, ax.scene.viewport[].widths[1])

	@savefig "density_xy_if"

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

# ╔═╡ 8df71dd5-ba0b-4918-9dc0-791c6c3bad5f
let 
	fig = Figure()
	ax = Axis(fig[1,1], yscale=log10, limits=(nothing, (1e-8, 1)),
		xlabel=L"\epsilon", ylabel="count")

	bins = 200
	stephist!(lguys.calc_ϵ(snap_i), bins=bins, weights=snap_i.weights,
		label = "initial"
	)
	es = lguys.calc_ϵ(snap_f)
	filt = es .> 0
	es = es[filt]
	stephist!(es, bins=bins, weights = snap_f.weights[filt],
		label = "final"
	)

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

# ╔═╡ 44e065d8-e85b-4972-adc7-340a2af22a54
prof_i = lguys.StellarProfile3D(snap_i)

# ╔═╡ 11d8c251-28ae-412c-8a96-646e19adef56
prof_f = lguys.StellarProfile3D(snap_f)

# ╔═╡ 6e34b91c-c336-4538-a961-60833d37f070
function v_rad_hist(snap, bins=100)

	mass = snap.weights

	v_rad = lguys.calc_v_rad(snap)
	logr = log10.(lguys.calc_r(snap))
	x_bins, v_bins, err = lguys.histogram(logr, bins, weights=v_rad .* mass)
	_, counts, _ = lguys.histogram(logr, x_bins, weights=mass)


	return x_bins, v_bins ./ counts, err ./ counts
end

# ╔═╡ dfc045b0-2cf1-4397-b7e5-3b855767b71f
function v_rad_hist_dm(snap, bins=100)


	v_rad = lguys.calc_v_rad(snap)
	logr = log10.(lguys.calc_r(snap))
	x_bins, v_bins, err = lguys.histogram(logr, bins, weights=v_rad)
	_, counts, _ = lguys.histogram(logr, x_bins)


	return x_bins, v_bins ./ counts, err ./ counts
end

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

# ╔═╡ fbd46bd2-79d7-460e-b0ab-0e34a68a1f0a
gaussian(x, μ, σ) = 1/sqrt(2π)* 1/σ * exp(-(x-μ)^2/(2σ^2))

# ╔═╡ c2ccf9de-e3cd-4950-9a4d-6f425d261ccb
function calc_v_tan(snap)
	v = lguys.calc_v(snap)
	v_rad = lguys.calc_v_rad(snap)
	v_tan = @. sqrt(v^2 - v_rad^2)
	return vec(v_tan)
end

# ╔═╡ 76438e61-d227-41cf-b9ea-e658bc389772
let
	snap = snap_f
	logr_max = 1
	vmax = 30

	σv = calc_σv(snap, r_max=10 .^ logr_max)
	
	
	mass = snap.weights
	v_x = (snap.velocities[1, :] .- snap.v_cen[1]) * V2KMS
	logr = log10.(lguys.calc_r(snap))
	filt = logr .< logr_max

	h = DE.histogram(v_x[filt], 25,
		weights=mass[filt], normalization=:pdf, limits=(-vmax, vmax))


	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"v_{x} \quad [\textrm{km\,s^{-1}}]",
		ylabel="stellar density"
	)
	
	scatter!(lguys.midpoints(h.bins), h.values)

	x_model = LinRange(-vmax, vmax, 100)
	y_model = gaussian.(x_model, 0, σv)
	lines!(x_model, y_model, 
		color=COLORS[2], label="gaussian (σ = $(round(σv, digits=2)) km / s)"
	)

	axislegend()

	@savefig "sigma_v_x_final"
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

	ymax = log10(maximum(prof_i.rho))
	
	ax = Axis(fig[1,1], xlabel=L"\log r / \textrm{kpc}", ylabel = L"\log \rho_\star", 
		limits=((-1.9, 1), (ymax - 15, ymax)))

	lines!(prof_i.log_r, log10.(prof_i.rho), label="initial")
	lines!(prof_f.log_r, log10.(prof_f.rho), label="final")
	
	x = LinRange(-2, 1, 1000)
	r = 10 .^ x
	y = log10.(ρ_s.(r))

	lines!(x, y, label="analytic", linestyle=:dot)
	
	vlines!(log10.(r_b_kpc), linestyle=:dash, color=:black, label="break radius")
	
	axislegend(ax, position=:lb)
	
	@savefig "rho_3d_if"

	fig
end

# ╔═╡ 0c26f965-0381-4f79-a6ce-0772ae922b3f
let
	snap = snap_f
	
	mass = snap.weights
	v_rad = lguys.calc_v_rad(snap) 
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

# ╔═╡ eeed6799-b7f7-43d4-82bd-07f178b49f22
r_b_kpc

# ╔═╡ 57d5decd-8259-4ec2-87ac-44d28625cd7b
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"log $r$ / kpc",
		ylabel=L"mean 3D radial velocity in bin (km\,s$^{-1}$)",
		limits=((nothing, 1.5), (-20, 20))
	)
	
	x, y, ye = v_rad_hist(snap_f, 150)
	errscatter!(lguys.midpoints(x), y * V2KMS, yerr=ye*V2KMS)

	vlines!(log10.(r_b_kpc), linestyle=:dash, color=:black)
	text!(log10.(r_b_kpc), -20, text=L"r_b")

	@savefig "break_radius"
	fig
end

# ╔═╡ ab55323b-d402-48e0-903f-dee904488f91
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"log $r$ / kpc",
		ylabel=L"mean 3D radial velocity in bin (km\,s$^{-1}$)",
		limits=((nothing, 1.5), (-20, 20))
	)
	
	x, y, ye = v_rad_hist_dm(snap_f, 200)
	errscatter!(lguys.midpoints(x), y * V2KMS, yerr=ye*V2KMS)

	vlines!(log10.(r_b_kpc), linestyle=:dash, color=:black)
	text!(log10.(r_b_kpc), -20, text=L"r_b")

	fig
end

# ╔═╡ 54d0ee8e-52d6-4b8b-84a9-1ddc66659137
orbit_props["t_last_peri"]

# ╔═╡ d53669fc-84a1-4138-8445-5f31c3ec44a5
r_b_arcmin = lguys.kpc_to_arcmin(r_b_kpc, orbit_props["distance_f"])

# ╔═╡ 27184a9d-07c9-4b4d-ad17-74ed279ed4e3
r_b_arcmin / 60

# ╔═╡ 9b75409d-55f5-47c3-ab63-8168d31d3d54
md"""
# Evolutionary Properties
"""

# ╔═╡ e26b40a5-b69f-4c9c-baf3-0d8aea717a63
function get_M_h(out, rh, skip=1)

	idx = eachindex(out)[1:skip:end]
	M = Vector{Float64}(undef, length(idx))

	for (i, i_out) in enumerate(idx)
		r = calc_r(out[i_out])
		f = r .< rh
		m = sum(out[i_out].weights[f])
		M[i] = m
	end

	return out.times[idx], M
end

# ╔═╡ a985bc1a-6bf7-448a-b0d3-61bbac301819
sum(prof_f.mass_in_shell_err[.!isnan.(prof_f.mass_in_shell_err)])

# ╔═╡ ade76f31-e638-4061-bcbc-509d0512b162
r_h = LilGuys.calc_r_h(profile)

# ╔═╡ 9729d5e4-5e35-4432-82d0-ff53bd491e6e
#t_dm_h, M_dm_h = get_M_h(out, n_rh * r_h)

# ╔═╡ 16c7a61a-6ee4-43cb-8a46-136c8e3956dd
prof_f.bound_mass

# ╔═╡ df229124-325f-4f0c-9cf3-9930002b7767
profiles = lguys.read_structs_from_hdf5(joinpath(stars_dir_out, "stellar_profiles.hdf5"), lguys.StellarProfile)

# ╔═╡ 49052b15-7cfc-488c-974e-7b9c1872d082
profiles_3D = lguys.read_structs_from_hdf5(joinpath(stars_dir_out, "stellar_profiles_3d.hdf5"), lguys.StellarProfile3D)

# ╔═╡ 5b23eaca-c5ab-4907-ab92-7cb6bc0e0381
function get_M_h_fast(profs, rh)

	idx = eachindex(profs)
	M = Vector{Float64}(undef, length(idx))
	times = Vector{Float64}(undef, length(idx))
	M_err = Vector{Float64}(undef, length(idx))

	for i in idx
		prof = profiles_3D[i].second
		rs = prof.log_r
		Ms = cumsum(prof.mass_in_shell)
		M[i] = lguys.lerp(rs, Ms)(log10(rh))
		times[i] = prof.time
	end

	s = sortperm(times[idx])
	return times[idx][s], M[s]
end

# ╔═╡ 0f2af99f-355a-4477-aa5c-265dc509d458
get_M_h_fast(profiles_3D, 1 * r_h)

# ╔═╡ d97d7d5d-6f6d-4f4b-97fe-bb825b71ee80
let
	fig = Figure()

	n_rh = 10
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="stellar mass within $(round(r_b_kpc, digits=2)) kpc (break radius)",
		# yscale=log10,
		#yticks=[1, 0.1],
	)
	t_dm_h2, M_dm_h2 = get_M_h_fast(profiles_3D,  r_b_kpc)

	
	scatterlines!(t_dm_h2 * lguys.T2GYR, M_dm_h2)


	@savefig "M_in_rb"
	
	fig
end

# ╔═╡ cdbed68e-0da2-4648-a5d2-61b5b07fb4a2
let
	fig = Figure()

	n_rh = 10
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="stellar mass within $n_rh rh",
		# yscale=log10,
		#yticks=[1, 0.1],
	)
	t_dm_h2, M_dm_h2 = get_M_h_fast(profiles_3D, n_rh * r_h)

	
	scatterlines!(t_dm_h2 * lguys.T2GYR, M_dm_h2)


	@savefig "M_in_10rh"
	
	fig
end

# ╔═╡ 7431ea86-7009-45e8-a0e5-5a84dd0a57ed
let
	fig = Figure()

	n_rh = 5
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="stellar mass within $n_rh rh",
		# yscale=log10,
		#yticks=[1, 0.1],
	)
	t_dm_h2, M_dm_h2 = get_M_h_fast(profiles_3D, n_rh * r_h)

	
	scatterlines!(t_dm_h2 * lguys.T2GYR, M_dm_h2)


	@savefig "M_in_5rh"
	
	fig
end

# ╔═╡ c9fdee7a-9e7f-4941-8cd2-492cfd427509
let
	fig = Figure()

	n_rh = 3
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="stellar mass within $n_rh rh",
		# yscale=log10,
		#yticks=[1, 0.1],
	)
	t_dm_h2, M_dm_h2 = get_M_h_fast(profiles_3D, n_rh * r_h)

	
	scatterlines!(t_dm_h2 * lguys.T2GYR, M_dm_h2)


	@savefig "M_in_3rh"
	
	fig
end

# ╔═╡ 23ef95b2-c35f-4e66-896c-54a66be59e61
let
	fig = Figure()

	n_rh = 10
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="bound mass",
		# yscale=log10,
		#yticks=[1, 0.1],
	)


	times = [prof.time * T2GYR for (_, prof) in profiles_3D]

	boundmass = [prof.bound_mass  for (_, prof) in profiles_3D]

	scatter!(times, boundmass)

	@savefig "boundmass"
	fig
end

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


	times = [prof.time * T2GYR for (_, prof) in profiles_3D]

	vs = [prof.sigma_vx * V2KMS for (_, prof) in profiles_3D]

	scatter!(times, vs)

	@savefig "vcirc_time"
	fig
end

# ╔═╡ 60ffd004-9770-4b96-8845-0d694cfbced8
profiles_3D[1].second

# ╔═╡ 0112797c-ecb8-4ea7-8a66-365f9fbe952e
let
	fig = Figure()
	
	ax = Axis(fig[1,1],
		limits = (-2.5, 2, -20, 3),
		xlabel = L"log $r$ / kpc",
		ylabel = L"log $\Sigma_\star$",
	)

	colorrange = (prof_i.time, prof_f.time) .* T2GYR
	for (_, prof) in profiles
		lines!(prof.log_r, prof.log_Sigma, color=prof.time* T2GYR, colorrange=colorrange)
	end


	Colorbar(fig[1,2], label="time / Gyr", colorrange=colorrange)

	@savefig "Sigma_all_times"

	fig
end

# ╔═╡ b0b2948a-0299-4b90-83be-4bdcf2a75eaf
let
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel = L"log $r$ / kpc",
		ylabel = L"log $\rho_\star$",
		limits = (-2.5, 2, -15, 5)
	)

	colorrange = (prof_i.time, prof_f.time) .* T2GYR
	for (_, prof) in profiles_3D
		lines!(prof.log_r, log10.(prof.rho), color=prof.time* T2GYR, colorrange=colorrange)
	end


	Colorbar(fig[1,2], label="time / Gyr", colorrange=colorrange)

	@savefig "rho_all_times"

	fig
end

# ╔═╡ a955cee5-609e-46bd-8212-51e18c0d1e86
properties = Dict(
	"break_radius_kpc" => r_b_kpc,
	"break_radius_arcmin" => r_b_arcmin,
)

# ╔═╡ a3311d9e-430b-49af-ac7c-460aa49e19bc
open(joinpath(stars_dir_out, "properties.toml"), "w") do f
	TOML.print(f, properties)
end

# ╔═╡ Cell order:
# ╟─377284f2-dcee-44d3-9a04-728605cea92a
# ╠═340ffbbe-17bd-11ef-35c6-63505bb128b7
# ╠═faeaf38d-8c06-4646-8179-57ffb05f720e
# ╠═d401ec4b-048e-4aae-85a8-f7f0d8e44a79
# ╠═f0d2b68a-fae2-4486-a434-a8816e400e84
# ╟─b3a16249-b8d9-4a6b-9294-cd654a17dc17
# ╠═2106bfe1-a53d-4ef8-9111-e191a8056351
# ╠═f0d74eaa-81e9-4b04-9765-24a0935b1430
# ╠═9c42eb0a-029d-46f7-afb0-de03f82c5889
# ╠═08aa0f76-3d74-45b5-b9e9-6abbf6350910
# ╠═7c64a7c7-232b-47b6-973a-62e1ac21df3a
# ╠═64350409-6bae-4e1f-be11-b2ec7d48d1f1
# ╠═1b5c00d2-9df6-4a9c-ae32-05abcbf0e41a
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
# ╟─5ee4f95d-0587-44ab-b543-9b7039d545e6
# ╠═77479cd4-513c-4603-9aa0-1acd964c403a
# ╠═7bc2c15c-33f7-43f3-a47f-ca39ffc22071
# ╠═865b194d-9efa-41c2-8184-63bcfcd90f6f
# ╠═8df71dd5-ba0b-4918-9dc0-791c6c3bad5f
# ╟─c6ce37f4-09f6-43c2-9ee3-286d9e6baf7f
# ╠═a1138e64-5fa5-4a0e-aeef-487ee78a7adc
# ╠═04e315d5-8684-4206-bc8f-a9483d721b2a
# ╠═91daae57-94dc-4a8e-b981-75f1406e0768
# ╠═44e065d8-e85b-4972-adc7-340a2af22a54
# ╠═11d8c251-28ae-412c-8a96-646e19adef56
# ╟─1866c280-89c3-4a71-9dbf-50b583360145
# ╠═0c26f965-0381-4f79-a6ce-0772ae922b3f
# ╠═6e34b91c-c336-4538-a961-60833d37f070
# ╠═dfc045b0-2cf1-4397-b7e5-3b855767b71f
# ╠═27184a9d-07c9-4b4d-ad17-74ed279ed4e3
# ╠═eeed6799-b7f7-43d4-82bd-07f178b49f22
# ╠═57d5decd-8259-4ec2-87ac-44d28625cd7b
# ╠═ab55323b-d402-48e0-903f-dee904488f91
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
# ╠═fbd46bd2-79d7-460e-b0ab-0e34a68a1f0a
# ╠═c2ccf9de-e3cd-4950-9a4d-6f425d261ccb
# ╠═76438e61-d227-41cf-b9ea-e658bc389772
# ╠═b9a4b2b7-be95-4ccb-ad74-9b761abfae8a
# ╠═13a87549-1318-494c-9147-3f71095bf2ef
# ╠═54d0ee8e-52d6-4b8b-84a9-1ddc66659137
# ╠═d53669fc-84a1-4138-8445-5f31c3ec44a5
# ╟─9b75409d-55f5-47c3-ab63-8168d31d3d54
# ╠═e26b40a5-b69f-4c9c-baf3-0d8aea717a63
# ╠═5b23eaca-c5ab-4907-ab92-7cb6bc0e0381
# ╠═a985bc1a-6bf7-448a-b0d3-61bbac301819
# ╠═ade76f31-e638-4061-bcbc-509d0512b162
# ╠═9729d5e4-5e35-4432-82d0-ff53bd491e6e
# ╠═0f2af99f-355a-4477-aa5c-265dc509d458
# ╠═16c7a61a-6ee4-43cb-8a46-136c8e3956dd
# ╠═d97d7d5d-6f6d-4f4b-97fe-bb825b71ee80
# ╠═cdbed68e-0da2-4648-a5d2-61b5b07fb4a2
# ╠═7431ea86-7009-45e8-a0e5-5a84dd0a57ed
# ╠═c9fdee7a-9e7f-4941-8cd2-492cfd427509
# ╠═23ef95b2-c35f-4e66-896c-54a66be59e61
# ╠═d42795d0-bd69-4c2c-be5b-e27e85199ee3
# ╠═df229124-325f-4f0c-9cf3-9930002b7767
# ╠═49052b15-7cfc-488c-974e-7b9c1872d082
# ╠═60ffd004-9770-4b96-8845-0d694cfbced8
# ╠═0112797c-ecb8-4ea7-8a66-365f9fbe952e
# ╠═b0b2948a-0299-4b90-83be-4bdcf2a75eaf
# ╠═a955cee5-609e-46bd-8212-51e18c0d1e86
# ╠═a3311d9e-430b-49af-ac7c-460aa49e19bc
