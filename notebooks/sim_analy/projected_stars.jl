### A Pluto.jl notebook ###
# v0.19.45

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

# ╔═╡ 33a75908-3d98-4006-a8ef-833d9a161b01
using KernelDensity

# ╔═╡ 377284f2-dcee-44d3-9a04-728605cea92a
md"""
Given a stellar probability file, calculates initial-final density profiles, 
and projects stars onto the sky
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
model_dir = "/arc7/home/dboyea/sculptor/orbits/orbit1/"

# ╔═╡ 29988108-b02c-418c-a720-5766f47c39ff
starsname = "exp2d_rs0.05_stars.fits"

# ╔═╡ 64350409-6bae-4e1f-be11-b2ec7d48d1f1
fig_dir = joinpath(dirname(model_dir),  "figures"); mkpath(fig_dir)

# ╔═╡ 44dec2f8-c149-461d-b586-56b73a97c0a2
obs_today_filename = "../../models/sculptor/mc_orbits/orbit1.toml"

# ╔═╡ 396cd0a8-1d73-44dd-89db-3243fb9e8ac4
md"""
# File loading
"""

# ╔═╡ b9d43c14-8b04-4c38-a35c-4aa29ca59635
md"""
## Misc
"""

# ╔═╡ 436a5be3-f597-4fc4-80a8-dc5af302ad66
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ 84dc77f7-14b3-4a2e-a556-c025d7df0095
stars = lguys.load_fits(joinpath(model_dir, "stars", starsname))

# ╔═╡ 2612e3a2-6e6e-494e-b140-720dd2db6ec2
obs_today_file = TOML.parsefile(obs_today_filename)

# ╔═╡ 672bdffc-7d61-4727-82cf-c819ebf4aa99
obs_today_icrs = lguys.ICRS(;
	ra=obs_today_file["ra"], dec=obs_today_file["dec"],
	distance=obs_today_file["distance"],
	pmra=obs_today_file["pmra"],
	pmdec=obs_today_file["pmdec"],
	radial_velocity=obs_today_file["radial_velocity"],
)

# ╔═╡ 2d61ba92-d587-42bd-a1e3-1b03e1c6c884
obs_today = lguys.transform(lguys.GSR, obs_today_icrs)

# ╔═╡ 9e1a868f-05a0-46d5-a554-7f59a762ec51
obs_today_err = lguys.ICRS(;
	ra=0, dec=0,
	distance=obs_today_file["distance_err"],
	pmra=obs_today_file["pmra_err"],
	pmdec=obs_today_file["pmdec_err"],
	radial_velocity=obs_today_file["radial_velocity_err"],
)

# ╔═╡ 240077b7-dc20-4cfa-9ada-d3aedcf36a75
sky_orbit = lguys.load_fits(joinpath(model_dir, "stars", "sky_orbit.fits"))

# ╔═╡ 9b5e75e7-171d-40e9-9148-718bb498c56d
idx_f = orbit_props["idx_f"]

# ╔═╡ 638b644f-a1e2-4d7f-a3b7-2d543d556729
md"""
# Plots
"""

# ╔═╡ 3655adc4-44a7-4746-9ea9-036f2cac43f3
sort(stars.r_ell / 60)

# ╔═╡ 20e742e9-a909-4be0-822b-f4ca6015b8aa
let
	fig, ax = FigAxis(aspect=1,
		limits=(-2, 2, -2, 2),
		xlabel=L"\xi",
		ylabel=L"\eta",
		xgridvisible=false,
		ygridvisible=false
	)
	
	scatter!(stars.xi, stars.eta, 
		alpha=0.1, color=:black, markersize=3)
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

# ╔═╡ e904f104-2d01-45f0-a6f1-2040131d8780
function ra_dec_axis(ddeg=5; kwargs...)
	fig = Figure(;kwargs...)
	
	dy = ddeg
	obs_c = stars[1, :]
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
	x = stars.ra
	y = stars.dec


	hi = Arya.histogram2d(x, y, bins, weights=stars.weights, limits=limits)
	areas = diff(hi.xbins) .* (diff(hi.ybins)')
	hi.values ./= areas
	
		
	h = heatmap!(hi, colorscale=log10, colorrange=(1e-12, maximum(hi.values)))
	errscatter!([obs_today.ra], [obs_today.dec], color=COLORS[3], size=10)

	# idx = idx_f - 20: idx_f + 20
	# lines!(sky_orbit.ra[idx], sky_orbit.dec[idx])
	
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

	h = Arya.histogram2d(stars.xi, stars.eta, bins, weights=stars.weights)

	p = heatmap!(h.xbins, h.ybins, h.values, 
		colorscale=log10, colorrange=(1e-20, maximum(h.values)))

	Colorbar(fig[1, 2], p)

	fig
end

# ╔═╡ b4778d19-cb91-4f0f-97fa-4ef69448f849
save = CairoMakie.save

# ╔═╡ e3fdb5b0-acf1-4ee1-bd3f-56fbfd60f646
abspath(fig_dir)

# ╔═╡ 7fa19be2-4014-4df0-821e-4c509eca4f28
let
	dr = 10
	fig, ax = xi_eta_axis(dr, dr)

	bins = LinRange(-dr, dr, 100)

	h = mean_2d(stars,  stars.radial_velocity, bins=bins)

	p = heatmap!(h.xbins, h.ybins, h.values, 
	colormap=:redsblues)

	Colorbar(fig[1, 2], p)

	fig
end

# ╔═╡ 96307998-07a0-45bf-bf10-cd14bfcfe20a
let
	dr = 2
	fig, ax = xi_eta_axis(dr, dr)

	bins = LinRange(-dr, dr, 20)

	h = mean_2d(stars,  stars.radial_velocity, bins=bins)

	p = heatmap!(h.xbins, h.ybins, h.values, 
	colormap=:redsblues, # colorrange=(50, 80)
	)

	Colorbar(fig[1, 2], p)

	fig
end

# ╔═╡ 35c52910-089c-4507-956a-2b0649507495
filt_cen = stars.r_ell .< 2

# ╔═╡ 2ef43371-bfae-4a65-8fa9-d1ab5ade32f1


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

	filt = stars.r_ell .< r_cut
	scatter!(stars.ra[filt], stars.dec[filt], 
				alpha=0.1, color=:black, markersize=3)

	fig
end

# ╔═╡ 9d74ffbb-4c31-4062-9434-7755f53e4da0
let
	dr = 0.1
	r_max = r_cut
	
	limits = (stars.pmra[1] .+ (-dr, dr), stars.pmdec[1] .+ (-dr, dr))
	
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel=L"\tilde{\mu}_\delta / \textrm{mas\,yr^{-1}}",
		limits=limits,
		aspect=DataAspect(),
		title="today",
	)

	filt = stars.r_ell .< r_max

	x = stars.pmra[filt]
	y = stars.pmdec[filt]
	
	h = Arya.hist2d!(ax, x, y, weights=stars.weights[filt], bins=100, 
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
	
	limits = (stars.pmra[1] .+ (-dr, dr), nothing)
	
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel=L"\tilde{\mu}_{{\alpha}\!*} / \textrm{mas\,yr^{-1}}",
		ylabel="density",
	limits=limits)

	filt = stars.r_ell .< r_cut

	x = stars.pmra[filt]
	y = stars.pmdec[filt]
	
	h = DensityEstimators.histogram(x, bins, weights=stars.weights[filt], normalization=:pdf)

	
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


	
	x = stars.distance
	y = stars.radial_velocity
	
	h = Arya.hist2d!(ax, x, y, weights=stars.weights, bins=100,
		colorscale=log10, colorrange=(1e-10, nothing)
	)

	
	Colorbar(fig[1, 2], h)
	
	fig
end

# ╔═╡ 4eae6b35-e71d-47da-bde9-d67b30ce143a
hist(stars.radial_velocity, weights=stars.weights, bins=100)

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



	filt = stars.xi .^ 2 .+ stars.eta .^ 2 .< r_max ^ 2
	x = stars.distance[filt]
	y = stars.radial_velocity[filt]
	
	h = Arya.hist2d!(ax, x, y, weights=stars.weights[filt], bins=100,
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

	hist2d!(log10.(stars.r_ell), stars.radial_velocity, weights=stars.weights,
		colorscale=log10,
		colorrange=(1e-8, 0.3),
		bins=100
	)
	
	fig
end

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

	x = log10.(stars.r_ell)
	prob = stars.weights
	v = stars.radial_velocity

	filt = @. !isnan(x)
	filt .&= x .> -0.5
	filt .&= x .< 3

	x = x[filt]
	prob = prob[filt]
	v = v[filt]

	bins = 20
	
	r_bins = DensityEstimators.make_bins(x, (-0.5, 2.3), DensityEstimators.bins_equal_number, n=bins)

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
σv = calc_σv(stars.r_ell, stars.radial_velocity, stars.weights, r_max=r_cut)

# ╔═╡ cfccf1a1-22ac-4fb8-8b91-c7af98ad3c4d
lguys.arcmin_to_kpc(240, 86)

# ╔═╡ fbd46bd2-79d7-460e-b0ab-0e34a68a1f0a
gaussian(x, μ, σ) = 1/sqrt(2π)* 1/σ * exp(-(x-μ)^2/(2σ^2))

# ╔═╡ ed097e29-d2dc-4fa7-94e5-483b380600cc
let
	
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="RV (km/s)",
		ylabel="density",
		title="within $r_cut arcmin"
	)
	
	filt = stars.r_ell .< r_cut
	
	rv = stars.radial_velocity[filt]
	mass = stars.weights[filt]

	μ = DensityEstimators.mean(rv, DensityEstimators.sb.weights(mass))
	
	h = DensityEstimators.histogram(rv, weights=mass, normalization=:pdf)

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

# ╔═╡ 02183d68-bf7a-4774-9b48-f50712eeb552
r_b_arcmin = lguys.kpc_to_arcmin(r_b_kpc, obs_today.distance)

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
# ╠═64350409-6bae-4e1f-be11-b2ec7d48d1f1
# ╠═44dec2f8-c149-461d-b586-56b73a97c0a2
# ╟─396cd0a8-1d73-44dd-89db-3243fb9e8ac4
# ╠═b9d43c14-8b04-4c38-a35c-4aa29ca59635
# ╠═672bdffc-7d61-4727-82cf-c819ebf4aa99
# ╠═9e1a868f-05a0-46d5-a554-7f59a762ec51
# ╠═2d61ba92-d587-42bd-a1e3-1b03e1c6c884
# ╠═436a5be3-f597-4fc4-80a8-dc5af302ad66
# ╠═84dc77f7-14b3-4a2e-a556-c025d7df0095
# ╠═2612e3a2-6e6e-494e-b140-720dd2db6ec2
# ╠═240077b7-dc20-4cfa-9ada-d3aedcf36a75
# ╠═9b5e75e7-171d-40e9-9148-718bb498c56d
# ╟─638b644f-a1e2-4d7f-a3b7-2d543d556729
# ╠═9fed63d6-c139-4d28-b00c-37dc1b8dc004
# ╠═3655adc4-44a7-4746-9ea9-036f2cac43f3
# ╠═20e742e9-a909-4be0-822b-f4ca6015b8aa
# ╠═d7fece88-3327-4435-ab61-b45ff62b3b2e
# ╠═e904f104-2d01-45f0-a6f1-2040131d8780
# ╠═6ebfff07-c43f-4d4d-8604-9fd4f1de5d25
# ╠═8ad01781-8b5d-4d57-a0b5-7a445fb09b5b
# ╠═33a75908-3d98-4006-a8ef-833d9a161b01
# ╠═edf68b42-4fe9-4e14-b7ed-739e89a1541a
# ╠═b4778d19-cb91-4f0f-97fa-4ef69448f849
# ╠═e3fdb5b0-acf1-4ee1-bd3f-56fbfd60f646
# ╠═7fa19be2-4014-4df0-821e-4c509eca4f28
# ╠═96307998-07a0-45bf-bf10-cd14bfcfe20a
# ╠═35c52910-089c-4507-956a-2b0649507495
# ╠═2ef43371-bfae-4a65-8fa9-d1ab5ade32f1
# ╠═1c242116-e66d-453b-ad62-b6a65cdbe284
# ╠═b3e68e32-c058-467d-b214-aab6a4cd1e19
# ╠═9d74ffbb-4c31-4062-9434-7755f53e4da0
# ╠═975a2008-cf02-4442-9ee9-0b1bbb20889d
# ╠═25ccc096-1d27-43ce-8fe1-5a3a9d7cd54e
# ╠═4eae6b35-e71d-47da-bde9-d67b30ce143a
# ╠═6fe6deb4-ae44-4ca0-9617-95841fdaf791
# ╟─e23e85f9-7667-4f6e-8af6-2516fa292e2b
# ╠═0dd09c1e-67c0-4f23-bd12-41cbef62e4de
# ╠═5355bcc3-2494-473d-9da7-fc38930b8ee7
# ╠═ca228ca7-ff7c-4f73-8fa7-28707c61d8e5
# ╠═ed097e29-d2dc-4fa7-94e5-483b380600cc
# ╠═d32757d2-dc08-488c-9ddd-d3eefefa2db7
# ╠═b4f7cbb0-fb5b-4dc5-b66e-679a8e5b630d
# ╠═422839f0-6da4-46b9-8689-2dd13b03188b
# ╠═194ff30a-31a7-44bc-ac54-722a629457fc
# ╠═cfccf1a1-22ac-4fb8-8b91-c7af98ad3c4d
# ╠═fbd46bd2-79d7-460e-b0ab-0e34a68a1f0a
# ╠═c2ccf9de-e3cd-4950-9a4d-6f425d261ccb
# ╟─b9a4b2b7-be95-4ccb-ad74-9b761abfae8a
# ╠═54d0ee8e-52d6-4b8b-84a9-1ddc66659137
# ╠═750185c9-a317-4978-aa04-486e5bfb7a63
# ╠═02183d68-bf7a-4774-9b48-f50712eeb552