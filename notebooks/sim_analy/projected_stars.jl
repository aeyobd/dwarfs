### A Pluto.jl notebook ###
# v0.20.8

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

# ╔═╡ 7d82571f-4d80-412a-9813-1eee1f916e56
using PyFITS

# ╔═╡ acbb16c7-07ff-4793-b9b0-173678c93e2e
using PlutoUI

# ╔═╡ 33a75908-3d98-4006-a8ef-833d9a161b01
using KernelDensity

# ╔═╡ 377284f2-dcee-44d3-9a04-728605cea92a
md"""
Makes some basic plots of a projected stellar file
"""

# ╔═╡ 145490bd-d87e-4153-865b-342506af8f5a
md"""
# Setup
"""

# ╔═╡ faeaf38d-8c06-4646-8179-57ffb05f720e
import DensityEstimators

# ╔═╡ f0d2b68a-fae2-4486-a434-a8816e400e84
import TOML

# ╔═╡ b80bb3ad-966a-4a9b-9ec5-6e6fab872a66
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

# ╔═╡ 1de340a6-da4e-47a6-b499-ced8c469a70d
@bind inputs confirm(notebook_inputs(;
	galaxyname = TextField(default="ursa_minor"),
	modelname = TextField(60, default="1e6_v37_r5.0/orbit"),
	starsname = TextField(default="exp2d_rs0.13"),
))

# ╔═╡ 0a73bf88-3f46-4864-97f5-41705ea6913d
model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", inputs.galaxyname, inputs.modelname)

# ╔═╡ 29988108-b02c-418c-a720-5766f47c39ff
starsname = inputs.starsname

# ╔═╡ 8b149ffe-cbfc-4026-b08f-130eab7c7093
starsdir = joinpath(dirname(model_dir), "stars", starsname)

# ╔═╡ 64350409-6bae-4e1f-be11-b2ec7d48d1f1
FIGDIR = joinpath(starsdir, "figures")

# ╔═╡ d401ec4b-048e-4aae-85a8-f7f0d8e44a79
begin 
	using LilGuys
	FIGDIR
	FIGSUFFIX = ".projected_stars"
end

# ╔═╡ dadf0dd2-f4cb-4657-837a-f99d5dcdfe93
if isdir(starsdir) && .!isdir(FIGDIR)
	mkdir(FIGDIR)
end

# ╔═╡ 44dec2f8-c149-461d-b586-56b73a97c0a2
obs_today_filename = ENV["DWARFS_ROOT"] * "/observations/$(inputs.galaxyname)/observed_properties.toml"

# ╔═╡ 396cd0a8-1d73-44dd-89db-3243fb9e8ac4
md"""
# File loading
"""

# ╔═╡ 436a5be3-f597-4fc4-80a8-dc5af302ad66
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ 84dc77f7-14b3-4a2e-a556-c025d7df0095
stars = read_fits(joinpath(model_dir, "stars", starsname, "final.fits"))

# ╔═╡ 1e3d61ea-16d0-4bbf-a19e-03e54f0f5a0f
stars_i = read_fits(joinpath(model_dir, "stars", starsname, "initial.fits"))

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
sky_orbit = read_fits(joinpath(model_dir, "skyorbit.fits"))

# ╔═╡ 9b5e75e7-171d-40e9-9148-718bb498c56d
idx_f = orbit_props["idx_f"]

# ╔═╡ 8c8a3bad-3f8c-4bf3-8892-ee6b930b3f97
if "r_ell" ∉ names(stars) 
	ra0, dec0 = lguys.spherical_mean(stars.ra, stars.dec)
	stars[:, :xi], stars[:, :eta] = lguys.to_tangent(stars.ra, stars.dec, ra0, dec0)
	stars[!, :r_ell] = @. sqrt(stars.xi^2 + stars.eta^2)
end

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
		xlabel=L"$\xi$ / degree",
		ylabel=L"$\eta$ / degree",
		xgridvisible=false,
		ygridvisible=false
	)
	
	scatter!(stars.xi, stars.eta, 
		alpha=1, color=:black, markersize=1)
	fig
end

# ╔═╡ d7fece88-3327-4435-ab61-b45ff62b3b2e
function mean_2d(obs_df, values; bins=100, centre=false, limits=nothing)
	if centre
		val_mean = lguys.mean(values, obs_df.weights)
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
		ygridvisible=false,
		xreversed=true
		
	)

	return fig, ax
end

# ╔═╡ 6ebfff07-c43f-4d4d-8604-9fd4f1de5d25
function xi_eta_axis(dx=10, dy=5; kwargs...)
	fig = Figure(;kwargs...)
	

	limits = ((-dx, dx), (-dy, dy))

	ax = Axis(fig[1,1],
		xlabel=L"$\xi$ / degrees",
		ylabel=L"$\eta$ / degrees",
		limits=limits,
		aspect = DataAspect(),
		xgridvisible=false,
		ygridvisible=false,
		xreversed=true,
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
	
		
	h = heatmap!(hi, colorscale=log10, colorrange=(1e-6, maximum(hi.values)))
	errorscatter!([obs_today.ra], [obs_today.dec], color=COLORS[3], size=10)

	# idx = idx_f - 20: idx_f + 20
	# lines!(sky_orbit.ra[idx], sky_orbit.dec[idx])
	
	Colorbar(fig[1, 2], h,
		label="stellar density"
	)
	fig
end

# ╔═╡ 9a5e143b-2295-4db9-a945-642fd6adeef0
let 
	fig, ax = ra_dec_axis(10)

	bins = 100
	limits = ax.limits.val
	x = stars.ra
	y = stars.dec


	hi = Arya.histogram2d(x, y, bins, weights=stars.weights, limits=limits)
	areas = diff(hi.xbins) .* (diff(hi.ybins)')
	hi.values ./= areas
	
		
	h = heatmap!(hi, colorscale=log10, colorrange=(1e-5*maximum(hi.values), maximum(hi.values)), colormap=Reverse(:greys))

	# idx = idx_f - 20: idx_f + 20
	# lines!(sky_orbit.ra[idx], sky_orbit.dec[idx])
	
	Colorbar(fig[1, 2], h,
		label="stellar density"
	)
	fig
end

# ╔═╡ 7b4284a3-65bb-4338-ac9c-5ffd51ca8e80
let
	dr = 3
	fig, ax = xi_eta_axis(dr, dr)

	bins = LinRange(-dr, dr, 100)

	h = Arya.histogram2d(stars.xi, stars.eta, bins, weights=stars.weights)

	p = heatmap!(h.xbins, h.ybins, h.values, 
		colorscale=log10, colorrange=(1e-10*maximum(h.values), maximum(h.values)))

	Colorbar(fig[1, 2], p)

	fig
end

# ╔═╡ edf68b42-4fe9-4e14-b7ed-739e89a1541a
let
	dr = 10
	fig, ax = xi_eta_axis(dr, dr)

	bins = LinRange(-dr, dr, 100)

	h = Arya.histogram2d(stars.xi, stars.eta, bins, weights=stars.weights, normalization=:pdf)

	p = heatmap!(h.xbins, h.ybins, h.values, 
		colorscale=log10, colorrange=(1e-7*maximum(h.values), maximum(h.values)))

	Colorbar(fig[1, 2], p, label="fractional stellar density")

	fig
end

# ╔═╡ b4778d19-cb91-4f0f-97fa-4ef69448f849
save = CairoMakie.save

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

# ╔═╡ 922c8d7d-1685-4c78-80ab-a9d4b6b425ce
sky_orbit.radial_velocity[idx_f]

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

# ╔═╡ 1c7b234b-c440-42c5-a601-b3e515ba31af
let
	dr = 2
	fig, ax = xi_eta_axis(dr, dr)

	bins = LinRange(-dr, dr, 20)

	h = mean_2d(stars,  stars.radial_velocity_gsr, bins=bins, centre=true)

	p = heatmap!(h.xbins, h.ybins, h.values, 
	colormap=:redsblues, # colorrange=(50, 80)
	)

	Colorbar(fig[1, 2], p)

	fig
end

# ╔═╡ bae97026-90ac-401f-80a2-7d53f9b3594b
hist(stars.radial_velocity_gsr, weights=stars.weights, bins=100)

# ╔═╡ 312387f2-84b3-432a-85b1-2fcb8baa073c
stars

# ╔═╡ 3341a6c7-e92f-4e68-94d9-073095217a1c
stars.radial_velocity_gsr[1]

# ╔═╡ 35c52910-089c-4507-956a-2b0649507495
filt_cen = stars.r_ell .< 2

# ╔═╡ 2ef43371-bfae-4a65-8fa9-d1ab5ade32f1


# ╔═╡ 1c242116-e66d-453b-ad62-b6a65cdbe284
sky_orbit[idx_f, :]

# ╔═╡ b3e68e32-c058-467d-b214-aab6a4cd1e19
r_cut = 60 # arcmin

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
	
	errorscatter!([obs_today.pmra], [obs_today.pmdec], xerr=[obs_today_err.pmra], yerror=[obs_today_err.pmra], color=COLORS[3])
	
	scatter!([sky_orbit.pmra[idx_f]], 
		[sky_orbit.pmdec[idx_f]], markersize=10)

	idx = idx_f - 20: idx_f
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
	
	h = DensityEstimators.histogram(x, bins, weights=stars.weights[filt], normalization=:pdf, errors=:weighted)

	
	barplot!(h)
	
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
		ylabel = L"${v}_\textrm{rad}$ / km s$^{-1}$",
		limits=limits,
		title="today",
	)



	filt = stars.xi .^ 2 .+ stars.eta .^ 2 .< r_max ^ 2
	x = stars.distance[filt]
	y = stars.radial_velocity[filt]
	
	h = Arya.hist2d!(ax, x, y, weights=stars.weights[filt], bins=100,
		colorscale=log10, colorrange=(1e-10, nothing)
	)
	
	errorscatter!([obs_today.distance], [obs_today.radial_velocity], xerror=[obs_today_err.distance], yerror=[obs_today_err.radial_velocity], color=COLORS[3])

	scatter!([sky_orbit.distance[idx_f]], 
		[sky_orbit.radial_velocity[idx_f]], markersize=10)
	
	idx = idx_f - 20: idx_f
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
	dr = 40

	
	fig, ax = FigAxis(
		limits = (-1, 2.5, nothing, nothing),
		xlabel = "log r",
		ylabel="radial velocity"
		
	)

	bins = LinRange(sky_orbit.radial_velocity[idx_f] - dr, sky_orbit.radial_velocity[idx_f] + dr, 100)

	hist2d!(log10.(stars.r_ell), stars.radial_velocity, weights=stars.weights,
		colorscale=log10,
		colorrange=(1e-8, 0.3),
		bins=(100, bins)
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

	errorscatter!(lguys.midpoints(r_bins), σs, yerror=err)

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
lguys.arcmin2kpc(240, 86)

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
	
	h = DensityEstimators.histogram(rv, 60, weights=mass, normalization=:pdf, errors=:weighted)

	errorscatter!(lguys.midpoints(h.bins), h.values, yerror=h.err)

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

# ╔═╡ 9eb85c51-cf51-4fb4-b95b-7ba74fad9f8d
 lguys.kpc2arcmin(calc_rb(9.6, orbit_props["t_last_peri"]), obs_today.distance) # kpc

# ╔═╡ 02183d68-bf7a-4774-9b48-f50712eeb552
r_b_arcmin = lguys.kpc2arcmin(r_b_kpc, obs_today.distance)

# ╔═╡ 3b3f3219-fae4-4f3c-be82-27759b7fece7
@info "r break = $r_b_arcmin"

# ╔═╡ a1578e05-dc0e-4ed3-b291-af59254e5eab
md"""
# Quick density profile 
"""

# ╔═╡ 5fcd4b8c-7eb6-4d30-885d-25b857f1263c
nan_filt = isfinite.(stars.r_ell) .& (stars.r_ell .> 0)

# ╔═╡ a8688258-e818-4b7a-b238-5629d99413ed
prof = lguys.StellarDensityProfile(stars.r_ell[nan_filt], weights=stars.weights[nan_filt], normalization=:none,  bins=lguys.Interface.bins_both(stars.r_ell[nan_filt], nothing, bin_width=0.05, num_per_bin=1))

# ╔═╡ 0fb0dfd7-8536-4d8f-b368-d5ca3dc74783
prof_i = lguys.StellarDensityProfile(stars_i.r_ell, weights=stars_i.weights, normalization=:none, bins=lguys.Interface.bins_both(stars_i.r_ell, nothing, bin_width=0.05, num_per_bin=1))

# ╔═╡ 1fde438a-ad46-4b60-bc9f-fddc533d9cdb
prof_expected = lguys.StellarDensityProfile(ENV["DWARFS_ROOT"] * "/observations/$(inputs.galaxyname)/density_profiles/jax_2c_eqw_profile.toml") |> lguys.filter_empty_bins

# ╔═╡ b5fd1bcd-b554-48d3-8472-024cb0bd0792
lguys.GalactocentricFrame().d

# ╔═╡ f433b4f0-3214-469f-909a-3cc66d693263
function get_normalization(prof)
	return -lguys.mean(prof.log_Sigma[1:5])
end

# ╔═╡ 901f5ba3-3dab-41ee-ba4b-50f2bf6ffeff
let 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=L"\log\ R \ /\ \textrm{arcmin}",
		ylabel = L"\log\ \Sigma\ /\ \textrm{fraction\ arcmin^{-2}}",
		limits=((0, 2.5), (-6, nothing))
	)

	errorscatter!(prof_expected.log_R, prof_expected.log_Sigma .+ get_normalization(prof_expected),
		yerror=error_interval.(prof_expected.log_Sigma),
		label="J+24",
		color=:black
	)


	dy = get_normalization(prof)
	
	lines!(prof_i.log_R, prof_i.log_Sigma .+ dy, 
			label="initial")

	lines!(prof.log_R, prof.log_Sigma .+ dy, 
			label="final")

	vlines!(log10(r_b_arcmin), color=:grey, label="break radius", linestyle=:dot)

	axislegend(position=:lb)
	@savefig "density2d_vs_observations" fig
	fig
end

# ╔═╡ 35b870a6-7934-4729-ba92-ab62acc4139d
let 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=L"\log r \ /\ \textrm{arcmin}",
		ylabel = L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
		limits=((nothing, 2.5), (-10, nothing))
	)

	errscatter!(prof_expected.log_r, prof_expected.log_Sigma .+ get_normalization(prof_expected),
		yerr=prof_expected.log_Sigma_err,
		label="J+24",
		color=:black
	)


	dy = get_normalization(prof)
	
	lines!(prof_i.log_r, prof_i.log_Sigma .+ dy, 
			label="initial")

	lines!(prof.log_r, prof.log_Sigma .+ dy, 
			label="final")

	vlines!(log10(r_b_arcmin), color=:grey, label="break radius")

	axislegend(position=:lb)
	@savefig "density2d_vs_observations_nolim" fig
	fig
end

# ╔═╡ c901bf9d-46bf-4cd7-af11-227d477663b4
r_b_arcmin

# ╔═╡ 54912dae-b4d9-49b3-80ef-af6f277cd73b
idx_f

# ╔═╡ b6dff4d6-a89b-4b46-858a-5d490c47eeb7
let 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=L" r \ / \textrm{kpc}",
		ylabel = L"\Gamma",
		limits=((0, 120), (-12, 2))
	)


	

	errscatter!(10 .^ prof.log_r, prof.Gamma, yerr=prof.Gamma_err, 
			label="model")

	errscatter!(10 .^ prof_expected.log_r, prof_expected.Gamma,
		yerr=prof_expected.Gamma_err,
		label="J+24",
		color=:black
	)
	
	vlines!((r_b_arcmin), color=:grey, label="break radius")
	fig
end

# ╔═╡ Cell order:
# ╟─377284f2-dcee-44d3-9a04-728605cea92a
# ╟─1de340a6-da4e-47a6-b499-ced8c469a70d
# ╟─145490bd-d87e-4153-865b-342506af8f5a
# ╠═340ffbbe-17bd-11ef-35c6-63505bb128b7
# ╠═7d82571f-4d80-412a-9813-1eee1f916e56
# ╠═faeaf38d-8c06-4646-8179-57ffb05f720e
# ╠═d401ec4b-048e-4aae-85a8-f7f0d8e44a79
# ╠═f0d2b68a-fae2-4486-a434-a8816e400e84
# ╠═acbb16c7-07ff-4793-b9b0-173678c93e2e
# ╠═33a75908-3d98-4006-a8ef-833d9a161b01
# ╠═b80bb3ad-966a-4a9b-9ec5-6e6fab872a66
# ╠═0a73bf88-3f46-4864-97f5-41705ea6913d
# ╠═29988108-b02c-418c-a720-5766f47c39ff
# ╠═8b149ffe-cbfc-4026-b08f-130eab7c7093
# ╠═64350409-6bae-4e1f-be11-b2ec7d48d1f1
# ╠═dadf0dd2-f4cb-4657-837a-f99d5dcdfe93
# ╠═44dec2f8-c149-461d-b586-56b73a97c0a2
# ╟─396cd0a8-1d73-44dd-89db-3243fb9e8ac4
# ╠═672bdffc-7d61-4727-82cf-c819ebf4aa99
# ╠═9e1a868f-05a0-46d5-a554-7f59a762ec51
# ╠═2d61ba92-d587-42bd-a1e3-1b03e1c6c884
# ╠═436a5be3-f597-4fc4-80a8-dc5af302ad66
# ╠═84dc77f7-14b3-4a2e-a556-c025d7df0095
# ╠═1e3d61ea-16d0-4bbf-a19e-03e54f0f5a0f
# ╠═2612e3a2-6e6e-494e-b140-720dd2db6ec2
# ╠═240077b7-dc20-4cfa-9ada-d3aedcf36a75
# ╠═9b5e75e7-171d-40e9-9148-718bb498c56d
# ╠═8c8a3bad-3f8c-4bf3-8892-ee6b930b3f97
# ╟─638b644f-a1e2-4d7f-a3b7-2d543d556729
# ╠═9fed63d6-c139-4d28-b00c-37dc1b8dc004
# ╠═3655adc4-44a7-4746-9ea9-036f2cac43f3
# ╠═20e742e9-a909-4be0-822b-f4ca6015b8aa
# ╠═d7fece88-3327-4435-ab61-b45ff62b3b2e
# ╠═e904f104-2d01-45f0-a6f1-2040131d8780
# ╠═6ebfff07-c43f-4d4d-8604-9fd4f1de5d25
# ╠═8ad01781-8b5d-4d57-a0b5-7a445fb09b5b
# ╠═9a5e143b-2295-4db9-a945-642fd6adeef0
# ╠═7b4284a3-65bb-4338-ac9c-5ffd51ca8e80
# ╠═edf68b42-4fe9-4e14-b7ed-739e89a1541a
# ╠═b4778d19-cb91-4f0f-97fa-4ef69448f849
# ╠═7fa19be2-4014-4df0-821e-4c509eca4f28
# ╠═922c8d7d-1685-4c78-80ab-a9d4b6b425ce
# ╠═96307998-07a0-45bf-bf10-cd14bfcfe20a
# ╠═1c7b234b-c440-42c5-a601-b3e515ba31af
# ╠═bae97026-90ac-401f-80a2-7d53f9b3594b
# ╠═312387f2-84b3-432a-85b1-2fcb8baa073c
# ╠═3341a6c7-e92f-4e68-94d9-073095217a1c
# ╠═35c52910-089c-4507-956a-2b0649507495
# ╠═2ef43371-bfae-4a65-8fa9-d1ab5ade32f1
# ╠═1c242116-e66d-453b-ad62-b6a65cdbe284
# ╠═b3e68e32-c058-467d-b214-aab6a4cd1e19
# ╠═9d74ffbb-4c31-4062-9434-7755f53e4da0
# ╠═975a2008-cf02-4442-9ee9-0b1bbb20889d
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
# ╠═9eb85c51-cf51-4fb4-b95b-7ba74fad9f8d
# ╠═02183d68-bf7a-4774-9b48-f50712eeb552
# ╠═3b3f3219-fae4-4f3c-be82-27759b7fece7
# ╟─a1578e05-dc0e-4ed3-b291-af59254e5eab
# ╠═5fcd4b8c-7eb6-4d30-885d-25b857f1263c
# ╠═a8688258-e818-4b7a-b238-5629d99413ed
# ╠═0fb0dfd7-8536-4d8f-b368-d5ca3dc74783
# ╠═1fde438a-ad46-4b60-bc9f-fddc533d9cdb
# ╠═b5fd1bcd-b554-48d3-8472-024cb0bd0792
# ╠═f433b4f0-3214-469f-909a-3cc66d693263
# ╠═901f5ba3-3dab-41ee-ba4b-50f2bf6ffeff
# ╠═35b870a6-7934-4729-ba92-ab62acc4139d
# ╠═c901bf9d-46bf-4cd7-af11-227d477663b4
# ╠═54912dae-b4d9-49b3-80ef-af6f277cd73b
# ╠═b6dff4d6-a89b-4b46-858a-5d490c47eeb7
