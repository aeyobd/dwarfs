### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ f653a17a-8062-4898-9d0c-168a13dfbf6f
using OrderedCollections

# ╔═╡ 9821e893-507b-45c3-b833-6de6c7445ee1
using PyFITS

# ╔═╡ 7d1e9712-1f02-4713-a1f1-a61cebe6d741
using DataFrames, Distributions

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl"); scale_theme_element!(:linewidth, 0.5)

# ╔═╡ ed8867fb-ca90-4577-905d-0ad93ec37097
import TOML

# ╔═╡ 325ff1ea-ddfa-4661-9dd0-22866afca67d
import LinearAlgebra: eigen, diagm

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ be9e982a-904d-497c-949a-1fd265fcb67a
"""
	get_slope(f, x, h)

Calculate the slope of the function f at x with stepsize h
"""
function get_slope(f, x0, h=1e-8x0)
	return (f(x0 + h) - f(x0)) / h
end

# ╔═╡ e5d68f88-0863-4824-89ae-db211b3de68a
module Rapha
	include(joinpath(ENV["DWARFS_ROOT"], "utils/rapha_utils.jl"))
end

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
md"""
# Setup
"""

# ╔═╡ c515f500-27c5-403d-b50f-f0f25c7850cd
σ_fattahi = 0.035

# ╔═╡ d10fb6a6-6bd1-4748-9c7f-b69cccabe176
σ_ludlow = 0.10

# ╔═╡ 22629f89-3d2a-457d-b80a-7df412b9c791
halos_boo3 = OrderedDict(
	"mean" => NFW(v_circ_max = 22 / V2KMS, r_circ_max = 3.9),
	"selected" => NFW(v_circ_max = 35 / V2KMS, r_circ_max = 3.0),
	)

# ╔═╡ f6178d1d-3687-4538-80e7-4f2efe93ab54
obs_props_boo3 = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", "bootes3", "observed_properties.toml"))

# ╔═╡ 567d0e17-9749-405c-9c51-f3d58629d4b6


# ╔═╡ b4bd8f82-ec5d-490b-8be2-08a70f311e24
⊕(x, y) = sqrt(x^2 + y^2)

# ╔═╡ 1e43350b-ebae-42a5-ad9e-778da7087732
dlog_Ms = 2

# ╔═╡ f5794260-7b09-492b-a7e3-0f9559ce0d35
function derived_properties(obs_props; dlog_Ms=0)
	log_Ms = log10(obs_props["M_star"]) + dlog_Ms
	log_Ms_err = obs_props["M_star_err"] / obs_props["M_star"] / log(10)

	v_0 =  LilGuys.find_zero(x -> log10(LilGuys.M_s_from_vel_fattahi(x)) + 10 - log_Ms, [0.01, 0.3])
	log_v_0  = log10(v_0)
	m_fattahi = get_slope(x -> log10(LilGuys.M_s_from_vel_fattahi(10^x)), log_v_0)
	log_r_0 = log10(LilGuys.Ludlow.solve_rmax(v_0))

	
	σy_fattahi = (log_Ms_err / m_fattahi) ⊕ σ_fattahi

	m_ludlow = 1/  get_slope(x -> log10(LilGuys.Ludlow.solve_rmax(10 ^ x)),log_v_0)


	σx_ludlow = 1/2 * (log10(LilGuys.Ludlow.solve_rmax(v_0, -σ_ludlow)) - log10(LilGuys.Ludlow.solve_rmax(v_0, σ_ludlow)))
	σy_ludlow = σx_ludlow * m_ludlow
	
	Σxx = (σy_fattahi^2 + σy_ludlow^2) / m_ludlow^2
	Σyy = (σy_fattahi^2) 
	Σxy = (σy_fattahi^2) / m_ludlow
	Σ = [
		Σxx Σxy
		Σxy Σyy]
	μ = [log_r_0, log_v_0 + log10(V2KMS)]

	R_h = LilGuys.arcmin2kpc(obs_props["R_h"], obs_props["distance"]) 


	return (;
		log_Ms = log_Ms,
		v_0 = v_0,
		log_v_0 = log_v_0,
		log_v_0_err = σy_fattahi,
		R_h = R_h,
	)
end

# ╔═╡ 4f4dd26d-3b93-49fe-9d9a-472f504e365d
der_props_boo3 = derived_properties(obs_props_boo3)

# ╔═╡ cc422906-8ec9-43ec-9c34-3ccad2509e74
der_props_boo3_100Ms = derived_properties(obs_props_boo3, dlog_Ms=2)

# ╔═╡ e2a2e7b0-fe92-4abc-9b52-842e07138bcc
10 ^(Measurement(der_props_boo3.log_v_0, der_props_boo3.log_v_0_err)) * V2KMS

# ╔═╡ 27425538-2104-46fe-8fc9-e6c28c0c6041
α_exp = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 0ec48582-0a73-4b57-87bc-24c19f79209a
R2r = LilGuys.r_h(LilGuys.Exp2D()) / LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 1f74ca1b-ac53-425c-8f89-1afe361525cd
function solve_v_from_dispersion(r0, σv, R_h)
	halo(v) = NFW(r_circ_max=r0, v_circ_max=v)
	stars = LilGuys.Exp2D(R_s=R_h/α_exp)

	f(v) = LilGuys.σv_star_mean(halo(v), stars) - σv/V2KMS
	v0 = LilGuys.find_zero(f, [0.01, 0.5])
	return v0
end

# ╔═╡ 38712f63-8f52-4b1c-8da1-788101798f55
function solve_v_from_dispersion_simple(r0, σv, R_h)
	f(v) = LilGuys.v_circ(NFW(r_circ_max=r0, v_circ_max=v), R_h*R2r)/sqrt(3) - σv/V2KMS
	
	return LilGuys.find_zero(f, [0.01, 0.5])
end

# ╔═╡ dc0add50-5003-4e6e-a814-32ca84aaf187
function σv_wolf(halo, R_h)
	return LilGuys.v_circ(halo, R_h*R2r)/sqrt(3)
end

# ╔═╡ a0e7532b-2aaf-4fb3-af00-1d0912787c35
function plot_error_ellipses!(μ, Σ; kwargs...)
	λ, V = eigen(Σ)
	t = LinRange(0, 2π, 10000)

	Λ = diagm(λ)

	for sigma in [2.3, 6.17, 11.8]
		xy = [μ .+ V * sqrt.(sigma * Λ) * [cos(tt), sin(tt)] for tt in t]
	
		lines!(10 .^ first.(xy), 10 .^ last.(xy); kwargs...)
	end
end

# ╔═╡ 21c57d80-3dda-412f-a912-0248a2f6d600
smallfontsize = @lift 0.8 * $(theme(:fontsize))

# ╔═╡ e9ccfc08-12c7-43af-a372-7f8c04017462
function plot_fattahi!(log_v_0, log_v_0_err; label="Fattahi+18", shade=false)
	
	v_0 = V2KMS * 10^log_v_0
	v_l, v_h = V2KMS .* 10 .^ ((log_v_0 - log_v_0_err), (log_v_0 + log_v_0_err))
	color = COLORS[3]

	if shade
		hspan!(v_l, v_h, xmax=0.2, color=(color, 0.1))
	end
	hlines!(v_0, color=(color, 0.5), xmax=0.2)
	text!(0.1, v_0, text=label, color=color, fontsize=smallfontsize)
end

# ╔═╡ 4ff87ef7-2132-4be9-9b71-15510fe3345a
function log_derivative(f, x0; h=0.001)
	y0 = f(x0)

	return (log10(f(x0 * 10^h)) - log10(y0)) / h
end

# ╔═╡ 5106a242-4fe5-47d8-94f7-45c12053c882
function plot_sigma_v!(σv; R_h, x0=15, rf=1, vmax = solve_v_from_dispersion_simple, kwargs...)
	color = COLORS[2]
	
	x = LinRange(0.1, 30, 300)
	y = vmax.(x, σv, R_h) * V2KMS
	lines!(x, y; color=color, kwargs...)

	y0 = vmax(x0, σv, R_h) * V2KMS
	
	slope = log_derivative(x->vmax(x, σv, R_h) * V2KMS, x0)
	θ = @lift atan(slope * $rf)

	if σv == 8
		label = L"$\sigma_\textrm{v} = %$σv$\,km\,s$^{-1}$"
	else
		label = string(σv)
	end
	
	text!(x0, y0, text=label, rotation=θ, color=color, fontsize=smallfontsize, align=(:center, :bottom))
end

# ╔═╡ 12177747-6d65-42d4-be9f-c4a86d37e65e
function plot_mass_concentration!(; rf=1)
	color = COLORS[8]

	vs = LinRange(10, 90, 1000) / V2KMS
	rs = LilGuys.Ludlow.solve_rmax.(vs)
	xl = (LilGuys.Ludlow.solve_rmax.(vs, -σ_ludlow))
	xh =  (LilGuys.Ludlow.solve_rmax.(vs, σ_ludlow))
	
	lines!(rs, vs*V2KMS, color=color)
	band!(vs*V2KMS, xl, xh, color=(color, 0.1), direction=:y)

	y0 = 12
	x0 = LilGuys.Ludlow.solve_rmax(y0 / V2KMS)
	slope = 1/log_derivative(x->LilGuys.Ludlow.solve_rmax(x / V2KMS), y0)
	θ = @lift atan(slope * $rf)
	
	text!(x0, y0, text="Ludlow+16", color=color, rotation =θ)
end

# ╔═╡ 2f97e7db-b2d5-44af-9c03-72e78f3edcde
function rotation_factor(ax, log=false)
	if log
		limits = ax.finallimits
		x1 = @lift $limits.origin[1]
		x2 = @lift $x1 + $limits.widths[1]
		y1 = @lift $limits.origin[2]
		y2 = @lift $y1 + $limits.widths[2]

		Δx_data = @lift log10($x2) - log10($x1)
		Δy_data = @lift log10($y2) - log10($y1)
	else
		Δx_data = @lift $(ax.finallimits).widths[2]
		Δy_data = @lift $(ax.finallimits).widths[2]
	end
	Δx_screen = @lift $(ax.scene.viewport).widths[1]
	Δy_screen = @lift $(ax.scene.viewport).widths[2]
	correction_factor = @lift $Δx_data / $Δx_screen * $Δy_screen / $Δy_data
	return correction_factor
end

# ╔═╡ 4bf130aa-7cc4-470b-bbe6-bdc86f831683
LilGuys.mean_density(NFW(r_circ_max=1, v_circ_max=0.001), der_props_boo3.R_h*R2r)

# ╔═╡ c66ed3b6-4d1c-43ad-ac22-eced49adc811
function tidal_track_axis(gs)
	ax = Axis(gs,
		ylabel=L"$\,\textrm{v}_\textrm{max}$ / km\,s$^{-1}$",
		xlabel=L"$\,r_\textrm{max}$ / kpc",
		xscale=log10,
		yscale=log10,
		limits=(0.1, 30, 10, 80),
		xminorticks = [0.1:0.1:1; 1:10; 10:10:30],
			  xticks = [0.1, 1, 10],
	)
end

# ╔═╡ 7d36f094-70bc-4114-b920-3655f73dc75f
function plot_halo_constraints(gs, der_props, halos)
	ax = tidal_track_axis(gs)

	rf = rotation_factor(ax, true)


	plot_mass_concentration!(rf=rf)
	R_h = der_props.R_h
	plot_sigma_v!(7.7; x0 = 1, rf=rf, R_h=R_h)
	plot_sigma_v!(7.7-1.5; x0 = 1, rf=rf, R_h=R_h)
	plot_sigma_v!(7.7+2; x0 = 1, rf=rf, R_h=R_h)


	plot_fattahi!(der_props.log_v_0, der_props.log_v_0_err)

	
	# for (label, halo) in halos
	# 	y = LilGuys.v_circ_max(halo) * V2KMS
	# 	x = LilGuys.r_circ_max(halo) 
	# 	scatter!((x), (y), label=string(label))
	# end

	

	ax
end

# ╔═╡ 391d2ecc-424e-4e81-adc8-0d4b03ed9a40
function plot_tidal_track!(v0, r0)
	x, y = LilGuys.EN21_tidal_track(r0, v0/V2KMS, x_min=0.01)

	lines!(x, y*V2KMS, color=COLORS[5], alpha=0.5)
end

# ╔═╡ 8f8eb577-5af9-4d08-ac50-77130a8fa880
import Agama

# ╔═╡ 5e1388e9-1674-4fbc-865c-d1eb9ee74445
pot = Agama.Potential(file = joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))

# ╔═╡ 14ceb32c-19bb-44ec-b795-c807fad60525
function plot_jacobi!(pericentre)
	ρ_mean = Agama.enclosed_mass(pot, pericentre) / (4π/3 * pericentre^3)


	r = logrange(0.1, 30, 100)

	v = [LilGuys.find_zero(v -> LilGuys.mean_density(NFW(r_circ_max=rr, v_circ_max=v, c=100), der_props_boo3.R_h*R2r) - 3ρ_mean, [1e-2, 1]) for rr in r]

	lines!(r, v*V2KMS, color=COLORS[1])
	idx = argmin(abs.(r .- 3))
	@info r[idx], v[idx] * V2KMS

	text!(r[idx], v[idx]*V2KMS, text="peri = $pericentre", color=COLORS[1])
	
end

# ╔═╡ 7ecee6bf-69b3-4570-aa70-8491b96a46b3
Agama.enclosed_mass(pot, 14) / (4π/3 * 14^3)

# ╔═╡ ec8c28f1-b6d6-4616-a6b3-302d8dae1d6d
(9.7/V2KMS)^2/ (4π/3 * der_props_boo3.R_h^2)

# ╔═╡ a12a2cfb-d94c-44ec-9a27-782e0120ac97
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 # yscale=log10,
			 # yticks=[1, 10, 100],
			 xlabel = "pericentre / kpc",
			 ylabel = L"Jacob $\sigma_\text{v,\ \star}$ (km/s)")

	x = LinRange(0.5, 30, 1000)

	ρ_mean = @. Agama.enclosed_mass(pot, x) / (4π/3 * x^3)
	r_h = der_props_boo3.R_h * R2r

	v = @. sqrt(3ρ_mean * 4π/3*r_h^2)
	lines!(x, v*V2KMS / sqrt(3))
	ylims!(6, 30)
	fig
end


# ╔═╡ 8b287edc-bcbc-4ce2-b355-1d174bc9dc7a
let
	x = 20
	ρ_mean = @. Agama.enclosed_mass(pot, x) / (4π/3 * x^3)
	r_h = der_props_boo3.R_h * R2r

	v = @. sqrt(3ρ_mean * 4π/3*r_h^2)
	v * V2KMS / sqrt(3)
end

# ╔═╡ d800f480-3847-43fd-a458-aaaece1d9dcd
LilGuys.Ludlow.solve_rmax(30/V2KMS, 0.0)

# ╔═╡ a15dbc88-8dbf-48d4-8a72-1881151c8e26
let
	fig = Figure(size=(4, 3) .* 72)

	plot_halo_constraints(fig[1,1], der_props_boo3, halos_boo3)


	for (v, δc) in [(30, 0), (30, 3*σ_ludlow), (25, 0.3), (30, 0.15)]
		r = LilGuys.Ludlow.solve_rmax(v/V2KMS, δc)
		plot_tidal_track!(v, r)
	end

	annotation!(0, 30, der_props_boo3.R_h, 10, text=L"R_h")
	plot_fattahi!(der_props_boo3_100Ms.log_v_0, der_props_boo3_100Ms.log_v_0_err, label="100x M⋆", shade=true)
	# plot_jacobi!(10)

	# plot_jacobi!(18)
	ylims!(10, 40)

	plot_sigma_v!(15; x0 = 1, R_h=der_props_boo3.R_h)

	@savefig "initial_halo_selection"
	fig
end

# ╔═╡ 273d1aae-1b1e-42b8-85e4-efdcbcd04d06
orbit_props = read_fits(joinpath(ENV["DWARFS_ROOT"], "orbits/bootes3/EP2020/orbital_properties.fits"))

# ╔═╡ 079fae18-072c-4a1c-b18b-016421d3efbe
md"""
# MCMC samples
"""

# ╔═╡ 56095498-3113-46b5-afb9-2e5e2d501ff1
α_exp_cusp = LilGuys.r_circ_max(LilGuys.ExpCusp())

# ╔═╡ 1aa9e317-2f53-43bc-af7d-9294da42b607
samples = let
	df = DataFrame()


	N = 100_000
	vmax = rand(Uniform(20, 35), N) / V2KMS

	lc = abs.(randn(N)) * 0.1
	rmax = [LilGuys.Ludlow.solve_rmax(v, c) for (v, c) in zip(vmax, lc)]

	df[!, :v_circ_max] = vmax
	df[!, :r_circ_max] = rmax
	df[!, :delta_log_c] = lc

	peri = rand(LogNormal(7, 0.3), N)

	df[!, :peri] = orbit_props.pericentre[1:N]
	df[!, :t_end] = 10 * rand(N) ./ T2GYR
	df[!, :period] = -orbit_props.period[1:N]
	df[!, :n_peri] = ceil.(df.t_end ./ df.period)

	vr_end = [Rapha.rapha_final_halo(rmax[i], vmax[i], orbit_props.pericentre[i], orbit_props.apocentre[i], df.n_peri[i])
					for i in 1:N]

	df[!, :v_circ_end] = last.(vr_end)
	df[!, :r_circ_end] = first.(vr_end)



	r_s = df.r_circ_end ./ α_exp_cusp
	v0 = LilGuys.v_circ_max.(LilGuys.ExpCusp.(1, r_s))
	M = @. (df.v_circ_end/v0)^2
	halos_final_samples = LilGuys.ExpCusp.(M, r_s)

	
	df[!, :M] = [LilGuys.mass(LilGuys.TruncNFW(r_circ_max=rmax[i], v_circ_max=vmax[i], trunc=20, xi=3)) for i in eachindex(rmax)]
	df[!, :M_end] = M
	df[!, :r_s_end] = r_s
	df[!, :sigma_v] = σv_wolf.(halos_final_samples, der_props_boo3.R_h)

	
	df
end

# ╔═╡ fa0b0f0b-4efe-43f4-b879-572bd3197dbb
let
	fig = Figure(size=(4, 3) .* 72)

	plot_halo_constraints(fig[1,1], der_props_boo3, halos_boo3)


	annotation!(0, 30, der_props_boo3.R_h, 10, text=L"R_h")
	plot_fattahi!(der_props_boo3_100Ms.log_v_0, der_props_boo3_100Ms.log_v_0_err, label="100x M⋆", shade=true)

	ylims!(10, 40)

	plot_sigma_v!(15; x0 = 1, R_h=der_props_boo3.R_h)

	scatter!(samples.r_circ_max, samples.v_circ_max*V2KMS, markersize=1)


	fig
end

# ╔═╡ a09669c6-72fd-418d-957c-ffa2f6434c07
let
	fig = Figure(size=(4, 3) .* 72)

	plot_halo_constraints(fig[1,1], der_props_boo3, halos_boo3)


	annotation!(0, 30, der_props_boo3.R_h, 10, text=L"R_h")
	plot_fattahi!(der_props_boo3_100Ms.log_v_0, der_props_boo3_100Ms.log_v_0_err, label="100x M⋆", shade=true)

	ylims!(10, 40)

	plot_sigma_v!(15; x0 = 1, R_h=der_props_boo3.R_h)

	scatter!(samples.r_circ_max, samples.v_circ_max*V2KMS, markersize=1)

	p = scatter!(samples.r_circ_end, samples.v_circ_end*V2KMS, markersize=1, color=log10.(samples.peri), colorrange=(0, 1.3))

	Colorbar(fig[1,2], p, label="log peri")

	@savefig "initial_halo_selection"
	fig
end

# ╔═╡ d8e7b34f-4938-4ea1-a4c2-cecad0f3e5e5
filt_samples = obs_props_boo3["sigma_v"] - obs_props_boo3["sigma_v_em"] .< samples.sigma_v * V2KMS .< obs_props_boo3["sigma_v"] + obs_props_boo3["sigma_v_ep"]

# ╔═╡ 8faefce1-f32a-4119-8e1c-9dcc13293da6
let
	fig = Figure(size=(4, 3) .* 72)

	plot_halo_constraints(fig[1,1], der_props_boo3, halos_boo3)


	annotation!(0, 30, der_props_boo3.R_h, 10, text=L"R_h")
	plot_fattahi!(der_props_boo3_100Ms.log_v_0, der_props_boo3_100Ms.log_v_0_err, label="100x M⋆", shade=true)

	ylims!(10, 40)

	plot_sigma_v!(15; x0 = 1, R_h=der_props_boo3.R_h)

	scatter!(samples.r_circ_max[filt_samples], samples.v_circ_max[filt_samples]*V2KMS, markersize=1)

	p = scatter!(samples.r_circ_end[filt_samples], samples.v_circ_end[filt_samples]*V2KMS, markersize=1, color=log10.(samples.peri[filt_samples]), colorrange=(0, 1.3))

	Colorbar(fig[1,2], p, label="log peri")

	@savefig "initial_halo_selection"
	fig
end

# ╔═╡ 990172f0-6303-4cdb-8637-845dd60ce2e0
LilGuys.mean(filt_samples)

# ╔═╡ 6158b971-e5fb-4ef7-84e6-e0cae6072686
hist(samples.n_peri)

# ╔═╡ c263e88e-58db-4273-943e-593e43070a02
hist(samples.n_peri[samples.v_circ_end*V2KMS .> 10])

# ╔═╡ 5433f056-2a84-4439-85f1-1f4786162f5d
hist(samples.t_end[samples.v_circ_end*V2KMS .> 10])

# ╔═╡ cca4b599-298e-4033-b9aa-36f9c5ba85e2
hist(samples.v_circ_max[filt_samples] * V2KMS)

# ╔═╡ ffcd3315-5e32-47cb-b831-7a19cb1c2162
hist(samples.v_circ_max*V2KMS)

# ╔═╡ bcf8c8ae-cf45-4269-b6cf-c197147ff6ba
hist(samples.r_circ_max[filt_samples])

# ╔═╡ bd0d8191-841c-4dcb-abf0-6f9f1728a761
hist(samples.r_circ_max)

# ╔═╡ a136ee90-aad1-4bca-a40e-554ba921e90e

let
	fig = Figure()

	ax = Axis(fig[1,1], xlabel="r max", ylabel = "sigma v stars end / km/s")
		
	scatter!(samples.r_circ_max, samples.sigma_v * V2KMS, markersize=1, alpha=0.1, color=samples.n_peri,)


	ax = Axis(fig[1,2], xlabel = "v max")

	scatter!(samples.v_circ_max * V2KMS, samples.sigma_v * V2KMS, markersize=1, alpha=0.1, color=samples.peri, colorrange=(0, 20))

	fig
end

# ╔═╡ f79760a9-8fb7-4249-994f-29f54f5c4c3c
let
	filt = filt_samples	

	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "pericentre",
			 ylabel= "# pericentres",
			 yticks=0:7, 
			 yminorticksvisible=false)

	
	p = scatter!(samples.peri[filt] ./ samples.n_peri[filt] .^ 0.0, samples.n_peri[filt] .+ 0.1*randn(sum(filt)), markersize=1, alpha=1, color=samples.v_circ_max[filt] * V2KMS)

	for i in 1:maximum(samples.n_peri[filt])
		filt_n = (samples.n_peri .== i ) .& filt
		peri_med = median(samples.peri[filt_n])
		peri_low, peri_high = quantile(samples.peri[filt_n], [0.16, 0.84])
		errorscatter!([peri_med], [i], xerror =[(peri_med-peri_low, peri_high-peri_med)], color=:black, alpha=0.5)
		@info peri_med, peri_low, peri_high
	end
	
	

	Colorbar(fig[1,2], p, label="vmax")
	fig
end

# ╔═╡ 0555a067-17ac-43a4-82b8-440f0de638a6
let
	filt = filt_samples	

	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "pericentre",
			 ylabel= "# pericentres",
			 yticks=0:7, 
			 yminorticksvisible=false)

	
	p = scatter!(samples.peri[filt] ./ samples.n_peri[filt] .^ 0.0, samples.n_peri[filt] .+ 0.1*randn(sum(filt)), markersize=1, alpha=1, color=samples.delta_log_c[filt] * 10)

	for i in 1:maximum(samples.n_peri[filt])
		filt_n = (samples.n_peri .== i ) .& filt
		peri_med = median(samples.peri[filt_n])
		peri_low, peri_high = quantile(samples.peri[filt_n], [0.16, 0.84])
		errorscatter!([peri_med], [i], xerror =[(peri_med-peri_low, peri_high-peri_med)], color=:black, alpha=0.5)
		@info peri_med, peri_low, peri_high
	end
	
	

	Colorbar(fig[1,2], p, label="concentration (sigma)")
	fig
end

# ╔═╡ 933e9188-ec1f-4ae8-81ef-7e51adff262a
scatter(samples.r_circ_max[filt_samples], samples.v_circ_max[filt_samples]*V2KMS, alpha=1, markersize=1, color=samples.n_peri[filt_samples])

# ╔═╡ 3426557d-074d-471d-b8b8-f18be83a5d78
scatter(samples.v_circ_max[filt_samples], samples.peri[filt_samples], alpha=0.5, markersize=1)

# ╔═╡ 70a701b3-68e5-45d8-a7ba-c8bff151e8a6
scatter(samples.delta_log_c[filt_samples], samples.peri[filt_samples], alpha=0.5, markersize=1)

# ╔═╡ aa5f86a8-b35d-4849-af2c-b9d68fb53d77
let
	filt = filt_samples	

	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "pericentre",
			 ylabel= "# pericentres",
			 yticks=0:7, 
			 yminorticksvisible=false,
			 limits=(0, 30, nothing, nothing)
			  
)

	
	p = scatter!(samples.peri[filt] ./ samples.n_peri[filt] .^ 0.0, samples.n_peri[filt] .+ 0.1*randn(sum(filt)), markersize=1, alpha=1, color=samples.sigma_v[filt]*V2KMS)

	for i in 1:maximum(samples.n_peri[filt])
		filt_n = (samples.n_peri .== i ) .& filt
		peri_med = median(samples.peri[filt_n])
		peri_low, peri_high = quantile(samples.peri[filt_n], [0.16, 0.84])
		errorscatter!([peri_med], [i], xerror =[(peri_med-peri_low, peri_high-peri_med)], color=:black, alpha=0.5)
		@info peri_med, peri_low, peri_high
	end
	
	

	Colorbar(fig[1,2], p, label="σv / (km/s)")
	fig
end

# ╔═╡ 42f779b7-a1d1-4ed9-bd41-42a59143bc01
let
	f = scatter(samples.peri, samples.n_peri .+ 0.1*randn(sum(samples.v_circ_end*V2KMS .> 0)), markersize=1, alpha=0.05)

		lines!([0, 15] .+ 2.5, [0, 15] .* 0.5)

	f
end

# ╔═╡ 0bda0872-7566-4199-a85e-673e241783dc
let
	
	f = scatter(samples.peri[filt_samples], samples.t_end[filt_samples]*T2GYR, markersize=1)

	lines!([0, 10], [0, 10])
	f
end

# ╔═╡ 697b83ce-1787-484a-aa8a-12b91d112d88
let
	fig = Figure(size=(4, 3) .* 72)

	plot_halo_constraints(fig[1,1], der_props_boo3, halos_boo3)


	annotation!(0, 30, der_props_boo3.R_h, 10, text=L"R_h")
	plot_fattahi!(der_props_boo3_100Ms.log_v_0, der_props_boo3_100Ms.log_v_0_err, label="100x M⋆", shade=true)

	ylims!(10, 40)

	plot_sigma_v!(15; x0 = 1, R_h=der_props_boo3.R_h)

	scatter!(samples.r_circ_max, samples.v_circ_max*V2KMS, markersize=3)

	p = scatter!(samples.r_circ_end, samples.v_circ_end*V2KMS, markersize=3, color=(samples.n_peri))

	Colorbar(fig[1,2], p, label="# peri")

	@savefig "initial_halo_selection"
	fig
end

# ╔═╡ 03a263f1-b2f2-43a9-96c8-7895eac5632b
	scatter(samples.r_circ_end, samples.v_circ_end*V2KMS)


# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═ed8867fb-ca90-4577-905d-0ad93ec37097
# ╠═f653a17a-8062-4898-9d0c-168a13dfbf6f
# ╠═325ff1ea-ddfa-4661-9dd0-22866afca67d
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═be9e982a-904d-497c-949a-1fd265fcb67a
# ╠═e5d68f88-0863-4824-89ae-db211b3de68a
# ╟─3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═c515f500-27c5-403d-b50f-f0f25c7850cd
# ╠═d10fb6a6-6bd1-4748-9c7f-b69cccabe176
# ╠═22629f89-3d2a-457d-b80a-7df412b9c791
# ╠═f6178d1d-3687-4538-80e7-4f2efe93ab54
# ╠═567d0e17-9749-405c-9c51-f3d58629d4b6
# ╠═b4bd8f82-ec5d-490b-8be2-08a70f311e24
# ╠═1e43350b-ebae-42a5-ad9e-778da7087732
# ╠═f5794260-7b09-492b-a7e3-0f9559ce0d35
# ╠═4f4dd26d-3b93-49fe-9d9a-472f504e365d
# ╠═cc422906-8ec9-43ec-9c34-3ccad2509e74
# ╠═e2a2e7b0-fe92-4abc-9b52-842e07138bcc
# ╠═27425538-2104-46fe-8fc9-e6c28c0c6041
# ╠═0ec48582-0a73-4b57-87bc-24c19f79209a
# ╠═1f74ca1b-ac53-425c-8f89-1afe361525cd
# ╠═38712f63-8f52-4b1c-8da1-788101798f55
# ╠═dc0add50-5003-4e6e-a814-32ca84aaf187
# ╠═a0e7532b-2aaf-4fb3-af00-1d0912787c35
# ╠═5106a242-4fe5-47d8-94f7-45c12053c882
# ╠═21c57d80-3dda-412f-a912-0248a2f6d600
# ╠═e9ccfc08-12c7-43af-a372-7f8c04017462
# ╠═4ff87ef7-2132-4be9-9b71-15510fe3345a
# ╠═12177747-6d65-42d4-be9f-c4a86d37e65e
# ╠═2f97e7db-b2d5-44af-9c03-72e78f3edcde
# ╠═14ceb32c-19bb-44ec-b795-c807fad60525
# ╠═4bf130aa-7cc4-470b-bbe6-bdc86f831683
# ╠═c66ed3b6-4d1c-43ad-ac22-eced49adc811
# ╠═7d36f094-70bc-4114-b920-3655f73dc75f
# ╠═391d2ecc-424e-4e81-adc8-0d4b03ed9a40
# ╠═8f8eb577-5af9-4d08-ac50-77130a8fa880
# ╠═5e1388e9-1674-4fbc-865c-d1eb9ee74445
# ╠═7ecee6bf-69b3-4570-aa70-8491b96a46b3
# ╠═ec8c28f1-b6d6-4616-a6b3-302d8dae1d6d
# ╠═a12a2cfb-d94c-44ec-9a27-782e0120ac97
# ╠═8b287edc-bcbc-4ce2-b355-1d174bc9dc7a
# ╠═d800f480-3847-43fd-a458-aaaece1d9dcd
# ╠═a15dbc88-8dbf-48d4-8a72-1881151c8e26
# ╠═9821e893-507b-45c3-b833-6de6c7445ee1
# ╠═273d1aae-1b1e-42b8-85e4-efdcbcd04d06
# ╠═079fae18-072c-4a1c-b18b-016421d3efbe
# ╠═7d1e9712-1f02-4713-a1f1-a61cebe6d741
# ╠═1aa9e317-2f53-43bc-af7d-9294da42b607
# ╠═fa0b0f0b-4efe-43f4-b879-572bd3197dbb
# ╠═a09669c6-72fd-418d-957c-ffa2f6434c07
# ╠═8faefce1-f32a-4119-8e1c-9dcc13293da6
# ╠═d8e7b34f-4938-4ea1-a4c2-cecad0f3e5e5
# ╠═990172f0-6303-4cdb-8637-845dd60ce2e0
# ╠═6158b971-e5fb-4ef7-84e6-e0cae6072686
# ╠═c263e88e-58db-4273-943e-593e43070a02
# ╠═5433f056-2a84-4439-85f1-1f4786162f5d
# ╠═cca4b599-298e-4033-b9aa-36f9c5ba85e2
# ╠═ffcd3315-5e32-47cb-b831-7a19cb1c2162
# ╠═bcf8c8ae-cf45-4269-b6cf-c197147ff6ba
# ╠═bd0d8191-841c-4dcb-abf0-6f9f1728a761
# ╠═a136ee90-aad1-4bca-a40e-554ba921e90e
# ╠═f79760a9-8fb7-4249-994f-29f54f5c4c3c
# ╠═0555a067-17ac-43a4-82b8-440f0de638a6
# ╠═933e9188-ec1f-4ae8-81ef-7e51adff262a
# ╠═3426557d-074d-471d-b8b8-f18be83a5d78
# ╠═70a701b3-68e5-45d8-a7ba-c8bff151e8a6
# ╠═aa5f86a8-b35d-4849-af2c-b9d68fb53d77
# ╠═42f779b7-a1d1-4ed9-bd41-42a59143bc01
# ╠═0bda0872-7566-4199-a85e-673e241783dc
# ╠═56095498-3113-46b5-afb9-2e5e2d501ff1
# ╠═697b83ce-1787-484a-aa8a-12b91d112d88
# ╠═03a263f1-b2f2-43a9-96c8-7895eac5632b
