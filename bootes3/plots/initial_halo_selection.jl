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

# ╔═╡ 14c49a3a-efb2-449d-9efa-e7559768fb36
using PyFITS

# ╔═╡ 4cc61756-3ac3-4d10-9d81-38bacea93ede
using DataFrames, Distributions

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl"); scale_theme_element!(:linewidth, 0.5)

# ╔═╡ ed8867fb-ca90-4577-905d-0ad93ec37097
import TOML

# ╔═╡ 84ad44e5-f045-4166-b34a-a87929c6db25
module Rapha
	include(joinpath(ENV["DWARFS_ROOT"], "utils/rapha_utils.jl"))
end

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

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

# ╔═╡ b4bd8f82-ec5d-490b-8be2-08a70f311e24
⊕(x, y) = sqrt(x^2 + y^2)

# ╔═╡ 1e43350b-ebae-42a5-ad9e-778da7087732
dlog_Ms = 2

# ╔═╡ be9e982a-904d-497c-949a-1fd265fcb67a
"""
	get_slope(f, x, h)

Calculate the slope of the function f at x with stepsize h
"""
function get_slope(f, x0, h=1e-8x0)
	return (f(x0 + h) - f(x0)) / h
end

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
der_props_boo3_100Ms = derived_properties(obs_props_boo3, dlog_Ms=1)

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

# ╔═╡ 21c57d80-3dda-412f-a912-0248a2f6d600
smallfontsize = @lift 0.8 * $(theme(:fontsize))

# ╔═╡ e9ccfc08-12c7-43af-a372-7f8c04017462
function plot_fattahi!(log_v_0, log_v_0_err; label="Fattahi+18", shade=false)
	
	v_0 = V2KMS * 10^log_v_0
	v_l, v_h = V2KMS .* 10 .^ ((log_v_0 - log_v_0_err), (log_v_0 + log_v_0_err))
	v_ll, v_hh = V2KMS .* 10 .^ ((log_v_0 - 2log_v_0_err), (log_v_0 + 2log_v_0_err))
	color = COLORS[3]

	if shade
		hspan!(v_l, v_h, xmax=1, color=(color, 0.1))
		hspan!(v_ll, v_hh, xmax=1, color=(color, 0.1))
	end
	hlines!(v_0, color=(color, 0.5), xmax=1)
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
	xll = (LilGuys.Ludlow.solve_rmax.(vs, -2σ_ludlow))
	xhh =  (LilGuys.Ludlow.solve_rmax.(vs, 2σ_ludlow))
	
	lines!(rs, vs*V2KMS, color=color)
	band!(vs*V2KMS, xl, xh, color=(color, 0.1), direction=:y)
	band!(vs*V2KMS, xll, xhh, color=(color, 0.1), direction=:y)

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

# ╔═╡ 7d36f094-70bc-4114-b920-3655f73dc75f
function plot_halo_constraints(gs, der_props, halos)
	ax = Axis(gs,
		ylabel=L"$\,\textrm{v}_\textrm{max}$ / km\,s$^{-1}$",
		xlabel=L"$\,r_\textrm{max}$ / kpc",
		xscale=log10,
		yscale=log10,
		limits=(0.1, 30, 20, 80),
		xminorticks = [0.1:0.1:1; 1:10; 10:10:30],
			  xticks = [0.1, 1, 10],
		yminorticks = 20:1:80
	)

	rf = rotation_factor(ax, true)


	plot_mass_concentration!(rf=rf)
	R_h = der_props.R_h
	# plot_sigma_v!(7.7; x0 = 1, rf=rf, R_h=R_h)
	# plot_sigma_v!(7.7-1.5; x0 = 1, rf=rf, R_h=R_h)
	# plot_sigma_v!(7.7+2; x0 = 1, rf=rf, R_h=R_h)


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

# ╔═╡ b11526de-cd9d-44b5-98b5-e444cbeee825
let
	fig = Figure(size=(4, 3) .* 72)

	plot_halo_constraints(fig[1,1], der_props_boo3, halos_boo3)


	annotation!(0, 30, der_props_boo3.R_h, 15, text=L"R_h")
	plot_fattahi!(der_props_boo3_100Ms.log_v_0, der_props_boo3_100Ms.log_v_0_err, label=L"Fattahi+18, $10\times M_\star$", shade=true)

	ylims!(10, 40)



	scatter!(1.0, 30)

	scatter!(2.2, 30)
	scatter!(5.9, 30)

	
	@savefig "initial_halo_selection"
	fig
end

# ╔═╡ e8b3b478-9112-4c44-952e-c5bc1f91c148
orbit_props = read_fits(joinpath(ENV["DWARFS_ROOT"], "orbits/bootes3/EP2020/orbital_properties.fits"))

# ╔═╡ e8b5d7d3-2c9c-4b00-b429-207e307ffae2
function σv_wolf(halo, R_h)
	return LilGuys.v_circ(halo, R_h*R2r)/sqrt(3)
end

# ╔═╡ 435205e2-f0dc-4657-ac5a-68ef69b3e31f


# ╔═╡ ab13530c-9f82-42ac-b926-9fd83ba826b8
Rapha.rapha_final_halo(0.1, 10/V2KMS, 12, 100, 8)

# ╔═╡ 60d08767-3f2b-428c-a55e-50342a34191d
α_exp_cusp = LilGuys.r_circ_max(LilGuys.ExpCusp())

# ╔═╡ b25bff4a-39eb-48a9-99b8-9ff9e2a75974
function solve_vmax(rmax, σv, peri, apo, n_peri)


	function f(log_vmax)
		vmax = 10 ^ log_vmax
		r_end, v_end = Rapha.rapha_final_halo(rmax, vmax, peri, apo, n_peri)

		r_s = r_end ./ α_exp_cusp
		v0 = LilGuys.v_circ_max.(LilGuys.ExpCusp.(1, r_s))
		M = @. (v_end/v0)^2
		halo_f = LilGuys.ExpCusp.(M, r_s)

		return σv_wolf(halo_f, der_props_boo3.R_h) .- σv
	end


	local log_vmax
	try
		log_vmax =  LilGuys.find_zero(f, 0)
	catch e
		return NaN
	end

	10 ^ log_vmax
end

# ╔═╡ 90520423-cbf8-4c89-ad3f-2e78a2924bdd
function solve_rmax(vmax, σv, peri, apo, n_peri)


	function f(log_rmax)
		rmax = 10 ^ log_rmax
		r_end, v_end = Rapha.rapha_final_halo(rmax, vmax, peri, apo, n_peri)

		r_s = r_end ./ α_exp_cusp
		v0 = LilGuys.v_circ_max.(LilGuys.ExpCusp.(1, r_s))
		M = @. (v_end/v0)^2
		halo_f = LilGuys.ExpCusp.(M, r_s)

		return σv_wolf(halo_f, der_props_boo3.R_h) .- σv
	end


	local log_rmax
	try
		log_rmax =  LilGuys.find_zero(f, 0)
	catch e
		return NaN
	end

	10 ^ log_rmax
end

# ╔═╡ 65d79bf7-3569-4b41-b798-0877bcecf9b5
solve_rmax(10/V2KMS, 7/V2KMS, 12, 100, 8)

# ╔═╡ c5599cf2-beb0-46ba-b07e-5a9b8c4d0d75
solve_rmax(12/V2KMS, 7/V2KMS, 12, 100, 8)

# ╔═╡ 80ae3ce4-9a4e-4414-b0b1-4cd3136ef75f
let
	vmax, σv, peri, apo, n_peri = 20/V2KMS, 7/V2KMS, 4, 100, 8
	function f(log_rmax)
		rmax = 10 ^ log_rmax
		r_end, v_end = Rapha.rapha_final_halo(rmax, vmax, peri, apo, n_peri)

		r_s = r_end ./ α_exp_cusp
		v0 = LilGuys.v_circ_max.(LilGuys.ExpCusp.(1, r_s))
		M = @. (v_end/v0)^2
		halo_f = LilGuys.ExpCusp.(M, r_s)

		return σv_wolf(halo_f, der_props_boo3.R_h) .- σv
	end


	x = LinRange(-2, 2, 1000)
	lines(x, f.(x))
end

# ╔═╡ 9a2665d0-4995-4d13-93b2-7ea920e26972
let
	x = 0.9
	x / LilGuys.v_circ_EN21(x)
end

# ╔═╡ 3cb30f2d-eb1c-475a-b8a8-372d0880145f
lines(LinRange(0.5, 1, 100), x -> x/LilGuys.v_circ_EN21(x))

# ╔═╡ 70634c78-1ccf-47c0-b08c-cc6b713a6da2
function plot_vmax_solution!(σv, peri, apo, n_peri; kwargs...)

	rmax = logrange(0.1, 10, 1000) 
	vmax = solve_vmax.(rmax, σv, peri, apo, n_peri)

	lines!(rmax, vmax*V2KMS; kwargs...)


	x0 = 1.2
	y0 = solve_vmax.(x0, σv, peri, apo, n_peri) * V2KMS
	rf = rotation_factor(Makie.current_axis(), true)

	slope = log_derivative(x->solve_vmax.(x, σv, peri, apo, n_peri) * V2KMS, x0)
	θ = @lift atan(slope * $rf)


	text!(x0, y0, text="$peri kpc", rotation = θ)

	
end

# ╔═╡ 34cf4448-0656-43d1-bdc3-083ac6c69ccf
let
	fig = Figure(size=(4, 2.5) .* 72)

	plot_halo_constraints(fig[1,1], der_props_boo3, halos_boo3)


	annotation!(0, 30, der_props_boo3.R_h, 15, text=L"R_h")
	plot_fattahi!(der_props_boo3_100Ms.log_v_0, der_props_boo3_100Ms.log_v_0_err, label=L"Fattahi+18, $10\times M_\star$", shade=true)

	ylims!(10, 40)
	xlims!(0.3, 10)


	plot_vmax_solution!(7/V2KMS, 18, 150, 5)
	plot_vmax_solution!(5/V2KMS, 18, 150, 5)
	plot_vmax_solution!(9/V2KMS, 18, 150, 5)


	fig
end

# ╔═╡ a15dbc88-8dbf-48d4-8a72-1881151c8e26
let
	fig = Figure(size=(4, 2.5) .* 72)

	plot_halo_constraints(fig[1,1], der_props_boo3, halos_boo3)


	annotation!(0, 30, der_props_boo3.R_h, 15, text=L"R_h")
	plot_fattahi!(der_props_boo3_100Ms.log_v_0, der_props_boo3_100Ms.log_v_0_err, label=L"Fattahi+18, $10\times M_\star$", shade=true)

	ylims!(10, 40)
	xlims!(0.3, 10)


	plot_vmax_solution!(7/V2KMS, 18, 150, 5)
	plot_vmax_solution!(7/V2KMS, 12, 120, 7)
	plot_vmax_solution!(7/V2KMS, 7, 98, 8)

	# plot_vmax_solution!(5/V2KMS, 18, 150, 2)
	# plot_vmax_solution!(5/V2KMS, 12, 120, 2)
	# plot_vmax_solution!(5/V2KMS, 7, 98, 2)

	scatter!(3.0, 30)
	scatter!(2.2, 30)
	scatter!(1.0, 30)



	plot_halo_constraints(fig[1,2], der_props_boo3, halos_boo3)


	annotation!(0, 30, der_props_boo3.R_h, 15, text=L"R_h")
	plot_fattahi!(der_props_boo3_100Ms.log_v_0, der_props_boo3_100Ms.log_v_0_err, label=L"Fattahi+18, $10\times M_\star$", shade=true)

	ylims!(10, 40)
	xlims!(0.3, 10)


	plot_vmax_solution!(7/V2KMS, 7, 98, 2)
	plot_vmax_solution!(7/V2KMS, 4, 98, 2)
	plot_vmax_solution!(7/V2KMS, 1.5, 98, 2)


	scatter!(3.0, 30)

	hideydecorations!(ticks=false, minorticks=false)

	
	@savefig "initial_halo_selection"
	fig
end

# ╔═╡ 9f1cbc12-8592-4f14-b03f-4a046b7417f2
function plot_rmax_solution(σv, peri, apo, n_peri; kwargs...)

	vmax = LinRange(10, 40, 1000) ./ V2KMS
	rmax = solve_rmax.(vmax, σv, peri, apo, n_peri)

	lines!(rmax, vmax*V2KMS; kwargs...)


	y0 = 20
	x0 = solve_rmax.(y0/V2KMS, σv, peri, apo, n_peri)
	rf = rotation_factor(Makie.current_axis(), true)

	slope = 1/log_derivative(y->solve_rmax.(y/V2KMS, σv, peri, apo, n_peri), y0)
	θ = @lift atan(slope * $rf)


	text!(x0, y0, text="$peri kpc", rotation = θ)

	
end

# ╔═╡ 295979bb-350a-4908-bd93-ce08964ac540


# ╔═╡ bf26d49e-9461-4710-bc86-10653f0319a6
let
	fig = Figure(size=(4, 3) .* 72)

	plot_halo_constraints(fig[1,1], der_props_boo3, halos_boo3)


	annotation!(0, 30, der_props_boo3.R_h, 10, text=L"R_h")
	plot_fattahi!(der_props_boo3_100Ms.log_v_0, der_props_boo3_100Ms.log_v_0_err, label="100x M⋆", shade=true)

	ylims!(10, 40)

	plot_sigma_v!(15; x0 = 1, R_h=der_props_boo3.R_h)



	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═ed8867fb-ca90-4577-905d-0ad93ec37097
# ╠═f653a17a-8062-4898-9d0c-168a13dfbf6f
# ╠═14c49a3a-efb2-449d-9efa-e7559768fb36
# ╠═4cc61756-3ac3-4d10-9d81-38bacea93ede
# ╠═84ad44e5-f045-4166-b34a-a87929c6db25
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╟─3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═c515f500-27c5-403d-b50f-f0f25c7850cd
# ╠═d10fb6a6-6bd1-4748-9c7f-b69cccabe176
# ╠═22629f89-3d2a-457d-b80a-7df412b9c791
# ╠═f6178d1d-3687-4538-80e7-4f2efe93ab54
# ╠═b4bd8f82-ec5d-490b-8be2-08a70f311e24
# ╠═1e43350b-ebae-42a5-ad9e-778da7087732
# ╠═be9e982a-904d-497c-949a-1fd265fcb67a
# ╠═f5794260-7b09-492b-a7e3-0f9559ce0d35
# ╠═4f4dd26d-3b93-49fe-9d9a-472f504e365d
# ╠═cc422906-8ec9-43ec-9c34-3ccad2509e74
# ╠═e2a2e7b0-fe92-4abc-9b52-842e07138bcc
# ╠═27425538-2104-46fe-8fc9-e6c28c0c6041
# ╠═0ec48582-0a73-4b57-87bc-24c19f79209a
# ╠═1f74ca1b-ac53-425c-8f89-1afe361525cd
# ╠═38712f63-8f52-4b1c-8da1-788101798f55
# ╠═5106a242-4fe5-47d8-94f7-45c12053c882
# ╠═21c57d80-3dda-412f-a912-0248a2f6d600
# ╠═e9ccfc08-12c7-43af-a372-7f8c04017462
# ╠═4ff87ef7-2132-4be9-9b71-15510fe3345a
# ╠═12177747-6d65-42d4-be9f-c4a86d37e65e
# ╠═2f97e7db-b2d5-44af-9c03-72e78f3edcde
# ╠═14ceb32c-19bb-44ec-b795-c807fad60525
# ╠═4bf130aa-7cc4-470b-bbe6-bdc86f831683
# ╠═7d36f094-70bc-4114-b920-3655f73dc75f
# ╠═391d2ecc-424e-4e81-adc8-0d4b03ed9a40
# ╠═8f8eb577-5af9-4d08-ac50-77130a8fa880
# ╠═5e1388e9-1674-4fbc-865c-d1eb9ee74445
# ╠═7ecee6bf-69b3-4570-aa70-8491b96a46b3
# ╠═ec8c28f1-b6d6-4616-a6b3-302d8dae1d6d
# ╠═a12a2cfb-d94c-44ec-9a27-782e0120ac97
# ╠═8b287edc-bcbc-4ce2-b355-1d174bc9dc7a
# ╠═d800f480-3847-43fd-a458-aaaece1d9dcd
# ╠═34cf4448-0656-43d1-bdc3-083ac6c69ccf
# ╠═a15dbc88-8dbf-48d4-8a72-1881151c8e26
# ╠═b11526de-cd9d-44b5-98b5-e444cbeee825
# ╠═e8b3b478-9112-4c44-952e-c5bc1f91c148
# ╠═e8b5d7d3-2c9c-4b00-b429-207e307ffae2
# ╠═b25bff4a-39eb-48a9-99b8-9ff9e2a75974
# ╠═90520423-cbf8-4c89-ad3f-2e78a2924bdd
# ╠═435205e2-f0dc-4657-ac5a-68ef69b3e31f
# ╠═ab13530c-9f82-42ac-b926-9fd83ba826b8
# ╠═65d79bf7-3569-4b41-b798-0877bcecf9b5
# ╠═80ae3ce4-9a4e-4414-b0b1-4cd3136ef75f
# ╠═c5599cf2-beb0-46ba-b07e-5a9b8c4d0d75
# ╠═60d08767-3f2b-428c-a55e-50342a34191d
# ╠═9a2665d0-4995-4d13-93b2-7ea920e26972
# ╠═3cb30f2d-eb1c-475a-b8a8-372d0880145f
# ╠═70634c78-1ccf-47c0-b08c-cc6b713a6da2
# ╠═9f1cbc12-8592-4f14-b03f-4a046b7417f2
# ╠═295979bb-350a-4908-bd93-ce08964ac540
# ╠═bf26d49e-9461-4710-bc86-10653f0319a6
