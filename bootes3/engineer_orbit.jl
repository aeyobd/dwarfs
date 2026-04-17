### A Pluto.jl notebook ###
# v0.20.24

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

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ 85fa9a66-7354-4cba-b6a9-dd88e86194ef
using Distributions

# ╔═╡ 7da8d805-3505-45f4-8205-c317d2415b9a
using DataFrames

# ╔═╡ d4e32b48-c95f-4cc2-910e-2a5f4bcc81bc
using PlutoUI

# ╔═╡ f653a17a-8062-4898-9d0c-168a13dfbf6f
using OrderedCollections

# ╔═╡ 9821e893-507b-45c3-b833-6de6c7445ee1
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl"); scale_theme_element!(:linewidth, 0.5)

# ╔═╡ 91100369-bd73-4bf1-a459-c569f9d976bc
@bind vmax_kms confirm(NumberField(10:0.1:40, default=30))

# ╔═╡ a3835522-a1b8-44de-bfe5-68ee780f8278
@bind delta_logc NumberField(0:0.1:5, default=2)

# ╔═╡ 69205eb1-d9c8-4213-88ec-759424bf922a
@bind pericentre confirm(NumberField(0:1:40, default=7))

# ╔═╡ 224d25b1-46c8-4d1e-b98a-862ddce604a4
dwarfs_dir = ENV["DWARFS_ROOT"]

# ╔═╡ 57e2cf1b-58c3-4135-9dd5-cc28d79ee9b6
include(joinpath(dwarfs_dir, "utils/pluto_utils.jl"))

# ╔═╡ fb5afd29-d1e8-4ec6-bbec-42362b8ddf49
include(joinpath(dwarfs_dir, "utils/rapha_utils.jl"))

# ╔═╡ ed8867fb-ca90-4577-905d-0ad93ec37097
import TOML

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

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
md"""
# Setup
"""

# ╔═╡ c515f500-27c5-403d-b50f-f0f25c7850cd
σ_fattahi = 0.035

# ╔═╡ d10fb6a6-6bd1-4748-9c7f-b69cccabe176
σ_ludlow = 0.10

# ╔═╡ 5e4d1dbb-96b5-4b80-8b28-244223bc7912
vmax = vmax_kms / V2KMS

# ╔═╡ ea9f5ca5-a77b-45a8-96cd-3ae6bfcb2284
rmax = LilGuys.Ludlow.solve_rmax(vmax, delta_logc * σ_ludlow)

# ╔═╡ 9c76311e-a16d-4b4d-af0b-0c9e05a77fc5
halo_boo3 = NFW(v_circ_max = vmax, r_circ_max = rmax)

# ╔═╡ f6178d1d-3687-4538-80e7-4f2efe93ab54
obs_props_boo3 = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", "bootes3", "observed_properties.toml"))

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

# ╔═╡ 8f8eb577-5af9-4d08-ac50-77130a8fa880
import Agama

# ╔═╡ 5e1388e9-1674-4fbc-865c-d1eb9ee74445
pot = Agama.Potential(file = joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))

# ╔═╡ 4946fb19-4b0e-41cb-b08f-4dce34162764
md"""
## Plot utils
"""

# ╔═╡ 21c57d80-3dda-412f-a912-0248a2f6d600
smallfontsize = @lift 0.8 * $(theme(:fontsize))

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

# ╔═╡ 4bf130aa-7cc4-470b-bbe6-bdc86f831683
LilGuys.mean_density(NFW(r_circ_max=1, v_circ_max=0.001), der_props_boo3.R_h*R2r)

# ╔═╡ 7d36f094-70bc-4114-b920-3655f73dc75f
function plot_halo_constraints(gs, der_props, halo)
	ax = Axis(gs,
		ylabel=L"$\,\textrm{v}_\textrm{max}$ / km\,s$^{-1}$",
		xlabel=L"$\,r_\textrm{max}$ / kpc",
		xscale=log10,
		yscale=log10,
		limits=(0.1, 30, 10, 80),
		xminorticks = [0.1:0.1:1; 1:10; 10:10:30],
		 xticks = [0.1, 1, 10],
		 yminorticks = [0.1:0.1:10; 10:1:80]
	)

	rf = rotation_factor(ax, true)


	plot_mass_concentration!(rf=rf)
	R_h = der_props.R_h
	plot_sigma_v!(7.7; x0 = 1, rf=rf, R_h=R_h)
	plot_sigma_v!(7.7-1.5; x0 = 1, rf=rf, R_h=R_h)
	plot_sigma_v!(7.7+2; x0 = 1, rf=rf, R_h=R_h)



	
	y = LilGuys.v_circ_max(halo) * V2KMS
	x = LilGuys.r_circ_max(halo) 
	scatter!((x), (y),)


	ax
end

# ╔═╡ 391d2ecc-424e-4e81-adc8-0d4b03ed9a40
function plot_tidal_track!(v0, r0)
	x, y = LilGuys.EN21_tidal_track(r0, v0, x_min=0.01)

	lines!(x, y*V2KMS, color=COLORS[5], alpha=0.5)
end

# ╔═╡ 49cf2db8-f534-4b42-8fde-90c24fa9a122
md"""
# Tidal track machinery
"""

# ╔═╡ 273d1aae-1b1e-42b8-85e4-efdcbcd04d06
orbit_props = read_fits(joinpath(ENV["DWARFS_ROOT"], "orbits/bootes3/EP2020/orbital_properties.fits"))

# ╔═╡ fd8a09ed-e136-47d7-95bf-16a305055cba
quantile(orbit_props.pericentre, cdf(Normal(), [-2, -1, 0, 1, 2, 3, 4, 5]))

# ╔═╡ 26e66dca-2a4e-442d-ba8b-4775c67ba10e
V0 = Agama.circular_velocity(pot, pericentre)

# ╔═╡ 427b60f6-31ab-4435-a73d-d7e2f312e6cb
 f_ecc_e(x) = (2*x / (x + 1))^3.2

# ╔═╡ 655df646-ac13-4bf3-83af-789197230f26
n_scale = 1.5

# ╔═╡ 56095498-3113-46b5-afb9-2e5e2d501ff1


# ╔═╡ 70de09ee-46ef-474e-b899-5c771601e45a
idx_orbit = argmin(abs.(orbit_props.pericentre .- pericentre))

# ╔═╡ 295f7c6a-d07f-4379-a6a8-605f4959f254
orbit = orbit_props[idx_orbit, :]

# ╔═╡ 576d5649-f7e1-41d3-874f-10fbdd830b4a
apocentre = orbit.apocentre

# ╔═╡ 259f7ed4-9f43-4c9f-90d9-35e44a881da8
n_peri_max = round(Int, -10 / T2GYR / orbit.period)

# ╔═╡ daee353a-2283-4194-b7bf-83a64e59786d
@bind n_peri NumberField(0:1:n_peri_max, default=n_peri_max)

# ╔═╡ bb75d21c-2e17-4781-89a0-4d074d8e4095
time_end = -n_peri * orbit.period

# ╔═╡ 9bca0ddb-7b5c-49b6-bf9b-bd81fcff70e1
time_end * T2GYR

# ╔═╡ 9960527a-e5fd-4368-9cfc-aaeb1a284795
let
	fig = Figure(size=(4, 3) .* 72)

	plot_halo_constraints(fig[1,1], der_props_boo3, halo_boo3)

	annotation!(0, 30, der_props_boo3.R_h, 10, text=L"R_h")

	ylims!(10, 40)

	for n in 0:n_peri_max
		r, v = rapha_final_halo(rmax, vmax, pericentre, apocentre, n)
		scatter!(r, v * V2KMS, color=COLORS[1],)
		text!(r, v*V2KMS, text="$n", align=(:left, :center), offset=(smallfontsize[], 0))
	end


	@savefig "initial_halo_selection"
	fig
end

# ╔═╡ ec2871bf-3ec9-4df8-adc1-5556584e1e12
rmax_end, vmax_end = rapha_final_halo(rmax, vmax, pericentre, apocentre, n_peri)

# ╔═╡ a2a49d4d-8f60-4479-bffb-7dd0d34eea08
rmax_end, vmax_end * V2KMS

# ╔═╡ 3ff176d9-c355-4470-bb7f-97dd9c47f0fe


# ╔═╡ 8f7c4d78-3437-4937-9111-779bc71cbc0a
md"""
# Samples
"""

# ╔═╡ 98334af7-fe6e-48f0-958a-9d60adb7ae95
samples = let
	df = DataFrame()


	N = 100_000

	df[!, :peri] = orbit_props.pericentre[1:N]
	df[!, :t_end] = 10 * rand(N) ./ T2GYR
	df[!, :period] = -orbit_props.period[1:N]
	df[!, :n_peri] = ceil.(df.t_end ./ df.period)

	vr_end = [rapha_final_halo(rmax, vmax, orbit_props.pericentre[i], orbit_props.apocentre[i], df.n_peri[i])
					for i in 1:N]

	df[!, :v_circ_end] = last.(vr_end)
	df[!, :r_circ_end] = first.(vr_end)


	r_s = df.r_circ_end ./ α_exp_cusp
	v0 = LilGuys.v_circ_max.(LilGuys.ExpCusp.(1, r_s))
	M = @. (df.v_circ_end/v0)^2
	halos_final_samples = LilGuys.ExpCusp.(M, r_s)

	
	df[!, :M_end] = M
	df[!, :r_s_end] = r_s
	df[!, :sigma_v] = σv_wolf.(halos_final_samples, der_props_boo3.R_h)

	
	df
end

# ╔═╡ 315ff7a6-4768-47fe-bd92-79843b1b2a64
md"""
# Plots
"""

# ╔═╡ a15dbc88-8dbf-48d4-8a72-1881151c8e26
let
	fig = Figure(size=(4, 3) .* 72)

	plot_halo_constraints(fig[1,1], der_props_boo3, halo_boo3)

	annotation!(0, 30, der_props_boo3.R_h, 10, text=L"R_h")

	plot_tidal_track!(vmax, rmax)
	ylims!(10, 40)

	scatter!(rmax_end, vmax_end * V2KMS)


	@savefig "initial_halo_selection"
	fig
end

# ╔═╡ 1a779e71-47f6-4bdc-a1ad-109ec45f2956
N_samples = size(samples, 1)

# ╔═╡ 7f141f24-c911-4fc6-9191-83a132451b46
filt_samples = obs_props_boo3["sigma_v"] - obs_props_boo3["sigma_v_em"] .< samples.sigma_v * V2KMS .< obs_props_boo3["sigma_v"] + obs_props_boo3["sigma_v_ep"]

# ╔═╡ e7e05e13-b65c-4b9d-ba84-94a6cf9b1387
let
	df = samples[filt_samples, :]

	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "pericentre",
			 ylabel= "# pericentres",
			 yticks=0:7, 
			 yminorticksvisible=false,)

	
	p = scatter!(df.peri, df.n_peri .+ 0.1*randn(size(df, 1)),
				 markersize=.5, alpha=1)

	for i in 1:maximum(df.n_peri)
		filt_n = (df.n_peri .== i ) 
		peri_med = median(df.peri[filt_n])
		peri_low, peri_high = quantile(df.peri[filt_n], [0.16, 0.84])
		errorscatter!([peri_med], [i], xerror =[(peri_med-peri_low, peri_high-peri_med)], color=:black, alpha=0.5)
		@info peri_med, peri_low, peri_high
	end

	text!(2, 7, text="rmax=$(round(rmax, digits=1)), vmax=$(round(Int, vmax_kms))", align=(:left, :top))
	
	

	xlims!(0, 30)
	ylims!(0.5, 7.5)

	fig
end

# ╔═╡ a079e891-3642-4120-9c5c-68962fe761e5
let
	df = samples[filt_samples, :]

	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "pericentre",
			 ylabel= "time / Gyr",
			 yminorticksvisible=false,)

	
	p = scatter!(df.peri, df.t_end * T2GYR,
				 markersize=.5, alpha=1)

	text!(0.05, 0.95, text="rmax=$(round(rmax, digits=1)), vmax=$(round(Int, vmax_kms))", align=(:left, :top), space=:relative)
	
	

	xlims!(0, 30)

	fig
end

# ╔═╡ Cell order:
# ╠═91100369-bd73-4bf1-a459-c569f9d976bc
# ╠═a3835522-a1b8-44de-bfe5-68ee780f8278
# ╠═69205eb1-d9c8-4213-88ec-759424bf922a
# ╠═daee353a-2283-4194-b7bf-83a64e59786d
# ╠═fd8a09ed-e136-47d7-95bf-16a305055cba
# ╠═bb75d21c-2e17-4781-89a0-4d074d8e4095
# ╠═9bca0ddb-7b5c-49b6-bf9b-bd81fcff70e1
# ╠═a2a49d4d-8f60-4479-bffb-7dd0d34eea08
# ╠═ea9f5ca5-a77b-45a8-96cd-3ae6bfcb2284
# ╠═9960527a-e5fd-4368-9cfc-aaeb1a284795
# ╠═85fa9a66-7354-4cba-b6a9-dd88e86194ef
# ╠═7da8d805-3505-45f4-8205-c317d2415b9a
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═224d25b1-46c8-4d1e-b98a-862ddce604a4
# ╠═57e2cf1b-58c3-4135-9dd5-cc28d79ee9b6
# ╠═fb5afd29-d1e8-4ec6-bbec-42362b8ddf49
# ╠═d4e32b48-c95f-4cc2-910e-2a5f4bcc81bc
# ╠═ed8867fb-ca90-4577-905d-0ad93ec37097
# ╠═f653a17a-8062-4898-9d0c-168a13dfbf6f
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═be9e982a-904d-497c-949a-1fd265fcb67a
# ╟─3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═c515f500-27c5-403d-b50f-f0f25c7850cd
# ╠═d10fb6a6-6bd1-4748-9c7f-b69cccabe176
# ╠═9c76311e-a16d-4b4d-af0b-0c9e05a77fc5
# ╠═5e4d1dbb-96b5-4b80-8b28-244223bc7912
# ╠═f6178d1d-3687-4538-80e7-4f2efe93ab54
# ╠═4f4dd26d-3b93-49fe-9d9a-472f504e365d
# ╠═f5794260-7b09-492b-a7e3-0f9559ce0d35
# ╠═1f74ca1b-ac53-425c-8f89-1afe361525cd
# ╠═38712f63-8f52-4b1c-8da1-788101798f55
# ╠═8f8eb577-5af9-4d08-ac50-77130a8fa880
# ╠═5e1388e9-1674-4fbc-865c-d1eb9ee74445
# ╟─4946fb19-4b0e-41cb-b08f-4dce34162764
# ╠═21c57d80-3dda-412f-a912-0248a2f6d600
# ╠═5106a242-4fe5-47d8-94f7-45c12053c882
# ╠═4ff87ef7-2132-4be9-9b71-15510fe3345a
# ╠═12177747-6d65-42d4-be9f-c4a86d37e65e
# ╠═2f97e7db-b2d5-44af-9c03-72e78f3edcde
# ╠═14ceb32c-19bb-44ec-b795-c807fad60525
# ╠═4bf130aa-7cc4-470b-bbe6-bdc86f831683
# ╠═7d36f094-70bc-4114-b920-3655f73dc75f
# ╠═391d2ecc-424e-4e81-adc8-0d4b03ed9a40
# ╟─49cf2db8-f534-4b42-8fde-90c24fa9a122
# ╠═9821e893-507b-45c3-b833-6de6c7445ee1
# ╠═273d1aae-1b1e-42b8-85e4-efdcbcd04d06
# ╠═26e66dca-2a4e-442d-ba8b-4775c67ba10e
# ╠═427b60f6-31ab-4435-a73d-d7e2f312e6cb
# ╠═655df646-ac13-4bf3-83af-789197230f26
# ╠═56095498-3113-46b5-afb9-2e5e2d501ff1
# ╠═70de09ee-46ef-474e-b899-5c771601e45a
# ╠═295f7c6a-d07f-4379-a6a8-605f4959f254
# ╠═576d5649-f7e1-41d3-874f-10fbdd830b4a
# ╠═259f7ed4-9f43-4c9f-90d9-35e44a881da8
# ╠═ec2871bf-3ec9-4df8-adc1-5556584e1e12
# ╠═3ff176d9-c355-4470-bb7f-97dd9c47f0fe
# ╠═8f7c4d78-3437-4937-9111-779bc71cbc0a
# ╠═98334af7-fe6e-48f0-958a-9d60adb7ae95
# ╟─315ff7a6-4768-47fe-bd92-79843b1b2a64
# ╠═a15dbc88-8dbf-48d4-8a72-1881151c8e26
# ╠═1a779e71-47f6-4bdc-a1ad-109ec45f2956
# ╠═7f141f24-c911-4fc6-9191-83a132451b46
# ╠═e7e05e13-b65c-4b9d-ba84-94a6cf9b1387
# ╠═a079e891-3642-4120-9c5c-68962fe761e5
