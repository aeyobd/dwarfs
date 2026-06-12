### A Pluto.jl notebook ###
# v1.0.1

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

# ╔═╡ aee2f674-d84c-4d6e-9536-c2d2cb1048ca
include("./utils.jl")

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
	"fiducial" => NFW(v_circ_max = 30 / V2KMS, r_circ_max = 2.2),
	"mean" => NFW(v_circ_max = 22 / V2KMS, r_circ_max = 3.9),
	"compact" => NFW(v_circ_max = 30 / V2KMS, r_circ_max = 1.0),
	"big" => NFW(v_circ_max = 30 / V2KMS, r_circ_max = 3),
	)

# ╔═╡ f6178d1d-3687-4538-80e7-4f2efe93ab54
obs_props_boo3 = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", "bootes3", "observed_properties.toml"))

# ╔═╡ 4ebadd1d-4b5e-4de9-987f-4cdb6b985952
md"""
# Utils
"""

# ╔═╡ b4bd8f82-ec5d-490b-8be2-08a70f311e24
⊕(x, y) = sqrt(x^2 + y^2)

# ╔═╡ 27425538-2104-46fe-8fc9-e6c28c0c6041
α_exp = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 0ec48582-0a73-4b57-87bc-24c19f79209a
R2r = LilGuys.r_h(LilGuys.Exp2D()) / LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 60d08767-3f2b-428c-a55e-50342a34191d
α_exp_cusp = LilGuys.r_circ_max(LilGuys.ExpCusp())

# ╔═╡ 21c57d80-3dda-412f-a912-0248a2f6d600
smallfontsize = @lift 0.8 * $(theme(:fontsize))

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
	
	σy_fattahi = (log_Ms_err / m_fattahi) ⊕ σ_fattahi

	return (;
		log_Ms = log_Ms,
		v_0 = v_0,
		log_v_0 = log_v_0,
		log_v_0_err = σy_fattahi,
		R_h = obs_props["R_h_kpc"],
	)
end

# ╔═╡ 4f4dd26d-3b93-49fe-9d9a-472f504e365d
der_props_boo3 = derived_properties(obs_props_boo3)

# ╔═╡ cc422906-8ec9-43ec-9c34-3ccad2509e74
der_props_boo3_10Ms = derived_properties(obs_props_boo3, dlog_Ms=1)

# ╔═╡ e2a2e7b0-fe92-4abc-9b52-842e07138bcc
10 ^(Measurement(der_props_boo3.log_v_0, der_props_boo3.log_v_0_err)) * V2KMS

# ╔═╡ 86e8d0fc-95b2-405d-87ce-dce0c8666641
σv_obs = obs_props_boo3["sigma_v"]

# ╔═╡ e9ccfc08-12c7-43af-a372-7f8c04017462
function plot_fattahi!(log_v_0, log_v_0_err; label="Fattahi+18", shade=false)
	
	v_0 = V2KMS * 10^log_v_0
	v_l, v_h = V2KMS .* 10 .^ ((log_v_0 - log_v_0_err), (log_v_0 + log_v_0_err))
	v_ll, v_hh = V2KMS .* 10 .^ ((log_v_0 - 2log_v_0_err), (log_v_0 + 2log_v_0_err))
	color = COLORS[3]

	if shade
		hspan!(v_l, v_h, xmax=1, color=(color, 0.1))
		# hspan!(v_ll, v_hh, xmax=1, color=(color, 0.1))
	end
	hlines!(v_0, color=(color, 0.5), xmax=1)
	text!(0.3, v_0, text=label, color=color, fontsize=smallfontsize)
end

# ╔═╡ 12177747-6d65-42d4-be9f-c4a86d37e65e
function plot_mass_concentration!()
	color = COLORS[8]

	vs = LinRange(10, 90, 1000) / V2KMS
	rs = LilGuys.Ludlow.solve_rmax.(vs)
	xl = (LilGuys.Ludlow.solve_rmax.(vs, -σ_ludlow))
	xh =  (LilGuys.Ludlow.solve_rmax.(vs, σ_ludlow))
	xll = (LilGuys.Ludlow.solve_rmax.(vs, -2σ_ludlow))
	xhh =  (LilGuys.Ludlow.solve_rmax.(vs, 2σ_ludlow))
	
	lines!(rs, vs*V2KMS, color=color)
	band!(vs*V2KMS, xl, xh, color=(color, 0.1), direction=:y)
	# band!(vs*V2KMS, xll, xhh, color=(color, 0.1), direction=:y)

	y0 = 24
	x0 = LilGuys.Ludlow.solve_rmax(y0 / V2KMS)	
	text_along_line_log!(rs, vs*V2KMS, x0,
						 text="Ludlow+16", color=color)
end

# ╔═╡ 3095f45a-5a66-45e7-a26b-27bf9b37c5d9
der_props_boo3_10Ms.log_v_0 .- der_props_boo3.log_v_0

# ╔═╡ 278c521a-062a-4760-a1e4-309ed291ee68
md"""
## Solutions given velocity
"""

# ╔═╡ e8b5d7d3-2c9c-4b00-b429-207e307ffae2
function σv_wolf(halo, R_h)
	return LilGuys.v_circ(halo, R_h*R2r)/sqrt(3)
end

# ╔═╡ f597d368-2f3d-4817-8131-6f933eeed76e
function predict_σv(rmax, vmax, orbit_props)
	r_end, v_end = Rapha.rapha_final_halo(rmax, vmax, orbit_props["pericentre"], orbit_props["apocentre"], orbit_props["n_peris"])

	r_s = r_end ./ α_exp_cusp
	v0 = LilGuys.v_circ_max.(LilGuys.ExpCusp.(1, r_s))
	M = @. (v_end/v0)^2
	halo_f = LilGuys.ExpCusp.(M, r_s)

	return σv_wolf(halo_f, der_props_boo3.R_h)
end

# ╔═╡ b25bff4a-39eb-48a9-99b8-9ff9e2a75974
function solve_vmax(rmax, σv, orbit_props)


	function f(log_vmax)
		vmax = 10 ^ log_vmax
		σv_pred = predict_σv(rmax, vmax, orbit_props)
		return σv_pred - σv
	end


	log_vmax =  LilGuys.find_zero(f, 0)
	10 ^ log_vmax
end

# ╔═╡ 70634c78-1ccf-47c0-b08c-cc6b713a6da2
function plot_vmax_solution!(σv, orbit_props; label="", kwargs...)
	
	rmax = logrange(0.1, 10, 1000) 
	vmax = solve_vmax.(rmax, σv, [orbit_props])
	l = lines!(rmax, vmax*V2KMS; kwargs...)
	
	x0 = 1.2
	text_along_line_log!(
		rmax, vmax*V2KMS, x0,
		text=label, color=l.color[], fontsize=smallfontsize
	)
	
end

# ╔═╡ 80013bbb-c03b-4ae9-8be5-4687fe58b7e3
md"""
# Plots
"""

# ╔═╡ be5a1a3d-8a2c-409e-b116-8f446abb3513
function load_orbit_props(modelname)
	orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "orbits/bootes3/EP2020_special_cases", modelname * ".toml"))
end

# ╔═╡ 5b35f67b-dbbe-4f4d-9b2b-a4903d844862
orbit_props = OrderedDict(
	"7kpc" => load_orbit_props("orbit_mean"),
	"12kpc" => load_orbit_props("orbit_peri_12"),
	"18kpc" => load_orbit_props("orbit_peri_18"),
	"26kpc" => load_orbit_props("orbit_peri_26"),
)

# ╔═╡ 74db2d15-d0d7-4b01-b398-74bafddf6876
styles = OrderedDict(
	"fiducial" => (;color=:black, marker=:star5),
	"mean" => (; color=COLORS[1], marker=:circle),
	"compact" => (; color=COLORS[4], marker=:utriangle),
	"big" => (; color=COLORS[2], marker=:rect)
)

# ╔═╡ 7d36f094-70bc-4114-b920-3655f73dc75f
function plot_halo_constraints(gs, der_props, halos)
	ax = Axis(gs,
		ylabel=L"$\,\textrm{v}_\textrm{max}$ / km\,s$^{-1}$",
		xlabel=L"$\,r_\textrm{max}$ / kpc",
		xscale=log10,
		yscale=log10,
		limits=(0.3, 10, 14, 40),
		xminorticks = [0.1:0.1:1; 1:10; 10:10:30],
			  xticks = [0.1, 1, 10],
		yminorticks = 10:2:80
	)


	plot_mass_concentration!()
	R_h = der_props.R_h

	
	annotation!(0, 30, der_props_boo3.R_h, 15, text=L"R_h", color=:black)
	
	plot_fattahi!(der_props_boo3_10Ms.log_v_0, der_props_boo3_10Ms.log_v_0_err, label="Fattahi+18", shade=true)
	arrows2d!([0.5], [V2KMS * 10^der_props_boo3.log_v_0], [1e-3], [V2KMS  * (10^der_props_boo3_10Ms.log_v_0 .- 10^der_props_boo3.log_v_0)], color=COLORS[3], minshaftlength=0, shaftwidth=1, tipwidth=5 / sqrt(3) * 2, tiplength=5,)
	text!(0.5, V2KMS * 10^der_props_boo3.log_v_0, align=(:center, :top), text=L"10\times M_\star", color=COLORS[3], fontsize=smallfontsize)
	
	for (label, halo) in halos
		y = LilGuys.v_circ_max(halo) * V2KMS
		x = LilGuys.r_circ_max(halo) 
		scatter!((x), (y), label=string(label); styles[label]...)
	end

	

	ax
end

# ╔═╡ 10423004-6c87-4bed-add3-80291ccccd67
# turn off scatter strokewidth cycling for arrows
theme(:Scatter)[:cycle] = Cycle([:marker] => :marker)

# ╔═╡ a15dbc88-8dbf-48d4-8a72-1881151c8e26
let
	fig = Figure(size=(6, 3) .* 72)

	ax = plot_halo_constraints(fig[1,1], der_props_boo3, halos_boo3)

	plot_vmax_solution!(σv_obs/V2KMS, orbit_props["26kpc"], label="26 kpc")
	plot_vmax_solution!(σv_obs/V2KMS, orbit_props["18kpc"], label="18 kpc")
	plot_vmax_solution!(σv_obs/V2KMS, orbit_props["12kpc"], label="12 kpc", color=:black)
	plot_vmax_solution!(σv_obs/V2KMS, orbit_props["7kpc"], label="7 kpc")

	ax.title[] = "pericentre"


	ax2 = plot_halo_constraints(fig[1,2], der_props_boo3, halos_boo3)

	df = copy(orbit_props["12kpc"])
	plot_vmax_solution!(σv_obs/V2KMS, df, label="mean", color=:black)
	plot_vmax_solution!((σv_obs - obs_props_boo3["sigma_v_em"])/V2KMS, df, label=L"low $\sigma_\text{v}$")
	plot_vmax_solution!((σv_obs + obs_props_boo3["sigma_v_ep"])/V2KMS, df, label=L"high $\sigma_\text{v}$")

	ax2.title[] = "velocity dispersion"
	hideydecorations!(ticks=false, minorticks=false)


	ax3 = plot_halo_constraints(fig[1,3], der_props_boo3, halos_boo3)

	df["n_peris"] = 1
	plot_vmax_solution!(σv_obs/V2KMS, df, label="1 peri")
	df["n_peris"] = 2
	plot_vmax_solution!(σv_obs/V2KMS, df, label="2 peri")
	df["n_peris"] = 4
	plot_vmax_solution!(σv_obs/V2KMS, df, label="4 peri")
	df["n_peris"] = 7
	plot_vmax_solution!(σv_obs/V2KMS, df, label="7 peri", color=:black)
	
	hideydecorations!(ticks=false, minorticks=false)
	ax3.title[] = "numer pericentres"

	
	Legend(fig[2, :], ax, tellwidth=false, tellheight=true, nbanks=4)
	@savefig "initial_halo_selection"
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
# ╠═aee2f674-d84c-4d6e-9536-c2d2cb1048ca
# ╟─3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═c515f500-27c5-403d-b50f-f0f25c7850cd
# ╠═d10fb6a6-6bd1-4748-9c7f-b69cccabe176
# ╠═22629f89-3d2a-457d-b80a-7df412b9c791
# ╠═f6178d1d-3687-4538-80e7-4f2efe93ab54
# ╟─4ebadd1d-4b5e-4de9-987f-4cdb6b985952
# ╠═b4bd8f82-ec5d-490b-8be2-08a70f311e24
# ╠═27425538-2104-46fe-8fc9-e6c28c0c6041
# ╠═0ec48582-0a73-4b57-87bc-24c19f79209a
# ╠═60d08767-3f2b-428c-a55e-50342a34191d
# ╠═21c57d80-3dda-412f-a912-0248a2f6d600
# ╠═be9e982a-904d-497c-949a-1fd265fcb67a
# ╠═f5794260-7b09-492b-a7e3-0f9559ce0d35
# ╠═4f4dd26d-3b93-49fe-9d9a-472f504e365d
# ╠═cc422906-8ec9-43ec-9c34-3ccad2509e74
# ╠═e2a2e7b0-fe92-4abc-9b52-842e07138bcc
# ╠═86e8d0fc-95b2-405d-87ce-dce0c8666641
# ╠═e9ccfc08-12c7-43af-a372-7f8c04017462
# ╠═12177747-6d65-42d4-be9f-c4a86d37e65e
# ╠═7d36f094-70bc-4114-b920-3655f73dc75f
# ╠═3095f45a-5a66-45e7-a26b-27bf9b37c5d9
# ╠═278c521a-062a-4760-a1e4-309ed291ee68
# ╠═e8b5d7d3-2c9c-4b00-b429-207e307ffae2
# ╠═70634c78-1ccf-47c0-b08c-cc6b713a6da2
# ╠═f597d368-2f3d-4817-8131-6f933eeed76e
# ╠═b25bff4a-39eb-48a9-99b8-9ff9e2a75974
# ╟─80013bbb-c03b-4ae9-8be5-4687fe58b7e3
# ╠═be5a1a3d-8a2c-409e-b116-8f446abb3513
# ╠═5b35f67b-dbbe-4f4d-9b2b-a4903d844862
# ╠═74db2d15-d0d7-4b01-b398-74bafddf6876
# ╠═10423004-6c87-4bed-add3-80291ccccd67
# ╠═a15dbc88-8dbf-48d4-8a72-1881151c8e26
