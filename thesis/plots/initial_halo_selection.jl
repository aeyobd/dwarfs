### A Pluto.jl notebook ###
# v0.20.17

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

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ eb4a5051-43b6-4afa-bf42-5a74125cf60a
galaxy = "ursa_minor"

# ╔═╡ ed8867fb-ca90-4577-905d-0ad93ec37097
import TOML

# ╔═╡ 325ff1ea-ddfa-4661-9dd0-22866afca67d
import LinearAlgebra: eigen, diagm

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ b7a4e67b-a3c9-44c2-acf3-b61f60d387ad
scale_theme_element!(:linewidth, 0.5)

# ╔═╡ 48289b64-08a8-4b26-9506-1843b3c69009
scale_theme_element!(:markersize, 2)

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

# ╔═╡ 22629f89-3d2a-457d-b80a-7df412b9c791
halos_scl = OrderedDict(
		"Scl: fiducial" => NFW(v_circ_max = 31 / V2KMS, r_circ_max = 3.2),
		"Scl: small" => NFW(v_circ_max=25/V2KMS, r_circ_max=2.5),
	)

# ╔═╡ 61aa6484-ae3c-4478-a426-51e78ae02ed7
halos_umi = OrderedDict(
		"UMi: fiducial" =>  NFW(v_circ_max = 38 / V2KMS, r_circ_max=4.0)
	)

# ╔═╡ f6178d1d-3687-4538-80e7-4f2efe93ab54
obs_props_scl = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", "sculptor", "observed_properties.toml"))

# ╔═╡ 847d94fa-ce1e-4a3a-8f15-51f47b53336f
obs_props_umi = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", "ursa_minor", "observed_properties.toml"))

# ╔═╡ b4bd8f82-ec5d-490b-8be2-08a70f311e24
⊕(x, y) = sqrt(x^2 + y^2)

# ╔═╡ f5794260-7b09-492b-a7e3-0f9559ce0d35
function derived_properties(obs_props)
	log_Ms = log10(LilGuys.mag_to_L(obs_props["Mv"]) * obs_props["M_L_s"])
	log_Ms_err = LilGuys.get_uncertainty(obs_props, "M_L_s") / obs_props["M_L_s"] / log(10) ⊕ LilGuys.get_uncertainty(obs_props, "Mv") * 2/5

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
der_props_scl = derived_properties(obs_props_scl)

# ╔═╡ b76167d8-0c5a-47c1-9cdb-4ca6a16e0651
der_props_umi = derived_properties(obs_props_umi)

# ╔═╡ 27425538-2104-46fe-8fc9-e6c28c0c6041
α_exp = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 1f74ca1b-ac53-425c-8f89-1afe361525cd
function solve_v_from_dispersion(r0, σv, R_h)
	halo(v) = NFW(r_circ_max=r0, v_circ_max=v)
	stars = LilGuys.Exp2D(R_s=R_h/α_exp)

	f(v) = LilGuys.σv_star_mean(halo(v), stars) - σv/V2KMS
	v0 = LilGuys.find_zero(f, [0.01, 0.5])
	return v0
end

# ╔═╡ 273ccadb-04c6-485e-90e5-20cb4935f58d
solve_v_from_dispersion(1, 2, 0.24)

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
function plot_fattahi!(log_v_0, log_v_0_err; label="Fattahi+18")
	
	v_0 = V2KMS * 10^log_v_0
	v_l, v_h = V2KMS .* 10 .^ ((log_v_0 - log_v_0_err), (log_v_0 + log_v_0_err))
	color = COLORS[7]

	hspan!(v_l, v_h, color=(color, 0.1))
	hlines!(v_0, color=(color, 0.5))
	text!(1.1, v_0, text=label, color=color, fontsize=smallfontsize)
end

# ╔═╡ 4ff87ef7-2132-4be9-9b71-15510fe3345a
function log_derivative(f, x0; h=0.001)
	y0 = f(x0)

	return (log10(f(x0 * 10^h)) - log10(y0)) / h
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

# ╔═╡ 53967222-f243-4529-9544-caf2d68215c2
md"""
# Sanity checks
"""

# ╔═╡ 80005476-d236-4a62-bfc0-2e7d8ca914ff
stars = LilGuys.Exp2D(R_s = R_h/α_exp)

# ╔═╡ e0e418d2-16f2-4bbe-92db-ccdc6e89669e
[label => LilGuys.σv_star_mean(halo, stars) * V2KMS for (label, halo) in halos_scl]

# ╔═╡ f29e7c86-9fe2-4f9b-9d03-9a029773123e
R2r = LilGuys.r_h(LilGuys.Exp2D()) / LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 38712f63-8f52-4b1c-8da1-788101798f55
function solve_v_from_dispersion_simple(r0, σv, R_h)
	f(v) = LilGuys.v_circ(NFW(r_circ_max=r0, v_circ_max=v), R_h*R2r)/sqrt(3) - σv/V2KMS
	
	return LilGuys.find_zero(f, [0.01, 0.5])
end

# ╔═╡ 5555b64e-c234-4b32-b8dd-965e3842671c
solve_v_from_dispersion_simple(1, 2, 0.24)

# ╔═╡ 5106a242-4fe5-47d8-94f7-45c12053c882
function plot_sigma_v!(σv; R_h, x0=15, rf=1, vmax = solve_v_from_dispersion_simple, kwargs...)
	color = COLORS[5]
	
	x = LinRange(1, 30, 100)
	y = vmax.(x, σv, R_h) * V2KMS
	lines!(x, y; color=color, kwargs...)

	y0 = vmax(x0, σv, R_h) * V2KMS
	
	slope = log_derivative(x->vmax(x, σv, R_h) * V2KMS, x0)
	θ = @lift atan(slope * $rf)
	
	text!(x0, y0, text=L"$\sigma_\textrm{v} = %$σv\,$km\,s$^{-1}$", rotation=θ, color=color, fontsize=smallfontsize)
end

# ╔═╡ 7d36f094-70bc-4114-b920-3655f73dc75f
function plot_halo_constraints(gs, der_props, halos)
	ax = Axis(gs,
		ylabel=L"$\,\textrm{v}_\textrm{max}$ / km\,s$^{-1}$",
		xlabel=L"$\,r_\textrm{max}$ / kpc",
		xscale=log10,
		yscale=log10,
		limits=(1, 30, 12, 80),
		xminorticks = [1:10; 10:10:30]
	)

	rf = rotation_factor(ax, true)


	plot_mass_concentration!(rf=rf)
	R_h = der_props.R_h
	plot_sigma_v!(12; x0=11.8, rf=rf, linestyle=:dot, R_h=R_h)
	plot_sigma_v!(10; x0=13, rf=rf, linestyle=:dash, R_h=R_h)
	plot_sigma_v!(8; rf=rf, R_h=R_h)


	plot_fattahi!(der_props.log_v_0, der_props.log_v_0_err, label="Fattahi+18")

	
	for (label, halo) in halos
		y = LilGuys.v_circ_max(halo) * V2KMS
		x = LilGuys.r_circ_max(halo) 
		scatter!((x), (y), label=string(label))
	end

	
	axislegend(position=:rb)

	ax
end

# ╔═╡ a15dbc88-8dbf-48d4-8a72-1881151c8e26
let
	fig = Figure(size=(5, 3) .* 72)

	ax1 = plot_halo_constraints(fig[1,1], der_props_scl, halos_scl)
	ax2 = plot_halo_constraints(fig[1,2], der_props_umi, halos_umi)

	linkaxes!(ax1, ax2)
	hideydecorations!(ax2, ticks=false, minorticks=false)
	colgap!(fig.layout, 0.0)
	
	@savefig "initial_halo_selection"
	fig
end

# ╔═╡ 0ca4187b-d63c-48e4-947c-c9a48973323d
[LilGuys.v_circ(halo, R_h * R2r) * V2KMS / sqrt(3) for (label, halo) in halos_scl]

# ╔═╡ 51075ee8-2cd7-4643-bf45-24a8a60d4eec
let

	fig = Figure()
	ax = Axis(fig[1,1])

	for (label, halo) in merge(halos_scl, halos_umi)
		x = LinRange(-2, 2, 10_000)
		y = @. LilGuys.potential(halo, 10^x) 

		lines!(x, y)
	end

	fig

end

# ╔═╡ 9fb0dfe4-5d49-4b05-899e-cb4138dc7aac
md"""
These values are from the MCMC halo mass estimator notebook in `analysis/sculptor` 
"""

# ╔═╡ f178ab03-06ba-4e83-9be9-3924543bd288
Σ_scl = [0.019467512114446032 0.0019034735328083907; 0.0019034735328083907 0.001516776581639168]

# ╔═╡ 28f7a671-8b64-43ed-b079-95cd7360a73d
μ_scl = 
[0.780442
1.49899]

# ╔═╡ Cell order:
# ╠═eb4a5051-43b6-4afa-bf42-5a74125cf60a
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═ed8867fb-ca90-4577-905d-0ad93ec37097
# ╠═f653a17a-8062-4898-9d0c-168a13dfbf6f
# ╠═325ff1ea-ddfa-4661-9dd0-22866afca67d
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═b7a4e67b-a3c9-44c2-acf3-b61f60d387ad
# ╠═48289b64-08a8-4b26-9506-1843b3c69009
# ╟─be9e982a-904d-497c-949a-1fd265fcb67a
# ╟─3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═c515f500-27c5-403d-b50f-f0f25c7850cd
# ╠═d10fb6a6-6bd1-4748-9c7f-b69cccabe176
# ╠═22629f89-3d2a-457d-b80a-7df412b9c791
# ╠═61aa6484-ae3c-4478-a426-51e78ae02ed7
# ╠═f6178d1d-3687-4538-80e7-4f2efe93ab54
# ╠═847d94fa-ce1e-4a3a-8f15-51f47b53336f
# ╠═b4bd8f82-ec5d-490b-8be2-08a70f311e24
# ╠═f5794260-7b09-492b-a7e3-0f9559ce0d35
# ╠═4f4dd26d-3b93-49fe-9d9a-472f504e365d
# ╠═b76167d8-0c5a-47c1-9cdb-4ca6a16e0651
# ╠═27425538-2104-46fe-8fc9-e6c28c0c6041
# ╠═1f74ca1b-ac53-425c-8f89-1afe361525cd
# ╠═38712f63-8f52-4b1c-8da1-788101798f55
# ╠═273ccadb-04c6-485e-90e5-20cb4935f58d
# ╠═5555b64e-c234-4b32-b8dd-965e3842671c
# ╠═a0e7532b-2aaf-4fb3-af00-1d0912787c35
# ╠═5106a242-4fe5-47d8-94f7-45c12053c882
# ╠═21c57d80-3dda-412f-a912-0248a2f6d600
# ╠═e9ccfc08-12c7-43af-a372-7f8c04017462
# ╠═4ff87ef7-2132-4be9-9b71-15510fe3345a
# ╠═12177747-6d65-42d4-be9f-c4a86d37e65e
# ╠═2f97e7db-b2d5-44af-9c03-72e78f3edcde
# ╠═7d36f094-70bc-4114-b920-3655f73dc75f
# ╠═a15dbc88-8dbf-48d4-8a72-1881151c8e26
# ╟─53967222-f243-4529-9544-caf2d68215c2
# ╠═0ca4187b-d63c-48e4-947c-c9a48973323d
# ╠═e0e418d2-16f2-4bbe-92db-ccdc6e89669e
# ╠═80005476-d236-4a62-bfc0-2e7d8ca914ff
# ╠═f29e7c86-9fe2-4f9b-9d03-9a029773123e
# ╠═51075ee8-2cd7-4643-bf45-24a8a60d4eec
# ╟─9fb0dfe4-5d49-4b05-899e-cb4138dc7aac
# ╠═f178ab03-06ba-4e83-9be9-3924543bd288
# ╠═28f7a671-8b64-43ed-b079-95cd7360a73d
