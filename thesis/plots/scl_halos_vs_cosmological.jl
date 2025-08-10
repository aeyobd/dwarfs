### A Pluto.jl notebook ###
# v0.20.13

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
halos =
	if galaxy == "sculptor" 
		OrderedDict(
			:compact => NFW(v_circ_max = 31 / V2KMS, r_circ_max = 3.2),
			:lmc => NFW(v_circ_max = 31 / V2KMS, r_circ_max = 4.2),
			:small => NFW(v_circ_max=25/V2KMS, r_circ_max=2.5)
		)
	elseif galaxy == "ursa_minor"

		OrderedDict(
			:fiducial => NFW(v_circ_max = 38 / V2KMS, r_circ_max=4.0)
		)
	end

# ╔═╡ f6178d1d-3687-4538-80e7-4f2efe93ab54
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxy, "observed_properties.toml"))

# ╔═╡ b4bd8f82-ec5d-490b-8be2-08a70f311e24
⊕(x, y) = sqrt(x^2 + y^2)

# ╔═╡ 407d62b9-8e5c-4052-b281-e74b0a3290c4
log_Ms_0 = log10(LilGuys.mag_to_L(obs_props["Mv"]) * obs_props["M_L_s"])

# ╔═╡ c8454dc7-8aaf-40d8-9231-2e6e4f5f2186
log_Ms_0_err = LilGuys.get_uncertainty(obs_props, "M_L_s") / obs_props["M_L_s"] / log(10) ⊕ LilGuys.get_uncertainty(obs_props, "Mv") * 2/5

# ╔═╡ d1636c10-37d3-4c92-ad9c-e58a19016469
v_0 = LilGuys.find_zero(x -> log10(LilGuys.M_s_from_vel_fattahi(x)) + 10 - log_Ms_0, [0.01, 0.3])

# ╔═╡ b575403d-b2d3-4c81-9957-8f384c9811ac
log_v_0 = log10(v_0)

# ╔═╡ 0020549b-40f7-4542-8cd8-aad0e047a689
m_ludlow = 1/  get_slope(x -> log10(LilGuys.Ludlow.solve_rmax(10 ^ x)),log_v_0)

# ╔═╡ f61ca53d-ef76-47ce-bd2f-015645048d44
m_fattahi = get_slope(x -> log10(LilGuys.M_s_from_vel_fattahi(10^x)), log_v_0)

# ╔═╡ 25e17fae-4ae6-4a3d-99c9-646e3dfb574d
σy_fattahi = (log_Ms_0_err / m_fattahi) ⊕ σ_fattahi

# ╔═╡ 816c93e1-3f6a-4b76-9003-2622e8d1f9fc
σx_ludlow = 1/2 * (log10(LilGuys.Ludlow.solve_rmax(v_0, -σ_ludlow)) - log10(LilGuys.Ludlow.solve_rmax(v_0, σ_ludlow)))

# ╔═╡ 17dd04ed-2c60-488e-9b73-a8adea7dc9ff
σy_ludlow = σx_ludlow * m_ludlow

# ╔═╡ 73ab0963-1d9a-4488-be1c-5c061b2e0e7e
Σxx = (σy_fattahi^2 + σy_ludlow^2) / m_ludlow^2

# ╔═╡ a8e71e9a-cb26-4058-948b-0d27b60c5b48
Σyy = (σy_fattahi^2) 

# ╔═╡ 1e7d486e-7d0f-4f51-9175-a1a378774f1e
Σxy = (σy_fattahi^2) / m_ludlow

# ╔═╡ e29ca4df-c319-4e8d-9dfa-93ac7bf28ee8
Σ = [
	Σxx Σxy
	Σxy Σyy]

# ╔═╡ 4f4dd26d-3b93-49fe-9d9a-472f504e365d
log_r_0 = log10(LilGuys.Ludlow.solve_rmax(v_0))

# ╔═╡ 23959023-9412-4899-a546-60fecdf6fe6d
μ = [log_r_0, log_v_0 + log10(V2KMS)]

# ╔═╡ 2967e69c-8732-49eb-bdb2-91161993a5bb
σv_min = 9

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

# ╔═╡ 104e4af5-cafc-4de0-82b4-e3726bc7fc6e
R_s = LilGuys.arcmin2kpc(obs_props["R_h"], obs_props["distance"]) / LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ b01d7319-2f22-49a5-9303-2ef28b24b798
stars = LilGuys.Exp2D(R_s=R_s)

# ╔═╡ 1f74ca1b-ac53-425c-8f89-1afe361525cd
function solve_v_from_dispersion(r0, stars=stars)
	return LilGuys.find_zero(x -> LilGuys.σv_star_mean(NFW(r_circ_max=r0, v_circ_max=x), stars) - σv_min/V2KMS, [0.01, 0.5])
end

# ╔═╡ b013d48a-7d33-4630-a694-60c7562f3d52
solve_v_from_dispersion(10)

# ╔═╡ 3909983c-1a96-4fd5-9705-95cb4e9a184c
galaxy_short = Dict(
	"sculptor" => "scl",
	"ursa_minor" => "umi",
)[galaxy]

# ╔═╡ a15dbc88-8dbf-48d4-8a72-1881151c8e26
let
	fig = Figure(size=(5, 3) .* 72)
	ax = Axis(fig[1,1],
		ylabel=L"$\,v_\textrm{max}$ / km\,s$^{-1}$",
		xlabel=L"$\,r_\textrm{max}$ / kpc",
		xscale=log10,
		yscale=log10,
		limits=(1, 30, 12, 80),
		xminorticks = [1:10; 10:10:30]
	)

	vs = LinRange(10, 90, 1000) / V2KMS
	

	# ludlow constraint
	color = COLORS[8]
	rs = LilGuys.Ludlow.solve_rmax.(vs)
	xl = (LilGuys.Ludlow.solve_rmax.(vs, -σ_ludlow))
	xh =  (LilGuys.Ludlow.solve_rmax.(vs, σ_ludlow))
	lines!(rs, vs*V2KMS, color=color)
	band!(vs*V2KMS, xl, xh, color=(color, 0.1), direction=:y)

	text!(LilGuys.Ludlow.solve_rmax(50 / V2KMS), 50, text="Ludlow+16", color=color, rotation = π/5)

	color = COLORS[7]

	hspan!(10^(μ[2] - sqrt(Σ[2,2])), 10^(μ[2] + sqrt(Σ[2,2])), color=(color, 0.1))
	hlines!(10^(μ[2]), color=(color, 0.5))
	text!(1.1, v_0 * V2KMS, text="Fattahi+18", color=color)



	# velocity dispersion constraint
	color = COLORS[5]
	x = LinRange(1, 30, 100)
	y = solve_v_from_dispersion.(x) * V2KMS
	lines!(x, y, color=color)
	band!(x, y, y .+ 100, color=(color, 0.1))
	text!(x[70], y[70], text=L"\sigma_v > %$σv_min", rotation=0.3, color=color)
	
	for (label, halo) in halos
		y = LilGuys.v_circ_max(halo) * V2KMS
		x = LilGuys.r_circ_max(halo) 
		scatter!((x), (y), label=string(label))
	end

	plot_error_ellipses!(μ, Σ, color=(COLORS[4]))
	
	axislegend(position=:rb, title="halo")

	@savefig "$(galaxy_short)_initial_halos"
	
	fig
end

# ╔═╡ 5d0a1520-e8e0-4b90-b83b-f8605ed0cbb5
COLORS[8]

# ╔═╡ 679a471a-da30-4ef4-b8f7-60a69c8d8d94
R_h = 0.15

# ╔═╡ 0ca4187b-d63c-48e4-947c-c9a48973323d
[LilGuys.v_circ(halo, R_h) * V2KMS / sqrt(2) for (label, halo) in halos]

# ╔═╡ e0e418d2-16f2-4bbe-92db-ccdc6e89669e
[label => LilGuys.σv_star_mean(halo, stars) * V2KMS for (label, halo) in halos]

# ╔═╡ a085b24c-9ac6-41c8-811e-8e972eef0efb
let

	fig = Figure()
	ax = Axis(fig[1,1])

	for (label, halo) in halos
		x = LinRange(-2, 2, 10_000)
		y = @. LilGuys.v_circ(halo, 10^x) * V2KMS

		lines!(x, y, label=string(label))
	end

	axislegend(position=:lt)

	fig

end

# ╔═╡ 51075ee8-2cd7-4643-bf45-24a8a60d4eec
let

	fig = Figure()
	ax = Axis(fig[1,1])

	for (label, halo) in halos
		x = LinRange(-2, 2, 10_000)
		y = @. LilGuys.potential(halo, 10^x) 

		lines!(x, y)
	end

	fig

end

# ╔═╡ 53967222-f243-4529-9544-caf2d68215c2
md"""
# Sanity checks
"""

# ╔═╡ b289a718-e6ae-40c6-825d-f9706b87d9f1
let
	fig = Figure()
	ax = Axis(fig[1,1])

	v = LinRange(20, 50, 1000)  ./ V2KMS
	x = log10.(v)
	y  = log10.(LilGuys.M_s_from_vel_fattahi.(v))


	lines!(x, y)

	y_approx = log10.(LilGuys.M_s_from_vel_fattahi(v_0)) .+ m_fattahi *(x .- log10(v_0)) 

	lines!(x, y_approx)
	fig
end

# ╔═╡ e7a2cd33-e964-4861-97a3-e1e6ae1b80e0
let
	fig = Figure()
	ax = Axis(fig[1,1])

	v = LinRange(20, 50, 1000)  ./ V2KMS
	x = log10.(v)
	y  = log10.(LilGuys.Ludlow.solve_rmax.(v))


	lines!(x, y)
	x0 = log_v_0
	y0 = log10.(LilGuys.Ludlow.solve_rmax(10^x0))
	y_approx = y0 .+ 1/m_ludlow * (x .- x0) 

	lines!(x, y_approx)
	scatter!(x0, y0)
	fig
end

# ╔═╡ 9fb0dfe4-5d49-4b05-899e-cb4138dc7aac
md"""
These values are from the MCMC halo mass estimator notebook in `analysis/sculptor` 
"""

# ╔═╡ f178ab03-06ba-4e83-9be9-3924543bd288
Σ_scl = [0.019467512114446032 0.0019034735328083907; 0.0019034735328083907 0.001516776581639168]

# ╔═╡ b77ca8f4-f300-412a-8250-2fbe8846b93c
Σ ./ Σ_scl

# ╔═╡ 28f7a671-8b64-43ed-b079-95cd7360a73d
μ_scl = 
[0.780442
1.49899]

# ╔═╡ b3994aa4-8e12-4e4b-8d7a-4646ba3785f7
μ - μ_scl

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
# ╠═f6178d1d-3687-4538-80e7-4f2efe93ab54
# ╠═b4bd8f82-ec5d-490b-8be2-08a70f311e24
# ╠═407d62b9-8e5c-4052-b281-e74b0a3290c4
# ╠═c8454dc7-8aaf-40d8-9231-2e6e4f5f2186
# ╠═d1636c10-37d3-4c92-ad9c-e58a19016469
# ╠═b575403d-b2d3-4c81-9957-8f384c9811ac
# ╠═25e17fae-4ae6-4a3d-99c9-646e3dfb574d
# ╠═0020549b-40f7-4542-8cd8-aad0e047a689
# ╠═f61ca53d-ef76-47ce-bd2f-015645048d44
# ╠═816c93e1-3f6a-4b76-9003-2622e8d1f9fc
# ╠═17dd04ed-2c60-488e-9b73-a8adea7dc9ff
# ╠═73ab0963-1d9a-4488-be1c-5c061b2e0e7e
# ╠═a8e71e9a-cb26-4058-948b-0d27b60c5b48
# ╠═1e7d486e-7d0f-4f51-9175-a1a378774f1e
# ╠═e29ca4df-c319-4e8d-9dfa-93ac7bf28ee8
# ╠═4f4dd26d-3b93-49fe-9d9a-472f504e365d
# ╠═23959023-9412-4899-a546-60fecdf6fe6d
# ╠═2967e69c-8732-49eb-bdb2-91161993a5bb
# ╠═1f74ca1b-ac53-425c-8f89-1afe361525cd
# ╠═b013d48a-7d33-4630-a694-60c7562f3d52
# ╠═a0e7532b-2aaf-4fb3-af00-1d0912787c35
# ╠═b01d7319-2f22-49a5-9303-2ef28b24b798
# ╠═104e4af5-cafc-4de0-82b4-e3726bc7fc6e
# ╠═3909983c-1a96-4fd5-9705-95cb4e9a184c
# ╠═a15dbc88-8dbf-48d4-8a72-1881151c8e26
# ╠═5d0a1520-e8e0-4b90-b83b-f8605ed0cbb5
# ╠═679a471a-da30-4ef4-b8f7-60a69c8d8d94
# ╠═0ca4187b-d63c-48e4-947c-c9a48973323d
# ╠═e0e418d2-16f2-4bbe-92db-ccdc6e89669e
# ╠═a085b24c-9ac6-41c8-811e-8e972eef0efb
# ╠═51075ee8-2cd7-4643-bf45-24a8a60d4eec
# ╠═53967222-f243-4529-9544-caf2d68215c2
# ╠═b289a718-e6ae-40c6-825d-f9706b87d9f1
# ╠═e7a2cd33-e964-4861-97a3-e1e6ae1b80e0
# ╠═9fb0dfe4-5d49-4b05-899e-cb4138dc7aac
# ╠═f178ab03-06ba-4e83-9be9-3924543bd288
# ╠═b77ca8f4-f300-412a-8250-2fbe8846b93c
# ╠═28f7a671-8b64-43ed-b079-95cd7360a73d
# ╠═b3994aa4-8e12-4e4b-8d7a-4646ba3785f7
