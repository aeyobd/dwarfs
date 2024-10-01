### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ d324513e-4268-11ef-1cdd-c712289bd296
begin
	import Pkg; Pkg.activate()
	using CairoMakie

	using LilGuys
	using Arya
end

# ╔═╡ 66870b4e-33b6-4c0d-9585-21bd51b3f3f8
using QuadGK

# ╔═╡ 39f7d95c-1785-4708-a930-c92b8d93c174
md"""
# Explore potentials


This notebook is simply to make some little comparison plots for two different halos to explore the effects on the density and velocity profiles.
"""

# ╔═╡ 57474fd0-6ee8-4a76-b018-9df5ad087367
V0 = 31 / V2KMS

# ╔═╡ d5703c11-08c7-4cd5-9c5e-54c93416da44
LilGuys.Ludlow.solve_rmax(V0)

# ╔═╡ 76e8a69e-6d6a-4a67-b320-e77980428d70
R0 = 5.9

# ╔═╡ 45cc3ae0-0c2f-415e-ade8-9415ae4795d8
V = 31/ V2KMS

# ╔═╡ 53af1c45-d189-4b87-a28c-f1d18b3f3746
n_sigma_R = -2

# ╔═╡ 4e846290-ffeb-4b5d-b6f4-ac8831c1a8be
R = LilGuys.Ludlow.solve_rmax(V) * 10 ^ (0.135 * n_sigma_R)

# ╔═╡ 03769904-cc37-4847-8cad-fb96f9f611dd
LilGuys.Ludlow.solve_rmax(V, -0.1 * n_sigma_R)

# ╔═╡ 40a60cce-0ec8-405f-b8df-78e8b3ee0f8b
LilGuys.Ludlow.solve_rmax(V)

# ╔═╡ 1c7e2102-6905-42c3-9f3a-237e99c3c1a2
0.44 * R / R0

# ╔═╡ 1e641adf-1147-4606-bfb4-693a81892031
p0 = NFW(; v_circ_max=V0, r_circ_max=R0)

# ╔═╡ c9a564e4-35ad-413e-b336-1eb7936d7cd3
p1 = NFW(; v_circ_max=V, r_circ_max=R)

# ╔═╡ c50cde53-613c-427d-8d98-142f018f9559
R_s = 0.13

# ╔═╡ 4f226d3d-9e73-43b2-95b2-f9a3f6ac5f34
stellar_profile = LilGuys.Exp2D(R_s=R_s)

# ╔═╡ 428004de-036a-4d44-8178-255f777421fc
function calc_σv_star_mean(p; stellar_profile=stellar_profile, R_max=Inf, R = 10 .^ LinRange(-4, 5, 1000))
	integrand(r) = calc_ρ(stellar_profile, r) * calc_M(p, r) * LilGuys.G / r^2
	
	weighted_sigma(r) = quadgk(integrand, r, R_max)[1] * 4π * r^2
	mass(r) = calc_ρ(stellar_profile, r) * 4π * r^2
	sigmas = weighted_sigma.(R)
	
	sqrt(sum(sigmas) / sum(mass.(R)))
end

# ╔═╡ 75116f74-8d56-4282-b0cd-49784ac14285
calc_σv_star_mean(p0) * V2KMS

# ╔═╡ 6def8955-87f8-487e-b91e-29b225909f14
calc_σv_star_mean(p1) * V2KMS

# ╔═╡ 9c9253d8-cfbe-4229-ba7c-fd3696fe9680
function calc_σv_star(p, r; stellar_profile=stellar_profile, R_max=Inf)
	integrand(r) = calc_ρ(stellar_profile, r) * calc_M(p, r) * LilGuys.G / r^2

	i, δi = quadgk(integrand, r, R_max)
	if maximum(δi ./ i) > 1e-5
		@warn "relative integration error is large $(δi / i)"
	end
	
	ρ0 = calc_ρ(stellar_profile, r)
	sqrt(i ./ ρ0)
end

# ╔═╡ f0058eab-cd2e-4f27-b1a6-e8e74644d575
begin 
	println(calc_v_circ(p0, 0.2) * V2KMS / √3 / 1.2)
	println(calc_v_circ(p1, 0.15) * V2KMS / √3 / 1.2)
end

# ╔═╡ 4c6d27c3-cc04-4512-bfb6-8f05d3b5ba40
begin 
	#battaglia suggest M(r < 1.6) ≈ 0.034 ± 0.007
	println(calc_M(p0, 1.6))
	println(calc_M(p1, 1.6))
end

# ╔═╡ 13a22e15-2ef3-47a4-8b9d-f93228aa09db
r_model = 10 .^ LinRange(-3, 2, 1000)

# ╔═╡ c3e1e662-b169-4ea9-98d6-eab65b164be6
let
	fig, ax = FigAxis(
		xscale=log10,
		yscale=log10,
		xlabel="r",
		ylabel="rho"
	)

	lines!(r_model, calc_ρ.(p0, r_model))
	lines!(r_model, calc_ρ.(p1, r_model))
	#lines!(r_model, calc_ρ.(p1, r_model) ./ V^2)

	fig
end

# ╔═╡ cdcdb9d4-ea44-4c2b-af34-c9a548a35e70
let
	fig, ax = FigAxis(
		xscale=log10,
		#yscale=log10,
		xlabel="r",
		ylabel="v circ"
	)

	lines!(r_model, LilGuys.calc_v_circ.(p0, r_model) *V2KMS)
	lines!(r_model, LilGuys.calc_v_circ.(p1, r_model) * V2KMS)
	vlines!(0.1)

	fig
end

# ╔═╡ e3bf65f5-a59f-4696-940a-fde472cb080a
0.33 * V2KMS

# ╔═╡ e4328247-467a-417a-9cd9-f8b6d1bb9985
let
	fig, ax = FigAxis(
		xscale=log10,
		yscale=log10,
		xlabel="r",
		ylabel="M in"
	)

	lines!(r_model, LilGuys.calc_M.(p0, r_model))
	lines!(r_model, LilGuys.calc_M.(p1, r_model))

	fig
end

# ╔═╡ 9f0ebcb4-de65-47a6-be00-f5bae00f72ac
let
	fig, ax = FigAxis(
		xscale=log10,
		yscale=log10,
		xlabel="r",
		ylabel="tc"
	)

	lines!(r_model, r_model ./ LilGuys.calc_v_circ.(p0, r_model) * T2GYR)
	lines!(r_model, r_model ./ LilGuys.calc_v_circ.(p1, r_model) * T2GYR)

	fig
end

# ╔═╡ Cell order:
# ╟─39f7d95c-1785-4708-a930-c92b8d93c174
# ╠═d324513e-4268-11ef-1cdd-c712289bd296
# ╠═57474fd0-6ee8-4a76-b018-9df5ad087367
# ╠═d5703c11-08c7-4cd5-9c5e-54c93416da44
# ╠═76e8a69e-6d6a-4a67-b320-e77980428d70
# ╠═45cc3ae0-0c2f-415e-ade8-9415ae4795d8
# ╠═4e846290-ffeb-4b5d-b6f4-ac8831c1a8be
# ╠═53af1c45-d189-4b87-a28c-f1d18b3f3746
# ╠═03769904-cc37-4847-8cad-fb96f9f611dd
# ╠═40a60cce-0ec8-405f-b8df-78e8b3ee0f8b
# ╠═1c7e2102-6905-42c3-9f3a-237e99c3c1a2
# ╠═1e641adf-1147-4606-bfb4-693a81892031
# ╠═c9a564e4-35ad-413e-b336-1eb7936d7cd3
# ╠═4f226d3d-9e73-43b2-95b2-f9a3f6ac5f34
# ╠═c50cde53-613c-427d-8d98-142f018f9559
# ╠═66870b4e-33b6-4c0d-9585-21bd51b3f3f8
# ╠═75116f74-8d56-4282-b0cd-49784ac14285
# ╠═6def8955-87f8-487e-b91e-29b225909f14
# ╠═428004de-036a-4d44-8178-255f777421fc
# ╠═9c9253d8-cfbe-4229-ba7c-fd3696fe9680
# ╠═f0058eab-cd2e-4f27-b1a6-e8e74644d575
# ╠═4c6d27c3-cc04-4512-bfb6-8f05d3b5ba40
# ╠═13a22e15-2ef3-47a4-8b9d-f93228aa09db
# ╠═c3e1e662-b169-4ea9-98d6-eab65b164be6
# ╠═cdcdb9d4-ea44-4c2b-af34-c9a548a35e70
# ╠═e3bf65f5-a59f-4696-940a-fde472cb080a
# ╠═e4328247-467a-417a-9cd9-f8b6d1bb9985
# ╠═9f0ebcb4-de65-47a6-be00-f5bae00f72ac
