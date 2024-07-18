### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ d324513e-4268-11ef-1cdd-c712289bd296
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya
	using CairoMakie
end

# ╔═╡ 57474fd0-6ee8-4a76-b018-9df5ad087367
V0 = 31.5 / V2KMS

# ╔═╡ d5703c11-08c7-4cd5-9c5e-54c93416da44
LilGuys.Ludlow.solve_rmax(V0)

# ╔═╡ 76e8a69e-6d6a-4a67-b320-e77980428d70
R0 = 5.938

# ╔═╡ 45cc3ae0-0c2f-415e-ade8-9415ae4795d8
V = 70/ V2KMS

# ╔═╡ 53af1c45-d189-4b87-a28c-f1d18b3f3746
n_sigma_R = -3

# ╔═╡ 4e846290-ffeb-4b5d-b6f4-ac8831c1a8be
R = LilGuys.Ludlow.solve_rmax(V) * 10 ^ (0.135 * n_sigma_R)

# ╔═╡ 03769904-cc37-4847-8cad-fb96f9f611dd
LilGuys.Ludlow.solve_rmax(V, -0.1 * n_sigma_R)

# ╔═╡ 1c7e2102-6905-42c3-9f3a-237e99c3c1a2
0.44 * R / R0

# ╔═╡ 1e641adf-1147-4606-bfb4-693a81892031
p0 = NFW(; v_circ_max=V0, r_circ_max=R0)

# ╔═╡ c9a564e4-35ad-413e-b336-1eb7936d7cd3
p1 = NFW(; v_circ_max=V, r_circ_max=R)

# ╔═╡ c50cde53-613c-427d-8d98-142f018f9559
r_h = 0.15

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
		yscale=log10,
		xlabel="r",
		ylabel="v circ"
	)

	lines!(r_model, LilGuys.calc_v_circ.(p0, r_model))
	lines!(r_model, LilGuys.calc_v_circ.(p1, r_model))

	fig
end

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
# ╠═d324513e-4268-11ef-1cdd-c712289bd296
# ╠═57474fd0-6ee8-4a76-b018-9df5ad087367
# ╠═d5703c11-08c7-4cd5-9c5e-54c93416da44
# ╠═76e8a69e-6d6a-4a67-b320-e77980428d70
# ╠═45cc3ae0-0c2f-415e-ade8-9415ae4795d8
# ╠═53af1c45-d189-4b87-a28c-f1d18b3f3746
# ╠═4e846290-ffeb-4b5d-b6f4-ac8831c1a8be
# ╠═03769904-cc37-4847-8cad-fb96f9f611dd
# ╠═1c7e2102-6905-42c3-9f3a-237e99c3c1a2
# ╠═1e641adf-1147-4606-bfb4-693a81892031
# ╠═c9a564e4-35ad-413e-b336-1eb7936d7cd3
# ╠═c50cde53-613c-427d-8d98-142f018f9559
# ╠═f0058eab-cd2e-4f27-b1a6-e8e74644d575
# ╠═4c6d27c3-cc04-4512-bfb6-8f05d3b5ba40
# ╠═13a22e15-2ef3-47a4-8b9d-f93228aa09db
# ╠═c3e1e662-b169-4ea9-98d6-eab65b164be6
# ╠═cdcdb9d4-ea44-4c2b-af34-c9a548a35e70
# ╠═e4328247-467a-417a-9cd9-f8b6d1bb9985
# ╠═9f0ebcb4-de65-47a6-be00-f5bae00f72ac
