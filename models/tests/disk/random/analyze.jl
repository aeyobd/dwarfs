### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 71e3f40f-a912-4bcb-aa30-313f8b5dae9e
begin 
	import Pkg
	Pkg.activate()
	import LilGuys as lguys
	using Plots

	pwd()
end

# ╔═╡ 7fde36fd-8f27-45c3-b255-989e394998c2
out = lguys.Output("out")

# ╔═╡ 8b30ba3c-8441-432c-a893-8e38db3599f5
begin 
	id = 4
	positions = lguys.extract(out, :positions, id)
	velocities = lguys.extract(out, :velocities, id)
	Φs = lguys.extract(out, :Φs_ext, id)
end

# ╔═╡ 05b392b8-43eb-4f97-81d5-985513497908
lguys.plot_xyz(positions)

# ╔═╡ eb154a33-67da-4bb4-819d-1a4fefdbc66d
begin 
	rs = lguys.calc_r(positions)
	Φ_nfw(r, M, r_s, c) = -lguys.G* M * 1/r * log(1 + r/r_s) / (log(1+c) - c/(1+c))
	vs = lguys.calc_r(velocities)
end

# ╔═╡ f9dd10c3-40cd-492c-9d62-0fa505ee065b
plot(Φ_nfw.(rs, 10, 10, 10) ./ Φs)

# ╔═╡ e798f63f-c5ad-46c7-a07e-478f95635367
begin
	plot(out.times, lguys.calc_E_spec.(Φs, vs))
	ylabel!("total energy")
end

# ╔═╡ dd0f4746-5873-433c-918c-97ad34f9c770
begin 
	Ls = lguys.calc_L_spec(positions, velocities)

	plot(transpose(Ls .- Ls[:, 1]))
end

# ╔═╡ 45c20cb4-fbc6-4e5b-987f-080b7a66fd9e
md"""
Check that the peri and apocenter are as expected...
"""

# ╔═╡ 3a106dcd-7e52-4f58-a9ab-bb88ba7137ec
begin
	calc_Φ(r) = Φ_nfw(r, 10, 10, 10)
	L = lguys.calc_r(Ls[:, 1])
	E = 1/2 * vs[1] + Φs[1]

end

# ╔═╡ 49ed11e0-9cfb-4309-b64b-ad7b4118280d
begin 
	peri_obs = minimum(rs)
	apo_obs = maximum(rs)
	peri_exp = find_zero(f, peri_obs)
	apo_exp = find_zero(f, apo_obs)

	println(peri_obs/peri_exp)
	println(apo_obs/apo_exp)

	T_exp = quadgk(g, peri_exp + 1e-4, apo_exp - 1e-4)
end

# ╔═╡ 367a74bc-806f-4e3a-99fc-833bb85a9e6a
begin 
	xs_1 = LinRange(0.7*peri_exp, 1.3*apo_exp, 1000)
	plot(xs_1, f.(xs_1))
	vline!([apo_exp, peri_exp])
	hline!([0])
end

# ╔═╡ ab00e3c3-fbeb-4e9b-882f-f5c425e99e88
begin
	plot(out.times, rs)
	hline!([apo_exp, peri_exp])
	xlabel!("time")
	ylabel!("r")
end

# ╔═╡ f32c522f-a979-40c3-a64f-97a7053b95b9
begin
	plot(out.times, 1/2*vs.^2, label="T")
	plot!(out.times, Φs, label="V")
	plot!(out.times, Φs + 1/2*vs.^2, label="total")
	xlabel!("time")
	ylabel!("energy")
end

# ╔═╡ bfab46bf-970f-4027-8a3d-7f69a0f736af
plot(out.times, positions[3, :])

# ╔═╡ Cell order:
# ╠═71e3f40f-a912-4bcb-aa30-313f8b5dae9e
# ╠═7fde36fd-8f27-45c3-b255-989e394998c2
# ╠═8b30ba3c-8441-432c-a893-8e38db3599f5
# ╠═05b392b8-43eb-4f97-81d5-985513497908
# ╠═eb154a33-67da-4bb4-819d-1a4fefdbc66d
# ╠═f9dd10c3-40cd-492c-9d62-0fa505ee065b
# ╠═e798f63f-c5ad-46c7-a07e-478f95635367
# ╠═dd0f4746-5873-433c-918c-97ad34f9c770
# ╟─45c20cb4-fbc6-4e5b-987f-080b7a66fd9e
# ╠═3a106dcd-7e52-4f58-a9ab-bb88ba7137ec
# ╠═367a74bc-806f-4e3a-99fc-833bb85a9e6a
# ╠═49ed11e0-9cfb-4309-b64b-ad7b4118280d
# ╠═ab00e3c3-fbeb-4e9b-882f-f5c425e99e88
# ╠═f32c522f-a979-40c3-a64f-97a7053b95b9
# ╠═bfab46bf-970f-4027-8a3d-7f69a0f736af
