### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 71e3f40f-a912-4bcb-aa30-313f8b5dae9e
begin 
	import Pkg
	Pkg.activate()
	import LilGuys as lguys
	using Plots; plotly()
	using QuadGK
	using Roots

	pwd()
end

# ╔═╡ 7fde36fd-8f27-45c3-b255-989e394998c2
begin 
	out = lguys.Output("out")
	M = 115.0
	c = 9.545
	r_s = 20.2
	calc_Φ(r) = -lguys.G* M * 1/r * log(1 + r/r_s) / (log(1+c) - c/(1+c))
end

# ╔═╡ 8b30ba3c-8441-432c-a893-8e38db3599f5
begin 
	id = 45
	positions = lguys.extract(out, :positions, id)
	velocities = lguys.extract(out, :velocities, id)
	accelerations = lguys.extract(out, :accelerations, id)

	Φs = lguys.extract(out, :Φs_ext, id)
end

# ╔═╡ 05b392b8-43eb-4f97-81d5-985513497908
lguys.plot_xyz(positions)

# ╔═╡ 85cc1e79-2556-46cf-823f-485ca908b0b5
lguys.plot_xyz(velocities)

# ╔═╡ 6a376f08-b2f7-4fc4-bf95-bce4c978afbe
lguys.plot_xyz(accelerations)

# ╔═╡ eb154a33-67da-4bb4-819d-1a4fefdbc66d
begin 
	rs = lguys.calc_r(positions)
	vs = lguys.calc_r(velocities)
end

# ╔═╡ f9dd10c3-40cd-492c-9d62-0fa505ee065b
plot(calc_Φ.(rs) ./ Φs)

# ╔═╡ e798f63f-c5ad-46c7-a07e-478f95635367
begin
	plot(out.times, lguys.calc_E_spec.(Φs, vs))
	ylabel!("total energy")
end

# ╔═╡ dd0f4746-5873-433c-918c-97ad34f9c770
begin 
	Ls = lguys.calc_L_spec(positions, velocities)

	L_rel_err = (Ls .- Ls[:, 1]) ./ Ls[:, 1]
	plot(transpose(L_rel_err))
end

# ╔═╡ 45c20cb4-fbc6-4e5b-987f-080b7a66fd9e
md"""
Check that the peri and apocenter are as expected...
"""

# ╔═╡ 3a106dcd-7e52-4f58-a9ab-bb88ba7137ec
begin
	L = lguys.calc_r(Ls[:, 1])
	E = 1/2 * vs[1]^2 + Φs[1]
	
	f(r) =  r^-2 + 2*(calc_Φ(r) - E) / L^2
	g(r) = 2 / sqrt(2*(E - calc_Φ(r)) - L^2/r^2)
end

# ╔═╡ 49ed11e0-9cfb-4309-b64b-ad7b4118280d
begin 
	peri_obs = minimum(rs)
	apo_obs = maximum(rs)
	peri_exp = find_zero(f, peri_obs)
	apo_exp = find_zero(f, apo_obs)

	println("rel peri error ", peri_obs/peri_exp)
	println("rel apo error ", apo_obs/apo_exp)

	ϵ= 1e-8
	T_exp, int_err = quadgk(g, peri_exp + ϵ, apo_exp - ϵ)
	println("integration relerr ", int_err)
end

# ╔═╡ 367a74bc-806f-4e3a-99fc-833bb85a9e6a
begin 
	xs_1 = LinRange(0.7*peri_exp, 1.3*apo_exp, 1000)
	plot(xs_1, f.(xs_1))
	vline!([apo_exp, peri_exp])
	hline!([0])
end

# ╔═╡ 4e35ac67-fe15-44f0-8a67-3e552660d672
begin
	xs_2 = LinRange(peri_exp+ϵ, apo_exp-ϵ, 1000)
	plot(xs_2, g.(xs_2))
	ylims!(0, 20)
end

# ╔═╡ ab00e3c3-fbeb-4e9b-882f-f5c425e99e88
begin
	t_0 = out.times[argmax(rs[out.times .< T_exp])]
	N_periods = floor(Int, (out.times[end] - t_0) / T_exp)
	
	plot(out.times, rs)
	hline!([apo_exp, peri_exp])
	vline!(collect(0:N_periods) .* T_exp .+ t_0)
	xlabel!("time")
	ylabel!("r")
end

# ╔═╡ 3506bc94-695e-49d5-aaca-7ba2448565a8
collect(1:3) * T_exp[1]

# ╔═╡ f32c522f-a979-40c3-a64f-97a7053b95b9
begin
	plot(out.times, 1/2*vs.^2, label="T")
	plot!(out.times, Φs, label="V")
	plot!(out.times, Φs + 1/2*vs.^2, label="total")
	xlabel!("time")
	ylabel!("energy")
end

# ╔═╡ 972cb86d-b688-4a66-ad7b-babcfd0f7291
scatter(out.times, lguys.calc_r(accelerations))

# ╔═╡ Cell order:
# ╠═71e3f40f-a912-4bcb-aa30-313f8b5dae9e
# ╠═7fde36fd-8f27-45c3-b255-989e394998c2
# ╠═8b30ba3c-8441-432c-a893-8e38db3599f5
# ╠═05b392b8-43eb-4f97-81d5-985513497908
# ╠═85cc1e79-2556-46cf-823f-485ca908b0b5
# ╠═6a376f08-b2f7-4fc4-bf95-bce4c978afbe
# ╠═eb154a33-67da-4bb4-819d-1a4fefdbc66d
# ╠═f9dd10c3-40cd-492c-9d62-0fa505ee065b
# ╠═e798f63f-c5ad-46c7-a07e-478f95635367
# ╠═dd0f4746-5873-433c-918c-97ad34f9c770
# ╟─45c20cb4-fbc6-4e5b-987f-080b7a66fd9e
# ╠═3a106dcd-7e52-4f58-a9ab-bb88ba7137ec
# ╠═367a74bc-806f-4e3a-99fc-833bb85a9e6a
# ╠═4e35ac67-fe15-44f0-8a67-3e552660d672
# ╠═49ed11e0-9cfb-4309-b64b-ad7b4118280d
# ╠═ab00e3c3-fbeb-4e9b-882f-f5c425e99e88
# ╠═3506bc94-695e-49d5-aaca-7ba2448565a8
# ╠═f32c522f-a979-40c3-a64f-97a7053b95b9
# ╠═972cb86d-b688-4a66-ad7b-babcfd0f7291
