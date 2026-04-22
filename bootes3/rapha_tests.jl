### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 156c47a4-39b1-11f1-861a-df30869de806
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using CairoMakie, Arya
end

# ╔═╡ 35e39c2d-cac3-4686-8ac4-37027f5b56a2
vmax_i = 31 / V2KMS

# ╔═╡ 95f28817-ec88-4e40-bf1d-f499f0217451
rmax_i = 3.2

# ╔═╡ 64233a5a-cf7a-4770-b154-b6e2fa759cbe
peri = 20

# ╔═╡ ccfa0f7b-d942-4c31-b657-9e2214e2db0b
apo = 100

# ╔═╡ d6b145b6-4be1-488e-9572-aff43e17e6a0
V0 = 220/V2KMS

# ╔═╡ 9f80a0a7-d83b-4a27-8038-dc375a82dc15
CairoMakie.activate!(type=:png)

# ╔═╡ fddd99d4-fd2f-4f90-bf04-8f8d20c39c48
dwarfs_dir = ENV["DWARFS_ROOT"]

# ╔═╡ cc264eb9-65df-4fc0-9be1-26142ecd2690
module Rapha
	include(joinpath(ENV["DWARFS_ROOT"], "utils/rapha_utils.jl"))
end

# ╔═╡ 67580761-e143-4824-b829-64d34496675a
import Agama

# ╔═╡ a829e274-961f-456d-89e6-fc5f1014af1e
pot = Agama.Potential(type="Logarithmic", v0=220/V2KMS, scaleRadius=0)

# ╔═╡ 8f7c919f-b7d4-4a07-9061-cd2247af4962
md"""
The below cell sanity checks that the potential is isothermal (`v_circ` is constant)
"""

# ╔═╡ 4f64e911-e57c-4756-bdda-834cb4948d93
@assert all(Agama.circular_velocity(pot, [1e-10, 0.1, 1, 100, 10000]) ./ V0 .≈ 1)

# ╔═╡ 54f2797a-ed34-41d2-b1c5-e23b7eb48f6f
md"""
# Tidal tracks
"""

# ╔═╡ 929b599c-38a3-48ad-9a24-41e85e35ff3e
md"""
# Tidal prediction
"""

# ╔═╡ 6b617252-7929-4421-bcad-d05c64bfc388
T_peri = 2π * peri / V0

# ╔═╡ 39aad1b7-080c-4d3e-89ad-12dfdd2ba27e
md"""
Reproduces Fig. 5
"""

# ╔═╡ 23ef0e79-15cd-4c23-80c5-8a10cd2341a2
Φ(r) = Agama.potential(pot, [0, 0, r])

# ╔═╡ 7d071e1f-9d49-48a2-b568-5772b662b8ad
function peri_apo_to_E_L(peri, apo)
	L = sqrt(
		(Φ(peri) - Φ(apo)) / 
			(1/2 * (1/apo^2 - 1/peri^2))
	)

	E = Φ(peri) + L^2 / 2peri^2

	@assert Φ(apo) + L^2 / 2apo^2 ≈ E
	return E, L
end

# ╔═╡ a9bbb390-fbc9-475b-bbc1-1d557ad9c232
function orbital_period(peri, apo)
	E, L = peri_apo_to_E_L(peri, apo)
	return 2*LilGuys.integrate(
		r -> 1 / sqrt(2*(E - Φ(r)) - L^2 / r^2),
		peri * (1 + 1e-12), apo*(1 - 1e-12)
	)
end

# ╔═╡ 6702b75a-dfba-4894-8729-10af81b0034c
function make_track(; rmax_i=rmax_i, vmax_i=vmax_i, peri=peri, apo=apo, t_rel_max=20, V0=V0)

	f_ecc = Rapha.ecc_factor(peri, apo)

	t_rel = LinRange(0, 20, 1000) |> collect

	
	rvmax = Rapha.rapha_final_halo.(rmax_i, vmax_i, peri, apo, t_rel .* f_ecc, V0=V0)


	T_peri = 2π * peri / V0

	
	rmax = first.(rvmax)
	vmax = last.(rvmax)

	t_max = @. rmax * 2π / vmax

	LilGuys.DataFrame(
		r_max = rmax,
		v_max = vmax,
		t_max = t_max,
		t_rel = t_rel,
		time = orbital_period(peri, apo) * t_rel .* f_ecc,
		M_max = vmax .^2 .* rmax, 
	)
end

# ╔═╡ 4db3abf0-aa4a-4928-9806-687ee1a09a54
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 limits=(0, 20, -0.7, -0.1))

	vmax_i = 0.255

	scatter!(0, log10(2π * rmax_i /  vmax_i/ T_peri))

	for apo_rel in [1, 5, 10, 20]
		df = make_track(apo=apo*apo_rel, vmax_i=vmax_i)

		filt = df.time .< 200 / T2GYR

		lines!(df.t_rel[filt], log10.(df.t_max[filt] ./ T_peri), linewidth = 1 + log10(apo_rel))
	end

	fig
end

# ╔═╡ c76e56e8-1425-4aaf-bb10-fafd70a148b4
T_orb = orbital_period(peri, apo)

# ╔═╡ bef630c1-444e-4d8a-b7bd-5ec3cbe3da59
2π * peri / Agama.circular_velocity(pot, apo)

# ╔═╡ 1733ed18-3466-458f-bfd2-3cf34ec69540
E, L = peri_apo_to_E_L(peri, apo)

# ╔═╡ f397349c-b23e-42d0-a7d1-933f25f3e4c6
orbit = LilGuys.agama_orbit(pot, Galactocentric([0, 0, apo], [0, L/apo * V2KMS, 0]), timerange=(0, 10*T_orb))

# ╔═╡ 74f23d66-3e3b-47aa-ba75-83b6d46b8595
orbit.pericenter, peri

# ╔═╡ d63ad974-c5a2-44b5-b24c-7d7a8b5d46f4
orbit.apocenter, apo

# ╔═╡ 520ae3b5-5a01-425e-887b-3b210954b123
let
	fig = Figure()
	ax = Axis(fig[1,1])

	lines!(LilGuys.times(orbit) / T_orb, radii(orbit))

	fig
end

# ╔═╡ 99df26cf-1707-4c9a-815f-d367652fa746
E0 = Φ(apo) + 1/2 * (L/apo)^2

# ╔═╡ fa14719c-4462-4385-a423-ab63f6e3e0b5
Rapha.ecc_factor(1, 1)

# ╔═╡ fae22d07-0f9f-4e04-aa22-8878e92e2855
Rapha.ecc_factor(1, 5)

# ╔═╡ 1781ca2e-b151-453c-be5d-c1e964c6411d
Rapha.ecc_factor(1, 10)

# ╔═╡ bb96ffc3-6caf-423d-93d0-26f1dd228ef3
Rapha.ecc_factor(1, 20)

# ╔═╡ f2c04684-d489-46ac-a2e7-deafcffecfd5
md"""
# Reproduces Fig. 4
"""

# ╔═╡ 7d9b6f14-c2e3-4d73-869d-43f94f6cfe94
let
	fig = Figure()
	ax = Axis(fig[1,1])

	for log_tmax in LinRange(0.2, 0.6, 10)
		df = make_track(vmax_i = 2π*rmax_i / (1/T2GYR * 10 .^ log_tmax ), peri=80, apo=1.05*80)

		filt = df.M_max .> 0.003 * df.M_max[1]

		lines!(df.t_rel[filt], log10.(df.t_max[filt] * T2GYR ))
	end

	# hlines!(log10(T_peri / 4 * T2GYR))
	fig
end

# ╔═╡ 1ffb8d96-6b61-4cba-9c73-c6188bdcfafe
LilGuys.v_circ_max(NFW(M200=3.5, c=11))

# ╔═╡ f49840ed-dc5f-4900-a9c4-f6bab47fa721
0.2769600010857212^2 * 13.604511226243485

# ╔═╡ 97fd14cd-d84d-4390-9a76-531ee5b8cf13
LilGuys.find_zero(v ->  v^2 * LilGuys.Ludlow.solve_rmax(v) - 1, 0.3)

# ╔═╡ 70e9185f-a6f5-4c20-bf78-f28b7d68c832
LilGuys.Ludlow.solve_rmax(0.27803169630265384)

# ╔═╡ f73db2cd-6c3c-462a-95bb-a590842e97dc
md"""
Reproduces Fig. 13 of EN2021
"""

# ╔═╡ 4dfa1123-6b47-47b7-be1e-8a8dee0e4779
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 aspect = 1, 
			 limits=(-2, 1.5, 0, 2.5))

	rmax_i=2.209
	vmax_i = 0.06620374917628208
	
	track = make_track(rmax_i=rmax_i, vmax_i=vmax_i,  peri=20, apo=20.01)

	lines!(log10.(track.r_max), log10.(track.v_max * V2KMS),)

	for i in 0:20
		r, v = Rapha.rapha_final_halo.(rmax_i, vmax_i, 20, 20.01, i, V0=V0)
		scatter!(log10(r), log10(v*V2KMS), color=COLORS[1], marker=:circle, strokewidth=0)
	end


	rmax_i=12.936338802835191
	vmax_i = 0.27803169630265384
	
	track = make_track(rmax_i=rmax_i, vmax_i=vmax_i,  peri=20, apo=20.01)

	lines!(log10.(track.r_max), log10.(track.v_max * V2KMS),)

	for i in 0:20
		r, v = Rapha.rapha_final_halo.(rmax_i, vmax_i, 20, 20.01, i, V0=V0)
		scatter!(log10(r), log10(v*V2KMS), color=COLORS[2], marker=:circle, strokewidth=0)
	end


	
	fig
end



# ╔═╡ a2865557-31d4-4bf0-9149-87b9a7b85fbd
md"""
# Extra
"""

# ╔═╡ 32b0e8fb-5b79-4e3f-97b2-837252ffa7dd


# ╔═╡ b3724f59-d9ab-45d9-a0fb-de0217a4671c
md"""
## Does the `r_cut` work with the tidal tracks?
"""

# ╔═╡ 01296f19-9c0b-41ff-8961-ea53f14059ea
LilGuys.r_circ_max(Rapha.TidalNFW(r_circ_max=1, v_circ_max=1, r_cut=0.3))

# ╔═╡ 4c331187-8192-47ec-aea3-1060ec1f33ea
LilGuys.v_circ_max(Rapha.TidalNFW(r_circ_max=1, v_circ_max=1, r_cut=0.3))

# ╔═╡ 3fc349a5-5c33-4489-9c9b-54fd01f8b801
let
	fig = Figure()
	ax = Axis(fig[1,1])

	x = LinRange(-2, 2, 1000)

	nfw_i = NFW(r_circ_max=1, v_circ_max=1)
	y = log10.(LilGuys.v_circ.(nfw_i, 10 .^ x))
	lines!(x, y, color=:black, linewidth=2)

	for M_rel in logrange(1, 1e-2, 30)
		h = Rapha.final_halo(1, 1, M_rel)
	
		y = log10.(LilGuys.v_circ.(h, 10 .^ x))
		lines!(x, y)
		# r_max = LilGuys.solve_r_circ_max(h)
		# scatter!(r_max, LilGuys.v_circ(h, r_max))
		i = argmax(y)
		scatter!(x[i], y[i])
		r_rel, v_rel = Rapha.r_v_rel(M_rel)
		@info r_rel, v_rel
		scatter!(log10(r_rel), log10(v_rel), color=:black, marker=:circle, strokewidth=0)
	end

	

	fig
end

# ╔═╡ Cell order:
# ╠═35e39c2d-cac3-4686-8ac4-37027f5b56a2
# ╠═95f28817-ec88-4e40-bf1d-f499f0217451
# ╠═64233a5a-cf7a-4770-b154-b6e2fa759cbe
# ╠═ccfa0f7b-d942-4c31-b657-9e2214e2db0b
# ╠═d6b145b6-4be1-488e-9572-aff43e17e6a0
# ╠═156c47a4-39b1-11f1-861a-df30869de806
# ╠═9f80a0a7-d83b-4a27-8038-dc375a82dc15
# ╠═fddd99d4-fd2f-4f90-bf04-8f8d20c39c48
# ╠═cc264eb9-65df-4fc0-9be1-26142ecd2690
# ╠═67580761-e143-4824-b829-64d34496675a
# ╠═a829e274-961f-456d-89e6-fc5f1014af1e
# ╠═8f7c919f-b7d4-4a07-9061-cd2247af4962
# ╠═4f64e911-e57c-4756-bdda-834cb4948d93
# ╟─54f2797a-ed34-41d2-b1c5-e23b7eb48f6f
# ╟─929b599c-38a3-48ad-9a24-41e85e35ff3e
# ╠═6702b75a-dfba-4894-8729-10af81b0034c
# ╠═6b617252-7929-4421-bcad-d05c64bfc388
# ╟─39aad1b7-080c-4d3e-89ad-12dfdd2ba27e
# ╠═4db3abf0-aa4a-4928-9806-687ee1a09a54
# ╠═23ef0e79-15cd-4c23-80c5-8a10cd2341a2
# ╠═a9bbb390-fbc9-475b-bbc1-1d557ad9c232
# ╠═7d071e1f-9d49-48a2-b568-5772b662b8ad
# ╠═c76e56e8-1425-4aaf-bb10-fafd70a148b4
# ╠═f397349c-b23e-42d0-a7d1-933f25f3e4c6
# ╠═74f23d66-3e3b-47aa-ba75-83b6d46b8595
# ╠═d63ad974-c5a2-44b5-b24c-7d7a8b5d46f4
# ╠═99df26cf-1707-4c9a-815f-d367652fa746
# ╠═520ae3b5-5a01-425e-887b-3b210954b123
# ╠═bef630c1-444e-4d8a-b7bd-5ec3cbe3da59
# ╠═1733ed18-3466-458f-bfd2-3cf34ec69540
# ╠═fa14719c-4462-4385-a423-ab63f6e3e0b5
# ╠═fae22d07-0f9f-4e04-aa22-8878e92e2855
# ╠═1781ca2e-b151-453c-be5d-c1e964c6411d
# ╠═bb96ffc3-6caf-423d-93d0-26f1dd228ef3
# ╠═f2c04684-d489-46ac-a2e7-deafcffecfd5
# ╠═7d9b6f14-c2e3-4d73-869d-43f94f6cfe94
# ╠═1ffb8d96-6b61-4cba-9c73-c6188bdcfafe
# ╠═f49840ed-dc5f-4900-a9c4-f6bab47fa721
# ╠═97fd14cd-d84d-4390-9a76-531ee5b8cf13
# ╠═70e9185f-a6f5-4c20-bf78-f28b7d68c832
# ╟─f73db2cd-6c3c-462a-95bb-a590842e97dc
# ╠═4dfa1123-6b47-47b7-be1e-8a8dee0e4779
# ╠═a2865557-31d4-4bf0-9149-87b9a7b85fbd
# ╠═32b0e8fb-5b79-4e3f-97b2-837252ffa7dd
# ╠═b3724f59-d9ab-45d9-a0fb-de0217a4671c
# ╠═01296f19-9c0b-41ff-8961-ea53f14059ea
# ╠═4c331187-8192-47ec-aea3-1060ec1f33ea
# ╠═3fc349a5-5c33-4489-9c9b-54fd01f8b801
