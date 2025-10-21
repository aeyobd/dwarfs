### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 9991bd2a-77c4-11f0-2bf0-2fbe36425e63
begin 
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya
	using CairoMakie
end

# ╔═╡ 4f471846-ae8b-419e-8cd5-920032f19c82
import LilGuys: R200, M200, find_zero

# ╔═╡ b88a3b23-6f3e-487b-848c-24fab0e8a286
a1, a2, a3 = [2.9, 0.614, 0.995] #[2.521, 0.729, 0.988]

# ╔═╡ c8f78a5b-fbf5-4b0f-849a-eec18e8a4639
function c_mah(c)
	return ((c / a1)^(1/a3)  - 1) / a2
end

# ╔═╡ 5ff465a3-846e-4fac-b2c3-bbe92a5a068e
function c_enclosed_density(c)
	return a1 * (1 + a2 * c)^a3
end

# ╔═╡ cc36e1fc-cd68-45c8-a1dd-4a75e19e3057
ρ_c(z) = LilGuys.ρ_crit * (1+z)^3 

# ╔═╡ d690b0c4-3f9a-4c46-b356-f439982adbba
ρ_0 = ρ_c(0.)

# ╔═╡ 74437045-5ee0-4c80-b8b9-90adf1a06f61
function scale_density(halo)
	return LilGuys.mean_density(halo, halo.r_s)
end

# ╔═╡ f433807e-620a-4191-b94e-ebe9f68d3a21
function scale_point_density(halo)
	return LilGuys.density(halo, halo.r_s)
end

# ╔═╡ 8325d421-a7a7-41ee-a0e3-3ff1b0df0c22
function mass_given_density(halo, ρ)
	log_r = LilGuys.find_zero(lr -> LilGuys.mean_density(halo, 10 ^ lr) - ρ, 
			log10(halo.r_s))
	return LilGuys.mass(halo, 10 ^ log_r)
end

# ╔═╡ cd58027b-ff2d-45d6-9b37-08d95eaf40ab
function MAH(nfw; c_scale=1)
	c = c_mah(nfw.c) * c_scale

	M = find_zero(M -> mass_given_density(NFW(M200=M, c=c), 200ρ_0) - M200(nfw), M200(nfw))
	
	return NFW(M200=M, c=c)
end

# ╔═╡ ac6e60e7-b9b1-4c98-9144-061d25429c40
md"""
## Calculations and checks
"""

# ╔═╡ 4c59a8e8-1b59-4817-802a-79db329b1eb3
A_ρ = 853

# ╔═╡ fbafc7d5-cfde-4a79-ba98-77157d951cc9
ρ_c(0)

# ╔═╡ 72e2f6a1-1214-4cc5-bb7b-ca4148134f16
A_ρ

# ╔═╡ e684cc0e-089d-46c1-b779-9d9437eba071
log10(ρ_c(4) / ρ_c(0))

# ╔═╡ 87182b31-fd19-48f4-abaf-90d2cd5f2bb4
md"""
# Input
"""

# ╔═╡ 423f83f7-3530-4689-9aa9-0fbc50df0e54
nfw = NFW(M200 = 1e5, c=5)

# ╔═╡ 6fa3783a-14b3-4943-a5d8-a04a48effd43
scale_density(nfw), scale_point_density(nfw)

# ╔═╡ b88a0324-6a65-4bd3-8ee3-05e4fbd04a44
scale_density(nfw)

# ╔═╡ 17276a12-4365-487a-bcfe-28c56c4f6cb4
1/A_ρ * scale_density(nfw)

# ╔═╡ de3d6614-0f58-4110-bcc0-ce61e2f57fae
nfw_ah = MAH(nfw; c_scale=1)

# ╔═╡ 6cc3c980-a208-4484-a75f-05efd461c73f
mass_given_density(nfw_ah, ρ_c(51.5))

# ╔═╡ 311ec430-fcf6-4fe1-b167-e548f5061ef9
M200(nfw), mass_given_density(nfw_ah, ρ_c(0.)), M200(nfw_ah)

# ╔═╡ 096ed683-3b2a-4a13-a7ec-f1300904d0e4
mass_given_density(nfw_ah, ρ_c(0.))/ mass_given_density(nfw_ah, ρ_c(2.))

# ╔═╡ 2710b8d2-2339-46bf-9951-8b45b58d8d71
@assert LilGuys.concentration(nfw_ah) ≈ nfw_ah.c

# ╔═╡ fdb40c9f-59a7-4403-94ae-cc7b09a11083
c_mah(nfw.c), c_enclosed_density(nfw.c)

# ╔═╡ 0065b475-ed69-405e-8a83-08026535b9ab
@assert mass_given_density(nfw, 200ρ_c(0)) ≈ LilGuys.M200(nfw)

# ╔═╡ a8a26ba9-d0cc-437d-afaa-66ee16d483b9
@assert mass_given_density(nfw_ah, 200ρ_c(0)) ≈ LilGuys.M200(nfw)

# ╔═╡ 7144c3df-bde4-46a6-ba1c-744894691eb5
@assert nfw_ah.c ≈ c_mah(nfw.c)

# ╔═╡ 2a222c23-5873-46de-b900-2abc8ab4d4a9
LilGuys.Ludlow.c_ludlow(M200(nfw), 0.)

# ╔═╡ c157c90e-aa07-43d0-bfcb-9a2d780dc05b
md"""
# Reproduced plots
"""

# ╔═╡ dde42a8a-e641-4eaf-a513-411e1fcf096e
let
	fig = Figure()
	Ax = Axis(fig[1,1],
			 xlabel = L"\log c_\textrm{MAH}",
			 ylabel = L"\log c_\textrm{NFW}",
			 
		)


	x = LinRange(-0.5, 0.6, 1000)
	y = @. log10(c_enclosed_density(10^x))
	x2 = @.log10(c_mah(10^y))
	lines!(x, y)
	lines!(x2, y, linestyle=:dash)

	scatter!(log10(nfw_ah.c), log10(nfw.c))
	fig
end

# ╔═╡ d438ef61-a565-401c-9e90-1c92e9e1d4d3
function plot_M_ρ!(nfw)
	t = 10 .^ LinRange(-3.5, 1.5, 1000) .*  LilGuys.R200(nfw)

	M = LilGuys.mass.(nfw, t)
	ρ = LilGuys.mean_density.(nfw, t)

	lines!(log10.(ρ), log10.(M))
	
end

# ╔═╡ e3acc2f0-a8dc-4158-b0d6-b0f1c31245e9
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xlabel = "log ρ",
			 ylabel = "log M"
			 )


	plot_M_ρ!(nfw_ah)

	fig
end

# ╔═╡ 689370d3-88a7-4f20-b4c7-63f7a835198e
scale_mass(nfw) = nfw.M_s * LilGuys.A_NFW(1)

# ╔═╡ 8fad5d11-a441-4cb5-a205-43893e8152ed
z_coll = LilGuys.find_zero(z -> scale_mass(nfw) - mass_given_density(nfw_ah, 200ρ_c(z)), 1)

# ╔═╡ f0e9cd0c-c220-4f09-a461-bda303c6e1c5
ρ_coll = ρ_c(z_coll)

# ╔═╡ 0872d86a-0a42-4323-913c-6f6e337e1f8b
scale_density(nfw) / ρ_coll 

# ╔═╡ 5ba06492-f3d7-48c5-bf95-52933cffd9a1
ρ_coll / ρ_0

# ╔═╡ 0130caca-b7a7-48cc-8935-29d67f6a2308
let
	fig = Figure()
	Ax = Axis(fig[1,1],
			 xlabel = L"\log \rho_\textrm{crit} / \rho_0",
			 ylabel = L"\log \bar\rho / \rho_0"
			 
		)


	x = LinRange(0, 2.5, 1000)
	ρ = 10 .^ x .* ρ_0

	ρ_m = A_ρ * ρ


	y = @. log10(ρ_m / ρ_0)
	lines!(x, y)

	scatter!(log10(ρ_coll / ρ_0), log10(scale_density(nfw) / ρ_0))

	fig
end

# ╔═╡ 001525fd-6c28-454f-a70b-14354471a177
scale_mass(nfw)

# ╔═╡ 246ce2c3-bd41-4fc4-8387-b6d768e97716
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			xlabel = "ρ c / ρ 0",
			  ylabel = "M"
			 )


	t = 10 .^ LinRange(-2.5, 0, 1000) .*  LilGuys.R200(nfw_ah)

	ρ = LilGuys.mean_density.(nfw_ah,t )

	x = log10.(ρ ./ 200ρ_0)
	y = log10.(LilGuys.mass.(nfw_ah, t)  ./ M200(nfw))
	lines!(x, y)


	scatter!(log10(1/A_ρ * scale_density(nfw) ./ ρ_0), log10(scale_mass(nfw) / M200(nfw)))

	scatter!(0, 0)
	fig
end

# ╔═╡ f7e4babc-ea6f-4cc3-bc31-b9837e248663
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xlabel = "-log 1 + z",
			 ylabel = "log M",
			#xreversed = true,
			 # limits=(0, nothing, nothing, nothing)
			 )


	scatter!(-log10(1+z_coll), log10.(scale_mass(nfw)))



	x = LinRange(-0.8, 0, 1000)
	a = 10 .^ x
	z = @.  1 / a - 1
	ρs = @. ρ_c(z)

	Ms = mass_given_density.(nfw_ah, 200ρs)

	@info Ms[end]
	lines!(x, log10.(Ms))

	fig
end

# ╔═╡ 0c46be62-efb1-4b71-b071-f7bf102363ad
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xlabel = "1 + z",
			 ylabel = "log M",
			#xreversed = true,
			 # limits=(0, nothing, nothing, nothing)
			 )


	scatter!(z_coll, log10.(scale_mass(nfw)))



	x = LinRange(-1.0, 0, 1000)
	a = 10 .^ x
	z = @.  1 / a - 1
	ρs = @.  ρ_c(z)

	Ms = mass_given_density.(nfw_ah, 200ρs)

	lines!(z, log10.(Ms))

	fig
end

# ╔═╡ 102236d3-f2b3-4174-b9ac-837deb9ca083
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			#xreversed = true,
			 # limits=(0, nothing, nothing, nothing)
			 )


	t = 10 .^ LinRange(-5, 1, 1000) .*  LilGuys.R200(nfw_ah)

	ρ = LilGuys.mean_density.(nfw_ah,t )

	lines!(log10.(t), log10.(ρ ./ ρ_coll .* t.^2))
	hlines!(ρ_coll)

	fig
end

# ╔═╡ bf6262c6-e858-4ecd-bfa6-91587d90704e
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			#xreversed = true,
			 # limits=(0, nothing, nothing, nothing)
			 )


	t = 10 .^ LinRange(-4, 0, 1000) .*  LilGuys.R200(nfw_ah)

	ρ = LilGuys.mean_density.(nfw,t )

	x = log10.(ρ ./ scale_density(nfw))
	y = log10.(LilGuys.mass.(nfw, t) ./ nfw.M_s)
	lines!(x, y)
	fig
end

# ╔═╡ d964b115-2183-40cc-bb00-f55d13aeadca
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			#xreversed = true,
			 # limits=(0, nothing, nothing, nothing)
			 )


	t = 10 .^ LinRange(-4, -1.5, 1000) .*  LilGuys.R200(nfw_ah)

	ρ = LilGuys.mean_density.(nfw_ah,t )

	x = log10.(ρ ./ scale_density(nfw))
	y = log10.(LilGuys.mass.(nfw_ah, t) ./ nfw_ah.M_s)
	lines!(x, y)



	ρ = LilGuys.mean_density.(nfw_ah,t )

	y2 = log10.(mass_given_density.(nfw_ah, ρ) ./ nfw_ah.M_s)
	
	lines!(x, y2)
	fig
end

# ╔═╡ 63d5ded7-8860-4d94-af48-84a45550c3aa
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			#xreversed = true,
			 # limits=(0, nothing, nothing, nothing)
			 )


	x = LinRange(0, 5, 100)
	ρ = ρ_c(0) * 10 .^ x

	y = mass_given_density.(nfw_ah, ρ)

	lines!(x, log10.(y ./ M200(nfw)))
	fig
end

# ╔═╡ a5e38f8c-c838-4fb0-8915-4acbb6b7f184
md"""
# Density evolution
"""

# ╔═╡ 662868fe-4a0d-40c2-9744-39efe781ee76
LilGuys.ρ_crit

# ╔═╡ 5c1e028c-be21-4431-9b63-135d926c3822
c_0 = nfw.c

# ╔═╡ 2beba66d-9a40-42f3-bf4d-0ec09985ac08
c_from_z(z) = LilGuys.Ludlow.c_ludlow(mass_given_density(nfw_ah, ρ_c(z)), z)

# ╔═╡ b005ca96-a516-4a9a-9afc-10c5be92843a
c_from_z(0)

# ╔═╡ ae9270cf-b8cf-4684-b559-59c1a43b324f
halo_at_z(z) = LilGuys.NFW(M200=mass_given_density(nfw_ah, ρ_c(z)), c=c_from_z(z), z=z)

# ╔═╡ 768cd933-9dbe-4ae7-b812-169e8e991050
LilGuys.Ludlow.c_ludlow(1, 7.7)

# ╔═╡ 44c39128-3f06-4cbf-8276-9e63649dfe19
halo_at_z(7)

# ╔═╡ 83c1fcad-ab1a-4fb2-8397-46ee4e8ca802
mass_given_density(nfw_ah, ρ_c(9))

# ╔═╡ 0963eb03-2c6d-45b9-b34d-344a0bc5b545
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "z", ylabel = "c")


	z = LinRange(0, 7, 100)
	c = c_from_z.(z)

	lines!(z, c)
	fig
end

# ╔═╡ e1dc5fee-600a-4a81-87cd-9db46ca223e0
let
	fig = Figure()
	ax = Axis(fig[1,1])


	zs = 0:7

	for z in zs
		h = halo_at_z(z)

		x = LinRange(-2, log10(R200(h, z)), 1000) 
		r = 10 .^ x
		
		v = LilGuys.v_circ.(h, r)

		y = log10.(v .* V2KMS)
		lines!(x, y, color=z, colorrange=extrema(zs))


		scatter!(log10(LilGuys.r_circ_max(h)), log10(LilGuys.v_circ_max(h)*V2KMS), 
				color=z, colorrange=extrema(zs))
		r_vir = R200(h, z)
		m_vir = mass(h, r_vir)
		ρ_vir = m_vir / (4π/3 * r_vir^3)
		@info m_vir, mass_given_density(nfw_ah, ρ_c(z))
		@info ρ_vir, 200ρ_c(z)
		@info LilGuys.v_circ_max(h)
	end

	Colorbar(fig[1,2], colorrange=(extrema(zs)))

	fig
end

# ╔═╡ faa9c7d6-2d5e-45d7-8cb0-3e47a025f080
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			  xlabel = "log radius / kpc",
			  ylabel = "log density"
			 )




	zs = 0:7

	for z in zs
		h = halo_at_z(z)

		x = LinRange(-2, log10(R200(h, z)), 1000) 
		r = 10 .^ x
		
		v = LilGuys.density.(h, r)

		y = log10.(v)
		lines!(x, y, color=z, colorrange=extrema(zs))
		#hlines!(log10(200ρ_c(z)), color=z, colorrange=(extrema(zs)))
		scatter!(log10(h.r_s), log10(LilGuys.density(h, h.r_s)), color=z, colorrange=extrema(zs))

	end

	Colorbar(fig[1,2], colorrange=(extrema(zs)))

	fig
end

# ╔═╡ 2e6f6253-8408-4bb1-bd5e-36c4a47a814f
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			  xlabel = "log M200",
			  ylabel = "log c"
			 )




	zs = [0, 1, 2, 3, 5, 7]

	for z in zs
	
		x = LinRange(-15, 5, 1000) 
		M = 10 .^ x		
		c = LilGuys.Ludlow.c_ludlow.(M, z)
		y = log10.(c)
		lines!(x, y, color=z, colorrange=extrema(zs))

	end

	Colorbar(fig[1,2], colorrange=(extrema(zs)))

	fig
end

# ╔═╡ 990c126b-83b4-4f3f-a2d6-d775d2174f1c
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			  xlabel = "log radius / R200",
			  ylabel = "log density"
			 )




	zs = 0:7

	for z in zs
		h = halo_at_z(z)

		x = LinRange(-2, log10(R200(h, z)), 1000) 
		r = 10 .^ x
		
		v = LilGuys.density.(h, r)

		y = log10.(v)

		scale = R200(h)
		
		lines!(x .- log10(scale), y .+ 3log10(scale), color=z, colorrange=extrema(zs))
		#hlines!(log10(200ρ_c(z)), color=z, colorrange=(extrema(zs)))
		#

	end

	Colorbar(fig[1,2], colorrange=(extrema(zs)))

	fig
end

# ╔═╡ 9e1b080d-ec12-4345-b56c-5f9c0e233fdc
ρ_c(3) * 4^-3

# ╔═╡ 292b53c4-c1e2-4572-99ea-65842248852e
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			  xlabel = "log comoving radius / ckpc",
			  ylabel = "log comoving density"
			 )




	zs = 0:7

	for z in zs
		h = halo_at_z(z)

		x = LinRange(-2, log10(R200(h, z)), 1000) 
		r = 10 .^ x
		
		v = LilGuys.density.(h, r)

		y = log10.(v)
		lines!(x .+ log10(1+z), y .- 3log10(1+z), color=z, colorrange=extrema(zs))
		#hlines!(log10(200ρ_c(z)), color=z, colorrange=(extrema(zs)))
		scatter!(log10(h.r_s) + log10(1+z), log10(LilGuys.density(h, h.r_s) ) - 3log10(1+z), color=z, colorrange=extrema(zs))

	end

	Colorbar(fig[1,2], colorrange=(extrema(zs)))

	fig
end

# ╔═╡ 36e585f8-fee4-4d7e-aa40-ad808b3382f2
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			  xlabel = "log radius / pkpc",
			  ylabel = "potential"
			 )




	zs = 0:7

	for z in zs
		h = halo_at_z(z)

		x = LinRange(-2, log10(R200(h, z)), 1000) 
		r = 10 .^ x
		
		v = LilGuys.potential.(h, r)

		lines!(x, v, color=z, colorrange=extrema(zs))


	end

	Colorbar(fig[1,2], colorrange=(extrema(zs)))

	fig
end

# ╔═╡ 0d3b11e3-ce38-47be-8ecd-9041bea887ad
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			  xlabel = "log comoving radius / pkpc",
			  ylabel = "comoving potential"
			 )




	zs = 0:7

	for z in zs
		h = halo_at_z(z)

		x = LinRange(-2, log10(R200(h, z)), 1000) 
		r = 10 .^ x
		
		v = LilGuys.potential.(h, r)

		lines!(x .- log10(1+z), v ./ (1+z), color=z, colorrange=extrema(zs))


	end

	Colorbar(fig[1,2], colorrange=(extrema(zs)))

	fig
end

# ╔═╡ Cell order:
# ╠═9991bd2a-77c4-11f0-2bf0-2fbe36425e63
# ╠═4f471846-ae8b-419e-8cd5-920032f19c82
# ╠═b88a3b23-6f3e-487b-848c-24fab0e8a286
# ╠═c8f78a5b-fbf5-4b0f-849a-eec18e8a4639
# ╠═5ff465a3-846e-4fac-b2c3-bbe92a5a068e
# ╠═cc36e1fc-cd68-45c8-a1dd-4a75e19e3057
# ╠═d690b0c4-3f9a-4c46-b356-f439982adbba
# ╠═74437045-5ee0-4c80-b8b9-90adf1a06f61
# ╠═f433807e-620a-4191-b94e-ebe9f68d3a21
# ╠═8325d421-a7a7-41ee-a0e3-3ff1b0df0c22
# ╠═cd58027b-ff2d-45d6-9b37-08d95eaf40ab
# ╟─ac6e60e7-b9b1-4c98-9144-061d25429c40
# ╠═4c59a8e8-1b59-4817-802a-79db329b1eb3
# ╠═6fa3783a-14b3-4943-a5d8-a04a48effd43
# ╠═fbafc7d5-cfde-4a79-ba98-77157d951cc9
# ╠═8fad5d11-a441-4cb5-a205-43893e8152ed
# ╠═f0e9cd0c-c220-4f09-a461-bda303c6e1c5
# ╠═0872d86a-0a42-4323-913c-6f6e337e1f8b
# ╠═72e2f6a1-1214-4cc5-bb7b-ca4148134f16
# ╠═5ba06492-f3d7-48c5-bf95-52933cffd9a1
# ╠═b88a0324-6a65-4bd3-8ee3-05e4fbd04a44
# ╠═17276a12-4365-487a-bcfe-28c56c4f6cb4
# ╠═e684cc0e-089d-46c1-b779-9d9437eba071
# ╠═001525fd-6c28-454f-a70b-14354471a177
# ╠═6cc3c980-a208-4484-a75f-05efd461c73f
# ╠═311ec430-fcf6-4fe1-b167-e548f5061ef9
# ╠═096ed683-3b2a-4a13-a7ec-f1300904d0e4
# ╠═de3d6614-0f58-4110-bcc0-ce61e2f57fae
# ╠═fdb40c9f-59a7-4403-94ae-cc7b09a11083
# ╠═2710b8d2-2339-46bf-9951-8b45b58d8d71
# ╠═0065b475-ed69-405e-8a83-08026535b9ab
# ╠═a8a26ba9-d0cc-437d-afaa-66ee16d483b9
# ╠═7144c3df-bde4-46a6-ba1c-744894691eb5
# ╟─87182b31-fd19-48f4-abaf-90d2cd5f2bb4
# ╠═423f83f7-3530-4689-9aa9-0fbc50df0e54
# ╠═2a222c23-5873-46de-b900-2abc8ab4d4a9
# ╟─c157c90e-aa07-43d0-bfcb-9a2d780dc05b
# ╠═246ce2c3-bd41-4fc4-8387-b6d768e97716
# ╠═0130caca-b7a7-48cc-8935-29d67f6a2308
# ╠═dde42a8a-e641-4eaf-a513-411e1fcf096e
# ╠═f7e4babc-ea6f-4cc3-bc31-b9837e248663
# ╠═0c46be62-efb1-4b71-b071-f7bf102363ad
# ╠═d438ef61-a565-401c-9e90-1c92e9e1d4d3
# ╠═e3acc2f0-a8dc-4158-b0d6-b0f1c31245e9
# ╠═689370d3-88a7-4f20-b4c7-63f7a835198e
# ╠═102236d3-f2b3-4174-b9ac-837deb9ca083
# ╠═bf6262c6-e858-4ecd-bfa6-91587d90704e
# ╠═d964b115-2183-40cc-bb00-f55d13aeadca
# ╠═63d5ded7-8860-4d94-af48-84a45550c3aa
# ╟─a5e38f8c-c838-4fb0-8915-4acbb6b7f184
# ╠═662868fe-4a0d-40c2-9744-39efe781ee76
# ╠═5c1e028c-be21-4431-9b63-135d926c3822
# ╠═2beba66d-9a40-42f3-bf4d-0ec09985ac08
# ╠═b005ca96-a516-4a9a-9afc-10c5be92843a
# ╠═ae9270cf-b8cf-4684-b559-59c1a43b324f
# ╠═768cd933-9dbe-4ae7-b812-169e8e991050
# ╠═44c39128-3f06-4cbf-8276-9e63649dfe19
# ╠═83c1fcad-ab1a-4fb2-8397-46ee4e8ca802
# ╠═0963eb03-2c6d-45b9-b34d-344a0bc5b545
# ╠═e1dc5fee-600a-4a81-87cd-9db46ca223e0
# ╠═faa9c7d6-2d5e-45d7-8cb0-3e47a025f080
# ╠═2e6f6253-8408-4bb1-bd5e-36c4a47a814f
# ╠═990c126b-83b4-4f3f-a2d6-d775d2174f1c
# ╠═9e1b080d-ec12-4345-b56c-5f9c0e233fdc
# ╠═292b53c4-c1e2-4572-99ea-65842248852e
# ╠═36e585f8-fee4-4d7e-aa40-ad808b3382f2
# ╠═0d3b11e3-ce38-47be-8ecd-9041bea887ad
