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

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
md"""
# Setup
"""

# ╔═╡ b7a4e67b-a3c9-44c2-acf3-b61f60d387ad
scale_theme_element!(:linewidth, 0.5)

# ╔═╡ 48289b64-08a8-4b26-9506-1843b3c69009
scale_theme_element!(:markersize, 2)

# ╔═╡ 7587767d-b2f7-4a18-b4a3-8ca410078dd8
10^(1.5 - 0.04)

# ╔═╡ 22629f89-3d2a-457d-b80a-7df412b9c791
halos = OrderedDict(
	:mean => NFW(v_circ_max = 31 / V2KMS, r_circ_max = 5.9),
	:compact => NFW(v_circ_max = 31 / V2KMS, r_circ_max = 3.2),
	:middle => NFW(v_circ_max = 31 / V2KMS, r_circ_max = 4.2),
	:smol => NFW(v_circ_max=25/V2KMS, r_circ_max=2.5)
)

# ╔═╡ 0dc48789-f708-4a4d-9731-6dbfe94b0913
begin
	M200_mean = 10 .^ LinRange(-2, 1.5, 1000)
	c_mean2 = LilGuys.Ludlow.c_ludlow.(M200_mean, 0)
	halo_mean = [NFW(M200=M200_mean[i], c=c_mean2[i]) for i in eachindex(M200_mean)]
	Vc_mean = v_circ_max.(halo_mean)
	Rc_mean = r_circ_max.(halo_mean)
	Ms_mean = LilGuys.M_s_from_vel_fattahi.(Vc_mean)
end

# ╔═╡ e4c5705b-4547-4679-a669-f51a3217f578
log10(31.5)

# ╔═╡ f178ab03-06ba-4e83-9be9-3924543bd288
Σ = [0.019467512114446032 0.0019034735328083907; 0.0019034735328083907 0.001516776581639168]

# ╔═╡ 28f7a671-8b64-43ed-b079-95cd7360a73d
μ = 
[0.780442
1.49899]

# ╔═╡ e470d921-fec4-436c-9cf4-f7d74d142e3c


# ╔═╡ 325ff1ea-ddfa-4661-9dd0-22866afca67d
import LinearAlgebra: eigen, diagm

# ╔═╡ 1c0b9737-217a-4870-bc70-1a5959fe2b01
λ, V = eigen(Σ)

# ╔═╡ f0e7895a-e7c5-4db2-bfdc-716121211612
function ellipse!(sigma; kwargs...)
	t = LinRange(0, 2π, 10000)

	Λ = diagm(λ)
	xy = [μ .+ V * sqrt.(sigma * Λ) * [cos(tt), sin(tt)] for tt in t]

	lines!(10 .^ first.(xy), 10 .^ last.(xy); kwargs...)
end

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

	v = (Vc_mean * V2KMS)
	lines!((Rc_mean), v, color=:grey)
	
	xl = (LilGuys.Ludlow.solve_rmax.(Vc_mean, 0.09))
	xh =  (LilGuys.Ludlow.solve_rmax.(Vc_mean, -0.09))
	x = [xl; reverse(xh)]
	y = [v; reverse(v)]
	poly!(x, y, color=(:grey,0.2))
	scatter!(x, y, markersize=1)

	hspan!(10^(μ[2] - sqrt(Σ[2,2])), 10^(μ[2] + sqrt(Σ[2,2])), color=(COLORS[2], 0.1))
	hlines!(10^(μ[2]), color=(COLORS[2], 0.5))

	
	
	for (label, halo) in halos
		y = LilGuys.v_circ_max(halo) * V2KMS
		x = LilGuys.r_circ_max(halo) 
		scatter!((x), (y), label=string(label))
	end


		

	ellipse!(2.3, color=:green, alpha=0.1)
	ellipse!(6.17, color=:green, alpha=0.1)
	ellipse!(11.8, color=:green, alpha=0.1)

	
	axislegend(position=:lt, title="halo")

	fig
end

# ╔═╡ b01d7319-2f22-49a5-9303-2ef28b24b798
stars = LilGuys.Exp2D(R_s=0.15)

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

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═f653a17a-8062-4898-9d0c-168a13dfbf6f
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═b7a4e67b-a3c9-44c2-acf3-b61f60d387ad
# ╠═48289b64-08a8-4b26-9506-1843b3c69009
# ╠═7587767d-b2f7-4a18-b4a3-8ca410078dd8
# ╠═22629f89-3d2a-457d-b80a-7df412b9c791
# ╠═0dc48789-f708-4a4d-9731-6dbfe94b0913
# ╠═e4c5705b-4547-4679-a669-f51a3217f578
# ╠═f178ab03-06ba-4e83-9be9-3924543bd288
# ╠═28f7a671-8b64-43ed-b079-95cd7360a73d
# ╠═e470d921-fec4-436c-9cf4-f7d74d142e3c
# ╠═325ff1ea-ddfa-4661-9dd0-22866afca67d
# ╠═1c0b9737-217a-4870-bc70-1a5959fe2b01
# ╠═f0e7895a-e7c5-4db2-bfdc-716121211612
# ╠═a15dbc88-8dbf-48d4-8a72-1881151c8e26
# ╠═b01d7319-2f22-49a5-9303-2ef28b24b798
# ╠═e0e418d2-16f2-4bbe-92db-ccdc6e89669e
# ╠═a085b24c-9ac6-41c8-811e-8e972eef0efb
# ╠═51075ee8-2cd7-4643-bf45-24a8a60d4eec
