### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# ╔═╡ 28ed6182-698c-11f0-3368-f5b5af3ea05c
begin
	import Pkg; Pkg.activate()


	using LilGuys
	using CairoMakie, Arya

end

# ╔═╡ ec95bcc8-2c24-40d7-80f1-32e83c4c17aa
include("./paper_style.jl")

# ╔═╡ 7efef2ac-d002-4a86-a8c3-e9e7036de668
FIGDIR = "./figures"

# ╔═╡ 4671142c-e341-4927-bfe7-809b41daf251
CairoMakie.activate!(type=:png)

# ╔═╡ 903cd96e-0768-4edb-9126-4b6eafeb4ca9
prof_exp = LilGuys.Exp2D()

# ╔═╡ cb7dc629-fbc9-4e08-ba56-2074009ee874
prof_sersic_n0_8 = LilGuys.Sersic(n=0.8)

# ╔═╡ d8694866-324a-4bac-8856-77a130559a8f
prof_sersic_n1_2 = LilGuys.Sersic(n=1.2)

# ╔═╡ 5fb7529d-0b1a-4fa1-bbaf-ed13f2f01027
prof_plummer = LilGuys.Plummer()

# ╔═╡ 645b8702-c0df-4152-b6c3-7e21f11ab305
import SpecialFunctions: gamma

# ╔═╡ 929dfb76-3580-49e2-bf92-098e131f4c01
Base.@kwdef struct GeneralizedPlummer <: LilGuys.SphericalProfile
	r_s::Float64
	M::Float64
	n::Float64

	_k::Float64 = gamma(n / 2) / (gamma((n-3)/2) * π^(3/2))
	
end

# ╔═╡ 6e612d6d-4605-47de-af5f-c38d022650d8
LilGuys.density(prof::GeneralizedPlummer, r::Real) = prof._k * prof.M * prof.r_s^(prof.n - 3) / (prof.r_s^2 + r^2)^(prof.n/2)

# ╔═╡ cd2f1f91-2f9d-476d-8f84-6e17abf6df66
outer_slope(n) = n - 1

# ╔═╡ f2725829-5ca7-43a5-b43e-1dfddbac1f7e
LilGuys.R_h(prof::LilGuys.Sersic) = prof.R_h

# ╔═╡ 15ed2abd-a8b4-4bc4-94ab-48276475e9a8
steep_plummer = GeneralizedPlummer(M=1, r_s=1, n=9)

# ╔═╡ c0e9018f-5ac6-4f1c-b800-bb7348a4f5f7
profiles = [
	"Exponential" => prof_exp,
	"Sérsic (𝑛=0.8)" => prof_sersic_n0_8,
	"Sérsic (𝑛=1.2)" => prof_sersic_n1_2,
	"Plummer" => prof_plummer,
	"Steeper Plummer" => steep_plummer,
	"King (𝑐=10)" => LilGuys.KingProfile(R_c=1, R_t=10, k=1)
]

# ╔═╡ 76a2ba8f-81c4-405e-a9f9-b99c34541a5e
LilGuys.find_zero(R -> LilGuys.mass_2D(GeneralizedPlummer(M=1, r_s=1, n=9), R) .- 1/2, 1)

# ╔═╡ 12ab6677-b944-4ac7-941c-5699598e2ba2
LilGuys.surface_density(steep_plummer, 0.0001)

# ╔═╡ c1aaca51-d80d-4c9f-a3b5-86d38944a348


# ╔═╡ 1c4d83c9-a6db-4e11-a0e4-18bef5894427
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		ylabel = "log Sigma"
			 )
	ax2 = Axis(fig[2,1], 
		xlabel = "log R / Rh",
		ylabel = "Gamma"
			 )

	x = LinRange(-1, 1, 1000)

	R = 10 .^ x

	for (label, prof) in profiles

		R_h = LilGuys.R_h(prof)
		y = @. log10(LilGuys.surface_density(prof, R * R_h))

		y .-= log10(LilGuys.surface_density(prof, R_h))

		lines!(ax, x, y, label=label)

		Gamma = LilGuys.gradient(y, x)

		lines!(ax2, x, Gamma, label=label)

	end

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ ff84653f-4f4b-4631-bb73-fcb3ce03aef6
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		ylabel = "log Sigma"
			 )
	ax2 = Axis(fig[2,1], 
		xlabel = "R / Rh",
		ylabel = "Gamma"
			 )

	x = LinRange(-1, 1.0, 1000)

	R = 10 .^ x

	for (label, prof) in profiles

		R_h = LilGuys.R_h(prof)
		y = @. log10(LilGuys.surface_density(prof, R * R_h))

		y .-= log10(LilGuys.surface_density(prof, R_h))

		lines!(ax, R, y, label=label)

		Gamma = LilGuys.gradient(y, x)

		lines!(ax2, R, Gamma, label=label)

	end

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 9a087611-eb84-4fff-b945-df6585445b1e
@savefig "analytic_profile_comparison" let
	fig = Figure(figure_padding=8)
	ax = Axis(fig[1,1], 
		ylabel = L"\log\, \Sigma"
			 )
	ax2 = Axis(fig[2,1], 
		xlabel = "R / Rh",
		ylabel = L"\log\left( \Sigma - \Sigma_\textrm{exp}\right)",
			   limits=(-1, 1, -2, 2)
			 )

	x = LinRange(-1, 1.0, 2000)

	R = 10 .^ x
	prof_ref = prof_exp

	for (label, prof) in profiles

		R_h = LilGuys.R_h(prof)
		y = @. log10(LilGuys.surface_density(prof, R * R_h))

		y .-= log10(LilGuys.surface_density(prof, R_h))

		R_h_ref = LilGuys.R_h(prof_ref)
		y_ref = @. log10(LilGuys.surface_density(prof_ref, R * R_h_ref))
		y_ref .-= log10(LilGuys.surface_density(prof_ref, R_h_ref))

		lines!(ax, x, y, label=label, linewidth=2)

		lines!(ax2, x, y .- y_ref, label=label, linewidth=2)

	end

	hidexdecorations!(ax, ticks=false, minorticks=false)
	rowgap!(fig.layout, 0)
	rowsize!(fig.layout, 2, Relative(1/4))


	linkxaxes!(ax, ax2)
	axislegend(ax, position=:lb)
	fig
end

# ╔═╡ 2546b28b-72c6-4a8a-9851-49592c3e5eb1
pwd()

# ╔═╡ Cell order:
# ╠═28ed6182-698c-11f0-3368-f5b5af3ea05c
# ╠═7efef2ac-d002-4a86-a8c3-e9e7036de668
# ╠═ec95bcc8-2c24-40d7-80f1-32e83c4c17aa
# ╠═4671142c-e341-4927-bfe7-809b41daf251
# ╠═903cd96e-0768-4edb-9126-4b6eafeb4ca9
# ╠═cb7dc629-fbc9-4e08-ba56-2074009ee874
# ╠═d8694866-324a-4bac-8856-77a130559a8f
# ╠═5fb7529d-0b1a-4fa1-bbaf-ed13f2f01027
# ╠═c0e9018f-5ac6-4f1c-b800-bb7348a4f5f7
# ╠═645b8702-c0df-4152-b6c3-7e21f11ab305
# ╠═929dfb76-3580-49e2-bf92-098e131f4c01
# ╠═6e612d6d-4605-47de-af5f-c38d022650d8
# ╠═cd2f1f91-2f9d-476d-8f84-6e17abf6df66
# ╠═f2725829-5ca7-43a5-b43e-1dfddbac1f7e
# ╠═15ed2abd-a8b4-4bc4-94ab-48276475e9a8
# ╠═76a2ba8f-81c4-405e-a9f9-b99c34541a5e
# ╠═12ab6677-b944-4ac7-941c-5699598e2ba2
# ╠═c1aaca51-d80d-4c9f-a3b5-86d38944a348
# ╠═1c4d83c9-a6db-4e11-a0e4-18bef5894427
# ╠═ff84653f-4f4b-4631-bb73-fcb3ce03aef6
# ╠═9a087611-eb84-4fff-b945-df6585445b1e
# ╠═2546b28b-72c6-4a8a-9851-49592c3e5eb1
