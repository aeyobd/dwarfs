### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 62c47b3c-3699-11f0-16dd-5d8818928d6b
begin 
	import Pkg; Pkg.activate()

	using Agama
	using CairoMakie
	using Arya
end

# ╔═╡ 12111900-614e-4a81-8cbf-67c2a9a470e9
using OrderedCollections

# ╔═╡ 1a510a83-6109-4021-8bd8-5a2dc9f4a890
import LilGuys

# ╔═╡ 20f425a2-5f1a-4097-b522-f5d25cf1f772
h = LilGuys.NFW(v_circ_max=31/LilGuys.V2KMS, r_circ_max  = 4.2)

# ╔═╡ 4798201b-6a8a-44d9-bd0d-856c84830b2e
ρ_s = h.M_s / (4π * h.r_s^3)

# ╔═╡ 405535c6-bb35-46a2-8342-3d8eee7c2bd7
Rvir = LilGuys.R200(h)

# ╔═╡ 44d2b7c1-60ac-402d-ad2a-cc65c54332e7
Rvir / h.r_s

# ╔═╡ 09cdace6-b0c7-4db8-abfd-6cea3f0604ef
nfw = Agama.Potential(type="spheroid", alpha=1, beta=3, gamma=1, densityNorm=ρ_s, scaleRadius=h.r_s)

# ╔═╡ 4442735e-112f-4874-8d62-5ca7970f2587
nfw_exptrunc = Agama.Potential(type="spheroid", alpha=1, beta=3, gamma=1, densityNorm=ρ_s, scaleRadius=h.r_s,
	outerCutoffRadius=Rvir, cutoffStrength=1
)

# ╔═╡ 14e61f2b-71a9-4219-94e2-f6c41b1380d2
nfw_cubeexptrunc = Agama.Potential(type="spheroid", alpha=1, beta=3, gamma=1, densityNorm=ρ_s, scaleRadius=h.r_s,
	outerCutoffRadius=20h.r_s, cutoffStrength=3
)

# ╔═╡ 8ac25b73-3052-4669-ba42-b7eadd25b013
nfw_longexp = Agama.Potential(type="spheroid", alpha=1, beta=3, gamma=1, densityNorm=ρ_s, scaleRadius=h.r_s,
	outerCutoffRadius=100h.r_s, cutoffStrength=1
)

# ╔═╡ 88faf52b-984e-4c3b-aace-51aa02951420
nfw_sharp = Agama.Potential(type="spheroid", alpha=1, beta=3, gamma=1, densityNorm=ρ_s, scaleRadius=h.r_s,
	outerCutoffRadius=Rvir, cutoffStrength=10
)

# ╔═╡ a469c963-8f5a-4973-ab0c-914bebe53465
potentials = OrderedDict(
	"NFW" => nfw,
	"exp" => nfw_exptrunc,
	"cube" => nfw_cubeexptrunc,
	"long exp" => nfw_longexp,
	"sharp" => nfw_sharp,
)

# ╔═╡ 3c01ae79-56ed-4c52-8330-6202cbfb37c4
CairoMakie.activate!(type=:png)

# ╔═╡ d3e1ca46-9a25-4c08-b0e0-3c528b2cb9a2
let
	fig = Figure()
	ax = Axis(fig[1,1],
		ylabel = "log rho"
	)

	x = LinRange(-2, log10(2Rvir), 1000)

	for (label, pot) in potentials
		r = 10 .^ x
		y =  log10.(max.(0, Agama.density(pot, r' .* [1, 0, 0])))

		lines!(x, y, label=label)
	end
	hidexdecorations!(ticks=false, minorticks=false)

	axislegend(position=:lb)
	vlines!(log10(Rvir))
	vlines!(log10(h.r_s))

	ax_res = Axis(fig[2, 1], ylabel="res", xlabel="log r",
				  limits=(nothing, (-0.2, 0.05))
		)

	r = 10 .^ x
	y2 = log10.(Agama.density(nfw, r' .* [1, 0, 0]))

	for (label, pot) in potentials
		y =  log10.(max.(1e-10, Agama.density(pot, r' .* [1, 0, 0])))
		lines!(x, y .- y2, label=label)
	end

	linkxaxes!(ax, ax_res)
	rowsize!(fig.layout, 2, Relative(1/4))
	
	fig
end

# ╔═╡ d1536ac1-65fa-4c70-b252-ed620ad9f47b
10^-0.04

# ╔═╡ d9a9a0fd-46ee-413b-983f-1295f220d0f2
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log r",
		ylabel = "v circ"
	)

	x = LinRange(-2, log10(2Rvir), 1000)

	for (label, pot) in potentials
		r = 10 .^ x
		y =  Agama.circular_velocity(pot, r)

		lines!(x, y .* LilGuys.V2KMS, label=label)
	end
	hidexdecorations!(ticks=false, minorticks=false)

	axislegend(position=:lb)
	vlines!(log10(Rvir))
	vlines!(log10(h.r_s))

	ax_res = Axis(fig[2, 1], ylabel="res", xlabel="log r",
				  limits=(nothing, (-0.1, 0.05))
				 )

	r = 10 .^ x
	y2 = Agama.circular_velocity(nfw, r)

	for (label, pot) in potentials
		y =  Agama.circular_velocity(pot, r)
		lines!(x, log10.(y) .- log10.(y2), label=label)
	end

	linkxaxes!(ax, ax_res)
	rowsize!(fig.layout, 2, Relative(1/4))
	
	fig
end

# ╔═╡ a634a872-d30c-4440-8581-779a25fe59c3
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log r",
		ylabel = "enclosed mass"
	)

	x = LinRange(-2, log10(20Rvir), 1000)

	for (label, pot) in potentials
		r = 10 .^ x
		y =  Agama.enclosed_mass(pot, r)

		lines!(x, y, label=label)
	end

	axislegend(position=:lb)
	vlines!(log10(Rvir))
	vlines!(log10(h.r_s))

	
	fig
end

# ╔═╡ d7c6b9d8-0b62-41c9-a930-9189cbaf7e93
for (label, pot) in potentials
	println(label, "\t", enclosed_mass(pot, 1000Rvir))
end

# ╔═╡ 7e607d22-0c8b-4dd0-ad16-c72d8386c049
for (label, pot) in potentials
	println(label, "\t", enclosed_mass(pot, 1Rvir))
end

# ╔═╡ 596e747f-1ddb-4fd3-a358-054842575a98
20/1.1

# ╔═╡ 13073a02-4ae6-41b6-8514-4a8393c1cb7c
for (label, pot) in potentials
	println(label, "\t", enclosed_mass(pot, Rvir/2) / enclosed_mass(nfw, Rvir/2))
end

# ╔═╡ 1c37d00d-27a5-41db-816d-a58b8dac7671
for (label, pot) in potentials
	mtot = enclosed_mass(pot, 1000Rvir)
	m_in = enclosed_mass(pot, Rvir)
	println(label, "\t",  (m_in)/mtot)
end

# ╔═╡ Cell order:
# ╠═62c47b3c-3699-11f0-16dd-5d8818928d6b
# ╠═12111900-614e-4a81-8cbf-67c2a9a470e9
# ╠═1a510a83-6109-4021-8bd8-5a2dc9f4a890
# ╠═20f425a2-5f1a-4097-b522-f5d25cf1f772
# ╠═4798201b-6a8a-44d9-bd0d-856c84830b2e
# ╠═405535c6-bb35-46a2-8342-3d8eee7c2bd7
# ╠═44d2b7c1-60ac-402d-ad2a-cc65c54332e7
# ╠═09cdace6-b0c7-4db8-abfd-6cea3f0604ef
# ╠═4442735e-112f-4874-8d62-5ca7970f2587
# ╠═14e61f2b-71a9-4219-94e2-f6c41b1380d2
# ╠═8ac25b73-3052-4669-ba42-b7eadd25b013
# ╠═88faf52b-984e-4c3b-aace-51aa02951420
# ╠═a469c963-8f5a-4973-ab0c-914bebe53465
# ╠═3c01ae79-56ed-4c52-8330-6202cbfb37c4
# ╠═d3e1ca46-9a25-4c08-b0e0-3c528b2cb9a2
# ╠═d1536ac1-65fa-4c70-b252-ed620ad9f47b
# ╠═d9a9a0fd-46ee-413b-983f-1295f220d0f2
# ╠═a634a872-d30c-4440-8581-779a25fe59c3
# ╠═d7c6b9d8-0b62-41c9-a930-9189cbaf7e93
# ╠═7e607d22-0c8b-4dd0-ad16-c72d8386c049
# ╠═596e747f-1ddb-4fd3-a358-054842575a98
# ╠═13073a02-4ae6-41b6-8514-4a8393c1cb7c
# ╠═1c37d00d-27a5-41db-816d-a58b8dac7671
