### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 41c8b302-32a3-11f0-2a7a-d5be92a0c4b4
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Agama
	using CairoMakie
	using Arya
end

# ╔═╡ 182303d1-663b-436e-b246-21ecec1613b6
md"""
Plots the circular velocity curve of the given potential.
"""

# ╔═╡ ba8f9b70-a324-478d-8376-4c5ff32a92cc
begin
	potential_file = "hunter+2024/MWPotentialHunter24_axi.ini"
	#potential_file = "portail+2017.ini"
end

# ╔═╡ 9f440b9d-6087-4017-83f6-95505b5189f5
units = :ones

# ╔═╡ 3abeb382-fbed-445d-89d9-7feb696c1c9f
Φ = Agama.Potential(file = joinpath(ENV["DWARFS_ROOT"], "agama/potentials", potential_file))

# ╔═╡ 12abf249-3a2c-4165-9b7e-688eff82eee6
Φ_ref = Agama.Potential(file = joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))

# ╔═╡ 879504bc-4a4c-4b83-a12e-aa136d6445c8
if units == :ones
	M_scale = 1 #4.301476e-7
	V_scale = 0.004822
	R_scale = 1.0
end

# ╔═╡ 85ac46cd-f140-4020-9efa-607afe66aab9
T_scale = R_scale / V_scale

# ╔═╡ b87153e3-70d6-40ec-831a-0cd512fe663f
A_scale =  V_scale / T_scale * M_scale

# ╔═╡ e4ee804c-9524-4809-bf93-079effe31af1
Pot_scale = A_scale  * R_scale * M_scale

# ╔═╡ 770287e7-8fab-473b-8a41-867668697ab4
CairoMakie.activate!(type=:png)

# ╔═╡ 5b9af774-f356-46a0-aa96-96633ecd4b8a
let
	fig = Figure()
	ax = Axis(fig[1,1])

	x = LinRange(0, 100, 1000)

	m = Φ._py.enclosedMass(x) |> py2vec

	vs = V_scale .* sqrt.(m ./ x)
	lines!(x, vs * V2KMS)

	m = Φ_ref._py.enclosedMass(x) |> py2vec
	vs = sqrt.(m ./ x)
	lines!(x, vs * V2KMS)
	
	fig

end

# ╔═╡ d36b5a2e-de28-4c7b-b7cb-ea5daa26471d
x_vec = [1, 1, 1]

# ╔═╡ c5ea9e03-9916-45ef-b56f-569716aa2803


# ╔═╡ 80c256fd-1be3-45a7-bc6c-5e73ec33c2f8
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 limits=(0, 100, -5, 0))

	x = LinRange(0, 100, 1000)

	a = Agama.calc_acceleration(Φ, x' .* x_vec)

	lines!(x, log10.(radii(a) .* A_scale))


	a = Agama.calc_acceleration(Φ_ref, x' .* x_vec)

	lines!(x, log10.(radii(a)))

	
	fig

end

# ╔═╡ 745c3a87-b62a-4a38-81fd-acd4623a71b8
let
	fig = Figure()
	ax = Axis(fig[1,1])

	x = LinRange(0, 100, 1000)

	a = Agama.calc_Φ(Φ, x' .* x_vec)

	lines!(x, a .* Pot_scale)


	a = Agama.calc_Φ(Φ_ref, x' .* x_vec)

	lines!(x, (a))

	
	fig

end

# ╔═╡ 29ab5fc8-f7ed-4a67-b54f-3f72c104e901
md"""
# Project surface density
"""

# ╔═╡ 86b9a237-8271-4d5d-9c99-30fae9f6cd00
function plot_surface_density(Φ;
		N = 300,
		Rmax = 20,
		Σ_min = 1e-5,
		time = 0,
		β = 0,
	)
	
	fig = Figure()
	ax = Axis(fig[1,1], 
			  xlabel = "x' / kpc", ylabel = "y' / kpc", aspect=DataAspect())
	
	x = LinRange(-Rmax, Rmax, N)
	y = LinRange(-Rmax, Rmax, N)

	pos = [repeat(x, N) repeat(y, inner=N)]

	density = Φ._py.projectedDensity(pos, beta=β, t=time / T2GYR / T_scale) |> py2vec

	density = reshape(density, N, N)
	p = image!(extrema(x), extrema(y), log10.(max.(density, Σ_min)))

	Colorbar(fig[1,2], p, label=L"\log\, \Sigma")
	
	fig
end

# ╔═╡ 14699a06-5158-43df-8ee0-a9fc740d88e3
plot_surface_density(Φ)

# ╔═╡ 22368d3a-109e-4988-9f86-fd11990e467f
plot_surface_density(Φ, time=0.01)

# ╔═╡ d8f51813-9317-4454-aa3a-c0a8d05a82f2
plot_surface_density(Φ, β = π/2)

# ╔═╡ Cell order:
# ╠═182303d1-663b-436e-b246-21ecec1613b6
# ╠═ba8f9b70-a324-478d-8376-4c5ff32a92cc
# ╠═9f440b9d-6087-4017-83f6-95505b5189f5
# ╠═41c8b302-32a3-11f0-2a7a-d5be92a0c4b4
# ╠═3abeb382-fbed-445d-89d9-7feb696c1c9f
# ╠═12abf249-3a2c-4165-9b7e-688eff82eee6
# ╠═879504bc-4a4c-4b83-a12e-aa136d6445c8
# ╠═85ac46cd-f140-4020-9efa-607afe66aab9
# ╠═b87153e3-70d6-40ec-831a-0cd512fe663f
# ╠═e4ee804c-9524-4809-bf93-079effe31af1
# ╠═770287e7-8fab-473b-8a41-867668697ab4
# ╠═5b9af774-f356-46a0-aa96-96633ecd4b8a
# ╠═d36b5a2e-de28-4c7b-b7cb-ea5daa26471d
# ╠═c5ea9e03-9916-45ef-b56f-569716aa2803
# ╠═80c256fd-1be3-45a7-bc6c-5e73ec33c2f8
# ╠═745c3a87-b62a-4a38-81fd-acd4623a71b8
# ╠═29ab5fc8-f7ed-4a67-b54f-3f72c104e901
# ╠═86b9a237-8271-4d5d-9c99-30fae9f6cd00
# ╠═14699a06-5158-43df-8ee0-a9fc740d88e3
# ╠═22368d3a-109e-4988-9f86-fd11990e467f
# ╠═d8f51813-9317-4454-aa3a-c0a8d05a82f2
