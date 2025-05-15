### A Pluto.jl notebook ###
# v0.20.8

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

# ╔═╡ 25f6474b-2ee7-4d2e-9ce5-fdf782adfa72
using Agama

# ╔═╡ dea35f32-ba19-474f-a7b9-fdaa5e89fc43
include("./style.jl")

# ╔═╡ 2a7db353-9975-41d0-8a03-d91622d549d4
potential_dir = joinpath(ENV["DWARFS_ROOT"], "agama", "potentials")

# ╔═╡ 5cbabcc8-6473-4ac1-8084-80da24a84617
Φ_ep = Agama.Potential(file = joinpath(potential_dir, "EP2020.ini"))

# ╔═╡ babbbe3b-13b7-4450-98be-f807f21706d8
Φ_bulge, Φ_thin, Φ_thick, Φ_nfw = Φ_ep._py

# ╔═╡ 86961f65-fa3e-4899-81d5-94c6bec3b8de
function plot_v_circ!(pot; x=LinRange(0, 2.5, 1000), kwargs...)
	r = 10 .^ x

	m = pot.enclosedMass(r) |> Agama.py2vec
	v = @. sqrt( m / r) * V2KMS

	lines!(x, v; kwargs...)
	
end

# ╔═╡ cf74dd35-70f3-4a88-b1c0-580d0ffe69aa
CairoMakie.activate!(type=:png)

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "v_circ_potential" let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = L"$\log\, R$ / kpc",
		ylabel = L"$v_\textrm{circ}$ / km\,s$^{-1}$",
		yscale=log10, 
		yticks = [10, 20, 100, 200],
		yminorticks = [10:10:100; 200:100:1000;],
		limits=(0, 2.5, nothing, nothing)
	)

	plot_v_circ!(Φ_ep._py, label="total", linestyle=:solid, color=:black, linewidth=8)

	plot_v_circ!(Φ_bulge, label="bulge", linestyle=:dot, color=COLORS[1], linewidth=6)
	plot_v_circ!(Φ_thin + Φ_thick, label = "disk", linestyle=:dashdot,  color=COLORS[2], linewidth=6)
	plot_v_circ!(Φ_nfw, label="halo", linestyle=:dash,  color=COLORS[3], linewidth=6)


	axislegend(position=:lb, margin=(19, 19, 19, 19))
	#Legend(fig[1,2], ax)
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═25f6474b-2ee7-4d2e-9ce5-fdf782adfa72
# ╠═dea35f32-ba19-474f-a7b9-fdaa5e89fc43
# ╠═2a7db353-9975-41d0-8a03-d91622d549d4
# ╠═5cbabcc8-6473-4ac1-8084-80da24a84617
# ╠═babbbe3b-13b7-4450-98be-f807f21706d8
# ╠═86961f65-fa3e-4899-81d5-94c6bec3b8de
# ╠═cf74dd35-70f3-4a88-b1c0-580d0ffe69aa
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
