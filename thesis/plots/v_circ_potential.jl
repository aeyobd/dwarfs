### A Pluto.jl notebook ###
# v0.20.17

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

# ╔═╡ 874f376a-c16b-411c-a260-d925e25d1277
using CSV, DataFrames

# ╔═╡ 25f6474b-2ee7-4d2e-9ce5-fdf782adfa72
using Agama

# ╔═╡ dea35f32-ba19-474f-a7b9-fdaa5e89fc43
include("./paper_style.jl")

# ╔═╡ 2a7db353-9975-41d0-8a03-d91622d549d4
potential_dir = joinpath(ENV["DWARFS_ROOT"], "agama", "potentials")

# ╔═╡ 5cbabcc8-6473-4ac1-8084-80da24a84617
Φ_ep = Agama.Potential(file = joinpath(potential_dir, "EP2020.ini"))

# ╔═╡ babbbe3b-13b7-4450-98be-f807f21706d8
Φ_bulge, Φ_thin, Φ_thick, Φ_nfw = Φ_ep._py

# ╔═╡ 4c0c9322-186e-4ba3-a98e-81e3c3c1b037
function v_circ(pot, r)
	m = pot.enclosedMass(r) |> Agama.py2f
	v = @. sqrt( m / r) * V2KMS
	return v
end

# ╔═╡ 86961f65-fa3e-4899-81d5-94c6bec3b8de
function plot_v_circ!(pot; x=LinRange(0, 2.5, 1000), kwargs...)
	r = 10 .^ x

	m = pot.enclosedMass(r) |> Agama.py2vec
	v = @. sqrt( m / r) * V2KMS

	lines!(r, v; kwargs...)
	
end

# ╔═╡ cf74dd35-70f3-4a88-b1c0-580d0ffe69aa
CairoMakie.activate!(type=:png)

# ╔═╡ 0b3729e3-184a-4eed-91df-67a8fee539b9
theme(:fontsize)

# ╔═╡ 33c4f27a-9bf7-404f-9acd-03fc2aa16827
eilers19 = CSV.read("eilers+19.dat", DataFrame, delim=' ')

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "v_circ_potential" let
	fig = Figure(figure_padding=12)

	ax = Axis(fig[1,1], 
		xlabel = "radius / kpc",
		ylabel = L"circular velocity / km\,s$^{-1}$",
		limits=(0, 200, nothing, nothing)
	)

	plot_v_circ!(Φ_ep._py, label="total", linestyle=:solid, color=:black, linewidth=2*theme(:linewidth)[])

	plot_v_circ!(Φ_bulge, label="bulge", linestyle=:dot, color=COLORS[1])
	plot_v_circ!(Φ_thin + Φ_thick, label = "disk", linestyle=:dashdot,  color=COLORS[2])
	plot_v_circ!(Φ_nfw, label="halo", linestyle=:dash,  color=COLORS[3])

	
	smallfontsize=theme(:fontsize)[]
	dy = 2
	text!(100, 191 + dy, text="total", fontsize=smallfontsize,  rotation=-0π/15, font=:bold)
	text!(100, 180 - dy, text="halo", fontsize=smallfontsize, color=COLORS[3], align=(:left, :top))
	text!(100, 57 + dy, text="disk", fontsize=smallfontsize, color=COLORS[2], rotation=0)
	text!(100, 29 - dy, text="bulge", fontsize=smallfontsize, color=COLORS[1], align=(:left, :top), rotation=0)

	#Legend(fig[1,2], ax)
	fig
end

# ╔═╡ e83f16b2-3d84-4a67-90f3-ffe9b439a863
v_circ(Φ_ep._py, 100.)

# ╔═╡ d01934f7-ab69-4715-b68c-e53709c57dee
v_circ(Φ_bulge, 100.)

# ╔═╡ d9699355-d786-436b-8811-ec8b8caff091
v_circ(Φ_thin + Φ_thick, 100.)

# ╔═╡ c4bc926f-b600-47a6-8df7-070ee57adca7
v_circ(Φ_nfw, 100.)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═874f376a-c16b-411c-a260-d925e25d1277
# ╠═25f6474b-2ee7-4d2e-9ce5-fdf782adfa72
# ╠═dea35f32-ba19-474f-a7b9-fdaa5e89fc43
# ╠═2a7db353-9975-41d0-8a03-d91622d549d4
# ╠═5cbabcc8-6473-4ac1-8084-80da24a84617
# ╠═babbbe3b-13b7-4450-98be-f807f21706d8
# ╠═4c0c9322-186e-4ba3-a98e-81e3c3c1b037
# ╠═86961f65-fa3e-4899-81d5-94c6bec3b8de
# ╠═cf74dd35-70f3-4a88-b1c0-580d0ffe69aa
# ╠═0b3729e3-184a-4eed-91df-67a8fee539b9
# ╠═33c4f27a-9bf7-404f-9acd-03fc2aa16827
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═e83f16b2-3d84-4a67-90f3-ffe9b439a863
# ╠═d01934f7-ab69-4715-b68c-e53709c57dee
# ╠═d9699355-d786-436b-8811-ec8b8caff091
# ╠═c4bc926f-b600-47a6-8df7-070ee57adca7
