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

	lines!(r, v; kwargs...)
	
end

# ╔═╡ a2e7a222-cca9-47a3-901e-69895193a1db
function plot_v_circ_agama!(pot; x=LinRange(0, 2.5, 1000), kwargs...)
	r = 10 .^ x

	m = pot.enclosedMass(r) |> Agama.py2vec
	v = @. sqrt( m / r) 

	lines!(r, v; kwargs...)
	
end

# ╔═╡ cf74dd35-70f3-4a88-b1c0-580d0ffe69aa
CairoMakie.activate!(type=:png)

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "v_circ_potential" let
	fig = Figure(size=(8*72, 8*72))

	ax = Axis(fig[1,1], 
		xlabel = L"$R$ / kpc",
		ylabel = L"$v_\textrm{circ}$ / km\,s$^{-1}$",
		limits=(0, 100, nothing, nothing),
	)
	lw = theme(:linewidth)

	plot_v_circ!(Φ_ep._py, label="total", linestyle=:solid, color=:black, linewidth=lw,)

	plot_v_circ!(Φ_bulge, label="bulge", linestyle=:dash, color=COLORS[1], linewidth=lw)
	plot_v_circ!(Φ_thin + Φ_thick, label = "disk", linestyle=:dashdot,  color=COLORS[2], linewidth=lw)
	plot_v_circ!(Φ_nfw, label="halo", linestyle=:dot,  color=COLORS[3], linewidth=lw)


	smallfontsize=0.8*theme(:fontsize)[]
	text!(50, 212, text="total", fontsize=smallfontsize,  rotation=-π/15)
	text!(50, 188, text="halo", fontsize=smallfontsize, color=COLORS[3], align=(:left, :top))
	text!(50, 80, text="disk", fontsize=smallfontsize, color=COLORS[2], rotation=-π/12)
	text!(50, 41, text="bulge", fontsize=smallfontsize, color=COLORS[1], align=(:left, :top), rotation=-π/24)

	resize_to_layout!(fig)
	fig
end

# ╔═╡ d1022570-eecb-49ab-8429-ca01bc3e7116
Agama.circular_velocity(Agama.Potential(Φ_nfw), 50) * V2KMS

# ╔═╡ 3de23e3a-b613-4330-9a08-5544dbbb911d
Φ_v24_mw_halo = Agama.Potential(file = joinpath(potential_dir, "vasiliev24/L3M11/potential_mw_halo.ini"))

# ╔═╡ 023ae372-9d5d-48c0-9f51-54e014e10223
Φ_v24_mw = Agama.Potential(file = joinpath(potential_dir, "vasiliev24/L3M11/potential_mw_init.ini"))

# ╔═╡ 8cb923e5-e7f0-448e-ab4a-41d6d1328175
Φ_v24_mw_stars = Agama.Potential(file = joinpath(potential_dir, "vasiliev24/L3M11/potential_stars.ini"))

# ╔═╡ 8b8b7053-2ce2-4fcb-8964-8e5795c35cd1
readdir("/cosma/home/durham/dc-boye1/data/dwarfs/agama/potentials/vasiliev24/L3M11")

# ╔═╡ c4917e59-1502-4a06-86eb-a31d008739a0
Φ_v24_lmc = Agama.Potential(file = joinpath(potential_dir, "vasiliev24/L3M11/potential_lmc_init.ini"))

# ╔═╡ f6bec4d8-7a6a-4b2f-be11-1658a3e39cdd
@savefig "v_circ_lmc" let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = L"$R$ / kpc",
		ylabel = L"$v_\textrm{circ}$ / km\,s$^{-1}$",
		limits=(0, 100, nothing, nothing),
	)
	lw = theme(:linewidth)

	plot_v_circ_agama!(Φ_v24_mw._py, label="total", linestyle=:solid, color=:black, linewidth=lw)

	plot_v_circ_agama!(Φ_v24_mw_stars._py, label="stars", linestyle=:dash, color=COLORS[1], linewidth=lw)
	plot_v_circ_agama!(Φ_v24_mw_halo._py, label="halo", linestyle=:dot,  color=COLORS[3], linewidth=lw)

	plot_v_circ_agama!(Φ_v24_lmc._py, label="LMC", linestyle=:dashdot,  color=COLORS[4], linewidth=lw)

	Legend(fig[1,2], ax, patchsize=(4/3*theme(:fontsize)[], theme(:fontsize)[]), margin=zeros(4))

	resize_to_layout!(fig)
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
# ╠═a2e7a222-cca9-47a3-901e-69895193a1db
# ╠═cf74dd35-70f3-4a88-b1c0-580d0ffe69aa
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═d1022570-eecb-49ab-8429-ca01bc3e7116
# ╠═3de23e3a-b613-4330-9a08-5544dbbb911d
# ╠═023ae372-9d5d-48c0-9f51-54e014e10223
# ╠═8cb923e5-e7f0-448e-ab4a-41d6d1328175
# ╠═8b8b7053-2ce2-4fcb-8964-8e5795c35cd1
# ╠═c4917e59-1502-4a06-86eb-a31d008739a0
# ╠═f6bec4d8-7a6a-4b2f-be11-1658a3e39cdd
