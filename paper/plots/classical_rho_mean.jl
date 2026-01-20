### A Pluto.jl notebook ###
# v0.20.21

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

# ╔═╡ df4f20bf-97d9-408a-859e-e4070edcd0ef
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 3abb0de9-81b7-4fe8-b457-92ae4b2e78e3
import CSV

# ╔═╡ 7a199993-0622-49cc-867b-8df8504447fe
import DataFrames: DataFrame

# ╔═╡ 11cc55a3-d166-4bf3-b147-1c789d690f90
import Agama

# ╔═╡ b74a4535-1e12-422e-a524-135db4d8a9f7
import TOML

# ╔═╡ 191b6df6-0fa4-4393-b69d-2aefeb3f9373
import StatsBase: quantile, median

# ╔═╡ 20c338e3-3d10-41f4-b8ee-d0eda4e755bd
CairoMakie.activate!(type=:png)

# ╔═╡ a4c87a6e-976d-4af6-868e-09cb85e3d424
module Utils
	include("utils.jl")
end

# ╔═╡ 50839b67-9514-47b9-add2-0c84d05f12da
function get_obs_props(galaxyname)
	return TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))
end

# ╔═╡ 43dad299-d8d2-4146-8320-c90a62f3a3f0
md"""
# Data loading
"""

# ╔═╡ 0dbc4770-fd1d-4a57-9825-3576900dd7a3
function get_peri(galaxyname, potname="EP2020", colname="pericentre")
	modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, potname)
	props = read_fits(joinpath(modeldir, "orbital_properties.fits"))

	m = median((props[!, colname]))
	l, h = quantile(props[!, colname], [0.16, 0.84])
	LilGuys.Measurement(m, m-l, h-m)
end

# ╔═╡ ce0a518e-1cc9-4b3c-a9a6-b465830ddcea
α_3d_2d = LilGuys.r_h(LilGuys.Exp2D()) / LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 475816ca-baa6-4b51-a4b1-aeeec13ff105
function get_rho_mean_h(galaxyname)
	props = get_obs_props(galaxyname) |> LilGuys.collapse_errors
	σv = props["sigma_v"] / LilGuys.V2KMS
	r_h = LilGuys.arcmin2kpc(props["R_h"], props["distance"]) * α_3d_2d
	M_h_over_r_h = 3 / LilGuys.G * σv^2 # times r_h normally

	return M_h_over_r_h / (4π/3 * r_h^2) * LilGuys.M2MSUN
end

# ╔═╡ aadbdcc6-1695-431d-92de-c08e86e1c1f0
classicals = [
	"sagittarius",
	"fornax",
	"leo1", 
	"sculptor", 
	"antlia2",
	"leo2",
	"carina", 
	"draco",
	"ursa_minor", 
	"canes_venatici1",
	"sextans1",
	"crater2"
]

# ╔═╡ 0c08a41a-b21f-4168-82fe-b0ff3da6111e
classical_labels = Dict(
	"sagittarius" => "Sagittarius",
	"fornax" => "Fornax",
	"leo1" => "Leo I", 
	"sculptor" => "Sculptor", 
	"antlia2" => "Antlia2",
	"leo2" => "Leo2",
	"carina" => "Carina", 
	"draco" => "Draco",
	"ursa_minor" => "Ursa Minor", 
	"canes_venatici1" => "Canes Venatici I",
	"sextans1" => "Sextans I",
	"crater2" => "Crater II"
)

# ╔═╡ ef6af06e-492d-44dd-b68e-cf71b3b07347
pericentres = [get_peri(galaxyname) for galaxyname in classicals]

# ╔═╡ 7b09fdb3-6edb-4eaf-994b-c2f508112c02
# pericentres_lmc = [get_peri(galaxyname, "vasiliev24_L3M11") for galaxyname in classicals]

# ╔═╡ 0a237127-8ebd-4400-b6b5-57f4cee8febc
ρ_bar = [get_rho_mean_h(galaxyname) for galaxyname in classicals]

# ╔═╡ f3f2f788-ce44-4390-b7c5-fdcd862aff49
pot_mw = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))

# ╔═╡ 1ded763d-8f25-44ab-be6f-935f6d26e8c4
pot_lmc = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/vasiliev24/L3M11/potential_lmc_init.ini"))

# ╔═╡ 99c6b883-f7ed-46e7-855b-23e234e1fae4
pot_mw_lmc = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/vasiliev24/L3M11/potential_mw_init.ini"))

# ╔═╡ 1aea04a4-a1bf-4e04-b73e-2fb9afd052f5
function get_MW_mean_density(radii=LinRange(7, 150, 10000))
	Ms = Agama.enclosed_mass(pot_mw, radii)
	return radii, Ms ./ (4π/3 * radii .^3) * M2MSUN
end

# ╔═╡ ec7d88d3-7503-4c18-9088-04319c6f2e98
smallfontsize=0.8 * theme(:fontsize)[]

# ╔═╡ 697e49e7-0d4a-4109-95a1-e7602d7cfd70
let
	fig = Figure()
	ax = Axis(fig[1, 1],
			 xlabel = "pericentre / kpc",
			 ylabel = L"$\bar{\rho}_h$ / $\textrm{M}_\odot\,\textrm{kpc}^{-3}$", yscale=log10, 
			  xscale = log10,
			  yticks=[1e6, 1e7, 1e8],
			 limits =(10, 150, 3e5, 3e8)
			 )


	x, y = get_MW_mean_density()
	lines!(x, 3y, color="black")

	x = middle.(pericentres)
	y = middle.(ρ_bar)
	xerr = error_interval.(pericentres)
	yerr = error_interval.(ρ_bar)
	scatter!(x, y)
	errorscatter!(x, y, xerror=xerr, yerror=yerr)
	annotation!(x, y, text=[classical_labels[galaxyname] for galaxyname in classicals], fontsize=smallfontsize)


	@savefig "rho_mean_pericentre"
	fig
end

# ╔═╡ 6c106a00-445b-4e0f-9168-66a9aeb35767
smalllinewidth=theme(:linewidth)[]/2

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═3abb0de9-81b7-4fe8-b457-92ae4b2e78e3
# ╠═7a199993-0622-49cc-867b-8df8504447fe
# ╠═11cc55a3-d166-4bf3-b147-1c789d690f90
# ╠═b74a4535-1e12-422e-a524-135db4d8a9f7
# ╠═df4f20bf-97d9-408a-859e-e4070edcd0ef
# ╠═191b6df6-0fa4-4393-b69d-2aefeb3f9373
# ╠═20c338e3-3d10-41f4-b8ee-d0eda4e755bd
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═a4c87a6e-976d-4af6-868e-09cb85e3d424
# ╠═50839b67-9514-47b9-add2-0c84d05f12da
# ╟─43dad299-d8d2-4146-8320-c90a62f3a3f0
# ╠═0dbc4770-fd1d-4a57-9825-3576900dd7a3
# ╠═475816ca-baa6-4b51-a4b1-aeeec13ff105
# ╠═ce0a518e-1cc9-4b3c-a9a6-b465830ddcea
# ╠═aadbdcc6-1695-431d-92de-c08e86e1c1f0
# ╠═0c08a41a-b21f-4168-82fe-b0ff3da6111e
# ╠═ef6af06e-492d-44dd-b68e-cf71b3b07347
# ╠═7b09fdb3-6edb-4eaf-994b-c2f508112c02
# ╠═0a237127-8ebd-4400-b6b5-57f4cee8febc
# ╠═f3f2f788-ce44-4390-b7c5-fdcd862aff49
# ╠═1ded763d-8f25-44ab-be6f-935f6d26e8c4
# ╠═99c6b883-f7ed-46e7-855b-23e234e1fae4
# ╠═1aea04a4-a1bf-4e04-b73e-2fb9afd052f5
# ╠═697e49e7-0d4a-4109-95a1-e7602d7cfd70
# ╠═ec7d88d3-7503-4c18-9088-04319c6f2e98
# ╠═6c106a00-445b-4e0f-9168-66a9aeb35767
