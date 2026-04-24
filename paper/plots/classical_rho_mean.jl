### A Pluto.jl notebook ###
# v0.20.24

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
# Calculations
"""

# ╔═╡ 98d6c4aa-002d-44b1-9a0b-896f5c068a01
md"""
## Pericentres
"""

# ╔═╡ 1a41df8c-40b9-41d6-b33d-a09f050093c0
function calc_peri(pot, galaxyname, units=Agama.AgamaUnits(), time_max=-1060)
	ic = get_obs_props(galaxyname) |> ICRS
	orbit = LilGuys.agama_orbit(pot, ic, timerange=(0, time_max), agama_units=units, N=10_000)
	return LilGuys.pericenter(orbit)
end

# ╔═╡ f3f2f788-ce44-4390-b7c5-fdcd862aff49
pot_mw = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))

# ╔═╡ 1058cf63-0fa5-4b72-ac74-ca59a173715c
pot_lmc_evolving = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/vasiliev24/L3M11/potential.ini"))

# ╔═╡ ce0a518e-1cc9-4b3c-a9a6-b465830ddcea
α_3d_2d = LilGuys.r_h(LilGuys.Exp2D()) / LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 08b220e8-cfe0-478f-85cd-d73b85f47204
md"""
### Mean densities
"""

# ╔═╡ 1aea04a4-a1bf-4e04-b73e-2fb9afd052f5
function get_MW_mean_density(pot=pot_mw, radii=LinRange(7, 150, 10000); units=Agama.AgamaUnits())
	Ms = Agama.enclosed_mass(pot, radii, units)
	return radii, Ms ./ (4π/3 * radii .^3) * M2MSUN
end

# ╔═╡ dfada36d-3769-462b-a4df-a48069df3b45
function sample_key(dict, key, N)
	if key * "_err" in keys(dict)
		δx = dict[key * "_err"] 
	else
		δx = max(dict[key * "_em"], dict[key * "_ep"]) 
	end

	return dict[key] .+ δx * randn(N)
end

# ╔═╡ f5ef70ba-4471-4e40-88ee-fa8c200fc2ca
function get_rho_mean_h_mc(galaxyname, N=10_000)	
	props = get_obs_props(galaxyname)
	
	σv = sample_key(props, "sigma_v", N) / V2KMS
	R_h = sample_key(props, "R_h", N)
	dm = sample_key(props, "distance_modulus", N)
	dist = LilGuys.dm2kpc.(dm)
	
	r_h = LilGuys.arcmin2kpc.(R_h, dist) * α_3d_2d
	M_h = @. 3*r_h ./ LilGuys.G * σv^2 

	ρ = @. M_h / (4π/3 * r_h^3) * LilGuys.M2MSUN

	m = median(ρ)
	l, h = quantile(ρ, [0.16, 0.84])
	

	return LilGuys.Measurement(m, m-l, h-m)
end

# ╔═╡ 19d388e6-bd12-42c9-ac59-29c7a47d3f26
md"""
# Plotting
"""

# ╔═╡ ec7d88d3-7503-4c18-9088-04319c6f2e98
smallfontsize=0.8 * theme(:fontsize)[]

# ╔═╡ 6c106a00-445b-4e0f-9168-66a9aeb35767
smalllinewidth=theme(:linewidth)[]/2

# ╔═╡ 7e063f49-0b32-4aff-9a0d-4cac8251fbb3
markers = theme(:palette)[:marker][]

# ╔═╡ 0c08a41a-b21f-4168-82fe-b0ff3da6111e
classical_labels = Dict(
	"sagittarius" => "Sagittarius",
	"fornax" => "Fornax",
	"leo1" => "Leo I", 
	"sculptor" => "Sculptor", 
	"antlia2" => "Antlia II",
	"leo2" => "Leo II",
	"carina" => "Carina", 
	"draco" => "Draco",
	"ursa_minor" => "Ursa Minor", 
	"canes_venatici1" => "Canes Venatici I",
	"sextans1" => "Sextans I",
	"crater2" => "Crater II"
)

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

# ╔═╡ 98c02ad0-ee51-477b-94d5-6c1525b5960e
pericentres = [calc_peri(pot_mw, galaxyname) for galaxyname in classicals]

# ╔═╡ 253342ea-e551-4cc3-bc4f-e9294949e933
pericentres_v24 = [calc_peri(pot_lmc_evolving, galaxyname, Agama.VASILIEV_UNITS) for galaxyname in classicals]

# ╔═╡ 0a237127-8ebd-4400-b6b5-57f4cee8febc
ρ_bar = [get_rho_mean_h_mc(galaxyname) for galaxyname in classicals]

# ╔═╡ 6deb29ed-e2f0-4d54-ae61-f8b8bb0d700c
aligns = Dict(
	"carina" => (:left, :center),
	"sculptor" => (:left, :center),
	"sagittarius" => (:left, :center),
	"leo2" => (:left, :top),
	"draco" => (:left, :bottom),
	"fornax" => (:left, :bottom)
)

# ╔═╡ 34d39ebc-cc90-427b-a492-a556ec0455be
colors = [
	1
	3
	4
	1
	5
	2
	1
	3
	5
	4
	2
	3
	]

# ╔═╡ d9606d05-c279-43c6-8662-8eced5bdfc64
offsets = Dict(
	"leo1" => (-15, 0),
	"draco" => (2, 2),
	"leo2" => (0, -6),
	"carina" => (3, 0),
	"sculptor" => (3, 0),
	"sagittarius" => (7, 0),
	"crater2" => (-3, -8),
	"antlia2" => (-2, 9),
	"fornax" => (-2, 2),

)

# ╔═╡ d11de3cb-bbfb-446a-b384-06c9f254cd9b
light_grey = RGBf(0.9, 0.9, 0.9)

# ╔═╡ 7875cdf7-5cdb-4ef4-8a9c-6f4d057705d0
function plot_mw_disruption_region!()
	x, y = get_MW_mean_density()
	lines!(x, 3y, color="black", label=L"3\bar{\rho}_\textrm{MW}")
	text!(15, 2e7, text="likely tidal disruption", fontsize=smallfontsize, rotation=-0.6)
	band!(x, y*1e-3, 3y, color=light_grey)
end

# ╔═╡ 697e49e7-0d4a-4109-95a1-e7602d7cfd70
let
	fig = Figure()
	ax = Axis(fig[1, 1],
			 xlabel = "pericentre / kpc",
			 ylabel = L"$\bar{\rho}_h$ / $\textrm{M}_\odot\,\textrm{kpc}^{-3}$", yscale=log10, 
			  xscale = log10,
			  xticks = [10, 30, 100, 130],
			  xminorticks = 10:10:150,
			  yticks=[1e6, 1e7, 1e8],
			 limits =(10, 150, 3e5, 3e8)
			 )

	plot_mw_disruption_region!()

	for i in eachindex(classicals)
		x1 = (pericentres[i])
		x2 = (pericentres_v24[i])
		y = middle(ρ_bar[i])

	
		color = COLORS[colors[i]]
		marker = :circle
		fillcolor = y < 1e7 ? light_grey : "white"

		lines!([x1, x2], 
			   fill(y, 2), linewidth=smalllinewidth, color=color
			  )
		
		errorscatter!([x1, x2], [y, y], yerror=fill(error_interval(ρ_bar[i]), 2),
				 color=color, marker=marker, strokewidth=0)

		# fill 
		scatter!([x2], [y], 
				 color=fillcolor, strokecolor=color, strokewidth=smalllinewidth, marker=marker)

		# label
		offset = get(offsets, classicals[i], (-5,0))
		align = get(aligns, classicals[i], (:right, :center) )
		text!((x1), y, text=classical_labels[classicals[i]], color=color, align=align, fontsize=smallfontsize * 0.8, offset=offset)

	end
	

	# create legend
	scatter!([NaN], [NaN], color=:black, marker=:circle, label="MW-only")
	scatter!([NaN], [NaN], color=:white, strokecolor=:black, strokewidth=smalllinewidth, marker=:circle, label="MW+LMC")
	axislegend(position=:rb)

	@savefig "rho_mean_pericentre"
	fig
end

# ╔═╡ f0b7f0ae-e19c-4d97-8175-9365ba597b1a
md"""
# Sanity checks, etc.
"""

# ╔═╡ 475816ca-baa6-4b51-a4b1-aeeec13ff105
function get_rho_mean_h(galaxyname)
	props = get_obs_props(galaxyname) |> LilGuys.collapse_errors
	σv = props["sigma_v"] / LilGuys.V2KMS
	r_h = LilGuys.arcmin2kpc(props["R_h"], props["distance"]) * α_3d_2d
	M_h_over_r_h = 3 * σv^2 / LilGuys.G  # times r_h normally

	return M_h_over_r_h / (4π/3 * r_h^2) * LilGuys.M2MSUN
end

# ╔═╡ 7d0ea3e7-100d-4a39-9a2a-c44f39279bb2
let
	ic = get_obs_props("sculptor") |> ICRS
	println(ic)
	gc = LilGuys.transform(Galactocentric, ic)
	orbit = LilGuys.agama_orbit(pot_lmc_evolving, ic, timerange=(0, -1060), agama_units=Agama.VASILIEV_UNITS, N=10_000)


	lines(orbit.positions[2, :], orbit.positions[3, :])
end

# ╔═╡ 65d80746-7e16-41ab-8db5-feb7b3f23d8c
let

	fig = Figure()
	ax = Axis(fig[1,1])


	for gal in classicals
		ic = get_obs_props(gal) |> ICRS
		orbit = LilGuys.agama_orbit(pot_mw, ic, timerange=(0, -1200), N=10_000)
	
	
		lines!(orbit.times, radii(orbit))
	end
	fig
end

# ╔═╡ 68f8ff23-f5a6-4dfd-8446-085d32fa3e62
let

	fig = Figure()
	ax = Axis(fig[1,1])


	for gal in classicals
		ic = get_obs_props(gal) |> ICRS
		orbit = LilGuys.agama_orbit(pot_lmc_evolving, ic, timerange=(0, -1200), N=10_000, agama_units=Agama.VASILIEV_UNITS)
	
	
		lines!(orbit.times, radii(orbit))
	end
	fig
end

# ╔═╡ c8a8d5d8-1d50-4bac-8f4f-30d064ef158b
function get_rel_err(dict, key)
	if key * "_err" in keys(dict)
		return dict[key * "_err"] / dict[key]
	else
		return max(dict[key * "_em"], dict[key * "_ep"]) / dict[key]
	end
end

# ╔═╡ 72264878-8f46-43b7-8795-9986fd277b1e
⊕(a, b) = sqrt(a^2 + b^2)

# ╔═╡ 354a9bba-93c0-4c02-8cdc-9fb20063780a
function get_rho_mean_h_uncert(galaxyname)
	props = get_obs_props(galaxyname)
	σv = props["sigma_v"] / LilGuys.V2KMS
	r_h = LilGuys.arcmin2kpc(props["R_h"], props["distance"]) * α_3d_2d
	M_h_over_r_h = 3 / LilGuys.G * σv^2 # times r_h normally

	R_h_err = get_rel_err(props, "R_h")
	d_err = get_rel_err(props, "distance")
	σv_err = get_rel_err(props, "sigma_v")
	rel_err = 2R_h_err ⊕ 2d_err ⊕ 2σv_err

	val =  M_h_over_r_h / (4π/3 * r_h^2) * LilGuys.M2MSUN

	return Measurement(val, val * (rel_err))
end

# ╔═╡ 0dbc4770-fd1d-4a57-9825-3576900dd7a3
function get_peri(galaxyname, potname="EP2020", colname="pericentre")
	modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, potname)

	filename = joinpath(modeldir, "orbital_properties.fits")
	if !isfile(filename)
		return Measurement(NaN, NaN)
	end
	props = read_fits(filename)

	m = median((props[!, colname]))
	l, h = quantile(props[!, colname], [0.16, 0.84])
	LilGuys.Measurement(m, m-l, h-m)
end

# ╔═╡ ef6af06e-492d-44dd-b68e-cf71b3b07347
pericentres2 = [get_peri(galaxyname) for galaxyname in classicals]

# ╔═╡ 3a4ee436-3a33-468e-9ede-d3ed45356432
for name in classicals
	println(name)
	println(get_rho_mean_h(name))
	println(get_rho_mean_h_uncert(name))
	println(get_rho_mean_h_mc(name))
end

# ╔═╡ 7e79a3b3-b231-47bb-a984-f35049174f66
pericentres_v24_2 = [get_peri(galaxyname, "vasiliev24_L3M11") for galaxyname in classicals]

# ╔═╡ 122dcae8-39eb-4fa7-863a-002ad44faef1
for i in eachindex(classicals)
	println(classicals[i], "\t", pericentres[i], "\t", pericentres2[i])
end

# ╔═╡ a0095620-9286-472e-8236-ff8f9b8ad008
for i in eachindex(classicals)
	println(classicals[i], "\t", pericentres_v24[i], "\t", pericentres_v24_2[i])
end

# ╔═╡ 7533ef4c-38f1-4c12-9953-3a7d4c0cefe1
scatter(pericentres, pericentres2)

# ╔═╡ 08d717ef-ba69-41dc-b5d2-f1273a8f30d3
scatter(pericentres_v24, middle.(pericentres_v24_2))

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
# ╟─98d6c4aa-002d-44b1-9a0b-896f5c068a01
# ╠═1a41df8c-40b9-41d6-b33d-a09f050093c0
# ╠═f3f2f788-ce44-4390-b7c5-fdcd862aff49
# ╠═1058cf63-0fa5-4b72-ac74-ca59a173715c
# ╠═ce0a518e-1cc9-4b3c-a9a6-b465830ddcea
# ╠═98c02ad0-ee51-477b-94d5-6c1525b5960e
# ╠═253342ea-e551-4cc3-bc4f-e9294949e933
# ╟─08b220e8-cfe0-478f-85cd-d73b85f47204
# ╠═1aea04a4-a1bf-4e04-b73e-2fb9afd052f5
# ╠═dfada36d-3769-462b-a4df-a48069df3b45
# ╠═f5ef70ba-4471-4e40-88ee-fa8c200fc2ca
# ╠═0a237127-8ebd-4400-b6b5-57f4cee8febc
# ╟─19d388e6-bd12-42c9-ac59-29c7a47d3f26
# ╠═ec7d88d3-7503-4c18-9088-04319c6f2e98
# ╠═6c106a00-445b-4e0f-9168-66a9aeb35767
# ╠═7e063f49-0b32-4aff-9a0d-4cac8251fbb3
# ╠═0c08a41a-b21f-4168-82fe-b0ff3da6111e
# ╠═aadbdcc6-1695-431d-92de-c08e86e1c1f0
# ╠═6deb29ed-e2f0-4d54-ae61-f8b8bb0d700c
# ╠═34d39ebc-cc90-427b-a492-a556ec0455be
# ╠═d9606d05-c279-43c6-8662-8eced5bdfc64
# ╠═d11de3cb-bbfb-446a-b384-06c9f254cd9b
# ╠═7875cdf7-5cdb-4ef4-8a9c-6f4d057705d0
# ╠═697e49e7-0d4a-4109-95a1-e7602d7cfd70
# ╟─f0b7f0ae-e19c-4d97-8175-9365ba597b1a
# ╠═ef6af06e-492d-44dd-b68e-cf71b3b07347
# ╠═475816ca-baa6-4b51-a4b1-aeeec13ff105
# ╠═7d0ea3e7-100d-4a39-9a2a-c44f39279bb2
# ╠═65d80746-7e16-41ab-8db5-feb7b3f23d8c
# ╠═68f8ff23-f5a6-4dfd-8446-085d32fa3e62
# ╠═c8a8d5d8-1d50-4bac-8f4f-30d064ef158b
# ╠═72264878-8f46-43b7-8795-9986fd277b1e
# ╠═354a9bba-93c0-4c02-8cdc-9fb20063780a
# ╠═0dbc4770-fd1d-4a57-9825-3576900dd7a3
# ╠═3a4ee436-3a33-468e-9ede-d3ed45356432
# ╠═7e79a3b3-b231-47bb-a984-f35049174f66
# ╠═122dcae8-39eb-4fa7-863a-002ad44faef1
# ╠═a0095620-9286-472e-8236-ff8f9b8ad008
# ╠═7533ef4c-38f1-4c12-9953-3a7d4c0cefe1
# ╠═08d717ef-ba69-41dc-b5d2-f1273a8f30d3
