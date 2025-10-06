### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

	using PyFITS
end

# ╔═╡ 1482481d-a5f3-48c2-a4d2-1353afe7fd72
using OrderedCollections

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ c1b20bd0-bc18-4c53-b687-5f6b792fdcd0
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ 8081e889-c1d4-4e50-b8e1-91aa8e82adf1
function modeldir(modelname)
	joinpath(ENV["DWARFS_ROOT"], "analysis",modelname,)
end

# ╔═╡ f49c9325-4319-4435-ac33-a8980e1c4f7f
function get_profiles(modelname)
	mass_profs = LilGuys.read_ordered_structs(
		joinpath(modeldir(modelname), "profiles_densities.hdf5"),
		LilGuys.DensityProfile
	)
	idx_f = TOML.parsefile(joinpath(modeldir(modelname), "orbital_properties.toml"))["idx_f"]
						   
	return mass_profs[1].second, mass_profs[idx_f].second
end

# ╔═╡ b50a1346-19e1-49a7-b723-c7fb001c9d20
get_profiles("sculptor/1e6_v49_r7.1_c0.1/orbit_smallperi")

# ╔═╡ cb38048e-99d9-431c-ae90-de8817e02ff7
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log radius / kpc", ylabel = "log density",
			  limits=(-1.5, 2, -8, -1.0)
			 )





	profs = get_profiles("sculptor/1e6_new_v43_r7/orbit_smallperi")

	lines!(profs[1].log_r, LilGuys.log_densities(profs[1]), color=COLORS[1], linestyle=:dot, label="cuspy initial")
	lines!(profs[2].log_r, LilGuys.log_densities(profs[2]), color=COLORS[1], label="cuspy final")

	profs = get_profiles("sculptor/1e6_v49_r7.1_c0.1/orbit_smallperi")

	lines!(profs[1].log_r, LilGuys.log_densities(profs[1]), linestyle=:dot, color=COLORS[2], label="cored initial")
	lines!(profs[2].log_r, LilGuys.log_densities(profs[2]), color=COLORS[2], label="cored final")

	axislegend(position=:lb)


	@savefig "cored_density_i_f"

	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═c1b20bd0-bc18-4c53-b687-5f6b792fdcd0
# ╠═1482481d-a5f3-48c2-a4d2-1353afe7fd72
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═8081e889-c1d4-4e50-b8e1-91aa8e82adf1
# ╠═f49c9325-4319-4435-ac33-a8980e1c4f7f
# ╠═b50a1346-19e1-49a7-b723-c7fb001c9d20
# ╠═cb38048e-99d9-431c-ae90-de8817e02ff7
