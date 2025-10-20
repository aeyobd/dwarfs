### A Pluto.jl notebook ###
# v0.20.19

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

# ╔═╡ 6161cd69-594b-4eb9-83d9-12a682db7c88
function get_scalars(modelname)
	galaxy, model, stars = modelname
	df = read_fits(joinpath(modeldir(joinpath(galaxy, model)), "profiles_scalars.fits")
	)
	idx_f = TOML.parsefile(joinpath(modeldir(joinpath(galaxy, model)), "orbital_properties.toml"))["idx_f"]

	if idx_f < size(df, 1)
		df.time .-= df.time[idx_f]
	else
		@warn "idx_f > number calculated profiles for $galaxy/$model"
		df.time .-= df.time[end]
	end
	df
end

# ╔═╡ f49c9325-4319-4435-ac33-a8980e1c4f7f
function get_profiles(modelname)
	mass_profs = LilGuys.read_ordered_structs(
		joinpath(modeldir(modelname), "profiles.hdf5"),
		LilGuys.MassProfile
	)
	return mass_profs[1], mass_profs[end]
end

# ╔═╡ 7d8f9954-0bbc-4fc0-ba21-aad54b3aaee4
function compare_boundmass(modelnames; relative=false, legend_position=:rb)
	
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "simulation time / Gyr",
		ylabel = L"$\textrm{v}_\textrm{max}$ / km\,s$^{-1}$",
	)

	
	for (label, modelname) in modelnames
		df = get_scalars(modelname)

		y = log10.(df.v_circ_max * V2KMS)
		if relative
			y .-=  log10(df.v_circ_max[1])
		end
		
		lines!(df.time * T2GYR, y, label=label)
	end




	ax2 = Axis(fig[1,2], 
		xlabel = L"$r_\textrm{max}$ / kpc",
		ylabel = L"$\textrm{v}_\textrm{max}$ / km\,s$^{-1}$",
	)

	
	for (label, modelname) in modelnames
		df = get_scalars(modelname)

		y = log10.(df.v_circ_max * V2KMS)
		if relative
			y .-=  log10(df.v_circ_max[1])
		end
		
		lines!(df.r_circ_max, y, label=label)
	end
	axislegend(position=legend_position, patchsize=(24, 12), backgroundcolor=(:white, 0.8))


	linkyaxes!(ax, ax2)
	hideydecorations!(ax2, ticks=false, minorticks=false)

	rowsize!(fig.layout, 1, Aspect(1, 3/3))
	resize_to_layout!()
	fig
end

# ╔═╡ 93ba7c53-81c3-40e8-9171-626c04c834ff
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ 483b7a71-3e64-41d8-b77d-7f26d5bd33da


# ╔═╡ faff05dd-31b7-47bf-b6da-32d94ca80514
@savefig "scl_mw_halo_boundmass" compare_boundmass(OrderedDict(
	"fiducial" => modelnames["scl_smallperi"],
	"heavier halo" => modelnames["heavier"],
	"heavier, new orbit" => modelnames["heavier_new_orbit"],
	"lighter halo" => modelnames["lighter"],
),
				 )

# ╔═╡ f99c5233-a05c-46e9-b7ac-6f8076377419
@savefig "scl_mw_structure_boundmass" compare_boundmass(OrderedDict(
	"heavy NFW" => modelnames["heavier"],
	"cored" => modelnames["cored"],
	"anisotropic" =>  modelnames["anisotropic"],
	"oblate" => modelnames["oblate"],
),
				 )

# ╔═╡ ffad6fd8-a767-459e-aa3d-e9fd6741a47f
@savefig "scl_orbits_boundmass"  compare_boundmass(OrderedDict(
	"fiducial" => modelnames["scl_smallperi"],
	"mean orbit" => modelnames["scl_mean"],
	"MW-impact" => modelnames["mw_impact"],

))

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═c1b20bd0-bc18-4c53-b687-5f6b792fdcd0
# ╠═1482481d-a5f3-48c2-a4d2-1353afe7fd72
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═8081e889-c1d4-4e50-b8e1-91aa8e82adf1
# ╠═6161cd69-594b-4eb9-83d9-12a682db7c88
# ╠═f49c9325-4319-4435-ac33-a8980e1c4f7f
# ╠═7d8f9954-0bbc-4fc0-ba21-aad54b3aaee4
# ╠═93ba7c53-81c3-40e8-9171-626c04c834ff
# ╠═483b7a71-3e64-41d8-b77d-7f26d5bd33da
# ╠═faff05dd-31b7-47bf-b6da-32d94ca80514
# ╠═f99c5233-a05c-46e9-b7ac-6f8076377419
# ╠═ffad6fd8-a767-459e-aa3d-e9fd6741a47f
