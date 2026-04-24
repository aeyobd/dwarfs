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

	using PyFITS

	import TOML
	using OrderedCollections
end

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 3400e874-28f2-496d-9493-1624de9e7fdc
theme(:Lines)[:cycle][] = Cycle([[:color] => :color, [:linestyle] => :linestyle], true)

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ 15ee589c-cb4d-4796-b296-544b09d2b59f
md"""
# Functions
"""

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

	@assert idx_f < size(df, 1)
	df.time .-= df.time[idx_f]

	df
end

# ╔═╡ f06150d9-2c43-4414-98de-d657c873e4cc
function compare_boundmass(modelnames; legend_position=:rb)
	
	fig = Figure()

	ax_boundmass_time = Axis(fig[1,1], 
		xlabel = "simulation time / Gyr",
		ylabel = L"log bound mass / $10^{10}\,\textrm{M}_\odot$",
	)
	
	for (label, modelname) in modelnames
		df = get_scalars(modelname)

		y = df.bound_mass .|> log10	
		lines!(df.time * T2GYR, y, label=label, joinstyle=:miter)
	end


	ax_rmax_vmax = Axis(fig[1,2], 
		xlabel = L"$r_\textrm{max}$ / kpc",
		xscale = log10,
		xticks=[1; 2:1:8],
		xminorticks=[1:0.1:2; 2:0.5:8],
		yticks=10:5:60,
		yminorticks=10:1:60,
		ylabel = L"$\textrm{v}_\textrm{max}$ / km\,s$^{-1}$",
		yscale = log10
	)

	
	for (label, modelname) in modelnames
		df = get_scalars(modelname)
		y = (df.v_circ_max * V2KMS)
		
		lines!(df.r_circ_max, y, label=label, joinstyle=:miter)
	end

	
	Legend(fig[2, :], ax_boundmass_time, tellwidth=false, tellheight=true, nbanks=2, patchsize=(20, 5))

	resize_to_layout!()
	fig
end

# ╔═╡ 93ba7c53-81c3-40e8-9171-626c04c834ff
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ f99c5233-a05c-46e9-b7ac-6f8076377419
@savefig "extra_tidal_tracks" compare_boundmass(
	OrderedDict(
		"fiducial" => modelnames["scl_smallperi"],
		"heavy NFW" => modelnames["heavier"],
		"cored" => modelnames["cored"],
		"anisotropic" =>  modelnames["anisotropic"],
	),
)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═3400e874-28f2-496d-9493-1624de9e7fdc
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╟─15ee589c-cb4d-4796-b296-544b09d2b59f
# ╠═8081e889-c1d4-4e50-b8e1-91aa8e82adf1
# ╠═6161cd69-594b-4eb9-83d9-12a682db7c88
# ╠═f06150d9-2c43-4414-98de-d657c873e4cc
# ╠═93ba7c53-81c3-40e8-9171-626c04c834ff
# ╠═f99c5233-a05c-46e9-b7ac-6f8076377419
