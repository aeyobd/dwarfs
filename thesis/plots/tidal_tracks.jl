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

	using PyFITS
end

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ 8081e889-c1d4-4e50-b8e1-91aa8e82adf1
function modeldir(modelname)
	joinpath(ENV["DWARFS_ROOT"], "analysis", "sculptor",modelname,)
end

# ╔═╡ 6161cd69-594b-4eb9-83d9-12a682db7c88
function get_scalars(modelname)
	read_fits(joinpath(modeldir(modelname), "profiles_scalars.fits")
	)
end

# ╔═╡ f49c9325-4319-4435-ac33-a8980e1c4f7f
function get_profiles(modelname)
	mass_profs = LilGuys.read_ordered_structs(
		joinpath(modeldir(modelname), "profiles.hdf5"),
		LilGuys.MassProfile
	)
	return mass_profs[1], mass_profs[end]
end

# ╔═╡ d8ef459b-bd41-4f71-bdf8-d7795f16cd58
df = get_scalars("1e7_new_v31_r3.2/orbit_smallperi")

# ╔═╡ 415f249c-de75-454f-8ab8-be3bec7e4222
profs = get_profiles("1e7_new_v31_r3.2/orbit_smallperi")

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "scl_tidal_track" let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = L"r_\textrm{circ}",
		ylabel = L"\textrm{v}_\textrm{circ}",
		xscale=log10,
		yscale=log10,
		xticks = [0.1, 1, 10],
		title = "Sculptor Dark Matter: smallperi orbit"
	)


	x, y = LilGuys.EN21_tidal_track(df.r_circ_max[1], df.v_circ_max[1], x_min=0.3)

	
	scatter!(df.r_circ_max, df.v_circ_max * V2KMS, label="max", markersize=theme(:markersize)[]/2)

	lines!(radii(profs[1].second), LilGuys.circular_velocity(profs[1].second) * V2KMS, label=L"initial ($t=-9.2$\,Gyr)")
	lines!(radii(profs[end].second), LilGuys.circular_velocity(profs[end].second) * V2KMS, label=L"final ($t=0.0$\,Gyr)")
	lines!(x, y * V2KMS, color=:black,linestyle=:dot, label="EN21 tidal track")

	axislegend(position=:lt, patchsize=(24, 12), backgroundcolor=(:white, 0.5))
	xlims!(0.1, 20)
	ylims!(10, 33)

	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═8081e889-c1d4-4e50-b8e1-91aa8e82adf1
# ╠═6161cd69-594b-4eb9-83d9-12a682db7c88
# ╠═f49c9325-4319-4435-ac33-a8980e1c4f7f
# ╠═d8ef459b-bd41-4f71-bdf8-d7795f16cd58
# ╠═415f249c-de75-454f-8ab8-be3bec7e4222
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
