### A Pluto.jl notebook ###
# v0.20.15

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

# ╔═╡ c2c320b8-3b14-4e49-9048-92a546a6b275
using HDF5

# ╔═╡ a1606881-5138-4c90-b66f-34bcffd563eb
using OrderedCollections

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./style.jl"); scale_theme_element!(:linewidth, 1/2)

# ╔═╡ 45f1d142-cce5-4d81-830c-b61efb849e77
import TOML

# ╔═╡ 1dc687cf-f56d-4e28-92f5-f0218fd5ff43
galaxyname = "sculptor"

# ╔═╡ 432dc7d3-c390-48cb-9e51-d4d0ff21620c
Nmax = 100

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png, px_per_unit=1)

# ╔═╡ 2d94b40a-2a38-48d5-9b7e-6b493d504aec


# ╔═╡ cbc79255-9254-4e01-9353-c87e62620ec8
module OrbitUtils
	include(joinpath(ENV["DWARFS_ROOT"], "orbits/orbit_utils.jl"))
end

# ╔═╡ f7d1d155-eaab-4589-b5e7-5efac2b8fb5f
module Utils 
	include("./utils.jl")
end

# ╔═╡ afb27a46-7b9f-4d61-8d22-2bc2566723e1
import Agama

# ╔═╡ b9b99de2-1c01-4562-8295-75b770e232dc
md"""
# Utils
"""

# ╔═╡ bd1e330a-a7b7-4273-87ba-7100cb6d5af8
orbits_dir = joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname)

# ╔═╡ 14f862cc-8e53-4d10-9d73-579ed1d26cd9
function read_traj(modelname)
	structs = LilGuys.read_ordered_structs(joinpath(orbits_dir, modelname, "orbits.hdf5"), LilGuys.Orbit)

	filt = 1:min(Nmax, length(structs))
	idxs = first.(structs)[filt]
	orbits = last.(structs)[filt]
	return orbits
end

# ╔═╡ 70c6eab7-0f69-41b9-8c09-11ab0aa2de46
orbits_fiducial = read_traj("EP2020")

# ╔═╡ 93ca6a6b-5f41-40d1-8661-31d0b3c5aa58
orbits_w_lmc = read_traj("vasiliev24_L3M11")

# ╔═╡ 18b45a96-ddcd-4a9d-88e2-4dde209676cf
read_lmc_orbit(filename) = LilGuys.resample(OrbitUtils.get_lmc_orbit(joinpath(orbits_dir, filename)), orbits_w_lmc[1].times)

# ╔═╡ af5cc193-033b-48d2-b2f5-266f7c54c8d6
orbit_lmc = read_lmc_orbit("vasiliev24_L3M11")

# ╔═╡ bc7f29c3-0368-48e9-b74d-728ed072df44
orbits_no_lmc = read_traj("vasiliev24_M11")

# ╔═╡ 70d81f4a-8e01-45b1-a695-234d0249cca9
function plot_x_y_traj!(orbits; thin=1, x_direction=2, y_direction=3, alpha=0.1, color=:black, t_min = -2/T2GYR, kwargs...)

    for i in 1:thin:length(orbits)
		filt = orbits[i].times .> t_min 

        x = orbits[i].positions[x_direction, filt]
        y = orbits[i].positions[y_direction, filt]
        
        lines!(x, y; rasterize=true, alpha=alpha, color=color, kwargs...)
    
    end
end

# ╔═╡ 6b0210cc-353d-48d4-bd59-f81fdd422d55
function get_nbody_orbit(model)
	
	orbit = Orbit(joinpath(ENV["DWARFS_ROOT"], "analysis/$model/centres.hdf5"))

	idx_f = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis/$model/orbital_properties.toml"))["idx_f"]

	return orbit[1:idx_f]
end

# ╔═╡ 38f593ae-590b-4f41-908d-75b17e4c4a88
md"""
# References
"""

# ╔═╡ 9d78acb1-5439-4338-aaff-2860f7f30b77
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/observed_properties.toml"))

# ╔═╡ 27f32059-5af5-4b32-8b58-619b378acb19
icrs0 = ICRS(obs_props)

# ╔═╡ bfac7dd2-264c-448d-8755-725ebff4df95
gc0 = LilGuys.transform(LilGuys.Galactocentric, icrs0)

# ╔═╡ f4431bdf-b495-42ec-8aac-b7cf07f5746a
pos_0 = LilGuys.position(gc0)

# ╔═╡ 06212962-dc5c-4a61-a490-2f418a91ac4f
pot = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama", "potentials", "EP2020.ini"))

# ╔═╡ a4a79f8c-cfca-4345-bbe7-15b6e455c476
pot_stars = Agama.Potential(pot._py[0] + pot._py[1] + pot._py[2])

# ╔═╡ 79a75204-7727-431b-a370-1e6e2c2898ce
xz_iso = Utils.integrate_isodensity(pot_stars, x_direction=1, s_scale=3e-4, h_scale=1e-4)

# ╔═╡ d9ecd245-da11-4fcb-9044-0089abec59c0
yz_iso = Utils.integrate_isodensity(pot_stars, s_scale=3e-4)

# ╔═╡ efae8808-1768-40c1-a11d-253402cdfe43
function compare_x_y_traj(trajectories; r_max=300, kwargs...)
    fig = Figure()
  	ax2 = Axis(fig[1, 1], xlabel="x / kpc", ylabel="z / kpc",
        xgridvisible=false, ygridvisible=false, 
			  limits=(-1, 1, -1, 1) .* r_max,
			  #aspect= DataAspect()
    )

    for (i, (label, traj)) in enumerate(trajectories)
        plot_x_y_traj!(traj, x_direction=1, label=label=>(; alpha=0.5), color=COLORS[i]; kwargs...)
    end

	poly!(xz_iso..., color=COLORS[8])
	Utils.plot_sun!(x_direction=1)

	
    ax = Axis(fig[1, 2], xlabel="y / kpc", ylabel="z / kpc",
        xgridvisible=false, ygridvisible=false, 
			  limits=(-1, 1, -1, 1) .* r_max,
			  #aspect = DataAspect()
    )
    
    for (i, (label, traj)) in enumerate(trajectories)
        plot_x_y_traj!(traj, label=label=>(; alpha=0.5), color=COLORS[i]; kwargs...)
    end
	
	poly!(yz_iso..., color=COLORS[8])
	Utils.plot_sun!()
	
	hideydecorations!(ticks=false, minorticks=false)

	colsize!(fig.layout, 1, Aspect(1, 1.0))
	colsize!(fig.layout, 2, Aspect(1, 1.0))
	linkyaxes!(ax, ax2)

	resize_to_layout!(fig)
  
    fig
end

# ╔═╡ a1c44523-e31e-4b50-bb00-db80ab080004
md"""
# Plots
"""

# ╔═╡ 83da41b8-ffb5-4bef-8d1e-175094d69d41
orbit_nbody_smallperi = get_nbody_orbit("sculptor/1e7_new_v31_r3.2/orbit_smallperi")

# ╔═╡ f62f0c8b-bfa3-4bfe-8c02-b3d232f6fef6
lw_thick = @lift 2*$(theme(:linewidth))

# ╔═╡ c50bea4a-1f2e-4867-8f4c-8dc194e6ec4f
@savefig "scl_mc_orbits_xy_empty" let
	fig = compare_x_y_traj(Dict("fiducial" => []), t_min=-5/T2GYR, r_max=120, thin=2,)
	ax = fig.content[1]



	
	scatter!(pos_0[2], pos_0[3], markersize=20, color=COLORS[1], strokewidth=theme(:linewidth), strokecolor=:black)
	
	text!(pos_0[2], pos_0[3], text="today", align=(:left, :center), offset=(24, 0), fontsize=0.8*theme(:fontsize)[])


	Makie.current_axis!(fig.content[1])
	scatter!(pos_0[1], pos_0[3], markersize=20, color=COLORS[1], strokewidth=theme(:linewidth), strokecolor=:black)

	o = orbit_nbody_smallperi.positions[:, orbit_nbody_smallperi.times .- orbit_nbody_smallperi.times[end] .> -5/T2GYR]
	
	#lines!(o[1, :], o[3, :], color=:black, linewidth=lw_thick)

	Makie.current_axis!(fig.content[2])

	#lines!(o[2, :], o[3, :], color=:black, linewidth=lw_thick )


	resize_to_layout!(fig)
	fig
end

# ╔═╡ 702a673d-ea9d-44f7-8255-854aa5b27623
@savefig "scl_mc_orbits_xy" let
	fig = compare_x_y_traj(Dict("fiducial" => orbits_fiducial), t_min=-5/T2GYR, r_max=120, thin=2,)
	ax = fig.content[1]



	
	scatter!(pos_0[2], pos_0[3], markersize=20, color=COLORS[1], strokewidth=theme(:linewidth), strokecolor=:black)
	
	text!(pos_0[2], pos_0[3], text="today", align=(:left, :center), offset=(24, 0), fontsize=0.8*theme(:fontsize)[])


	Makie.current_axis!(fig.content[1])
	scatter!(pos_0[1], pos_0[3], markersize=20, color=COLORS[1], strokewidth=theme(:linewidth), strokecolor=:black)

	o = orbit_nbody_smallperi.positions[:, orbit_nbody_smallperi.times .- orbit_nbody_smallperi.times[end] .> -5/T2GYR]
	
	#lines!(o[1, :], o[3, :], color=:black, linewidth=lw_thick)

	Makie.current_axis!(fig.content[2])

	#lines!(o[2, :], o[3, :], color=:black, linewidth=lw_thick )


	resize_to_layout!(fig)
	fig
end

# ╔═╡ fb3f61f2-f425-40b4-9983-dfc7528df6ab
@savefig "scl_mc_orbits_xy_w_act" let
	fig = compare_x_y_traj(Dict("fiducial" => orbits_fiducial), t_min=-5/T2GYR, r_max=120, thin=2,)
	ax = fig.content[1]



	
	scatter!(pos_0[2], pos_0[3], markersize=20, color=COLORS[1], strokewidth=theme(:linewidth), strokecolor=:black)
	
	text!(pos_0[2], pos_0[3], text="today", align=(:left, :center), offset=(24, 0), fontsize=0.8*theme(:fontsize)[])


	Makie.current_axis!(fig.content[1])
	scatter!(pos_0[1], pos_0[3], markersize=20, color=COLORS[1], strokewidth=theme(:linewidth), strokecolor=:black)

	o = orbit_nbody_smallperi.positions[:, orbit_nbody_smallperi.times .- orbit_nbody_smallperi.times[end] .> -5/T2GYR]
	
	lines!(o[1, :], o[3, :], color=:black, linewidth=lw_thick)

	Makie.current_axis!(fig.content[2])

	lines!(o[2, :], o[3, :], color=:black, linewidth=lw_thick )


	resize_to_layout!(fig)
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═45f1d142-cce5-4d81-830c-b61efb849e77
# ╠═1dc687cf-f56d-4e28-92f5-f0218fd5ff43
# ╠═432dc7d3-c390-48cb-9e51-d4d0ff21620c
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═c2c320b8-3b14-4e49-9048-92a546a6b275
# ╠═a1606881-5138-4c90-b66f-34bcffd563eb
# ╠═2d94b40a-2a38-48d5-9b7e-6b493d504aec
# ╠═cbc79255-9254-4e01-9353-c87e62620ec8
# ╠═f7d1d155-eaab-4589-b5e7-5efac2b8fb5f
# ╠═afb27a46-7b9f-4d61-8d22-2bc2566723e1
# ╟─b9b99de2-1c01-4562-8295-75b770e232dc
# ╠═bd1e330a-a7b7-4273-87ba-7100cb6d5af8
# ╠═14f862cc-8e53-4d10-9d73-579ed1d26cd9
# ╠═18b45a96-ddcd-4a9d-88e2-4dde209676cf
# ╠═af5cc193-033b-48d2-b2f5-266f7c54c8d6
# ╠═70c6eab7-0f69-41b9-8c09-11ab0aa2de46
# ╠═93ca6a6b-5f41-40d1-8661-31d0b3c5aa58
# ╠═bc7f29c3-0368-48e9-b74d-728ed072df44
# ╠═efae8808-1768-40c1-a11d-253402cdfe43
# ╠═70d81f4a-8e01-45b1-a695-234d0249cca9
# ╠═6b0210cc-353d-48d4-bd59-f81fdd422d55
# ╟─38f593ae-590b-4f41-908d-75b17e4c4a88
# ╠═9d78acb1-5439-4338-aaff-2860f7f30b77
# ╠═27f32059-5af5-4b32-8b58-619b378acb19
# ╠═bfac7dd2-264c-448d-8755-725ebff4df95
# ╠═f4431bdf-b495-42ec-8aac-b7cf07f5746a
# ╠═06212962-dc5c-4a61-a490-2f418a91ac4f
# ╠═a4a79f8c-cfca-4345-bbe7-15b6e455c476
# ╠═79a75204-7727-431b-a370-1e6e2c2898ce
# ╠═d9ecd245-da11-4fcb-9044-0089abec59c0
# ╟─a1c44523-e31e-4b50-bb00-db80ab080004
# ╠═83da41b8-ffb5-4bef-8d1e-175094d69d41
# ╠═f62f0c8b-bfa3-4bfe-8c02-b3d232f6fef6
# ╠═c50bea4a-1f2e-4867-8f4c-8dc194e6ec4f
# ╠═702a673d-ea9d-44f7-8255-854aa5b27623
# ╠═fb3f61f2-f425-40b4-9983-dfc7528df6ab
