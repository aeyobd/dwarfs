### A Pluto.jl notebook ###
# v0.20.20

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
include("./paper_style.jl"); scale_theme_element!(:linewidth, 1/2)

# ╔═╡ 45f1d142-cce5-4d81-830c-b61efb849e77
import TOML

# ╔═╡ afb27a46-7b9f-4d61-8d22-2bc2566723e1
import Agama

# ╔═╡ cbc79255-9254-4e01-9353-c87e62620ec8
module OrbitUtils
	include(joinpath(ENV["DWARFS_ROOT"], "orbits/orbit_utils.jl"))
end

# ╔═╡ f7d1d155-eaab-4589-b5e7-5efac2b8fb5f
module Utils 
	include("./utils.jl")
end

# ╔═╡ 1dc687cf-f56d-4e28-92f5-f0218fd5ff43
galaxyname = "sculptor"

# ╔═╡ 432dc7d3-c390-48cb-9e51-d4d0ff21620c
Nmax = 100

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png, px_per_unit=4)

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

# ╔═╡ 70d81f4a-8e01-45b1-a695-234d0249cca9
function plot_x_y_traj!(orbits; thin=1, x_direction=2, y_direction=3, alpha=0.1, color=:black, t_min = -5/T2GYR, kwargs...)

    for i in 1:thin:length(orbits)
		filt = orbits[i].times .> t_min 

        x = orbits[i].positions[x_direction, filt]
        y = orbits[i].positions[y_direction, filt]
        
        lines!(x, y; alpha=alpha, color=color, kwargs...)
    
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

# ╔═╡ 8091ec68-c602-4d96-a418-d040aee84e41
pot_v24 = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama", "potentials/vasiliev24/L3M11/potential_stars.ini"))

# ╔═╡ 1ebe70d8-34b0-48d0-ad60-341223fb55fe
xz_iso_v24 = Utils.integrate_isodensity(pot_v24, s_scale=3e-4, x_direction=1)

# ╔═╡ d75f20c3-4a75-4a0b-b0fd-825c10113be8
yz_iso_v24 = Utils.integrate_isodensity(pot_v24, s_scale=3e-4)

# ╔═╡ 6d17d26a-9946-4656-be7b-13d84002c3fa
md"""
# Data loading
"""

# ╔═╡ 70c6eab7-0f69-41b9-8c09-11ab0aa2de46
orbits_fiducial = read_traj("EP2020")

# ╔═╡ 93ca6a6b-5f41-40d1-8661-31d0b3c5aa58
orbits_w_lmc = read_traj("vasiliev24_L3M11_9Gyr")

# ╔═╡ 18b45a96-ddcd-4a9d-88e2-4dde209676cf
read_lmc_orbit(filename) = LilGuys.resample(OrbitUtils.get_lmc_orbit(joinpath(orbits_dir, filename)), orbits_w_lmc[1].times)

# ╔═╡ af5cc193-033b-48d2-b2f5-266f7c54c8d6
orbit_lmc = read_lmc_orbit("vasiliev24_L3M11")

# ╔═╡ 4529b1cc-5046-4491-8932-4c3fbb1f1e48
pos_0_lmc = orbit_lmc.positions[:, 1]

# ╔═╡ bc7f29c3-0368-48e9-b74d-728ed072df44
orbits_no_lmc = read_traj("vasiliev24_M11")

# ╔═╡ 83da41b8-ffb5-4bef-8d1e-175094d69d41
orbit_nbody_lmc = get_nbody_orbit("sculptor/1e7_new_v25_r2.5/smallperilmc")

# ╔═╡ a1c44523-e31e-4b50-bb00-db80ab080004
md"""
# Plots
"""

# ╔═╡ f62f0c8b-bfa3-4bfe-8c02-b3d232f6fef6
lw_thick = @lift 2*$(theme(:linewidth))

# ╔═╡ ff89d7db-c954-48b8-b87a-446ccfb2d79b
trajectories = OrderedDict(

	"Sculptor: MW-only" => orbits_no_lmc,
	 "Sculptor: MW & LMC" => orbits_w_lmc,
)

# ╔═╡ eac841e7-5015-48cc-898b-f0232c8d0d72
orbits_w_lmc[1].times[end], orbits_no_lmc[1].times[end]

# ╔═╡ 92fb5b1f-5dc6-4766-b69d-f8644861d4cd
md"""
# Trajectories with LMC
"""

# ╔═╡ 0f3018a2-42e2-4140-8e9e-f094276b76f2
orbits_m_lmc = orbits_w_lmc .- [orbit_lmc]

# ╔═╡ 751ae520-59a2-471c-9247-90fba399cb42
Arya.update_fontsize!(12)

# ╔═╡ 6f1a1958-708a-4aac-935d-0e5a87091a98
function plot_traj_lmc!(; x_direction = 2, t_min=-5/T2GYR)
	times = orbit_lmc.times
	times = times[times .> t_min]
	plot_x_y_traj!([orbit_lmc], label="LMC", color=:black, linewidth=lw_thick, t_min=t_min, alpha=1, x_direction=x_direction)
end


# ╔═╡ f4bf5ad4-7c8e-456b-91b0-45f1a6baedf9
function plot_orbits_w_lmc(trajectories; lmc_visible=false, act_visible=false)
	fig = Figure(size=(6.5, 4) .* 72)
	r_max = 300
	t_min = -5/T2GYR
	
	o = orbit_nbody_lmc

  	ax2 = Axis(fig[1, 1], xlabel="x / kpc", ylabel="z / kpc",
        xgridvisible=false, ygridvisible=false, 
			  limits=(-1, 1, -1, 1) .* r_max,
			   aspect= DataAspect(),
    )

    for (i, (label, traj)) in enumerate(trajectories)
        plot_x_y_traj!(traj, x_direction=1, label=label=>(; alpha=0.5), 	
					   color=COLORS[i]; t_min=t_min)
    end

	
	if lmc_visible
		scatter!(pos_0_lmc[1], pos_0_lmc[3], markersize=8, color=:black)
		plot_traj_lmc!(x_direction=1, t_min=t_min)
	end
	
	scatter!(pos_0[1], pos_0[3], markersize=8, color=COLORS[1], strokewidth=lw_thick[]/2, strokecolor=:black)

	if act_visible
		lines!(o.positions[1, :], o.positions[3, :], linewidth=lw_thick, color=COLORS[3])
	end

	# poly!(xz_iso_v24..., color=COLORS[8])
	Utils.plot_sun!(x_direction=1)

	
    ax = Axis(fig[1, 2], xlabel="y / kpc", ylabel="z / kpc",
        xgridvisible=false, ygridvisible=false, 
        aspect=DataAspect(),
			  limits=(-1, 1, -1, 1) .* r_max,
				  xticks = -200:100:300

    )
    
    for (i, (label, traj)) in enumerate(trajectories)
        plot_x_y_traj!(traj, label=label=>(; alpha=0.5), color=COLORS[i], t_min=-5/T2GYR;)
    end

	if lmc_visible
		plot_traj_lmc!(t_min=t_min)
		scatter!(pos_0_lmc[2], pos_0_lmc[3], markersize=8, color=:black)
	end
	
	scatter!(pos_0[2], pos_0[3], markersize=8, color=COLORS[1], strokewidth=lw_thick[]/2, strokecolor=:black)

	
	text!(pos_0[2], pos_0[3], text="today", align=(:left, :center), offset=(5, 0), fontsize=0.8*theme(:fontsize)[])


	text!(0, 0, text="Sun", align=(:left, :center), offset=(5, 0), fontsize=0.8*theme(:fontsize)[])


	if act_visible
		lines!(o.positions[2, :], o.positions[3, :], linewidth=lw_thick, color=COLORS[3])
	end

	# poly!(yz_iso_v24..., color=COLORS[8])
	Utils.plot_sun!()
	
	hideydecorations!(ax, ticks=false, minorticks=false)

	rowsize!(fig.layout, 1, Aspect(1, 1.0))
	linkyaxes!(ax, ax2)
	

	axislegend(ax, merge=true, unique=true, position=:lb, )

	resize_to_layout!(fig)
	fig
end

# ╔═╡ 2683e73a-e9d9-4f79-a101-8f306d36e428
@savefig "scl_lmc_mc_orbits_xy" plot_orbits_w_lmc(trajectories, lmc_visible=true)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═45f1d142-cce5-4d81-830c-b61efb849e77
# ╠═c2c320b8-3b14-4e49-9048-92a546a6b275
# ╠═afb27a46-7b9f-4d61-8d22-2bc2566723e1
# ╠═a1606881-5138-4c90-b66f-34bcffd563eb
# ╠═cbc79255-9254-4e01-9353-c87e62620ec8
# ╠═f7d1d155-eaab-4589-b5e7-5efac2b8fb5f
# ╠═1dc687cf-f56d-4e28-92f5-f0218fd5ff43
# ╠═432dc7d3-c390-48cb-9e51-d4d0ff21620c
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╟─b9b99de2-1c01-4562-8295-75b770e232dc
# ╠═bd1e330a-a7b7-4273-87ba-7100cb6d5af8
# ╠═14f862cc-8e53-4d10-9d73-579ed1d26cd9
# ╠═18b45a96-ddcd-4a9d-88e2-4dde209676cf
# ╠═70d81f4a-8e01-45b1-a695-234d0249cca9
# ╠═6b0210cc-353d-48d4-bd59-f81fdd422d55
# ╟─38f593ae-590b-4f41-908d-75b17e4c4a88
# ╠═9d78acb1-5439-4338-aaff-2860f7f30b77
# ╠═27f32059-5af5-4b32-8b58-619b378acb19
# ╠═bfac7dd2-264c-448d-8755-725ebff4df95
# ╠═f4431bdf-b495-42ec-8aac-b7cf07f5746a
# ╠═8091ec68-c602-4d96-a418-d040aee84e41
# ╠═1ebe70d8-34b0-48d0-ad60-341223fb55fe
# ╠═d75f20c3-4a75-4a0b-b0fd-825c10113be8
# ╠═4529b1cc-5046-4491-8932-4c3fbb1f1e48
# ╟─6d17d26a-9946-4656-be7b-13d84002c3fa
# ╠═af5cc193-033b-48d2-b2f5-266f7c54c8d6
# ╠═70c6eab7-0f69-41b9-8c09-11ab0aa2de46
# ╠═93ca6a6b-5f41-40d1-8661-31d0b3c5aa58
# ╠═bc7f29c3-0368-48e9-b74d-728ed072df44
# ╠═83da41b8-ffb5-4bef-8d1e-175094d69d41
# ╟─a1c44523-e31e-4b50-bb00-db80ab080004
# ╠═f62f0c8b-bfa3-4bfe-8c02-b3d232f6fef6
# ╠═ff89d7db-c954-48b8-b87a-446ccfb2d79b
# ╠═eac841e7-5015-48cc-898b-f0232c8d0d72
# ╟─92fb5b1f-5dc6-4766-b69d-f8644861d4cd
# ╠═0f3018a2-42e2-4140-8e9e-f094276b76f2
# ╠═751ae520-59a2-471c-9247-90fba399cb42
# ╠═2683e73a-e9d9-4f79-a101-8f306d36e428
# ╠═f4bf5ad4-7c8e-456b-91b0-45f1a6baedf9
# ╠═6f1a1958-708a-4aac-935d-0e5a87091a98
