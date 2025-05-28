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

# ╔═╡ c2c320b8-3b14-4e49-9048-92a546a6b275
using HDF5

# ╔═╡ a1606881-5138-4c90-b66f-34bcffd563eb
using OrderedCollections

# ╔═╡ 2d94b40a-2a38-48d5-9b7e-6b493d504aec
using CSV, DataFrames

# ╔═╡ 60c30dac-1831-4ee7-ac8f-a797ed6bcbf8
using ForwardDiff

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./style.jl")

# ╔═╡ 08c7e69b-2c64-48fa-a12f-4d1c8d60a0d2
include("./utils.jl")

# ╔═╡ 15ef3ea9-fe17-44d1-924e-c18722ddf35a
import TOML

# ╔═╡ 7f77b542-8f03-4474-8e7f-f6353f97128b
let
	fig = Figure()
	ax = Axis(fig[1,1])
	x = sind.(0:1.05:360)
	y = cosd.(0:1.05:360)
	lines!(x, y)
	xlims!(-30, 30)
	ylims!(-30, 30)

	fig

end

# ╔═╡ b9b99de2-1c01-4562-8295-75b770e232dc
md"""
# Utils
"""

# ╔═╡ 68d1a6c9-128c-4abf-93e0-3f3e36431ccb
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor/observed_properties.toml"))

# ╔═╡ 940c54c4-8331-4c26-a8e8-8ef94486bb3b
icrs0 = ICRS(obs_props)

# ╔═╡ 9aca37e5-6d53-48d7-b4db-b84770e66a4a
gc0 = LilGuys.transform(LilGuys.Galactocentric, icrs0)

# ╔═╡ 37e8a3c3-03e0-45e6-890d-8071c3343c64
pos_0 = LilGuys.position(gc0)

# ╔═╡ 14f862cc-8e53-4d10-9d73-579ed1d26cd9
function read_traj(name)
    local positions, velocities, times
    
    h5open(joinpath(ENV["DWARFS_ROOT"], "analysis/ursa_minor/mc_orbits/$name/trajectory.hdf5"), "r") do f
        positions = f["positions"][:, :, :]
        velocities = f["velocities"][:, :, :]
        times = -f["times"][:]
    end

    return positions, velocities, times
end

# ╔═╡ 18b45a96-ddcd-4a9d-88e2-4dde209676cf
function read_lmc_traj(name)    
	df = CSV.read(joinpath(ENV["DWARFS_ROOT"], "analysis/ursa_minor/mc_orbits/$name/lmc_traj.csv"), DataFrame)

	positions = [df.x df.y df.z]'
	velocities = [df.v_x df.v_y df.v_z]'
	times = -df.time
    return reshape(positions, (3, 1, :)), reshape(velocities, (3, 1, :)), times
end

# ╔═╡ af5cc193-033b-48d2-b2f5-266f7c54c8d6
traj_lmc = read_lmc_traj("vasiliev24_L3M11")

# ╔═╡ 70c6eab7-0f69-41b9-8c09-11ab0aa2de46
traj_fiducial = read_traj("EP2020")

# ╔═╡ 93ca6a6b-5f41-40d1-8661-31d0b3c5aa58
traj_w_lmc = read_traj("vasiliev24_L3M11")

# ╔═╡ bc7f29c3-0368-48e9-b74d-728ed072df44
traj_nolmc = read_traj("vasiliev24_M11")

# ╔═╡ 088dceb9-fd77-4c95-9040-5f30f3fe1a53
function plot_r_t_traj!(traj; alpha=0.01, thin=1, color=:black, kwargs...)
    positions, velocities, times = traj
    for i in 1:thin:size(positions, 2)
        x = times * T2GYR
        y = radii(positions[:, i, :])
        lines!(x, y; rasterize=true, alpha=alpha, color=color, kwargs...)
    
    end
end

# ╔═╡ 70d81f4a-8e01-45b1-a695-234d0249cca9
function plot_x_y_traj!(traj; thin=1, x_direction=2, y_direction=3, alpha=0.01, color=:black, t_min = -2/T2GYR, kwargs...)
    positions, velocities, times = traj

	filt = times .> t_min 
    for i in 1:thin:size(positions, 2)
        x = positions[x_direction, i, filt]
        y = positions[y_direction, i, filt]
        
        lines!(x, y; rasterize=true, alpha=alpha, color=color, kwargs...)
    
    end
end

# ╔═╡ 1014f681-8014-4df2-bd07-bce4f3348056
function compare_r_t_traj(trajectories; limits = (nothing, nothing), colors=COLORS, kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="time / Gyr", ylabel = "distance / kpc",
        xgridvisible=false, ygridvisible=false,
			  limits=limits
    )

    for (i, (label, traj)) in enumerate(trajectories)
        plot_r_t_traj!(traj, label=label =>(; alpha=0.5, ), color=colors[i]; kwargs...)
    end
    
    fig
end

# ╔═╡ ff89d7db-c954-48b8-b87a-446ccfb2d79b
trajectories = OrderedDict(
	#"fiducial" => traj_fiducial,
	# "nolmc" => traj_nolmc,
	"no LMC" => traj_nolmc,
	 "w/ LMC" => traj_w_lmc,
)

# ╔═╡ 859df394-8b33-410a-9979-ef5fdac72369
CairoMakie.activate!(type=:png)

# ╔═╡ 07d86e69-6c9f-4b0a-a1da-a63aad3f9b0f
import Agama

# ╔═╡ 0e02e7b9-294d-4849-a16c-d4719c07b626
pot = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama", "potentials", "EP2020.ini"))

# ╔═╡ fc35faac-7df4-42ef-8ff2-6455ae99c340
pot_stars = Agama.Potential(pot._py[0] + pot._py[1] + pot._py[2])

# ╔═╡ efae8808-1768-40c1-a11d-253402cdfe43
function compare_x_y_traj(trajectories; r_max=300, kwargs...)
    fig = Figure()
  	ax2 = Axis(fig[1, 1], xlabel="x / kpc", ylabel="z / kpc",
        xgridvisible=false, ygridvisible=false, 
			  limits=(-1, 1., -1, 1) .* r_max,
    )

    for (i, (label, traj)) in enumerate(trajectories)
        plot_x_y_traj!(traj, x_direction=1, label=label=>(; alpha=0.5), color=COLORS[i]; kwargs...)
    end

	plot_isodensity!(pot_stars, color=:black, x_direction=1)
	plot_sun!(x_direction=1)

	
    ax = Axis(fig[1, 2], xlabel="y / kpc", ylabel="z / kpc",
        xgridvisible=false, ygridvisible=false, 
			  limits=(-1, 1, -1, 1) .* r_max,
			  aspect = DataAspect()
    )
    
    for (i, (label, traj)) in enumerate(trajectories)
        plot_x_y_traj!(traj, label=label=>(; alpha=0.5), color=COLORS[i]; kwargs...)
    end
	hideydecorations!(ticks=false, minorticks=false)

	plot_isodensity!(pot_stars, color=:black)
	plot_sun!()
	
	colsize!(fig.layout, 1, Aspect(1, 1.0))
	colsize!(fig.layout, 2, Aspect(1, 1.0))
	linkyaxes!(ax, ax2)




    fig
end

# ╔═╡ 044b353c-6d00-47f3-9def-7f9fa9e361fa
let
	fig = Figure()

	ax = Axis(fig[1,1], aspect=DataAspect())
	plot_isodensity!(pot_stars, plot_reflection=true)
	xlims!(-15, 15)


	fig
end

# ╔═╡ 9d64514b-0f59-4cd3-9f17-5dfabaacb707
let
	fig = Figure()

	ax = Axis(fig[1,1], aspect=DataAspect())
	plot_isodensity!(pot_stars, x_direction=1, y_direction=2, plot_reflection=true)

plot_sun!(x_direction=1, y_direction=2)
	fig
end

# ╔═╡ a1c44523-e31e-4b50-bb00-db80ab080004
md"""
# Plots
"""

# ╔═╡ ce0314a1-a07e-4608-9ee4-a2db61dc6043
function get_nbody_orbit(model)
	positions, times = HDF5.h5open(joinpath(ENV["DWARFS_ROOT"], "analysis/$model/centres.hdf5")) do f
		return f["positions"][:, :], f["times"][:]
	end

	idx_f = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis/$model/orbital_properties.toml"))["idx_f"]

	return positions[:, 1:idx_f], times[1:idx_f]
end
						

# ╔═╡ 61e3ae56-873c-48a4-8058-34f792660df8
orbit_nbody = get_nbody_orbit("ursa_minor/1e6_v38_r4.0/orbit_smallperi.3")

# ╔═╡ 43ef72f8-7d7e-4976-82ce-74831bf0cbeb
@savefig "umi_mc_orbits_xy_act" let

	fig = compare_x_y_traj(Dict("fiducial" => traj_fiducial), t_min=-5/T2GYR, r_max=120, thin=2,)

	scatter!(pos_0[2], pos_0[3],  color=:black)
	ax = fig.content[1]


	
	text!(pos_0[2], pos_0[3], text="today", align=(:left, :center), offset=(24, 0), fontsize=0.8*theme(:fontsize)[])

	o = orbit_nbody[1][:, orbit_nbody[2] .- orbit_nbody[2][end] .> -5/T2GYR]
	
	lines!(ax, o[1, :], o[3, :])

	Makie.current_axis!(fig.content[2])

	lines!(o[2, :], o[3, :])
	scatter!(ax, pos_0[1], pos_0[3],  color=:black)


	CairoMakie.current_axis!(ax)
	plot_sun!(x_direction=1)

	resize_to_layout!(fig)

	fig
end

# ╔═╡ 1d3f999b-b25d-49dc-b1ca-7dfdbf1d3d59
let
	fig = compare_r_t_traj(["" => traj_fiducial], thin=3)
	lines!((orbit_nbody[2] .- orbit_nbody[2][end]) * T2GYR, radii(orbit_nbody[1]))

	fig
end

# ╔═╡ 0c9eb9da-16c2-4d49-a798-7cd3faa5ec09
@savefig "umi_mc_orbits_w_lmc_yz" let
	fig = Figure()

	r_max=150
	
    ax = Axis(fig[1, 1], xlabel="y / kpc", ylabel="z / kpc",
        xgridvisible=false, ygridvisible=false, 
			  limits=(-1, 1, -1, 1) .* r_max,
			  aspect = DataAspect()
    )
    
	  for (i, (label, traj)) in enumerate(trajectories)
        plot_x_y_traj!(traj, label=label=>(; alpha=0.5), color=COLORS[i], t_min=-5/T2GYR)
    end

	plot_x_y_traj!(traj_lmc, label="LMC", alpha=1, t_min=-5/T2GYR)
	axislegend(position=:lb, merge=true, unique=true, margin=(19,19,19,19))
	plot_isodensity!(pot_stars)
	plot_sun!()



	fig
end

# ╔═╡ 515bbb85-29e6-4fd9-a031-98fc45e0f426
@savefig "umi_mc_orbits_w_lmc" let
	fig = compare_x_y_traj(trajectories, t_min=-10/T2GYR, thin=2 )
	
	#ax.xticks = collect(-100:100:100)
	ax = fig.content[1]

	plot_x_y_traj!(traj_lmc, label="LMC", alpha=1, t_min=-10/T2GYR)
	axislegend(position=:rb, merge=true, unique=true, margin=(19,19,19,19))

	CairoMakie.current_axis!(fig.content[1])

	plot_x_y_traj!(traj_lmc, label="LMC", alpha=1, x_direction=1, t_min=-10/T2GYR)


	resize_to_layout!(fig)

	fig
end

# ╔═╡ c1062412-37ff-488e-a423-67c5adaa8e66
theme(:Legend)

# ╔═╡ 34e480a6-8e9e-4278-80a2-c788ccb01327
function subtract_lmc(traj, traj_lmc)
	@assert all(traj[3] .≈ traj_lmc[3])
	pos_scl_lmc = traj[1] .- reshape(traj_lmc[1], (3, 1, :))
	vel_lmc = traj[2] .- reshape(traj_lmc[2], (3, 1, :))

	return pos_scl_lmc, vel_lmc, traj[3]
end

# ╔═╡ 081e6a08-92e3-4f9d-bc8e-7ff80ca74693
trajectories_radii = OrderedDict(
	"no lmc" => traj_nolmc,
	 "lmc" => traj_w_lmc,

)

# ╔═╡ a876edb6-d0c4-4420-95b2-09c90f33a4b5
let 
	fig = compare_r_t_traj(trajectories_radii, limits=(-10, 0, 0, 300), colors=[COLORS[1], COLORS[2], :black], thin=2)

	axislegend(position=:rt, merge=true, unique=true)
	
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═15ef3ea9-fe17-44d1-924e-c18722ddf35a
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═08c7e69b-2c64-48fa-a12f-4d1c8d60a0d2
# ╠═7f77b542-8f03-4474-8e7f-f6353f97128b
# ╠═c2c320b8-3b14-4e49-9048-92a546a6b275
# ╠═a1606881-5138-4c90-b66f-34bcffd563eb
# ╠═2d94b40a-2a38-48d5-9b7e-6b493d504aec
# ╟─b9b99de2-1c01-4562-8295-75b770e232dc
# ╠═68d1a6c9-128c-4abf-93e0-3f3e36431ccb
# ╠═940c54c4-8331-4c26-a8e8-8ef94486bb3b
# ╠═9aca37e5-6d53-48d7-b4db-b84770e66a4a
# ╠═37e8a3c3-03e0-45e6-890d-8071c3343c64
# ╠═14f862cc-8e53-4d10-9d73-579ed1d26cd9
# ╠═18b45a96-ddcd-4a9d-88e2-4dde209676cf
# ╠═af5cc193-033b-48d2-b2f5-266f7c54c8d6
# ╠═70c6eab7-0f69-41b9-8c09-11ab0aa2de46
# ╠═93ca6a6b-5f41-40d1-8661-31d0b3c5aa58
# ╠═bc7f29c3-0368-48e9-b74d-728ed072df44
# ╠═efae8808-1768-40c1-a11d-253402cdfe43
# ╠═088dceb9-fd77-4c95-9040-5f30f3fe1a53
# ╠═70d81f4a-8e01-45b1-a695-234d0249cca9
# ╠═1014f681-8014-4df2-bd07-bce4f3348056
# ╠═ff89d7db-c954-48b8-b87a-446ccfb2d79b
# ╠═859df394-8b33-410a-9979-ef5fdac72369
# ╠═07d86e69-6c9f-4b0a-a1da-a63aad3f9b0f
# ╠═0e02e7b9-294d-4849-a16c-d4719c07b626
# ╠═60c30dac-1831-4ee7-ac8f-a797ed6bcbf8
# ╠═fc35faac-7df4-42ef-8ff2-6455ae99c340
# ╠═044b353c-6d00-47f3-9def-7f9fa9e361fa
# ╠═9d64514b-0f59-4cd3-9f17-5dfabaacb707
# ╟─a1c44523-e31e-4b50-bb00-db80ab080004
# ╠═ce0314a1-a07e-4608-9ee4-a2db61dc6043
# ╠═61e3ae56-873c-48a4-8058-34f792660df8
# ╠═43ef72f8-7d7e-4976-82ce-74831bf0cbeb
# ╠═1d3f999b-b25d-49dc-b1ca-7dfdbf1d3d59
# ╠═0c9eb9da-16c2-4d49-a798-7cd3faa5ec09
# ╠═515bbb85-29e6-4fd9-a031-98fc45e0f426
# ╠═c1062412-37ff-488e-a423-67c5adaa8e66
# ╠═34e480a6-8e9e-4278-80a2-c788ccb01327
# ╠═081e6a08-92e3-4f9d-bc8e-7ff80ca74693
# ╠═a876edb6-d0c4-4420-95b2-09c90f33a4b5
