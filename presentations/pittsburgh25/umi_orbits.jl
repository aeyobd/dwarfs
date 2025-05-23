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

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./style.jl")

# ╔═╡ 15ef3ea9-fe17-44d1-924e-c18722ddf35a
import TOML

# ╔═╡ 20bc3e9a-ea85-46f3-9a56-1fae435a8575
set_theme!(theme_arya())

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

# ╔═╡ efae8808-1768-40c1-a11d-253402cdfe43
function compare_x_y_traj(trajectories; r_max=300, kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="y / kpc", ylabel="z / kpc",
        xgridvisible=false, ygridvisible=false, 
        aspect=DataAspect(),
			  limits=(-1, 1, -1, 1) .* r_max,
    )
    
    for (i, (label, traj)) in enumerate(trajectories)
        plot_x_y_traj!(traj, label=label=>(; alpha=0.5), color=COLORS[i]; kwargs...)
    end
        
    fig
end

# ╔═╡ 1014f681-8014-4df2-bd07-bce4f3348056
function compare_r_t_traj(trajectories; limits = (nothing, nothing), colors=COLORS, kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="time / Gyr", ylabel = "distance / kpc",
        xgridvisible=false, ygridvisible=false,
			  limits=limits
    )

    for (i, (label, traj)) in enumerate(trajectories)
        plot_r_t_traj!(traj, label=label, color=colors[i]; kwargs...)
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

# ╔═╡ a1c44523-e31e-4b50-bb00-db80ab080004
md"""
# Plots
"""

# ╔═╡ fb3f61f2-f425-40b4-9983-dfc7528df6ab
@savefig "umi_mc_orbits_xy" let

	fig = compare_x_y_traj(Dict("fiducial" => traj_fiducial), t_min=-5/T2GYR, r_max=120, thin=2,)

	scatter!(pos_0[2], pos_0[3], markersize=20, color=:black)
	
	text!(pos_0[2], pos_0[3], text="today", align=(:left, :center), offset=(24, 0), fontsize=0.8*theme(:fontsize)[])


	fig
end

# ╔═╡ 1d3f999b-b25d-49dc-b1ca-7dfdbf1d3d59
compare_r_t_traj(["" => traj_fiducial], thin=3)

# ╔═╡ 0c9eb9da-16c2-4d49-a798-7cd3faa5ec09
let
	fig = compare_x_y_traj(trajectories, t_min=-5/T2GYR, thin=2 )

	plot_x_y_traj!(traj_lmc, label="LMC", alpha=1, t_min=-5/T2GYR)

	axislegend(position=:lb, merge=true, unique=true, margin=(19,19,19,19))
	fig
end

# ╔═╡ 515bbb85-29e6-4fd9-a031-98fc45e0f426
let
	fig = compare_x_y_traj(trajectories, t_min=-10/T2GYR, thin=2 )

	plot_x_y_traj!(traj_lmc, label="LMC", alpha=1, t_min=-10/T2GYR)

	axislegend(position=:lt, merge=true, unique=true, margin=(19,19,19,19))
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

# ╔═╡ f48b494b-8b65-47d8-aa8f-3ea37e8a48ed
let 
	fig = compare_r_t_traj(["" => traj_fiducial], limits=(-5, 0, 0, 100), colors=[COLORS[1], COLORS[2], :black], thin=2)

	vlines!(-0.12, color=:black)

	fig
end

# ╔═╡ a876edb6-d0c4-4420-95b2-09c90f33a4b5
let 
	fig = compare_r_t_traj(trajectories_radii, limits=(-4, 0, 0, 300), colors=[COLORS[1], COLORS[2], :black], thin=2)

	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═15ef3ea9-fe17-44d1-924e-c18722ddf35a
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═20bc3e9a-ea85-46f3-9a56-1fae435a8575
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
# ╟─a1c44523-e31e-4b50-bb00-db80ab080004
# ╠═fb3f61f2-f425-40b4-9983-dfc7528df6ab
# ╠═1d3f999b-b25d-49dc-b1ca-7dfdbf1d3d59
# ╠═0c9eb9da-16c2-4d49-a798-7cd3faa5ec09
# ╠═515bbb85-29e6-4fd9-a031-98fc45e0f426
# ╠═c1062412-37ff-488e-a423-67c5adaa8e66
# ╠═34e480a6-8e9e-4278-80a2-c788ccb01327
# ╠═081e6a08-92e3-4f9d-bc8e-7ff80ca74693
# ╠═f48b494b-8b65-47d8-aa8f-3ea37e8a48ed
# ╠═a876edb6-d0c4-4420-95b2-09c90f33a4b5
