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

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ 1e5451e3-87ec-4c04-ae1a-87a2f62c2891
import TOML

# ╔═╡ 76cf78d3-5f15-48f1-927f-37266bd49865
let 
	fig = Figure()
	ax = Axis(fig[1,1])
scatter!([1], [1], label="hi")
axislegend(margin=(20,20,20,20))
fig
end

# ╔═╡ a2d26cab-b562-4d4c-88f7-5f5b8ee49116
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/observed_properties.toml"))

# ╔═╡ e095772c-7c5e-4dec-b29a-4cad97683a7d
icrs0 = ICRS(obs_props)

# ╔═╡ b9f1f612-efb3-4b8a-9b3e-1d7d61e08726
gc0 = LilGuys.transform(LilGuys.Galactocentric, icrs0)

# ╔═╡ 71e3de8e-91e4-4214-9602-89572ce7fa70
pos_0 = LilGuys.position(gc0)

# ╔═╡ b9b99de2-1c01-4562-8295-75b770e232dc
md"""
# Utils
"""

# ╔═╡ 14f862cc-8e53-4d10-9d73-579ed1d26cd9
function read_traj(name)
    local positions, velocities, times
    
    h5open(joinpath(ENV["DWARFS_ROOT"], "analysis/sculptor/mc_orbits/$name/trajectory.hdf5"), "r") do f
        positions = f["positions"][:, :, :]
        velocities = f["velocities"][:, :, :]
        times = -f["times"][:]
    end

    return positions, velocities, times
end

# ╔═╡ 18b45a96-ddcd-4a9d-88e2-4dde209676cf
function read_lmc_traj(name)    
	df = CSV.read(joinpath(ENV["DWARFS_ROOT"], "analysis/sculptor/mc_orbits/$name/lmc_traj.csv"), DataFrame)

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
        
        lines!(x, y; rasterize=true, alpha=alpha, color=color, linewidth=3, kwargs...)
    
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
        plot_r_t_traj!(traj, label=label => (; alpha=1.0), color=colors[i]; kwargs...)
    end
    
    fig
end

# ╔═╡ ff89d7db-c954-48b8-b87a-446ccfb2d79b
trajectories = OrderedDict(
	#"fiducial" => traj_fiducial,
	# "nolmc" => traj_nolmc,
	"MW only" => traj_nolmc,
	 "MW & LMC" => traj_w_lmc,
)

# ╔═╡ a1c44523-e31e-4b50-bb00-db80ab080004
md"""
# Plots
"""

# ╔═╡ c4761049-4e08-4fce-85f9-6970f6530e4b
theme(:markersize)

# ╔═╡ fb3f61f2-f425-40b4-9983-dfc7528df6ab
@savefig "scl_mc_orbits_xy" let

	fig = compare_x_y_traj(Dict("fiducial" => traj_fiducial), t_min=-5/T2GYR, r_max=120, thin=2,)

	scatter!(pos_0[2], pos_0[3], markersize=20, color=:black)
	
	text!(pos_0[2], pos_0[3], text="today", align=(:left, :center), offset=(24, 0), fontsize=0.8*theme(:fontsize)[])


	fig
end

# ╔═╡ 0272f30a-49f5-405c-9bcb-4c03f32d697d
let
	fig = Figure()
	ax = Axis(fig[1,1])
		plot_x_y_traj!(traj_lmc, label="LMC", color=-filter(x->x > -4/T2GYR, Vector(traj_lmc[3])), colormap=:greys, linewidth=6, t_min=-4/T2GYR, alpha=1)

	fig
end

# ╔═╡ b2709a10-2bf6-4f60-8167-ea56c67f218a
filter(x->x > -4/T2GYR, Vector(traj_lmc[3]))

# ╔═╡ 56b4119d-3d74-483e-9340-c98512cda640
theme(:size)

# ╔═╡ 0c9eb9da-16c2-4d49-a798-7cd3faa5ec09
@savefig "scl_lmc_mc_orbits_xy" let
    fig = Figure(size=(1000, 645))
	
    ax = Axis(fig[1, 1], xlabel="y / kpc", ylabel="z / kpc",
        xgridvisible=false, ygridvisible=false, 
        aspect=DataAspect(),
			  limits=(-1, 1, -1, 1) .* 300,
    )
    
    for (i, (label, traj)) in enumerate(trajectories)
        plot_x_y_traj!(traj, label=label=>(; alpha=0.5), color=COLORS[i], t_min=-5/T2GYR;)
    end
        
		plot_x_y_traj!(traj_lmc, label="LMC", color=-filter(x->x > -5/T2GYR, Vector(traj_lmc[3])), colormap=:greys, linewidth=6, t_min=-5/T2GYR, alpha=1)
	scatter!(pos_0_lmc[2], pos_0_lmc[3], markersize=20, color=:black)

	scatter!(pos_0[2], pos_0[3], markersize=20, color=:black)
	text!(pos_0[2], pos_0[3], text="today", align=(:left, :center), offset=(12, 0), fontsize=0.8*theme(:fontsize)[])

	Legend(fig[1,2], Makie.current_axis(), merge=true, unique=true)
	fig
end

# ╔═╡ 865f71c1-0cf8-4a63-b417-a2aed5c5a3ae
traj_lmc[3]

# ╔═╡ bb635c9c-facf-4fbf-b9c9-2445d017c046
pos_0_lmc = traj_lmc[1][:, 1, 1]

# ╔═╡ dfecdf9a-b28c-4561-991e-d555c4fdcb39
theme(:linewidth)

# ╔═╡ 34e480a6-8e9e-4278-80a2-c788ccb01327
function subtract_lmc(traj, traj_lmc)
	@assert all(traj[3] .≈ traj_lmc[3])
	pos_scl_lmc = traj[1] .- reshape(traj_lmc[1], (3, 1, :))
	vel_lmc = traj[2] .- reshape(traj_lmc[2], (3, 1, :))

	return pos_scl_lmc, vel_lmc, traj[3]
end

# ╔═╡ 0f3018a2-42e2-4140-8e9e-f094276b76f2
traj_scl_lmc = subtract_lmc(traj_w_lmc, traj_lmc)

# ╔═╡ 081e6a08-92e3-4f9d-bc8e-7ff80ca74693
trajectories_radii = OrderedDict(
	"MW only" => traj_nolmc,
	 "MW & LMC" => traj_w_lmc,
	#"scl - lmc" => traj_scl_lmc

)

# ╔═╡ aa0abdb1-325d-4dd6-b91a-adc4a758bb2a
theme(:Legend)

# ╔═╡ 702369bd-7e44-4309-9d52-2b02be7422fa
theme(:Legend)

# ╔═╡ baf713c2-fa6e-4c26-9f7d-9025c05b58fb
R_b_obs = LilGuys.arcmin2kpc(30, obs_props["distance"])

# ╔═╡ 3a561930-fadd-4cfc-8980-cd5b35d5e2ed
 R_b_obs / ( obs_props["sigma_v"] / V2KMS) * T2GYR / 0.55

# ╔═╡ daab8d6c-5c0b-462d-af91-1d26848eda15
t_break = 0.130

# ╔═╡ a876edb6-d0c4-4420-95b2-09c90f33a4b5
let 
	#compare_r_t_traj(trajectories_radii, , , thin=2)
	colors=[COLORS[1], COLORS[2], :black]
	limits=(-5, 0, 0, 300)
	
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="time / Gyr", ylabel = "distance to MW / kpc",
        xgridvisible=false, ygridvisible=false,
			  limits=limits,
			  xticks = -5:0,
			  xminorticks = -5:0.2:0,
			  yticks = 0:100:300,
			  yminorticks = 0:20:300,

    )

    for (i, (label, traj)) in enumerate(trajectories_radii)
        plot_r_t_traj!(traj, label=label => (; alpha=1.0), color=colors[i]; thin=2)
    end
    
    

	axislegend(merge=true, unique=true, margin=(19, 19, 19, 19))

		vlines!(-t_break, color=:black)

	fig
end

# ╔═╡ 05498dce-6458-4254-bb51-5eb016a44de5
let 
	#compare_r_t_traj(trajectories_radii, , , thin=2)
	colors=[COLORS[1], COLORS[2], :black]
	limits=(-5, 0, 0, 300)
	
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="time / Gyr", ylabel = "distance to LMC / kpc",
        xgridvisible=false, ygridvisible=false,
			  limits=limits
    )

    plot_r_t_traj!(traj_scl_lmc, color=COLORS[2]; thin=2)
    
	vlines!(-t_break, color=:black)
	text!(-t_break, 210, text="observed break", align=(:center, :bottom), rotation=π/2, fontsize=0.8*theme(:fontsize)[])

	fig
end

# ╔═╡ 4c227fe5-dbb3-4f92-8eb5-2f0de0c6c209
LilGuys.break_radius(0.130/T2GYR, obs_props["sigma_v"] / V2KMS)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═c2c320b8-3b14-4e49-9048-92a546a6b275
# ╠═a1606881-5138-4c90-b66f-34bcffd563eb
# ╠═2d94b40a-2a38-48d5-9b7e-6b493d504aec
# ╠═1e5451e3-87ec-4c04-ae1a-87a2f62c2891
# ╠═76cf78d3-5f15-48f1-927f-37266bd49865
# ╠═a2d26cab-b562-4d4c-88f7-5f5b8ee49116
# ╠═e095772c-7c5e-4dec-b29a-4cad97683a7d
# ╠═b9f1f612-efb3-4b8a-9b3e-1d7d61e08726
# ╠═71e3de8e-91e4-4214-9602-89572ce7fa70
# ╟─b9b99de2-1c01-4562-8295-75b770e232dc
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
# ╠═c4761049-4e08-4fce-85f9-6970f6530e4b
# ╠═fb3f61f2-f425-40b4-9983-dfc7528df6ab
# ╠═0272f30a-49f5-405c-9bcb-4c03f32d697d
# ╠═b2709a10-2bf6-4f60-8167-ea56c67f218a
# ╠═56b4119d-3d74-483e-9340-c98512cda640
# ╠═0c9eb9da-16c2-4d49-a798-7cd3faa5ec09
# ╠═865f71c1-0cf8-4a63-b417-a2aed5c5a3ae
# ╠═bb635c9c-facf-4fbf-b9c9-2445d017c046
# ╠═dfecdf9a-b28c-4561-991e-d555c4fdcb39
# ╠═34e480a6-8e9e-4278-80a2-c788ccb01327
# ╠═0f3018a2-42e2-4140-8e9e-f094276b76f2
# ╠═081e6a08-92e3-4f9d-bc8e-7ff80ca74693
# ╠═aa0abdb1-325d-4dd6-b91a-adc4a758bb2a
# ╠═a876edb6-d0c4-4420-95b2-09c90f33a4b5
# ╠═702369bd-7e44-4309-9d52-2b02be7422fa
# ╠═05498dce-6458-4254-bb51-5eb016a44de5
# ╠═baf713c2-fa6e-4c26-9f7d-9025c05b58fb
# ╠═3a561930-fadd-4cfc-8980-cd5b35d5e2ed
# ╠═daab8d6c-5c0b-462d-af91-1d26848eda15
# ╠═4c227fe5-dbb3-4f92-8eb5-2f0de0c6c209
