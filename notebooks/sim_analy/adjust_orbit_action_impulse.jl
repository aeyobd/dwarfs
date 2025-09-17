### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ b18d1734-e40e-11ef-1007-314371eb1a54
begin
	using Pkg;Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ 0b6b20b0-bba8-4ed3-b8a5-d8361bc3674f
using PlutoUI

# ╔═╡ 4bc528f9-4ac0-407e-88a6-5e606117a3d3
using CSV, DataFrames

# ╔═╡ e13b1c88-75c4-498a-87c5-6d11c88fa569
using OrderedCollections

# ╔═╡ 481edcca-51a8-4f11-88ec-84939b223bdf
include(ENV["DWARFS_ROOT"] * "/utils/agama_utils.jl")

# ╔═╡ 9bfa1465-2ca9-495c-bab6-5110072f02ed
md"""
# Adjusting Orbits post impact

This notebook is specifically designed to solve the impact orbit deviation problem

For my so-called MW-impact model, the galaxy passes right through the inner milky way. This notebook uses the assumption that the kinematic deviations of the orbit from a point trajectory happen instantaneously at pericentre. As a result, we can adjust the orbit by accounting for the action loss in the orbit prior to the impact. 

In detail, this notebook calculates the deviations in action space from the desired orbit assuming the initial, axisymmetric potential. Actions are more challenging to interpret in the evolving potential. 
"""

# ╔═╡ 9c404056-2980-4620-9e4f-459157533c77
units = Agama.VASILIEV_UNITS

# ╔═╡ e9c9d1d9-42e6-4732-b425-477ec0bd7589
md"""
# Setup
"""

# ╔═╡ 34e54152-cebb-4364-b1b7-37e0dd4f10c8
import LinearAlgebra: ⋅, ×

# ╔═╡ f0f7b2db-b926-46d5-a9bc-baa9cff0b6cf
import TOML

# ╔═╡ bab15d74-24e2-4bec-a9c7-39163297b360
CairoMakie.activate!(type=:png, px_per_unit=2)

# ╔═╡ 598be49e-82ef-4820-8996-ea2b753f4275
function notebook_inputs(; kwargs...)
	return PlutoUI.combine() do Child
		
		user_inputs = [
			md""" $(string(name)): $(
				Child(name, obj)
			)"""
			
			for (name, obj) in kwargs
		]
		
		md"""
		#### Inputs
		$(user_inputs)
		"""
	end
end

# ╔═╡ 5006f134-648a-4318-9dc3-5de383ac4d0e
@bind inputs confirm(notebook_inputs(;
	galaxyname = TextField(default="sculptor"),
	haloname = TextField(default="1e6_new_v31_r3.2"),
	orbitname = TextField(default="L3M11_9Gyr_smallperi.x"),
	potname = TextField(default = "simulation/agama_potential.ini")
))

# ╔═╡ fba1564a-c12d-4657-86fe-d4d771be1037
galaxyname = inputs.galaxyname

# ╔═╡ 0b6b2020-282b-4384-ac58-8665364cf86e
haloname = inputs.haloname

# ╔═╡ 06fbb3db-26ae-4380-935c-883a02c6e119
orbitname =  inputs.orbitname

# ╔═╡ c8e2c46d-917c-4fa8-a841-ea26430d17d3
modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, orbitname) * "/"

# ╔═╡ fa2ccfd9-59ee-4484-88ae-997afa83ec86
FIGDIR = joinpath(modeldir, "figures")

# ╔═╡ b5302ffe-76e1-424d-9e08-1656981b068e
md"""
# Data loading
"""

# ╔═╡ 5c8ac63e-3661-455b-a434-7415981b9d89
md"""
## Load in the potentials
We also need a `potential_mw_init.ini` in the analysis directory which is the chosen axisymetric static potential
"""

# ╔═╡ de8bdfef-efca-4d2a-a9cd-a79ccc931b3f
pot = Agama.Potential(file=joinpath(modeldir  * inputs.potname))

# ╔═╡ afeb4045-01be-42e1-9998-36df9b197bb9
pot_static = Agama.Potential(file=joinpath(modeldir  * "/potential_mw_init.ini"))

# ╔═╡ 7bcdef79-7391-4f35-902b-62d2dda1d858
Φ_in(x, t=0) = Agama.potential(pot, x,units, t=t)

# ╔═╡ f9c28bef-7beb-42c7-a708-c32dc1d43924
md"""
## Load in the orbits
- `orbit_old` is the expected point particle orbit. We want the final simulation to end up in the same place
- `orbit_nbody` is the result from the nbody model
- `times` comes from orbit nbody for convenience
"""

# ╔═╡ 7ada86e5-bc6e-4157-b9cd-a7cf154d25d0
orbit_nbody = Orbit(joinpath(modeldir, "centres.hdf5"))

# ╔═╡ 7fe3e7a2-3ad9-4f37-8a9c-0a50c8aa0c01
times = orbit_nbody.times

# ╔═╡ 087ae8e4-47ea-4c16-8f70-b38168834268
orbit_old = LilGuys.resample(Orbit(modeldir * "orbit_original.csv"), times)

# ╔═╡ c1acdfae-e69d-4c62-9bc7-1a0e1f0f6492
orbit_point = LilGuys.resample(Orbit(modeldir * "simulation/orbit.csv"), times)

# ╔═╡ acf9bb4c-c03f-455f-bd12-65766f55bfaa
md"""
# Observational reference
"""

# ╔═╡ f7659682-7f06-4618-aa70-4f2c7d3ee838
obs_props_file = modeldir * "simulation/orbit.toml"

# ╔═╡ 64d8c003-8b5d-4eb4-adb3-ccfae58f08d0
obs_props = TOML.parsefile(obs_props_file)

# ╔═╡ 51a6bc2e-7e1c-4dad-8a88-40015dc3a1f6
gc_today = LilGuys.transform(Galactocentric, ICRS(obs_props))

# ╔═╡ 0706ee1c-74a6-4a5d-907f-8884895b0361
E_obs = Φ_in(LilGuys.position(gc_today)) .+ 1/2 * radii(LilGuys.velocity(gc_today)/V2KMS) .^ 2

# ╔═╡ 7131d149-c2b6-4ea4-a022-954102bf0870
md"""
## For the orbits
"""

# ╔═╡ 5cd10512-5472-4cd5-b9ef-2b5eb0f4bb0c
md"""
# Inputs for analysis
"""

# ╔═╡ c7aa2829-9048-42b0-9348-3e36714e8e37
md"""
time range: $(@bind time_range confirm(NumberField(0:10:1000, default=100)))
"""

# ╔═╡ 6a22bbb7-f070-4cc5-9328-9e21cac65487
md"""
time smoothing: $(@bind window confirm(NumberField(0:1:100, default=0)))
"""

# ╔═╡ d7ca78c0-3898-4bd8-b28f-115644bcf314
md"""
# Action differences
"""

# ╔═╡ cc180e47-3dca-4209-9ac6-1aa4031faca4
function get_actions(pot, orbit; units=units)
	return get_actions(pot, orbit.positions, orbit.velocities, units=units)
end

# ╔═╡ d83bda34-f614-459d-8fb5-283a11b371f1
function get_actions(pot, pos, vel; units=units)
	af = Agama.ActionFinder(pot)

	return Agama.actions_angles(af, pos, vel, units)[[1, 2]]
end

# ╔═╡ 3d33f8f9-16a4-47ac-be22-d35f09a3ee7e
time_shift = -(orbit_nbody.times[argmin(radii(orbit_nbody))] - orbit_old.times[argmin(radii(orbit_old))])

# ╔═╡ 1e1f37a8-c3d0-4f78-b542-61b8e0814a9c
orbit_old_resampled = LilGuys.resample(orbit_old, orbit_nbody.times .+ time_shift)

# ╔═╡ ae7444ec-b934-4b98-966d-1b55760cc592
orbit_point_resampled = LilGuys.resample(orbit_point, orbit_nbody.times .+ time_shift)

# ╔═╡ 8b36e2c4-9073-4af7-a9f4-7ea391e12e56
orbit_nbody.times[argmin(radii(orbit_nbody))] - orbit_old_resampled.times[argmin(radii(orbit_old_resampled))]

# ╔═╡ 65afad37-79b3-408d-87ca-492bd55d4eca
J_old, Theta_old = get_actions(pot_static, orbit_old_resampled)

# ╔═╡ 09d3bbb1-31f6-42b0-9362-4a1defe4958a
J_point, Theta_point = get_actions(pot_static, orbit_point_resampled)

# ╔═╡ 5c4c0767-6c17-43e4-a9cf-823fde8a4228
J_nbody, Theta_nbody = get_actions(pot_static, orbit_nbody)

# ╔═╡ 01060309-2889-4055-8d53-b127cb88d3ea
t_encounter = orbit_nbody.times[argmin(radii(orbit_nbody))]

# ╔═╡ fa257945-41dd-481b-b719-ffffe7b87d67
md"""
### Sanity check
just check that the new orbit has a similar pericentric time
"""

# ╔═╡ cc8d1ced-99a0-41af-8fe8-f4c5bd8eef91
orbit_old.times[argmin(radii(orbit_old))]

# ╔═╡ 14667772-b5cd-4d79-9278-1451bf1bdeaa
md"""
## Plots
"""

# ╔═╡ 81ff5445-40fb-41e2-b4f8-d0bba8e702a4
idx_before = argmin(abs.(orbit_nbody.times .- t_encounter .+ time_range/2))

# ╔═╡ ec74d9e2-b2ae-447d-b9c6-5d55dd464d28
pos_i_point = LilGuys.positions(orbit_point)[:, idx_before]

# ╔═╡ fcb9b299-d660-48ba-afec-93c52a40d930
vel_i_point = LilGuys.velocities(orbit_point)[:, idx_before]

# ╔═╡ 1e4c0c1e-6809-4083-aec1-d3d3a2374068
idx_after = argmin(abs.(orbit_nbody.times .- t_encounter .- time_range/2))

# ╔═╡ 32e71749-6fd2-4640-838f-ece1ee354c6b
let
	fig = Figure(size=(6*72, 5*72))

	local ax
	for i in 1:3
		coord = ["R", "z", "\\phi"][i]
		
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"\Delta J_%$coord"
		)


		lines!(orbit_nbody.times, J_nbody[i, :] .- J_old[i, :], label = "nbody - desired")
		
		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end


		vlines!(orbit_nbody.times[[idx_before, idx_after]])
		hlines!(0, color=:black)
	end
	#ylims!(-1, 1)
	linkyaxes!(fig.content...)
	axislegend(fig.content[1])
	fig
end


# ╔═╡ 4dff152a-4426-4e47-addc-4ea3bf9157f0
LilGuys.mean(J_old[:, idx_after:idx_after+window] - J_nbody[:, idx_after:idx_after+window], dims=2)

# ╔═╡ e57d04f4-3908-4316-ad61-407bd1909b9b
dJ_suggested =  J_old[:, idx_after] - J_nbody[:, idx_after]

# ╔═╡ 77b485d5-da49-437f-8117-52a57f276090
LilGuys.mean(Theta_old[:, idx_after:idx_after+window] - Theta_nbody[:, idx_after:idx_after+window], dims=2)

# ╔═╡ 71c78bfc-5ba6-437f-b6dc-462db6840b08
dθ_suggested = Theta_old[:, idx_after] - Theta_nbody[:, idx_after] 

# ╔═╡ 94323325-0d5c-43b8-8793-e28ebad685d9
@bind change_in_act_angles confirm(notebook_inputs(;
	Jr = NumberField(-60:0.01:60, default=dJ_suggested[1]),
	Jz = NumberField(-30:0.01:30, default=dJ_suggested[2]),
	Jϕ = NumberField(-30:0.01:30, default=dJ_suggested[3]),
	Θr = NumberField(-3.15:0.001:3.15, default=dθ_suggested[1]),
	Θz = NumberField(-3.15:0.001:3.15, default=dθ_suggested[2]),
	Θϕ = NumberField(-3.15:0.001:3.15, default=dθ_suggested[3]),
))

# ╔═╡ 4d2d9e31-54eb-45cb-b228-9375a829c5a9
dJ = [change_in_act_angles[:Jr], change_in_act_angles[:Jz], change_in_act_angles[:Jϕ]]

# ╔═╡ 5ad6cf47-3ea5-41e1-81f0-16e4c13fa3a7
dθ = [change_in_act_angles[:Θr], change_in_act_angles[:Θz], change_in_act_angles[:Θϕ]]

# ╔═╡ a4ead8dc-c983-4a1e-a6fc-c5dc6d77adbf
let
	fig = Figure(size=(6*72, 5*72))

	local ax
	for i in 1:3
		coord = ["x", "y", "z"][i]
		
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"J %$coord"
		)


		scatter!(orbit_nbody.times, J_nbody[i, :] .- J_old[i, :], label="point orbit old")
		
		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end
		xlims!(t_encounter - time_range, t_encounter + time_range)
		vlines!(orbit_nbody.times[[idx_before, idx_after]])
		hlines!(0, color=:black)
	end

	ylims!(-30, 30)
	linkyaxes!(fig.content...)

	Legend(fig[2, 2], ax)
	fig
end


# ╔═╡ e5fef40f-fa04-4156-9e4b-56fcf935ebba
let
	fig = Figure(size=(6*72, 5*72))

	local ax
	for i in 1:3
		coord = ["x", "y", "z"][i]
		
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"\theta  %$coord"
		)


		scatter!(orbit_nbody.times, Theta_nbody[i, :] .- Theta_old[i, :], label="point orbit old")
		
		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end
		xlims!(t_encounter - time_range, t_encounter + time_range)
		vlines!(orbit_nbody.times[[idx_before, idx_after]])
		hlines!(0, color=:black)

	end

	ylims!(-0.3, 0.3)
	linkyaxes!(fig.content...)

	Legend(fig[2, 2], ax)
	fig
end


# ╔═╡ e0f1db72-8af1-4312-885e-85f74ce04ce7
md"""
# Saving output
"""

# ╔═╡ 0aea9fe2-65ef-42f3-bf58-fa3b81dbba8f
@bind write_ic PlutoUI.Button()

# ╔═╡ e26737e1-39a4-4633-9801-5fb39d652044
md"""
# Reverse actions
"""

# ╔═╡ b0e530c1-ae4b-4c9d-bbee-55d651baf2fb
md"""
Unfortunantly the action mapper class does not always work.
As a result, we just solve for the new position/velocities resulting
in the desired actions.
"""

# ╔═╡ 623661ac-9790-41f0-b9c4-926b50a7f834
import Optim

# ╔═╡ 25368b02-885d-4123-afe4-6766232c9412
af = Agama.ActionFinder(pot_static)

# ╔═╡ a62fb28e-3836-4215-a752-db81948273f1
act_i_point, ang_i_point = get_actions(pot_static, pos_i_point, vel_i_point)

# ╔═╡ 60d591bc-c59c-4a7b-9602-543aac15aef0
act_i_new = act_i_point .+ dJ

# ╔═╡ 4fdc6b72-653f-48db-bb3e-607914fd43f1
ang_i_new = ang_i_point .+ dθ 

# ╔═╡ 7ee29efa-81e7-499d-93ac-11929946939b
function f_to_min(x)
	act, ang, _ = Agama.actions_angles(af, x[1:3], x[4:6], units)
	if any(isnan.(act))
		return Inf
	end

	
	return radii(act, act_i_new) ./ radii(act_i_new) + radii(ang, ang_i_new) 
end

# ╔═╡ 4f235023-42c1-401c-b2eb-3d42952935dc
idx_guess = argmin(f_to_min.(eachcol(vcat(orbit_point.positions, orbit_point.velocities))))

# ╔═╡ d4ccdeb2-055c-4956-9157-0d89fb41c566
minimized = Optim.optimize(f_to_min, vcat(pos_i_point, vel_i_point), Optim.LBFGS(), 
						  Optim.Options(g_tol=1e-3, iterations=1_000))

# ╔═╡ 66c82cb5-d010-418c-96f0-6b11a78546f4
pos_new, vel_new = minimized.minimizer[1:3], minimized.minimizer[4:6]

# ╔═╡ 7e640a85-597a-4852-ad2e-582321199d38
orbit_new_reverse = LilGuys.agama_orbit(pot, LilGuys.Galactocentric(pos_new, vel_new*V2KMS), agama_units=units, timerange=(orbit_nbody.times[idx_before], orbit_nbody.times[1]))

# ╔═╡ b3d27aba-cd7d-4ffc-93ea-e3676e9c2faa
orbit_new = LilGuys.agama_orbit(pot, 
	LilGuys.Galactocentric(orbit_new_reverse.positions[:, end], orbit_new_reverse.velocities[:, end]*V2KMS),
	timerange = (orbit_nbody.times[1], orbit_nbody.times[end]),
	agama_units = units
	)

# ╔═╡ 3e8e7fcf-0c42-4ac3-a76e-dfab44d101db
LilGuys.plot_xyz(LilGuys.positions.((orbit_nbody, orbit_old, orbit_point, orbit_new))..., labels=["nbody", "desired", "point", "new"])

# ╔═╡ 9ea9076a-2d51-41ca-ac10-c71a7e379c79
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "time", ylabel = "galcen radius")

	lines!(orbit_nbody.times, radii(orbit_nbody), label="nbody")
	lines!(orbit_old.times, radii(orbit_old), label="desired")
	lines!(orbit_old.times, radii(orbit_point), label="point")
	lines!(orbit_new.times, radii(orbit_new), label="new")

	vlines!(orbit_nbody.times[[idx_before, idx_after]])
	ylims!(0, nothing)
	Legend(fig[1,2], ax)

	fig
end

# ╔═╡ 14a9600c-2495-4dea-b475-7d7f2ec032c1
[orbit_old, orbit_point, orbit_new]

# ╔═╡ 04de7674-3013-446d-9889-c0749c007f64
J_new, Theta_new = get_actions(pot_static, orbit_new)

# ╔═╡ 108911f0-ca4e-42e4-bf44-c350bebec552
orbit_new.times[argmin(radii(orbit_new))]

# ╔═╡ 876ca278-e9e7-44ad-bbdf-b3f93b300959
let
	fig = Figure(size=(6*72, 5*72))

	local ax
	for i in 1:3
		coord = ["R", "z", "phi"][i]
		
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"J %$coord"
		)


		lines!(orbit_nbody.times, J_nbody[i, :] , label="nbody")
		lines!(orbit_nbody.times, J_old[i, :] , label="point orbit old")

		lines!(orbit_new.times, J_new[i, :], label="point orbit new")
		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

	end
	#ylims!(-1, 1)
	linkyaxes!(fig.content...)

	vlines!(t_encounter)
	Legend(fig[2, 2], ax)
	fig
end


# ╔═╡ 963bd378-2603-45b1-8623-979aad5c2538
df_new = LilGuys.to_frame(orbit_new)

# ╔═╡ ee0a9cd2-5174-42c6-be6b-4273f886be2e
act_new, ang_new = get_actions(pot_static, orbit_new)

# ╔═╡ 4fda883d-7173-49c0-a56e-6572ebc80ba4
let write_ic
	@info "writing  to $modeldir/next_orbit.*"
	
	CSV.write("$modeldir/next_orbit.csv", df_new)
	open("$modeldir/next_orbit.toml", "w") do f
		d = OrderedDict(
			"lastmodel" => orbitname,
			"x_i" => orbit_new.positions[1, 1],
			"y_i" => orbit_new.positions[2, 1],
			"z_i" => orbit_new.positions[3, 1],
			"v_x_i" => orbit_new.velocities[1, 1] * V2KMS,
			"v_y_i" => orbit_new.velocities[2, 1] * V2KMS,
			"v_z_i" => orbit_new.velocities[3, 1] * V2KMS,
			"dJ_R" => dJ[1],
			"dJ_z" => dJ[2],
			"dJ_phi" => dJ[3],
			"dtheta_R" => dθ[1],
			"dtheta_z" => dθ[2],
			"dtheta_phi" => dθ[3],
			"J_R" => act_new[1, 1],
			"J_z" => act_new[2, 1],
			"J_phi" => act_new[3, 1],
			"theta_R" => ang_new[1, 1],
			"theta_z" => ang_new[2, 1],
			"theta_phi" => ang_new[3, 1],
			"time_adjust" => orbit_nbody.times[idx_before],
			"time_measure" => orbit_nbody.times[idx_after],
		)

		for k in ["ra", "dec", "pmra", "pmdec", "distance", "radial_velocity"]
			d[k] = obs_props[k]
			if k * "_err" ∈ keys(obs_props)
				d[k*"_err"] = obs_props[k * "_err"]
			elseif k * "_em" ∈ keys(obs_props)
				d[k*"_em"] = obs_props[k * "_em"]
				d[k*"_ep"] = obs_props[k * "_ep"]
			else
				@warn("error associated with $k missing")
			end
		end
		
		TOML.print(f, d)
		d
	end
end

# ╔═╡ d2eecec1-e01c-4427-874e-c9cbbfa8e30a
act_i_recovered, ang_i_recovered = get_actions(pot_static, minimized.minimizer[1:3], minimized.minimizer[4:6])

# ╔═╡ 983ddeaa-1654-4c4b-9798-4d9b77372a7b
act_i_new - act_i_recovered, ang_i_new - ang_i_recovered

# ╔═╡ 5fdb6129-0496-4f55-a22c-ff5a7ec5352e
act_i_new, ang_i_new

# ╔═╡ 33b006bd-d00a-4f08-833a-524970ed50f3
act_i_point, ang_i_point

# ╔═╡ 5a63da10-6f23-452a-b820-acec63c259dd
act_i_new - act_i_point, ang_i_new - ang_i_point

# ╔═╡ 3d9e7313-3060-46d3-8a7d-c72f062bc827
act_i_recovered - act_i_point, ang_i_recovered - ang_i_point

# ╔═╡ 72970a0d-1347-42ea-907e-cc17692f2b7c
md"""
# Additional points
"""

# ╔═╡ 706a0753-018a-48f7-8c77-23af747141fd
@savefig "delta_xyz_time" let
	fig = Figure(size=(6*72, 5*72))

	local ax
	for i in 1:3
		coord = ["x", "y", "z"][i]
		
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"\Delta %$coord"
		)

		lines!(times, 0*orbit_nbody.positions[i, :], label="nbody")

		lines!(orbit_old_resampled.times, orbit_old_resampled.positions[i, :] .- orbit_nbody.positions[i, :], label="point orbit old")
		
	
		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

		vlines!(orbit_nbody.times[[idx_before, idx_after]])

	end

	Legend(fig[2, 2], ax)
	fig
end


# ╔═╡ d69c4ec1-c596-4f9e-97a6-1b9fbb0e6feb
idx_future = idx_after #+ 20

# ╔═╡ d8e81768-ed46-4012-be13-194dde6ccdaa
ic_future = Galactocentric(orbit_nbody.positions[:, idx_future], orbit_nbody.velocities[:, idx_future] .* V2KMS)

# ╔═╡ 4ed2282e-a16f-4c3a-83b1-61898045fac9
orbit_nbody.times

# ╔═╡ 6238247a-e9ff-44d5-9e83-7f57d056a2e1
orbit_future = LilGuys.agama_orbit(pot,  ic_future, agama_units=units, timerange=(orbit_nbody.times[idx_future]+0, orbit_nbody.times[end]))

# ╔═╡ fd9769e7-9a3e-411a-b9af-fafaf8bf04b9
LilGuys.plot_xyz(LilGuys.positions.((orbit_nbody, orbit_future, orbit_old))..., labels=["nbody", "future", "point"])

# ╔═╡ ed6a519b-f7d4-4ff1-859b-d704667287e8
orbit_future_resampled = LilGuys.resample(orbit_future, orbit_nbody.times)

# ╔═╡ 35b06869-b449-4193-b3b2-cff7bfc06083
J_future, Theta_future = get_actions(pot_static, orbit_future_resampled)

# ╔═╡ b1439b38-6dad-487c-b78b-32736c8aa560
@savefig "xyz_time" let
	fig = Figure(size=(6*72, 5*72))

	local ax
	for i in 1:3
		coord = ["x", "y", "z"][i]
		
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"%$coord"
		)

		lines!(times, orbit_nbody.positions[i, :], label="nbody")

		lines!(orbit_old.times, orbit_old.positions[i, :], label="point orbit old")
		if @isdefined orbit_future
			lines!(orbit_future.times, orbit_future.positions[i, :], label="point orbit future")
		end

		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

	end

	Legend(fig[2, 2], ax)
	fig
end


# ╔═╡ d3b3f0ec-bf69-47c5-ba45-8291acad8c5b
@savefig "delta_xyz_time" let
	fig = Figure(size=(6*72, 5*72))

	local ax
	for i in 1:3
		coord = ["x", "y", "z"][i]
		
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"\Delta %$coord"
		)

		lines!(times, 0*orbit_nbody.positions[i, :], label="nbody")
		
		lines!(orbit_future.times, orbit_future.positions[i, :] .- LilGuys.resample(orbit_nbody, orbit_future.times).positions[i, :], label="point orbit future")

		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

		vlines!(orbit_nbody.times[[idx_before, idx_after, idx_future]])

	end

	linkaxes!(fig.content...)

	Legend(fig[2, 2], ax)
	fig
end


# ╔═╡ fcac8470-c0a7-477e-84f8-3a294dbacd6d
md"""
# Energy and Angular Momentum
This is a nice double check on the method
"""

# ╔═╡ d64c3ce0-7a58-4bea-9159-f63be1141954
E_old = Φ_in(orbit_old.positions, orbit_old.times) .+ 1/2 * speeds(orbit_old) .^2

# ╔═╡ 071df488-b34a-4636-9677-930c8a0f9a56
E_point = Φ_in(orbit_point.positions, orbit_point.times) .+ 1/2 * speeds(orbit_point) .^2

# ╔═╡ 8f7b9a24-d612-4680-9d85-336d93e2662d
E_new = Φ_in(orbit_new.positions, orbit_new.times) .+ 1/2 * radii(orbit_new.velocities) .^2

# ╔═╡ b1f28df3-000d-4f4e-bc5e-32a11116395b
E_new[1] .- E_old[1]

# ╔═╡ ff040d49-f2e7-4042-9931-e94773e1e09e
E_nbody = Φ_in(orbit_nbody.positions, orbit_nbody.times) .+ 1/2 * radii(orbit_nbody.velocities) .^2

# ╔═╡ f18c492b-e67a-4e5f-9544-b7cff2172abc
E_nbody_f = LilGuys.mean(E_nbody[end-window:end])

# ╔═╡ 75f560d3-c29b-4557-b368-ada4ff2611b8
dE_suggested = (E_obs) .- E_nbody_f

# ╔═╡ 302cb462-8413-4cfc-ad75-dfd44b681694
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = "energy"
	)
	hlines!(E_obs)
	
	lines!(times, E_nbody)
	lines!(orbit_old.times, E_old)
	lines!(orbit_point.times, E_point)
	if @isdefined orbit_new
		lines!(orbit_new.times, E_new )
	end



	fig
end

# ╔═╡ 4d60bb0a-34bf-4146-8185-b636e492c425
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = "energy"
	)
	
	lines!(times, E_nbody .- E_old)
	hlines!(0, color=:black)
	
		vlines!(orbit_nbody.times[[idx_before, idx_after]])

	fig
end

# ╔═╡ c98ee04d-adb3-4630-bb0c-45fff32913f7
L_new = LilGuys.angular_momenta(orbit_new)

# ╔═╡ 740c8828-f1f8-44ab-8a38-18e49e4ca62f
L_old = LilGuys.angular_momenta(orbit_old)

# ╔═╡ 2f0414b7-e4ba-45cd-97e8-2473dd85ed1f
L_nbody = LilGuys.angular_momenta(orbit_nbody)

# ╔═╡ 81da9ee9-f7b9-47b8-8fef-951103f32d8e
let
	fig = Figure()

	for i in 1:3
		coord = ["x", "y", "z"][i]
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"L_%$coord",
		)

		lines!(times, L_nbody[i, :])
		lines!(orbit_old.times, L_old[i, :])
		if @isdefined orbit_new
			lines!(orbit_new.times, L_new[i, :])
		end

		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

		L_mean = LilGuys.mean(L_nbody[i, end-window:end])
		lines!(times[end-window:end], fill(L_mean, window+1), color=COLORS[4], linestyle=:solid, linewidth=2)

	end

	fig
end

# ╔═╡ d5a79624-ea87-49fd-86c5-09581b8b1200
let
	fig = Figure()

	for i in 1:3
		coord = ["x", "y", "z"][i]
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"\Delta L_%$coord",
				limits=(-2000, 0, -5, 5)
		)

		lines!(times, L_nbody[i, :] .-  L_old[i, :])

		if i < 3
			hidexdecorations!(ax)
		end

		vlines!(orbit_nbody.times[[idx_before, idx_after]])

		hlines!(0, color=:black)
	end

	linkaxes!(fig.content...)
	
	fig
end

# ╔═╡ 2f1441f5-25e0-41e0-b7f2-4298c4ee6bdc
let
	fig = Figure()

	for i in 1:3
		coord = ["x", "y", "z"][i]
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"\Delta L_%$coord",
				limits=(orbit_nbody.times[[idx_before, idx_after]]..., -2, 2)
		)

		lines!(times, L_nbody[i, :] .-  L_old[i, :])

		if i < 3
			hidexdecorations!(ax)
		end

	
	end

	linkaxes!(fig.content...)
	
	fig
end

# ╔═╡ 5b816336-e545-4bc2-abef-27d8dac642ca
let
	fig = Figure()

	ax = Axis(fig[1, 1],
		xlabel = "time",
		ylabel = L"L",
	)

	lines!(times, radii(L_nbody) .- radii(L_old))


	fig
end

# ╔═╡ Cell order:
# ╟─9bfa1465-2ca9-495c-bab6-5110072f02ed
# ╠═5006f134-648a-4318-9dc3-5de383ac4d0e
# ╠═9c404056-2980-4620-9e4f-459157533c77
# ╟─e9c9d1d9-42e6-4732-b425-477ec0bd7589
# ╠═b18d1734-e40e-11ef-1007-314371eb1a54
# ╠═0b6b20b0-bba8-4ed3-b8a5-d8361bc3674f
# ╠═34e54152-cebb-4364-b1b7-37e0dd4f10c8
# ╠═f0f7b2db-b926-46d5-a9bc-baa9cff0b6cf
# ╠═4bc528f9-4ac0-407e-88a6-5e606117a3d3
# ╠═bab15d74-24e2-4bec-a9c7-39163297b360
# ╠═481edcca-51a8-4f11-88ec-84939b223bdf
# ╠═598be49e-82ef-4820-8996-ea2b753f4275
# ╠═fba1564a-c12d-4657-86fe-d4d771be1037
# ╠═0b6b2020-282b-4384-ac58-8665364cf86e
# ╠═06fbb3db-26ae-4380-935c-883a02c6e119
# ╠═c8e2c46d-917c-4fa8-a841-ea26430d17d3
# ╠═fa2ccfd9-59ee-4484-88ae-997afa83ec86
# ╟─b5302ffe-76e1-424d-9e08-1656981b068e
# ╟─5c8ac63e-3661-455b-a434-7415981b9d89
# ╠═de8bdfef-efca-4d2a-a9cd-a79ccc931b3f
# ╠═afeb4045-01be-42e1-9998-36df9b197bb9
# ╠═7bcdef79-7391-4f35-902b-62d2dda1d858
# ╟─f9c28bef-7beb-42c7-a708-c32dc1d43924
# ╠═087ae8e4-47ea-4c16-8f70-b38168834268
# ╠═c1acdfae-e69d-4c62-9bc7-1a0e1f0f6492
# ╠═7ada86e5-bc6e-4157-b9cd-a7cf154d25d0
# ╠═7fe3e7a2-3ad9-4f37-8a9c-0a50c8aa0c01
# ╟─acf9bb4c-c03f-455f-bd12-65766f55bfaa
# ╠═f7659682-7f06-4618-aa70-4f2c7d3ee838
# ╠═64d8c003-8b5d-4eb4-adb3-ccfae58f08d0
# ╠═51a6bc2e-7e1c-4dad-8a88-40015dc3a1f6
# ╠═0706ee1c-74a6-4a5d-907f-8884895b0361
# ╟─7131d149-c2b6-4ea4-a022-954102bf0870
# ╠═3e8e7fcf-0c42-4ac3-a76e-dfab44d101db
# ╠═fd9769e7-9a3e-411a-b9af-fafaf8bf04b9
# ╠═9ea9076a-2d51-41ca-ac10-c71a7e379c79
# ╠═32e71749-6fd2-4640-838f-ece1ee354c6b
# ╟─5cd10512-5472-4cd5-b9ef-2b5eb0f4bb0c
# ╠═c7aa2829-9048-42b0-9348-3e36714e8e37
# ╠═6a22bbb7-f070-4cc5-9328-9e21cac65487
# ╠═94323325-0d5c-43b8-8793-e28ebad685d9
# ╠═14a9600c-2495-4dea-b475-7d7f2ec032c1
# ╟─d7ca78c0-3898-4bd8-b28f-115644bcf314
# ╠═ec74d9e2-b2ae-447d-b9c6-5d55dd464d28
# ╠═fcb9b299-d660-48ba-afec-93c52a40d930
# ╠═4d2d9e31-54eb-45cb-b228-9375a829c5a9
# ╠═5ad6cf47-3ea5-41e1-81f0-16e4c13fa3a7
# ╠═7e640a85-597a-4852-ad2e-582321199d38
# ╠═b3d27aba-cd7d-4ffc-93ea-e3676e9c2faa
# ╠═cc180e47-3dca-4209-9ac6-1aa4031faca4
# ╠═d83bda34-f614-459d-8fb5-283a11b371f1
# ╠═1e1f37a8-c3d0-4f78-b542-61b8e0814a9c
# ╠═ae7444ec-b934-4b98-966d-1b55760cc592
# ╠═3d33f8f9-16a4-47ac-be22-d35f09a3ee7e
# ╠═8b36e2c4-9073-4af7-a9f4-7ea391e12e56
# ╠═35b06869-b449-4193-b3b2-cff7bfc06083
# ╠═ed6a519b-f7d4-4ff1-859b-d704667287e8
# ╠═65afad37-79b3-408d-87ca-492bd55d4eca
# ╠═09d3bbb1-31f6-42b0-9362-4a1defe4958a
# ╠═5c4c0767-6c17-43e4-a9cf-823fde8a4228
# ╠═04de7674-3013-446d-9889-c0749c007f64
# ╠═01060309-2889-4055-8d53-b127cb88d3ea
# ╟─fa257945-41dd-481b-b719-ffffe7b87d67
# ╠═108911f0-ca4e-42e4-bf44-c350bebec552
# ╠═cc8d1ced-99a0-41af-8fe8-f4c5bd8eef91
# ╟─14667772-b5cd-4d79-9278-1451bf1bdeaa
# ╠═876ca278-e9e7-44ad-bbdf-b3f93b300959
# ╠═81ff5445-40fb-41e2-b4f8-d0bba8e702a4
# ╠═1e4c0c1e-6809-4083-aec1-d3d3a2374068
# ╠═4dff152a-4426-4e47-addc-4ea3bf9157f0
# ╠═e57d04f4-3908-4316-ad61-407bd1909b9b
# ╠═77b485d5-da49-437f-8117-52a57f276090
# ╠═71c78bfc-5ba6-437f-b6dc-462db6840b08
# ╠═a4ead8dc-c983-4a1e-a6fc-c5dc6d77adbf
# ╠═e5fef40f-fa04-4156-9e4b-56fcf935ebba
# ╟─e0f1db72-8af1-4312-885e-85f74ce04ce7
# ╠═e13b1c88-75c4-498a-87c5-6d11c88fa569
# ╠═0aea9fe2-65ef-42f3-bf58-fa3b81dbba8f
# ╠═963bd378-2603-45b1-8623-979aad5c2538
# ╠═ee0a9cd2-5174-42c6-be6b-4273f886be2e
# ╠═4fda883d-7173-49c0-a56e-6572ebc80ba4
# ╟─e26737e1-39a4-4633-9801-5fb39d652044
# ╟─b0e530c1-ae4b-4c9d-bbee-55d651baf2fb
# ╠═623661ac-9790-41f0-b9c4-926b50a7f834
# ╠═25368b02-885d-4123-afe4-6766232c9412
# ╠═7ee29efa-81e7-499d-93ac-11929946939b
# ╠═a62fb28e-3836-4215-a752-db81948273f1
# ╠═60d591bc-c59c-4a7b-9602-543aac15aef0
# ╠═4fdc6b72-653f-48db-bb3e-607914fd43f1
# ╠═4f235023-42c1-401c-b2eb-3d42952935dc
# ╠═d4ccdeb2-055c-4956-9157-0d89fb41c566
# ╠═66c82cb5-d010-418c-96f0-6b11a78546f4
# ╠═d2eecec1-e01c-4427-874e-c9cbbfa8e30a
# ╠═983ddeaa-1654-4c4b-9798-4d9b77372a7b
# ╠═5fdb6129-0496-4f55-a22c-ff5a7ec5352e
# ╠═33b006bd-d00a-4f08-833a-524970ed50f3
# ╠═5a63da10-6f23-452a-b820-acec63c259dd
# ╠═3d9e7313-3060-46d3-8a7d-c72f062bc827
# ╠═72970a0d-1347-42ea-907e-cc17692f2b7c
# ╠═b1439b38-6dad-487c-b78b-32736c8aa560
# ╠═706a0753-018a-48f7-8c77-23af747141fd
# ╠═d69c4ec1-c596-4f9e-97a6-1b9fbb0e6feb
# ╠═d8e81768-ed46-4012-be13-194dde6ccdaa
# ╠═4ed2282e-a16f-4c3a-83b1-61898045fac9
# ╠═6238247a-e9ff-44d5-9e83-7f57d056a2e1
# ╠═d3b3f0ec-bf69-47c5-ba45-8291acad8c5b
# ╠═fcac8470-c0a7-477e-84f8-3a294dbacd6d
# ╠═d64c3ce0-7a58-4bea-9159-f63be1141954
# ╠═071df488-b34a-4636-9677-930c8a0f9a56
# ╠═8f7b9a24-d612-4680-9d85-336d93e2662d
# ╠═75f560d3-c29b-4557-b368-ada4ff2611b8
# ╠═b1f28df3-000d-4f4e-bc5e-32a11116395b
# ╠═ff040d49-f2e7-4042-9931-e94773e1e09e
# ╠═f18c492b-e67a-4e5f-9544-b7cff2172abc
# ╠═302cb462-8413-4cfc-ad75-dfd44b681694
# ╠═4d60bb0a-34bf-4146-8185-b636e492c425
# ╠═c98ee04d-adb3-4630-bb0c-45fff32913f7
# ╠═740c8828-f1f8-44ab-8a38-18e49e4ca62f
# ╠═2f0414b7-e4ba-45cd-97e8-2473dd85ed1f
# ╠═81da9ee9-f7b9-47b8-8fef-951103f32d8e
# ╠═d5a79624-ea87-49fd-86c5-09581b8b1200
# ╠═2f1441f5-25e0-41e0-b7f2-4298c4ee6bdc
# ╠═5b816336-e545-4bc2-abef-27d8dac642ca
