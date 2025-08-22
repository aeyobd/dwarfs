### A Pluto.jl notebook ###
# v0.20.15

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

# ╔═╡ b35b362e-f275-4d70-a703-f398a08b7a47
using Measurements

# ╔═╡ a47e8c62-2b01-4be4-b91a-20b08f592160
using PythonCall

# ╔═╡ 4bc528f9-4ac0-407e-88a6-5e606117a3d3
using CSV, DataFrames

# ╔═╡ 0b6b20b0-bba8-4ed3-b8a5-d8361bc3674f
using PlutoUI

# ╔═╡ 4ad5c3b3-d4c4-4257-bac0-39930ddae499
using HDF5

# ╔═╡ e13b1c88-75c4-498a-87c5-6d11c88fa569
using OrderedCollections

# ╔═╡ 481edcca-51a8-4f11-88ec-84939b223bdf
include(ENV["DWARFS_ROOT"] * "/utils/agama_utils.jl")

# ╔═╡ 9bfa1465-2ca9-495c-bab6-5110072f02ed
md"""
# Adjusting Orbits

The goal of this notebook is to update the initial position and velocity of an n-body simulation to better agree with the desired present-day kinematics. There are two different frameworks for facilitating this transformation. One is to use angular momentum (which a shift in initial angular momentum should be approximantly conserved) and a second is to use actions. 

The angular momentum framework is potential independent and is more physically understandable. However, there is complexity in defining a useful frame of reference and especially in adjusting positions to better line up with the observed positions. 
"""

# ╔═╡ 0c21b10f-ad86-4894-8959-721742d2a2c1
md"""
I prefer the actions/angles framwork. For a test-particle orbit, the actions are conserved and action angles linearly increase with time. Thus, the deviations from the test-particle orbit in the N-body simulation should (to first order) be easily corrected by equivialnt shifts in the action (angle) initial conditions. The challenge with this framework is action (angles) are more nontrival quantities to calculate and depend on our knowledge of the potential. I use `agama` to compute these quantities, and in a fixed (azimuthal) potential, this framework should be a good approximation. 
"""

# ╔═╡ 9c404056-2980-4620-9e4f-459157533c77
units = Agama.VASILIEV_UNITS

# ╔═╡ 9dad1bb5-4d04-4dd8-88a9-335f87688044
window = 10

# ╔═╡ 7372c7f3-3ea0-413a-81af-6aa93a04bd62
pos_err_scale = 10

# ╔═╡ 056cf4aa-cb61-4cce-9429-22366e9f7e97
vel_err_scale = 2

# ╔═╡ e9c9d1d9-42e6-4732-b425-477ec0bd7589
md"""
# Setup
"""

# ╔═╡ 34e54152-cebb-4364-b1b7-37e0dd4f10c8
import LinearAlgebra: ⋅, ×

# ╔═╡ f0f7b2db-b926-46d5-a9bc-baa9cff0b6cf
import TOML

# ╔═╡ bab15d74-24e2-4bec-a9c7-39163297b360
CairoMakie.activate!(type=:png)

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
	galaxyname = TextField(default="ursa_minor"),
	haloname = TextField(default="1e6_new_v38_r4.0"),
	orbitname = TextField(default="orbit_"),
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

# ╔═╡ e279753d-6be6-4e47-a72a-dc0983acef5e
idx_f = TOML.parsefile(modeldir * "orbital_properties.toml")["idx_f"]

# ╔═╡ 9cbf9f1e-1a61-48ae-9fd8-d4253a98f98b
out = Output(modeldir)

# ╔═╡ de8bdfef-efca-4d2a-a9cd-a79ccc931b3f
pot = Agama.Potential(file=joinpath(modeldir  * inputs.potname))

# ╔═╡ 7bcdef79-7391-4f35-902b-62d2dda1d858
Φ_in(x) = Agama.potential(pot, x)

# ╔═╡ 087ae8e4-47ea-4c16-8f70-b38168834268
orbit_old = Orbit(modeldir * "simulation/orbit.csv")

# ╔═╡ 17c59968-27bf-4ab6-ac27-26db28bdbcbc
h5open(joinpath(modeldir, "centres.hdf5")) do centres
	global x_cen, v_cen, x_cen_err, v_cen_err, times
	
	x_cen = centres["positions"][:, :] 
	v_cen = centres["velocities"][:, :]
	x_cen_err = centres["position_errs"][:] .* pos_err_scale
	v_cen_err = centres["velocity_errs"][:] .* vel_err_scale
	times = centres["times"][:]
end

# ╔═╡ acf9bb4c-c03f-455f-bd12-65766f55bfaa
md"""
# Observational reference
"""

# ╔═╡ f7659682-7f06-4618-aa70-4f2c7d3ee838
obs_props_file = modeldir * "simulation/orbit.toml"

# ╔═╡ eee089cf-fa34-454d-afc6-b195cbb80ed9
r_h = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxyname/observed_properties.toml"))["r_h"] / 60

# ╔═╡ 5906da72-48ce-44f6-b4d7-2738c2df64b5
icrs = LilGuys.coord_from_file(obs_props_file)

# ╔═╡ 1cc54a8f-44df-4f35-aa01-3244bdeae52e
obs_props = let
	df = TOML.parsefile(obs_props_file)
	dm = LilGuys.kpc2dm(df["distance"] ± df["distance_err"])
	df["distance_modulus"] = Measurements.value(dm)
	df["distance_modulus_err"] = Measurements.uncertainty(dm)

	df["ra_err"] = 0
	df["dec_err"] = 0
	df
end

# ╔═╡ 18b29878-b8e3-44e5-8843-ccde9c6bb273
function plot_meas!(meas)
	hlines!(Measurements.value(meas), color=:black, linestyle=:solid)

	y = Measurements.value(meas)
	ye = Measurements.uncertainty(meas)

	hspan!(y-ye, y+ye, color=:black, alpha=0.1)

end

# ╔═╡ 9344b6b4-0b06-4638-bfe0-1c5b3d024d3c
md"""
# New initial conditions
"""

# ╔═╡ ec74d9e2-b2ae-447d-b9c6-5d55dd464d28
pos_i_old = LilGuys.positions(orbit_old)[:, 1]

# ╔═╡ fcb9b299-d660-48ba-afec-93c52a40d930
vel_i_old = LilGuys.velocities(orbit_old)[:, 1]

# ╔═╡ e30a0824-bf8f-4f1f-a8f8-9b9f3bb48a3e
print(pos_i_old,", ", vel_i_old)

# ╔═╡ 297aef42-be47-480d-9a07-0d7c7c6897e8
Agama.length_scale(units) * Agama.velocity_scale(units)

# ╔═╡ cb35d50e-bee5-4aa3-946f-f98954fb0c4e
act_0, ang_0 = Agama.agama[].ActionFinder(pot._py)(vcat([16.13, 92.47, 39.63] ./ 5, [-2.49, -42.78, 86.10] ./ 5 / V2KMS / Agama.velocity_scale(units)), angles=true)

# ╔═╡ 40f9a908-0b90-4657-9e21-11acb2afb93b
ang_0

# ╔═╡ e830b46c-5848-4cfa-8212-16ef88d5a6ed
Agama.agama[].ActionMapper(pot._py, tol=1e-15)(vcat(py2vec(act_0), py2vec(ang_0)))

# ╔═╡ 9c9a1dee-881d-4e7d-b9b3-321afde406da
L_i_old = pos_i_old × vel_i_old

# ╔═╡ a2c949eb-21e2-4647-a53f-b245b92ea2a7
md"""
### Validation check
All of the below vectors should be very similar. This is a validation test of the from_actions method.
"""

# ╔═╡ 19803b8e-df6a-41cd-a2dd-dc4f04f95a46
x_cen[:, 1], v_cen[:, 1]

# ╔═╡ d84f975b-dcb6-4d6f-aaa9-497036cf7edf
pos_i_old, vel_i_old

# ╔═╡ 5bfc3f07-853e-45ea-984a-07a76c372a18
md"""
# Actions, Orbits, and validation
"""

# ╔═╡ 95149181-b6c1-429e-be8e-f46d1a2616ef
function get_actions(pot, pos, vel; units=units)
	af = Agama.ActionFinder(pot)

	return Agama.actions_angles(af, pos, vel, units)[[1, 2]]
end

# ╔═╡ 2294328a-67a5-4d49-bae5-dee3c71e17bf
function get_actions(af, pos, vel, pos_err, vel_err; units=units, N=1000)

	N = size(pos, 2)
	actions = Matrix{Float64}(undef, 3, N)
	angles = Matrix{Float64}(undef, 3, N)

	actions_err = Matrix{Float64}(undef, 3, N)
	angles_err = Matrix{Float64}(undef, 3, N)

	for i in 1:size(pos, 2)
		pos_samples = pos[:, i] .+ pos_err[i] * randn(3, N)
		vel_samples = vel[:, i] .+ vel_err[i] * randn(3, N)
		act, ang, freq_py = Agama.actions_angles(af, pos_samples, vel_samples, units)
		
		act_m = LilGuys.mean.(eachrow(act))
		act_err = LilGuys.std.(eachrow(act))
		ang_m = LilGuys.mean.(eachrow(ang))
		ang_err = LilGuys.std.(eachrow(ang))
		actions[:, i] = act_m
		actions_err[:, i] = act_err
		angles[:, i] = ang_m
		angles_err[:, i] = ang_err

	end

	return actions, angles, actions_err, angles_err
end

# ╔═╡ 95255ae4-dc4c-43da-9f19-5fb8d9817e3f
get_actions(pot, pos_i_old, vel_i_old)

# ╔═╡ 4445f7b7-20a5-4a01-b5b2-ff82aa7929c4
get_actions(pot, pos_i_old, vel_i_old)

# ╔═╡ 3ab820ad-1e95-4533-aa70-ad56b11b549a
am = Agama.ActionMapper(pot)

# ╔═╡ 6b2fa1a5-6031-4c37-92ac-6a6088570dba
"""
	shift_actions(Potential, pos_i, vel_i; dJ, dθ)

Shifts the actions and angles of the initial conditions by the specified delta vectors
"""
function shift_actions(Φ, pos_i, vel_i; dJ, dθ)
	act, ang = get_actions(Φ, pos_i, vel_i)

	return Agama.from_actions(am, act .+ dJ, ang .+ dθ, units)
end

# ╔═╡ b3fd0308-2056-4d16-8ff2-bb4d17aaa131
vec.(shift_actions(pot, pos_i_old, vel_i_old, dJ=zeros(3), dθ=zeros(3)))

# ╔═╡ 68669a9d-625a-4b23-99ef-0f235f64454e
Agama.from_actions(am, get_actions(pot, pos_i_old, vel_i_old)...)

# ╔═╡ 373cff92-1e12-4832-9427-47bc3922ba99
Agama.from_actions(am, [2, 0, 0], [3.14, 1.85, 4.97], units)

# ╔═╡ 0cc89288-a915-49a5-a7f3-10abbc02b184
af = Agama.ActionFinder(pot)

# ╔═╡ eec1b26c-d028-4305-b2ff-f5772868dfed
md"""
## Observed properties
"""

# ╔═╡ d2ab7faf-27ed-43f3-b915-374b244149d3
md"""
To properly estimate the uncertainties on the actions, we just MC sample over the uncertainties on the observed properties (rv, distances, pms, etc.).
Then we take the mean and standard deviaton for the action distributions
"""

# ╔═╡ 8af92209-eafb-4da6-b641-a9c3c2ea1080
df_coord_samples = let
	coords = LilGuys.rand_coords(obs_props, 3000)
	gcs = LilGuys.transform.(Galactocentric, coords)
	pos = hcat(LilGuys.position.(gcs)...)
	vel = hcat(LilGuys.velocity.(gcs)...) ./ V2KMS
	
	acts, angs = get_actions(pot, pos, vel)

	df_coord_samples = LilGuys.to_frame(coords)
	df_coord_samples[!, :J_R] = acts[1, :]
	df_coord_samples[!, :J_z] = acts[2, :]
	df_coord_samples[!, :J_phi] = acts[3, :]

	df_coord_samples[!, :Theta_R] = angs[1, :]
	df_coord_samples[!, :Theta_z] = angs[2, :]
	df_coord_samples[!, :Theta_phi] = angs[3, :]


	L = LilGuys.angular_momenta(pos, vel)
	E = Φ_in(pos) .+ 1/2 * radii(vel) .^ 2


	df_coord_samples[:, :L_x] = L[1, :]
	df_coord_samples[:, :L_y] = L[2, :]
	df_coord_samples[:, :L_z] = L[3, :]

	df_coord_samples[:, :E] = E
	df_coord_samples
end

# ╔═╡ 0706ee1c-74a6-4a5d-907f-8884895b0361
E_obs = LilGuys.mean(df_coord_samples.E) ± LilGuys.std(df_coord_samples.E)

# ╔═╡ bbfbf56b-eb98-45f1-bafc-4b04b1f0eda1
L_obs = [LilGuys.mean(x) ± LilGuys.std(x) for x in eachcol(df_coord_samples[:, [:L_x, :L_y, :L_z]])]

# ╔═╡ c133bc01-0493-4ccf-b584-6f34722e9c91
J_obs = [LilGuys.mean(x) ± LilGuys.std(x) for x in eachcol(df_coord_samples[:, [:J_R, :J_z, :J_phi]])]

# ╔═╡ 6d0d7acb-cf73-4ad8-ad74-1b84b343bc10
hist(df_coord_samples.J_phi)

# ╔═╡ 7131d149-c2b6-4ea4-a022-954102bf0870
md"""
## For the orbits
"""

# ╔═╡ 8f3f6ff6-eaf9-4ee4-8a66-0edbe7785162
J_exp, Theta_exp = get_actions(pot, orbit_old.positions, orbit_old.velocities)

# ╔═╡ 95359526-d9e0-4b6f-8f65-6394ef433e35
Theta_obs = [LilGuys.mean(x) ± LilGuys.std(x) for x in eachcol(df_coord_samples[:, [:Theta_R, :Theta_z, :Theta_phi]])]

# ╔═╡ 8303b321-1b95-4c7c-8c4f-d4c6f4e1b344
act_nbody, ang_nbody, act_err_nbody, ang_nbody_err = get_actions(af, x_cen, v_cen, x_cen_err, v_cen_err)

# ╔═╡ a1caf3e7-5b10-4872-a980-5b4fa8448f65
md"""
Note: actions are more important to adjust on initial iterations. Angle evolution changes depending on the actions, so this is best adjusted at the last step.
"""

# ╔═╡ 571ae542-406c-4e62-9b8f-af34e777748f
dθ_suggested = [Theta_obs[i] .- ang_nbody[i, idx_f] for i in 1:3]

# ╔═╡ 8ebf47fd-acd6-48d1-95df-d83a602e75f4
J_f_mean = [LilGuys.mean((act_nbody .± act_err_nbody)[i, idx_f-window:idx_f]) for i in 1:3]

# ╔═╡ b266eafd-debd-4942-a210-71102ba172ff
dJ_suggested =  Measurements.value.(J_obs) .- J_f_mean

# ╔═╡ ae8f7bbd-46ee-4fcb-900a-557e5bb11088
@bind change_in_act_angles confirm(notebook_inputs(;
	Jr = NumberField(-3:0.01:3, default=dJ_suggested[1]),
	Jz = NumberField(-3:0.01:3, default=dJ_suggested[2]),
	Jϕ = NumberField(-3:0.01:3, default=dJ_suggested[3]),
	Θr = NumberField(-3.15:0.01:3.15, default=dθ_suggested[1]),
	Θz = NumberField(-3.15:0.01:3.15, default=dθ_suggested[2]),
	Θϕ = NumberField(-3.15:0.01:3.15, default=dθ_suggested[3]),
))

# ╔═╡ 9da81ee9-5638-48dc-949e-603874b3f627
dJ = [change_in_act_angles[:Jr], change_in_act_angles[:Jz], change_in_act_angles[:Jϕ]]

# ╔═╡ 0fef057f-2f68-447d-9926-8704536cc7e5
dθ = [change_in_act_angles[:Θr], change_in_act_angles[:Θz], change_in_act_angles[:Θϕ]]

# ╔═╡ 9facd961-6c18-4b4a-9e84-9a4db31600c1
pos_new, vel_new = vec.(shift_actions(pot, pos_i_old, vel_i_old, dJ=dJ, dθ=dθ))

# ╔═╡ 80556944-3d5f-4e4b-b3e3-68b092d006cd
gc_new = LilGuys.Galactocentric(pos_new, vel_new*V2KMS)

# ╔═╡ d5fda85b-45b8-442b-b604-4d45236301d4
orbit = let
	orbit_a = Agama.orbit(pot, pos_new, vel_new, units, timerange=(0, 10 / T2GYR), )

	orbit = Orbit(positions=orbit_a.positions, velocities=orbit_a.velocities, times=orbit_a.times)
	orbit
end

# ╔═╡ bb0195a0-d06d-4f72-8646-31ef848e987a
x_new = orbit.positions

# ╔═╡ af07bdd5-4baa-481d-a973-a4c119af6c5f
v_new = orbit.velocities

# ╔═╡ 3fa86c92-e6bd-41b6-bd64-cfd33f74b229
act_new, ang_new = get_actions(pot, x_new, v_new)

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

		lines!(out.times .- out.times[idx_f], x_cen[i, :], label="nbody")

		lines!(orbit_old.times .- orbit_old.times[end], orbit_old.positions[i, :], label="point orbit old")
		if @isdefined orbit
			lines!(orbit.times .- orbit.times[end], orbit.positions[i, :], label="point orbit new")
		end
		plot_meas!(J_obs[i])
		
		
		#lines!(out.times[[idx_f-window, idx_f]] .- out.times[idx_f], Measurements.value.([J_f_mean[i], J_f_mean[i]]), linestyle=:solid, linewidth=2, label="adopted mean")
		
		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

	end

	Legend(fig[2, 2], ax)
	fig
end


# ╔═╡ b39cbb71-b98f-4a6b-ba01-55d5b2bb2190
get_actions(pot, pos_new, vel_new)

# ╔═╡ b8f2475a-3588-4611-99ac-de156f25b853
@savefig "actions_adjustment" let
	fig = Figure(size=(6*72, 5*72))

	local ax
	for i in 1:3
		coord = ["r", "z", "ϕ"][i]
		
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"J_%$coord"
		)

		lines!(out.times .- out.times[idx_f], act_nbody[i, :], label="nbody")
		band!(out.times .- out.times[idx_f], act_nbody[i, :] .- act_err_nbody[i, :], act_nbody[i, :] .+ act_err_nbody[i, :], alpha=0.5)

		lines!(orbit_old.times .- orbit_old.times[end], J_exp[i, :], label="point orbit old")
		if @isdefined orbit
			lines!(orbit.times .- orbit.times[end], act_new[i, :], label="point orbit new")
		end
		plot_meas!(J_obs[i])
		
		
		lines!(out.times[[idx_f-window, idx_f]] .- out.times[idx_f], Measurements.value.([J_f_mean[i], J_f_mean[i]]), linestyle=:solid, linewidth=2, label="adopted mean")
		
		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

	end

	Legend(fig[2, 2], ax)
	fig
end


# ╔═╡ d5be7108-a6c6-4e2a-98b9-562dbd0dc999
ang_nbody_err[:, idx_f]

# ╔═╡ 706a0753-018a-48f7-8c77-23af747141fd
@savefig "action_angles_adjustment" let
	fig = Figure(size=(5*72, 5*72))
	local ax

	for i in 1:3
		coord = ["r", "z", "ϕ"][i]

		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"\theta_%$coord",
		)

			plot_meas!(Theta_obs[i])

		idx =idx_f-5:idx_f+0
		lines!(out.times[idx] .- out.times[idx_f], ang_nbody[i, idx], label="nbody")

		
		lines!(orbit_old.times .- orbit_old.times[end], Theta_exp[i, :], label="point orbit old")
		if @isdefined orbit
			lines!(orbit.times .- orbit.times[end], ang_new[i, :], label="point orbit new")
		end

		xlims!(-50, 10)
		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end
	
	end

	Legend(fig[2,2], ax)

	fig
end


# ╔═╡ e7aba352-84a1-417c-b25b-7d0503030e6a
idx_f

# ╔═╡ 2af93853-0df1-45a5-ace6-75ad532e49dd
@savefig "action_angles_evolution"  let
	fig = Figure(size=(5*72, 5*72))

	local ax
	for i in 1:3
		coord = ["r", "z", "ϕ"][i]

		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"\theta_%$coord",
		)
		
		lines!(out.times , ang_nbody[i, :], label="nbody")
		lines!(orbit_old.times, Measurements.value.(Theta_exp[i, :]), label="initial point")
		lines!(orbit.times, ang_new[i, :], label="final point")

		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

	end
	Legend(fig[2, 2], ax)

	fig
end


# ╔═╡ a419c394-d597-4f6f-812f-cd0a0bc7ea0e
LilGuys.plot_xyz(orbit_old.positions, x_new, labels=["old point orbit", "new point orbit"])

# ╔═╡ 9ea9076a-2d51-41ca-ac10-c71a7e379c79
let
	fig = Figure(size=(6, 4) .* Arya.UNITS_PER_INCH)
	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = L"r",
	)
	lines!(out.times, radii(x_cen), label = "nbody")

	lines!(orbit_old.times, radii(orbit_old), label="initial point")
	if @isdefined orbit
		lines!(orbit.times, radii(orbit), label="next test")
	end

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 5ced8710-8187-482c-b00c-29f21f6084d5
md"""
# Energy and Angular Momentum
This is a nice double check on the method
"""

# ╔═╡ 4473d1b0-fea7-4d87-8837-1966124f1cf0
E_old = Φ_in(orbit_old.positions) .+ 1/2 * speeds(orbit_old) .^2

# ╔═╡ e2f1e41d-4996-4f81-8dda-35cf6ba0fbba
E_new = Φ_in(x_new) .+ 1/2 * radii(v_new) .^2

# ╔═╡ b4247610-8250-4e28-bb30-79dac76cc1bb
E_new[1] .- E_old[1]

# ╔═╡ c2b025c6-0bd9-49be-9848-3f79e31e6140
E_nbody = Φ_in(x_cen) .+ 1/2 * radii(v_cen) .^2

# ╔═╡ 588295b0-b629-4526-8069-58a17ad9a32b
E_nbody_f = LilGuys.mean(E_nbody[idx_f-window:idx_f])

# ╔═╡ 7c3d6c08-a345-4418-91a4-a14e401c50ea
dE_suggested = Measurements.value(E_obs) .- E_nbody_f

# ╔═╡ da3f27fc-0f2c-4e03-bc36-ccbb5eee1ea2
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = "energy"
	)
	plot_meas!(E_obs)
	
	lines!(out.times, E_nbody)
	lines!(orbit_old.times, E_old)
	if @isdefined orbit
		lines!(orbit.times, E_new )
	end


	lines!(out.times[idx_f-window:idx_f], fill(E_nbody_f, window+1), color=:black, linestyle=:solid, linewidth=2)

	fig
end


# ╔═╡ f8e73044-8d07-4b09-81ee-54f7f7278496
L_new = LilGuys.angular_momenta(x_new, v_new)

# ╔═╡ 471fae8a-4998-4662-87e9-d4277f5f8604
L_old = LilGuys.angular_momenta(orbit_old)

# ╔═╡ 74a73f14-426b-45b6-afdd-55bc4a919b34
L_nbody = LilGuys.angular_momenta(x_cen, v_cen)

# ╔═╡ 0f936517-1924-44db-ba03-2c7fd589afa3
let
	fig = Figure()

	for i in 1:3
		coord = ["x", "y", "z"][i]
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"L_%$coord",
		)

		lines!(out.times, L_nbody[i, :])
		lines!(orbit_old.times, L_old[i, :])
		if @isdefined orbit
			lines!(orbit.times, L_new[i, :])
		end
		plot_meas!(L_obs[i])

		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

		L_mean = LilGuys.mean(L_nbody[i, idx_f-window:idx_f])
		lines!(out.times[idx_f-window:idx_f], fill(L_mean, window+1), color=COLORS[4], linestyle=:solid, linewidth=2)

	end

	fig
end


# ╔═╡ e0f1db72-8af1-4312-885e-85f74ce04ce7
md"""
# Saving output
"""

# ╔═╡ 0aea9fe2-65ef-42f3-bf58-fa3b81dbba8f
@bind write_ic PlutoUI.Button()

# ╔═╡ 963bd378-2603-45b1-8623-979aad5c2538
df_new = LilGuys.to_frame(orbit)

# ╔═╡ 4fda883d-7173-49c0-a56e-6572ebc80ba4
let write_ic
	@info "writing  to $modeldir/next_orbit.*"
	
	CSV.write("$modeldir/next_orbit.csv", df_new)
	open("$modeldir/next_orbit.toml", "w") do f
		d = OrderedDict(
			"lastmodel" => orbitname,
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
			"x_i" => orbit.positions[1, 1],
			"y_i" => orbit.positions[2, 1],
			"z_i" => orbit.positions[3, 1],
			"v_x_i" => orbit.velocities[1, 1] * V2KMS,
			"v_y_i" => orbit.velocities[2, 1] * V2KMS,
			"v_z_i" => orbit.velocities[3, 1] * V2KMS,
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

# ╔═╡ Cell order:
# ╟─9bfa1465-2ca9-495c-bab6-5110072f02ed
# ╟─0c21b10f-ad86-4894-8959-721742d2a2c1
# ╠═5006f134-648a-4318-9dc3-5de383ac4d0e
# ╠═9c404056-2980-4620-9e4f-459157533c77
# ╠═9dad1bb5-4d04-4dd8-88a9-335f87688044
# ╠═7372c7f3-3ea0-413a-81af-6aa93a04bd62
# ╠═056cf4aa-cb61-4cce-9429-22366e9f7e97
# ╟─e9c9d1d9-42e6-4732-b425-477ec0bd7589
# ╠═b18d1734-e40e-11ef-1007-314371eb1a54
# ╠═34e54152-cebb-4364-b1b7-37e0dd4f10c8
# ╠═b35b362e-f275-4d70-a703-f398a08b7a47
# ╠═a47e8c62-2b01-4be4-b91a-20b08f592160
# ╠═f0f7b2db-b926-46d5-a9bc-baa9cff0b6cf
# ╠═4bc528f9-4ac0-407e-88a6-5e606117a3d3
# ╠═0b6b20b0-bba8-4ed3-b8a5-d8361bc3674f
# ╠═4ad5c3b3-d4c4-4257-bac0-39930ddae499
# ╠═bab15d74-24e2-4bec-a9c7-39163297b360
# ╠═481edcca-51a8-4f11-88ec-84939b223bdf
# ╠═598be49e-82ef-4820-8996-ea2b753f4275
# ╠═fba1564a-c12d-4657-86fe-d4d771be1037
# ╠═0b6b2020-282b-4384-ac58-8665364cf86e
# ╠═06fbb3db-26ae-4380-935c-883a02c6e119
# ╠═c8e2c46d-917c-4fa8-a841-ea26430d17d3
# ╠═fa2ccfd9-59ee-4484-88ae-997afa83ec86
# ╟─b5302ffe-76e1-424d-9e08-1656981b068e
# ╠═e279753d-6be6-4e47-a72a-dc0983acef5e
# ╠═9cbf9f1e-1a61-48ae-9fd8-d4253a98f98b
# ╠═de8bdfef-efca-4d2a-a9cd-a79ccc931b3f
# ╠═7bcdef79-7391-4f35-902b-62d2dda1d858
# ╠═087ae8e4-47ea-4c16-8f70-b38168834268
# ╠═17c59968-27bf-4ab6-ac27-26db28bdbcbc
# ╟─acf9bb4c-c03f-455f-bd12-65766f55bfaa
# ╠═f7659682-7f06-4618-aa70-4f2c7d3ee838
# ╠═95255ae4-dc4c-43da-9f19-5fb8d9817e3f
# ╠═e30a0824-bf8f-4f1f-a8f8-9b9f3bb48a3e
# ╠═eee089cf-fa34-454d-afc6-b195cbb80ed9
# ╠═5906da72-48ce-44f6-b4d7-2738c2df64b5
# ╠═1cc54a8f-44df-4f35-aa01-3244bdeae52e
# ╠═18b29878-b8e3-44e5-8843-ccde9c6bb273
# ╟─9344b6b4-0b06-4638-bfe0-1c5b3d024d3c
# ╠═6b2fa1a5-6031-4c37-92ac-6a6088570dba
# ╠═ec74d9e2-b2ae-447d-b9c6-5d55dd464d28
# ╠═fcb9b299-d660-48ba-afec-93c52a40d930
# ╠═68669a9d-625a-4b23-99ef-0f235f64454e
# ╠═4445f7b7-20a5-4a01-b5b2-ff82aa7929c4
# ╠═373cff92-1e12-4832-9427-47bc3922ba99
# ╠═297aef42-be47-480d-9a07-0d7c7c6897e8
# ╠═cb35d50e-bee5-4aa3-946f-f98954fb0c4e
# ╠═40f9a908-0b90-4657-9e21-11acb2afb93b
# ╠═e830b46c-5848-4cfa-8212-16ef88d5a6ed
# ╠═9c9a1dee-881d-4e7d-b9b3-321afde406da
# ╠═9facd961-6c18-4b4a-9e84-9a4db31600c1
# ╟─a2c949eb-21e2-4647-a53f-b245b92ea2a7
# ╠═19803b8e-df6a-41cd-a2dd-dc4f04f95a46
# ╠═d84f975b-dcb6-4d6f-aaa9-497036cf7edf
# ╠═b3fd0308-2056-4d16-8ff2-bb4d17aaa131
# ╟─5bfc3f07-853e-45ea-984a-07a76c372a18
# ╠═95149181-b6c1-429e-be8e-f46d1a2616ef
# ╠═2294328a-67a5-4d49-bae5-dee3c71e17bf
# ╠═3ab820ad-1e95-4533-aa70-ad56b11b549a
# ╠═0cc89288-a915-49a5-a7f3-10abbc02b184
# ╟─eec1b26c-d028-4305-b2ff-f5772868dfed
# ╠═d2ab7faf-27ed-43f3-b915-374b244149d3
# ╠═8af92209-eafb-4da6-b641-a9c3c2ea1080
# ╠═0706ee1c-74a6-4a5d-907f-8884895b0361
# ╠═bbfbf56b-eb98-45f1-bafc-4b04b1f0eda1
# ╠═c133bc01-0493-4ccf-b584-6f34722e9c91
# ╠═6d0d7acb-cf73-4ad8-ad74-1b84b343bc10
# ╟─7131d149-c2b6-4ea4-a022-954102bf0870
# ╠═8f3f6ff6-eaf9-4ee4-8a66-0edbe7785162
# ╠═95359526-d9e0-4b6f-8f65-6394ef433e35
# ╠═80556944-3d5f-4e4b-b3e3-68b092d006cd
# ╠═d5fda85b-45b8-442b-b604-4d45236301d4
# ╠═bb0195a0-d06d-4f72-8646-31ef848e987a
# ╠═af07bdd5-4baa-481d-a973-a4c119af6c5f
# ╠═8303b321-1b95-4c7c-8c4f-d4c6f4e1b344
# ╠═3fa86c92-e6bd-41b6-bd64-cfd33f74b229
# ╠═b39cbb71-b98f-4a6b-ba01-55d5b2bb2190
# ╠═9da81ee9-5638-48dc-949e-603874b3f627
# ╠═0fef057f-2f68-447d-9926-8704536cc7e5
# ╟─a1caf3e7-5b10-4872-a980-5b4fa8448f65
# ╠═b266eafd-debd-4942-a210-71102ba172ff
# ╠═571ae542-406c-4e62-9b8f-af34e777748f
# ╠═ae8f7bbd-46ee-4fcb-900a-557e5bb11088
# ╠═b8f2475a-3588-4611-99ac-de156f25b853
# ╠═b1439b38-6dad-487c-b78b-32736c8aa560
# ╠═8ebf47fd-acd6-48d1-95df-d83a602e75f4
# ╠═d5be7108-a6c6-4e2a-98b9-562dbd0dc999
# ╠═706a0753-018a-48f7-8c77-23af747141fd
# ╠═e7aba352-84a1-417c-b25b-7d0503030e6a
# ╠═2af93853-0df1-45a5-ace6-75ad532e49dd
# ╠═a419c394-d597-4f6f-812f-cd0a0bc7ea0e
# ╠═9ea9076a-2d51-41ca-ac10-c71a7e379c79
# ╟─5ced8710-8187-482c-b00c-29f21f6084d5
# ╠═4473d1b0-fea7-4d87-8837-1966124f1cf0
# ╠═e2f1e41d-4996-4f81-8dda-35cf6ba0fbba
# ╠═7c3d6c08-a345-4418-91a4-a14e401c50ea
# ╠═b4247610-8250-4e28-bb30-79dac76cc1bb
# ╠═c2b025c6-0bd9-49be-9848-3f79e31e6140
# ╠═588295b0-b629-4526-8069-58a17ad9a32b
# ╠═da3f27fc-0f2c-4e03-bc36-ccbb5eee1ea2
# ╠═f8e73044-8d07-4b09-81ee-54f7f7278496
# ╠═471fae8a-4998-4662-87e9-d4277f5f8604
# ╠═74a73f14-426b-45b6-afdd-55bc4a919b34
# ╠═0f936517-1924-44db-ba03-2c7fd589afa3
# ╟─e0f1db72-8af1-4312-885e-85f74ce04ce7
# ╠═e13b1c88-75c4-498a-87c5-6d11c88fa569
# ╠═0aea9fe2-65ef-42f3-bf58-fa3b81dbba8f
# ╠═963bd378-2603-45b1-8623-979aad5c2538
# ╠═4fda883d-7173-49c0-a56e-6572ebc80ba4
