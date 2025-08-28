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

# ╔═╡ cb4ae829-de51-4ab4-b98c-3e17f0c487d2
using Printf

# ╔═╡ 481edcca-51a8-4f11-88ec-84939b223bdf
include(ENV["DWARFS_ROOT"] * "/utils/agama_utils.jl")

# ╔═╡ 9bfa1465-2ca9-495c-bab6-5110072f02ed
md"""
# Adjusting Orbits

The goal of this notebook is to update the initial position and velocity of an n-body simulation to better agree with the desired present-day kinematics. There are two different frameworks for facilitating this transformation. One is to use angular momentum (which a shift in initial angular momentum should be approximantly conserved) and a second is to use actions. 

The angular momentum framework is potential independent and is more physically understandable. However, there is complexity in defining a useful frame of reference and especially in adjusting positions to better line up with the observed positions. 
"""

# ╔═╡ 0cd25191-b443-4c13-9aa3-5af0c9e9693d
md"""
Similar to @vasiliev2024, we create a reference frame centred at the current orbit to use to adjust the orbits. This is essentially a 6D space.

If we definne the instantaneous ``\vec{L}`` value as the ``z`` axis in a cylendrical reference frame, then we can reparameterize the position and velocity in a more physical manner. However, we furthermore rotate the frame such that the ``\hat \phi'`` direction is along the current orbital velocity. 

Velocity
- ``\Delta v_\phi'`` Change in along-orbit velocity
- ``\Delta v_{R'}``. Change in perpendicular-orbit velocity.
- ``\Delta v_{z}`` Out of plane change in velocity.
Position
- ``\Delta R'`` Current orbital impact parameter shift. 
- ``\Delta z`` Out of plane shift.
- ``\Delta \phi' ``Along orbit shift. Corresponds to a change in the initial time. 

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

# ╔═╡ 12536adf-4683-43d8-920e-03040f115cfe
import LinearAlgebra: normalize

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
	galaxyname = TextField(default="sculptor"),
	haloname = TextField(default="1e6_new_v31_r3.2"),
	orbitname = TextField(default="L3M11_9Gyr_smallperi_example"),
	potname = TextField(default = "simulation/agama_potential.ini")
))

# ╔═╡ ae8f7bbd-46ee-4fcb-900a-557e5bb11088
@bind change_in_positions confirm(notebook_inputs(;
	dx = NumberField(-1:0.01:1, default=0.02),
	dv = NumberField(-1:0.01:1, default=0.05),
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

# ╔═╡ 9cbf9f1e-1a61-48ae-9fd8-d4253a98f98b
out = Output(modeldir)

# ╔═╡ de8bdfef-efca-4d2a-a9cd-a79ccc931b3f
pot = Agama.Potential(file=joinpath(modeldir  * inputs.potname))

# ╔═╡ 7bcdef79-7391-4f35-902b-62d2dda1d858
Φ_in(x, t) = Agama.potential(pot, x,units, t=t)

# ╔═╡ 087ae8e4-47ea-4c16-8f70-b38168834268
orbit_old = LilGuys.resample(Orbit(modeldir * "simulation/orbit.csv"), out.times)

# ╔═╡ 17c59968-27bf-4ab6-ac27-26db28bdbcbc
h5open(joinpath(modeldir, "centres.hdf5")) do centres
	global x_cen, v_cen, x_cen_err, v_cen_err, times
	
	x_cen = centres["positions"][:, :] 
	v_cen = centres["velocities"][:, :]
	x_cen_err = centres["position_errs"][:] .* pos_err_scale
	v_cen_err = centres["velocity_errs"][:] .* vel_err_scale
	times = centres["times"][:]
end

# ╔═╡ e279753d-6be6-4e47-a72a-dc0983acef5e
idx_f = length(times) #85 #TOML.parsefile(modeldir * "orbital_properties.toml")["idx_f"]

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

# ╔═╡ 9c9a1dee-881d-4e7d-b9b3-321afde406da
L_i_old = pos_i_old × vel_i_old

# ╔═╡ 015b8fff-9da2-4449-8aae-8f892f662c16
function orbit_unit_vectors(pos_old, vel_old)
	z_hat = normalize(pos_old × vel_old)
	R_hat = normalize(pos_old)
	phi_hat = R_hat × z_hat
	
	return R_hat, z_hat, phi_hat
end

# ╔═╡ ab9a394e-1a46-44ed-a4ea-984d5a20290a
Rh, zh, phih  = orbit_unit_vectors(pos_i_old, vel_i_old)

# ╔═╡ 7f919011-6346-4bcc-9aba-10d199dd92a3
Rh, pos_i_old

# ╔═╡ 7f8f5e17-9024-4e40-87cf-f467b990894c
function adjust_orbit(pos_old, vel_old, dpos, dvel)
	dR, dz, dphi = dpos
	dvR, dvz, dvphi = dvel

	Rhat, zhat, phihat = orbit_unit_vectors(pos_old, vel_old)

	d_pos = dR * Rhat + dz * zhat + dphi * phihat
	d_vel = dvR * Rhat + dvz * zhat + dvphi * phihat
	return pos_old + d_pos, vel_old + d_vel
end

# ╔═╡ 7131d149-c2b6-4ea4-a022-954102bf0870
md"""
## For the orbits
"""

# ╔═╡ e3e3af5a-c531-4d83-b5a1-6ff28026ba5e
dx = change_in_positions.dx * radii(pos_i_old)

# ╔═╡ bc86bf99-80e9-42ed-a664-b6775609f4db
dv = change_in_positions.dv * radii(vel_i_old)

# ╔═╡ 10ed9395-1b5e-4a8c-b263-2d0a16269d7a
dv * V2KMS

# ╔═╡ 0a6bbb9f-6d4a-4420-9777-85d06f418e6f
function make_orbit(delta_pos, delta_vel)
	pos_new, vel_new = adjust_orbit(pos_i_old, vel_i_old, delta_pos, delta_vel)

	orbit_a = Agama.orbit(pot, pos_new, vel_new, units, timerange=(orbit_old.times[1], orbit_old.times[end]), )

	orbit = Orbit(positions=orbit_a.positions, velocities=orbit_a.velocities, times=orbit_a.times)
	orbit
end

# ╔═╡ 33782762-cc6e-41b4-a1ee-4545bbe03489
orbits = [
	"R" => make_orbit(    [dx, 0, 0], zeros(3)),
	"R_m" => make_orbit(  [-dx, 0, 0], zeros(3)),
	"z" => make_orbit(    [0, dx, 0], zeros(3)),
	"z_m" => make_orbit(  [0, -dx, 0], zeros(3)),
	"phi" => make_orbit(  [0, 0, dx], zeros(3)),
	"phi_m" => make_orbit([0, 0, -dx], zeros(3)),
	"v_R" => make_orbit(zeros(3),     [dv, 0, 0]),
	"v_R_m" => make_orbit(zeros(3),   [-dv, 0, 0]),
	"v_z" => make_orbit(zeros(3),     [0, dv, 0]),
	"v_z_m" => make_orbit(zeros(3),   [0, -dv, 0]),
	"v_phi" => make_orbit(zeros(3),   [0, 0, dv]),
	"v_phi_m" => make_orbit(zeros(3), [0, 0, -dv]),
]

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

		lines!(out.times, x_cen[i, :], label="nbody")

		lines!(orbit_old.times, orbit_old.positions[i, :], label="point orbit old")
		for (label, orbit) in orbits
			lines!(orbit.times, orbit.positions[i, :], label="point orbit $label")
		end
		
		
		#lines!(out.times[[idx_f-window, idx_f]] .- out.times[idx_f], Measurements.value.([J_f_mean[i], J_f_mean[i]]), linestyle=:solid, linewidth=2, label="adopted mean")
		
		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

	end

	Legend(fig[2, 2], ax)
	fig
end


# ╔═╡ a419c394-d597-4f6f-812f-cd0a0bc7ea0e
LilGuys.plot_xyz(LilGuys.positions.(last.(orbits))..., labels=first.(orbits))

# ╔═╡ 9ea9076a-2d51-41ca-ac10-c71a7e379c79
let
	fig = Figure(size=(6, 4) .* Arya.UNITS_PER_INCH)
	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = L"r",
	)
	lines!(out.times, radii(x_cen), label = "nbody")

	lines!(orbit_old.times, radii(orbit_old), label="initial point")
	for (label, orbit) in orbits
		lines!(orbit.times, radii(orbit), label="new $label")
	end

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 86c6ac08-2705-4b63-a711-207e636d504c
let
	fig = Figure(size=(6, 4) .* Arya.UNITS_PER_INCH)
	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = "distance from old orbit",
	)

	for (label, orbit) in orbits
		orbit = LilGuys.resample(orbit, orbit_old.times)
		lines!(orbit.times, radii(orbit.positions, orbit_old.positions), label="new $label")
	end

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 41b934aa-bff4-4a88-84b6-827dff39c2c7
let
	fig = Figure(size=(6, 4) .* Arya.UNITS_PER_INCH)
	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = "distance from mean orbit",
	)

	for (label, orbit) in orbits
		lines!(orbit.times, radii(orbit.positions, orbits[1].second.positions), label="new $label")
	end

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 14cd6367-1e69-4ab7-b744-dc1cbe5c138a
let
	fig = Figure(size=(6, 4) .* Arya.UNITS_PER_INCH)
	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = "velocity distance from mean orbit",
	)

	for (label, orbit) in orbits
		lines!(orbit.times, V2KMS * radii(orbit.velocities, orbits[1].second.velocities), label="new $label")
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
E_old = Φ_in(orbit_old.positions, orbit_old.times) .+ 1/2 * speeds(orbit_old) .^2

# ╔═╡ c2b025c6-0bd9-49be-9848-3f79e31e6140
E_nbody = Φ_in(x_cen, times) .+ 1/2 * radii(v_cen) .^2

# ╔═╡ 7c3d6c08-a345-4418-91a4-a14e401c50ea
dE_suggested = (E_old .- E_nbody)[idx_f]

# ╔═╡ b4247610-8250-4e28-bb30-79dac76cc1bb
(E_old .- E_nbody)[idx_f]

# ╔═╡ 588295b0-b629-4526-8069-58a17ad9a32b
E_nbody_f = LilGuys.mean(E_nbody[idx_f-window:idx_f])

# ╔═╡ da3f27fc-0f2c-4e03-bc36-ccbb5eee1ea2
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = "energy",
	)
	
	lines!(out.times, E_nbody)
	lines!(orbit_old.times, E_old)
	for (label, orbit) in orbits
		 E_new = Φ_in(orbit.positions, orbit.times) .+ 1/2 * radii(orbit.velocities) .^2
		lines!(orbit.times, E_new, label=label )
	end


	lines!(out.times[idx_f-window:idx_f], fill(E_nbody_f, window+1), color=:black, linestyle=:solid, linewidth=2)

		Legend(fig[1, 2], ax)

	
	fig
end


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
		for (label, orbit) in orbits
			L_new = LilGuys.angular_momenta(orbit.positions, orbit.velocities)
			lines!(orbit.times, L_new[i, :])
		end

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

# ╔═╡ 99594bcf-77bc-453b-8145-1669e190f09d
function project(v)
	return [ v ⋅ Rh, v ⋅ zh, v ⋅ phih]
end

# ╔═╡ d49d92db-33c4-4cab-a3c1-f16bfd70ebee
for (label, orbit) in orbits
	x_i = orbit.positions[:, 1]
	v_i = orbit.velocities[:, 1]

	pos_cyl = project(x_i)
	vel_cyl = project(v_i)
	pos_cyl_old = project(pos_i_old)
	vel_cyl_old = project(vel_i_old)
	println(label)
	@printf "x_i\t%8.2f %8.2f %8.2f\n" x_i...
	@printf "x_o\t%8.2f %8.2f %8.2f\n" pos_i_old...
	@printf "v_i\t%8.2f %8.2f %8.2f\n" V2KMS * v_i...
	@printf "v_o\t%8.2f %8.2f %8.2f\n" V2KMS * vel_i_old...
	@printf "Rzϕ = \t%8.2f %8.2f %8.2f\n" pos_cyl...
	@printf "ΔRzϕ = \t%8.2f %8.2f %8.2f\n" (pos_cyl - pos_cyl_old)...
	@printf "Δvzϕ = \t%8.2f %8.2f %8.2f\n" V2KMS*(vel_cyl - vel_cyl_old)...

	println()
end

# ╔═╡ 4fda883d-7173-49c0-a56e-6572ebc80ba4
for (label, orbit) in orbits
	@info "writing  to $modeldir/orbit_$label.*"
	
	CSV.write("$modeldir/orbit_$label.csv", LilGuys.to_frame(orbit))
	open("$modeldir/orbit_$label.toml", "w") do f
		d = OrderedDict(
			"lastmodel" => orbitname,
			"x_i" => orbit.positions[1, 1],
			"y_i" => orbit.positions[2, 1],
			"z_i" => orbit.positions[3, 1],
			"v_x_i" => orbit.velocities[1, 1] * V2KMS,
			"v_y_i" => orbit.velocities[2, 1] * V2KMS,
			"v_z_i" => orbit.velocities[3, 1] * V2KMS,
			"dx" => dx,
			"dv" => dv,
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
# ╟─0cd25191-b443-4c13-9aa3-5af0c9e9693d
# ╠═5006f134-648a-4318-9dc3-5de383ac4d0e
# ╠═ae8f7bbd-46ee-4fcb-900a-557e5bb11088
# ╠═e279753d-6be6-4e47-a72a-dc0983acef5e
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
# ╠═12536adf-4683-43d8-920e-03040f115cfe
# ╠═bab15d74-24e2-4bec-a9c7-39163297b360
# ╠═481edcca-51a8-4f11-88ec-84939b223bdf
# ╠═598be49e-82ef-4820-8996-ea2b753f4275
# ╠═fba1564a-c12d-4657-86fe-d4d771be1037
# ╠═0b6b2020-282b-4384-ac58-8665364cf86e
# ╠═06fbb3db-26ae-4380-935c-883a02c6e119
# ╠═c8e2c46d-917c-4fa8-a841-ea26430d17d3
# ╠═fa2ccfd9-59ee-4484-88ae-997afa83ec86
# ╟─b5302ffe-76e1-424d-9e08-1656981b068e
# ╠═9cbf9f1e-1a61-48ae-9fd8-d4253a98f98b
# ╠═de8bdfef-efca-4d2a-a9cd-a79ccc931b3f
# ╠═7bcdef79-7391-4f35-902b-62d2dda1d858
# ╠═087ae8e4-47ea-4c16-8f70-b38168834268
# ╠═17c59968-27bf-4ab6-ac27-26db28bdbcbc
# ╟─acf9bb4c-c03f-455f-bd12-65766f55bfaa
# ╠═f7659682-7f06-4618-aa70-4f2c7d3ee838
# ╠═e30a0824-bf8f-4f1f-a8f8-9b9f3bb48a3e
# ╠═eee089cf-fa34-454d-afc6-b195cbb80ed9
# ╠═5906da72-48ce-44f6-b4d7-2738c2df64b5
# ╠═1cc54a8f-44df-4f35-aa01-3244bdeae52e
# ╟─9344b6b4-0b06-4638-bfe0-1c5b3d024d3c
# ╠═ec74d9e2-b2ae-447d-b9c6-5d55dd464d28
# ╠═fcb9b299-d660-48ba-afec-93c52a40d930
# ╠═9c9a1dee-881d-4e7d-b9b3-321afde406da
# ╠═ab9a394e-1a46-44ed-a4ea-984d5a20290a
# ╠═7f919011-6346-4bcc-9aba-10d199dd92a3
# ╠═015b8fff-9da2-4449-8aae-8f892f662c16
# ╠═7f8f5e17-9024-4e40-87cf-f467b990894c
# ╟─7131d149-c2b6-4ea4-a022-954102bf0870
# ╠═e3e3af5a-c531-4d83-b5a1-6ff28026ba5e
# ╠═bc86bf99-80e9-42ed-a664-b6775609f4db
# ╠═10ed9395-1b5e-4a8c-b263-2d0a16269d7a
# ╠═33782762-cc6e-41b4-a1ee-4545bbe03489
# ╠═0a6bbb9f-6d4a-4420-9777-85d06f418e6f
# ╠═b1439b38-6dad-487c-b78b-32736c8aa560
# ╠═a419c394-d597-4f6f-812f-cd0a0bc7ea0e
# ╠═9ea9076a-2d51-41ca-ac10-c71a7e379c79
# ╠═86c6ac08-2705-4b63-a711-207e636d504c
# ╠═41b934aa-bff4-4a88-84b6-827dff39c2c7
# ╠═14cd6367-1e69-4ab7-b744-dc1cbe5c138a
# ╟─5ced8710-8187-482c-b00c-29f21f6084d5
# ╠═4473d1b0-fea7-4d87-8837-1966124f1cf0
# ╠═7c3d6c08-a345-4418-91a4-a14e401c50ea
# ╠═b4247610-8250-4e28-bb30-79dac76cc1bb
# ╠═c2b025c6-0bd9-49be-9848-3f79e31e6140
# ╠═588295b0-b629-4526-8069-58a17ad9a32b
# ╠═da3f27fc-0f2c-4e03-bc36-ccbb5eee1ea2
# ╠═471fae8a-4998-4662-87e9-d4277f5f8604
# ╠═74a73f14-426b-45b6-afdd-55bc4a919b34
# ╠═0f936517-1924-44db-ba03-2c7fd589afa3
# ╟─e0f1db72-8af1-4312-885e-85f74ce04ce7
# ╠═e13b1c88-75c4-498a-87c5-6d11c88fa569
# ╠═cb4ae829-de51-4ab4-b98c-3e17f0c487d2
# ╠═99594bcf-77bc-453b-8145-1669e190f09d
# ╠═d49d92db-33c4-4cab-a3c1-f16bfd70ebee
# ╠═4fda883d-7173-49c0-a56e-6572ebc80ba4
