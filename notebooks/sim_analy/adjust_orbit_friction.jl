### A Pluto.jl notebook ###
# v0.20.17

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
	galaxyname = TextField(default="sculptor"),
	haloname = TextField(default="1e6_new_v31_r3.2"),
	orbitname = TextField(default="L3M11_9Gyr_smallperi"),
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

# ╔═╡ afeb4045-01be-42e1-9998-36df9b197bb9
pot_static = Agama.Potential(file=joinpath(modeldir  * "/potential_mw_init.ini"))

# ╔═╡ 7bcdef79-7391-4f35-902b-62d2dda1d858
Φ_in(x, t=0) = Agama.potential(pot, x,units, t=t)

# ╔═╡ 087ae8e4-47ea-4c16-8f70-b38168834268
orbit_old = LilGuys.resample(Orbit(modeldir * "simulation/orbit.csv"), out.times)

# ╔═╡ 7ada86e5-bc6e-4157-b9cd-a7cf154d25d0
orbit_nbody = Orbit(joinpath(modeldir, "centres.hdf5"))

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

	df_coord_samples = LilGuys.to_frame(coords)

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

# ╔═╡ 7131d149-c2b6-4ea4-a022-954102bf0870
md"""
## For the orbits
"""

# ╔═╡ 80556944-3d5f-4e4b-b3e3-68b092d006cd
gc_new = LilGuys.Galactocentric(orbit_old.positions[:, end], orbit_old.velocities[:, end]*V2KMS)

# ╔═╡ a1caf3e7-5b10-4872-a980-5b4fa8448f65
md"""
Note: actions are more important to adjust on initial iterations. Angle evolution changes depending on the actions, so this is best adjusted at the last step.
"""

# ╔═╡ 706a0753-018a-48f7-8c77-23af747141fd


# ╔═╡ 5ced8710-8187-482c-b00c-29f21f6084d5
md"""
# Energy and Angular Momentum
This is a nice double check on the method
"""

# ╔═╡ 4473d1b0-fea7-4d87-8837-1966124f1cf0
E_old = Φ_in(orbit_old.positions, orbit_old.times) .+ 1/2 * speeds(orbit_old) .^2

# ╔═╡ c2b025c6-0bd9-49be-9848-3f79e31e6140
E_nbody = Φ_in(x_cen, times) .+ 1/2 * radii(v_cen) .^2

# ╔═╡ 588295b0-b629-4526-8069-58a17ad9a32b
E_nbody_f = LilGuys.mean(E_nbody[idx_f-window:idx_f])

# ╔═╡ 7c3d6c08-a345-4418-91a4-a14e401c50ea
dE_suggested = Measurements.value(E_obs) .- E_nbody_f

# ╔═╡ b0c77794-b94f-4daf-9631-529d6fd9e2e3
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = "energy"
	)
	
	lines!(out.times, E_nbody .- E_old)
	

	fig
end


# ╔═╡ 471fae8a-4998-4662-87e9-d4277f5f8604
L_old = LilGuys.angular_momenta(orbit_old)

# ╔═╡ 74a73f14-426b-45b6-afdd-55bc4a919b34
L_nbody = LilGuys.angular_momenta(x_cen, v_cen)

# ╔═╡ 1a90fe1c-c0cd-40c0-a7cd-27cce5034810
let
	fig = Figure()

	for i in 1:3
		coord = ["x", "y", "z"][i]
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"\Delta L_%$coord",
				limits=(-2000, 1000, -5, 5)
		)

		lines!(out.times, L_nbody[i, :] .-  L_old[i, :])

		if i < 3
			hidexdecorations!(ax)
		end

	
	end

	linkaxes!(fig.content...)
	
	fig
end


# ╔═╡ 65bcb37f-32a6-4e01-9454-0bb815cd93f9
let
	fig = Figure()

	ax = Axis(fig[1, 1],
		xlabel = "time",
		ylabel = L"L",
	)

	lines!(out.times, radii(L_nbody) .- radii(L_old))


	# plot_meas!(L_obs[i])

	# 	L_mean = LilGuys.mean(L_nbody[i, idx_f-window:idx_f])
	# 	lines!(out.times[idx_f-window:idx_f], fill(L_mean, window+1), color=COLORS[4], linestyle=:solid, linewidth=2)

	# end

	fig
end


# ╔═╡ e0f1db72-8af1-4312-885e-85f74ce04ce7
md"""
# Saving output
"""

# ╔═╡ 0aea9fe2-65ef-42f3-bf58-fa3b81dbba8f
@bind write_ic PlutoUI.Button()

# ╔═╡ 45dbf4df-edc2-4705-a796-2528accc1449
md"""
# Dynamical friction model
"""

# ╔═╡ 914a5864-7225-4384-8941-d497bac408d2
function get_sigma_v(pot)
	gm = Agama.GalaxyModel(pot, Agama.DistributionFunction(pot, type="QuasiSpherical"))


	N = 100
	Rs = 10 .^ LinRange(-1, 3, N)

	sigmas = gm._py.moments(Agama.mat2py([zeros(N) Rs zeros(N)]'), dens=false) |> Agama.py2mat

	sigma3 = sigmas[1, :] .+ sigmas[2, :] + sigmas[3, :]

	sigma1 = sigma3 ./ 3
	return LilGuys.lerp(Rs, sqrt.(sigma1))
end

# ╔═╡ 6a3113d9-80c3-4e3e-b0ee-40f164d97dac
σv = get_sigma_v(pot_static)

# ╔═╡ 4424837b-9749-4ac5-a856-8b7dbf717622


# ╔═╡ 264c3eb9-4d72-4a2e-a0bf-dab32dc8454e
import Optim

# ╔═╡ bb1c934e-9c27-483c-b73d-4fbd71fc0cfa
function make_orbit_with_fric( coord=Galactocentric(orbit_old.positions[:, 1], orbit_old.velocities[:, 1] .* V2KMS), 
							   ; timerange=(orbit_old.times[1], 0), r_s=3, M=0.2, σv=r->150/V2KMS, logΛ=nothing, kwargs...)

	if isnothing(logΛ)
		Λ = nothing
	else
		Λ = exp(logΛ)
	end
	
	dyn_fric = LilGuys.ChandrashakarDynamicalFriction(;
		r_s=r_s, σv=σv, M=M, ρ = x->Agama.density(pot_static, x, units), Λ=Λ, kwargs...)

	f_int = (pos, vel, t) -> (LilGuys.acceleration(dyn_fric, pos, vel) + Agama.acceleration(pot, pos, units,t=t))


	orbit_dyn_fric = LilGuys.leapfrog(f_int, coord, timerange=timerange) 

end

# ╔═╡ 39140e52-67b7-4004-a24d-9bdd633a42d5
σv(5)

# ╔═╡ 6aa12dc1-da84-47f2-b0ed-1cd3e9617b89


# ╔═╡ 4df3341b-754e-420e-9486-59dac6131509
obs_props["t_i"]

# ╔═╡ 2afa98f0-9dce-4571-b558-d5a6920271f0
function f_diff_final(params; kwargs...)
	M, logΛ = params
	o = make_orbit_with_fric(; M=M, logΛ=logΛ, σv=150, kwargs...)
	return radii(x_cen[:, end], o.positions[:, end])
end

# ╔═╡ e77b1014-81bc-485c-835c-d74153469e8a
function f_diff_final_free(params; kwargs...)
	M, r_s = params
	o = make_orbit_with_fric(; M=M, r_s=r_s, kwargs...)
	return radii(x_cen[:, end], o.positions[:, end])
end

# ╔═╡ a4b67dc3-4709-4359-9168-d37d7c845526
f_diff_final([0.11,4.0])

# ╔═╡ 271edaeb-3300-4102-8fba-59ac11decc59
soln = Optim.optimize(f_diff_final, [0.1, 3] )

# ╔═╡ 7bd0b3fb-04f4-480c-8a08-05219376cdbc
soln.minimizer

# ╔═╡ 7ae87fe8-c3e2-4bc2-ac91-a48fc5d9e0eb
soln2 = Optim.optimize(f_diff_final_free, [0.01, 0.1], [10, 10], [0.1, 2],  Optim.SAMIN(), Optim.Options(iterations=10^2))

# ╔═╡ 748d47ca-0f98-4d55-af19-ee5f9fc97c5b
soln2.minimizer

# ╔═╡ 2c911eed-4c1b-419d-ab05-aa6db135562e
orbit_old.times

# ╔═╡ 3d86a0bf-c3c5-436d-a4f8-16116140240b
@bind orbital_params confirm(notebook_inputs(;
	r = NumberField(0:0.01:10, default=10^soln2.minimizer[2]),
	M = NumberField(0:0.01:10, default=soln2.minimizer[1]),
	σv = NumberField(0:1:300, default=150),
))

# ╔═╡ d5fda85b-45b8-442b-b604-4d45236301d4
orbit = make_orbit_with_fric(gc_new, 
	M=orbital_params.M, r_s=orbital_params.r,timerange=(0, orbit_old.times[1]))


# ╔═╡ bb0195a0-d06d-4f72-8646-31ef848e987a
x_new = orbit.positions

# ╔═╡ af07bdd5-4baa-481d-a973-a4c119af6c5f
v_new = orbit.velocities

# ╔═╡ f8e73044-8d07-4b09-81ee-54f7f7278496
L_new = LilGuys.angular_momenta(x_new, v_new)

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

		lines!(orbit_old.times, orbit_old.positions[i, :], label="point orbit old")
		if @isdefined orbit
			lines!(orbit.times, orbit.positions[i, :], label="point orbit new")
		end

		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

	end

	Legend(fig[2, 2], ax)
	fig
end


# ╔═╡ a419c394-d597-4f6f-812f-cd0a0bc7ea0e
LilGuys.plot_xyz(orbit_old.positions, orbit.positions, labels=["old point orbit", "new point orbit"])

# ╔═╡ 9ea9076a-2d51-41ca-ac10-c71a7e379c79
let
	fig = Figure()
	ax = Axis(fig[1,1])


	lines!(orbit_old.times, radii(orbit_old))
	lines!(orbit.times, radii(orbit))
	lines!(times, radii(x_cen))

	fig
end

# ╔═╡ 77387a88-7f36-4453-a3b5-df1277497dbc
orbit.pericenter

# ╔═╡ e2f1e41d-4996-4f81-8dda-35cf6ba0fbba
E_new = Φ_in(x_new, orbit.times) .+ 1/2 * radii(v_new) .^2

# ╔═╡ b4247610-8250-4e28-bb30-79dac76cc1bb
E_new[1] .- E_old[1]

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


# ╔═╡ 8f28e702-b5e0-4972-90f9-0df943d7b2a4
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = "velocity"
	)
	plot_meas!(E_obs)
	
	lines!(out.times, radii(v_cen) .- speeds(orbit_old))
	if @isdefined orbit
		lines!(orbit.times, E_new )
	end


	lines!(out.times[idx_f-window:idx_f], fill(E_nbody_f, window+1), color=:black, linestyle=:solid, linewidth=2)

	fig
end


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


# ╔═╡ 963bd378-2603-45b1-8623-979aad5c2538
df_new = LilGuys.to_frame(orbit)

# ╔═╡ 4fda883d-7173-49c0-a56e-6572ebc80ba4
let write_ic
	@info "writing  to $modeldir/next_orbit.*"
	
	CSV.write("$modeldir/next_orbit.csv", df_new)
	open("$modeldir/next_orbit.toml", "w") do f
		d = OrderedDict(
			"lastmodel" => orbitname,
			"x_i" => orbit.positions[1, 1],
			"y_i" => orbit.positions[2, 1],
			"z_i" => orbit.positions[3, 1],
			"v_x_i" => orbit.velocities[1, 1] * V2KMS,
			"v_y_i" => orbit.velocities[2, 1] * V2KMS,
			"v_z_i" => orbit.velocities[3, 1] * V2KMS,
			"M_sat_dynf" => orbital_params.M, 
			"r_sat_dynf" => orbital_params.r, 
			"sigma_v_dynf" => orbital_params.σv,
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

# ╔═╡ c7cf21c6-ab6c-4720-9bf9-d0d9768e2132
orbit_dyn_fric = make_orbit_with_fric(M=orbital_params.M, r_s=orbital_params.r)

# ╔═╡ aa59896a-80f2-4ec5-9842-b8afe42f0133
orbit_dyn_fric.velocities

# ╔═╡ d63af59f-9c48-4d75-a4b0-34926b30d366
radii(orbit_dyn_fric.positions[:, end], x_cen[:, end])

# ╔═╡ cb1545ef-da60-4a8f-885f-5ed3030a916c
LilGuys.plot_xyz(x_cen, orbit_dyn_fric.positions, orbit_old.positions)

# ╔═╡ 9788fab5-b00d-450b-80dc-9265899a824d
md"""
# Additional acceleration

We can calculate the acceleration due to the galactic potential, recovering the missing friction force not accounted for in the orbit
"""

# ╔═╡ d88e2a01-7ac8-4dbb-b699-f99225948365
orbit_interp = Agama.agama[].Spline(times, x_cen[1, :], v_cen[1, :]), Agama.agama[].Spline(times, x_cen[2, :], v_cen[2, :]), Agama.agama[].Spline(times, x_cen[3, :], v_cen[3, :])

# ╔═╡ 25c9260c-d740-477f-8217-c204592d87df
orbit_interp[1](5)

# ╔═╡ 48d583ce-a9dd-411f-83df-14b8f5e8c5b3
accs = Agama.acceleration(pot, x_cen, units, t=times)

# ╔═╡ 0d0c04c8-a618-4fd6-b331-7e1f649b8f30
function predict_next_position(pos, vel, time, time_next)
	o = LilGuys.agama_orbit(pot, Galactocentric(pos, vel .* V2KMS), timerange=(time, time_next), agama_units=units)
	return o.positions[:, end], o.velocities[:, end]
end

# ╔═╡ 55b64f3b-05c0-48ed-a941-5fe4f6482e2e
function predict_next_positions(x_cen, v_cen, times)
	x_pred = similar(x_cen)
	v_pred = similar(v_cen)


	for i in eachindex(times)[1:end-1]
		x_new, v_new = predict_next_position(x_cen[:, i], v_cen[:, i], times[i], times[i+1])

		x_pred[:, i+1] .= x_new
		v_pred[:, i+1] .= v_new
	end

	x_pred[:, 1] = x_cen[:, 1]
	v_pred[:, 1] = v_cen[:, 1]

	return x_pred, v_pred
end

# ╔═╡ b20f4af4-cee7-415f-9f07-366b0e16bd0e
x_pred, v_pred = predict_next_positions(x_cen, v_cen, times)

# ╔═╡ 607b091e-3069-40d3-a43e-68d15c060535
dx_pred = x_pred .- x_cen

# ╔═╡ e84842b4-2aa6-4c74-8212-fc7ce0604060
x_pred

# ╔═╡ 0874a3a9-dfbc-44b4-94e1-ab60b18b2758
dv_pred = v_cen .- v_pred

# ╔═╡ 3b60257e-089e-4ccb-9123-6c4b10a6a478
acc_res = dv_pred[:, 2:end] ./ diff(times)'

# ╔═╡ b68a0b1b-407a-4834-bef2-f97a9a22c5e1
plot(log10.(radii(x_cen)[2:end]), log10.(radii(acc_res)))

# ╔═╡ ac498cfd-7ce6-4344-be27-55057eb838d1
plot(log10.(radii(v_cen)[2:end]), log10.(radii(acc_res)))

# ╔═╡ b4ed176b-c2ce-4bd0-88ed-7b9ba1b491e0
acc_proj = [acc_res[:, i] ⋅ v_cen[:, i+1] for i in eachindex(times)[1:end-1]]

# ╔═╡ 453622ee-8d14-416d-8b1a-fb90f916838f
plot(log10.(radii(acc_res)), log10.(abs.(acc_proj) ./ radii(v_cen)[2:end] ./ radii(acc_res)))

# ╔═╡ d7ca78c0-3898-4bd8-b28f-115644bcf314
md"""
# Action differences
"""

# ╔═╡ cc180e47-3dca-4209-9ac6-1aa4031faca4
function get_actions(pot, orbit; units=units)
	af = Agama.ActionFinder(pot)

	return Agama.actions_angles(af, orbit.positions, orbit.velocities, units)[[1, 2]]
end

# ╔═╡ f72939d8-87c7-4fb3-9660-c11ac77e389f
am = Agama.ActionMapper(pot_static)

# ╔═╡ 25368b02-885d-4123-afe4-6766232c9412
af = Agama.ActionFinder(pot_static)

# ╔═╡ 3d33f8f9-16a4-47ac-be22-d35f09a3ee7e
time_shift = -(orbit_nbody.times[argmin(radii(orbit_nbody))] - orbit_old.times[argmin(radii(orbit_old))])

# ╔═╡ 1e1f37a8-c3d0-4f78-b542-61b8e0814a9c
orbit_old_resampled = LilGuys.resample(orbit_old, orbit_nbody.times .+ time_shift)

# ╔═╡ 8b36e2c4-9073-4af7-a9f4-7ea391e12e56
orbit_nbody.times[argmin(radii(orbit_nbody))] - orbit_old_resampled.times[argmin(radii(orbit_old_resampled))]

# ╔═╡ 35b06869-b449-4193-b3b2-cff7bfc06083
J_old, Theta_old = get_actions(pot_static, orbit_old_resampled)

# ╔═╡ 5c4c0767-6c17-43e4-a9cf-823fde8a4228
J_nbody, Theta_nbody = get_actions(pot_static, orbit_nbody)

# ╔═╡ 01060309-2889-4055-8d53-b127cb88d3ea
t_encounter = orbit_nbody.times[argmin(radii(orbit_nbody))]

# ╔═╡ 14667772-b5cd-4d79-9278-1451bf1bdeaa
md"""
## Plots
"""

# ╔═╡ 32e71749-6fd2-4640-838f-ece1ee354c6b
let
	fig = Figure(size=(6*72, 5*72))

	local ax
	for i in 1:3
		coord = ["x", "y", "z"][i]
		
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"J %$coord"
		)


		lines!(orbit_nbody.times, J_nbody[i, :] .- J_old[i, :], label="point orbit old")
		
		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

	end
	ylims!(-1, 1)
	linkyaxes!(fig.content...)

	vlines!(t_encounter)
	Legend(fig[2, 2], ax)
	fig
end


# ╔═╡ c7aa2829-9048-42b0-9348-3e36714e8e37
time_range = 100

# ╔═╡ 81ff5445-40fb-41e2-b4f8-d0bba8e702a4
idx_before = argmin(abs.(orbit_nbody.times .- t_encounter .- time_range/2))

# ╔═╡ 1e4c0c1e-6809-4083-aec1-d3d3a2374068
idx_after = argmin(abs.(orbit_nbody.times .- t_encounter .+ time_range/2))

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

	end

	ylims!(-0.2, 0.2)
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
			ylabel = L"theta J %$coord"
		)


		scatter!(orbit_nbody.times, Theta_nbody[i, :] .- Theta_old[i, :], label="point orbit old")
		
		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end
		xlims!(t_encounter - time_range, t_encounter + time_range)

	end

	ylims!(-0.01, 0.01)
	linkyaxes!(fig.content...)

	Legend(fig[2, 2], ax)
	fig
end


# ╔═╡ Cell order:
# ╟─9bfa1465-2ca9-495c-bab6-5110072f02ed
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
# ╠═afeb4045-01be-42e1-9998-36df9b197bb9
# ╠═7bcdef79-7391-4f35-902b-62d2dda1d858
# ╠═087ae8e4-47ea-4c16-8f70-b38168834268
# ╠═7ada86e5-bc6e-4157-b9cd-a7cf154d25d0
# ╠═17c59968-27bf-4ab6-ac27-26db28bdbcbc
# ╟─acf9bb4c-c03f-455f-bd12-65766f55bfaa
# ╠═f7659682-7f06-4618-aa70-4f2c7d3ee838
# ╠═e30a0824-bf8f-4f1f-a8f8-9b9f3bb48a3e
# ╠═eee089cf-fa34-454d-afc6-b195cbb80ed9
# ╠═5906da72-48ce-44f6-b4d7-2738c2df64b5
# ╠═1cc54a8f-44df-4f35-aa01-3244bdeae52e
# ╠═18b29878-b8e3-44e5-8843-ccde9c6bb273
# ╟─9344b6b4-0b06-4638-bfe0-1c5b3d024d3c
# ╠═ec74d9e2-b2ae-447d-b9c6-5d55dd464d28
# ╠═fcb9b299-d660-48ba-afec-93c52a40d930
# ╠═9c9a1dee-881d-4e7d-b9b3-321afde406da
# ╟─a2c949eb-21e2-4647-a53f-b245b92ea2a7
# ╠═19803b8e-df6a-41cd-a2dd-dc4f04f95a46
# ╠═d84f975b-dcb6-4d6f-aaa9-497036cf7edf
# ╟─eec1b26c-d028-4305-b2ff-f5772868dfed
# ╠═d2ab7faf-27ed-43f3-b915-374b244149d3
# ╠═8af92209-eafb-4da6-b641-a9c3c2ea1080
# ╠═0706ee1c-74a6-4a5d-907f-8884895b0361
# ╠═bbfbf56b-eb98-45f1-bafc-4b04b1f0eda1
# ╟─7131d149-c2b6-4ea4-a022-954102bf0870
# ╠═80556944-3d5f-4e4b-b3e3-68b092d006cd
# ╠═d5fda85b-45b8-442b-b604-4d45236301d4
# ╠═bb0195a0-d06d-4f72-8646-31ef848e987a
# ╠═af07bdd5-4baa-481d-a973-a4c119af6c5f
# ╟─a1caf3e7-5b10-4872-a980-5b4fa8448f65
# ╠═b1439b38-6dad-487c-b78b-32736c8aa560
# ╠═706a0753-018a-48f7-8c77-23af747141fd
# ╠═a419c394-d597-4f6f-812f-cd0a0bc7ea0e
# ╠═9ea9076a-2d51-41ca-ac10-c71a7e379c79
# ╠═77387a88-7f36-4453-a3b5-df1277497dbc
# ╟─5ced8710-8187-482c-b00c-29f21f6084d5
# ╠═4473d1b0-fea7-4d87-8837-1966124f1cf0
# ╠═e2f1e41d-4996-4f81-8dda-35cf6ba0fbba
# ╠═7c3d6c08-a345-4418-91a4-a14e401c50ea
# ╠═b4247610-8250-4e28-bb30-79dac76cc1bb
# ╠═c2b025c6-0bd9-49be-9848-3f79e31e6140
# ╠═588295b0-b629-4526-8069-58a17ad9a32b
# ╠═da3f27fc-0f2c-4e03-bc36-ccbb5eee1ea2
# ╠═b0c77794-b94f-4daf-9631-529d6fd9e2e3
# ╠═8f28e702-b5e0-4972-90f9-0df943d7b2a4
# ╠═f8e73044-8d07-4b09-81ee-54f7f7278496
# ╠═471fae8a-4998-4662-87e9-d4277f5f8604
# ╠═74a73f14-426b-45b6-afdd-55bc4a919b34
# ╠═0f936517-1924-44db-ba03-2c7fd589afa3
# ╠═1a90fe1c-c0cd-40c0-a7cd-27cce5034810
# ╠═65bcb37f-32a6-4e01-9454-0bb815cd93f9
# ╟─e0f1db72-8af1-4312-885e-85f74ce04ce7
# ╠═e13b1c88-75c4-498a-87c5-6d11c88fa569
# ╠═0aea9fe2-65ef-42f3-bf58-fa3b81dbba8f
# ╠═963bd378-2603-45b1-8623-979aad5c2538
# ╠═4fda883d-7173-49c0-a56e-6572ebc80ba4
# ╟─45dbf4df-edc2-4705-a796-2528accc1449
# ╠═914a5864-7225-4384-8941-d497bac408d2
# ╠═6a3113d9-80c3-4e3e-b0ee-40f164d97dac
# ╠═4424837b-9749-4ac5-a856-8b7dbf717622
# ╠═264c3eb9-4d72-4a2e-a0bf-dab32dc8454e
# ╠═bb1c934e-9c27-483c-b73d-4fbd71fc0cfa
# ╠═39140e52-67b7-4004-a24d-9bdd633a42d5
# ╠═6aa12dc1-da84-47f2-b0ed-1cd3e9617b89
# ╠═4df3341b-754e-420e-9486-59dac6131509
# ╠═2afa98f0-9dce-4571-b558-d5a6920271f0
# ╠═e77b1014-81bc-485c-835c-d74153469e8a
# ╠═a4b67dc3-4709-4359-9168-d37d7c845526
# ╠═271edaeb-3300-4102-8fba-59ac11decc59
# ╠═7bd0b3fb-04f4-480c-8a08-05219376cdbc
# ╠═7ae87fe8-c3e2-4bc2-ac91-a48fc5d9e0eb
# ╠═748d47ca-0f98-4d55-af19-ee5f9fc97c5b
# ╠═2c911eed-4c1b-419d-ab05-aa6db135562e
# ╠═aa59896a-80f2-4ec5-9842-b8afe42f0133
# ╠═c7cf21c6-ab6c-4720-9bf9-d0d9768e2132
# ╠═d63af59f-9c48-4d75-a4b0-34926b30d366
# ╠═cb1545ef-da60-4a8f-885f-5ed3030a916c
# ╠═3d86a0bf-c3c5-436d-a4f8-16116140240b
# ╟─9788fab5-b00d-450b-80dc-9265899a824d
# ╠═d88e2a01-7ac8-4dbb-b699-f99225948365
# ╠═25c9260c-d740-477f-8217-c204592d87df
# ╠═48d583ce-a9dd-411f-83df-14b8f5e8c5b3
# ╠═0d0c04c8-a618-4fd6-b331-7e1f649b8f30
# ╠═55b64f3b-05c0-48ed-a941-5fe4f6482e2e
# ╠═b20f4af4-cee7-415f-9f07-366b0e16bd0e
# ╠═607b091e-3069-40d3-a43e-68d15c060535
# ╠═e84842b4-2aa6-4c74-8212-fc7ce0604060
# ╠═0874a3a9-dfbc-44b4-94e1-ab60b18b2758
# ╠═3b60257e-089e-4ccb-9123-6c4b10a6a478
# ╠═b68a0b1b-407a-4834-bef2-f97a9a22c5e1
# ╠═ac498cfd-7ce6-4344-be27-55057eb838d1
# ╠═453622ee-8d14-416d-8b1a-fb90f916838f
# ╠═b4ed176b-c2ce-4bd0-88ed-7b9ba1b491e0
# ╠═d7ca78c0-3898-4bd8-b28f-115644bcf314
# ╠═cc180e47-3dca-4209-9ac6-1aa4031faca4
# ╠═f72939d8-87c7-4fb3-9660-c11ac77e389f
# ╠═25368b02-885d-4123-afe4-6766232c9412
# ╠═1e1f37a8-c3d0-4f78-b542-61b8e0814a9c
# ╠═3d33f8f9-16a4-47ac-be22-d35f09a3ee7e
# ╠═8b36e2c4-9073-4af7-a9f4-7ea391e12e56
# ╠═35b06869-b449-4193-b3b2-cff7bfc06083
# ╠═5c4c0767-6c17-43e4-a9cf-823fde8a4228
# ╠═01060309-2889-4055-8d53-b127cb88d3ea
# ╠═14667772-b5cd-4d79-9278-1451bf1bdeaa
# ╠═32e71749-6fd2-4640-838f-ece1ee354c6b
# ╠═81ff5445-40fb-41e2-b4f8-d0bba8e702a4
# ╠═1e4c0c1e-6809-4083-aec1-d3d3a2374068
# ╠═c7aa2829-9048-42b0-9348-3e36714e8e37
# ╠═a4ead8dc-c983-4a1e-a6fc-c5dc6d77adbf
# ╠═e5fef40f-fa04-4156-9e4b-56fcf935ebba
