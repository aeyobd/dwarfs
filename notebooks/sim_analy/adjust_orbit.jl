### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
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

# ╔═╡ e9c9d1d9-42e6-4732-b425-477ec0bd7589
md"""
# Setup
"""

# ╔═╡ 34e54152-cebb-4364-b1b7-37e0dd4f10c8
import LinearAlgebra: ⋅, ×

# ╔═╡ f0f7b2db-b926-46d5-a9bc-baa9cff0b6cf
import TOML

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
	haloname = TextField(default="1e6_v37_r5.0"),
	orbitname = TextField(default="orbit_"),
))

# ╔═╡ fba1564a-c12d-4657-86fe-d4d771be1037
galaxyname = inputs.galaxyname

# ╔═╡ 0b6b2020-282b-4384-ac58-8665364cf86e
haloname = inputs.haloname

# ╔═╡ 06fbb3db-26ae-4380-935c-883a02c6e119
orbitname =  inputs.orbitname

# ╔═╡ c8e2c46d-917c-4fa8-a841-ea26430d17d3
modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, orbitname) * "/"

# ╔═╡ b5302ffe-76e1-424d-9e08-1656981b068e
md"""
# Data loading
"""

# ╔═╡ 9cbf9f1e-1a61-48ae-9fd8-d4253a98f98b
out = Output(modeldir)

# ╔═╡ de8bdfef-efca-4d2a-a9cd-a79ccc931b3f
pot = agama.Potential(modeldir * "simulation/agama_potential.ini")

# ╔═╡ e279753d-6be6-4e47-a72a-dc0983acef5e
idx_f = TOML.parsefile(modeldir * "orbital_properties.toml")["idx_f"]

# ╔═╡ fd1adb47-9a55-49c9-b1ca-5009d91b0408
begin 
	orbit_exp = CSV.read(modeldir * "simulation/orbit.csv", DataFrame)
	orbit_exp.t .-= orbit_exp.t[1]
	orbit_exp[!, :R] = orbit_exp.x .⊕ orbit_exp.y
	orbit_exp[!, :v_R] = @. 1/orbit_exp.R *( orbit_exp.v_x * orbit_exp.x + orbit_exp.v_y*orbit_exp.y)
	orbit_exp[!, :v_phi] = @. 1/orbit_exp.R *( orbit_exp.v_x * orbit_exp.y - orbit_exp.v_y*orbit_exp.x)
	
end

# ╔═╡ be68781b-d2f9-4f76-8120-aa9f99c4a0f5
x_cen_exp = [orbit_exp.x orbit_exp.y orbit_exp.z]'

# ╔═╡ 17a0c0f3-37dc-4301-ac0c-3eaebdbdecdd
v_cen_exp = [orbit_exp.v_x orbit_exp.v_y orbit_exp.v_z]'

# ╔═╡ 199a0b4c-b828-4070-8628-50c27d8664ee
x_cen = out.x_cen

# ╔═╡ 4c8b7c67-c4c3-4616-82b7-546df559599d
v_cen = out.v_cen

# ╔═╡ 11c39da0-cf98-4f3f-a08b-ec9075135adb
times = out.times

# ╔═╡ acf9bb4c-c03f-455f-bd12-65766f55bfaa
md"""
# Observational reference
"""

# ╔═╡ f7659682-7f06-4618-aa70-4f2c7d3ee838
obs_props_file = modeldir * "simulation/orbit.toml"

# ╔═╡ eee089cf-fa34-454d-afc6-b195cbb80ed9
r_h = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxyname/observed_properties.toml"))["rh"] / 60

# ╔═╡ 5906da72-48ce-44f6-b4d7-2738c2df64b5
icrs = LilGuys.coord_from_file(obs_props_file)

# ╔═╡ 1cc54a8f-44df-4f35-aa01-3244bdeae52e
obs_props = TOML.parsefile(obs_props_file)

# ╔═╡ 040e03a3-ce59-4ee2-bcdb-1b78168c1794
icrs_u = LilGuys.ICRS(
	ra = icrs.ra ± r_h,
	dec = icrs.dec ± r_h,
	pmra = icrs.pmra ± obs_props["pmra_err"],
	pmdec = icrs.pmdec ± obs_props["pmdec_err"],
	distance = icrs.distance ± obs_props["distance_err"],
	radial_velocity = icrs.radial_velocity ± obs_props["radial_velocity_err"],
)

# ╔═╡ 4fba414e-ea58-46b5-8e88-22282690df13
gc = LilGuys.transform(Galactocentric, icrs_u)

# ╔═╡ 1258d03a-d4c4-4047-97cc-0fe718f2107b
gc2 = LilGuys.transform(Galactocentric, icrs)

# ╔═╡ 9344b6b4-0b06-4638-bfe0-1c5b3d024d3c
md"""
# New initial conditions
"""

# ╔═╡ 9a4c1b07-8b98-4d71-938e-736e0dcafc9b
orbit_in = calc_orbit(gc2, pot, time=10/T2GYR)

# ╔═╡ ec74d9e2-b2ae-447d-b9c6-5d55dd464d28
pos_i = [orbit_exp.x[1], orbit_exp.y[1], orbit_exp.z[1]]

# ╔═╡ fcb9b299-d660-48ba-afec-93c52a40d930
vel_i = [orbit_exp.v_x[1], orbit_exp.v_y[1], orbit_exp.v_z[1]]

# ╔═╡ 9c9a1dee-881d-4e7d-b9b3-321afde406da
L_in = pos_i × vel_i

# ╔═╡ 7bcdef79-7391-4f35-902b-62d2dda1d858
Φ_in(x::AbstractVector{<:Real}) = pyconvert(Float64, pot.potential(x))

# ╔═╡ 4f3330e4-455d-43ee-8cd1-244e3288393e
Φ_in(x::AbstractMatrix{<:Real}) = py2vec(pot.potential(np.array(x')))

# ╔═╡ 19803b8e-df6a-41cd-a2dd-dc4f04f95a46
x_cen[:, 1], v_cen[:, 1]

# ╔═╡ 9d49a44d-7147-4d3b-b848-aedfe1216219
pos_i, vel_i

# ╔═╡ e34091d8-c9f8-40f0-abb5-1d6b12b002cd
md"""
## Angular momentum method (old)
"""

# ╔═╡ 50b9fe74-de86-495c-8de7-f1e7e9b1bde5
E_in = Φ_in(pos_i) + 1/2 * calc_r(vel_i)^2

# ╔═╡ 083eb8df-ed6f-4f01-a5f9-5029c6fd8e82
function get_peri_apo(Φ, pos_i, vel_i, time=10/T2GYR)
	gc = LilGuys.Galactocentric(pos_i, vel_i*V2KMS)
	orbit = calc_orbit(gc, Φ, time=time)
	return extrema(calc_r(orbit.position))
end

# ╔═╡ ac9008eb-02b4-40b2-83fa-486fb5f3676e
function calc_peri_apo_simple(Φ, L, E, x0=[1,0,0])
	apo = LilGuys.find_zero(x -> -E + Φ_in(x0 * x) + calc_r(L)^2 / (2*x^2), 200)
	peri = LilGuys.find_zero(x -> -E + Φ_in(x0 * x) + calc_r(L)^2 / (2*x^2), 2)

	return peri, apo
end

# ╔═╡ c3080823-05de-4dd4-8a37-8117e97e4561
"""
	shift_initial_conditions()

Older method to shift initial conditions based on angular momentum. 
Abandoned since shifting positions and velocities to agree with shifts in phase
is much more nontrivial than with action angles.
"""
function shift_initial_conditions(Φ, pos_i, vel_i; dE=0, dlogr=0, dLx=0, dLy=0, dLz=0)
	r_i = LilGuys.calc_r(pos_i)
	r_hat = pos_i / r_i
	
	L_i = LilGuys.calc_L_spec(pos_i, vel_i)

	E_i = Φ(pos_i) + 1/2 * calc_r(vel_i)^2

	r_apo_i = calc_peri_apo_simple(Φ, L_i, E_i, r_hat)[2]

	E_new = E_i + dE
	L_new = [L_i[1] + dLx, L_i[2] + dLy, L_i[3] + dLz]
	r_apo_new = calc_peri_apo_simple(Φ, L_new, E_new, r_hat)[2]

	r_new = r_i * r_apo_new / r_apo_i * 10^dlogr
	pos_new = r_hat * r_new
	
	v_new = sqrt(2E_new - 2*Φ(pos_new))		

	vel_new_sign = sign(pos_i ⋅ vel_i)
	vel_new = v_new * (L_new × r_hat / r_new/v_new + vel_new_sign*r_hat * sqrt(1-calc_r(L_new/r_new/v_new)^2))

	println(L_new)
	println(pos_new × vel_new)
	println(E_new, "\t", Φ(pos_new) + 1/2*calc_r(vel_new)^2)
	println("r", r_new, "\t", calc_r(pos_new), "\t", r_i)
	println("apo new\t", r_apo_new, "\t", r_apo_i)

	return pos_new, vel_new
end

# ╔═╡ 2bb30ad3-2ddc-47bb-8356-8e7d647d5623
calc_peri_apo_simple(pot, L_in, E_in, pos_i / calc_r(pos_i))

# ╔═╡ 0b36f6e2-5e0e-43b1-aa0b-8de785748576
get_peri_apo(pot, pos_i, vel_i)

# ╔═╡ ecd27c0a-47a2-4fc3-8ade-7499fc36a496
shift_initial_conditions(Φ_in, pos_i, vel_i)

# ╔═╡ 5bfc3f07-853e-45ea-984a-07a76c372a18
md"""
# Actions, Orbits, and validation
"""

# ╔═╡ 95149181-b6c1-429e-be8e-f46d1a2616ef
function get_actions(pot, pos, vel)
	af = agama.ActionFinder(pot)

	act_py, ang_py, freq_py = af(np.array(vcat(pos, vel)'), angles=true, frequencies=true)

	return py2mat(act_py), py2mat(ang_py)
end

# ╔═╡ 3ab820ad-1e95-4533-aa70-ad56b11b549a
am = agama.ActionMapper(pot)

# ╔═╡ 3b247e3e-6aab-4c3d-9904-94123ccb00b7
function from_actions(pot, act, ang)
	am = agama.ActionMapper(pot)

	xv_py = am(np.array(vcat(act, ang)'))

	xv = py2mat(xv_py)
	return xv[1:3, :], xv[4:6, :]
end

# ╔═╡ 6b2fa1a5-6031-4c37-92ac-6a6088570dba
function shift_actions(Φ, pos_i, vel_i; dJ, dθ)
	act, ang = get_actions(Φ, pos_i, vel_i)
	act .+= dJ
	ang .+= dθ
	

	return from_actions(Φ, act, ang)
end

# ╔═╡ b3fd0308-2056-4d16-8ff2-bb4d17aaa131
shift_actions(pot, pos_i, vel_i, dJ=zeros(3), dθ=zeros(3))

# ╔═╡ 8303b321-1b95-4c7c-8c4f-d4c6f4e1b344
act_nbody, ang_nbody = get_actions(pot, x_cen, v_cen)

# ╔═╡ e2831eac-7c7a-4972-98f9-c6c500e91a41
x_obs = Measurements.value.(LilGuys.position_of(gc))

# ╔═╡ 7ebc25cb-c10a-4b38-b170-760aeeff38bf
v_obs = Measurements.value.(LilGuys.velocity_of(gc)) / V2KMS

# ╔═╡ 275c9856-a706-4926-b76f-f26b8b51203f
act_obs, ang_obs = get_actions(pot, x_obs, v_obs)

# ╔═╡ 97fb1541-59b6-4a42-adf3-f3084478ea04
act_exp, ang_exp = get_actions(pot, x_cen_exp, v_cen_exp)

# ╔═╡ 471fae8a-4998-4662-87e9-d4277f5f8604
L_exp = LilGuys.calc_L_spec(x_cen_exp, v_cen_exp)

# ╔═╡ 74a73f14-426b-45b6-afdd-55bc4a919b34
L_nbody = LilGuys.calc_L_spec(x_cen, v_cen)

# ╔═╡ b8f2475a-3588-4611-99ac-de156f25b853
let
	global dJ_suggested

	dJ_suggested = zeros(3)
	fig = Figure()

	for i in 1:3
		coord = ["r", "z", "ϕ"][i]
		
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"J_%$coord"
		)

		lines!(out.times .- out.times[idx_f], act_nbody[i, :])
		lines!(orbit_exp.t .- orbit_exp.t[end], act_exp[i, :])
		#lines!(orbit.time .- orbit.time[end], act_new[i, :])
		hlines!(act_obs[i], color="black")

		didx = 5
		acts = act_nbody[i, idx_f-didx:idx_f]
		act = LilGuys.mean(acts)
		
		lines!(out.times[[idx_f-didx, idx_f]] .- out.times[idx_f], [act, act])
		
		dJ_suggested[i] = act_obs[i] - act
		println("J", coord, "\t", dJ_suggested[i], "\t ± ", LilGuys.std(acts))

		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

	end

	fig
end


# ╔═╡ 2ceeecba-d375-4de9-86c4-a7278c3d016f
Number

# ╔═╡ 571ae542-406c-4e62-9b8f-af34e777748f
dθ_suggested = [ang_obs[i] .- ang_nbody[i, idx_f] for i in 1:3]

# ╔═╡ ae8f7bbd-46ee-4fcb-900a-557e5bb11088
@bind change_in_act_angles confirm(notebook_inputs(;
	Jr = NumberField(-1:0.01:1, default=dJ_suggested[1]),
	Jz = NumberField(-1:0.01:1, default=dJ_suggested[2]),
	Jϕ = NumberField(-1:0.01:1, default=dJ_suggested[3]),
	Θr = NumberField(-π:0.01:π, default=dθ_suggested[1]),
	Θz = NumberField(-π:0.01:π, default=dθ_suggested[2]),
	Θϕ = NumberField(-π:0.01:π, default=dθ_suggested[3]),
))

# ╔═╡ 9da81ee9-5638-48dc-949e-603874b3f627
dJ = [change_in_act_angles[:Jr], change_in_act_angles[:Jz], change_in_act_angles[:Jϕ]]

# ╔═╡ 0fef057f-2f68-447d-9926-8704536cc7e5
dθ = [change_in_act_angles[:Θr], change_in_act_angles[:Θz], change_in_act_angles[:Θϕ]]

# ╔═╡ 9facd961-6c18-4b4a-9e84-9a4db31600c1
begin 
	pos_new, vel_new = shift_actions(pot, pos_i, vel_i, dJ=dJ, dθ=dθ)
	#shift_initial_conditions(Φ_in, pos_i, vel_i, dE=0.01, dLz=-1.0, dLy=+0.0)
	pos_new = vec(pos_new)
	vel_new = vec(vel_new)
	pos_new, vel_new
end

# ╔═╡ 80556944-3d5f-4e4b-b3e3-68b092d006cd
gc_new = LilGuys.Galactocentric(pos_new, vel_new*V2KMS)

# ╔═╡ d5fda85b-45b8-442b-b604-4d45236301d4
orbit = calc_orbit(gc_new, pot, time=10 / T2GYR)

# ╔═╡ bb0195a0-d06d-4f72-8646-31ef848e987a
x_new = orbit.position

# ╔═╡ a419c394-d597-4f6f-812f-cd0a0bc7ea0e
LilGuys.plot_xyz(x_cen_exp, x_new, labels=["old point orbit", "new point orbit"])

# ╔═╡ af07bdd5-4baa-481d-a973-a4c119af6c5f
v_new = orbit.velocity

# ╔═╡ 3fa86c92-e6bd-41b6-bd64-cfd33f74b229
act_new, ang_new = get_actions(pot, x_new, v_new)

# ╔═╡ f8e73044-8d07-4b09-81ee-54f7f7278496
L_new = LilGuys.calc_L_spec(x_new, v_new)

# ╔═╡ 9ea9076a-2d51-41ca-ac10-c71a7e379c79
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = L"r",
	)
	lines!(out.times, calc_r(x_cen), label = "nbody")

	lines!(orbit_exp.t, calc_r(x_cen_exp), label="initial point")
	lines!(orbit.time, calc_r(orbit.position), label="next test")

	axislegend()
	fig
end

# ╔═╡ 7a18885a-364f-419b-85cf-d43519f4a37c
minimum(calc_r(orbit.position))

# ╔═╡ 1b369070-1c75-4879-976c-4b80ff66ac67
let
	fig = Figure()

	for i in 1:3
		coord = ["r", "z", "ϕ"][i]
		
		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"J_%$coord"
		)

		lines!(out.times .- out.times[idx_f], act_nbody[i, :])
		lines!(orbit_exp.t .- orbit_exp.t[end], act_exp[i, :])
		lines!(orbit.time .- orbit.time[end], act_new[i, :])
		hlines!(act_obs[i], color="black")

		
		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

	end

	fig
end


# ╔═╡ b39cbb71-b98f-4a6b-ba01-55d5b2bb2190
get_actions(pot, pos_new, vel_new)

# ╔═╡ 706a0753-018a-48f7-8c77-23af747141fd
let
	fig = Figure()
	local ax

	for i in 1:3
		coord = ["r", "z", "ϕ"][i]

		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"\theta_%$coord",
		)

		idx =idx_f-5:idx_f+0
		lines!(out.times[idx] .- out.times[idx_f], ang_nbody[i, idx])
		hlines!(ang_obs[i], color="black", label="observed")

		lines!(orbit_exp.t .- orbit_exp.t[end], ang_exp[i, :], label="initial point")
		lines!(orbit.time .- orbit.time[end], ang_new[i, :], label="new point")

		xlims!(-50, 0)
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
let
	fig = Figure()

	local ax
	for i in 1:3
		coord = ["r", "z", "ϕ"][i]

		ax = Axis(fig[i, 1],
			xlabel = "time",
			ylabel = L"\theta_%$coord",
		)
		
		lines!(out.times , ang_nbody[i, :], label="nbody")
		lines!(orbit_exp.t, ang_exp[i, :], label="initial point")
		lines!(orbit.time, ang_new[i, :], label="final point")

		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end
		println(ang_nbody[i, 1], "\t", ang_new[i, 1], "\t", dθ[i])

	end
	Legend(fig[2, 2], ax)

	fig
end


# ╔═╡ 947aba67-2d02-40e7-9e9c-16a39a0617be
act_nbody

# ╔═╡ 4473d1b0-fea7-4d87-8837-1966124f1cf0
E_exp = Φ_in(x_cen_exp) .+ 1/2 * calc_r(v_cen_exp) .^2

# ╔═╡ e2f1e41d-4996-4f81-8dda-35cf6ba0fbba
E_new = Φ_in(x_new) .+ 1/2 * calc_r(v_new) .^2

# ╔═╡ c2b025c6-0bd9-49be-9848-3f79e31e6140
E_nbody = Φ_in(x_cen) .+ 1/2 * calc_r(v_cen) .^2

# ╔═╡ b7049ed6-183a-4920-8d89-ac396979ae2b
E_obs = Φ_in(Measurements.value.(LilGuys.position_of(gc))) .+ 1/2 * calc_r(Measurements.value.(LilGuys.velocity_of(gc) / V2KMS)) .^ 2

# ╔═╡ da3f27fc-0f2c-4e03-bc36-ccbb5eee1ea2
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "time",
		ylabel = "energy"
	)
	
	lines!(out.times, E_nbody)
	lines!(orbit_exp.t, E_exp)
	lines!(orbit.time, E_new)

	hlines!(E_obs, color="black")
	fig
end


# ╔═╡ 963bd378-2603-45b1-8623-979aad5c2538
df_new = DataFrame(
	t = orbit.time,
	x = x_new[1, :],
	y = x_new[2, :],
	z = x_new[3, :],
	v_x = v_new[1, :],
	v_y = v_new[2, :],
	v_z = v_new[3, :],
)

# ╔═╡ 0aea9fe2-65ef-42f3-bf58-fa3b81dbba8f
@bind write_ic PlutoUI.Button()

# ╔═╡ 4fda883d-7173-49c0-a56e-6572ebc80ba4
let write_ic
	@info "writing  to $modeldir/next_orbit..."
	
	CSV.write("$modeldir/next_orbit.csv", df_new)
	open("$modeldir/next_orbit_shifts.toml", "w") do f
		d = OrderedDict(
			"lastmodel" => orbitname,
			"dJ" => dJ,
			"dθ" => dθ,
			"Ji" => act_new[:, 1],
			"Θi" => ang_new[:, 1],
			"position_i" => orbit.position[:, 1],
			"velocity_i" => orbit.velocity[:, 1],
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

# ╔═╡ 6946c00b-2f1d-4d97-8434-ca93cc6cf1c9
obs_props

# ╔═╡ 96ef21f4-1915-46ed-9c82-94b3bba456ac
L_obs = LilGuys.calc_L_spec(LilGuys.position_of(gc), LilGuys.velocity_of(gc) / V2KMS)

# ╔═╡ 18b29878-b8e3-44e5-8843-ccde9c6bb273
function plot_meas!(meas)
	hlines!(Measurements.value(meas), color=:black, linestyle=:dot)

	y = Measurements.value(meas)
	ye = Measurements.uncertainty(meas)

	hspan!(y-ye, y+ye, color=(:black, 0.2))

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
		lines!(orbit_exp.t, L_exp[i, :])
		lines!(orbit.time, L_new[i, :])
		plot_meas!(L_obs[i])

		if i < 3
			hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
		end

	end

	fig
end


# ╔═╡ Cell order:
# ╟─9bfa1465-2ca9-495c-bab6-5110072f02ed
# ╟─0c21b10f-ad86-4894-8959-721742d2a2c1
# ╠═5006f134-648a-4318-9dc3-5de383ac4d0e
# ╠═e9c9d1d9-42e6-4732-b425-477ec0bd7589
# ╠═b18d1734-e40e-11ef-1007-314371eb1a54
# ╠═34e54152-cebb-4364-b1b7-37e0dd4f10c8
# ╠═b35b362e-f275-4d70-a703-f398a08b7a47
# ╠═a47e8c62-2b01-4be4-b91a-20b08f592160
# ╠═f0f7b2db-b926-46d5-a9bc-baa9cff0b6cf
# ╠═4bc528f9-4ac0-407e-88a6-5e606117a3d3
# ╠═0b6b20b0-bba8-4ed3-b8a5-d8361bc3674f
# ╠═481edcca-51a8-4f11-88ec-84939b223bdf
# ╠═598be49e-82ef-4820-8996-ea2b753f4275
# ╠═fba1564a-c12d-4657-86fe-d4d771be1037
# ╠═0b6b2020-282b-4384-ac58-8665364cf86e
# ╠═06fbb3db-26ae-4380-935c-883a02c6e119
# ╠═c8e2c46d-917c-4fa8-a841-ea26430d17d3
# ╟─b5302ffe-76e1-424d-9e08-1656981b068e
# ╠═9cbf9f1e-1a61-48ae-9fd8-d4253a98f98b
# ╠═de8bdfef-efca-4d2a-a9cd-a79ccc931b3f
# ╠═e279753d-6be6-4e47-a72a-dc0983acef5e
# ╠═fd1adb47-9a55-49c9-b1ca-5009d91b0408
# ╠═be68781b-d2f9-4f76-8120-aa9f99c4a0f5
# ╠═17a0c0f3-37dc-4301-ac0c-3eaebdbdecdd
# ╠═199a0b4c-b828-4070-8628-50c27d8664ee
# ╠═4c8b7c67-c4c3-4616-82b7-546df559599d
# ╠═11c39da0-cf98-4f3f-a08b-ec9075135adb
# ╟─acf9bb4c-c03f-455f-bd12-65766f55bfaa
# ╠═f7659682-7f06-4618-aa70-4f2c7d3ee838
# ╠═eee089cf-fa34-454d-afc6-b195cbb80ed9
# ╠═5906da72-48ce-44f6-b4d7-2738c2df64b5
# ╠═1cc54a8f-44df-4f35-aa01-3244bdeae52e
# ╠═040e03a3-ce59-4ee2-bcdb-1b78168c1794
# ╠═4fba414e-ea58-46b5-8e88-22282690df13
# ╠═1258d03a-d4c4-4047-97cc-0fe718f2107b
# ╟─9344b6b4-0b06-4638-bfe0-1c5b3d024d3c
# ╠═9a4c1b07-8b98-4d71-938e-736e0dcafc9b
# ╠═9c9a1dee-881d-4e7d-b9b3-321afde406da
# ╠═6b2fa1a5-6031-4c37-92ac-6a6088570dba
# ╠═ec74d9e2-b2ae-447d-b9c6-5d55dd464d28
# ╠═fcb9b299-d660-48ba-afec-93c52a40d930
# ╠═7bcdef79-7391-4f35-902b-62d2dda1d858
# ╠═4f3330e4-455d-43ee-8cd1-244e3288393e
# ╠═9facd961-6c18-4b4a-9e84-9a4db31600c1
# ╠═19803b8e-df6a-41cd-a2dd-dc4f04f95a46
# ╠═9d49a44d-7147-4d3b-b848-aedfe1216219
# ╠═b3fd0308-2056-4d16-8ff2-bb4d17aaa131
# ╟─e34091d8-c9f8-40f0-abb5-1d6b12b002cd
# ╟─c3080823-05de-4dd4-8a37-8117e97e4561
# ╠═50b9fe74-de86-495c-8de7-f1e7e9b1bde5
# ╠═083eb8df-ed6f-4f01-a5f9-5029c6fd8e82
# ╠═ac9008eb-02b4-40b2-83fa-486fb5f3676e
# ╠═2bb30ad3-2ddc-47bb-8356-8e7d647d5623
# ╠═0b36f6e2-5e0e-43b1-aa0b-8de785748576
# ╠═ecd27c0a-47a2-4fc3-8ade-7499fc36a496
# ╟─5bfc3f07-853e-45ea-984a-07a76c372a18
# ╠═95149181-b6c1-429e-be8e-f46d1a2616ef
# ╠═3ab820ad-1e95-4533-aa70-ad56b11b549a
# ╠═3b247e3e-6aab-4c3d-9904-94123ccb00b7
# ╠═bb0195a0-d06d-4f72-8646-31ef848e987a
# ╠═af07bdd5-4baa-481d-a973-a4c119af6c5f
# ╠═8303b321-1b95-4c7c-8c4f-d4c6f4e1b344
# ╠═e2831eac-7c7a-4972-98f9-c6c500e91a41
# ╠═7ebc25cb-c10a-4b38-b170-760aeeff38bf
# ╠═275c9856-a706-4926-b76f-f26b8b51203f
# ╠═97fb1541-59b6-4a42-adf3-f3084478ea04
# ╠═3fa86c92-e6bd-41b6-bd64-cfd33f74b229
# ╠═80556944-3d5f-4e4b-b3e3-68b092d006cd
# ╠═d5fda85b-45b8-442b-b604-4d45236301d4
# ╠═b39cbb71-b98f-4a6b-ba01-55d5b2bb2190
# ╠═a419c394-d597-4f6f-812f-cd0a0bc7ea0e
# ╠═9ea9076a-2d51-41ca-ac10-c71a7e379c79
# ╠═7a18885a-364f-419b-85cf-d43519f4a37c
# ╠═f8e73044-8d07-4b09-81ee-54f7f7278496
# ╠═471fae8a-4998-4662-87e9-d4277f5f8604
# ╠═74a73f14-426b-45b6-afdd-55bc4a919b34
# ╠═0f936517-1924-44db-ba03-2c7fd589afa3
# ╠═b8f2475a-3588-4611-99ac-de156f25b853
# ╟─1b369070-1c75-4879-976c-4b80ff66ac67
# ╠═9da81ee9-5638-48dc-949e-603874b3f627
# ╠═0fef057f-2f68-447d-9926-8704536cc7e5
# ╠═2ceeecba-d375-4de9-86c4-a7278c3d016f
# ╠═ae8f7bbd-46ee-4fcb-900a-557e5bb11088
# ╠═571ae542-406c-4e62-9b8f-af34e777748f
# ╠═706a0753-018a-48f7-8c77-23af747141fd
# ╠═e7aba352-84a1-417c-b25b-7d0503030e6a
# ╠═2af93853-0df1-45a5-ace6-75ad532e49dd
# ╠═da3f27fc-0f2c-4e03-bc36-ccbb5eee1ea2
# ╠═947aba67-2d02-40e7-9e9c-16a39a0617be
# ╠═4473d1b0-fea7-4d87-8837-1966124f1cf0
# ╠═e2f1e41d-4996-4f81-8dda-35cf6ba0fbba
# ╠═c2b025c6-0bd9-49be-9848-3f79e31e6140
# ╠═b7049ed6-183a-4920-8d89-ac396979ae2b
# ╠═963bd378-2603-45b1-8623-979aad5c2538
# ╠═e13b1c88-75c4-498a-87c5-6d11c88fa569
# ╠═0aea9fe2-65ef-42f3-bf58-fa3b81dbba8f
# ╠═4fda883d-7173-49c0-a56e-6572ebc80ba4
# ╠═6946c00b-2f1d-4d97-8434-ca93cc6cf1c9
# ╠═96ef21f4-1915-46ed-9c82-94b3bba456ac
# ╠═18b29878-b8e3-44e5-8843-ccde9c6bb273
