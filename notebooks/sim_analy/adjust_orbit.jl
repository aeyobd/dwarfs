### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ b18d1734-e40e-11ef-1007-314371eb1a54
begin
	using Pkg;Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ 4bc528f9-4ac0-407e-88a6-5e606117a3d3
using CSV, DataFrames

# ╔═╡ c9786ddd-b2b0-401e-ab97-585fcbea31cb
begin
	using PythonCall
	np = pyimport("numpy")
end

# ╔═╡ 481edcca-51a8-4f11-88ec-84939b223bdf
include(ENV["DWARFS_ROOT"] * "/utils/agama_utils.jl")

# ╔═╡ b4c4733c-3c62-4480-9fc7-0511a20a0357
py2mat(x) = pyconvert(Matrix{Float64}, x)'

# ╔═╡ 60a4b629-ad00-436e-a9d4-a29201aa3e91
py2vec(x) = pyconvert(Vector{Float64}, x)

# ╔═╡ 8473ee55-78a6-45bc-9db5-ce4d9f636144
py2f(x) = pyconvert(Float64, x)

# ╔═╡ 1959e6a0-9ad0-4bb3-955f-662ae9024eb7
⊕(x::Real, y::Real) = sqrt(x^2 + y^2)

# ╔═╡ 121d4497-ee8c-40f0-a7a4-e04c6650cb25
galaxyname = "ursa_minor"

# ╔═╡ f5ec8d3e-d34d-481f-aec8-0113d7d9402a
modelname = "1e6_v37_r5.0"

# ╔═╡ a80b0905-7638-495b-aa97-a83478005492
orbitname = "orbit_mean"

# ╔═╡ c8e2c46d-917c-4fa8-a841-ea26430d17d3
modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, orbitname) * "/"

# ╔═╡ de8bdfef-efca-4d2a-a9cd-a79ccc931b3f
pot = agama.Potential(modeldir * "simulation/agama_potential.ini")

# ╔═╡ 9cbf9f1e-1a61-48ae-9fd8-d4253a98f98b
out = Output(modeldir)

# ╔═╡ 992b7c5d-02c0-4978-a063-26dd80256ad9


# ╔═╡ fd1adb47-9a55-49c9-b1ca-5009d91b0408
begin 
	orbit_exp = CSV.read(modeldir * "simulation/orbit.csv", DataFrame)
	orbit_exp.t .-= orbit_exp.t[1]
	orbit_exp[!, :R] = orbit_exp.x .⊕ orbit_exp.y
	orbit_exp[!, :Vr] = @. 1/orbit_exp.R *( orbit_exp.v_x * orbit_exp.x + orbit_exp.v_y*orbit_exp.y)
	orbit_exp[!, :Vphi] = @. 1/orbit_exp.R *( orbit_exp.v_x * orbit_exp.y - orbit_exp.v_y*orbit_exp.x)
	
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

# ╔═╡ 9344b6b4-0b06-4638-bfe0-1c5b3d024d3c
md"""
# New initial conditions
"""

# ╔═╡ 083eb8df-ed6f-4f01-a5f9-5029c6fd8e82
function get_peri_apo(Φ, pos_i, vel_i, time=10/T2GYR)
	gc = LilGuys.Galactocentric(pos_i, vel_i*V2KMS)
	orbit = calc_orbit(gc, Φ, time=time)
	return extrema(calc_r(orbit.position))
end

# ╔═╡ 34e54152-cebb-4364-b1b7-37e0dd4f10c8
import LinearAlgebra: ⋅, ×

# ╔═╡ ec74d9e2-b2ae-447d-b9c6-5d55dd464d28
pos_i = [orbit_exp.x[1], orbit_exp.y[1], orbit_exp.z[1]]

# ╔═╡ fcb9b299-d660-48ba-afec-93c52a40d930
vel_i = [orbit_exp.v_x[1], orbit_exp.v_y[1], orbit_exp.v_z[1]]

# ╔═╡ 0b36f6e2-5e0e-43b1-aa0b-8de785748576
get_peri_apo(pot, pos_i, vel_i)

# ╔═╡ 9c9a1dee-881d-4e7d-b9b3-321afde406da
L_in = pos_i × vel_i

# ╔═╡ 7bcdef79-7391-4f35-902b-62d2dda1d858
Φ_in(x) = pyconvert(Float64, pot.potential(x))

# ╔═╡ 50b9fe74-de86-495c-8de7-f1e7e9b1bde5
E_in = Φ_in(pos_i) + 1/2 * calc_r(vel_i)^2

# ╔═╡ ac9008eb-02b4-40b2-83fa-486fb5f3676e
function calc_peri_apo_simple(Φ, L, E, x0=[1,0,0])
	apo = LilGuys.find_zero(x -> -E + Φ_in(x0 * x) + calc_r(L)^2 / (2*x^2), 200)
	peri = LilGuys.find_zero(x -> -E + Φ_in(x0 * x) + calc_r(L)^2 / (2*x^2), 2)

	return peri, apo
end

# ╔═╡ 2bb30ad3-2ddc-47bb-8356-8e7d647d5623
calc_peri_apo_simple(pot, L_in, E_in, pos_i / calc_r(pos_i))

# ╔═╡ c3080823-05de-4dd4-8a37-8117e97e4561
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

# ╔═╡ 62f1bb5e-b0e9-4d9c-af99-be0df8f7dca1
LilGuys.find_zero(x -> -E_in + Φ_in([0,0,1] * x) + calc_r(L_in)^2 / (2*x^2), 80)

# ╔═╡ d37c6844-cc84-489f-8041-c1e4c093b4b6
LilGuys.find_zero(x -> -E_in + Φ_in([1,0,0] * x) + calc_r(L_in)^2 / (2*x^2), 2)

# ╔═╡ 9facd961-6c18-4b4a-9e84-9a4db31600c1
pos_new, vel_new = shift_initial_conditions(Φ_in, pos_i, vel_i, dE=0.05, dLx=-3)

# ╔═╡ 19803b8e-df6a-41cd-a2dd-dc4f04f95a46
x_cen[:, 1], v_cen[:, 1]

# ╔═╡ 9d49a44d-7147-4d3b-b848-aedfe1216219
pos_i, vel_i

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

# ╔═╡ 3b247e3e-6aab-4c3d-9904-94123ccb00b7
function from_actions(pot, act, ang)
	am = agama.ActionMapper(pot)

	xv_py, Ohm = am(np.array(vcat(act, ang)'))

	xv = py2mat(xv_py)
	return xv[1:3, :], xv[4:6, :]
end

# ╔═╡ 8303b321-1b95-4c7c-8c4f-d4c6f4e1b344
act_nbody, _ = get_actions(pot, x_cen, v_cen)

# ╔═╡ 97fb1541-59b6-4a42-adf3-f3084478ea04
act_exp, _ = get_actions(pot, x_cen_exp, v_cen_exp)

# ╔═╡ 3f958975-7abf-4564-9bac-432a87e66bcc
gc = LilGuys.Galactocentric(pos_i, vel_i*V2KMS)

# ╔═╡ 9a4c1b07-8b98-4d71-938e-736e0dcafc9b
orbit_in = calc_orbit(gc, pot, time=10/T2GYR)

# ╔═╡ 53ce58a7-8ddd-4c3e-96b8-180c06623ce0
calc_r(orbit_in.position)

# ╔═╡ 80556944-3d5f-4e4b-b3e3-68b092d006cd
gc_new = LilGuys.Galactocentric(pos_new, vel_new*V2KMS)

# ╔═╡ d5fda85b-45b8-442b-b604-4d45236301d4
orbit = calc_orbit(gc_new, pot, time=10 / T2GYR)

# ╔═╡ bb0195a0-d06d-4f72-8646-31ef848e987a
x_new = orbit.position

# ╔═╡ af07bdd5-4baa-481d-a973-a4c119af6c5f
v_new = orbit.velocity

# ╔═╡ 3fa86c92-e6bd-41b6-bd64-cfd33f74b229
act_new, _ = get_actions(pot, x_new, v_new)

# ╔═╡ b39cbb71-b98f-4a6b-ba01-55d5b2bb2190
get_actions(pot, pos_new, vel_new)

# ╔═╡ a419c394-d597-4f6f-812f-cd0a0bc7ea0e
LilGuys.plot_xyz(orbit_in.position, x_cen_exp, x_new)

# ╔═╡ 9ea9076a-2d51-41ca-ac10-c71a7e379c79
let
	fig = Figure()
	ax = Axis(fig[1,1])
	lines!(out.times, calc_r(x_cen))

	lines!(orbit_exp.t, calc_r(x_cen_exp))
	lines!(orbit.time, calc_r(orbit.position))

	fig
end

# ╔═╡ b8f2475a-3588-4611-99ac-de156f25b853
let
	fig = Figure()

	for i in 1:3
		ax = Axis(fig[i, 1])

		lines!(out.times, act_nbody[i, :])
		lines!(orbit_exp.t, act_exp[i, :])
		lines!(orbit.time, act_new[i, :])

	end

	fig
end


# ╔═╡ f8e73044-8d07-4b09-81ee-54f7f7278496
L_new = LilGuys.calc_L_spec(x_new, v_new)

# ╔═╡ 471fae8a-4998-4662-87e9-d4277f5f8604
L_exp = LilGuys.calc_L_spec(x_cen_exp, v_cen_exp)

# ╔═╡ 74a73f14-426b-45b6-afdd-55bc4a919b34
L_nbody = LilGuys.calc_L_spec(x_cen, v_cen)

# ╔═╡ 0f936517-1924-44db-ba03-2c7fd589afa3
let
	fig = Figure()

	for i in 1:3
		ax = Axis(fig[i, 1])

		lines!(out.times, L_nbody[i, :])
		lines!(orbit_exp.t, L_exp[i, :])
		lines!(orbit.time, L_new[i, :])

	end

	fig
end


# ╔═╡ 0aea9fe2-65ef-42f3-bf58-fa3b81dbba8f
write_ic = true

# ╔═╡ 061f913b-b842-41a3-aa9f-a918c86dab0a
new_model_dir = joinpath(modeldir, "../$orbitname.1")

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

# ╔═╡ 934e3300-61d6-48f9-af57-791b96bfb839
0.037 * sqrt(10)

# ╔═╡ 4fda883d-7173-49c0-a56e-6572ebc80ba4
if write_ic
	if isdir(modeldir)
		mkpath(new_model_dir)
	end
	CSV.write("$new_model_dir/new_orbit.csv", df_new)
end

# ╔═╡ Cell order:
# ╠═b18d1734-e40e-11ef-1007-314371eb1a54
# ╠═4bc528f9-4ac0-407e-88a6-5e606117a3d3
# ╠═c9786ddd-b2b0-401e-ab97-585fcbea31cb
# ╠═481edcca-51a8-4f11-88ec-84939b223bdf
# ╠═b4c4733c-3c62-4480-9fc7-0511a20a0357
# ╠═60a4b629-ad00-436e-a9d4-a29201aa3e91
# ╠═8473ee55-78a6-45bc-9db5-ce4d9f636144
# ╠═1959e6a0-9ad0-4bb3-955f-662ae9024eb7
# ╠═121d4497-ee8c-40f0-a7a4-e04c6650cb25
# ╠═f5ec8d3e-d34d-481f-aec8-0113d7d9402a
# ╠═a80b0905-7638-495b-aa97-a83478005492
# ╠═c8e2c46d-917c-4fa8-a841-ea26430d17d3
# ╠═de8bdfef-efca-4d2a-a9cd-a79ccc931b3f
# ╠═9cbf9f1e-1a61-48ae-9fd8-d4253a98f98b
# ╠═992b7c5d-02c0-4978-a063-26dd80256ad9
# ╠═fd1adb47-9a55-49c9-b1ca-5009d91b0408
# ╠═be68781b-d2f9-4f76-8120-aa9f99c4a0f5
# ╠═17a0c0f3-37dc-4301-ac0c-3eaebdbdecdd
# ╠═199a0b4c-b828-4070-8628-50c27d8664ee
# ╠═4c8b7c67-c4c3-4616-82b7-546df559599d
# ╠═11c39da0-cf98-4f3f-a08b-ec9075135adb
# ╟─9344b6b4-0b06-4638-bfe0-1c5b3d024d3c
# ╠═083eb8df-ed6f-4f01-a5f9-5029c6fd8e82
# ╠═9a4c1b07-8b98-4d71-938e-736e0dcafc9b
# ╠═53ce58a7-8ddd-4c3e-96b8-180c06623ce0
# ╠═0b36f6e2-5e0e-43b1-aa0b-8de785748576
# ╠═2bb30ad3-2ddc-47bb-8356-8e7d647d5623
# ╠═50b9fe74-de86-495c-8de7-f1e7e9b1bde5
# ╠═9c9a1dee-881d-4e7d-b9b3-321afde406da
# ╠═ac9008eb-02b4-40b2-83fa-486fb5f3676e
# ╠═62f1bb5e-b0e9-4d9c-af99-be0df8f7dca1
# ╠═d37c6844-cc84-489f-8041-c1e4c093b4b6
# ╠═c3080823-05de-4dd4-8a37-8117e97e4561
# ╠═34e54152-cebb-4364-b1b7-37e0dd4f10c8
# ╠═ec74d9e2-b2ae-447d-b9c6-5d55dd464d28
# ╠═fcb9b299-d660-48ba-afec-93c52a40d930
# ╠═7bcdef79-7391-4f35-902b-62d2dda1d858
# ╠═9facd961-6c18-4b4a-9e84-9a4db31600c1
# ╠═19803b8e-df6a-41cd-a2dd-dc4f04f95a46
# ╠═9d49a44d-7147-4d3b-b848-aedfe1216219
# ╠═ecd27c0a-47a2-4fc3-8ade-7499fc36a496
# ╟─5bfc3f07-853e-45ea-984a-07a76c372a18
# ╠═95149181-b6c1-429e-be8e-f46d1a2616ef
# ╠═3b247e3e-6aab-4c3d-9904-94123ccb00b7
# ╠═bb0195a0-d06d-4f72-8646-31ef848e987a
# ╠═af07bdd5-4baa-481d-a973-a4c119af6c5f
# ╠═8303b321-1b95-4c7c-8c4f-d4c6f4e1b344
# ╠═97fb1541-59b6-4a42-adf3-f3084478ea04
# ╠═3fa86c92-e6bd-41b6-bd64-cfd33f74b229
# ╠═3f958975-7abf-4564-9bac-432a87e66bcc
# ╠═80556944-3d5f-4e4b-b3e3-68b092d006cd
# ╠═d5fda85b-45b8-442b-b604-4d45236301d4
# ╠═b39cbb71-b98f-4a6b-ba01-55d5b2bb2190
# ╠═a419c394-d597-4f6f-812f-cd0a0bc7ea0e
# ╠═9ea9076a-2d51-41ca-ac10-c71a7e379c79
# ╠═b8f2475a-3588-4611-99ac-de156f25b853
# ╠═f8e73044-8d07-4b09-81ee-54f7f7278496
# ╠═471fae8a-4998-4662-87e9-d4277f5f8604
# ╠═74a73f14-426b-45b6-afdd-55bc4a919b34
# ╠═0f936517-1924-44db-ba03-2c7fd589afa3
# ╠═0aea9fe2-65ef-42f3-bf58-fa3b81dbba8f
# ╠═061f913b-b842-41a3-aa9f-a918c86dab0a
# ╠═963bd378-2603-45b1-8623-979aad5c2538
# ╠═934e3300-61d6-48f9-af57-791b96bfb839
# ╠═4fda883d-7173-49c0-a56e-6572ebc80ba4
