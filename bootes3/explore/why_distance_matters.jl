### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 92946976-38f5-11f1-813f-833fa39bd393
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya, CairoMakie

end

# ╔═╡ e3fedce4-eb08-4242-beaa-55a6b714f3c0
using Distributions

# ╔═╡ 4bbe30af-1059-45aa-a6be-0000b0c48a04
using LinearAlgebra: normalize

# ╔═╡ 9240885e-6f30-4f4b-8904-fe43f9297ef3
md"""
# Introduction
The aim of this notebook is to understand why the distance uncertainty in the initial coordinates primarily affects the pericentre for Bootes 3
"""

# ╔═╡ 6eb088f3-66b0-4f0b-ac84-e6b88c8d8ed1
import TOML

# ╔═╡ eeb6217f-6530-40e3-a313-12d2f6c07b60
import Agama

# ╔═╡ ffc0957b-ffd9-4594-b9e9-85aeaf2e48ec
CairoMakie.activate!(type=:png)

# ╔═╡ 5e0b2289-e57f-48cf-a701-10d66dba1925
md"""
# data
"""

# ╔═╡ ec309258-c5de-40ad-abf2-01d6808163d2
pot = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))

# ╔═╡ cb89b58b-48d3-4269-8158-ecbc81464e9d
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/bootes3/observed_properties.toml"))

# ╔═╡ 6931ed06-99a2-4a3f-ac63-7dc7f4d67ef7
icrs0 = ICRS(obs_props)

# ╔═╡ c366985d-07a9-4ca3-85a8-ed36298cd501
md"""
# Samples
"""

# ╔═╡ 5ad2202b-ad37-4c1c-9f0a-0a83da6fdc67
N_samples = 1000

# ╔═╡ 2b90d979-58a9-4f1f-9813-1193ace020ad
dm = obs_props["distance_modulus"] 

# ╔═╡ 7d688254-2008-43f2-903c-f6d160e7aa0f
dm_err = obs_props["distance_modulus_em"]

# ╔═╡ 1e3b7236-7c35-497d-9157-e68d15898535
dist_dm = Normal(dm, dm_err) #Uniform(dm - 3*dm_err, dm + 3dm_err)

# ╔═╡ 2b8b0ee9-0111-42f4-90ba-44cd7bc433e3
dist_samples = LilGuys.dm2kpc.(rand(dist_dm, N_samples)) |> sort

# ╔═╡ 402e3945-d039-4fe7-a6ba-571e9be908c6
coord_from_dist(dist) = ICRS(merge(obs_props, Dict("distance" => dist)))

# ╔═╡ df363bdc-ac69-4725-bde6-81c5377c99ca
coord_from_dist(10)

# ╔═╡ 807628e0-83af-431f-8014-a34387b883d5
orbits = LilGuys.agama_orbit(pot, coord_from_dist.(dist_samples), timerange=(0, -10/T2GYR))

# ╔═╡ f5fe9855-0e8e-4fb7-9f46-0aefb97a1aa5
md"""
# Analysis/calculations
"""

# ╔═╡ 2f7d73bf-8510-4d32-bedc-7e9ff84171d8
L0 = hcat([LilGuys.angular_momenta(o)[:, 1] for o in orbits]...)

# ╔═╡ 1e8efe39-25d6-4b32-be58-1e6b6f2e31be
v0 = [LilGuys.speeds(o)[1] for o in orbits]

# ╔═╡ ba736ef1-2c1d-4091-af50-9f8ba45961b0
x0 = hcat([LilGuys.positions(o)[:, 1] for o in orbits]...)

# ╔═╡ e061cbeb-7a7d-4765-873a-3bda33931af7
vel0 = hcat([LilGuys.velocities(o)[:, 1] for o in orbits]...)

# ╔═╡ 64fe2dce-a020-4d7f-a546-19753136dae6
ϕ0 = Agama.potential(pot, x0)

# ╔═╡ 53465bc5-dca7-44a2-b137-49876fe054a7
ϵ0 = @. ϕ0 + 1/2*v0^2

# ╔═╡ d63fcaa5-4e5f-4a8b-b714-3d026f0dff30
LilGuys.velocity(LilGuys.transform(Galactocentric, icrs0))/V2KMS |> radii

# ╔═╡ d78c9fe3-4070-4cac-93f1-5a92cfacfbcb
peris = LilGuys.pericenter.(orbits)

# ╔═╡ 2fe531fd-9724-4e4e-997e-ba308d2f7005
md"""
# Plots
"""

# ╔═╡ 2615f2e5-f7b6-453d-8471-8c4a7fa3b566
let
	fig = Figure()
	ax = Axis(fig[1,1])

	for o in orbits
		lines!(o.positions[1, 1:10], o.positions[2, 1:10])
	end

	fig

end

# ╔═╡ 0fff09cc-a1fc-43e8-bc8b-1cb2e32d78fb
LilGuys.plot_xyz(LilGuys.positions.(orbits)[rand(eachindex(orbits), 100)]..., color=:black, alpha=0.1, linewidth=1)

# ╔═╡ 9d8add14-2706-4197-b99e-4060afaa00c7
let
	f = LilGuys.plot_xyz([LilGuys.positions(o)[:, 1:5] for o in orbits[1:300:end]]...,  color=:black, alpha=0.5, linewidth=0.5)

	LilGuys.plot_xyz!(f.content, [-8.122, 0, 0], plot=:scatter)

	f
end

# ╔═╡ 27215e5e-52f1-4f57-8f15-982887cfc71e
orbits[1].times[1:5] * T2GYR

# ╔═╡ ff51107b-579c-4a9a-9975-698fa1bd3bd7
gsr_mean = LilGuys.transform(GSR, icrs0, )

# ╔═╡ 2de6de24-f317-4faa-bf5c-2dc93736761e
icrs_pmra = GSR(ra=icrs0.ra, dec=icrs0.dec, distance=icrs0.distance, pmra=0.1, pmdec=0, radial_velocity=0)

# ╔═╡ 1606b369-d217-45bd-80b4-c0328518a0dd
icrs_pmdec = GSR(ra=icrs0.ra, dec=icrs0.dec, distance=icrs0.distance, pmra=0, pmdec=0.1, radial_velocity=0)

# ╔═╡ 06ac66d4-dcf8-4376-8b2e-23b70e5c89a0
icrs_radial_vel = GSR(ra=icrs0.ra, dec=icrs0.dec, distance=icrs0.distance, pmra=0, pmdec=0, radial_velocity=1)

# ╔═╡ 0ca9c4c6-7ee6-42f2-b9a6-de7d3149a0e9
v_pmra = LilGuys.velocity(LilGuys.transform(Galactocentric, icrs_pmra)) |> normalize

# ╔═╡ 784d2f59-adf4-4564-9f29-9e2f1727ffa1
v_pmdec = LilGuys.velocity(LilGuys.transform(Galactocentric, icrs_pmdec)) |> normalize

# ╔═╡ dc6b7ad1-7bcd-4462-8e03-e83f935a6daf


# ╔═╡ a74c7575-70e7-4301-b05f-3b5dee5e0e0c
v_rv = LilGuys.velocity(LilGuys.transform(Galactocentric, icrs_radial_vel)) |> normalize

# ╔═╡ 069efc5e-99ca-4cc1-9c94-75c5e67f80e3
mean(x0, dims=2) |> normalize

# ╔═╡ 4333f006-4f9f-4d9a-8eca-d6d7d7974cdd
LilGuys.pm2kms(gsr_mean.pmra, gsr_mean.distance)

# ╔═╡ c3b20d5e-5186-47e8-b32b-ad2dc49b6693
gsr_mean.pmra, gsr_mean.pmdec

# ╔═╡ b9423067-a382-48eb-afe3-3a617d8ac475
LilGuys.pm2kms(gsr_mean.pmdec, gsr_mean.distance)

# ╔═╡ 784b9c5e-7705-499d-a593-c47f48ca0899
v_sun = LilGuys.velocity(LilGuys.transform(Galactocentric, ICRS(ra=0., dec=0, distance=0, pmra=0, pmdec=0, radial_velocity=0)))

# ╔═╡ b744f1db-f69b-4b19-a14d-a76c210ff305
gsr_mean.pmra

# ╔═╡ ae4f0dfd-aa6d-4f93-a0a1-f0c03e4d341c
gsr_mean.pmdec

# ╔═╡ 20f32802-7793-4135-8411-7db1db8bc2b3
let
	fig = Figure(size=(5, 3) .* 72)

	ax = Axis(fig[1,1], xlabel="helcen distance / kpc", 
			 ylabel = "galcen / km/s")

	lines!(dist_samples, vel0[1, :] .* V2KMS, label="vx")
	vel_simple = @. (
		v_pmra* LilGuys.pm2kms.(icrs0.pmra, dist_samples')
			+ v_pmdec* LilGuys.pm2kms.(icrs0.pmdec, dist_samples')
				+ v_rv * icrs0.radial_velocity
		.+ v_sun
	)
	# lines!(dist_samples, vel_simple[1, :])


	lines!(dist_samples, vel0[2, :] * V2KMS, label="vy")
	# lines!(dist_samples, vel_simple[2, :])




	lines!(dist_samples, vel0[3, :] * V2KMS, label="vz")
	# lines!(dist_samples, vel_simple[3, :])
	vel_0 = @. (
		v_rv * icrs0.radial_velocity
		.+ v_sun
	)
	
	scatter!(5, v_sun[1], marker=:star5, label="sun")
	scatter!(5, v_sun[2], marker=:star5)
	scatter!(5, v_sun[3], marker=:star5)


	hlines!(0, color=:black)
	for i in 1:3
		scatter!(0, vel_0[i], marker=:circle, color=COLORS[i], strokewidth=0)
	end
	
	axislegend(position=:rt)
	fig
end

# ╔═╡ 5fc1116c-1adc-471a-a46b-fb4e3364b8fe
coords_icrs = coord_from_dist.(dist_samples)

# ╔═╡ 644580f0-e228-4368-bc21-fa952e9b51ff
let
	fig = Figure()
	ax = Axis(fig[1,1])

	scatter!(dist_samples, peris, markersize=1)

	fig
end

# ╔═╡ 16e882ab-8122-46d4-94d4-25efd84af8cd
let
	fig = Figure()
	ax = Axis(fig[1,1])

	scatter!(dist_samples, radii(L0), markersize=1)

	fig
end

# ╔═╡ 2f6fac94-db27-45d4-a2ae-773a29aec8f7
let
	fig = Figure()
	ax = Axis(fig[1,1])

	scatter!(dist_samples, v0, markersize=1)

	fig
end

# ╔═╡ 7253deeb-773b-4f57-9b80-d1ac6e1e2391


# ╔═╡ d93b9694-0d25-44cd-a01a-4196c0ca4216
let
	fig = Figure()
	ax = Axis(fig[1,1])

	scatter!(dist_samples, ϵ0, markersize=1)

	fig
end

# ╔═╡ c9a2ca22-93a0-4469-9362-60bc7ea02f31
Ltot = radii(L0)

# ╔═╡ 986df19d-1a96-425f-8f50-73092a4b5d44
r0 = radii(x0)

# ╔═╡ a83b96b0-7dec-4745-9c63-1e9a5e29c6a5
ϕ_peri = Agama.potential(pot, peris' .* [0, 0, 1])

# ╔═╡ 2c000e19-d7c6-46fb-8d4f-2f8b8e4f10f4
E1 = @. ϵ0 - ϕ_peri  - 1/2*(Ltot/peris)^2

# ╔═╡ cc154b67-380a-4971-8b01-ac370c9a854e
hist(E1)

# ╔═╡ 72851dc5-5250-4ac8-a729-3cf182f504ef
pericentre(ϵ, L) = LilGuys.find_zero(r -> -ϵ + Agama.potential(pot, r*[0,0,1])
	+ 1/2 * (L/r)^2, 0.1)

# ╔═╡ 5547d7e0-c922-4dda-a645-de94278bc8ad
peris_pred = pericentre.(ϵ0, Ltot)

# ╔═╡ ad615220-8b74-489c-a770-fff8709dc358
orbit_mean = LilGuys.agama_orbit(pot, icrs0, timerange=(0, -10/T2GYR))

# ╔═╡ fb5fb33b-0b16-417d-aa0e-2e7aa5046572
ϕ_mean = Agama.potential(pot, orbit_mean.positions[:, 1])

# ╔═╡ 945a4f74-befe-4d81-9195-9d7d74a5cc07
v_mean = radii(orbit_mean.velocities[:, 1])

# ╔═╡ c59455ac-0198-4cab-9555-f51527b5b8ad
L_mean = LilGuys.angular_momenta(orbit_mean)[:, 1]

# ╔═╡ 19b55549-3a93-4409-a38a-9a4fe50f1000
ϵ_mean = ϕ_mean + 1/2 * v_mean^2

# ╔═╡ 7a51abe9-d8a8-43ec-8e11-d86b86aa3afe
peris_pred_L = pericentre.(ϵ0, radii(L_mean))

# ╔═╡ bd973afd-3c3f-4024-a8ee-04fe0de6f8e9
peris_pred_L_big = pericentre.(ϵ0, quantile(Ltot, 0.95))

# ╔═╡ 2ca0d9d3-eb57-49e1-80ca-39b7b230ca09
dist_mean = icrs0.distance

# ╔═╡ b87480a4-b10e-423f-9e04-d8b7b4381173
v_pmra* LilGuys.pm2kms.(icrs0.pmra, dist_mean) .+ v_pmdec* LilGuys.pm2kms.(icrs0.pmdec, dist_mean)

# ╔═╡ 7ed39b1e-f4c7-4ec1-8c85-b92c938da41b
peris_pred_L2 = pericentre.(ϵ0, radii(L_mean) * (dist_samples ./ dist_mean) .^ 2)

# ╔═╡ 8b1a593d-6d98-426f-b31d-1e8e3fc6f8fe


# ╔═╡ 4443c98e-d74b-4d2b-b3b6-13c5ee2ce1eb
calc_r_circ(ϵ) = LilGuys.find_zero(r -> -ϵ + Agama.potential(pot, [0, 0, r]) + 1/2 * Agama.circular_velocity(pot, r)^2/2, 10)

# ╔═╡ 87ebacc6-b671-4cb4-966b-330c97f362a2
function calc_L_max(ϵ)
	r_c = calc_r_circ(ϵ)
	v_c = Agama.circular_velocity(pot, r_c)

	return r_c * v_c
end

# ╔═╡ 5a7f4cd0-ba2e-4c67-84d2-6e3225fccc3f
peris_pred_ϵ_sub = pericentre.(ϵ_mean, min.(calc_L_max(ϵ_mean) , Ltot))

# ╔═╡ 47fe025b-8688-480e-96cd-41aead935699
L_rel_mean = radii(L_mean) / calc_L_max(ϵ_mean)

# ╔═╡ 3dc00654-006d-4792-a0d9-81fbabc05f7a
L_rel = Ltot ./ calc_L_max.(ϵ0)

# ╔═╡ e4b89156-f449-4ea6-86c6-dd605e7ab8ae
peris_pred_ϵ = pericentre.(ϵ_mean, L_rel ./ L_rel_mean .* radii(L_mean) )

# ╔═╡ 53f4019f-95d2-4331-a556-f00480a09ef0
ϵ0

# ╔═╡ 0b571cb3-dadc-4a7d-8501-ab7aca3627f7
L_naive =radii(L_mean) * (dist_samples ./ dist_mean) .^ 3

# ╔═╡ 8c80ccf0-d5b8-424b-a194-813c487f1fee
coords_gsr = LilGuys.transform.(GSR, coord_from_dist.(dist_samples))

# ╔═╡ 1d08792f-421d-42aa-9e64-7fc9a1b78793
let
	fig = Figure()

	ax = Axis(fig[1,1])

	lines!(dist_samples, [c.pmra for c in coords_gsr] )

	ax = Axis(fig[1,2])
	lines!(dist_samples, [c.pmdec for c in coords_gsr])


	fig
end

# ╔═╡ 4919ebfb-7064-460f-ab29-ba6c30f420cc
v_tan_2 = (Ltot ./ r0)

# ╔═╡ 5882b5e9-6a5f-4156-a582-981d56ed7e2b
scatter(dist_samples, L0[1, :])

# ╔═╡ 28c89637-3dd0-4bc7-9eb9-ac33499bc52a
scatter(dist_samples, L0[2, :])

# ╔═╡ aab8a240-617e-407d-8d76-af763479ea6e
scatter(dist_samples, L0[3, :])

# ╔═╡ e3bd7c81-9a06-4917-bcbc-3738cb4c1b1f
⊕(a, b) = sqrt(a^2 + b^2)

# ╔═╡ df919036-7a85-44b2-8964-0e17605860d9
v_tan = [LilGuys.pm2kms(c.pmra ⊕ c.pmdec, c.distance) for c in coords_gsr] ./ V2KMS

# ╔═╡ d03e507c-3187-4ed3-8ae8-475c100edeb5
scatter(dist_samples, v_tan ./ v_tan_2)

# ╔═╡ 59c5229f-17fc-4122-b0f9-a713ecc5781d
scatter(dist_samples, v_tan ./ dist_samples .^ 2)

# ╔═╡ 4a1a7516-3cbe-4d32-9a2a-c63a054d8d21
scatter(dist_samples, L_naive ./ Ltot)

# ╔═╡ a45113ab-9439-4e67-b358-e84a23788d92
calc_r_circ(-2.5)

# ╔═╡ 8cb41e63-9efa-4cb6-9f9d-687a01e82935
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xlabel = "helcen dist / kpc",
			 ylabel = "pericentre / kpc")
	
	
	scatter!(dist_samples, peris_pred_ϵ_sub, label="fixed energy")
	scatter!(dist_samples, peris_pred_L, label="fixed L")
	# scatter!(dist_samples, peris_pred_L_big, label="fixed L")

	scatter!(dist_samples, peris_pred)
	lines!(dist_samples,  peris, color=:black, label="actual orbits")

	axislegend(position=:lt)

	fig
end

# ╔═╡ Cell order:
# ╠═9240885e-6f30-4f4b-8904-fe43f9297ef3
# ╠═92946976-38f5-11f1-813f-833fa39bd393
# ╠═6eb088f3-66b0-4f0b-ac84-e6b88c8d8ed1
# ╠═eeb6217f-6530-40e3-a313-12d2f6c07b60
# ╠═ffc0957b-ffd9-4594-b9e9-85aeaf2e48ec
# ╟─5e0b2289-e57f-48cf-a701-10d66dba1925
# ╠═ec309258-c5de-40ad-abf2-01d6808163d2
# ╠═cb89b58b-48d3-4269-8158-ecbc81464e9d
# ╠═6931ed06-99a2-4a3f-ac63-7dc7f4d67ef7
# ╟─c366985d-07a9-4ca3-85a8-ed36298cd501
# ╠═e3fedce4-eb08-4242-beaa-55a6b714f3c0
# ╠═5ad2202b-ad37-4c1c-9f0a-0a83da6fdc67
# ╠═2b90d979-58a9-4f1f-9813-1193ace020ad
# ╠═7d688254-2008-43f2-903c-f6d160e7aa0f
# ╠═1e3b7236-7c35-497d-9157-e68d15898535
# ╠═2b8b0ee9-0111-42f4-90ba-44cd7bc433e3
# ╠═402e3945-d039-4fe7-a6ba-571e9be908c6
# ╠═df363bdc-ac69-4725-bde6-81c5377c99ca
# ╠═807628e0-83af-431f-8014-a34387b883d5
# ╟─f5fe9855-0e8e-4fb7-9f46-0aefb97a1aa5
# ╠═2f7d73bf-8510-4d32-bedc-7e9ff84171d8
# ╠═1e8efe39-25d6-4b32-be58-1e6b6f2e31be
# ╠═ba736ef1-2c1d-4091-af50-9f8ba45961b0
# ╠═e061cbeb-7a7d-4765-873a-3bda33931af7
# ╠═64fe2dce-a020-4d7f-a546-19753136dae6
# ╠═53465bc5-dca7-44a2-b137-49876fe054a7
# ╠═d63fcaa5-4e5f-4a8b-b714-3d026f0dff30
# ╠═d78c9fe3-4070-4cac-93f1-5a92cfacfbcb
# ╟─2fe531fd-9724-4e4e-997e-ba308d2f7005
# ╠═2615f2e5-f7b6-453d-8471-8c4a7fa3b566
# ╠═0fff09cc-a1fc-43e8-bc8b-1cb2e32d78fb
# ╠═9d8add14-2706-4197-b99e-4060afaa00c7
# ╠═27215e5e-52f1-4f57-8f15-982887cfc71e
# ╠═ff51107b-579c-4a9a-9975-698fa1bd3bd7
# ╠═2de6de24-f317-4faa-bf5c-2dc93736761e
# ╠═1606b369-d217-45bd-80b4-c0328518a0dd
# ╠═06ac66d4-dcf8-4376-8b2e-23b70e5c89a0
# ╠═4bbe30af-1059-45aa-a6be-0000b0c48a04
# ╠═0ca9c4c6-7ee6-42f2-b9a6-de7d3149a0e9
# ╠═784d2f59-adf4-4564-9f29-9e2f1727ffa1
# ╠═dc6b7ad1-7bcd-4462-8e03-e83f935a6daf
# ╠═a74c7575-70e7-4301-b05f-3b5dee5e0e0c
# ╠═069efc5e-99ca-4cc1-9c94-75c5e67f80e3
# ╠═4333f006-4f9f-4d9a-8eca-d6d7d7974cdd
# ╠═c3b20d5e-5186-47e8-b32b-ad2dc49b6693
# ╠═b9423067-a382-48eb-afe3-3a617d8ac475
# ╠═784b9c5e-7705-499d-a593-c47f48ca0899
# ╠═b744f1db-f69b-4b19-a14d-a76c210ff305
# ╠═ae4f0dfd-aa6d-4f93-a0a1-f0c03e4d341c
# ╠═b87480a4-b10e-423f-9e04-d8b7b4381173
# ╠═20f32802-7793-4135-8411-7db1db8bc2b3
# ╠═5fc1116c-1adc-471a-a46b-fb4e3364b8fe
# ╠═1d08792f-421d-42aa-9e64-7fc9a1b78793
# ╠═644580f0-e228-4368-bc21-fa952e9b51ff
# ╠═16e882ab-8122-46d4-94d4-25efd84af8cd
# ╠═2f6fac94-db27-45d4-a2ae-773a29aec8f7
# ╠═7253deeb-773b-4f57-9b80-d1ac6e1e2391
# ╠═d93b9694-0d25-44cd-a01a-4196c0ca4216
# ╠═c9a2ca22-93a0-4469-9362-60bc7ea02f31
# ╠═986df19d-1a96-425f-8f50-73092a4b5d44
# ╠═a83b96b0-7dec-4745-9c63-1e9a5e29c6a5
# ╠═2c000e19-d7c6-46fb-8d4f-2f8b8e4f10f4
# ╠═cc154b67-380a-4971-8b01-ac370c9a854e
# ╠═72851dc5-5250-4ac8-a729-3cf182f504ef
# ╠═5547d7e0-c922-4dda-a645-de94278bc8ad
# ╠═ad615220-8b74-489c-a770-fff8709dc358
# ╠═fb5fb33b-0b16-417d-aa0e-2e7aa5046572
# ╠═945a4f74-befe-4d81-9195-9d7d74a5cc07
# ╠═c59455ac-0198-4cab-9555-f51527b5b8ad
# ╠═19b55549-3a93-4409-a38a-9a4fe50f1000
# ╠═7a51abe9-d8a8-43ec-8e11-d86b86aa3afe
# ╠═bd973afd-3c3f-4024-a8ee-04fe0de6f8e9
# ╠═2ca0d9d3-eb57-49e1-80ca-39b7b230ca09
# ╠═7ed39b1e-f4c7-4ec1-8c85-b92c938da41b
# ╠═e4b89156-f449-4ea6-86c6-dd605e7ab8ae
# ╠═5a7f4cd0-ba2e-4c67-84d2-6e3225fccc3f
# ╠═47fe025b-8688-480e-96cd-41aead935699
# ╠═3dc00654-006d-4792-a0d9-81fbabc05f7a
# ╠═87ebacc6-b671-4cb4-966b-330c97f362a2
# ╠═8b1a593d-6d98-426f-b31d-1e8e3fc6f8fe
# ╠═4443c98e-d74b-4d2b-b3b6-13c5ee2ce1eb
# ╠═53f4019f-95d2-4331-a556-f00480a09ef0
# ╠═0b571cb3-dadc-4a7d-8501-ab7aca3627f7
# ╠═8c80ccf0-d5b8-424b-a194-813c487f1fee
# ╠═df919036-7a85-44b2-8964-0e17605860d9
# ╠═4919ebfb-7064-460f-ab29-ba6c30f420cc
# ╠═d03e507c-3187-4ed3-8ae8-475c100edeb5
# ╠═5882b5e9-6a5f-4156-a582-981d56ed7e2b
# ╠═28c89637-3dd0-4bc7-9eb9-ac33499bc52a
# ╠═aab8a240-617e-407d-8d76-af763479ea6e
# ╠═e3bd7c81-9a06-4917-bcbc-3738cb4c1b1f
# ╠═59c5229f-17fc-4122-b0f9-a713ecc5781d
# ╠═4a1a7516-3cbe-4d32-9a2a-c63a054d8d21
# ╠═a45113ab-9439-4e67-b358-e84a23788d92
# ╠═8cb41e63-9efa-4cb6-9f9d-687a01e82935
