### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ a4fc4b04-a211-4744-ba2e-a1da88a11e52
begin
	import Pkg; Pkg.activate()

	import Agama
	using Arya, CairoMakie

	using LilGuys
end

# ╔═╡ f8d52cf1-d6dc-4290-8434-fc9385e2b1ad
import TOML

# ╔═╡ 71002bea-e657-4d14-9d17-2564bbaff92e
Arya.update_figsize!(4)

# ╔═╡ dbe50390-6cd4-4287-8cea-9c050f8c19e1
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/bootes3/observed_properties.toml"))

# ╔═╡ 927389fb-ce7c-4f0b-ad17-884ee067ab34
N = 1000

# ╔═╡ 725770d5-03b4-40eb-a747-e7e9622ff308
coords_i = LilGuys.transform.(LilGuys.Galactocentric, LilGuys.rand_coords(obs_props, N))

# ╔═╡ 83e84e46-519a-43ef-b576-93d250fb934e
tmax = -5/T2GYR

# ╔═╡ ae7f9bca-fa3d-463c-ab13-49f0cb00e565
function get_potential(potname; kwargs...)
	Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/", potname * ".ini"); kwargs...)
end

# ╔═╡ 0345c9a9-20ad-499f-9715-158ba7817415
function plot_orbits(orbits::AbstractVector{<:Orbit}...; N=100)
	N = min(length(orbits[1]), N)

	fig = Figure()
	axes = LilGuys.axis_xyz(fig, limits=LilGuys.limits_xyz(orbits[1][1].positions))
	legend = false
	
	for (i, orbit) in enumerate(orbits)
		
		pos = (LilGuys.positions.(orbit)[1:N])

		LilGuys.plot_xyz!(axes, pos..., alpha=0.1, color=COLORS[i])
	end

	fig
end

# ╔═╡ 04152ab7-a6e3-4763-a464-533d68f7d8ac
function plot_orbits(orbits::Pair...; N=100)
	N = min(length(orbits[1].second), N)

	fig = Figure()
	axes = LilGuys.axis_xyz(fig, limits=LilGuys.limits_xyz(orbits[1].second[1].positions))
	legend = false
	
	for (i, (label, orbit)) in enumerate(orbits)
		pos = (LilGuys.positions.(orbit)[1:N])
		LilGuys.plot_xyz!(axes, pos..., alpha=0.1, color=COLORS[i], label=label=>(; alpha=0.5))
	end

	Legend(fig[1,2], axes[1], tellwidth=false, merge=true, unique=true)

	fig
end

# ╔═╡ 9a768080-29cf-4bd4-b040-7e58d448a041
function plot_peris(orbits::Pair...)

	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "pericentre / kpc")
	
	for (i, (label, orbit)) in enumerate(orbits)
		peris = minimum.(radii.(LilGuys.positions.(orbit)))
		stephist!(peris,  label=label)
	end

	Legend(fig[1,2], ax, tellwidth=false, merge=true, unique=true)

	fig
end

# ╔═╡ 545976f9-c846-46f4-b762-ea54af1985f6
function plot_apos(orbits::Pair...)

	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "apocentre / kpc")
	
	for (i, (label, orbit)) in enumerate(orbits)
		peris = maximum.(radii.(LilGuys.positions.(orbit)))
		stephist!(peris,  label=label)
	end

	Legend(fig[1,2], ax, tellwidth=false, merge=true, unique=true)

	fig
end

# ╔═╡ b1ee881c-3f3e-4a1c-8fcf-96a46cf779ce


# ╔═╡ e42de723-34c2-4a50-8e8c-6b0ee3d198a9
CairoMakie.activate!(type=:png)

# ╔═╡ f5455167-36d9-4560-a11a-cd957cb6e2be
function radii_axis(gs)
	ax = Axis(gs, xlabel = "time / Gyr", ylabel = "galcen radius / kpc")

end

# ╔═╡ 5cc6bfee-ef97-4749-8417-990ca3e17df0
function orbits(pot; kwargs...)
	LilGuys.agama_orbit(pot, coords_i; timerange=(0, tmax), N=1001, kwargs...)
end

# ╔═╡ c6ed2e5b-da5e-4206-9ad6-aeef24f5d64e
md"""
## EP2020
"""

# ╔═╡ 8cc8db0a-a126-4a7a-9191-48e4c80a8f4a
import LinearAlgebra

# ╔═╡ 7be90896-9b51-4ba0-b51c-31fd12c8831e
to_sym_mat(x) = [x[1] x[4] x[6] 
				x[4] x[2] x[5]
				x[6] x[5] x[3]
]

# ╔═╡ f3c28010-649d-4bcf-8bba-f8a24b8532a8
function max_tidal_force(pot, positions)
	T = Agama.stress(pot, positions)
	return  maximum.(LinearAlgebra.eigvals.(to_sym_mat.(eachcol(T))))
end

# ╔═╡ 5c612348-aef7-4a47-bf41-e33d21820150
pot_ep2020 = get_potential("EP2020")

# ╔═╡ 742d3997-528f-4da1-8d26-68e879dedb00
orbits_ep2020 = orbits(pot_ep2020)

# ╔═╡ fb7a5392-21cb-41f0-b32f-c61b00b72b82
function plot_radii_time!(orbits; color=COLORS[1], N=100, kwargs...)
	N =  min(length(orbits), N)
	pos = (LilGuys.positions.(orbits_ep2020)[1:N])
	for orbit in orbits[1:N]
		lines!(orbit.times * T2GYR, radii(orbit.positions); alpha=0.1, color=color, kwargs...)
	end
	
end

# ╔═╡ cedd460e-fe60-46ab-b6c6-47ae9494b0ea
function plot_radii(orbitss...; kwargs...)
    fig = Figure()
	ax = radii_axis(fig[1,1])
	legend = false
	for i in eachindex(orbitss)
		orbits = orbitss[i]
		if orbits isa AbstractVector{<:Orbit}
			label = ""
		else
			label, orbits = orbits
			legend = true
		end
		
		plot_radii_time!(orbits; color=COLORS[i], label=label=>(; alpha=0.5), kwargs...)
	end

	if legend
		Legend(fig[1,2], ax, merge=true, unique=true)
	end
	
	fig
end

# ╔═╡ 50331fc2-1e1e-48ca-8dce-3e4bf6f8c0d7
let
	fig = Figure()
	ax = Axis(fig[1,1])
		
	
	for i in 1:100
		T = max_tidal_force(pot_ep2020, orbits_ep2020[i].positions)
		lines!(radii(orbits_ep2020[i]), T, alpha=0.1, color=COLORS[1])
	end

	fig
end

# ╔═╡ 073023eb-c70e-4b40-aa1a-649a9dbd5c0a
let
	fig = Figure()
	ax = Axis(fig[1,1])
		
	
	for i in 1:100
		T = max_tidal_force(pot_ep2020, orbits_ep2020[i].positions)
		scatter!(minimum(radii(orbits_ep2020[i])), maximum(T), alpha=1, color=COLORS[1])
	end

	fig
end

# ╔═╡ 4fe32d31-d088-46d0-9d34-b53945f1af10
plot_orbits(orbits_ep2020)

# ╔═╡ 59c83d03-1c45-45db-bc13-27bf377d76d6
plot_radii(orbits_ep2020)

# ╔═╡ 75e2c575-020d-4c24-be9c-549cc98a48d2
md"""
# LMC
"""

# ╔═╡ 5249c3db-9abf-4e16-b5b9-7892e0a75fb6
md"""
## V21
"""

# ╔═╡ 2f620f98-4d3c-436e-a5cc-0292affdcee7
v21 = get_potential("vasiliev+21/potential_evolving")

# ╔═╡ c1a8be9a-203f-4c25-b80a-230132843224
v21_no = get_potential("vasiliev+21/potential_nolmc")

# ╔═╡ 19b2880d-992c-4d84-8695-bf53a105b50d
md"""
## V24
"""

# ╔═╡ ff217db6-dd5f-442c-8336-80006a304537
vasiliev_units = Agama.AgamaUnits(velocity_scale=1/V2KMS, )

# ╔═╡ f5c3d1f4-b4ff-440e-ab50-b8d101aaef0f
orbits_v21 = orbits(v21, agama_units=vasiliev_units)

# ╔═╡ 8b40ff21-f802-42ab-8271-6017218869a3
orbits_v21_no  = orbits(v21_no, agama_units=vasiliev_units)

# ╔═╡ bb6fb6c2-7624-4737-8dec-8132a602485a
plot_radii("MW" => orbits_v21_no, "MW + LMC" => orbits_v21)

# ╔═╡ 7c438819-0797-4161-aa9d-202fa38af673
plot_orbits("MW" => orbits_v21_no, "MW + LMC" => orbits_v21)

# ╔═╡ 2bb6cd6c-324c-47ee-b945-b426a30acc54
Agama.velocity_scale(vasiliev_units), Agama.length_scale(vasiliev_units), Agama.time_scale(vasiliev_units)

# ╔═╡ 790f8b06-3d42-4fef-836f-aee3a5804c68
tmax / Agama.time_scale(vasiliev_units)

# ╔═╡ d01e186f-243e-4273-8051-3a9c30c25583
v24_L3M11 = get_potential("vasiliev24/L3M11/potential")

# ╔═╡ f42b41fe-96a1-41bc-9cab-8dfbe0218057
v24_L3 = get_potential("vasiliev24/L3M11/potential_mw_init")

# ╔═╡ de57cd22-9937-48dd-a781-66431853ebd2
orbits_L3M11 = orbits(v24_L3M11, agama_units=vasiliev_units)

# ╔═╡ 911b180c-bbe1-4024-a1dd-c775effea7cc
Agama.time_scale(vasiliev_units), Agama.velocity_scale(vasiliev_units)

# ╔═╡ 16523ea3-d0f7-481a-a4f0-1d64b6c10d26
orbits_L3 = orbits(v24_L3, agama_units=vasiliev_units)

# ╔═╡ a58cd4e0-f445-4892-94fc-e8feeef18ba9
plot_orbits("MW"=>orbits_L3, "MW+LMC"=>orbits_L3M11)

# ╔═╡ 47c2bd5a-0a05-44aa-b6bb-41fdaebaa34a
v24_L2M10 = get_potential("vasiliev24/L2M10/potential")

# ╔═╡ 4e8d47f5-2de8-4945-8fc4-146892ffc4e2
orbits_L2M10 = orbits(v24_L2M10, agama_units=vasiliev_units)

# ╔═╡ d2ae04f7-d461-4238-88fc-b9ab79576d99
plot_orbits(orbits_L2M10)

# ╔═╡ 121cab37-3805-49c4-91a5-018b6c3f9fa5
plot_radii(orbits_ep2020, orbits_L3)

# ╔═╡ 4073ceeb-7bf7-4f07-ba2a-d72129a40eef
plot_radii("MW" => orbits_L3, "MW + LMC" => orbits_L3M11)

# ╔═╡ ad204806-782a-48f7-a28a-7617e365f39b
md"""
# Bar & Spiral arms
"""

# ╔═╡ 355e9d3f-8ded-4dcc-b836-e2d30bfa878f
md"""
## Sormani / Portali model
"""

# ╔═╡ 56d31ad9-9682-49bf-a5dc-eae995b7d361
md"""
This potential doesn't work well since it is only valid in the innermost galaxy :/.
"""

# ╔═╡ 237e5f17-2b48-4a27-98d3-c255715368cd
portali_units = Agama.AgamaUnits(velocity_scale=1/V2KMS)

# ╔═╡ f1e45a61-c447-4d50-aeb4-81d8d62fbfe2
4.301476e-6

# ╔═╡ 490edd58-75c8-4aa4-9aa1-5c55189580c4
pot_sormani = get_potential("sormani+2022/portali+2017")

# ╔═╡ eb61c6a2-889e-490f-9949-7e6564e00d73
pot_sormani_no = get_potential("sormani+2022/portali+2017_axi")

# ╔═╡ de475606-25cd-425e-baad-605902bb332e
orbits_nobar = orbits(pot_sormani_no, agama_units=portali_units)

# ╔═╡ 0a1f6962-0dba-46a0-b51a-c221988812a8
Agama.enclosed_mass(pot_sormani, 100)

# ╔═╡ 5cfe4631-bca6-4b5d-aad3-2788d09654ce
md"""
## Hunter 24
"""

# ╔═╡ 47e9180e-3ad0-4040-9ceb-95e002717e1b
pot_spiral = get_potential("hunter+2024/MWPotentialHunter24_rotspiral")

# ╔═╡ f00d04de-df07-4fe3-b24c-e2cbc0d6ad1a
Agama.enclosed_mass(pot_spiral, 100)

# ╔═╡ a1c10775-52f3-4dfb-9d53-348c8bed9690
pot_nospiral = get_potential("hunter+2024/MWPotentialHunter24_axi")

# ╔═╡ aab2cda3-b66c-4c79-8789-b3da7382c64c
pot_bar = get_potential("hunter+2024/MWPotentialHunter24_rotating")

# ╔═╡ 21c9d9b7-2009-45dd-b66c-06896d1a91ca
orbits_spiral = orbits(pot_spiral, agama_units=portali_units)

# ╔═╡ b3a3ed2c-7676-4c69-a24c-c5f8457435b6
orbits_nospiral = orbits(pot_nospiral, agama_units=portali_units)

# ╔═╡ 375566f1-4c5c-470e-b84d-7f0be18a4ff4
#=╠═╡
plot_radii(orbits_nobar, orbits_bar)
  ╠═╡ =#

# ╔═╡ f423d31e-282c-4610-8ba0-bbc56f10f544
#=╠═╡
plot_orbits("axisym" => orbits_nospiral, "bar" => orbits_bar, "bar+spiral" => orbits_spiral, )
  ╠═╡ =#

# ╔═╡ fd726be8-9053-4d44-b8c8-1d185e70e3e7
#=╠═╡
plot_radii("axisym" => orbits_nospiral, "bar" => orbits_bar, "bar+spiral" => orbits_spiral,)
  ╠═╡ =#

# ╔═╡ 9348378d-1712-4ba0-b583-1aca6f7784be
Agama.enclosed_mass(pot_spiral, 200)

# ╔═╡ 09cbd534-7082-4e0c-997e-34026be33004
Agama.enclosed_mass(v24_L3, 200)

# ╔═╡ 9042e75e-92eb-4361-9ef5-9912fef4ec5c
md"""
# Combined
"""

# ╔═╡ 80049420-d2e2-4f9d-8f20-910266d70d91
bar_angle = deg2rad(-25 )

# ╔═╡ a08a50d9-a488-477f-b811-3acadce90ea4
omega_bar = -37.5

# ╔═╡ c38417fa-1d36-43c3-b70d-08a6d0b1227c
pot_spiral_only = get_potential("hunter+2024/MWPotentialHunter24_spiral", rotation=[[0, bar_angle], [1, bar_angle + omega_bar]])

# ╔═╡ 1c327af7-abed-4568-92c6-e4203ee738cd
pot_bar_only = get_potential("hunter+2024/MWPotentialHunter24_bar", rotation=[[0, bar_angle], [1, bar_angle + omega_bar]])

# ╔═╡ e2d68e0f-aa9e-4895-9572-ca87ec247134
Agama.enclosed_mass(pot_nospiral, 100)

# ╔═╡ 8a0311a9-bbf2-4a08-9e78-008185618182
pot_lmc_dm = get_potential("vasiliev24/L3M11/potential_dm")

# ╔═╡ 66e2e47b-d077-461c-bc09-356336d6d466
pot_bary = get_potential("hunter+2024/MWPotentialHunter24_bary")

# ╔═╡ dc1a907c-08a5-4a18-bee7-0a9e8e8ca5b8
Agama.enclosed_mass(pot_bar_only, 100), Agama.enclosed_mass(pot_bary, 10), Agama.enclosed_mass(pot_spiral_only, 100)

# ╔═╡ 1fd94ae2-3ff4-4fa5-a6b9-ba82cbc9c729
pot_barspiral_lmc = pot_bar_only + pot_spiral_only + pot_bary + pot_lmc_dm

# ╔═╡ 52720a6d-9fde-4d75-8f6e-1536537490c1
Agama.circular_velocity(pot_barspiral_lmc, 8)

# ╔═╡ bffedda4-a70f-4807-9075-8ab086d210da
Agama.circular_velocity(pot_spiral, 8)

# ╔═╡ 87c9b1d5-e9e3-45cb-991d-1e2489ff3ae5
orbits_spiral_lmc = orbits(pot_barspiral_lmc, agama_units=vasiliev_units)

# ╔═╡ 998d2699-bf50-48b2-819f-ad63d49ebcf7
plot_radii("mw only" => orbits_L3, "mw+bar+spiral" => orbits_spiral,  "mw+lmc" => orbits_L3M11, "bar+spiral+lmc" => orbits_spiral_lmc, )

# ╔═╡ 590f4166-aad1-4ecc-b8b8-05b9e056df87
plot_orbits("mw only" => orbits_L3, "mw+bar+spiral" => orbits_spiral,  "mw+lmc" => orbits_L3M11, "bar+spiral+lmc" => orbits_spiral_lmc,)

# ╔═╡ d39f73ca-bca5-4dc4-a52d-97b12d1d261b
plot_peris("mw only" => orbits_L3, "mw+bar+spiral" => orbits_spiral,  "mw+lmc" => orbits_L3M11, "bar+spiral+lmc" => orbits_spiral_lmc,)

# ╔═╡ b12291ea-2090-4474-bf15-51f08788676b
plot_apos("mw only" => orbits_L3, "mw+bar+spiral" => orbits_spiral,  "mw+lmc" => orbits_L3M11, "bar+spiral+lmc" => orbits_spiral_lmc,)

# ╔═╡ 4c3be5af-5278-4c6f-9125-1de8f5494445
md"""
# Energies
"""

# ╔═╡ bd234d79-e8c9-4297-84e0-7997fb007d9b
function get_energies(pot, orbits; agama_units=Agama.current_units())

	Es = zeros(length(orbits))
	for (i, orbit) in enumerate(orbits)
		Φ = Agama.potential(pot, orbit.positions[:, 1], agama_units)
		V = LilGuys.radii(orbit.velocities[:, 1])

		Es[i] =  1/2 *V^2 + Φ
	end
	return Es
end

# ╔═╡ a17634af-c5fa-4cf5-8c8a-605336729c0f
function get_L(orbits; agama_units=Agama.current_units())

	Lz = zeros(3, length(orbits))
	for (i, orbit) in enumerate(orbits)
		Lz[:, i] = LilGuys.angular_momenta(orbit.positions[:, 1], orbit.velocities[:, 1])
	end
	return Lz
end

# ╔═╡ 3aff5898-4a6d-458e-abb1-fc576d95369e
orbits_ep2020[1].positions[1, :]

# ╔═╡ 0010c2cc-a7ed-4206-baf1-9e26b1f7f8e8
hist(get_energies(pot_ep2020, orbits_ep2020) * V2KMS^2 / 1e5)

# ╔═╡ 981b742c-75f4-45e6-99dc-8f74d1552fcf
hist(get_energies(v24_L3M11, orbits_L3M11, agama_units=vasiliev_units) * V2KMS^2 / 1e5, 
	axis=(; xlabel = "E / km^2 / s^2"))

# ╔═╡ f38a7281-d11e-4d4e-a48d-d7b54fad5d5c
hist(get_L(orbits_ep2020, agama_units=vasiliev_units)[3, :] * V2KMS^2, 
	axis=(; xlabel = "Lz / km^2 / s^2"))

# ╔═╡ 54e36834-5121-452c-bfd4-3bbcfcac0a62
# ╠═╡ disabled = true
#=╠═╡
orbits_bar = orbits(pot_sormani, agama_units=portali_units)
  ╠═╡ =#

# ╔═╡ 01596c58-082e-4d8b-a7ef-d3a85288f631
#=╠═╡
orbits_bar = orbits(pot_bar, agama_units=portali_units)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═a4fc4b04-a211-4744-ba2e-a1da88a11e52
# ╠═f8d52cf1-d6dc-4290-8434-fc9385e2b1ad
# ╠═71002bea-e657-4d14-9d17-2564bbaff92e
# ╠═dbe50390-6cd4-4287-8cea-9c050f8c19e1
# ╠═927389fb-ce7c-4f0b-ad17-884ee067ab34
# ╠═725770d5-03b4-40eb-a747-e7e9622ff308
# ╠═83e84e46-519a-43ef-b576-93d250fb934e
# ╠═ae7f9bca-fa3d-463c-ab13-49f0cb00e565
# ╠═0345c9a9-20ad-499f-9715-158ba7817415
# ╠═04152ab7-a6e3-4763-a464-533d68f7d8ac
# ╠═9a768080-29cf-4bd4-b040-7e58d448a041
# ╠═545976f9-c846-46f4-b762-ea54af1985f6
# ╠═b1ee881c-3f3e-4a1c-8fcf-96a46cf779ce
# ╠═e42de723-34c2-4a50-8e8c-6b0ee3d198a9
# ╠═fb7a5392-21cb-41f0-b32f-c61b00b72b82
# ╠═f5455167-36d9-4560-a11a-cd957cb6e2be
# ╠═cedd460e-fe60-46ab-b6c6-47ae9494b0ea
# ╠═5cc6bfee-ef97-4749-8417-990ca3e17df0
# ╠═c6ed2e5b-da5e-4206-9ad6-aeef24f5d64e
# ╠═8cc8db0a-a126-4a7a-9191-48e4c80a8f4a
# ╠═7be90896-9b51-4ba0-b51c-31fd12c8831e
# ╠═f3c28010-649d-4bcf-8bba-f8a24b8532a8
# ╠═50331fc2-1e1e-48ca-8dce-3e4bf6f8c0d7
# ╠═073023eb-c70e-4b40-aa1a-649a9dbd5c0a
# ╠═5c612348-aef7-4a47-bf41-e33d21820150
# ╠═742d3997-528f-4da1-8d26-68e879dedb00
# ╠═4fe32d31-d088-46d0-9d34-b53945f1af10
# ╠═59c83d03-1c45-45db-bc13-27bf377d76d6
# ╟─75e2c575-020d-4c24-be9c-549cc98a48d2
# ╟─5249c3db-9abf-4e16-b5b9-7892e0a75fb6
# ╠═2f620f98-4d3c-436e-a5cc-0292affdcee7
# ╠═c1a8be9a-203f-4c25-b80a-230132843224
# ╠═f5c3d1f4-b4ff-440e-ab50-b8d101aaef0f
# ╠═8b40ff21-f802-42ab-8271-6017218869a3
# ╠═bb6fb6c2-7624-4737-8dec-8132a602485a
# ╠═7c438819-0797-4161-aa9d-202fa38af673
# ╠═19b2880d-992c-4d84-8695-bf53a105b50d
# ╠═ff217db6-dd5f-442c-8336-80006a304537
# ╠═2bb6cd6c-324c-47ee-b945-b426a30acc54
# ╠═790f8b06-3d42-4fef-836f-aee3a5804c68
# ╠═d01e186f-243e-4273-8051-3a9c30c25583
# ╠═f42b41fe-96a1-41bc-9cab-8dfbe0218057
# ╠═de57cd22-9937-48dd-a781-66431853ebd2
# ╠═911b180c-bbe1-4024-a1dd-c775effea7cc
# ╠═16523ea3-d0f7-481a-a4f0-1d64b6c10d26
# ╠═a58cd4e0-f445-4892-94fc-e8feeef18ba9
# ╠═47c2bd5a-0a05-44aa-b6bb-41fdaebaa34a
# ╠═4e8d47f5-2de8-4945-8fc4-146892ffc4e2
# ╠═d2ae04f7-d461-4238-88fc-b9ab79576d99
# ╠═121cab37-3805-49c4-91a5-018b6c3f9fa5
# ╠═4073ceeb-7bf7-4f07-ba2a-d72129a40eef
# ╟─ad204806-782a-48f7-a28a-7617e365f39b
# ╟─355e9d3f-8ded-4dcc-b836-e2d30bfa878f
# ╟─56d31ad9-9682-49bf-a5dc-eae995b7d361
# ╠═237e5f17-2b48-4a27-98d3-c255715368cd
# ╠═f1e45a61-c447-4d50-aeb4-81d8d62fbfe2
# ╠═490edd58-75c8-4aa4-9aa1-5c55189580c4
# ╠═eb61c6a2-889e-490f-9949-7e6564e00d73
# ╠═54e36834-5121-452c-bfd4-3bbcfcac0a62
# ╠═de475606-25cd-425e-baad-605902bb332e
# ╠═0a1f6962-0dba-46a0-b51a-c221988812a8
# ╠═f00d04de-df07-4fe3-b24c-e2cbc0d6ad1a
# ╠═375566f1-4c5c-470e-b84d-7f0be18a4ff4
# ╟─5cfe4631-bca6-4b5d-aad3-2788d09654ce
# ╠═47e9180e-3ad0-4040-9ceb-95e002717e1b
# ╠═a1c10775-52f3-4dfb-9d53-348c8bed9690
# ╠═aab2cda3-b66c-4c79-8789-b3da7382c64c
# ╠═21c9d9b7-2009-45dd-b66c-06896d1a91ca
# ╠═b3a3ed2c-7676-4c69-a24c-c5f8457435b6
# ╠═01596c58-082e-4d8b-a7ef-d3a85288f631
# ╠═f423d31e-282c-4610-8ba0-bbc56f10f544
# ╠═fd726be8-9053-4d44-b8c8-1d185e70e3e7
# ╠═9348378d-1712-4ba0-b583-1aca6f7784be
# ╠═09cbd534-7082-4e0c-997e-34026be33004
# ╠═9042e75e-92eb-4361-9ef5-9912fef4ec5c
# ╠═80049420-d2e2-4f9d-8f20-910266d70d91
# ╠═a08a50d9-a488-477f-b811-3acadce90ea4
# ╠═c38417fa-1d36-43c3-b70d-08a6d0b1227c
# ╠═1c327af7-abed-4568-92c6-e4203ee738cd
# ╠═dc1a907c-08a5-4a18-bee7-0a9e8e8ca5b8
# ╠═e2d68e0f-aa9e-4895-9572-ca87ec247134
# ╠═8a0311a9-bbf2-4a08-9e78-008185618182
# ╠═1fd94ae2-3ff4-4fa5-a6b9-ba82cbc9c729
# ╠═66e2e47b-d077-461c-bc09-356336d6d466
# ╠═52720a6d-9fde-4d75-8f6e-1536537490c1
# ╠═bffedda4-a70f-4807-9075-8ab086d210da
# ╠═87c9b1d5-e9e3-45cb-991d-1e2489ff3ae5
# ╠═998d2699-bf50-48b2-819f-ad63d49ebcf7
# ╠═590f4166-aad1-4ecc-b8b8-05b9e056df87
# ╠═d39f73ca-bca5-4dc4-a52d-97b12d1d261b
# ╠═b12291ea-2090-4474-bf15-51f08788676b
# ╠═4c3be5af-5278-4c6f-9125-1de8f5494445
# ╠═bd234d79-e8c9-4297-84e0-7997fb007d9b
# ╠═a17634af-c5fa-4cf5-8c8a-605336729c0f
# ╠═3aff5898-4a6d-458e-abb1-fc576d95369e
# ╠═0010c2cc-a7ed-4206-baf1-9e26b1f7f8e8
# ╠═981b742c-75f4-45e6-99dc-8f74d1552fcf
# ╠═f38a7281-d11e-4d4e-a48d-d7b54fad5d5c
