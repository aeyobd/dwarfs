### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ a4fc4b04-a211-4744-ba2e-a1da88a11e52
begin
	import Pkg; Pkg.activate()

	import Agama
	using Arya, CairoMakie

	using LilGuys
end

# ╔═╡ 50edc139-077a-4e23-9c67-067fc133f16d
md"""
# Dynamical Friction
This notebook calculates the dynamical friction for dwarf spheroidals using a variety of different methods. 

Additionally, is a validation for galpy vs my agama + leapfrog formalism.

"""

# ╔═╡ 337d4229-53ff-4619-8a20-7f81cbd5fa4a
md"""
## Background
Dynamical friction is commonly approximated using the Chandrasekhar formalism

The help for dynamical friction below describes the current implementation.
"""

# ╔═╡ 058c1b2a-e461-4961-91fa-4ec0bfeb1b21
help(LilGuys.a_dyn_friction)

# ╔═╡ c23b938b-77d2-4267-9598-edb82e54f7bb
obs_i = ICRS(
	ra = 15.0183,
	dec = -33.7186,
	distance = 86,
	pmra = 0.14,
	pmdec = -0.151,
	radial_velocity = 111.2
)

# ╔═╡ 88310aa2-a470-4a58-abb9-6d4b9530a15b
begin 
	galaxy = "sculptor"
	M_tot = 0.15
	r_s = 1.5
	tmax = -10/T2GYR
end

# ╔═╡ 64695fe6-14b2-47f7-8c6f-1565f45ee34a
nfw = LilGuys.TruncNFW(M200=M_tot, r_s=r_s, trunc=20, xi=3)

# ╔═╡ 35cc230d-0f7b-4e10-9270-4244f29945dc
rhalf = LilGuys.r_h(nfw)

# ╔═╡ a702afe0-9c8a-46b5-aef2-dc5c907fc862
σv_0 = 150

# ╔═╡ 307b1eed-15a8-4d25-980c-6e7390c49510
md"""
Examples for Mtot and r_s are
- 0.15, 1.5 (sculptor)
- 0.29, 1.8 (ursa minor)
- 0.42, 3.2 (boo III?)
"""

# ╔═╡ a6c4bed5-f274-43d5-a18c-3222c3cda537
LilGuys.NFW(v_circ_max = 0.169, r_circ_max=6.9)

# ╔═╡ f8d52cf1-d6dc-4290-8434-fc9385e2b1ad
import TOML

# ╔═╡ 71002bea-e657-4d14-9d17-2564bbaff92e
Arya.update_figsize!(4)

# ╔═╡ dbe50390-6cd4-4287-8cea-9c050f8c19e1
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxy/observed_properties.toml"))

# ╔═╡ 927389fb-ce7c-4f0b-ad17-884ee067ab34
N = 1

# ╔═╡ 725770d5-03b4-40eb-a747-e7e9622ff308
coords_i = LilGuys.transform(LilGuys.Galactocentric, obs_i)

# ╔═╡ ae7f9bca-fa3d-463c-ab13-49f0cb00e565
function get_potential(potname; kwargs...)
	Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/", potname * ".ini"); kwargs...)
end

# ╔═╡ 0345c9a9-20ad-499f-9715-158ba7817415
function plot_orbits(orbits::AbstractVector{<:Orbit}...; N=100, alpha=1)
	N = min(length(orbits[1]), N)

	fig = Figure()
	axes = LilGuys.axis_xyz(fig, limits=LilGuys.limits_xyz(orbits[1][1].positions))
	legend = false
	
	for (i, orbit) in enumerate(orbits)
		
		pos = (LilGuys.positions.(orbit)[1:N])

		LilGuys.plot_xyz!(axes, pos..., alpha=alpha, color=COLORS[i])
	end

	fig
end

# ╔═╡ 04152ab7-a6e3-4763-a464-533d68f7d8ac
function plot_orbits(orbits::Pair...; N=100, alpha=1)
	N = min(length(orbits[1].second), N)

	fig = Figure()
	axes = LilGuys.axis_xyz(fig, limits=LilGuys.limits_xyz(orbits[1].second.positions))
	legend = false
	
	for (i, (label, orbit)) in enumerate(orbits)
		pos = (LilGuys.positions(orbit))
		LilGuys.plot_xyz!(axes, pos, alpha=alpha, color=COLORS[i], label=label=>(; alpha=max(0.5, alpha)))
	end

	Legend(fig[1,2], axes[1], tellwidth=false, merge=true, unique=true)

	fig
end

# ╔═╡ e42de723-34c2-4a50-8e8c-6b0ee3d198a9
CairoMakie.activate!(type=:png)

# ╔═╡ fb7a5392-21cb-41f0-b32f-c61b00b72b82
function plot_radii_time!(orbit; color=COLORS[1], N=100, kwargs...)
	lines!(orbit.times * T2GYR, radii(orbit.positions); color=color, kwargs...)
	
end

# ╔═╡ f5455167-36d9-4560-a11a-cd957cb6e2be
function radii_axis(gs)
	ax = Axis(gs, xlabel = "time / Gyr", ylabel = "galcen radius / kpc", limits=(nothing, nothing, 0, nothing))

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

# ╔═╡ 8cc8db0a-a126-4a7a-9191-48e4c80a8f4a
import LinearAlgebra

# ╔═╡ 7be90896-9b51-4ba0-b51c-31fd12c8831e
to_sym_mat(x) = [x[1] x[4] x[6] 
				x[4] x[2] x[5]
				x[6] x[5] x[3]
]

# ╔═╡ 2688df7b-075d-410a-b1ec-dc7c5a488d4b
md"""
# Comparison
"""

# ╔═╡ c6ed2e5b-da5e-4206-9ad6-aeef24f5d64e
md"""
## EP2020
"""

# ╔═╡ 5c612348-aef7-4a47-bf41-e33d21820150
pot_ep2020 = get_potential("vasiliev24/L3M11/potential")

# ╔═╡ c54baa68-26d8-4890-a5eb-0bde94a26133
pot_static = get_potential("vasiliev24/L3M11/potential_mw_init")

# ╔═╡ 464af32a-56ea-4707-87f6-c04d58ea459c
units = Agama.VASILIEV_UNITS

# ╔═╡ 5cc6bfee-ef97-4749-8417-990ca3e17df0
function orbits(pot; kwargs...)
	LilGuys.agama_orbit(pot, coords_i; agama_units=units, timerange=(0, tmax), N=1001, kwargs...)
end

# ╔═╡ f3c28010-649d-4bcf-8bba-f8a24b8532a8
function max_tidal_force(pot, positions)
	T = Agama.stress(pot, positions, units)
	return  maximum.(LinearAlgebra.eigvals.(to_sym_mat.(eachcol(T))))
end

# ╔═╡ 742d3997-528f-4da1-8d26-68e879dedb00
orbits_ep2020 = orbits(pot_ep2020)

# ╔═╡ bb4c9727-b8af-467b-b900-9c2a14231035
md"""
## Agama
"""

# ╔═╡ 3035f9a1-84f3-4ef7-84d5-3cb05b7d1668
function get_sigma_v(pot)
	gm = Agama.GalaxyModel(pot, Agama.DistributionFunction(pot, type="QuasiSpherical"))


	N = 100
	Rs = 10 .^ LinRange(-1, 3, N)

	sigmas = gm._py.moments(Agama.mat2py([zeros(N) Rs zeros(N)]'), dens=false) |> Agama.py2mat

	sigma3 = sigmas[1, :] .+ sigmas[2, :] + sigmas[3, :]

	sigma1 = sigma3 ./ 3
	return LilGuys.lerp(Rs, sqrt.(sigma1) .* Agama.velocity_scale(units))
end

# ╔═╡ b0887636-ad6c-43cc-b5a1-a76f678e8703
σv = get_sigma_v(pot_static)

# ╔═╡ 51ab5b01-ad58-42e1-b63c-3713bc2e8142
σv.([8, 100, 300]) .* V2KMS

# ╔═╡ c958051a-ae38-4284-8ecd-ce604ff5d23a
dyn_fric = LilGuys.ChandrashakarDynamicalFriction(r_s=r_s, σv=x->σv(radii(x)), M=M_tot, ρ = x->Agama.density(pot_static, x, units), Λ = exp(4))

# ╔═╡ 0c828c6d-76e3-42f8-9454-e787b3d8fac0
dyn_fric_moving = LilGuys.ChandrashakarDynamicalFriction(r_s=r_s, σv=x->σv(radii(x)), M=M_tot, ρ = x->Agama.density(pot_static, x, units))

# ╔═╡ 4b80bed0-65ef-46ff-b44b-086b0d563a81
LilGuys.density(pot::Agama.Potential, x) = Agama.density(pot, x)

# ╔═╡ 50f3c97d-4385-4cff-930c-a9d0c6cbe0c4
LilGuys.acceleration(pot::Agama.Potential, x) = Agama.acceleration(pot, x)

# ╔═╡ 22f5f57e-dcf0-4fda-b701-fd04b65a8b9f
f_tot = (pos, vel, t) -> (LilGuys.acceleration(dyn_fric, pos, vel) + Agama.acceleration(pot_ep2020, pos, units, t=t))

# ╔═╡ 3b681ed0-d7d5-4aab-a5fe-986f4b395825
orbits_dyn_fric = LilGuys.leapfrog(f_tot, coords_i, timerange=(0, tmax))

# ╔═╡ 2f3b22c5-aff6-4c30-8c23-77554463a6db
f_tot([100., 5, 0], [0., 0, 0.0000001], 0.0)

# ╔═╡ 08332449-06c2-4ced-b060-f08e1ca00971
f_tot_moving = (pos, vel, t) -> (LilGuys.acceleration(dyn_fric_moving, pos, vel) + Agama.acceleration(pot_ep2020, pos, units, t=t))

# ╔═╡ 45449bd6-140d-4009-82d7-54e387265957
orbits_dyn_fric_moving = LilGuys.leapfrog(f_tot_moving, coords_i, timerange=(0, tmax))

# ╔═╡ aaa1976b-2719-4a0d-b119-8f928b5a6003
orbits_list = [
	"no dyn fric" => orbits_ep2020, 
	"dyn fric" => orbits_dyn_fric, 
	"dyn fric moving" => orbits_dyn_fric_moving,

]

# ╔═╡ dbcc21e7-1b92-48b9-902e-e05b1436d5d6
plot_orbits(orbits_list...,)

# ╔═╡ 0632953f-c513-4e25-a14e-0e016ed180d3
plot_radii(orbits_list..., )

# ╔═╡ a6d1c4a2-3037-49e9-a209-e82f6eb734c4
all_orbits_list = [
	"no dyn fric" => orbits_ep2020, 
	"dyn fric" => orbits_dyn_fric, 
	"dyn fric moving" => orbits_dyn_fric_moving,
]

# ╔═╡ Cell order:
# ╟─50edc139-077a-4e23-9c67-067fc133f16d
# ╟─337d4229-53ff-4619-8a20-7f81cbd5fa4a
# ╠═058c1b2a-e461-4961-91fa-4ec0bfeb1b21
# ╠═c23b938b-77d2-4267-9598-edb82e54f7bb
# ╠═64695fe6-14b2-47f7-8c6f-1565f45ee34a
# ╠═35cc230d-0f7b-4e10-9270-4244f29945dc
# ╠═88310aa2-a470-4a58-abb9-6d4b9530a15b
# ╠═a702afe0-9c8a-46b5-aef2-dc5c907fc862
# ╠═307b1eed-15a8-4d25-980c-6e7390c49510
# ╠═a6c4bed5-f274-43d5-a18c-3222c3cda537
# ╠═a4fc4b04-a211-4744-ba2e-a1da88a11e52
# ╠═f8d52cf1-d6dc-4290-8434-fc9385e2b1ad
# ╠═71002bea-e657-4d14-9d17-2564bbaff92e
# ╠═dbe50390-6cd4-4287-8cea-9c050f8c19e1
# ╠═927389fb-ce7c-4f0b-ad17-884ee067ab34
# ╠═725770d5-03b4-40eb-a747-e7e9622ff308
# ╠═ae7f9bca-fa3d-463c-ab13-49f0cb00e565
# ╠═0345c9a9-20ad-499f-9715-158ba7817415
# ╠═04152ab7-a6e3-4763-a464-533d68f7d8ac
# ╠═e42de723-34c2-4a50-8e8c-6b0ee3d198a9
# ╠═fb7a5392-21cb-41f0-b32f-c61b00b72b82
# ╠═f5455167-36d9-4560-a11a-cd957cb6e2be
# ╠═cedd460e-fe60-46ab-b6c6-47ae9494b0ea
# ╠═5cc6bfee-ef97-4749-8417-990ca3e17df0
# ╠═8cc8db0a-a126-4a7a-9191-48e4c80a8f4a
# ╠═7be90896-9b51-4ba0-b51c-31fd12c8831e
# ╠═f3c28010-649d-4bcf-8bba-f8a24b8532a8
# ╟─2688df7b-075d-410a-b1ec-dc7c5a488d4b
# ╠═aaa1976b-2719-4a0d-b119-8f928b5a6003
# ╠═a6d1c4a2-3037-49e9-a209-e82f6eb734c4
# ╠═dbcc21e7-1b92-48b9-902e-e05b1436d5d6
# ╠═0632953f-c513-4e25-a14e-0e016ed180d3
# ╟─c6ed2e5b-da5e-4206-9ad6-aeef24f5d64e
# ╠═5c612348-aef7-4a47-bf41-e33d21820150
# ╠═c54baa68-26d8-4890-a5eb-0bde94a26133
# ╠═464af32a-56ea-4707-87f6-c04d58ea459c
# ╠═742d3997-528f-4da1-8d26-68e879dedb00
# ╟─bb4c9727-b8af-467b-b900-9c2a14231035
# ╠═3035f9a1-84f3-4ef7-84d5-3cb05b7d1668
# ╠═b0887636-ad6c-43cc-b5a1-a76f678e8703
# ╠═51ab5b01-ad58-42e1-b63c-3713bc2e8142
# ╠═c958051a-ae38-4284-8ecd-ce604ff5d23a
# ╠═0c828c6d-76e3-42f8-9454-e787b3d8fac0
# ╠═22f5f57e-dcf0-4fda-b701-fd04b65a8b9f
# ╠═08332449-06c2-4ced-b060-f08e1ca00971
# ╠═3b681ed0-d7d5-4aab-a5fe-986f4b395825
# ╠═45449bd6-140d-4009-82d7-54e387265957
# ╠═2f3b22c5-aff6-4c30-8c23-77554463a6db
# ╠═4b80bed0-65ef-46ff-b44b-086b0d563a81
# ╠═50f3c97d-4385-4cff-930c-a9d0c6cbe0c4
