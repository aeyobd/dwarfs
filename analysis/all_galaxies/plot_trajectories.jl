### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ 71389566-8f26-11f0-1805-853166b6d451
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
	
	using HDF5
	using DataFrames
	import CSV
end

# ╔═╡ 13eeedb1-814b-4682-8830-f7f86458118f
using Printf

# ╔═╡ ac9c37c2-a51a-4a10-81b3-2dfe0f7ab11e
import Agama

# ╔═╡ 1f097dca-e93f-4f08-9210-14a1356782f0
CairoMakie.activate!(type=:png)

# ╔═╡ e8cec7c0-9f08-47ed-a26e-6016fe382ffb
units = Agama.VASILIEV_UNITS

# ╔═╡ 903503e8-598a-441e-bf88-7d0a84d221ef
mw_halo = NFW(r_circ_max=43.7, v_circ_max=191/V2KMS)

# ╔═╡ 541a30dd-992c-4078-afa9-e1cb53dfc4a3
md"""
# Data loading
"""

# ╔═╡ 51141c7d-59e7-4df0-8928-1e83238fdcff
function load_trajectories(modelname)
	local pos, vel, t
	h5open(joinpath(modelname, "trajectory.hdf5")) do h5 
	
		pos = h5["positions"][:, :, :]
		vel = h5["velocities"][:, :, :]
		t = h5["times"][:]
	end

	return [Orbit(positions=pos[:, i, :], velocities=vel[:, i, :], times=t) for i in 1:size(pos, 2)]
end

# ╔═╡ 2ca5f1ae-37b5-4a96-b890-eae1ae4609c7
traj = load_trajectories("L3M11_hhr")

# ╔═╡ 0bcbb650-e54e-4a3e-86ff-7ab4b9504b6a
traj_alt = load_trajectories("L3M11_hr")

# ╔═╡ e8d209a7-fe07-4ca9-9151-84dd27a75023
traj_no = load_trajectories("M11")

# ╔═╡ bafed0b5-26b5-446c-942c-40634bc1102c
[LilGuys.positions(o)[:, 1]  for o in (traj .- traj_alt)]

# ╔═╡ cc25d8ef-b783-4826-92fa-b1bc0d3589e8
pot = Agama.Potential(file = "L3M11/simulation/agama_potential.ini")

# ╔═╡ 0dcf8e67-570f-4366-b36b-a96a3ecaf128
initial_df = CSV.read(joinpath(ENV["DWARFS_ROOT"], "simulations/all_galaxies/initial_galaxies.csv"), DataFrame)

# ╔═╡ fce7e21b-a2f7-4d8e-9be6-2ef9e3fe1805
icrs_0 = LilGuys.coords_from_df(initial_df)

# ╔═╡ 22ff5b38-edeb-4cd4-a2b7-efffaa638e43
gc_0 = LilGuys.transform.(Galactocentric, icrs_0)

# ╔═╡ 174ea3b9-f95e-40f1-a3ba-f970b6749b68
galaxynames = initial_df.galaxyname

# ╔═╡ 25b2f7b7-0474-437d-83b6-46f11be6b2d7
initial_df.masses[initial_df.galaxyname .== "sculptor"]

# ╔═╡ a83cbc6d-0035-47e0-be5b-df246954ddcf


# ╔═╡ bb151d8e-7329-48e9-be13-0578506d4a8f
md"""
# Comparisons
"""

# ╔═╡ c48273ad-8b2f-482f-851a-b19f5f394ad8
function get_galaxy(traj, galaxyname)
	if galaxyname ∉ galaxynames
		@error "galaxyname $galaxyname not found"
	end
	@assert sum(galaxynames .== galaxyname) == 1
	idx = findfirst(galaxynames .== galaxyname)

	return traj[idx]
end

# ╔═╡ 9a4600b4-d3f0-42c4-94bf-44c296701fd9
function compare_traj_galaxy(galaxyname)
	orbit1 = get_galaxy(traj, galaxyname)
	orbit2 = get_galaxy(traj_alt, galaxyname)
	orbit3 = get_galaxy(traj_no, galaxyname)


	LilGuys.plot_xyz(orbit1.positions, orbit2.positions, orbit3.positions)
end

# ╔═╡ e7b0890c-dcb6-4866-8dd8-49be259e5102
get_galaxy(traj, "sculptor")

# ╔═╡ 0eab06bb-7273-482d-a4ed-f2fa5ad42613
compare_traj_galaxy("sculptor")

# ╔═╡ 92fe77ec-23f9-49de-afe2-e8659c27085c
compare_traj_galaxy("ursa_minor")

# ╔═╡ ff763ea8-d733-46b4-a629-0c4a1f3f856a
function calc_relative_radii(traj, galaxyname)
	idx_no = galaxynames .!= galaxyname
	orbit_ref = get_galaxy(traj, galaxyname)
	
	radii_relative = []
	for orbit in traj
		rs = radii(orbit_ref.positions .- orbit.positions)
		push!(radii_relative, rs)

	end

	return radii_relative
end

# ╔═╡ 87e12931-2fbb-4ea3-9c39-dc3c6d2e1d12
function tidal_force(mass, radius)
	return mass ./ radius .^ 3
end

# ╔═╡ f3d4c4c3-9dba-4682-9a11-fa45749f1398
initial_df.v_circ_max

# ╔═╡ a55616cb-73c7-4594-b54a-cb665288d6c7
LilGuys.M200.([NFW(r_circ_max=initial_df.r_circ_max[i], v_circ_max=initial_df.v_circ_max[i])
			  for i in eachindex(galaxynames) ])

# ╔═╡ 139574b9-f2c2-41e0-85dd-59862a0c4d6b
initial_df.masses

# ╔═╡ d1b3606e-19f6-4e3d-82ac-befb6638b3fb
import LinearAlgebra

# ╔═╡ 71bad51a-17eb-445e-b1fc-736df5efa067
to_sym_matrix(v) = 
	 [
		v[1] v[4] v[6]
		v[4] v[2] v[5]
		v[6] v[5] v[3]
	]

# ╔═╡ dbe25f7c-8100-440b-ae96-ff415209d3ef
function derivative(f, r; dlr=1e-4)
	r1 = r * 10^dlr
	return @. (f(r1) - f(r)) / (r1 - r)
end

# ╔═╡ 4094009b-70bd-40db-bb89-a82c066c2ff6
derivative(x->LilGuys.potential(mw_halo, x), 1000)

# ╔═╡ 173803d4-5232-4fa6-890a-38383b77761b
LilGuys.mass(mw_halo, 1000) / 1000^2

# ╔═╡ d846151e-7c77-4d03-a79f-0092a6dea2f3
function tidal_force(h::NFW, r)
	return abs.(derivative(x->derivative(xx->LilGuys.potential(h, xx), x), r))
end

# ╔═╡ ec1d7f66-dcbd-4af8-b1ce-4cf705ffd84c
function calc_relative_tidal_force(traj, galaxyname)
	radii = calc_relative_radii(traj, galaxyname)

	forces = tidal_force.(initial_df.masses, radii)
end

# ╔═╡ d70eac87-387e-4e66-842f-c5833e522c9b
LilGuys.potential(mw_halo, 1)

# ╔═╡ 9755ecea-ea83-4913-9ac0-efe950210031
tidal_force(mw_halo, radii(traj[1].positions))

# ╔═╡ a4a55afa-f35d-4121-bf2c-5ca84fd0fe3c
LilGuys.acceleration

# ╔═╡ b633ea2b-3364-4de2-8382-c27fbfa59608
function calc_potential_tidal_force_naive(traj, galaxyname)
	orbit = get_galaxy(traj, galaxyname)
	ss = tidal_force(mw_halo, radii(orbit.positions))
	
end

# ╔═╡ 801b26ec-e9fc-4881-8999-06ebadb148cb
function calc_potential_tidal_force(traj, galaxyname)
	orbit = get_galaxy(traj, galaxyname)

	ss = Agama.stress(pot, orbit.positions, units)

	LinearAlgebra.eigmax.(to_sym_matrix.(eachcol(ss))) 
	
end

# ╔═╡ 4dbc3b4e-d75c-481f-9b71-4b7b6221faa9
tidal_scale(units) = Agama.acceleration_scale(units) / Agama.radius_scale(units)

# ╔═╡ feb6e5b3-17ab-4cf5-8dd7-8f8eab9e02fe
function plot_radii_relative(traj, galaxyname; n_max=10)
	rs = calc_relative_radii(traj, galaxyname)
	peris = minimum.(rs)

	idxs = sortperm(peris)[2:n_max+1]
	
	fig = Figure(size=(7*72, 3*72))
	ax = Axis(fig[1,1], 
			  ylabel="galaxy - $galaxyname distance",
			  xlabel = "lookback time",
			  limits=(nothing, nothing, 0, 100)
			 )

	for idx in idxs
		if galaxynames[idx] != galaxyname
			lines!(traj[idx].times * T2GYR, rs[idx], label="$(galaxynames[idx]), $(round(peris[idx], digits=1)) kpc")
		end

	end

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ d0522e05-b664-48ee-aebe-373e3c92519b
function _tidal_force(mass, r; order::Real)
	coeff = Dict(
		1 => 1,
		2 => 1,
		3 => 2,
		4 => 4,
	)[order]

	return coeff .* mass ./ r .^ order
end

# ╔═╡ b8e096eb-eab6-4cdd-8eb2-6476c99367a4
function calc_relative_force(traj, galaxyname; order=2)
	radii = calc_relative_radii(traj, galaxyname)

	halos = [NFW(r_circ_max=initial_df.r_circ_max[i], v_circ_max=initial_df.v_circ_max[i])
			  for i in eachindex(galaxynames) ]

	return [_tidal_force(mass.(halos[i], radii[i]), radii[i], order=order) for i in eachindex(galaxynames)]
end

# ╔═╡ 3b05d152-582b-40ef-b1f8-df198c59e8e8
function calc_potential_force_naive(traj, galaxyname; order=2)
	orbit = get_galaxy(traj, galaxyname)
	r = radii(orbit)

	return _tidal_force(mass.(mw_halo, r), r, order=order)
end

# ╔═╡ 355f5e29-c26a-41f1-9054-1e436e35f638
function plot_force_relative(traj, galaxyname; n_max=10, order=2)
	rs = calc_relative_force(traj, galaxyname; order=order)
	f_max = maximum.(rs)

	idxs = sortperm(f_max, order=Base.Order.ReverseOrdering())[2:n_max+1]
	
	fig = Figure(size=(7*72, 3*72))
	ax = Axis(fig[1,1], 
			  ylabel="force^($order) on $galaxyname from ...",
			  xlabel = "lookback time",
			  limits=(nothing, nothing, nothing, nothing),
			  yscale = log10,
			  yticks = Makie.automatic
			 )

	for idx in idxs
		if galaxynames[idx] != galaxyname
			lines!(traj[idx].times * T2GYR, rs[idx], label=
				   @sprintf "%s; %0.3e" galaxynames[idx] f_max[idx]
				  )
		end

	end


	lines!(traj[1].times * T2GYR, calc_potential_force_naive(traj, galaxyname; order=order), label="background potential", color=:black)
	

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ e7119f98-261f-42c7-9f10-e4bb84353f07
function plot_tides_relative(traj, galaxyname; n_max=10,)
	rs = calc_relative_tidal_force(traj, galaxyname; )
	f_max = maximum.(rs)

	idxs = sortperm(f_max, order=Base.Order.ReverseOrdering())[2:n_max+1]
	
	fig = Figure(size=(7*72, 3*72))
	ax = Axis(fig[1,1], 
			  ylabel="tidal stress on $galaxyname from ___",
			  xlabel = "lookback time",
			  limits=(nothing, nothing, nothing, nothing),
			  yscale = log10,
			  yticks = Makie.automatic
			 )

	for idx in idxs
		if galaxynames[idx] != galaxyname
			lines!(traj[idx].times * T2GYR, rs[idx], label=
				   @sprintf "%s; %0.3e" galaxynames[idx] f_max[idx]
				  )
		end

	end


	lines!(traj[1].times * T2GYR, calc_potential_tidal_force(traj, galaxyname), label="background potential", color=:black)
	

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 3dfaec55-f24e-4a7c-8dab-75ac4e6e602b
function plot_acceleration_relative(traj, galaxyname; n_max=10)
	rs = calc_relative_tidal_force(traj, galaxyname)
	f_max = maximum.(rs)

	idxs = sortperm(f_max, order=Base.Order.ReverseOrdering())[2:n_max+1]
	
	fig = Figure(size=(7*72, 3*72))
	ax = Axis(fig[1,1], 
			  ylabel="tidal force on $galaxyname from ...",
			  xlabel = "lookback time",
			  limits=(nothing, nothing, nothing, nothing),
			  yscale = log10,
			  yticks = Makie.automatic
			 )

	for idx in idxs
		if galaxynames[idx] != galaxyname
			lines!(traj[idx].times * T2GYR, rs[idx], label=
				   @sprintf "%s; %0.3e" galaxynames[idx] f_max[idx]
				  )
		end

	end


	lines!(traj[1].times * T2GYR, calc_potential_tidal_force(traj, galaxyname), label="background potential", color=:black)
	

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ a044ed30-63b3-4491-907d-8f1f29c0f36c
# plot_radii_relative(traj_alt, "sculptor")

# ╔═╡ 6946bf07-d73d-42ef-8664-268bd4ad00c5
plot_radii_relative(traj, "sculptor")

# ╔═╡ c22d1abe-95bf-4723-8ce7-aa0ff438edbc
plot_force_relative(traj, "sculptor")

# ╔═╡ 7e800951-d7df-408d-b0f5-97359e10ab48
plot_force_relative(traj, "sculptor", order=3)

# ╔═╡ cc4c5397-246c-4796-b3bc-7c6b76b844a2
plot_tides_relative(traj, "sculptor")

# ╔═╡ d5803bfe-8f48-4d51-b1ef-34b69112f717
plot_force_relative(traj, "sculptor", order=4)

# ╔═╡ 998f9879-be63-4690-90a4-3f14b1d639e1
md"""
# Ursa Minor
"""

# ╔═╡ 2a32dae8-e9bf-4784-b6bb-291e0c81bfc7
plot_radii_relative(traj, "ursa_minor")

# ╔═╡ 98e0e70b-0721-4c02-826e-d7949c0a21c4
# plot_radii_relative(traj_alt, "ursa_minor")

# ╔═╡ 5893d65c-fae8-483e-9069-e16319987b3f
plot_force_relative(traj, "ursa_minor")

# ╔═╡ 84f587a7-e824-4802-a7b7-af9d1d53555a
plot_force_relative(traj, "ursa_minor", order=3)

# ╔═╡ 9872a469-4e8e-4454-b5cf-58f09f9a4fd6
plot_tides_relative(traj, "ursa_minor")

# ╔═╡ 911dd315-9fb7-44fb-8aa0-0e9215665541
plot_force_relative(traj, "ursa_minor", order=4)

# ╔═╡ 513ae366-6ff1-4dba-84d0-a7b2fc16b7cb
md"""
# Sanity checks
"""

# ╔═╡ 85a364ec-567d-4be9-b270-6b58f2e475e7
init_errors_pos = radii.(LilGuys.position.(gc_0) .- [LilGuys.positions(o)[:, 1] for o in traj])

# ╔═╡ 71daf9d2-6169-40f9-8006-445f2848e812
init_errors_vel = radii.(LilGuys.velocity.(gc_0) ./ V2KMS .+ [LilGuys.velocities(o)[:, 1] for o in traj])

# ╔═╡ 9b594579-0c1e-4864-95b5-ef8cd649abba
@assert all(init_errors_vel .< 2e-5)

# ╔═╡ 67e97789-daf5-448a-9166-39fee71092c6
md"""
# How do derivatives depend on radius ?
"""

# ╔═╡ 0e6ee604-7b52-4cea-9872-c589d0df9dec
halo_scalefree = LilGuys.NFW(M_s=1, r_s=1)

# ╔═╡ f0da4649-c295-4928-ba04-8ad5ff4f3572
halo_agama_scalefree = Agama.Potential(type="NFW")

# ╔═╡ 3e6b041f-2635-4c38-b504-467d6fa4210e
let
	N = 10000
	log_r = LinRange(-2, 3, N)

	r = 10 .^ log_r
	xs = [r zeros(N) zeros(N)]'

	Mtot = LilGuys.mass.(halo_scalefree, r)
	
	
	fig = Figure(size=(3*72, 7*72))
	
	ax = Axis(fig[1, 1], ylabel = "potential")

	lines!(log_r, LilGuys.potential.(halo_scalefree, r))
	lines!(log_r, Agama.potential(halo_agama_scalefree, xs))
	lines!(log_r, -Mtot ./ r)

	ax = Axis(fig[2, 1], ylabel = "acceleration", yscale=log10, yticks=Makie.automatic)

	lines!(log_r, derivative(r->LilGuys.potential.(halo_scalefree, r), r))
	lines!(log_r, -Agama.acceleration(halo_agama_scalefree, xs)[1, :])

	lines!(log_r, Mtot ./ r .^2)
	ylims!(1e-10, 10)



	ax = Axis(fig[3, 1], ylabel = "tidal tensor", yscale=log10, yticks=Makie.automatic)

	lines!(log_r, -derivative(r -> derivative(r->LilGuys.potential.(halo_scalefree, r), r), r))

	s = Agama.stress(halo_agama_scalefree, xs)
	lines!(log_r, LinearAlgebra.eigmax.(to_sym_matrix.(eachcol(s))))
	lines!(log_r, 2 * Mtot ./ r .^3)
	ylims!(1e-10, 1)



	ax = Axis(fig[4, 1], ylabel = "octopole", yscale=log10, yticks=Makie.automatic)

	lines!(log_r, abs.(derivative(r -> derivative(r -> derivative(r->LilGuys.potential.(halo_scalefree, r), r), r), r)))

	lines!(log_r, 4 * Mtot ./ r .^4)
	ylims!(1e-10, 1)


	fig
end

# ╔═╡ Cell order:
# ╠═71389566-8f26-11f0-1805-853166b6d451
# ╠═13eeedb1-814b-4682-8830-f7f86458118f
# ╠═ac9c37c2-a51a-4a10-81b3-2dfe0f7ab11e
# ╠═1f097dca-e93f-4f08-9210-14a1356782f0
# ╠═e8cec7c0-9f08-47ed-a26e-6016fe382ffb
# ╠═903503e8-598a-441e-bf88-7d0a84d221ef
# ╟─541a30dd-992c-4078-afa9-e1cb53dfc4a3
# ╠═51141c7d-59e7-4df0-8928-1e83238fdcff
# ╠═2ca5f1ae-37b5-4a96-b890-eae1ae4609c7
# ╠═0bcbb650-e54e-4a3e-86ff-7ab4b9504b6a
# ╠═e8d209a7-fe07-4ca9-9151-84dd27a75023
# ╠═bafed0b5-26b5-446c-942c-40634bc1102c
# ╠═cc25d8ef-b783-4826-92fa-b1bc0d3589e8
# ╠═0dcf8e67-570f-4366-b36b-a96a3ecaf128
# ╠═fce7e21b-a2f7-4d8e-9be6-2ef9e3fe1805
# ╠═22ff5b38-edeb-4cd4-a2b7-efffaa638e43
# ╠═174ea3b9-f95e-40f1-a3ba-f970b6749b68
# ╠═25b2f7b7-0474-437d-83b6-46f11be6b2d7
# ╠═a83cbc6d-0035-47e0-be5b-df246954ddcf
# ╟─bb151d8e-7329-48e9-be13-0578506d4a8f
# ╠═c48273ad-8b2f-482f-851a-b19f5f394ad8
# ╠═9a4600b4-d3f0-42c4-94bf-44c296701fd9
# ╠═e7b0890c-dcb6-4866-8dd8-49be259e5102
# ╠═0eab06bb-7273-482d-a4ed-f2fa5ad42613
# ╠═92fe77ec-23f9-49de-afe2-e8659c27085c
# ╠═ff763ea8-d733-46b4-a629-0c4a1f3f856a
# ╠═87e12931-2fbb-4ea3-9c39-dc3c6d2e1d12
# ╠═ec1d7f66-dcbd-4af8-b1ce-4cf705ffd84c
# ╠═f3d4c4c3-9dba-4682-9a11-fa45749f1398
# ╠═b8e096eb-eab6-4cdd-8eb2-6476c99367a4
# ╠═a55616cb-73c7-4594-b54a-cb665288d6c7
# ╠═139574b9-f2c2-41e0-85dd-59862a0c4d6b
# ╠═d1b3606e-19f6-4e3d-82ac-befb6638b3fb
# ╠═71bad51a-17eb-445e-b1fc-736df5efa067
# ╠═dbe25f7c-8100-440b-ae96-ff415209d3ef
# ╠═4094009b-70bd-40db-bb89-a82c066c2ff6
# ╠═173803d4-5232-4fa6-890a-38383b77761b
# ╠═d846151e-7c77-4d03-a79f-0092a6dea2f3
# ╠═d70eac87-387e-4e66-842f-c5833e522c9b
# ╠═9755ecea-ea83-4913-9ac0-efe950210031
# ╠═a4a55afa-f35d-4121-bf2c-5ca84fd0fe3c
# ╠═b633ea2b-3364-4de2-8382-c27fbfa59608
# ╠═801b26ec-e9fc-4881-8999-06ebadb148cb
# ╠═4dbc3b4e-d75c-481f-9b71-4b7b6221faa9
# ╠═feb6e5b3-17ab-4cf5-8dd7-8f8eab9e02fe
# ╠═d0522e05-b664-48ee-aebe-373e3c92519b
# ╠═3b05d152-582b-40ef-b1f8-df198c59e8e8
# ╠═355f5e29-c26a-41f1-9054-1e436e35f638
# ╠═e7119f98-261f-42c7-9f10-e4bb84353f07
# ╠═3dfaec55-f24e-4a7c-8dab-75ac4e6e602b
# ╠═a044ed30-63b3-4491-907d-8f1f29c0f36c
# ╠═6946bf07-d73d-42ef-8664-268bd4ad00c5
# ╠═c22d1abe-95bf-4723-8ce7-aa0ff438edbc
# ╠═7e800951-d7df-408d-b0f5-97359e10ab48
# ╠═cc4c5397-246c-4796-b3bc-7c6b76b844a2
# ╠═d5803bfe-8f48-4d51-b1ef-34b69112f717
# ╟─998f9879-be63-4690-90a4-3f14b1d639e1
# ╠═2a32dae8-e9bf-4784-b6bb-291e0c81bfc7
# ╠═98e0e70b-0721-4c02-826e-d7949c0a21c4
# ╠═5893d65c-fae8-483e-9069-e16319987b3f
# ╠═84f587a7-e824-4802-a7b7-af9d1d53555a
# ╠═9872a469-4e8e-4454-b5cf-58f09f9a4fd6
# ╠═911dd315-9fb7-44fb-8aa0-0e9215665541
# ╟─513ae366-6ff1-4dba-84d0-a7b2fc16b7cb
# ╠═85a364ec-567d-4be9-b270-6b58f2e475e7
# ╠═71daf9d2-6169-40f9-8006-445f2848e812
# ╠═9b594579-0c1e-4864-95b5-ef8cd649abba
# ╟─67e97789-daf5-448a-9166-39fee71092c6
# ╠═0e6ee604-7b52-4cea-9872-c589d0df9dec
# ╠═f0da4649-c295-4928-ba04-8ad5ff4f3572
# ╠═3e6b041f-2635-4c38-b504-467d6fa4210e
