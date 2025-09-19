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

# ╔═╡ 9fd93eb0-3e63-4c59-bfd7-c1b9fd307a89
FIGDIR = "figures"

# ╔═╡ ac9c37c2-a51a-4a10-81b3-2dfe0f7ab11e
import Agama

# ╔═╡ 1f097dca-e93f-4f08-9210-14a1356782f0
CairoMakie.activate!(type=:png)

# ╔═╡ e8cec7c0-9f08-47ed-a26e-6016fe382ffb
units = Agama.VASILIEV_UNITS

# ╔═╡ 7d6ee0a5-8d1c-492b-9b7b-33f0bc3e2392
module NBody
	include("nbody_utils.jl")
end

# ╔═╡ 541a30dd-992c-4078-afa9-e1cb53dfc4a3
md"""
# Data loading
"""

# ╔═╡ 51141c7d-59e7-4df0-8928-1e83238fdcff
function load_trajectories(modelname)
	filenames = []
	trajectories = []
	for file in readdir(joinpath(modelname, "out"))
		orbits = last.(LilGuys.read_ordered_structs(joinpath(modelname, "out", file), LilGuys.Orbit))
		push!(trajectories, orbits)

		push!(filenames, file)
	end

	return trajectories
end

# ╔═╡ 2ca5f1ae-37b5-4a96-b890-eae1ae4609c7
traj = load_trajectories("L3M11")

# ╔═╡ 477a6ae1-c2a2-4293-b41e-7122ecad3333
traj_fric = load_trajectories("L3M11_dyn_fric")

# ╔═╡ cc25d8ef-b783-4826-92fa-b1bc0d3589e8
pot = NBody.get_potential("vasiliev24/L3M11/potential")

# ╔═╡ b952b5f7-35f5-462a-862a-0dfe5a33bcc6
pot_mw = NBody.get_potential("vasiliev24/L3M11/potential_mw_halo_evo") + NBody.get_potential("vasiliev24/L3M11/potential_stars")

# ╔═╡ a3e72a55-2456-4634-9008-dfdfc5aaff19
pot_lmc = NBody.get_potential("vasiliev24/L3M11/potential_lmc")

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
function compare_traj_galaxy(traj, galaxyname; color=COLORS[1], alpha=0.1)
	orbits = get_galaxy.(traj, [galaxyname])



	LilGuys.plot_xyz(LilGuys.positions.(orbits)..., color=color, linestyle=:solid, alpha=alpha)
end

# ╔═╡ 0eab06bb-7273-482d-a4ed-f2fa5ad42613
compare_traj_galaxy(traj, "sculptor")

# ╔═╡ 96eba240-aa23-4158-8110-7b222a825c11
compare_traj_galaxy(traj_fric, "sculptor")

# ╔═╡ 92fe77ec-23f9-49de-afe2-e8659c27085c
compare_traj_galaxy(traj, "ursa_minor")

# ╔═╡ 0e662dc3-d830-4cb9-aced-abfbd8e48adf
compare_traj_galaxy(traj_fric, "ursa_minor")

# ╔═╡ 26a167c2-2122-41d8-9a3b-82d18687f2fb
compare_traj_galaxy(traj, "bootes3")

# ╔═╡ 1a7b7f71-1c67-4060-a1f5-747940fb4de0
md"""
# Relative forces
"""

# ╔═╡ 6cff943e-6249-416b-9bef-30eece5f486c
traj

# ╔═╡ ff763ea8-d733-46b4-a629-0c4a1f3f856a
function calc_relative_radii(traj, galaxyname)	
	radii_relative = []
	for orbits in traj
		rs_orbit = []
		orbit_ref = orbits[galaxynames .== [galaxyname]][1]

		for orbit in orbits
			rs = radii(orbit_ref.positions .- orbit.positions)
			push!(rs_orbit, rs)
		end
		push!(radii_relative, rs_orbit)
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

# ╔═╡ 97d5d284-3e42-495a-8737-5829d93523ef
mw_halo = NFW(r_s=11.7, M200=110, delta=100)

# ╔═╡ 4094009b-70bd-40db-bb89-a82c066c2ff6
derivative(x->LilGuys.potential(mw_halo, x), 1000)

# ╔═╡ f420f018-0dda-42bb-8db4-3e53219d9648
LilGuys.M200(mw_halo)

# ╔═╡ 173803d4-5232-4fa6-890a-38383b77761b
LilGuys.mass(mw_halo, 1000) / 1000^2

# ╔═╡ d846151e-7c77-4d03-a79f-0092a6dea2f3
function tidal_force(h::NFW, r)
	return abs.(derivative(x->derivative(xx->LilGuys.potential(h, xx), x), r))
end

# ╔═╡ ec1d7f66-dcbd-4af8-b1ce-4cf705ffd84c
function calc_relative_tidal_force(traj, galaxyname)
	radii = calc_relative_radii(traj, galaxyname)

	forcess = []

	for idx in eachindex(traj)
		forces = tidal_force.(initial_df.masses, radii[idx])
		push!(forcess, forces)
	end

	return forcess
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
function calc_potential_tidal_force(pot, orbits)
	forces = []

	for orbit in orbits
		ss = Agama.stress(pot, orbit.positions, units)
	
		f = LinearAlgebra.eigmax.(to_sym_matrix.(eachcol(ss))) 
		push!(forces, f)
	end

	return forces
		
end

# ╔═╡ 4dbc3b4e-d75c-481f-9b71-4b7b6221faa9
tidal_scale(units) = Agama.acceleration_scale(units) / Agama.radius_scale(units)

# ╔═╡ 7f6d03e8-291d-4968-92e9-6b7a4af4224e
function plot_closest_orbits(traj, func, galaxyname; n_max=10, rev=false, legend=true)
	rs = func(traj, galaxyname)

	if rev
		agg = maximum
	else
		agg = minimum
	end
	
	peris = [10 ^ LilGuys.mean([agg(log10.(r[i])) for r in rs]) for i in eachindex(galaxynames)]
	peri_spread = [LilGuys.std([agg(log10.(r[i])) for r in rs]) for i in eachindex(galaxynames)]

	idxs = sortperm(peris)
	if rev
		idxs = reverse(idxs)
	end

	idxs = idxs[2:n_max+1]
	
	fig = Figure(size=(7*72, 3*72))
	ax = Axis(fig[1,1], 
			  ylabel="galaxy - $galaxyname",
			  xlabel = "lookback time",
			  limits=(nothing, nothing, nothing, nothing),
			  yscale=log10,
			  yticks=Makie.automatic,
			 )

	for (color, idx) in enumerate(reverse(idxs))
		if galaxynames[idx] != galaxyname
			for j in eachindex(rs)
				if j == 1
					label=@sprintf "%s: , %0.2e ± %0.2f dex " galaxynames[idx] peris[idx] peri_spread[idx]
				else
					label = nothing
				end
				
				lines!(traj[j][idx].times * T2GYR, rs[j][idx], label = label => (; alpha=0.5),
					  color = reverse(COLORS)[color], linestyle=:solid, alpha=0.1, )
			end
		end

	end

	if legend
		Legend(fig[1, 2], ax)
	end
	
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

# ╔═╡ bcef1260-2602-4843-901a-2c6125e0a2cb
function plot_mw_tides!(pot, traj, galaxyname; label="MW", color=:black)
	forces = calc_potential_tidal_force(pot, get_galaxy.(traj, [galaxyname]))

	f_max = maximum.([log10.(abs.(f)) for f in forces])

	for i in eachindex(traj)
		if i == 1
			label=@sprintf "%s: %0.2e ± %0.2f dex " label 10^LilGuys.mean(f_max) LilGuys.std(f_max)
		else
			label = nothing
		end
			
		lines!(traj[i][1].times * T2GYR, abs.(forces[i]), color=(color, 0.2), label=label)
	end
end

# ╔═╡ 2e3eacd5-e49a-4cb8-bd90-5a725fa2afec
function compare_tidal_forces(traj, galaxyname)

	fig = plot_closest_orbits(traj, calc_relative_tidal_force, galaxyname, rev=true, legend=false)

	plot_mw_tides!(pot_mw, traj, galaxyname, )
	plot_mw_tides!(pot_lmc, traj, galaxyname, label="LMC", color=COLORS[4])

	Legend(fig[1,2], fig.content[1])
	fig
end


# ╔═╡ 05db5836-70fa-454e-b465-329a26dbdec0
plot_closest_orbits(traj, calc_relative_radii, "sculptor")

# ╔═╡ b9e3b6c0-f373-4aec-8116-006701aa9dd3
compare_tidal_forces(traj, "sculptor")

# ╔═╡ 118314b2-fa8c-42d0-88ab-583769678eb6
compare_tidal_forces(traj_fric, "sculptor")

# ╔═╡ 998f9879-be63-4690-90a4-3f14b1d639e1
md"""
# Ursa Minor
"""

# ╔═╡ 2a32dae8-e9bf-4784-b6bb-291e0c81bfc7
plot_closest_orbits(traj, calc_relative_radii, "ursa_minor")

# ╔═╡ 98e0e70b-0721-4c02-826e-d7949c0a21c4
compare_tidal_forces(traj, "ursa_minor")

# ╔═╡ 6c7bbc07-9e24-48a3-b02e-594b8a4db9c7
compare_tidal_forces(traj_fric, "ursa_minor")

# ╔═╡ Cell order:
# ╠═71389566-8f26-11f0-1805-853166b6d451
# ╠═9fd93eb0-3e63-4c59-bfd7-c1b9fd307a89
# ╠═13eeedb1-814b-4682-8830-f7f86458118f
# ╠═ac9c37c2-a51a-4a10-81b3-2dfe0f7ab11e
# ╠═1f097dca-e93f-4f08-9210-14a1356782f0
# ╠═e8cec7c0-9f08-47ed-a26e-6016fe382ffb
# ╠═7d6ee0a5-8d1c-492b-9b7b-33f0bc3e2392
# ╟─541a30dd-992c-4078-afa9-e1cb53dfc4a3
# ╠═51141c7d-59e7-4df0-8928-1e83238fdcff
# ╠═2ca5f1ae-37b5-4a96-b890-eae1ae4609c7
# ╠═477a6ae1-c2a2-4293-b41e-7122ecad3333
# ╠═cc25d8ef-b783-4826-92fa-b1bc0d3589e8
# ╠═b952b5f7-35f5-462a-862a-0dfe5a33bcc6
# ╠═a3e72a55-2456-4634-9008-dfdfc5aaff19
# ╠═0dcf8e67-570f-4366-b36b-a96a3ecaf128
# ╠═fce7e21b-a2f7-4d8e-9be6-2ef9e3fe1805
# ╠═22ff5b38-edeb-4cd4-a2b7-efffaa638e43
# ╠═174ea3b9-f95e-40f1-a3ba-f970b6749b68
# ╠═25b2f7b7-0474-437d-83b6-46f11be6b2d7
# ╟─bb151d8e-7329-48e9-be13-0578506d4a8f
# ╠═c48273ad-8b2f-482f-851a-b19f5f394ad8
# ╠═9a4600b4-d3f0-42c4-94bf-44c296701fd9
# ╠═0eab06bb-7273-482d-a4ed-f2fa5ad42613
# ╠═96eba240-aa23-4158-8110-7b222a825c11
# ╠═92fe77ec-23f9-49de-afe2-e8659c27085c
# ╠═0e662dc3-d830-4cb9-aced-abfbd8e48adf
# ╠═26a167c2-2122-41d8-9a3b-82d18687f2fb
# ╠═1a7b7f71-1c67-4060-a1f5-747940fb4de0
# ╠═6cff943e-6249-416b-9bef-30eece5f486c
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
# ╠═97d5d284-3e42-495a-8737-5829d93523ef
# ╠═f420f018-0dda-42bb-8db4-3e53219d9648
# ╠═173803d4-5232-4fa6-890a-38383b77761b
# ╠═d846151e-7c77-4d03-a79f-0092a6dea2f3
# ╠═d70eac87-387e-4e66-842f-c5833e522c9b
# ╠═9755ecea-ea83-4913-9ac0-efe950210031
# ╠═a4a55afa-f35d-4121-bf2c-5ca84fd0fe3c
# ╠═b633ea2b-3364-4de2-8382-c27fbfa59608
# ╠═801b26ec-e9fc-4881-8999-06ebadb148cb
# ╠═4dbc3b4e-d75c-481f-9b71-4b7b6221faa9
# ╠═7f6d03e8-291d-4968-92e9-6b7a4af4224e
# ╠═d0522e05-b664-48ee-aebe-373e3c92519b
# ╠═3b05d152-582b-40ef-b1f8-df198c59e8e8
# ╠═bcef1260-2602-4843-901a-2c6125e0a2cb
# ╠═2e3eacd5-e49a-4cb8-bd90-5a725fa2afec
# ╠═05db5836-70fa-454e-b465-329a26dbdec0
# ╠═b9e3b6c0-f373-4aec-8116-006701aa9dd3
# ╠═118314b2-fa8c-42d0-88ab-583769678eb6
# ╟─998f9879-be63-4690-90a4-3f14b1d639e1
# ╠═2a32dae8-e9bf-4784-b6bb-291e0c81bfc7
# ╠═98e0e70b-0721-4c02-826e-d7949c0a21c4
# ╠═6c7bbc07-9e24-48a3-b02e-594b8a4db9c7
