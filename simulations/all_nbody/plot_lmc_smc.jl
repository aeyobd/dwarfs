### A Pluto.jl notebook ###
# v0.20.21

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
function load_trajectories(modelname, outdir="out")
	filenames = []
	trajectories = []
	for file in readdir(joinpath(modelname, outdir))
		if splitext(file)[end] == ".hdf5"
			orbits = last.(LilGuys.read_ordered_structs(joinpath(modelname, outdir, file), LilGuys.Orbit))
			push!(trajectories, orbits)
	
			push!(filenames, file)
		end
	end

	return trajectories
end

# ╔═╡ 2ca5f1ae-37b5-4a96-b890-eae1ae4609c7
traj_nosmc = load_trajectories("MW_LMC_SMC", "out_nosmc")

# ╔═╡ eb91aea8-79ae-4b75-b969-9abe72f7b7d3
traj = load_trajectories("MW_LMC_SMC")

# ╔═╡ 01dcda87-b5cf-45bc-abf5-fea731dd0d7d
traj_nofric = load_trajectories("MW_LMC_SMC", "out_nofric")

# ╔═╡ 0dcf8e67-570f-4366-b36b-a96a3ecaf128
initial_df = CSV.read(joinpath(ENV["DWARFS_ROOT"], "simulations/all_galaxies/initial_galaxies.csv"), DataFrame)

# ╔═╡ fce7e21b-a2f7-4d8e-9be6-2ef9e3fe1805
icrs_0 = LilGuys.coords_from_df(initial_df)

# ╔═╡ 22ff5b38-edeb-4cd4-a2b7-efffaa638e43
gc_0 = LilGuys.transform.(Galactocentric, icrs_0)

# ╔═╡ 25b2f7b7-0474-437d-83b6-46f11be6b2d7
initial_df.masses[initial_df.galaxyname .== "sculptor"]

# ╔═╡ bb151d8e-7329-48e9-be13-0578506d4a8f
md"""
# Comparisons
"""

# ╔═╡ 936c090a-dd4d-44e6-be5b-441f4039be3a
galaxynames = ["mw", "lmc", "smc", "sculptor", "ursa_minor"]

# ╔═╡ c48273ad-8b2f-482f-851a-b19f5f394ad8
function get_galaxy(traj, galaxyname, reference=nothing)
	if galaxyname ∉ galaxynames
		@error "galaxyname $galaxyname not found"
	end
	@assert sum(galaxynames .== [galaxyname]) == 1
	idx = findfirst(galaxynames .== [galaxyname])

	if !isnothing(reference)
		idx_ref = findfirst(galaxynames .== [reference])
		return 	traj[idx] - traj[idx_ref]
	end
		
	return traj[idx]
end

# ╔═╡ 9a4600b4-d3f0-42c4-94bf-44c296701fd9
function compare_traj_galaxy!(ax, traj, galaxyname, reference="mw"; color=COLORS[1], alpha=0.1, time_min=-10/T2GYR)
	orbits = get_galaxy.(traj, galaxyname, reference) 

	positions = []

	for orbit in orbits
		push!(positions, orbit.positions[:, orbit.times .> time_min])
	end

	LilGuys.plot_xyz!(ax, positions..., color=color, linestyle=:solid, alpha=alpha)
end

# ╔═╡ a8100604-2123-4ef0-ba88-3847ad485f57
function compare_traj_galaxy(traj, galaxyname, reference="mw"; color=COLORS[1], alpha=0.1, time_min=-10/T2GYR)
	orbits = get_galaxy.(traj, galaxyname, reference) 

	positions = []

	for orbit in orbits
		push!(positions, orbit.positions[:, orbit.times .> time_min])
	end

	LilGuys.plot_xyz(positions..., color=color, linestyle=:solid, alpha=alpha)
end

# ╔═╡ e7ba552c-fef6-4b4d-a270-6eceefdecd61
function compare_traj_galaxy_radii!(traj, traj2::Vector{<:AbstractVector}, galaxyname, reference; color=COLORS[1], alpha=0.1, kwargs...)
	orbits = get_galaxy.(traj, galaxyname, reference) 
	orbits2 = get_galaxy.(traj2, galaxyname, reference) 

	for i in eachindex(orbits)

		lines!(orbits[i].times, log10.(radii(orbits[i].positions, orbits2[i].positions));
			   color=color, linestyle=:solid, alpha=alpha, kwargs...)
	end

end

# ╔═╡ 8f13863e-b79c-4c45-840d-14636ec272be
function compare_traj_galaxy_radii!(traj, galaxyname::String, reference = "mw"; color=COLORS[1], alpha=0.1, kwargs...)
	orbits = get_galaxy.(traj, galaxyname, reference) 

	for orbit in orbits

		lines!(orbit.times, log10.(radii(orbit));
			   color=color, linestyle=:solid, alpha=alpha, kwargs...)
	end

end

# ╔═╡ 502ada8b-2d0d-4dbd-92ee-4891ad68cdd8
function compare_traj_galaxy_radii(args...;kwargs...)
	fig = Figure()
	ax = Axis(fig[1,1])
	compare_traj_galaxy_radii!(args...; kwargs...)

	fig
end

# ╔═╡ d597e65d-dfb2-41f5-8821-ad1f8f9614b9
compare_traj_galaxy(traj, "mw", nothing)

# ╔═╡ 0eab06bb-7273-482d-a4ed-f2fa5ad42613
compare_traj_galaxy(traj, "lmc")

# ╔═╡ 1507bbf5-e6a4-49a9-bb6a-2f079bb9669c
let
	f = compare_traj_galaxy(traj, "lmc", "sculptor", time_min=-0.3/T2GYR)
	compare_traj_galaxy!(f.content, traj, "smc", "sculptor", time_min=-0.3/T2GYR, color=COLORS[2])

	LilGuys.plot_xyz!(f.content, zeros(3, 1), plot=:scatter, color=:black)
	f
end

# ╔═╡ fa3492db-087c-413d-b438-dc09e5c1d80c
compare_traj_galaxy_radii(traj, "lmc")

# ╔═╡ 4fc6f578-b37b-4e1a-b8c8-fe0e9fdbbefd
compare_traj_galaxy_radii([(o1[[1,2,4,5]] .- o2) for (o1, o2) in zip(traj, traj_nosmc)], "lmc")

# ╔═╡ 2300e33c-1b4d-4751-ab3f-00afe02b5e28
compare_traj_galaxy_radii([(o1 .- o2) for (o1, o2) in zip(traj, traj_nofric)], "lmc")

# ╔═╡ 12d78b71-123a-4b32-88b0-6b41f564e91f
compare_traj_galaxy_radii(traj_nofric, "lmc")

# ╔═╡ b928f087-17ee-48fb-8446-c67d69e8288f
compare_traj_galaxy_radii(traj, "smc", "lmc")

# ╔═╡ 6d053852-c5f3-49fa-a110-d233ac6b06c4
let
	fig = compare_traj_galaxy_radii(traj, "sculptor", label="Scl_MW")
	compare_traj_galaxy_radii!(traj, "sculptor", "lmc", color=COLORS[2], label="Scl-LMC")
	compare_traj_galaxy_radii!(traj, "sculptor", "smc", color=COLORS[3], label="Scl-SMC")

	axislegend(merge=true)
	fig

end

# ╔═╡ 4a03701e-84fa-4818-a7e7-094f275ecbae
let
	fig = compare_traj_galaxy_radii(traj_nosmc, "smc", label="Scl_MW")
	compare_traj_galaxy_radii!(traj_nosmc, "smc", "lmc", color=COLORS[2], label="Scl-LMC")

	axislegend(merge=true)
	fig

end

# ╔═╡ 8c520387-1d0c-4fbd-8696-fb9acc1c3cfb
compare_traj_galaxy(traj_nosmc, "lmc")

# ╔═╡ e3ecb803-f33b-4f5a-8288-aa95411e3fd2
compare_traj_galaxy(traj, "smc")

# ╔═╡ f868386e-976f-4254-be16-0452bb2b0bd7
compare_traj_galaxy(traj, "smc", "lmc")

# ╔═╡ 598b692d-c083-4d66-bab2-a93af461957c
compare_traj_galaxy_radii(traj, "smc")

# ╔═╡ 92fe77ec-23f9-49de-afe2-e8659c27085c
compare_traj_galaxy(traj, "ursa_minor")

# ╔═╡ 02de13b9-fc64-4127-9beb-1321f6d8e22c
let
	fig = compare_traj_galaxy_radii(traj, "ursa_minor", label="UMi-MW")
	compare_traj_galaxy_radii!(traj, "ursa_minor", "lmc", color=COLORS[2], label="UMi-LMC")
	compare_traj_galaxy_radii!(traj, "ursa_minor", "smc", color=COLORS[3], label="UMi-SMC")

	axislegend(merge=true)
	fig

end

# ╔═╡ 95ebb8d1-5618-48df-882a-de47f56e150e
let
	fig = compare_traj_galaxy_radii(traj_nofric, "ursa_minor", label="UMi-MW")
	compare_traj_galaxy_radii!(traj_nofric, "ursa_minor", "lmc", color=COLORS[2], label="UMi-LMC")
	compare_traj_galaxy_radii!(traj_nofric, "ursa_minor", "smc", color=COLORS[3], label="UMi-SMC")

	axislegend(merge=true)
	fig

end

# ╔═╡ 3b9fd610-6a15-492c-b18f-e1e255a034c2
let
	fig = compare_traj_galaxy_radii(traj_nosmc, "sculptor", label="UMi-MW")
	compare_traj_galaxy_radii!(traj_nosmc, "sculptor", "lmc", color=COLORS[2], label="UMi-LMC")

	axislegend(merge=true)
	fig

end

# ╔═╡ 3af70d8e-e7bc-45d8-80cb-eb4a4c717221
hist(log10.(last.(radii.(get_galaxy.(traj, "ursa_minor", "lmc")))))

# ╔═╡ 95b9f993-a094-43f5-b99a-4bba5a2bd8c8
hist(log10.(last.(radii.(get_galaxy.(traj_nofric, "ursa_minor", "lmc")))))

# ╔═╡ 5527e3dc-952c-4b37-90f3-83f65d9d31aa
hist(log10.(last.(radii.(get_galaxy.(traj_nosmc, "sculptor", "lmc")))))

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
# ╠═eb91aea8-79ae-4b75-b969-9abe72f7b7d3
# ╠═01dcda87-b5cf-45bc-abf5-fea731dd0d7d
# ╠═0dcf8e67-570f-4366-b36b-a96a3ecaf128
# ╠═fce7e21b-a2f7-4d8e-9be6-2ef9e3fe1805
# ╠═22ff5b38-edeb-4cd4-a2b7-efffaa638e43
# ╠═25b2f7b7-0474-437d-83b6-46f11be6b2d7
# ╟─bb151d8e-7329-48e9-be13-0578506d4a8f
# ╠═936c090a-dd4d-44e6-be5b-441f4039be3a
# ╠═c48273ad-8b2f-482f-851a-b19f5f394ad8
# ╠═9a4600b4-d3f0-42c4-94bf-44c296701fd9
# ╠═a8100604-2123-4ef0-ba88-3847ad485f57
# ╠═502ada8b-2d0d-4dbd-92ee-4891ad68cdd8
# ╠═e7ba552c-fef6-4b4d-a270-6eceefdecd61
# ╠═8f13863e-b79c-4c45-840d-14636ec272be
# ╠═d597e65d-dfb2-41f5-8821-ad1f8f9614b9
# ╠═0eab06bb-7273-482d-a4ed-f2fa5ad42613
# ╠═1507bbf5-e6a4-49a9-bb6a-2f079bb9669c
# ╠═fa3492db-087c-413d-b438-dc09e5c1d80c
# ╠═4fc6f578-b37b-4e1a-b8c8-fe0e9fdbbefd
# ╠═2300e33c-1b4d-4751-ab3f-00afe02b5e28
# ╠═12d78b71-123a-4b32-88b0-6b41f564e91f
# ╠═b928f087-17ee-48fb-8446-c67d69e8288f
# ╠═6d053852-c5f3-49fa-a110-d233ac6b06c4
# ╠═4a03701e-84fa-4818-a7e7-094f275ecbae
# ╠═8c520387-1d0c-4fbd-8696-fb9acc1c3cfb
# ╠═e3ecb803-f33b-4f5a-8288-aa95411e3fd2
# ╠═f868386e-976f-4254-be16-0452bb2b0bd7
# ╠═598b692d-c083-4d66-bab2-a93af461957c
# ╠═92fe77ec-23f9-49de-afe2-e8659c27085c
# ╠═02de13b9-fc64-4127-9beb-1321f6d8e22c
# ╠═95ebb8d1-5618-48df-882a-de47f56e150e
# ╠═3b9fd610-6a15-492c-b18f-e1e255a034c2
# ╠═3af70d8e-e7bc-45d8-80cb-eb4a4c717221
# ╠═95b9f993-a094-43f5-b99a-4bba5a2bd8c8
# ╠═5527e3dc-952c-4b37-90f3-83f65d9d31aa
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
