### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 41f05fc0-8ce3-11f0-3c2a-d3ed3074bf54
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya
	using CairoMakie
end

# ╔═╡ e8b9ee9f-f1ba-42f1-8e34-1e850ccbfd77
using HDF5

# ╔═╡ c7856edd-9b38-47df-b2eb-e783d0aad32d
using CSV, DataFrames

# ╔═╡ 1a618a0b-4942-4d0d-998d-940f8adb5697
CairoMakie.activate!(type=:png)

# ╔═╡ 6e65cc1b-0d65-417f-b21a-6a60337b5f0f


# ╔═╡ 6bdb5577-b4a2-4402-9951-ac4aba8a4a0a


# ╔═╡ 33925b89-68fb-48eb-8172-c82f4e075fa7


# ╔═╡ bba14e40-dfb9-48ee-af59-d3e3cd65718d
function load_trajectories(modelname)
	local pos, vel, t
	h5open(joinpath(modelname, "trajectory.hdf5")) do h5 
	
		pos = h5["positions"][:, :, :]
		vel = h5["velocities"][:, :, :]
		t = h5["times"][:]
	end

	return [Orbit(positions=pos[:, i, :], velocities=vel[:, i, :], times=t) for i in 1:size(pos, 2)]
end

# ╔═╡ 729544e1-e52a-4657-9e53-a849cabd4900
traj = load_trajectories("L3M11")

# ╔═╡ 8cf441e3-e19f-49f0-8cae-c000f51f68a0
traj_alt = load_trajectories("L3M11_nointeract")

# ╔═╡ 60d970f9-27ad-4808-8675-716554ff8c7b
dr = [radii(a.positions .- b.positions) .^ 2 .+ radii(a.velocities .- b.velocities) .^ 2 for (a, b) in zip(traj, traj_alt)]

# ╔═╡ cb886ca9-84c8-40f8-a3b7-44ebcbbf294f
[rs[1] for rs in dr]

# ╔═╡ 5965d2f2-43d5-4095-9051-e532c388f5af
initial_df = CSV.read(joinpath(ENV["DWARFS_ROOT"], "simulations/all_galaxies/initial_galaxies.csv"), DataFrame)

# ╔═╡ e2d0bbc5-82ac-4ffd-8180-dd36a4a0e4a5
let
	fig = Figure()
	ax = Axis(fig[1,1], limits=(0, 1909, -10, 5), 
		xlabel = "log time", 
		ylabel = "log deviation"
			 )

	for i in eachindex(traj)
		lines!((traj[i].times .- traj[i].times[1]), log10.(dr[i]))
	end

	fig
end

# ╔═╡ b777bdfb-0587-4843-9555-293e1e8228ad
let
	fig = Figure()
	ax = Axis(fig[1,1], limits=(0, 190, -10,2), 
		xlabel = " time", 
		ylabel = "log deviation"
			 )

	for i in eachindex(traj)
		lines!((traj[i].times .- traj[i].times[1]), log10.(dr[i]))
	end

	fig
end

# ╔═╡ 30de90ab-6f86-4e7e-bc3a-8bb9000ea6f8
reverse(initial_df.galaxyname[sortperm(maximum.(dr))])

# ╔═╡ 6d216084-b58c-402b-9443-91a16a6fdf19
reverse(sortperm(maximum.(dr)))

# ╔═╡ f025ce8a-50e9-487a-bd8d-9abfc6f48a18
idx_ex = 24

# ╔═╡ 04c98f29-c606-494c-8449-3fc9edf741e0
initial_df.masses

# ╔═╡ da5f0835-2a44-4009-aaff-ad61f44c8bb9
LilGuys.plot_xyz(traj[idx_ex].positions, traj_alt[idx_ex].positions)

# ╔═╡ f228a24a-d3b9-4276-9ba3-92c36321606e
orbit_1 = traj[idx_ex]

# ╔═╡ 81395d03-35f0-4da8-829e-1cfdbd7fffbf
orbit_2 = traj_alt[idx_ex]

# ╔═╡ 0700aee3-7880-4caf-a05a-b0e4627a893b


# ╔═╡ 970fbf03-7597-48ba-a0aa-2c0527ce93fe
let
	fig = Figure()
	ax = Axis(fig[1,1], limits=(0, 1909, -15, 5))

	lines!((traj[idx_ex].times), log10.(dr[idx_ex]))
	lines!((traj[idx_ex].times), log10.(radii(orbit_1.positions .- orbit_2.positions)))
	lines!((traj[idx_ex].times), log10.(radii(orbit_1.velocities .- orbit_2.velocities)))
	fig
end

# ╔═╡ Cell order:
# ╠═41f05fc0-8ce3-11f0-3c2a-d3ed3074bf54
# ╠═e8b9ee9f-f1ba-42f1-8e34-1e850ccbfd77
# ╠═c7856edd-9b38-47df-b2eb-e783d0aad32d
# ╠═1a618a0b-4942-4d0d-998d-940f8adb5697
# ╠═6e65cc1b-0d65-417f-b21a-6a60337b5f0f
# ╠═6bdb5577-b4a2-4402-9951-ac4aba8a4a0a
# ╠═33925b89-68fb-48eb-8172-c82f4e075fa7
# ╠═bba14e40-dfb9-48ee-af59-d3e3cd65718d
# ╠═729544e1-e52a-4657-9e53-a849cabd4900
# ╠═8cf441e3-e19f-49f0-8cae-c000f51f68a0
# ╠═60d970f9-27ad-4808-8675-716554ff8c7b
# ╠═cb886ca9-84c8-40f8-a3b7-44ebcbbf294f
# ╠═5965d2f2-43d5-4095-9051-e532c388f5af
# ╠═e2d0bbc5-82ac-4ffd-8180-dd36a4a0e4a5
# ╠═b777bdfb-0587-4843-9555-293e1e8228ad
# ╠═30de90ab-6f86-4e7e-bc3a-8bb9000ea6f8
# ╠═6d216084-b58c-402b-9443-91a16a6fdf19
# ╠═f025ce8a-50e9-487a-bd8d-9abfc6f48a18
# ╠═04c98f29-c606-494c-8449-3fc9edf741e0
# ╠═da5f0835-2a44-4009-aaff-ad61f44c8bb9
# ╠═f228a24a-d3b9-4276-9ba3-92c36321606e
# ╠═81395d03-35f0-4da8-829e-1cfdbd7fffbf
# ╠═0700aee3-7880-4caf-a05a-b0e4627a893b
# ╠═970fbf03-7597-48ba-a0aa-2c0527ce93fe
