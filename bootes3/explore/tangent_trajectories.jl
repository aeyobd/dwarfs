### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 9550c0c2-6043-11f1-9d42-edfd720dde48
begin
	import Pkg; Pkg.activate()
	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ 8016aba1-27df-448d-a84e-448f70308359
import TOML

# ╔═╡ c34f96ce-d10f-4732-a176-5e0c2120435f
import Agama

# ╔═╡ d8010ec3-667d-4d26-ab8d-567b56c7c212
grillmair09_xi = [24.046071858346227, 20.027734800297978, 16.20445819723798, 12.44329837831643, 8.212939086585296, -0.0006876396768085158, 3.8821844020399965, -2.9607472351154662, -5.4944702309323254]


# ╔═╡ affdb8c0-ec56-4f6f-9125-d5fe5198c956
grillmair09_eta = [7.178929471322114, 6.553063084452138, 5.066975055716812, 3.9928468036547162, 2.5230717387475003, 0.18710127057309833, 1.505234699895075, -1.5705631418921513, -3.5235236001868717]

# ╔═╡ 81d10b2b-ed72-4535-bfd3-38dc64187d2b
orbit_file = joinpath(ENV["DWARFS_ROOT"], "orbits/bootes3/EP2020/orbits.hdf5")

# ╔═╡ 9a48a1a6-4d50-44d9-a89d-1c478fa8ae72
orbit_file_lmc = joinpath(ENV["DWARFS_ROOT"], "orbits/bootes3/vasiliev24_L3M11/orbits.hdf5")

# ╔═╡ e1200eb4-b05f-4cd8-b7e4-cafbf5bb019b
trajectories = last.(LilGuys.read_ordered_structs(orbit_file, LilGuys.Orbit))

# ╔═╡ c7c19b12-c4d9-4f75-8e9d-0b5b1836a2da
trajectories_lmc = last.(LilGuys.read_ordered_structs(orbit_file_lmc, LilGuys.Orbit))

# ╔═╡ 2103624c-e01b-4157-a530-42b551a1826e
pot = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))

# ╔═╡ 243d4a92-a978-4d1a-b21a-3065bb99d22a
pot_lmc = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials/vasiliev24/L3M11/potential.ini"))

# ╔═╡ cf2c9533-6077-41d9-9f34-8333acf32cbe
readdir(dirname(joinpath(ENV["DWARFS_ROOT"], "agama/potentials/vasiliev24/L3M11/potential.ini")))

# ╔═╡ 1a8f81b0-e15b-4c72-a100-a78d0baf343f
units_lmc = Agama.VASILIEV_UNITS

# ╔═╡ f635886b-d690-4777-b650-4790d07c9341
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/bootes3/observed_properties.toml"))

# ╔═╡ 252dc390-54e2-441b-8f3b-b0f08dec50e3
ra0, dec0 = obs_props["ra"], obs_props["dec"]

# ╔═╡ 63b9b50e-0254-4c74-8870-96b9fa5ff872
function to_tangent(orbit; idx_max=600)
	gc = [LilGuys.Galactocentric(orbit.positions[:, ii], V2KMS*orbit.velocities[:, ii]) for ii in 1:length(orbit)]
	icrs = LilGuys.transform.(ICRS, gc)
	ra = [c.ra for c in icrs]
	dec = [c.dec for c in icrs]
	return LilGuys.to_tangent(ra[1:idx_max], dec[1:idx_max], ra0, dec0)
end

# ╔═╡ ac871114-5a14-4143-b1d8-16068b1c1af7
to_tangent(trajectories[1])

# ╔═╡ 8798270a-ef44-4f93-8d1a-b1d68715c7eb
function future_orbit(orbit; pot=pot, units=Agama.AgamaUnits(), time_max=10/T2GYR)
	gc = Galactocentric(orbit.positions[:, 1], orbit.velocities[:, 1]*V2KMS)
	orbit_new = LilGuys.agama_orbit(pot, gc, timerange=(0, time_max), N=10000, agama_units=units)
end

# ╔═╡ b12003e8-14ec-48fc-a408-af0c86222a2f
let
	fig = Figure()
	ax = Axis(fig[1,1], limits=(-30, 30, -30, 30))

	o1 = trajectories[5]
	o2 = future_orbit(o1, time_max=3/T2GYR)
	o3 = future_orbit(o1, time_max=-3/T2GYR)

	for o in [o1, o2, o3]
		xi, eta = to_tangent(o)
		lines!(xi, eta, alpha=0.3)
	end

	
	fig
end

# ╔═╡ fbfc6f9b-e172-425d-896a-54a3d8c40a87
CairoMakie.activate!(type=:png)

# ╔═╡ 1b1418ca-0589-4f84-9063-6952947e9e01
let
	fig = Figure(size=(4, 2.5) .* 72)
	ax = Axis(
		fig[1,1],
		limits=(-15, 15, -15, 15),
		aspect=1, 
		xlabel=L"$\xi$ / deg",
		ylabel=L"$\eta$ / deg", 
		xreversed=true, 
		title="MW-only",
	)

	for i in 1:100
		orbit = trajectories[i]
		xi, eta = to_tangent(orbit)
		lines!(xi, eta, color=COLORS[1], alpha=0.1)
		
		orbit_future = future_orbit(trajectories[i])
		xi, eta = to_tangent(orbit_future)
		lines!(xi, eta, color=COLORS[1], alpha=0.1)

	end

	scatter!(grillmair09_xi, grillmair09_eta, color=:black)



	ax_lmc = Axis(
		fig[1,2],
		limits=(-15, 15, -15, 15),
		aspect=1, 
		xlabel=L"$\xi$ / deg",
		xreversed=true, 
		title="MW+LMC",
	)

	for i in 1:100
		orbit = trajectories_lmc[i]
		xi, eta = to_tangent(orbit)
		lines!(xi, eta, color=COLORS[2], alpha=0.1)
		
		orbit_future = future_orbit(trajectories_lmc[i], units=units_lmc, pot=pot_lmc)
		xi, eta = to_tangent(orbit_future)
		lines!(xi, eta, color=COLORS[2], alpha=0.1)

	end
	hideydecorations!(ax_lmc, ticks=false, minorticks=false)


	
	scatter!(grillmair09_xi, grillmair09_eta, color=:black)
	linkaxes!(ax, ax_lmc)


	@savefig "tangent_trajectories_w_wo_lmc"
	fig
end
	

# ╔═╡ Cell order:
# ╠═9550c0c2-6043-11f1-9d42-edfd720dde48
# ╠═8016aba1-27df-448d-a84e-448f70308359
# ╠═c34f96ce-d10f-4732-a176-5e0c2120435f
# ╠═d8010ec3-667d-4d26-ab8d-567b56c7c212
# ╠═affdb8c0-ec56-4f6f-9125-d5fe5198c956
# ╠═81d10b2b-ed72-4535-bfd3-38dc64187d2b
# ╠═9a48a1a6-4d50-44d9-a89d-1c478fa8ae72
# ╠═e1200eb4-b05f-4cd8-b7e4-cafbf5bb019b
# ╠═c7c19b12-c4d9-4f75-8e9d-0b5b1836a2da
# ╠═2103624c-e01b-4157-a530-42b551a1826e
# ╠═243d4a92-a978-4d1a-b21a-3065bb99d22a
# ╠═cf2c9533-6077-41d9-9f34-8333acf32cbe
# ╠═1a8f81b0-e15b-4c72-a100-a78d0baf343f
# ╠═f635886b-d690-4777-b650-4790d07c9341
# ╠═252dc390-54e2-441b-8f3b-b0f08dec50e3
# ╠═63b9b50e-0254-4c74-8870-96b9fa5ff872
# ╠═ac871114-5a14-4143-b1d8-16068b1c1af7
# ╠═8798270a-ef44-4f93-8d1a-b1d68715c7eb
# ╠═b12003e8-14ec-48fc-a408-af0c86222a2f
# ╠═fbfc6f9b-e172-425d-896a-54a3d8c40a87
# ╠═1b1418ca-0589-4f84-9063-6952947e9e01
