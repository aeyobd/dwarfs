### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./style.jl")

# ╔═╡ caf046d9-e475-4810-844e-b65e9c08c832
stars_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", "sculptor/1e7_V31_r3.2/", "stars")

# ╔═╡ 7a6d6345-fefb-4017-8230-3ff5bea55c83
snap = Snapshot(joinpath(stars_dir, "iso_paint.hdf5"))

# ╔═╡ a3a3a875-d0cc-4594-86c2-b242e859efef
df_probs = LilGuys.read_hdf5_table(joinpath(stars_dir, "exp2d_rs0.10", "probabilities_stars.hdf5"))

# ╔═╡ 01e34341-e806-44ef-9591-d226b345dec6
snap.weights = df_probs.probability[snap.index]

# ╔═╡ 95f811eb-b7c0-48e9-8b88-716cd667db3d
prof_dm = DensityProfile(snap)

# ╔═╡ 3edb5998-a200-4e94-9598-5b73c29def33
sum(densities(prof_dm) .* diff(π*LilGuys.radius_bins(prof_dm) .^3)) 

# ╔═╡ 5743880e-567f-4af5-acb5-e2ec420110aa
Mstar = 5e5 / M2MSUN

# ╔═╡ 098e33d5-22a5-496b-943b-b8c194cdcc62
prof_stars = DensityProfile(snap, snap.weights * Mstar)

# ╔═╡ f6d51dd1-91df-4a6b-85a8-cd9b26d2c9c9
sum(densities(prof_stars) .* diff(π*LilGuys.radius_bins(prof_stars) .^3)) 

# ╔═╡ 546e88aa-c3c5-41b6-84a7-ce733e8f309f


# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "initial_conditions" let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "log radius / kpc",
		ylabel = L"log $\rho$",
		limits = (-2, 2, -9.5, 0.5)
	)

	lines!(log_radii(prof_dm), log_densities(prof_dm), label="dark matter")
	lines!(log_radii(prof_stars), log_densities(prof_stars), label="stars")

	axislegend(margin=(19,19,19,19), position=:lb)
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═caf046d9-e475-4810-844e-b65e9c08c832
# ╠═7a6d6345-fefb-4017-8230-3ff5bea55c83
# ╠═a3a3a875-d0cc-4594-86c2-b242e859efef
# ╠═01e34341-e806-44ef-9591-d226b345dec6
# ╠═95f811eb-b7c0-48e9-8b88-716cd667db3d
# ╠═3edb5998-a200-4e94-9598-5b73c29def33
# ╠═f6d51dd1-91df-4a6b-85a8-cd9b26d2c9c9
# ╠═5743880e-567f-4af5-acb5-e2ec420110aa
# ╠═098e33d5-22a5-496b-943b-b8c194cdcc62
# ╠═546e88aa-c3c5-41b6-84a7-ce733e8f309f
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
