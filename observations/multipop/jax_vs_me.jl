### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 8a84cd7a-b75e-11f0-ac6e-5386c69b3015
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ 8d13fc1a-c60f-447f-b938-daf0ed3f8be0
using CSV, DataFrames

# ╔═╡ df00a157-8166-4bfe-9fb2-1e41baf70dd0
galaxyname = "ursa_minor"

# ╔═╡ a3767358-e88b-45e2-ac95-d0a635833860
import TOML

# ╔═╡ 30bc1fdb-d233-431f-b784-16baa272c0cc
α = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 678a0ff7-a901-498d-b739-e5b6b8d26a0e
function get_obs_props(galaxyname)
	TOML.parsefile(joinpath("..", galaxyname, "observed_properties.toml"))
end

# ╔═╡ 4230edae-3a81-4577-a48d-4c9a9f124744
function get_prof(galaxyname)
	prof = LilGuys.SurfaceDensityProfile(joinpath("..", galaxyname, "density_profiles/fiducial_profile.toml"))

	LilGuys.scale(prof, 1, 1/sum(prof.counts))

end

# ╔═╡ eada9dc1-640c-495d-bc6a-4553db4542a1
obs_prof = get_prof(galaxyname)

# ╔═╡ a17b43d0-4e62-455c-a2a9-d5243fd69ef2
function jax_M(R_s, Σ_0=1)
	return 2π * R_s^2 * Σ_0
end

# ╔═╡ 28cdcfb2-fe49-4fa7-9452-dc6448562569
function load_jax_model(galaxyname)
	fit = TOML.parsefile(joinpath("..", galaxyname, "j+24_fsat.toml"))["2c"]

	obs_props = get_obs_props(galaxyname)
	R_1 = obs_props["r_h"] / α

	M_1 = jax_M(R_1)

	R_2 = fit["r_s"] * 60
	M_2 = jax_M(R_2, fit["B"])
	prof = LilGuys.DoubleExp2D(M_1, R_1, M_2, R_2)

	LilGuys.scale(prof, sqrt(1 - obs_props["ellipticity"]), 1 / mass(prof))
end

# ╔═╡ 2d895f35-0bd5-4aba-94bb-0cfd61fdd98a
jax_model = load_jax_model(galaxyname)

# ╔═╡ d94c4a57-865c-4022-8b5d-5e4d5480dfce
function get_key(df, key)
	idx = only(findall(df.parameters .== [key]))
	return df.median[idx]
end

# ╔═╡ 7ec7132a-a405-4ab7-823c-e7a9a7328503
function load_mcmc_model(galaxyname)
	df = CSV.read(joinpath("..", galaxyname, "mcmc/summary.mcmc_2exp.csv"), DataFrame)

	R_1 = get_key(df, "R_s")
	f_b = get_key(df, "f_outer")
	M_1 = 1-f_b
	M_2 = f_b
	R_2 = get_key(df, "R_s_outer")

	prof = LilGuys.DoubleExp2D(M_1, R_1, M_2, R_2)

end

# ╔═╡ 5b37aea0-4082-4555-8849-43d67960c099
my_model = load_mcmc_model(galaxyname)

# ╔═╡ 7fc0c3a0-31a9-4dbc-abc2-7b917391e365
10 .^ obs_prof.log_R, 10 .^ obs_prof.log_Sigma

# ╔═╡ 35269808-441f-4c30-8861-d1b8efb5c94c
0.316 * 60


# ╔═╡ a2be855f-aa12-47ef-9d28-0191ada0866f
let
	fig = Figure()
	ax = Axis(fig[1,1], xscale=log10, yscale=log10, xticks=Makie.automatic, yticks=Makie.automatic)

	x = logrange(1, 100, 100)

	lines!(x, LilGuys.surface_density.(jax_model, x))
	lines!(x, LilGuys.surface_density.(my_model, x))

	scatter!(10 .^ obs_prof.log_R, 10 .^ obs_prof.log_Sigma)
	fig

end

# ╔═╡ Cell order:
# ╠═df00a157-8166-4bfe-9fb2-1e41baf70dd0
# ╠═8a84cd7a-b75e-11f0-ac6e-5386c69b3015
# ╠═a3767358-e88b-45e2-ac95-d0a635833860
# ╠═8d13fc1a-c60f-447f-b938-daf0ed3f8be0
# ╠═30bc1fdb-d233-431f-b784-16baa272c0cc
# ╠═678a0ff7-a901-498d-b739-e5b6b8d26a0e
# ╠═4230edae-3a81-4577-a48d-4c9a9f124744
# ╠═eada9dc1-640c-495d-bc6a-4553db4542a1
# ╠═a17b43d0-4e62-455c-a2a9-d5243fd69ef2
# ╠═28cdcfb2-fe49-4fa7-9452-dc6448562569
# ╠═2d895f35-0bd5-4aba-94bb-0cfd61fdd98a
# ╠═7ec7132a-a405-4ab7-823c-e7a9a7328503
# ╠═5b37aea0-4082-4555-8849-43d67960c099
# ╠═d94c4a57-865c-4022-8b5d-5e4d5480dfce
# ╠═7fc0c3a0-31a9-4dbc-abc2-7b917391e365
# ╠═35269808-441f-4c30-8861-d1b8efb5c94c
# ╠═a2be855f-aa12-47ef-9d28-0191ada0866f
