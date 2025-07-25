### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 5fa9aec4-f876-11ef-1eed-bd8da0982d2b
begin
	using Pkg; Pkg.activate()

	using CairoMakie
	using Arya

	using PlutoUI
end

# ╔═╡ 01240d81-5b6e-41b6-a2be-32b87a6568ef
include("./paper_style.jl")

# ╔═╡ 4b9e4fb8-f257-47c5-b1fa-0b4e9bdd042c
modelname = "1e6_v31_r3.2/orbit_15_100"

# ╔═╡ 6feb26d5-8b1a-4913-95a2-d57c237fdadb
idx_f = 98

# ╔═╡ 258525e6-7fa1-46b3-95f6-950c2b82aafc
idx_i = 1

# ╔═╡ ea3d8d00-11ee-4b60-9b6a-58a3e396a572
md"""
# Setup
"""

# ╔═╡ 221247cd-e447-400c-8491-b7248787f7a7
import TOML

# ╔═╡ d143090f-a2e1-4068-a3e7-6ef35635f758
galaxies = ["sculptor", "ursa_minor"]

# ╔═╡ 047193b4-dcf7-4c36-aafb-2179d8363f3c
function notebook_inputs(; kwargs...)
	return PlutoUI.combine() do Child
		
		user_inputs = [
			md""" $(string(name)): $(
				Child(name, obj)
			)"""
			
			for (name, obj) in kwargs
		]
		
		md"""
		#### Inputs
		$(user_inputs)
		"""
	end
end

# ╔═╡ 70ac654f-b45f-4b14-ad0f-01cf5f1e29e5
model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", "idealized", modelname)

# ╔═╡ 76fdb3ab-60a6-42ea-bb71-0204d37284ce
#=╠═╡
stars_dir_in = joinpath(model_dir, "../stars/$starsname")
  ╠═╡ =#

# ╔═╡ b1edbdfc-298b-40e2-b6c0-687a866fe180
#=╠═╡
stars_dir_out = joinpath(model_dir, "stars/$starsname")
  ╠═╡ =#

# ╔═╡ b51a4ddd-90fd-415f-8a48-e4f35aeee801
FIGDIR = "figures"

# ╔═╡ 81684bfc-3316-42c2-a2a9-d9db98858edd
using LilGuys; FIGDIR

# ╔═╡ 94f3bf4a-248a-466e-b2de-b90600ddfe2b


# ╔═╡ 3e947f36-69e2-433b-aaa7-4c8c681fd532
isdir("")

# ╔═╡ b8d70c2e-30bc-43f6-965f-c68399ba54d1
md"""
## Data Loading
"""

# ╔═╡ 81805051-a10d-40f8-b128-612a6dec304f
begin 
	obs_profs = LilGuys.StellarProfile[]
	
	for galaxy in galaxies
		filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxy, "density_profiles/jax_2c_profile.toml")

		prof = LilGuys.StellarProfile(filename)
		@assert prof.normalization == "mass"

		filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxy, "observed_properties.toml")
		obs_props = TOML.parsefile(filename)

		dist = obs_props["distance"]
		d_log_r = log10(LilGuys.arcmin_to_kpc(1, dist))

		prof.log_r .+= d_log_r
		prof.log_Sigma .-= 2*d_log_r
		
		push!(obs_profs, prof)
	end
end

# ╔═╡ 71089769-f037-4593-9d77-36c6cdf060ca


# ╔═╡ 529fb229-e633-4c73-a994-86e1130b5c9e
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

# ╔═╡ c51a8cfd-f8a4-495e-8319-d0c362b9a3aa
orbit_props["idx_peris"]

# ╔═╡ 74c6a109-7b1b-4a84-8744-feef39f6f2a6
orbit_props["idx_apos"]

# ╔═╡ 025a0ed4-90f9-469e-80ca-bb5258100eb9
#=╠═╡
df_probs = LilGuys.read_hdf5_table(stars_dir_in * "/probabilities_stars.hdf5")
  ╠═╡ =#

# ╔═╡ 6f8aea0a-91e6-40f4-863f-33a8cca3aee4
#=╠═╡
p_idx = df_probs.index; probabilities = df_probs.probability
  ╠═╡ =#

# ╔═╡ 932fd989-6856-42cf-a6e8-d5453af07c64
#=╠═╡
@assert isperm(p_idx)
  ╠═╡ =#

# ╔═╡ 5ea02b02-1ab2-4ed3-92d6-214a1772c01e
#=╠═╡
out = Output(model_dir, weights=probabilities)
  ╠═╡ =#

# ╔═╡ cc1607f9-083a-49d3-bc73-b274b55c4899
#=╠═╡
out.times[idx_f] * T2GYR
  ╠═╡ =#

# ╔═╡ abe2deee-32e3-43d5-8232-10aeb797bac1
md"""
# Analysis
"""

# ╔═╡ 9f09fb80-333d-4327-a8c6-fa74ffaa4e5d
#=╠═╡
σv = LilGuys.calc_σv_1d(out[idx_f])
  ╠═╡ =#

# ╔═╡ 52943af9-b6c5-4864-8aac-7c65dc2f056d
#=╠═╡
σv * V2KMS
  ╠═╡ =#

# ╔═╡ b30c5830-5bc0-479b-9d53-e614af04c6b4
#=╠═╡
Δt = out.times[idx_f] - filter(x->x < out.times[idx_f], orbit_props["t_peris"])[end]
  ╠═╡ =#

# ╔═╡ cfd02c43-1d4b-4651-9705-31068a0cc7a0
#=╠═╡
Δt * T2GYR
  ╠═╡ =#

# ╔═╡ ef5fa6b6-5111-4a2e-ab84-8452714ebb4a
#=╠═╡
r_b_kpc = LilGuys.calc_break_radius(σv, Δt)
  ╠═╡ =#

# ╔═╡ f4828aa8-011c-4c37-9437-0ecd95adaa68
#=╠═╡
prof_f = LilGuys.StellarProfile(out[idx_f])
  ╠═╡ =#

# ╔═╡ 73522902-8e61-427f-921f-1ef161ea9664
#=╠═╡
prof_i = LilGuys.StellarProfile(out[idx_i])
  ╠═╡ =#

# ╔═╡ 21ace95e-ba98-4728-a8f7-8303f0eac247
#=╠═╡
let
	fig = Figure()
	
	ax = Axis(fig[1,1],
		limits = (-2.5, 1, -5, 2),
		xlabel = L"log $r$ / kpc",
		ylabel = L"log $\Sigma_\star$",
	)


	lines!(prof_i.log_r, prof_i.log_Sigma, label="initial", linestyle=:dot)
	lines!(prof_f.log_r, prof_f.log_Sigma, label="final (illustrative)")

	vlines!(log10(r_b_kpc), color=:black, linestyle=:dash)
	text!(log10.(r_b_kpc) +0.025, -4, text=L"r_b")


	for i in eachindex(galaxies)
		galaxy = galaxies[i]
		prof = obs_profs[i]
		
		errorscatter!(prof.log_r, prof.log_Sigma, yerror=prof.log_Sigma_err, label=galaxy, color=COLORS[i+2])
	end

	LilGuys.hide_grid!(ax)

	axislegend(position=:lb)

	@savefig "scl_umi_vs_idealized" fig
end
  ╠═╡ =#

# ╔═╡ e6d91c8b-648d-416a-ab2e-5e6439f7d5b7
# ╠═╡ disabled = true
#=╠═╡
starsname = inputs.starsname
  ╠═╡ =#

# ╔═╡ 75c691d6-77ab-4715-a645-88352879313c
#=╠═╡
starsname = "exp2d_rs0.10"
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═4b9e4fb8-f257-47c5-b1fa-0b4e9bdd042c
# ╠═75c691d6-77ab-4715-a645-88352879313c
# ╠═6feb26d5-8b1a-4913-95a2-d57c237fdadb
# ╠═258525e6-7fa1-46b3-95f6-950c2b82aafc
# ╠═c51a8cfd-f8a4-495e-8319-d0c362b9a3aa
# ╠═74c6a109-7b1b-4a84-8744-feef39f6f2a6
# ╠═cc1607f9-083a-49d3-bc73-b274b55c4899
# ╟─ea3d8d00-11ee-4b60-9b6a-58a3e396a572
# ╠═5fa9aec4-f876-11ef-1eed-bd8da0982d2b
# ╠═221247cd-e447-400c-8491-b7248787f7a7
# ╠═81684bfc-3316-42c2-a2a9-d9db98858edd
# ╠═01240d81-5b6e-41b6-a2be-32b87a6568ef
# ╠═d143090f-a2e1-4068-a3e7-6ef35635f758
# ╠═047193b4-dcf7-4c36-aafb-2179d8363f3c
# ╠═70ac654f-b45f-4b14-ad0f-01cf5f1e29e5
# ╠═76fdb3ab-60a6-42ea-bb71-0204d37284ce
# ╠═b1edbdfc-298b-40e2-b6c0-687a866fe180
# ╠═e6d91c8b-648d-416a-ab2e-5e6439f7d5b7
# ╠═b51a4ddd-90fd-415f-8a48-e4f35aeee801
# ╠═94f3bf4a-248a-466e-b2de-b90600ddfe2b
# ╠═3e947f36-69e2-433b-aaa7-4c8c681fd532
# ╟─b8d70c2e-30bc-43f6-965f-c68399ba54d1
# ╠═81805051-a10d-40f8-b128-612a6dec304f
# ╠═71089769-f037-4593-9d77-36c6cdf060ca
# ╠═529fb229-e633-4c73-a994-86e1130b5c9e
# ╠═025a0ed4-90f9-469e-80ca-bb5258100eb9
# ╠═932fd989-6856-42cf-a6e8-d5453af07c64
# ╠═6f8aea0a-91e6-40f4-863f-33a8cca3aee4
# ╠═5ea02b02-1ab2-4ed3-92d6-214a1772c01e
# ╟─abe2deee-32e3-43d5-8232-10aeb797bac1
# ╠═9f09fb80-333d-4327-a8c6-fa74ffaa4e5d
# ╠═52943af9-b6c5-4864-8aac-7c65dc2f056d
# ╠═b30c5830-5bc0-479b-9d53-e614af04c6b4
# ╠═cfd02c43-1d4b-4651-9705-31068a0cc7a0
# ╠═ef5fa6b6-5111-4a2e-ab84-8452714ebb4a
# ╠═f4828aa8-011c-4c37-9437-0ecd95adaa68
# ╠═73522902-8e61-427f-921f-1ef161ea9664
# ╠═21ace95e-ba98-4728-a8f7-8303f0eac247
