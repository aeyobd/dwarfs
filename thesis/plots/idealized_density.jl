### A Pluto.jl notebook ###
# v0.20.18

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
modelname = "1e6/orbit_10_100"

# ╔═╡ 75c691d6-77ab-4715-a645-88352879313c
starsname = "exp2d_rs0.25"

# ╔═╡ 6feb26d5-8b1a-4913-95a2-d57c237fdadb
idx_f = 47

# ╔═╡ 258525e6-7fa1-46b3-95f6-950c2b82aafc
idx_i = 1

# ╔═╡ ea3d8d00-11ee-4b60-9b6a-58a3e396a572
md"""
# Setup
"""

# ╔═╡ 221247cd-e447-400c-8491-b7248787f7a7
import TOML

# ╔═╡ d143090f-a2e1-4068-a3e7-6ef35635f758
galaxies = ["sculptor", "ursa_minor", "fornax"]

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
stars_dir_in = joinpath(model_dir, "../stars/$starsname")

# ╔═╡ b1edbdfc-298b-40e2-b6c0-687a866fe180
stars_dir_out = joinpath(model_dir, "stars/$starsname")

# ╔═╡ b51a4ddd-90fd-415f-8a48-e4f35aeee801
FIGDIR = "figures"

# ╔═╡ 81684bfc-3316-42c2-a2a9-d9db98858edd
using LilGuys; FIGDIR

# ╔═╡ b8d70c2e-30bc-43f6-965f-c68399ba54d1
md"""
## Data Loading
"""

# ╔═╡ a0af7cb2-5ed6-4bd4-afc4-7d4481111127
α = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 470165c2-c675-4ad6-89e6-ec43eaab2301
function get_R_h(galaxy)
	props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxy, "observed_properties.toml"))
	return props["R_h"]
end

# ╔═╡ 81805051-a10d-40f8-b128-612a6dec304f
begin 
	obs_profs = LilGuys.SurfaceDensityProfile[]
		
	for galaxy in galaxies
		dir = joinpath(ENV["DWARFS_ROOT"], "observations", galaxy, "density_profiles")
		filepath = ""
		for filename in ["jax_2c_eqw_profile.toml", "jax_eqw_profile.toml", "jax_eqw_profile.toml"]
			if isfile(joinpath(dir, filename))
				filepath = joinpath(dir, filename)
				break
			end
		end
		

		prof = LilGuys.SurfaceDensityProfile(filepath) |> LilGuys.filter_empty_bins
		@assert prof.log_R_scale == 0.0
		R_h = get_R_h(galaxy)
		
		# density_fit = TOML.parsefile(filepath * "_inner_fits.toml")
		# R_h = 10 .^ density_fit["log_R_s_exp2d_inner"] * α

		# filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxy, "observed_properties.toml")
		# obs_props = TOML.parsefile(filename)
		Σ_h = 10 .^ LilGuys.lerp(prof.log_R, middle.(prof.log_Sigma))(log10(R_h))
		M_s = Σ_h * R_h .^ 2

    	prof = LilGuys.scale(prof, 1/R_h, 1/M_s)
	
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
df_probs = LilGuys.read_hdf5_table(stars_dir_in * "/probabilities_stars.hdf5")

# ╔═╡ 6f8aea0a-91e6-40f4-863f-33a8cca3aee4
p_idx = df_probs.index; probabilities = df_probs.probability

# ╔═╡ 932fd989-6856-42cf-a6e8-d5453af07c64
@assert isperm(p_idx)

# ╔═╡ 5ea02b02-1ab2-4ed3-92d6-214a1772c01e
out = Output(model_dir, weights=probabilities)

# ╔═╡ cc1607f9-083a-49d3-bc73-b274b55c4899
out.times[idx_f] * T2GYR

# ╔═╡ 3263692b-d918-4926-850b-8b127df24f20
T2GYR * (out.times[idx_f] - orbit_props["t_peris"][1])

# ╔═╡ abe2deee-32e3-43d5-8232-10aeb797bac1
md"""
# Analysis
"""

# ╔═╡ 9f09fb80-333d-4327-a8c6-fa74ffaa4e5d
σv = LilGuys.σv_1d(out[idx_f])

# ╔═╡ 52943af9-b6c5-4864-8aac-7c65dc2f056d
σv * V2KMS

# ╔═╡ b30c5830-5bc0-479b-9d53-e614af04c6b4
Δt = out.times[idx_f] - filter(x->x < out.times[idx_f], orbit_props["t_peris"])[end]

# ╔═╡ cfd02c43-1d4b-4651-9705-31068a0cc7a0
Δt * T2GYR

# ╔═╡ ef5fa6b6-5111-4a2e-ab84-8452714ebb4a
r_b_kpc = LilGuys.break_radius(σv, Δt)

# ╔═╡ f4828aa8-011c-4c37-9437-0ecd95adaa68
prof_f = LilGuys.SurfaceDensityProfile(out[idx_f])

# ╔═╡ 73522902-8e61-427f-921f-1ef161ea9664
prof_i = LilGuys.SurfaceDensityProfile(out[idx_i])

# ╔═╡ 80bb56d3-37f0-49fc-a261-4965032326fc
labels = Dict("sculptor"=> "Sculptor", "ursa_minor"=>"Ursa Minor", "fornax" => "Fornax")

# ╔═╡ b59b1601-a2fe-4481-a027-28b1e391c338
colors = Dict("sculptor" => COLORS[2], "ursa_minor" => COLORS[4], "fornax" => COLORS[3])

# ╔═╡ 2e3c0c9c-d7de-4ef8-aad1-972ed5957810
exp_profile = LilGuys.Sersic(n=1)

# ╔═╡ 00ed6c32-b602-4c8c-ae13-7b323c6232fd
fig_limits = (-1.0, 1.2, -5, 1)

# ╔═╡ 10a88d85-08b1-4401-bf9a-4873b9428676
@savefig "scl_umi_vs_fornax" let
	fig = Figure()
	
	ax = DensityAxis(fig[1,1])


	plot_exponential!(ax)

	plot_galaxy_densities!(ax)

	axislegend(position=:lb)


	fig
end

# ╔═╡ 3c2d3aa7-9937-46ad-a393-32c120b670b9
function DensityAxis(gs)
	ax = Axis(gs,
		xlabel = L"log $R$ / $R_h$",
		ylabel = L"log $\Sigma_\star\ /\ \Sigma_h$",
		limits = fig_limits
	)
end

# ╔═╡ 8d3ecf85-afbf-4b53-8669-b8390e8419db
function plot_exponential!(ax)
	prof = LilGuys.Sersic(R_h=1, n=1)

	x = LinRange(-1, 1, 1000)
	y = log10.(LilGuys.surface_density.(prof, 10 .^ x)) .- log10(LilGuys.surface_density(prof, 1))
	lines!(ax, x, y, color=:black, label="exponential")
end

# ╔═╡ 48b059bd-415e-47e5-b498-37e7b6ba21ad
function plot_galaxy_densities!(ax)
	for i in eachindex(galaxies)
		galaxy = galaxies[i]
		prof = obs_profs[i]

		errorscatter!(ax, prof.log_R, prof.log_Sigma, yerror=LilGuys.error_interval.(prof.log_Sigma), 
					  label=labels[galaxy], 
					  color=colors[galaxy],  marker=[:rect, :utriangle, :circle][i])
	end
end

# ╔═╡ 1d54de8d-c9be-42b6-85a1-235f7d8a08d0
@savefig "scl_umi_vs_idealized" let
	fig = Figure()
	
	ax = DensityAxis(fig[1,1])

	sim_scale = 0.45
	sim_R_scale = 0.3

	plot_galaxy_densities!(ax)



	lines!(prof_i.log_R .+ sim_R_scale, prof_i.log_Sigma .+ sim_scale, label="simulation initial", linestyle=:dot)
	lines!(prof_f.log_R .+ sim_R_scale, prof_f.log_Sigma .+ sim_scale, label="simulation final", color=:black)

	x = log10(r_b_kpc) + sim_R_scale
	y = LilGuys.lerp(LilGuys.log_radii(prof_f) .+ sim_R_scale, middle.(LilGuys.log_surface_density(prof_f)))(x) + sim_scale
	# dy = 0.25
	# h = 0.75
	# arrows2d!([x  + sim_R_scale], [y + sim_scale + h], [0], [-h], align=0.5, minshaftlength=0)

	annotation!(0, 50, x, y, text=L"r_\textrm{break}")


	x = LinRange(-2, 1, 1000)
	y = @. log10(LilGuys.surface_density.(exp_profile, 10^x))
	#lines!(x, y)

	axislegend(position=:lb)

	fig
end

# ╔═╡ Cell order:
# ╠═4b9e4fb8-f257-47c5-b1fa-0b4e9bdd042c
# ╠═75c691d6-77ab-4715-a645-88352879313c
# ╠═6feb26d5-8b1a-4913-95a2-d57c237fdadb
# ╠═258525e6-7fa1-46b3-95f6-950c2b82aafc
# ╠═c51a8cfd-f8a4-495e-8319-d0c362b9a3aa
# ╠═74c6a109-7b1b-4a84-8744-feef39f6f2a6
# ╠═cc1607f9-083a-49d3-bc73-b274b55c4899
# ╠═3263692b-d918-4926-850b-8b127df24f20
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
# ╠═b51a4ddd-90fd-415f-8a48-e4f35aeee801
# ╟─b8d70c2e-30bc-43f6-965f-c68399ba54d1
# ╠═a0af7cb2-5ed6-4bd4-afc4-7d4481111127
# ╠═470165c2-c675-4ad6-89e6-ec43eaab2301
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
# ╠═80bb56d3-37f0-49fc-a261-4965032326fc
# ╠═b59b1601-a2fe-4481-a027-28b1e391c338
# ╠═2e3c0c9c-d7de-4ef8-aad1-972ed5957810
# ╠═00ed6c32-b602-4c8c-ae13-7b323c6232fd
# ╠═10a88d85-08b1-4401-bf9a-4873b9428676
# ╠═3c2d3aa7-9937-46ad-a393-32c120b670b9
# ╠═8d3ecf85-afbf-4b53-8669-b8390e8419db
# ╠═48b059bd-415e-47e5-b498-37e7b6ba21ad
# ╠═1d54de8d-c9be-42b6-85a1-235f7d8a08d0
