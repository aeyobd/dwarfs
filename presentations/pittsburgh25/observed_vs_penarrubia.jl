### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 5fa9aec4-f876-11ef-1eed-bd8da0982d2b
begin
	using Pkg; Pkg.activate()

	using CairoMakie
	using Arya

end

# ╔═╡ 7c31c70f-9ade-407d-8db0-34ff94a43ffd
using CSV, DataFrames

# ╔═╡ 01240d81-5b6e-41b6-a2be-32b87a6568ef
include("./style.jl")

# ╔═╡ ea3d8d00-11ee-4b60-9b6a-58a3e396a572
md"""
# Setup
"""

# ╔═╡ 221247cd-e447-400c-8491-b7248787f7a7
import TOML

# ╔═╡ d143090f-a2e1-4068-a3e7-6ef35635f758
galaxies = ["sculptor", "ursa_minor", "fornax"]

# ╔═╡ b51a4ddd-90fd-415f-8a48-e4f35aeee801
FIGDIR = "figures"

# ╔═╡ 81684bfc-3316-42c2-a2a9-d9db98858edd
using LilGuys; FIGDIR

# ╔═╡ b8d70c2e-30bc-43f6-965f-c68399ba54d1
md"""
## Data Loading
"""

# ╔═╡ 9a7adce4-a829-42c4-9902-02b5d4e1891d
p08_notides = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/all/data", "penarrubia08_notides.csv"), DataFrame)

# ╔═╡ 39d1e497-a359-4d58-8f83-8aff1eeeb1f2
p08_tides = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/all/data", "penarrubia08_apo.csv"), DataFrame)

# ╔═╡ d3780086-bc13-4b35-a960-e7aa51f2433e
α = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 81805051-a10d-40f8-b128-612a6dec304f
begin 
	obs_profs = LilGuys.SurfaceDensityProfile[]
	
	for galaxy in galaxies
		dir = joinpath(ENV["DWARFS_ROOT"], "observations", galaxy, "density_profiles")
		filepath = ""
		for filename in ["jax_2c_eqw_profile.toml", "jax_eqw_profile.toml"]
			if isfile(joinpath(dir, filename))
				filepath = joinpath(dir, filename)
				break
			end
		end
		

		prof = LilGuys.SurfaceDensityProfile(filepath) |> LilGuys.filter_empty_bins
		@assert prof.log_R_scale == 0.0
		obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxy, "observed_properties.toml"))
		
		R_h = obs_props["R_h_inner"]

		# filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxy, "observed_properties.toml")
		# obs_props = TOML.parsefile(filename)
		Σ_h = 10 .^ LilGuys.lerp(prof.log_R, middle.(prof.log_Sigma))(log10(R_h))
		M_s = Σ_h * R_h .^ 2

    	prof = LilGuys.scale(prof, 1/R_h, 1/M_s)
	
		push!(obs_profs, prof)
	end
end

# ╔═╡ abe2deee-32e3-43d5-8232-10aeb797bac1
md"""
# Analysis
"""

# ╔═╡ 940c75b5-e19b-474f-8913-a489d770372a
r_scale = 0.7

# ╔═╡ b0e33ee3-ae02-400a-b468-33de44a44c06
m_scale = 5 * r_scale

# ╔═╡ 16697307-ef7b-49c7-9c96-9b9804f99dd3
theme(:size)

# ╔═╡ ef59af37-893f-4864-a01b-7c5b58afbaf0
CairoMakie.activate!(type=:png)

# ╔═╡ 0732a179-aca0-4bf4-8be7-be1e713ee18e
begin 
	prof_i = LilGuys.SurfaceDensityProfile(
		log_R=p08_notides.x,
		log_Sigma=p08_notides." y",
		R_units="", 
		log_R_bins=[]
	) 
	
	prof_f = LilGuys.SurfaceDensityProfile(
		log_R=p08_tides.x,
		log_Sigma=p08_tides." y",
		R_units="", 
		log_R_bins=[]
	)

	prof_i = LilGuys.scale(prof_i, r_scale, m_scale)

	prof_f = LilGuys.scale(prof_f, r_scale, m_scale)

end

# ╔═╡ 940f98f7-d1a2-4a37-9da5-aab82968561e
log_r_label = L"log radius / $R_h$"

# ╔═╡ c00be232-10ad-4905-a924-1e80b4da345b
log_sigma_label = L"log surface density / $\Sigma_{\star}(R_h)$"

# ╔═╡ 57d7cd00-9d14-43d9-86d4-2be8d03cfcba
let
	fig = Figure()
	
	ax = Axis(fig[1,1],
		limits = (-1.9, 1.2, -4.5, 1.5),
		xlabel = log_r_label,
		ylabel =log_sigma_label,
	)



	for i in eachindex(galaxies)
		galaxy = galaxies[i]
		prof = obs_profs[i]

		errorscatter!(prof.log_R, prof.log_Sigma, yerror=LilGuys.error_interval.(prof.log_Sigma), 
					  label=Dict("sculptor"=> "Sculptor", "ursa_minor"=>"Ursa Minor", "fornax" => "Fornax")[galaxy], 
					  color=[COLORS[4], COLORS[2], COLORS[1]][i])
	end



	axislegend(position=:lb, margin=(19,19,19,19))

	@savefig "scl_umi_vs_fornax" fig
end

# ╔═╡ d0c1d099-f07b-440e-8192-573b95c9354b


# ╔═╡ 2ab5af72-4667-42e5-8e3c-20a68952d2e6
let
	fig = Figure()
	
	ax = Axis(fig[1,1],
		limits = (-1.9, 1.2, -4.5, 1.5),
		xlabel = log_r_label,
		ylabel =log_sigma_label,
	)

	
	
	for i in eachindex(galaxies)
		galaxy = galaxies[i]
		prof = obs_profs[i]

		errorscatter!(prof.log_R, prof.log_Sigma, yerror=LilGuys.error_interval.(prof.log_Sigma), 
					  label=Dict("sculptor"=> "Sculptor", "ursa_minor"=>"Ursa Minor", "fornax" => "Fornax")[galaxy], 
					  color=[COLORS[4], COLORS[2], COLORS[1]][i])
	end

prof_ana = LilGuys.Exp2D(R_s = 1/LilGuys.R_h(LilGuys.Exp2D()))
	x = LinRange(-1.9, 1.2, 1000)
	y = @. log10(LilGuys.surface_density(prof_ana, 10^x))
	y .-= log10(LilGuys.surface_density(prof_ana, 1))

	lines!(x, y, color=:black)

	
	axislegend(position=:lb, margin=(19,19,19,19))

	@savefig "scl_umi_vs_exp" fig
end

# ╔═╡ 21ace95e-ba98-4728-a8f7-8303f0eac247
let
	fig = Figure()
	
	ax = Axis(fig[1,1],
		limits = (-1.9, 1.2, -4.5, 1.5),
		xlabel = log_r_label,
		ylabel =log_sigma_label,
	)



	for i in eachindex(galaxies)
		galaxy = galaxies[i]
		prof = obs_profs[i]

		errorscatter!(prof.log_R, prof.log_Sigma, yerror=LilGuys.error_interval.(prof.log_Sigma), 
					  label=Dict("sculptor"=> "Sculptor", "ursa_minor"=>"Ursa Minor", "fornax" => "Fornax")[galaxy], 
					  color=[COLORS[4], COLORS[2], COLORS[1]][i])
	end


	lines!(prof_i.log_R, prof_i.log_Sigma, label="no tides", linestyle=:dot, color=:black)
	lines!(prof_f.log_R, prof_f.log_Sigma, label="tides", linewidth=6, color=:black)

	axislegend(position=:lb, margin=(19,19,19,19))

	@savefig "scl_umi_vs_penarrubia" fig
end

# ╔═╡ Cell order:
# ╟─ea3d8d00-11ee-4b60-9b6a-58a3e396a572
# ╠═5fa9aec4-f876-11ef-1eed-bd8da0982d2b
# ╠═221247cd-e447-400c-8491-b7248787f7a7
# ╠═81684bfc-3316-42c2-a2a9-d9db98858edd
# ╠═01240d81-5b6e-41b6-a2be-32b87a6568ef
# ╠═d143090f-a2e1-4068-a3e7-6ef35635f758
# ╠═b51a4ddd-90fd-415f-8a48-e4f35aeee801
# ╟─b8d70c2e-30bc-43f6-965f-c68399ba54d1
# ╠═7c31c70f-9ade-407d-8db0-34ff94a43ffd
# ╠═9a7adce4-a829-42c4-9902-02b5d4e1891d
# ╠═39d1e497-a359-4d58-8f83-8aff1eeeb1f2
# ╠═d3780086-bc13-4b35-a960-e7aa51f2433e
# ╠═81805051-a10d-40f8-b128-612a6dec304f
# ╟─abe2deee-32e3-43d5-8232-10aeb797bac1
# ╠═940c75b5-e19b-474f-8913-a489d770372a
# ╠═b0e33ee3-ae02-400a-b468-33de44a44c06
# ╠═16697307-ef7b-49c7-9c96-9b9804f99dd3
# ╠═ef59af37-893f-4864-a01b-7c5b58afbaf0
# ╠═0732a179-aca0-4bf4-8be7-be1e713ee18e
# ╠═940f98f7-d1a2-4a37-9da5-aab82968561e
# ╠═c00be232-10ad-4905-a924-1e80b4da345b
# ╠═57d7cd00-9d14-43d9-86d4-2be8d03cfcba
# ╠═d0c1d099-f07b-440e-8192-573b95c9354b
# ╠═2ab5af72-4667-42e5-8e3c-20a68952d2e6
# ╠═21ace95e-ba98-4728-a8f7-8303f0eac247
