### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 2820df64-e001-11ef-08e2-d730b6110af3
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya
	import TOML
	using PlutoUI
end

# ╔═╡ 1de31271-dcf3-4513-9786-04822fe99294
md"""
The goal of this script is to compare the initial and final simulation density profiles for a given stars model against the observations
"""

# ╔═╡ e78d1592-f001-4710-9daa-826498258374
@bind galaxyname TextField(default="sculptor")

# ╔═╡ af83e902-3cdf-4fea-8a41-95daefe4e0df
@bind starsname TextField(default="exp2d_rs0.13")

# ╔═╡ 511e3e8d-8278-4e7b-998c-097b1513320e
@bind modelname TextField()

# ╔═╡ f01ad318-2ca3-430a-a121-7370b03085db
@bind run_program CheckBox()

# ╔═╡ ee38471c-043c-4571-bcf6-5fd93df84132
md"""
# Data loading
"""

# ╔═╡ a9ca639c-c7ab-4ba4-b47a-4a78ed0b2269
model_stars_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, "stars", starsname)

# ╔═╡ f3754857-3222-474f-9e9c-cb956a966125
figdir = model_stars_dir * "/figures"

# ╔═╡ 5d21d511-8e33-4c25-b9be-57dd8bf6998f
begin
	using LilGuys
	figdir
end

# ╔═╡ 9fecfed4-8aba-42f6-bcf2-1891e0b5f947
if run_program
	if isdir(model_stars_dir)
		mkpath(figdir)
	end
	
	prof_sim = LilGuys.StellarProfile(model_stars_dir *  "/final_profile.toml")
	prof_sim_i = LilGuys.StellarProfile(model_stars_dir *  "/initial_profile.toml")
	prof_obs = LilGuys.StellarProfile(
	joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "density_profiles", "fiducial_profile.toml")
)
end

# ╔═╡ 64aa553d-b024-4991-9265-e0ccbd57d6fb
md"""
# Plots
"""

# ╔═╡ 9a1dd286-f0d0-499a-a421-9e94a0d19bc0
norm_sim = prof_sim.log_Sigma[2]

# ╔═╡ ab47f03b-4eac-439d-869e-1ee6c2e22ebe
norm_obs = prof_obs.log_Sigma[2]

# ╔═╡ 54724d9f-0619-4374-abbe-fa600fc5cc92
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log r / arcmin",
		ylabel = L"\log \Sigma / \Sigma_0",
		limits = (-0.5, 2.5, -8, 1),
	)

	lines!(prof_sim.log_r, prof_sim.log_Sigma .- norm_sim, label = "simulation" )
	lines!(prof_sim_i.log_r, prof_sim_i.log_Sigma .- norm_sim, label = "simulation (initial)", color=COLORS[3], linestyle=:dot)
	
	errscatter!(prof_obs.log_r, prof_obs.log_Sigma .- norm_obs, yerr = prof_obs.log_Sigma_err,
		color=COLORS[2],
		label = "observation",
	)

	axislegend()

	@savefig "Sigma_vs_obs_i_f"
	fig
end

# ╔═╡ 2759953d-65c3-4645-a271-4996a83b62f1
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log r / arcmin",
		ylabel = L"\log \Sigma / \Sigma_0",
		limits = (-0.5, 2.5, -5, 0.3),
	)

	lines!(prof_sim.log_r, prof_sim.log_Sigma .- norm_sim, label = "simulation" )
	lines!(prof_sim_i.log_r, prof_sim_i.log_Sigma .- norm_sim, label = "simulation (initial)", color=COLORS[3], linestyle=:dot)
	
	errscatter!(prof_obs.log_r, prof_obs.log_Sigma .- norm_obs, yerr = prof_obs.log_Sigma_err,
		color=COLORS[2],
		label = "observation",
	)

	axislegend()

	fig
end

# ╔═╡ Cell order:
# ╠═1de31271-dcf3-4513-9786-04822fe99294
# ╠═2820df64-e001-11ef-08e2-d730b6110af3
# ╠═e78d1592-f001-4710-9daa-826498258374
# ╠═af83e902-3cdf-4fea-8a41-95daefe4e0df
# ╠═f3754857-3222-474f-9e9c-cb956a966125
# ╠═511e3e8d-8278-4e7b-998c-097b1513320e
# ╠═f01ad318-2ca3-430a-a121-7370b03085db
# ╠═5d21d511-8e33-4c25-b9be-57dd8bf6998f
# ╟─ee38471c-043c-4571-bcf6-5fd93df84132
# ╠═a9ca639c-c7ab-4ba4-b47a-4a78ed0b2269
# ╠═9fecfed4-8aba-42f6-bcf2-1891e0b5f947
# ╟─64aa553d-b024-4991-9265-e0ccbd57d6fb
# ╠═9a1dd286-f0d0-499a-a421-9e94a0d19bc0
# ╠═ab47f03b-4eac-439d-869e-1ee6c2e22ebe
# ╠═54724d9f-0619-4374-abbe-fa600fc5cc92
# ╠═2759953d-65c3-4645-a271-4996a83b62f1
