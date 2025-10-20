### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
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

# ╔═╡ 7caa0c82-b2ed-45e7-a1a9-beb0b20a2e65
using PyFITS

# ╔═╡ 1de31271-dcf3-4513-9786-04822fe99294
md"""
The goal of this script is to compare the initial and final simulation density profiles for a given stars model against the density profile and compare the velocity dispersions against observations. This enables a quick iteration to find the appropriate stellar profile and halo to match present-day stuctural properties.
"""

# ╔═╡ 0ec88057-7dab-45a3-a4c0-cf4cd7d86590
import HDF5

# ╔═╡ b47c355b-48d3-499f-9333-e469878eb4f2
lmc = false

# ╔═╡ ce947244-6184-48a3-9635-5c5ec7734d6b
md"""
# Setup
"""

# ╔═╡ 0d2ebb2b-0942-4f8b-8b4c-bf02ae60e20e
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

# ╔═╡ 2aba5629-9e59-489e-9831-45868365fbed
@bind inputs confirm(notebook_inputs(;
	galaxyname = TextField(default="sculptor"),
	haloname = TextField(default="1e7_new_v31_r3.2"),
	orbitname = TextField(default="orbit_"),
	starsname = TextField(default="exp2d_rs0.13"),
))

# ╔═╡ af83e902-3cdf-4fea-8a41-95daefe4e0df
starsname = inputs.starsname

# ╔═╡ a1e70d0d-f515-4233-9568-fb7261fe494e
modelname = joinpath(inputs.haloname, inputs.orbitname)

# ╔═╡ e78d1592-f001-4710-9daa-826498258374
galaxyname = inputs.galaxyname

# ╔═╡ ee38471c-043c-4571-bcf6-5fd93df84132
md"""
# Data loading
"""

# ╔═╡ a9ca639c-c7ab-4ba4-b47a-4a78ed0b2269
model_stars_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, "stars", starsname)

# ╔═╡ f3754857-3222-474f-9e9c-cb956a966125
FIGDIR = model_stars_dir * "/figures"

# ╔═╡ 5d21d511-8e33-4c25-b9be-57dd8bf6998f
begin
	using LilGuys
	FIGDIR
end

# ╔═╡ 9fecfed4-8aba-42f6-bcf2-1891e0b5f947
prof_obs = LilGuys.SurfaceDensityProfile(
joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "density_profiles", "fiducial_profile.toml")
) |> LilGuys.filter_empty_bins

# ╔═╡ e3fe8cb2-fc6f-4bbd-86e9-d1cd1c296c15
prof_sim_i = LilGuys.SurfaceDensityProfile(model_stars_dir *  "/initial_profile.toml") |> LilGuys.filter_empty_bins

# ╔═╡ 11a8ec0c-0565-4bed-9406-d4944d7ab333
prof_sim = LilGuys.SurfaceDensityProfile(model_stars_dir *  "/final_profile.toml") |> LilGuys.filter_empty_bins

# ╔═╡ 8358c354-29f3-4d62-940f-7085419a1970
if lmc 
	orbital_props_lmc = TOML.parsefile(joinpath(model_stars_dir, "../../orbital_properties_lmc.toml"))
end

# ╔═╡ 604a2e03-e3cc-4c48-80d3-04f10017170f
orbital_props = TOML.parsefile(joinpath(model_stars_dir, "../../orbital_properties.toml"))

# ╔═╡ 2226b7b1-fa78-4428-9cea-3f3a8b17db5f
r_J = let
	r_J = nothing
	if lmc
		filename = joinpath(model_stars_dir, "../../jacobi_lmc.toml")
		if isfile(filename)
			r_J = TOML.parsefile(filename)["r_J"]
		else
			@info "jacobi radius not calculated"
		end
	else
		filename = joinpath(model_stars_dir, "../../jacobi.toml")
		if isfile(filename)
			r_J = TOML.parsefile(filename)["r_J"]
		else
			@info "jacobi radius not calculated"
		end
	end
	r_J
end
	

# ╔═╡ f6df0a6f-ee08-49f1-92e6-24516397bda9
if lmc
	t_last_peri = orbital_props_lmc["t_last_peri"]
else
	t_last_peri = orbital_props["t_last_peri"]
end

# ╔═╡ a55ebe61-f364-4cde-9006-5e022226e110
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml")) |> LilGuys.collapse_errors

# ╔═╡ dc6e60d9-7203-48a5-b277-9555c5773851
times = HDF5.h5open(joinpath(model_stars_dir, "../../centres.hdf5")) do f
	return f["times"][:]
end

# ╔═╡ c4c8bb3c-7460-4946-b025-0f56f5502946
time_f = times[orbital_props["idx_f"]]

# ╔═╡ d0330154-d477-4199-9112-dd6ddd687ed6
time_f * T2GYR

# ╔═╡ 9a9ff46c-2237-4e01-a37d-374dff336f58
sigma_v = prof_sim.annotations["sigma_v"]

# ╔═╡ 9fc456ac-7ab9-4033-b682-d2cb5c386392
r_break = LilGuys.break_radius(t_last_peri/T2GYR, sigma_v)

# ╔═╡ 64aa553d-b024-4991-9265-e0ccbd57d6fb
md"""
# Plots
"""

# ╔═╡ a4eb7940-b145-4db9-b82c-94b46fc44035
md"""
## Profiles
"""

# ╔═╡ 940b085c-e78c-4800-8b3e-bdd13ed970d1
n_center = 13

# ╔═╡ b04a894b-33c2-4daf-9506-cb310ba0c218
R_b = LilGuys.kpc2arcmin(r_break, orbital_props["distance_f"])

# ╔═╡ ab47f03b-4eac-439d-869e-1ee6c2e22ebe
norm_obs = LilGuys.mean(prof_obs.log_Sigma[1:n_center])

# ╔═╡ 60cc712c-fea1-4d19-b97d-d14a6b275a77
CairoMakie.activate!(type=:png)

# ╔═╡ acc308e5-a3e5-46d7-9311-34f0c4ee7c54
ymin = -10

# ╔═╡ 7671fc47-15d4-4970-98a7-b0ad64f1114c
@bind dy confirm(NumberField(-10:0.05:10, default=0))

# ╔═╡ 9a1dd286-f0d0-499a-a421-9e94a0d19bc0
norm_sim = LilGuys.mean(prof_sim.log_Sigma[prof_sim.log_R .< prof_obs.log_R[n_center]]) - dy

# ╔═╡ 54724d9f-0619-4374-abbe-fa600fc5cc92
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log R / arcmin",
		ylabel = L"\log \Sigma / \Sigma_0",
			  limits=(-1, 2.5, ymin, 0.5)
	)

	lines!(log_radii(prof_sim_i), log_surface_density(prof_sim_i) .- norm_sim, label ="initial", color=COLORS[2], linestyle=:dot)
	
	lines!(log_radii(prof_sim), log_surface_density(prof_sim) .- norm_sim, label ="final", color=COLORS[2])

	prof = LilGuys.scale(prof_obs, 1.0, 10^-float(norm_obs))
	LilGuys.plot_log_Σ!(ax,prof, label ="data", markersize=3, color=:black, )

	axislegend(position=:lb)

	if @isdefined R_b
		arrows!([log10(R_b)], [-7], [0], [-0.5])
		text!([log10(R_b)], [-7], text=L"R_b", fontsize=8)
	end
	
	if !isnothing(r_J)
		arrows!([log10(r_J)], [-7], [0], [-0.5], color=COLORS[1])
		text!([log10(r_J)], [-7], text=L"R_j", fontsize=8, color=COLORS[1])
	end
	
	@savefig "Sigma_vs_obs_i_f"
	fig
end

# ╔═╡ 5876165d-6232-436b-aea2-a9636f1557b6


# ╔═╡ e86ff7fe-b424-4edf-bf23-264a327e3926
md"""
## Velocity Dispersions
"""

# ╔═╡ 41bcae3a-d1ea-4620-a2d0-aa9eddfe8554
readdir(model_stars_dir)

# ╔═╡ 7bc27716-2c02-47a3-9e04-b9d55a526000
scalars = read_fits(joinpath(model_stars_dir, "stellar_profiles_3d_scalars.fits"))

# ╔═╡ 250f4593-5ac1-472b-bdc1-95fd28ce2466
orbital_props["distance_f"]

# ╔═╡ b8796e05-b149-4aed-ad68-c61f56c41093
σ_obs = obs_props["sigma_v"]

# ╔═╡ 5205de59-6cda-47c1-a988-83b6d5f52b87
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "time", ylabel = L"$\sigma_v$ / km\,s$^{-1}$")

	hlines!(middle(σ_obs), color=:black)
	hspan!(σ_obs.middle - σ_obs.lower, σ_obs.middle + σ_obs.upper, alpha=0.2, color=:black)
	hlines!(middle(σ_obs) - 1, color=:black, linestyle=:dash)


	if @isdefined scalars
		lines!(scalars.time * T2GYR .- time_f*T2GYR, scalars.sigma_v * V2KMS)
	end
	scatter!(time_f * T2GYR, sigma_v)
	
	fig
	@savefig "sigma_v_evolution"
end

# ╔═╡ Cell order:
# ╟─1de31271-dcf3-4513-9786-04822fe99294
# ╠═2aba5629-9e59-489e-9831-45868365fbed
# ╠═7caa0c82-b2ed-45e7-a1a9-beb0b20a2e65
# ╠═0ec88057-7dab-45a3-a4c0-cf4cd7d86590
# ╠═b47c355b-48d3-499f-9333-e469878eb4f2
# ╟─ce947244-6184-48a3-9635-5c5ec7734d6b
# ╠═2820df64-e001-11ef-08e2-d730b6110af3
# ╟─0d2ebb2b-0942-4f8b-8b4c-bf02ae60e20e
# ╠═af83e902-3cdf-4fea-8a41-95daefe4e0df
# ╠═f3754857-3222-474f-9e9c-cb956a966125
# ╠═a1e70d0d-f515-4233-9568-fb7261fe494e
# ╠═e78d1592-f001-4710-9daa-826498258374
# ╠═5d21d511-8e33-4c25-b9be-57dd8bf6998f
# ╟─ee38471c-043c-4571-bcf6-5fd93df84132
# ╠═a9ca639c-c7ab-4ba4-b47a-4a78ed0b2269
# ╠═9fecfed4-8aba-42f6-bcf2-1891e0b5f947
# ╠═e3fe8cb2-fc6f-4bbd-86e9-d1cd1c296c15
# ╠═11a8ec0c-0565-4bed-9406-d4944d7ab333
# ╠═8358c354-29f3-4d62-940f-7085419a1970
# ╠═604a2e03-e3cc-4c48-80d3-04f10017170f
# ╠═2226b7b1-fa78-4428-9cea-3f3a8b17db5f
# ╠═f6df0a6f-ee08-49f1-92e6-24516397bda9
# ╠═a55ebe61-f364-4cde-9006-5e022226e110
# ╠═dc6e60d9-7203-48a5-b277-9555c5773851
# ╠═c4c8bb3c-7460-4946-b025-0f56f5502946
# ╠═d0330154-d477-4199-9112-dd6ddd687ed6
# ╠═9a9ff46c-2237-4e01-a37d-374dff336f58
# ╠═9fc456ac-7ab9-4033-b682-d2cb5c386392
# ╟─64aa553d-b024-4991-9265-e0ccbd57d6fb
# ╟─a4eb7940-b145-4db9-b82c-94b46fc44035
# ╠═940b085c-e78c-4800-8b3e-bdd13ed970d1
# ╠═b04a894b-33c2-4daf-9506-cb310ba0c218
# ╠═9a1dd286-f0d0-499a-a421-9e94a0d19bc0
# ╠═ab47f03b-4eac-439d-869e-1ee6c2e22ebe
# ╠═60cc712c-fea1-4d19-b97d-d14a6b275a77
# ╠═acc308e5-a3e5-46d7-9311-34f0c4ee7c54
# ╠═7671fc47-15d4-4970-98a7-b0ad64f1114c
# ╠═54724d9f-0619-4374-abbe-fa600fc5cc92
# ╠═5876165d-6232-436b-aea2-a9636f1557b6
# ╟─e86ff7fe-b424-4edf-bf23-264a327e3926
# ╠═41bcae3a-d1ea-4620-a2d0-aa9eddfe8554
# ╠═7bc27716-2c02-47a3-9e04-b9d55a526000
# ╠═250f4593-5ac1-472b-bdc1-95fd28ce2466
# ╠═5205de59-6cda-47c1-a988-83b6d5f52b87
# ╠═b8796e05-b149-4aed-ad68-c61f56c41093
