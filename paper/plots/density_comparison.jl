### A Pluto.jl notebook ###
# v0.20.19

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
include("./paper_style.jl")

# ╔═╡ 1f1780d7-c6ed-4a92-b263-f8e43e9fab68
md"""
Monolithic figure creation haha

"""

# ╔═╡ 3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:svg)

# ╔═╡ e93918f6-0c72-4dec-9cdf-97f185c0bceb
module ModelUtils
	include("model_utils.jl")
end

# ╔═╡ 7948bbd7-2dda-46c8-a850-d80944ba9096
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ 437ac44b-5440-4b3c-aa58-4bbd93be7174
function plot_prof!(prof; kwargs...)
	lines!(LilGuys.log_radii(prof), LilGuys.log_surface_density(prof); kwargs...)
end

# ╔═╡ 56365cfd-16e0-4567-a85c-f59372a611c7
function compare_exp_plummer(modelname_exp, modelname_plummer; limits=(-0.5, 2.5, -4, 2), lmc=false)

	galaxy = modelname_exp[1]
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log Radii / arcmin",
		ylabel = "log surface density",
		limits = limits,
			  title=Dict("sculptor"=>"Sculptor", "ursa_minor" => "Ursa Minor")[galaxy]
			 )
	ModelUtils.plot_expected_profile!(galaxy)

	prof_i, prof_f, _ = ModelUtils.load_stellar_profiles(modelname_exp...)

	plot_prof!(prof_i, color=COLORS[3], linewidth=1.5, linestyle=:dot, label="exponential final")
	plot_prof!(prof_f, color=COLORS[3], linewidth=1.5, linestyle=:solid, label="exponential initial")

	prof_i, prof_f, _ = ModelUtils.load_stellar_profiles(modelname_plummer...)

	plot_prof!(prof_i, color=COLORS[2], linestyle=:dot, label="Plummer final")
	plot_prof!(prof_f, color=COLORS[2], linestyle=:solid, label="Plummer initial")


	y_low = limits[3]
	r_b = ModelUtils.get_r_b(modelname_exp..., lmc=lmc)
	ModelUtils.plot_r_break_arrow!(r_b, y_low)

	r_j = ModelUtils.get_r_j(modelname_exp[1:2]..., lmc=lmc)
	ModelUtils.plot_r_jacobi_arrow!(r_j, y_low)

	R_h = ModelUtils.get_R_h(galaxy)
	ModelUtils.plot_R_h_arrow!(R_h, y_low)
	axislegend(position=:lb)
	
	fig
end

# ╔═╡ dbed5ec4-6030-4667-a87b-e8100deebe4a
@savefig "scl_density_i_f" compare_exp_plummer(modelnames["scl_lmc"], modelnames["scl_lmc_plummer"], lmc=true)

# ╔═╡ 7132c6c8-8b28-4d3a-9cac-b0bfcfa09cec
@savefig "umi_density_i_f" compare_exp_plummer(modelnames["umi_smallperi"], modelnames["umi_smallperi_plummer"], limits=(-0.5, 2.5, -4.5, 1.5))

# ╔═╡ Cell order:
# ╠═1f1780d7-c6ed-4a92-b263-f8e43e9fab68
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═e93918f6-0c72-4dec-9cdf-97f185c0bceb
# ╠═56365cfd-16e0-4567-a85c-f59372a611c7
# ╠═437ac44b-5440-4b3c-aa58-4bbd93be7174
# ╠═7948bbd7-2dda-46c8-a850-d80944ba9096
# ╠═dbed5ec4-6030-4667-a87b-e8100deebe4a
# ╠═7132c6c8-8b28-4d3a-9cac-b0bfcfa09cec
