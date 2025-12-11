### A Pluto.jl notebook ###
# v0.20.21

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

# ╔═╡ 3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:svg)

# ╔═╡ e93918f6-0c72-4dec-9cdf-97f185c0bceb
module ModelUtils
	include("model_utils.jl")
end

# ╔═╡ 437ac44b-5440-4b3c-aa58-4bbd93be7174
function plot_prof!(prof; kwargs...)
	lines!(LilGuys.log_radii(prof), LilGuys.log_surface_density(prof); kwargs...)
end

# ╔═╡ 56365cfd-16e0-4567-a85c-f59372a611c7
function compare_density(gs, modelname_exp; modelname_plummer=nothing, limits=(-0.5, 2.5, -6, 2.5), lmc=false, y_low_R_h=nothing, y_low=nothing, color=COLORS[3])

	galaxy = modelname_exp[1]
	ax = Axis(gs,
		xlabel = "log Radii / arcmin",
		# ylabel = "log surface density",
		limits = limits,
			  title=Dict("sculptor"=>"Sculptor", "ursa_minor" => "Ursa Minor")[galaxy]
			 )
	ModelUtils.plot_expected_profile!(galaxy)

	prof_i, prof_f, _ = ModelUtils.load_stellar_profiles(modelname_exp...)

	plot_prof!(prof_i, color=color, linewidth=1.5, linestyle=:dot, label="stars initial")
	plot_prof!(prof_f, color=color, linewidth=1.5, linestyle=:solid, label="stars final")


	if !isnothing(modelname_plummer)
		prof_i, prof_f, _ = ModelUtils.load_stellar_profiles(modelname_plummer...)
	
		plot_prof!(prof_i, color=COLORS[2], linestyle=:dot, label="2exp initial")
		plot_prof!(prof_f, color=COLORS[2], linestyle=:solid, label="2exp final")
	end


	limits!(limits...)
	if isnothing(y_low)
		y_low = limits[3]
	end
	r_b = ModelUtils.get_r_b(modelname_exp..., lmc=lmc)
	ModelUtils.plot_r_break_arrow!(r_b, y_low)

	r_j = ModelUtils.get_r_j(modelname_exp[1:2]..., lmc=lmc)
	ModelUtils.plot_r_jacobi_arrow!(r_j, y_low)

	R_h = ModelUtils.get_R_h(galaxy)
	if isnothing(y_low_R_h)
		y_low_R_h = limits[3]
	end
	ModelUtils.plot_R_h_arrow!(R_h, y_low_R_h)
	
	ax
end

# ╔═╡ 7948bbd7-2dda-46c8-a850-d80944ba9096
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ fc0a118f-f8b7-425e-aacf-c8c7630d61a5
function model_label!(text)
	text!(1, 1, space=:relative, text=text, align=(:right, :top), offset=(-4, -4), fontsize=8)
end

# ╔═╡ dbed5ec4-6030-4667-a87b-e8100deebe4a
@savefig "density_i_f" let
	fig = Figure(size=(3.5, 3.5) .* 72)


	ax_scl = compare_density(fig[1,1], modelnames["scl_smallperi"])
	hidexdecorations!(ticks=false, minorticks=false)
	model_label!("exponential")
	
	ax_scl_plummer = compare_density(fig[2,1], modelnames["scl_smallperi_plummer"], y_low_R_h=-3, y_low=-3)
	# hidexdecorations!(ticks=false, minorticks=false)
	model_label!("Plummer")
	ax_scl_plummer.title[] = ""


	ax_umi = compare_density(fig[1, 2], modelnames["umi_smallperi"], )
	hidedecorations!(ticks=false, minorticks=false)
	model_label!("exponential")

	ax_umi_plummer = compare_density(fig[2, 2], modelnames["umi_smallperi_plummer"], y_low=-4)
	hideydecorations!(ticks=false, minorticks=false)
	model_label!("Plummer")
	ax_umi_plummer.title[] = ""
	
	linkaxes!(ax_scl, ax_scl_plummer, ax_umi, ax_umi_plummer)
	ax_scl_plummer.xticks[] = -0.5:0.5:2

	axislegend(ax_scl_plummer, position=:lb)

	Label(fig[:, 0], L"log $\Sigma$ / stars arcmin$^{-2}$", rotation=π/2)
	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 0)

	rowsize!(fig.layout, 1, Aspect(1, 1))
	rowsize!(fig.layout, 2, Aspect(1, 1))

	resize_to_layout!()
	fig

end

# ╔═╡ 83f2045f-97f1-43ed-84c6-75518ed2beeb
@savefig "density_i_f_2exp" let
	fig = Figure(size=(3.5, 3.5) .* 72)


	ax_scl = compare_density(fig[1,1], modelnames["scl_smallperi"], modelname_plummer=modelnames["scl_smallperi_2exp"])
	hidexdecorations!(ticks=false, minorticks=false)
	


	ax_umi = compare_density(fig[1, 2], modelnames["umi_smallperi"], modelname_plummer=modelnames["umi_smallperi_2exp"])
	hidedecorations!(ticks=false, minorticks=false)

	
	linkaxes!(ax_scl, ax_umi)

	axislegend(ax_scl, position=:lb)

	Label(fig[:, 0], "log surface density", rotation=π/2)
	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 0)

	rowsize!(fig.layout, 1, Aspect(1, 1))

	resize_to_layout!()
	fig

end

# ╔═╡ d26db585-3e68-4870-99da-bb5d9182c923
@savefig "scl_lmc_density_i_f" let
	fig = Figure(size=(3.5, 3.5) .* 72)


	ax_scl = compare_density(fig[1,1], modelnames["scl_lmc"], y_low=-6, y_low_R_h=-3)
	model_label!("exponential")
	ax_scl.title[] = ""
	ax_scl.xticks[] = -0.5:0.5:2.0

	ax_scl_plummer = compare_density(fig[1,2], modelnames["scl_lmc_plummer"], )
	# hidexdecorations!(ticks=false, minorticks=false)
	model_label!("Plummer")
	hideydecorations!(ticks=false, minorticks=false)

	ax_scl_plummer.title[] = ""

	
	axislegend(ax_scl, position=:lb)


	Label(fig[:, 0], L"log $\Sigma$ / stars arcmin$^{-2}$", rotation=π/2, fontsize=10)
	Label(fig[0, :], "Sculptor: MW+LMC", font=:bold, fontsize=1.2 * theme(:fontsize)[])
	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 0)
	rowsize!(fig.layout, 1, Aspect(1, 1))
	fig

end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═e93918f6-0c72-4dec-9cdf-97f185c0bceb
# ╠═56365cfd-16e0-4567-a85c-f59372a611c7
# ╠═437ac44b-5440-4b3c-aa58-4bbd93be7174
# ╠═7948bbd7-2dda-46c8-a850-d80944ba9096
# ╠═fc0a118f-f8b7-425e-aacf-c8c7630d61a5
# ╠═dbed5ec4-6030-4667-a87b-e8100deebe4a
# ╠═83f2045f-97f1-43ed-84c6-75518ed2beeb
# ╠═d26db585-3e68-4870-99da-bb5d9182c923
