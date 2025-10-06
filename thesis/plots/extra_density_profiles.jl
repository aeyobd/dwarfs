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

# ╔═╡ e7ae2044-a31e-4a98-b496-3b1b3be056f0
using OrderedCollections

# ╔═╡ 53a18bc2-3bbd-4fe3-8a77-d9c49d2ddea9
include("paper_style.jl")

# ╔═╡ 0ac92e62-f8fd-4e17-8b2f-feecdab75497
import TOML

# ╔═╡ 322e01f3-32d2-4408-91ef-ab7777f0f925
module Utils
	include("gaia_utils.jl")
end

# ╔═╡ ad8414dd-9428-47db-b08d-65940091413e
CairoMakie.activate!(type=:png)

# ╔═╡ 54e449b3-6bf8-45d4-98da-5eb1c6e7e6ec
lw = theme(:linewidth)[]/2

# ╔═╡ b11f3588-9904-46ef-8f35-a7c1d623e020
log_r_ell_label = L"$\log\,R_\textrm{ell}$\,/\,arcmin"

# ╔═╡ d42a0cd3-cc8e-4a24-8887-7100f3927961
log_Σ_label = L"$\log\,\Sigma$\,/\,stars\ arcmin$^{-2}$"

# ╔═╡ 9c594ea7-8670-4235-8bde-cb2d670fe2c3
function get_R_h(galaxyname)
	obs_props = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/$galaxyname/observed_properties.toml")
	R_h = obs_props["R_h"]
end

# ╔═╡ 5d005025-8d6e-4d06-8b38-8aedd6189dfe
function load_profile(galaxyname, algname="jax")
	filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, 
		"density_profiles/$(algname)_profile.toml")
    prof = LilGuys.SurfaceDensityProfile(filename)

	prof
end

# ╔═╡ 4f854674-f062-4cb2-a026-f9bac6863990
function get_background(galaxyname, algname="best_eqw")
	filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, 
		"density_profiles/$(algname)_background.toml")
    prof = TOML.parsefile(filename)

	return prof["log_Sigma"], prof["log_Sigma_err"]
end

# ╔═╡ 78a7391b-511a-45b8-b4e1-52b0bb675058
# update_theme!(
# 	ErrorScatter=(; 
# 		cycle=Cycle([:color, :marker], covary=true)
# 		))

# ╔═╡ 7bfc460d-4ac2-4d02-aa3f-d93447fa71a6
function get_extrema(profiles)
	x_l = Inf
	x_h = -Inf
	y_l = Inf
	y_h = -Inf
	for (_, prof) in profiles
		log_bins = log_radii_bins(prof)
		if length(log_bins) > 0
			x_l = min(log_bins[1], x_l)
			x_h = max(log_bins[end], x_h)
		end

		ys = log_surface_density(prof)
		y_l = min(y_l, minimum(middle.(ys) .- lower_error.(ys)))
		y_h = max(y_h, maximum(float.(log_surface_density(prof))))
	end
	return x_l, x_h, y_l, y_h
end

# ╔═╡ 20b94af8-aaad-4832-8643-5d5d1f00e288
function compare_densities!(ax, profiles; 
		R_b_R_h=nothing, styles=nothing,
		prof_ana = nothing, y_R_b=nothing, R_h=nothing,
							y_min=-5
	)
	x_l, x_h, y_l, y_h = get_extrema(profiles)


	
	for (i, (name, prof)) in enumerate(profiles)
		if styles !== nothing
			kwargs = styles[name]
		else
			kwargs = (;)
		end

		filt = lower_error.(prof.log_Sigma) .< 1
		ys = prof.log_Sigma[filt]
		kwargs = Dict(k => kwargs[k] for k in propertynames(kwargs))
		if :label ∉ keys(kwargs)
			kwargs[:label] = name
		end
		scatter!(ax, prof.log_R, prof.log_Sigma; 
			#yerror=LilGuys.error_interval.(prof.log_Sigma), 
			kwargs...)
	end

	if prof_ana !== nothing
		x = LinRange(x_l, x_h, 1000)
	
		y = @. log10(surface_density(prof_ana, 10 .^ x))
		lines!(ax, x, y, color=:black, linewidth=lw)
	end

	
	# Plot annotations


	# if R_h !== nothing
	# 	x = log10(R_h)
	# 	y = y_l
	# 	annotation!(0, 36, x, y_min,linewidth=lw, text=L"R_h")
	# end


end

# ╔═╡ 3f13f074-e6b3-4b77-8080-bf8861b1e018
function compare_densities(profiles; kwargs...)
	
	fig = Figure()
	compare_densities(fig[1,1], profiles; kwargs...)

	fig

end

# ╔═╡ cc48ec77-c92c-40d9-8569-6550f6dde552
function compare_densities(gs, profiles; y_min=-5, kwargs...)
	ax = Axis(gs,  
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
		limits = (nothing, nothing, y_min, nothing)
	)

	compare_densities!(ax, profiles; y_min=y_min, kwargs...)

	ax
end

# ╔═╡ 183d7d89-b34d-4be8-9b7f-be039a6ce545
function scale!(dict, key, R_scale, M_scale)
	dict[key] = LilGuys.scale(dict[key], R_scale, M_scale)
end

# ╔═╡ 1160f544-3879-4892-a47f-fc88d3620b8a
md"""
# Sanity and validation plots
"""

# ╔═╡ 5917e45e-50bb-441f-9d12-6901025a81f3
profiles_scl_j24 = OrderedDict(
	"2-exp" => load_profile("sculptor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
	"1-exp" => load_profile("sculptor", "jax_eqw") |> LilGuys.filter_empty_bins,
	"simple" => load_profile("sculptor", "simple") |> LilGuys.filter_empty_bins,
	# "PM+CMD" => load_profile("sculptor", "jax_LLR_0_eqw") |> LilGuys.filter_empty_bins,

)

# ╔═╡ 66761f62-bfe3-4c5e-8630-17a8b63721e7
begin 
	profiles_scl_extra = OrderedDict(
		"2-exp" => load_profile("sculptor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
		"circ" => load_profile("sculptor", "jax_circ_eqw") |> LilGuys.filter_empty_bins,	

		"bright" => load_profile("sculptor", "jax_bright") |> LilGuys.filter_empty_bins,
		"DELVE" => load_profile("sculptor", "delve_new_sub"),
	)

	scale!(profiles_scl_extra, "bright", 1, 2)
	scale!(profiles_scl_extra, "DELVE", 1, 0.2)

	profiles_scl_extra
end

# ╔═╡ 8544dd0c-556f-45e2-b432-3258b257ad99
begin 
	profiles_umi_j24 = OrderedDict(
		"2-exp" => load_profile("ursa_minor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
		"1-exp" => load_profile("ursa_minor", "jax_eqw") |> LilGuys.filter_empty_bins,
		# "PM+CMD" => load_profile("ursa_minor", "jax_LLR_0_eqw") |> LilGuys.filter_empty_bins,
		"simple" => load_profile("ursa_minor", "simple") |> LilGuys.filter_empty_bins,


	)


	scale!(profiles_umi_j24, "simple", 1, 1.25)

end

# ╔═╡ ad8ab6ab-5a63-49f4-a60e-7b6130abfbad
begin 
	profiles_umi_extra = OrderedDict(
		"2-exp" => load_profile("ursa_minor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
		"circ" => load_profile("ursa_minor", "jax_circ_eqw") |> LilGuys.filter_empty_bins,
		"bright" => load_profile("ursa_minor", "jax_bright") |> LilGuys.filter_empty_bins,

		"UNIONS" => load_profile("ursa_minor", "unions_sub"),

	)

	scale!(profiles_umi_extra, "bright", 1, 2)
	scale!(profiles_umi_extra, "UNIONS", 1, 0.25)

	profiles_umi_extra
end

# ╔═╡ aeb9af8b-3129-4fbf-b2ef-970f6d76e605
r_limit_scl = 64.1

# ╔═╡ fcefddf3-1059-44fe-aa45-a9b53cf81e55
r_limit_umi = 86.4

# ╔═╡ 5609b085-5792-4d83-9395-81f977edf7cd
scale_theme_element!(:markersize, 2/3)

# ╔═╡ 9ab4431d-e64f-44c9-a641-f16a524f3c19
α = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 9c4c10a1-6161-4c2f-ba10-05509cabeb49
R_h_scl = get_R_h("sculptor")

# ╔═╡ 05eba4d6-f566-4e27-9bed-72a7520dfef0
ana_scl = LilGuys.Exp2D(R_s=R_h_scl/α, M=0.9sum(profiles_scl_extra["2-exp"].counts))

# ╔═╡ c691c460-2037-4559-9ef0-0c1619eafe3b
R_h_umi = get_R_h("ursa_minor")

# ╔═╡ 0278b08b-c7e5-48e0-a2fc-ebdfb46fb9ab
ana_umi = LilGuys.Exp2D(R_s=R_h_umi/α, M=0.8sum(profiles_umi_extra["2-exp"].counts))

# ╔═╡ a1109c5b-ea82-40ac-8194-c7dcc4a069ea
@savefig "density_methods_extra" let
	fig = Figure()

	
	ax1 = compare_densities(fig[1,1], profiles_scl_j24, R_h=R_h_scl, prof_ana=ana_scl)
	ax1.title = "Sculptor"
	axislegend(position=:lb)
	hidexdecorations!(ticks=false, minorticks=false)
	ax1.ylabel = ""

	vspan!(log10(r_limit_scl), 3, alpha=0.1, color=:black)
	text!((1.6), -1.6, text="2D exp", rotation=-0.95π/3, color=(:black), fontsize=theme(:fontsize)[] * 0.8, align=(:left, :top))

	
	ax2 = compare_densities(fig[2,1], profiles_scl_extra, R_h=R_h_scl, prof_ana=ana_scl)
	xlims!(-0.0, 2.1)
	axislegend(position=:lb)
	ax2.ylabel = ""
	vspan!(log10(r_limit_scl), 3, alpha=0.1, color=:black)

	Label(fig[:, 0], "log surface density", rotation=π/2)




	ax1_umi = compare_densities(fig[1,2], profiles_umi_j24, R_h=R_h_umi, prof_ana=ana_umi)
	axislegend(position=:lb)
	hidedecorations!(ticks=false, minorticks=false)
	ax1_umi.title = "Ursa Minor"
	vspan!(log10(r_limit_umi), 3, alpha=0.1, color=:black)

	text!(log10(r_limit_umi), -1, text="BG-limited", rotation=π/2, color=(:black, 0.8), fontsize=theme(:fontsize)[] * 0.8, align=(:left, :top))

	ax2_umi = compare_densities(fig[2,2], profiles_umi_extra, R_h=R_h_umi, prof_ana=ana_umi)
	axislegend(position=:lb)
	hideydecorations!(ticks=false, minorticks=false)
	vspan!(log10(r_limit_umi), 3, alpha=0.1, color=:black)

	
	linkaxes!(ax1, ax2)
	linkaxes!(ax1_umi, ax2_umi)
	linkyaxes!(ax1, ax1_umi, ax2, ax2_umi)
	xlims!(0.3, 2.3)
	ylims!(-4.5, 1.8)
	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 0)
	fig
end

# ╔═╡ e9f641fe-4ec6-495b-9634-a9be98395a4d
md"""
# Comparisons
"""

# ╔═╡ 16034a2d-574e-422a-978b-0add93e3b128
begin 
	profiles_umi_sanity = OrderedDict(
		"fiducial" => load_profile("ursa_minor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
		"unions" => load_profile("ursa_minor", "unions"),
		"m18" => load_profile("ursa_minor", "munoz+18"),
		"s23" => load_profile("ursa_minor", "sestito+23"),
	)

	# scale!(profiles_umi_sanity, "unions", 1, 0.2)
	scale!(profiles_umi_sanity, "m18", sqrt(1-0.55), 0.2*sqrt(0.44)^2)
	scale!(profiles_umi_sanity, "s23", sqrt(1-0.55), sqrt(1-0.55)^2)

end

# ╔═╡ 3d12e962-3f0a-447a-8af4-81ab0e7a34c3
compare_densities(profiles_umi_sanity, R_h=R_h_umi, prof_ana=ana_umi)

# ╔═╡ c29da90b-32b5-445b-9de4-03cba2349397
begin 
	profiles_scl_sanity = OrderedDict(
		"fiducial" => load_profile("sculptor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
		# "DELVE" => load_profile("sculptor", "delve_rgb_sub"),
		"m18" => load_profile("sculptor", "munoz+18"),
		"s23" => load_profile("sculptor", "sestito+23"),
	)

	scale!(profiles_scl_sanity, "m18", sqrt(1-0.37), 0.3*sqrt(1-0.37)^2)
	scale!(profiles_scl_sanity, "s23", sqrt(1-0.37), sqrt(1-0.37)^2)

	profiles_scl_j24
end

# ╔═╡ 523186ec-f7c8-4968-afd5-167ca5f06254
compare_densities(profiles_scl_sanity, R_h=R_h_scl, prof_ana=ana_scl)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═0ac92e62-f8fd-4e17-8b2f-feecdab75497
# ╠═322e01f3-32d2-4408-91ef-ab7777f0f925
# ╠═e7ae2044-a31e-4a98-b496-3b1b3be056f0
# ╠═ad8414dd-9428-47db-b08d-65940091413e
# ╠═53a18bc2-3bbd-4fe3-8a77-d9c49d2ddea9
# ╠═54e449b3-6bf8-45d4-98da-5eb1c6e7e6ec
# ╠═b11f3588-9904-46ef-8f35-a7c1d623e020
# ╠═d42a0cd3-cc8e-4a24-8887-7100f3927961
# ╠═9c594ea7-8670-4235-8bde-cb2d670fe2c3
# ╠═5d005025-8d6e-4d06-8b38-8aedd6189dfe
# ╠═4f854674-f062-4cb2-a026-f9bac6863990
# ╠═78a7391b-511a-45b8-b4e1-52b0bb675058
# ╠═7bfc460d-4ac2-4d02-aa3f-d93447fa71a6
# ╠═20b94af8-aaad-4832-8643-5d5d1f00e288
# ╠═3f13f074-e6b3-4b77-8080-bf8861b1e018
# ╠═cc48ec77-c92c-40d9-8569-6550f6dde552
# ╠═183d7d89-b34d-4be8-9b7f-be039a6ce545
# ╟─1160f544-3879-4892-a47f-fc88d3620b8a
# ╠═5917e45e-50bb-441f-9d12-6901025a81f3
# ╠═66761f62-bfe3-4c5e-8630-17a8b63721e7
# ╠═8544dd0c-556f-45e2-b432-3258b257ad99
# ╠═ad8ab6ab-5a63-49f4-a60e-7b6130abfbad
# ╠═aeb9af8b-3129-4fbf-b2ef-970f6d76e605
# ╠═fcefddf3-1059-44fe-aa45-a9b53cf81e55
# ╠═5609b085-5792-4d83-9395-81f977edf7cd
# ╠═9ab4431d-e64f-44c9-a641-f16a524f3c19
# ╠═9c4c10a1-6161-4c2f-ba10-05509cabeb49
# ╠═05eba4d6-f566-4e27-9bed-72a7520dfef0
# ╠═c691c460-2037-4559-9ef0-0c1619eafe3b
# ╠═0278b08b-c7e5-48e0-a2fc-ebdfb46fb9ab
# ╠═a1109c5b-ea82-40ac-8194-c7dcc4a069ea
# ╟─e9f641fe-4ec6-495b-9634-a9be98395a4d
# ╠═16034a2d-574e-422a-978b-0add93e3b128
# ╠═3d12e962-3f0a-447a-8af4-81ab0e7a34c3
# ╠═523186ec-f7c8-4968-afd5-167ca5f06254
# ╠═c29da90b-32b5-445b-9de4-03cba2349397
