### A Pluto.jl notebook ###
# v0.20.18

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

# ╔═╡ 5fc7e171-0fce-40f7-bc22-2aa02bb3694f
using CSV, DataFrames

# ╔═╡ 53a18bc2-3bbd-4fe3-8a77-d9c49d2ddea9
include("paper_style.jl")

# ╔═╡ 0ac92e62-f8fd-4e17-8b2f-feecdab75497
import TOML

# ╔═╡ ad8414dd-9428-47db-b08d-65940091413e
CairoMakie.activate!(type=:png)

# ╔═╡ 54e449b3-6bf8-45d4-98da-5eb1c6e7e6ec
lw = theme(:linewidth)[]/2

# ╔═╡ b11f3588-9904-46ef-8f35-a7c1d623e020
log_r_ell_label = L"$\log\,R_\textrm{ell}$\,/\,arcmin"

# ╔═╡ d42a0cd3-cc8e-4a24-8887-7100f3927961
log_Σ_label = L"$\log\,\Sigma$\,/\,stars\ arcmin$^{-2}$"

# ╔═╡ 389fefd6-60fd-4dd8-ba77-29e87b4ed846
obs_props_scl = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml")

# ╔═╡ 322e01f3-32d2-4408-91ef-ab7777f0f925
module Utils
	include("utils.jl")
end

# ╔═╡ fd74ebb7-625e-4396-a91a-37045a283a08
obs_props_umi = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/observed_properties.toml")

# ╔═╡ d50cc41b-9465-4481-86d9-c8cf221e1978
obs_props_fornax = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/fornax/observed_properties.toml")

# ╔═╡ 9ab4431d-e64f-44c9-a641-f16a524f3c19
α = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 9c594ea7-8670-4235-8bde-cb2d670fe2c3
function get_R_h(obs_props)	
	R_h = obs_props["R_h"]
end

# ╔═╡ 933c2e0a-db49-46fd-b3ae-a40be190022d
R_munoz_1 = 36.58

# ╔═╡ f3b219b4-7cc7-48a0-9675-91de7f40bb6d
function get_M_s(galaxyname)
	fit =  TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "density_profiles/jax_eqw_inner_fits.toml"))

	@info 10^fit["log_R_s_exp2d_inner"] * α
	return 10 ^ fit["log_M_exp2d_inner"]
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
update_theme!(
	ErrorScatter=(; 
		cycle=Cycle([:color, :marker], covary=true)
		))

# ╔═╡ c51c6f13-b5ed-4671-9f08-4adb4f7b38b7
import Statistics: median

# ╔═╡ ffaddffc-bbd0-4a08-9422-71bf872eb750
df_peris = CSV.read(ENV["DWARFS_ROOT"] * "/analysis/all_galaxies/vasiliev+21/properties_w_orbits.csv", DataFrame)


# ╔═╡ 50c692f7-93df-4edc-af0b-bb901a1bc907
function get_R_break(galaxyname, obs_props)
	dt = df_peris.t_last_peri[df_peris.galaxyname .== galaxyname] |> only
	dist = obs_props["distance"]
	sigma_v = obs_props["sigma_v"] / V2KMS
	r_b_kpc = LilGuys.break_radius(sigma_v, dt)
	R_b_arcmin = LilGuys.kpc2arcmin(r_b_kpc, dist)

	return R_b_arcmin
end

# ╔═╡ abdaea6d-504f-4708-8185-c63dcf90bad9
ana = LilGuys.Sersic(n=1)

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

	if prof_ana !== nothing
		x = LinRange(x_l, x_h, 1000)
	
		y = @. log10(surface_density(prof_ana, 10 .^ x))
		lines!(ax, x, y, label="Exp2D", color=:black, linewidth=lw)
	end
	
	for (name, prof) in profiles
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
		errorscatter!(ax, prof.log_R, prof.log_Sigma; 
			yerror=LilGuys.error_interval.(prof.log_Sigma), 
			kwargs...)
	end

	# Plot annotations


	if R_h !== nothing
		x = log10(R_h)
		y = y_l
		annotation!(ax, 0, 36, x, y_min, color=:black, linewidth=lw, text=L"R_h")
	end


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

# ╔═╡ 30a60b7d-630a-4fa3-ac57-e873655ad754
function compare_densities_residuals(profiles; 
		R_b_R_h=nothing, styles=nothing, reference=collect(keys(profiles))[begin],
		prof_ana = nothing, y_R_b=nothing, R_h=nothing,
	)
	
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
		limits=(nothing, nothing, -5, nothing)
	)

	compare_densities!(ax, profiles; R_b_R_h=R_b_R_h, styles=styles, 	
		prof_ana=prof_ana, y_R_b=y_R_b, R_h=R_h)

	
	# Residuals axis
	axislegend(position=:lb)
	ax2 = Axis(fig[2,1],
		xlabel = log_r_ell_label,
		ylabel = L"\log\,\Sigma / \Sigma_\textrm{Exp2D}",
		limits=(nothing, nothing, -2, 2)
	)

	hlines!(0, color=:black, linewidth=lw)

	for (key, prof) in profiles
		ym = log10.(surface_density.(prof_ana, radii(prof)))
		filt = isfinite.(prof.log_Sigma)

		x = prof.log_R[filt]
		y = LilGuys.middle.(prof.log_Sigma[filt]) .- ym[filt]
		ye = LilGuys.error_interval.(prof.log_Sigma)[filt]

		if styles !== nothing
			kwargs = styles[key]
		else
			kwargs = (;)
		end
		errorscatter!(x, y, yerror=ye; kwargs...)
	end


	linkxaxes!(ax, ax2)
	hidexdecorations!(ax, ticks=false, minorticks=false)
	
	rowgap!(fig.layout, 0.)
	rowsize!(fig.layout, 2, Relative(1/4))

	fig
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
	"fiducial" => load_profile("sculptor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
	"1-component" => load_profile("sculptor", "jax_eqw") |> LilGuys.filter_empty_bins,
	"circ" => load_profile("sculptor", "jax_circ_eqw") |> LilGuys.filter_empty_bins,	
	"MCMC" => load_profile("sculptor", "../mcmc/hist_fast"),
)

# ╔═╡ 66761f62-bfe3-4c5e-8630-17a8b63721e7
begin 
	profiles_scl_extra = OrderedDict(
		"fiducial" => load_profile("sculptor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
		"simple" => load_profile("sculptor", "simple") |> LilGuys.filter_empty_bins,
		"bright" => load_profile("sculptor", "jax_bright") |> LilGuys.filter_empty_bins,
		"DELVE" => load_profile("sculptor", "delve_rgb_sub"),
	)

	scale!(profiles_scl_extra, "bright", 1, 2)
	scale!(profiles_scl_extra, "DELVE", 1, 0.6)

	profiles_scl_extra
end

# ╔═╡ c29da90b-32b5-445b-9de4-03cba2349397
begin 
	profiles_scl_sanity = OrderedDict(
		"fiducial" => load_profile("sculptor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
		"DELVE" => load_profile("sculptor", "delve_rgb_sub"),
		"m18" => load_profile("sculptor", "munoz+18"),
		"s23" => load_profile("sculptor", "sestito+23"),
	)

	scale!(profiles_scl_sanity, "m18", sqrt(1-0.37), 0.3*sqrt(1-0.37)^2)
	scale!(profiles_scl_sanity, "s23", sqrt(1-0.37), sqrt(1-0.37)^2)

	profiles_scl_j24
end

# ╔═╡ 16034a2d-574e-422a-978b-0add93e3b128
begin 
	profiles_umi_sanity = OrderedDict(
		"fiducial" => load_profile("ursa_minor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
		"unions" => load_profile("ursa_minor", "unions_sub"),
		"m18" => load_profile("ursa_minor", "munoz+18"),
		"s23" => load_profile("ursa_minor", "sestito+23"),
	)

	scale!(profiles_umi_sanity, "unions", 1, 0.2)
	scale!(profiles_umi_sanity, "m18", sqrt(1-0.55), 0.2*sqrt(0.44)^2)
	scale!(profiles_umi_sanity, "s23", sqrt(1-0.55), sqrt(1-0.55)^2)

end

# ╔═╡ ad8ab6ab-5a63-49f4-a60e-7b6130abfbad
begin 
	profiles_umi_extra = OrderedDict(
		"fiducial" => load_profile("ursa_minor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
		"simple" => load_profile("ursa_minor", "simple") |> LilGuys.filter_empty_bins,
		"bright" => load_profile("ursa_minor", "jax_bright") |> LilGuys.filter_empty_bins,
		"UNIONS" => load_profile("ursa_minor", "unions_sub"),

	)

	scale!(profiles_umi_extra, "bright", 1, 2)
	scale!(profiles_umi_extra, "UNIONS", 1, 0.2)

	profiles_umi_extra
end

# ╔═╡ 8544dd0c-556f-45e2-b432-3258b257ad99
begin 
	profiles_umi_j24 = OrderedDict(
		"fiducial" => load_profile("ursa_minor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
		"1-component" => load_profile("ursa_minor", "jax_eqw") |> LilGuys.filter_empty_bins,
		"circ" => load_profile("ursa_minor", "jax_circ_eqw") |> LilGuys.filter_empty_bins,
		"MCMC" => load_profile("ursa_minor", "../mcmc/hist_fast"),

	)
end

# ╔═╡ e9f641fe-4ec6-495b-9634-a9be98395a4d
md"""
# Comparisons
"""

# ╔═╡ 27abef4c-62cf-4c82-a76f-38b735569328
begin 
	profiles_scl = OrderedDict(
		"all" => load_profile("sculptor", "best_eqw") |> LilGuys.filter_empty_bins,
		"BG subtracted" => load_profile("sculptor", "best_eqw_sub"),
		"CMD + PM" => load_profile("sculptor", "jax_LLR_0_eqw") |> LilGuys.filter_empty_bins,
		"probable members" => load_profile("sculptor", "fiducial") |> LilGuys.filter_empty_bins,
	)


	profiles_scl
end

# ╔═╡ 92db291e-afb7-490d-a599-d4bcb686dccb
get_extrema(profiles_scl)

# ╔═╡ ef3ad471-392d-435a-950b-1b34b2c3b469
grey = Utils.SELECTION_COLORS[1]

# ╔═╡ cc8c0231-d568-4ac8-a090-428fdea08696
plot_kwargs = Dict(
	"probable members" => (;
		color=Utils.SELECTION_COLORS[3],
		label = L"fiducial ($P_\textrm{sat} > 0.2$)",
		marker = :rect
	),
	"CMD + PM" => (;
		color = Utils.SELECTION_COLORS[2],
		marker = :utriangle,
		label = "CMD + PM"
	),
	"all" => (;
		color = grey,
		marker = :circle,
		label = "all",
	),
	"BG subtracted" => (;
		color = COLORS[2],
		marker = :pentagon,
		label = L"all $ – $ background"
	)
)

# ╔═╡ 0a029aaa-c23d-4126-bd21-323d978e1548
R_b_scl = get_R_break("sculptor", obs_props_scl)

# ╔═╡ 72c8e61d-03af-4bac-ae60-0f0d0ff19492
R_b_umi = get_R_break("ursa_minor", obs_props_umi)

# ╔═╡ b14b8b05-2d08-4e6b-bb85-dcd01038518f
R_b_fornax = get_R_break("fornax", obs_props_fornax)

# ╔═╡ 9c4c10a1-6161-4c2f-ba10-05509cabeb49
R_h_scl = get_R_h(obs_props_scl)

# ╔═╡ 05eba4d6-f566-4e27-9bed-72a7520dfef0
ana_scl = LilGuys.Exp2D(R_s=R_h_scl / α, M=get_M_s("sculptor"))

# ╔═╡ a1109c5b-ea82-40ac-8194-c7dcc4a069ea
@savefig "scl_density_methods_extra" let
	fig = Figure()

	
	ax2 = compare_densities(fig[1,1], profiles_scl_j24, R_h=R_h_scl, prof_ana=ana_scl)
	axislegend(position=:lb)
	hidexdecorations!(ticks=false, minorticks=false)

	ax1 = compare_densities(fig[2,1], profiles_scl_extra, R_h=R_h_scl, prof_ana=ana_scl)
	axislegend(position=:lb)


	linkaxes!(ax1, ax2)
	xlims!(-0.8, 2.3)
	ylims!(-5, 2)
	rowgap!(fig.layout, 0)
	fig
end

# ╔═╡ 523186ec-f7c8-4968-afd5-167ca5f06254
compare_densities(profiles_scl_sanity, R_h=R_h_scl, prof_ana=ana_scl)

# ╔═╡ 5772d3d1-0ff8-4a6a-a474-b96bc7d67347
R_h_umi = get_R_h(obs_props_umi)

# ╔═╡ 469acf5a-b7bf-4e53-bbef-6d57ab2988e9
ana_umi = LilGuys.Exp2D(R_s=R_h_umi / α, M=get_M_s("ursa_minor"))

# ╔═╡ 3d12e962-3f0a-447a-8af4-81ab0e7a34c3
compare_densities_residuals(profiles_umi_sanity, R_h=R_h_umi, prof_ana=ana_umi)

# ╔═╡ 5fa49bc9-a34d-4ae9-9cd3-4051281bc68f
@savefig "umi_density_methods_extra" let
	fig = Figure()

	
	ax1 = compare_densities(fig[1,1], profiles_umi_extra, R_h=R_h_umi, prof_ana=ana_umi)
	axislegend(position=:lb)
	hidexdecorations!(ticks=false, minorticks=false)

	ax2 = compare_densities(fig[2,1], profiles_umi_j24, R_h=R_h_umi, prof_ana=ana_umi)
	axislegend(position=:lb)

	linkaxes!(ax1, ax2)
	xlims!(-0.5, 2.3)
	rowgap!(fig.layout, 0)

	fig
end

# ╔═╡ 9002ba29-45ef-46dc-8632-e993c9ec2dc8
@savefig "umi_density_methods_j24" compare_densities_residuals(profiles_umi_j24, R_h=R_h_umi, prof_ana=ana_umi)

# ╔═╡ a2a77ae2-3a6b-46b2-bc5d-34a0e82fcef9
R_h_fornax = get_R_h(obs_props_fornax)

# ╔═╡ 0b33dfaa-055a-44b2-a463-1233047d624b
ana_fornax = LilGuys.Exp2D(R_s=R_h_fornax / α, M=get_M_s("fornax"))

# ╔═╡ fba7941e-8175-47df-8801-7d1f91211376
R_h_munoz1 = 0.49  # arcmin, munoz+2012

# ╔═╡ 7c8f5d6e-d53f-49a2-8b3a-239117c7c533
Utils.SELECTION_COLORS[[1, 2, 3]]

# ╔═╡ 58346276-5061-4d9a-8781-c1aba26d2c42
grey, COLORS[3], COLORS[1]

# ╔═╡ c070eb97-b7d8-45e8-ae95-0219da2cd3e6
begin 
	profiles_umi = OrderedDict(
		"all" => load_profile("ursa_minor", "best_eqw") |> LilGuys.filter_empty_bins,
		"BG subtracted" => load_profile("ursa_minor", "best_eqw_sub"),
		"CMD + PM" => load_profile("ursa_minor", "jax_LLR_0_eqw") |> LilGuys.filter_empty_bins,
		"probable members" => load_profile("ursa_minor", "fiducial") |> LilGuys.filter_empty_bins,
	)

	profiles_umi
end

# ╔═╡ 7c64d26a-bddd-4159-98d2-3a0f8d2125f5
begin 
	profiles_fornax = OrderedDict(
		"all" => load_profile("fornax", "best_eqw") |> LilGuys.filter_empty_bins,
		"BG subtracted" => load_profile("fornax", "best_eqw_sub") ,
		"CMD + PM" => load_profile("fornax", "jax_LLR_0_eqw") |> LilGuys.filter_empty_bins,
		"probable members" => load_profile("fornax", "fiducial") |> LilGuys.filter_empty_bins,
	)

end

# ╔═╡ 0459f053-b627-427e-a698-a49ea8e68b3d
@savefig "for_density_methods_j24" compare_densities_residuals(profiles_fornax, R_h=R_h_fornax, prof_ana=ana_fornax)

# ╔═╡ 19be0190-af5d-479b-b859-8f1c1e46bba1
@savefig "scl_umi_fnx_density_methods" let
	fig = Figure(size=(5.39, 5.39) .*72)
	ax = Axis(fig[1,1],
		#xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
	)

	y_min = -3
	
	hlines!(get_background("sculptor")[1], color=grey, linewidth=1)
	compare_densities!(ax, profiles_scl, styles=plot_kwargs, R_h=R_h_scl, y_R_b=-2, y_min=y_min)

	text!(0.8, 0.8, text="Sculptor", space=:relative)

	hidexdecorations!(ticks=false, minorticks=false)

	
	# Ursa Minor
	ax_umi = Axis(fig[2,1],
		#xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
	)

	hlines!(get_background("ursa_minor")[1], color=grey, linewidth=1)
	
	compare_densities!(ax_umi, profiles_umi, 
		styles=plot_kwargs,
		R_h=R_h_umi, y_min=y_min
	)

	text!(0.8, 0.8, text="Ursa Minor", space=:relative)

	#munoz 1
	lines!([log10(R_munoz_1 - 3R_h_munoz1), log10(R_munoz_1 + 3R_h_munoz1)], [-2.5, -2.5], color=COLORS[5])
	text!(ax_umi, log10(R_munoz_1), -2.5, text="Muñoz 1", color=COLORS[5], align=(:center, :bottom))
	
	axislegend(position=:lb)
	hidexdecorations!(ticks=false, minorticks=false)



	# Fornax
	ax_fnx = Axis(fig[3, 1],
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
	)
	
	hlines!(get_background("fornax")[1], color=grey, linewidth=1)

	compare_densities!(ax_fnx, profiles_fornax, styles=plot_kwargs, 	
		 R_h=R_h_fornax, y_min=y_min
	)
	
	text!(0.8, 0.8, text="Fornax", space=:relative)


	linkxaxes!(ax, ax_umi, ax_fnx)
	rowgap!(fig.layout, 0.)
	fig

end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═0ac92e62-f8fd-4e17-8b2f-feecdab75497
# ╠═ad8414dd-9428-47db-b08d-65940091413e
# ╠═53a18bc2-3bbd-4fe3-8a77-d9c49d2ddea9
# ╠═54e449b3-6bf8-45d4-98da-5eb1c6e7e6ec
# ╠═b11f3588-9904-46ef-8f35-a7c1d623e020
# ╠═d42a0cd3-cc8e-4a24-8887-7100f3927961
# ╠═389fefd6-60fd-4dd8-ba77-29e87b4ed846
# ╠═322e01f3-32d2-4408-91ef-ab7777f0f925
# ╠═fd74ebb7-625e-4396-a91a-37045a283a08
# ╠═d50cc41b-9465-4481-86d9-c8cf221e1978
# ╠═9ab4431d-e64f-44c9-a641-f16a524f3c19
# ╠═9c594ea7-8670-4235-8bde-cb2d670fe2c3
# ╠═933c2e0a-db49-46fd-b3ae-a40be190022d
# ╠═f3b219b4-7cc7-48a0-9675-91de7f40bb6d
# ╠═5d005025-8d6e-4d06-8b38-8aedd6189dfe
# ╠═4f854674-f062-4cb2-a026-f9bac6863990
# ╠═78a7391b-511a-45b8-b4e1-52b0bb675058
# ╠═e7ae2044-a31e-4a98-b496-3b1b3be056f0
# ╠═c51c6f13-b5ed-4671-9f08-4adb4f7b38b7
# ╠═5fc7e171-0fce-40f7-bc22-2aa02bb3694f
# ╠═ffaddffc-bbd0-4a08-9422-71bf872eb750
# ╠═50c692f7-93df-4edc-af0b-bb901a1bc907
# ╠═abdaea6d-504f-4708-8185-c63dcf90bad9
# ╠═7bfc460d-4ac2-4d02-aa3f-d93447fa71a6
# ╠═92db291e-afb7-490d-a599-d4bcb686dccb
# ╠═20b94af8-aaad-4832-8643-5d5d1f00e288
# ╠═3f13f074-e6b3-4b77-8080-bf8861b1e018
# ╠═cc48ec77-c92c-40d9-8569-6550f6dde552
# ╠═30a60b7d-630a-4fa3-ac57-e873655ad754
# ╠═183d7d89-b34d-4be8-9b7f-be039a6ce545
# ╟─1160f544-3879-4892-a47f-fc88d3620b8a
# ╠═5917e45e-50bb-441f-9d12-6901025a81f3
# ╠═66761f62-bfe3-4c5e-8630-17a8b63721e7
# ╠═c29da90b-32b5-445b-9de4-03cba2349397
# ╠═a1109c5b-ea82-40ac-8194-c7dcc4a069ea
# ╠═523186ec-f7c8-4968-afd5-167ca5f06254
# ╠═16034a2d-574e-422a-978b-0add93e3b128
# ╠═3d12e962-3f0a-447a-8af4-81ab0e7a34c3
# ╠═ad8ab6ab-5a63-49f4-a60e-7b6130abfbad
# ╠═8544dd0c-556f-45e2-b432-3258b257ad99
# ╠═5fa49bc9-a34d-4ae9-9cd3-4051281bc68f
# ╠═9002ba29-45ef-46dc-8632-e993c9ec2dc8
# ╠═0459f053-b627-427e-a698-a49ea8e68b3d
# ╟─e9f641fe-4ec6-495b-9634-a9be98395a4d
# ╠═27abef4c-62cf-4c82-a76f-38b735569328
# ╠═ef3ad471-392d-435a-950b-1b34b2c3b469
# ╠═cc8c0231-d568-4ac8-a090-428fdea08696
# ╠═0a029aaa-c23d-4126-bd21-323d978e1548
# ╠═72c8e61d-03af-4bac-ae60-0f0d0ff19492
# ╠═b14b8b05-2d08-4e6b-bb85-dcd01038518f
# ╠═9c4c10a1-6161-4c2f-ba10-05509cabeb49
# ╠═05eba4d6-f566-4e27-9bed-72a7520dfef0
# ╠═5772d3d1-0ff8-4a6a-a474-b96bc7d67347
# ╠═469acf5a-b7bf-4e53-bbef-6d57ab2988e9
# ╠═a2a77ae2-3a6b-46b2-bc5d-34a0e82fcef9
# ╠═0b33dfaa-055a-44b2-a463-1233047d624b
# ╠═fba7941e-8175-47df-8801-7d1f91211376
# ╠═7c8f5d6e-d53f-49a2-8b3a-239117c7c533
# ╠═58346276-5061-4d9a-8781-c1aba26d2c42
# ╠═19be0190-af5d-479b-b859-8f1c1e46bba1
# ╠═c070eb97-b7d8-45e8-ae95-0219da2cd3e6
# ╠═7c64d26a-bddd-4159-98d2-3a0f8d2125f5
