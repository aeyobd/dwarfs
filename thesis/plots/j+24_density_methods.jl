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

# ╔═╡ 5fc7e171-0fce-40f7-bc22-2aa02bb3694f
using CSV, DataFrames

# ╔═╡ 53a18bc2-3bbd-4fe3-8a77-d9c49d2ddea9
include("paper_style.jl")

# ╔═╡ 0ac92e62-f8fd-4e17-8b2f-feecdab75497
import TOML

# ╔═╡ ad8414dd-9428-47db-b08d-65940091413e
CairoMakie.activate!(type=:png)

# ╔═╡ c51c6f13-b5ed-4671-9f08-4adb4f7b38b7
import Statistics: median

# ╔═╡ 322e01f3-32d2-4408-91ef-ab7777f0f925
module Utils
	include("gaia_utils.jl")
end

# ╔═╡ 476200a0-8d64-4237-8b36-490cf2172e50
md"""
# Plot utils
"""

# ╔═╡ 54e449b3-6bf8-45d4-98da-5eb1c6e7e6ec
lw = theme(:linewidth)[]/2

# ╔═╡ b11f3588-9904-46ef-8f35-a7c1d623e020
log_r_ell_label = L"$\log\,R_\textrm{ell}$\,/\,arcmin"

# ╔═╡ d42a0cd3-cc8e-4a24-8887-7100f3927961
log_Σ_label = L"$\log\,\Sigma$\,/\,stars\ arcmin$^{-2}$"

# ╔═╡ 2e1e1605-a4f9-45f7-ab1b-56417d94237e
md"""
# Data inputs
"""

# ╔═╡ 389fefd6-60fd-4dd8-ba77-29e87b4ed846
obs_props_scl = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml")

# ╔═╡ fd74ebb7-625e-4396-a91a-37045a283a08
obs_props_umi = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/observed_properties.toml")

# ╔═╡ d50cc41b-9465-4481-86d9-c8cf221e1978
obs_props_fornax = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/fornax/observed_properties.toml")

# ╔═╡ 9ab4431d-e64f-44c9-a641-f16a524f3c19
α = LilGuys.R_h(LilGuys.Exp2D())

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
		y_R_b=nothing, R_h=nothing,
							y_min=-5
	)
	
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
	x_l, x_h, y_l, y_h = get_extrema(profiles)

	if R_h !== nothing
		x = log10(R_h)
		y = y_l
		annotation!(ax, 0, 36, x, y_min, color=:black, linewidth=lw, text=L"R_h")
	end


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

# ╔═╡ ef3ad471-392d-435a-950b-1b34b2c3b469
grey = Utils.SELECTION_COLORS[1]

# ╔═╡ 9c594ea7-8670-4235-8bde-cb2d670fe2c3
function get_R_h(obs_props)	
	R_h = obs_props["R_h"]
end

# ╔═╡ 9c4c10a1-6161-4c2f-ba10-05509cabeb49
R_h_scl = get_R_h(obs_props_scl)

# ╔═╡ 5772d3d1-0ff8-4a6a-a474-b96bc7d67347
R_h_umi = get_R_h(obs_props_umi)

# ╔═╡ a2a77ae2-3a6b-46b2-bc5d-34a0e82fcef9
R_h_fornax = get_R_h(obs_props_fornax)

# ╔═╡ fba7941e-8175-47df-8801-7d1f91211376
R_h_munoz1 = 0.49  # arcmin, munoz+2012

# ╔═╡ 933c2e0a-db49-46fd-b3ae-a40be190022d
R_munoz_1 = 36.58

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

# ╔═╡ 19be0190-af5d-479b-b859-8f1c1e46bba1
@savefig "scl_umi_fnx_density_methods" let
	fig = Figure(size=(5, 6) .*72)
	ax = Axis(fig[1,1],
		#xlabel = log_r_ell_label,
		#ylabel = log_Σ_label,
	)

	y_min = -3
	
	hlines!(get_background("sculptor")[1], color=grey, linewidth=1)
	compare_densities!(ax, profiles_scl, styles=plot_kwargs, R_h=R_h_scl, y_R_b=-2, y_min=y_min)

	text!(0.8, 0.8, text="Sculptor", space=:relative)

	hidexdecorations!(ticks=false, minorticks=false)

	
	# Ursa Minor
	ax_umi = Axis(fig[2,1],
		#xlabel = log_r_ell_label,
		#ylabel = log_Σ_label,
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
		#ylabel = log_Σ_label,
	)
	
	hlines!(get_background("fornax")[1], color=grey, linewidth=1)

	compare_densities!(ax_fnx, profiles_fornax, styles=plot_kwargs, 	
		 R_h=R_h_fornax, y_min=y_min
	)
	
	text!(0.8, 0.8, text="Fornax", space=:relative)


	linkxaxes!(ax, ax_umi, ax_fnx)
	rowgap!(fig.layout, 0.)

	Label(fig[:, 0], log_Σ_label, rotation=π/2)
	fig

end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═0ac92e62-f8fd-4e17-8b2f-feecdab75497
# ╠═ad8414dd-9428-47db-b08d-65940091413e
# ╠═53a18bc2-3bbd-4fe3-8a77-d9c49d2ddea9
# ╠═e7ae2044-a31e-4a98-b496-3b1b3be056f0
# ╠═c51c6f13-b5ed-4671-9f08-4adb4f7b38b7
# ╠═5fc7e171-0fce-40f7-bc22-2aa02bb3694f
# ╠═322e01f3-32d2-4408-91ef-ab7777f0f925
# ╠═476200a0-8d64-4237-8b36-490cf2172e50
# ╠═54e449b3-6bf8-45d4-98da-5eb1c6e7e6ec
# ╠═b11f3588-9904-46ef-8f35-a7c1d623e020
# ╠═d42a0cd3-cc8e-4a24-8887-7100f3927961
# ╟─2e1e1605-a4f9-45f7-ab1b-56417d94237e
# ╠═389fefd6-60fd-4dd8-ba77-29e87b4ed846
# ╠═fd74ebb7-625e-4396-a91a-37045a283a08
# ╠═d50cc41b-9465-4481-86d9-c8cf221e1978
# ╠═9ab4431d-e64f-44c9-a641-f16a524f3c19
# ╠═5d005025-8d6e-4d06-8b38-8aedd6189dfe
# ╠═4f854674-f062-4cb2-a026-f9bac6863990
# ╠═78a7391b-511a-45b8-b4e1-52b0bb675058
# ╠═7bfc460d-4ac2-4d02-aa3f-d93447fa71a6
# ╠═20b94af8-aaad-4832-8643-5d5d1f00e288
# ╟─e9f641fe-4ec6-495b-9634-a9be98395a4d
# ╠═27abef4c-62cf-4c82-a76f-38b735569328
# ╠═c070eb97-b7d8-45e8-ae95-0219da2cd3e6
# ╠═7c64d26a-bddd-4159-98d2-3a0f8d2125f5
# ╠═ef3ad471-392d-435a-950b-1b34b2c3b469
# ╠═9c594ea7-8670-4235-8bde-cb2d670fe2c3
# ╠═9c4c10a1-6161-4c2f-ba10-05509cabeb49
# ╠═5772d3d1-0ff8-4a6a-a474-b96bc7d67347
# ╠═a2a77ae2-3a6b-46b2-bc5d-34a0e82fcef9
# ╠═fba7941e-8175-47df-8801-7d1f91211376
# ╠═933c2e0a-db49-46fd-b3ae-a40be190022d
# ╠═cc8c0231-d568-4ac8-a090-428fdea08696
# ╠═19be0190-af5d-479b-b859-8f1c1e46bba1
