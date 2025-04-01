### A Pluto.jl notebook ###
# v0.20.5

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

# ╔═╡ 589d79d1-da7d-428b-abb2-0ae68ad49a5d
import PythonCall

# ╔═╡ 0ac92e62-f8fd-4e17-8b2f-feecdab75497
import TOML

# ╔═╡ ad8414dd-9428-47db-b08d-65940091413e
CairoMakie.activate!(type=:png)

# ╔═╡ b11f3588-9904-46ef-8f35-a7c1d623e020
log_r_ell_label = L"\log\,R\,/\,R_h"

# ╔═╡ d42a0cd3-cc8e-4a24-8887-7100f3927961
log_Σ_label = L"$\log\,\Sigma$\,/\,stars\ ${R_h}^{-2}$"

# ╔═╡ 389fefd6-60fd-4dd8-ba77-29e87b4ed846
obs_props_scl = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml")

# ╔═╡ fd74ebb7-625e-4396-a91a-37045a283a08
obs_props_umi = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/observed_properties.toml")

# ╔═╡ d50cc41b-9465-4481-86d9-c8cf221e1978
obs_props_fornax = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/fornax/observed_properties.toml")

# ╔═╡ 5d005025-8d6e-4d06-8b38-8aedd6189dfe
function load_profile(galaxyname, algname="jax")
	filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, 
		"density_profiles/$(algname)_profile.toml")
	obs_props = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/$galaxyname/observed_properties.toml")
	
    prof = LilGuys.StellarDensityProfile(filename)
	prof = LilGuys.scale(prof, 1/obs_props["r_h"], 1)

	prof
end

# ╔═╡ 78a7391b-511a-45b8-b4e1-52b0bb675058
update_theme!(
	ErrorScatter=(; 
		cycle=Cycle([:color, :marker], covary=true)
		))

# ╔═╡ c51c6f13-b5ed-4671-9f08-4adb4f7b38b7
import Statistics: median

# ╔═╡ 5917e45e-50bb-441f-9d12-6901025a81f3
begin 
	profiles = OrderedDict(
		"fiducial" => load_profile("sculptor", "jax_2c"),
		"circ" => load_profile("sculptor", "jax_circ"),
		"simple" => load_profile("sculptor", "simple"),
		"bright" => load_profile("sculptor", "jax_bright"),
		#"DELVE" => load_profile("delve"),
		"DELVE" => load_profile("sculptor", "delve_rgb_sub"),
	)

	profiles["bright"] = LilGuys.scale(profiles["bright"], 1, 2)
	#profiles["DELVE"] = LilGuys.scale(profiles["DELVE"], 1, 7000/1197)
	profiles["DELVE"] = LilGuys.scale(profiles["DELVE"], 1, 0.6)

	profiles
end

# ╔═╡ ffaddffc-bbd0-4a08-9422-71bf872eb750
df_peris = CSV.read(ENV["DWARFS_ROOT"] * "/analysis/all_galaxies/vasiliev+21/properties_w_orbits.csv", DataFrame)


# ╔═╡ 4be65b34-fbb9-46fe-9c8c-462abb5d0571
typeof((;))

# ╔═╡ 30a60b7d-630a-4fa3-ac57-e873655ad754
function compare_densities(profiles, r_h; styles=nothing)
	
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
	)

	for (name, prof) in profiles

		if styles !== nothing
			kwargs = styles[name]
		else
			kwargs = (;)
		end
		
		errorscatter!(prof.log_R, prof.log_Sigma; 
			yerror=LilGuys.error_interval.(prof.log_Sigma), 
			label=name, kwargs...)
	end

	#arrows!([log10(r_h)], [-2], [0], [-0.5], arrowsize=6)
	#text!([log10(r_h)], [-2], text=L"R_h")

	axislegend(position=:lb)
	ax2 = Axis(fig[2,1],
		xlabel = log_r_ell_label,
		ylabel = L"\log\Sigma / \Sigma_\textrm{%$(collect(keys(profiles))[begin])}",
		limits=(nothing, nothing, -1, 1)
	)

	hlines!(0, color=:black)

	prof_ref = profiles[collect(keys(profiles))[begin]]
	for (key, prof) in profiles
		ym = LilGuys.lerp(prof_ref.log_R, LilGuys.middle.(prof_ref.log_Sigma)).(prof.log_R)
		filt = prof.log_R .< prof_ref.log_R[end]
		filt .&= prof.log_R .> prof_ref.log_R[1]
		filt .&= isfinite.(ym)
		filt .&= isfinite.(prof.log_Sigma)

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

	
	
	# Legend(fig[1,1], fig.content[1], halign=0.1, valign=0.1, tellwidth=false)
	fig
end

# ╔═╡ d6509ab8-9923-443b-8b96-50e4f59498ab
@savefig "scl_density_methods_extra" compare_densities(profiles, obs_props_scl["r_h"])

# ╔═╡ c9e23e51-c2de-427f-8429-d620796e1baa
load_profile("sculptor", "best_eqw_sub") #|> LilGuys.filter_empty_bins

# ╔═╡ 27abef4c-62cf-4c82-a76f-38b735569328
begin 
	profiles_scl = OrderedDict(
		"probable members" => load_profile("sculptor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
		"CMD + PM" => load_profile("sculptor", "jax_LLR_0_eqw") |> LilGuys.filter_empty_bins,
		"all" => load_profile("sculptor", "best_eqw") |> LilGuys.filter_empty_bins,
		"BG subtracted" => load_profile("sculptor", "best_eqw_sub"),
	)

	# profiles_umi["bright"] = LilGuys.scale(profiles_umi["bright"], 1, 2)
	# #profiles["DELVE"] = LilGuys.scale(profiles["DELVE"], 1, 7000/1197)
	# profiles_umi["DELVE"] = LilGuys.scale(profiles_umi["DELVE"], 1, 0.6)

	profiles_scl
end

# ╔═╡ cc8c0231-d568-4ac8-a090-428fdea08696
plot_kwargs = Dict(
	"probable members" => (;
		color=COLORS[2]
	),
	"CMD + PM" => (;
		color = COLORS[1],
		marker = :utriangle
	),
	"all" => (;
		color = RGBf(0.7, 0.7, 0.7),
		marker = :rect
	),
	"BG subtracted" => (;
		color = COLORS[5],
		marker = :diamond
	)
)

# ╔═╡ dbeae6d9-1dd7-403b-80ac-82b20d047c3b
@savefig "scl_density_methods" compare_densities(profiles_scl, obs_props_scl["r_h"], styles=plot_kwargs)

# ╔═╡ c070eb97-b7d8-45e8-ae95-0219da2cd3e6
begin 
	profiles_umi = OrderedDict(
		"probable members" => load_profile("ursa_minor", "jax_2c_eqw") |> LilGuys.filter_empty_bins,
		"CMD + PM" => load_profile("ursa_minor", "jax_LLR_0_eqw") |> LilGuys.filter_empty_bins,
		"all" => load_profile("ursa_minor", "best_eqw") |> LilGuys.filter_empty_bins,
		"BG subtracted" => load_profile("ursa_minor", "best_eqw_sub"),
	)

	# profiles_umi["bright"] = LilGuys.scale(profiles_umi["bright"], 1, 2)
	# #profiles["DELVE"] = LilGuys.scale(profiles["DELVE"], 1, 7000/1197)
	# profiles_umi["DELVE"] = LilGuys.scale(profiles_umi["DELVE"], 1, 0.6)

	profiles_umi
end

# ╔═╡ 18b9cf25-527e-4bac-91ca-5d1fa69293c8
@savefig "umi_density_methods" compare_densities(profiles_umi, obs_props_umi["r_h"], styles=plot_kwargs)

# ╔═╡ 7c64d26a-bddd-4159-98d2-3a0f8d2125f5
begin 
	profiles_fornax = OrderedDict(
		"probable members" => load_profile("fornax", "jax_eqw") |> LilGuys.filter_empty_bins,
		"CMD + PM" => load_profile("fornax", "jax_LLR_0_eqw") |> LilGuys.filter_empty_bins,
		"all" => load_profile("fornax", "best_eqw") |> LilGuys.filter_empty_bins,
		"BG subtracted" => load_profile("fornax", "best_eqw_sub") ,
	)

	# profiles["bright"] = LilGuys.scale(profiles["bright"], 1, 2)
	# #profiles["DELVE"] = LilGuys.scale(profiles["DELVE"], 1, 7000/1197)
	# profiles["DELVE"] = LilGuys.scale(profiles["DELVE"], 1, 0.6)

	profiles
end

# ╔═╡ 188ce016-6011-45c3-9823-eafdc94f9617
@savefig "fornax_density_methods" compare_densities(profiles_fornax, obs_props_fornax["r_h"], styles=plot_kwargs)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═589d79d1-da7d-428b-abb2-0ae68ad49a5d
# ╠═0ac92e62-f8fd-4e17-8b2f-feecdab75497
# ╠═ad8414dd-9428-47db-b08d-65940091413e
# ╠═53a18bc2-3bbd-4fe3-8a77-d9c49d2ddea9
# ╠═b11f3588-9904-46ef-8f35-a7c1d623e020
# ╠═d42a0cd3-cc8e-4a24-8887-7100f3927961
# ╠═389fefd6-60fd-4dd8-ba77-29e87b4ed846
# ╠═fd74ebb7-625e-4396-a91a-37045a283a08
# ╠═d50cc41b-9465-4481-86d9-c8cf221e1978
# ╠═5d005025-8d6e-4d06-8b38-8aedd6189dfe
# ╠═78a7391b-511a-45b8-b4e1-52b0bb675058
# ╠═e7ae2044-a31e-4a98-b496-3b1b3be056f0
# ╠═c51c6f13-b5ed-4671-9f08-4adb4f7b38b7
# ╠═5917e45e-50bb-441f-9d12-6901025a81f3
# ╠═5fc7e171-0fce-40f7-bc22-2aa02bb3694f
# ╠═ffaddffc-bbd0-4a08-9422-71bf872eb750
# ╠═4be65b34-fbb9-46fe-9c8c-462abb5d0571
# ╠═30a60b7d-630a-4fa3-ac57-e873655ad754
# ╠═d6509ab8-9923-443b-8b96-50e4f59498ab
# ╠═c9e23e51-c2de-427f-8429-d620796e1baa
# ╠═27abef4c-62cf-4c82-a76f-38b735569328
# ╠═cc8c0231-d568-4ac8-a090-428fdea08696
# ╠═dbeae6d9-1dd7-403b-80ac-82b20d047c3b
# ╠═c070eb97-b7d8-45e8-ae95-0219da2cd3e6
# ╠═18b9cf25-527e-4bac-91ca-5d1fa69293c8
# ╠═7c64d26a-bddd-4159-98d2-3a0f8d2125f5
# ╠═188ce016-6011-45c3-9823-eafdc94f9617
