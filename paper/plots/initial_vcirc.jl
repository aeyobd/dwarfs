### A Pluto.jl notebook ###
# v0.20.24

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

# ╔═╡ b74a4535-1e12-422e-a524-135db4d8a9f7
import TOML

# ╔═╡ 20c338e3-3d10-41f4-b8ee-d0eda4e755bd
CairoMakie.activate!(type=:png)

# ╔═╡ a4c87a6e-976d-4af6-868e-09cb85e3d424
module Utils
	include("utils.jl")
end

# ╔═╡ cdb6446b-b129-4144-a90f-0a0d1c20e734
module ModelUtils
	include("model_utils.jl")
end

# ╔═╡ 43dad299-d8d2-4146-8320-c90a62f3a3f0
md"""
# Data loading
"""

# ╔═╡ d41f4781-5063-48c0-aade-fcd3980e19e7
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ 50839b67-9514-47b9-add2-0c84d05f12da
function get_obs_props(galaxyname)
	return TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))
end

# ╔═╡ 5ac62ad0-6413-4d5f-a48f-080590ac862e
function get_M_star(galaxyname)
	obs_props = get_obs_props(galaxyname)

	return LilGuys.mag_to_L(obs_props["Mv"]) * obs_props["M_L_s"] / M2MSUN
end

# ╔═╡ d66619f4-de41-45f4-9be1-4b938b9a959c
function load_snaps(galaxyname, modelname, starsname)
	haloname = joinpath(modelname, "..")
	stars_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, "stars")

	df_probs = LilGuys.read_hdf5_table(joinpath(stars_dir, starsname, "probabilities_stars.hdf5"))

	snap_i = Snapshot(joinpath(stars_dir, "iso_paint.hdf5"))
	LilGuys.add_stars!(snap_i, df_probs.probability)

	snap_f = Snapshot(joinpath(stars_dir, "iso_final.hdf5"))
	LilGuys.add_stars!(snap_f, df_probs.probability)

	return snap_i, snap_f
end

# ╔═╡ 890fa950-ca5e-4bec-9ef2-a2e5b3909909
function load_stars(galaxyname, modelname, starsname)

	snap_i, snap_f = load_snaps(galaxyname, modelname, starsname)
	Mstar = get_M_star(galaxyname)
	@info "stellar mass $Mstar"

	prof_s_i = MassProfile(snap_i, snap_i.weights * Mstar)
	prof_i = MassProfile(snap_i)
	prof_s_f = MassProfile(snap_f, snap_f.weights * Mstar)
	prof_f = MassProfile(snap_f)


	return (;dm_i = prof_i, dm_f = prof_f, stars_i=prof_s_i, stars_f=prof_s_f)
end

# ╔═╡ 8ac764c8-ca9e-400b-bde4-f3e9b1ebf6b9
scl_profs = load_stars(modelnames["scl_smallperi"]...)

# ╔═╡ edc6003e-e931-45a0-82c4-4aa394a3f060
umi_profs = load_stars(modelnames["umi_smallperi"]...)

# ╔═╡ dd6e4582-d610-4669-8fdf-b458bacb86b3
md"""
# Plot utils
"""

# ╔═╡ ec7d88d3-7503-4c18-9088-04319c6f2e98
smallfontsize=0.8 * theme(:fontsize)[]

# ╔═╡ 50751c0d-1b4e-4d56-bb6a-6ce6b2cfef3c
function plot_prof(prof; text=nothing, color=nothing, kwargs...)
	x = radii(prof)
	y = middle.(LilGuys.circular_velocity(prof) * V2KMS)
	lines!(x, y; joinstyle=:bevel, color=color, kwargs...)

	if !isnothing(text)
		Utils.text_along_line_log!(x, y, 0.05, text=text, color=color, align=(:left, :bottom), fontsize=smallfontsize)
	end

end

# ╔═╡ 7034614b-eb03-4333-bf69-c3fb64a5c0fc
function plot_i_f(profs; kwargs...)
	plot_prof(profs.dm_f;text="dark matter", kwargs...)
end

# ╔═╡ ee01807c-785e-43b8-a8bd-b093aacbf281
function plot_i_f_stars(profs; kwargs...)
	plot_prof(profs.stars_f; text="stars", kwargs...)
end

# ╔═╡ fb3dd9fe-2f1b-45fd-9d2d-3996859645b3
function vcirc_axis(gs)
	ax = Axis(gs, 
		xlabel = "radius / kpc",
		ylabel = L"circular velocity / km\,s$^{-1}$",
			  xscale = log10,
			  yscale = log10,
			  xticks = [0.1, 1, 10],
			  yticks = [1; 10:10:50],
			  yminorticks=[0.1:0.1:1; 1:10; 10:5:20; 20:5:50],
			  limits=(0.03, 100, 1, 40)
	)
end

# ╔═╡ 961f2851-8ebb-49ac-9459-5e6f543f80f4
function calc_r_J(halo, ρ_host)
	return LilGuys.find_zero(r -> 3ρ_host - LilGuys.mean_density(halo, r), 1)
end

# ╔═╡ 6c106a00-445b-4e0f-9168-66a9aeb35767
smalllinewidth=theme(:linewidth)[]/2

# ╔═╡ ebf4a089-0459-4241-b05d-69e1c514a4a5
function plot_sigma_v!(galaxy)
	obs_props = get_obs_props(galaxy) |> LilGuys.collapse_errors

	R_h = [LilGuys.arcmin2kpc(obs_props["R_h"], obs_props["distance"])]
	σv = [obs_props["sigma_v"]]
	
	errorscatter!(middle.(R_h), middle.(σv), 
				  xerror = error_interval.(R_h), 
				  yerror=error_interval.(σv), 
				  linewidth=smalllinewidth, color=COLORS[3])
	
	text!(middle.(R_h), middle.(σv), 
		  align=(:left, :center), 
		  offset=(theme(:fontsize)[]/2, 0), 
		  text=L"($\sigma_\textrm{v}$, $R_h$)",
		  color=COLORS[3], fontsize=smallfontsize)
end

# ╔═╡ 99675829-e4b2-4124-90c2-595e53290956
function plot_r_J(r_J)
	v0 = 10

	vlines!(r_J, color=:grey, linewidth=smalllinewidth,  linestyle=:dot)
	text!(r_J, v0, text=L"R_J", fontsize=smallfontsize, align=(:left, :center), color=:grey,)
end

# ╔═╡ e78a3af9-5e50-41e3-ab0b-4f8ba4a349f9
r_j_scl = ModelUtils.get_r_j(modelnames["scl_smallperi"][1:2]..., key="r_J_kpc")

# ╔═╡ d0675035-aa7f-44a1-9866-69c5567496fe
r_j_umi = ModelUtils.get_r_j(modelnames["umi_smallperi"][1:2]...,  key="r_J_kpc")

# ╔═╡ 7b2c0281-dba8-4e17-9e11-ef8acdc2eba6
@savefig "initial_velocity_nosigma" let
	fig = Figure(size=(10/3, 2) .*72)

	ax = vcirc_axis(fig[1,1])
	ax.title="Sculptor"
	
	plot_i_f(scl_profs, color=COLORS[5], label="Scl DM")
	plot_i_f_stars(scl_profs, 
				   color=COLORS[2], label="Scl stars",  linewidth=smalllinewidth)

	plot_r_J(r_j_scl)
	plot_sigma_v!("sculptor")


	xlims!(0.03, 30)
	ylims!(0.5, 50)

	
	ax_umi = vcirc_axis(fig[1,2])
	ax_umi.title = "Ursa Minor"
	
	plot_i_f(umi_profs, color=COLORS[5], label="UMi DM")
	plot_i_f_stars(umi_profs,
				   color=COLORS[2], label="UMi stars", linewidth=smalllinewidth)

	plot_r_J(r_j_umi)

	plot_sigma_v!("ursa_minor")
	hideydecorations!(ticks=false, minorticks=false)


	xlims!(0.03, 30)
	ylims!(0.5, 50)


	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═b74a4535-1e12-422e-a524-135db4d8a9f7
# ╠═20c338e3-3d10-41f4-b8ee-d0eda4e755bd
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═a4c87a6e-976d-4af6-868e-09cb85e3d424
# ╠═cdb6446b-b129-4144-a90f-0a0d1c20e734
# ╟─43dad299-d8d2-4146-8320-c90a62f3a3f0
# ╠═d41f4781-5063-48c0-aade-fcd3980e19e7
# ╠═50839b67-9514-47b9-add2-0c84d05f12da
# ╠═5ac62ad0-6413-4d5f-a48f-080590ac862e
# ╠═d66619f4-de41-45f4-9be1-4b938b9a959c
# ╠═890fa950-ca5e-4bec-9ef2-a2e5b3909909
# ╠═8ac764c8-ca9e-400b-bde4-f3e9b1ebf6b9
# ╠═edc6003e-e931-45a0-82c4-4aa394a3f060
# ╟─dd6e4582-d610-4669-8fdf-b458bacb86b3
# ╠═ec7d88d3-7503-4c18-9088-04319c6f2e98
# ╠═7034614b-eb03-4333-bf69-c3fb64a5c0fc
# ╠═ee01807c-785e-43b8-a8bd-b093aacbf281
# ╠═50751c0d-1b4e-4d56-bb6a-6ce6b2cfef3c
# ╠═ebf4a089-0459-4241-b05d-69e1c514a4a5
# ╠═fb3dd9fe-2f1b-45fd-9d2d-3996859645b3
# ╠═961f2851-8ebb-49ac-9459-5e6f543f80f4
# ╠═99675829-e4b2-4124-90c2-595e53290956
# ╠═6c106a00-445b-4e0f-9168-66a9aeb35767
# ╠═e78a3af9-5e50-41e3-ab0b-4f8ba4a349f9
# ╠═d0675035-aa7f-44a1-9866-69c5567496fe
# ╠═7b2c0281-dba8-4e17-9e11-ef8acdc2eba6
