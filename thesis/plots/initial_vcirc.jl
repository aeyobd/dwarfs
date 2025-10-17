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

# ╔═╡ df4f20bf-97d9-408a-859e-e4070edcd0ef
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 11cc55a3-d166-4bf3-b147-1c789d690f90
import Agama

# ╔═╡ b74a4535-1e12-422e-a524-135db4d8a9f7
import TOML

# ╔═╡ d41f4781-5063-48c0-aade-fcd3980e19e7
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ 20c338e3-3d10-41f4-b8ee-d0eda4e755bd
CairoMakie.activate!(type=:png)

# ╔═╡ 50839b67-9514-47b9-add2-0c84d05f12da
module Utils
	include("utils.jl")
end

# ╔═╡ 43dad299-d8d2-4146-8320-c90a62f3a3f0
md"""
# Data loading
"""

# ╔═╡ 0dbc4770-fd1d-4a57-9825-3576900dd7a3
function load_vcircs(galaxyname, haloname)
	modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, "stars")

	snap_i = Snapshot(joinpath(modeldir, "iso_initial.hdf5"))
	snap_f = Snapshot(joinpath(modeldir, "iso_final.hdf5"))

	prof_i = MassProfile(snap_i)
	prof_f = MassProfile(snap_f)

	return prof_i, prof_f
end

# ╔═╡ 890fa950-ca5e-4bec-9ef2-a2e5b3909909
function load_stars(galaxyname, modelname, starsname)
	haloname = joinpath(modelname, "..")
	stars_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, "stars")

	df_probs = LilGuys.read_hdf5_table(joinpath(stars_dir, starsname, "probabilities_stars.hdf5"))

	snap_i = Snapshot(joinpath(stars_dir, "iso_paint.hdf5"))
	LilGuys.add_stars!(snap_i, df_probs.probability)

	snap_f = Snapshot(joinpath(stars_dir, "iso_final.hdf5"))
	LilGuys.add_stars!(snap_f, df_probs.probability)

	Mstar = Utils.get_M_star(galaxyname)
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

# ╔═╡ 50751c0d-1b4e-4d56-bb6a-6ce6b2cfef3c
function plot_prof(prof; text=nothing, color=nothing, kwargs...)
	x = radii(prof)
	y = middle.(LilGuys.circular_velocity(prof) * V2KMS)
	lines!(x, y; color=color, kwargs...)

	if !isnothing(text)
		Utils.text_along_line_log!(x, y, 0.05, text=text, color=color, align=(:left, :bottom))
	end

end

# ╔═╡ 7034614b-eb03-4333-bf69-c3fb64a5c0fc
function plot_i_f(profs; kwargs...)
	plot_prof(profs.dm_i, text="dark matter", linestyle=:dot; kwargs...)
	plot_prof(profs.dm_f; kwargs...)
end

# ╔═╡ ee01807c-785e-43b8-a8bd-b093aacbf281
function plot_i_f_stars(profs; kwargs...)
	plot_prof(profs.stars_i, text="stars", linestyle=:dot; kwargs...)
	plot_prof(profs.stars_f; kwargs...)
end

# ╔═╡ c1d64f7b-b7c3-4968-b564-0dba5c119e09
α_3d = LilGuys.r_h(LilGuys.Exp2D()) / LilGuys.R_h(LilGuys.Exp2D()) 

# ╔═╡ fb3dd9fe-2f1b-45fd-9d2d-3996859645b3
function vcirc_axis(gs)
	ax = Axis(gs, 
		xlabel = "radius / kpc",
		ylabel = L"circular velocity / km\,s$^{-1}$",
			  xscale = log10,
			  yscale = log10,
			  xticks = [0.1, 1, 10],
			  yticks = [1; 10:10:40],
			  yminorticks=[1:10; 10:5:20; 20:5:40],
			  limits=(0.03, 100, 1, 40)
	)
end

# ╔═╡ 6c106a00-445b-4e0f-9168-66a9aeb35767
smalllinewidth=theme(:linewidth)[]/2

# ╔═╡ ebf4a089-0459-4241-b05d-69e1c514a4a5
function plot_sigma_v!(galaxy, R_s=0.10)
	obs_props = Utils.get_obs_props(galaxy) |> LilGuys.collapse_errors

	R_h = [LilGuys.arcmin2kpc(obs_props["R_h"], obs_props["distance"])]
	σv = [obs_props["sigma_v"]]
	
	errorscatter!(middle.(R_h), middle.(σv), 
				  xerror = error_interval.(R_h), 
				  yerror=error_interval.(σv), 
				  linewidth=smalllinewidth, markersize=1, color=COLORS[3])
	
	text!(middle.(R_h), middle.(σv), 
		  align=(:left, :center), 
		  offset=(theme(:fontsize)[]/2, 0), 
		  text=L"($\sigma_\textrm{v}$, $R_h$)",
		  color=COLORS[3])
end

# ╔═╡ 62407923-e1ba-445c-8726-2d72121a4ed2
function set_limits!(maxes)
	xlims!(0.03, maxes.r_circ_max * 10^1)
	ylims!(maxes.v_circ_max*V2KMS* 10^-1.9, maxes.v_circ_max*V2KMS * 10^0.1)

end

# ╔═╡ 1443385e-6d9b-4ddc-805c-2c12fb1c9067
scl_max = LilGuys.fit_v_r_circ_max(radii(scl_profs[1]), middle.(LilGuys.circular_velocity(scl_profs[1])))

# ╔═╡ 195582cf-83e3-44c8-9efe-b7d98f092c3c
umi_max = LilGuys.fit_v_r_circ_max(radii(umi_profs[1]), middle.(LilGuys.circular_velocity(umi_profs[1])))

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "initial_velocity" let
	fig = Figure(size=(6*72, 4*72))

	ax = vcirc_axis(fig[1,1])
	ax.title="Sculptor"
	
	plot_i_f(scl_profs, color=COLORS[1], label="Scl DM")
	plot_i_f_stars(scl_profs, 
				   color=COLORS[2], label="Scl stars",  linewidth=smalllinewidth)
	# plot_density_line!(get_mean_density("sculptor")...)
	# plot_r_J(NFW(r_circ_max=3.2, v_circ_max=31/V2KMS), get_mean_density("sculptor")[1])
	
	plot_sigma_v!("sculptor")
	set_limits!(scl_max)

	ax_umi = vcirc_axis(fig[1,2])
	ax_umi.title = "Ursa Minor"
	
	plot_i_f(umi_profs, color=COLORS[1], label="UMi DM")
	plot_i_f_stars(umi_profs,
				   color=COLORS[2], label="UMi stars", linewidth=smalllinewidth)

	# plot_density_line!(get_mean_density("ursa_minor")...)
	# plot_r_J(NFW(r_circ_max=4, v_circ_max=38/V2KMS), get_mean_density("ursa_minor")[1])

	plot_sigma_v!("ursa_minor")
	set_limits!(umi_max)

	hideydecorations!(ax_umi, ticks=false, minorticks=false)

	l1 = lines!([NaN], [NaN], linewidth=smalllinewidth, linestyle=:dot, label="initial")
	l2 = lines!([NaN], [NaN], linewidth=smalllinewidth, linestyle=:solid, label="final")
	linkyaxes!(ax, ax_umi)
	colgap!(fig.layout, 0)

	axislegend(ax, [l1, l2], ["initial conditions", "end of isolation"], position=:lb, backgroundcolor=(:white, 0.8))
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═11cc55a3-d166-4bf3-b147-1c789d690f90
# ╠═b74a4535-1e12-422e-a524-135db4d8a9f7
# ╠═df4f20bf-97d9-408a-859e-e4070edcd0ef
# ╠═d41f4781-5063-48c0-aade-fcd3980e19e7
# ╠═20c338e3-3d10-41f4-b8ee-d0eda4e755bd
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═50839b67-9514-47b9-add2-0c84d05f12da
# ╟─43dad299-d8d2-4146-8320-c90a62f3a3f0
# ╠═0dbc4770-fd1d-4a57-9825-3576900dd7a3
# ╠═890fa950-ca5e-4bec-9ef2-a2e5b3909909
# ╠═8ac764c8-ca9e-400b-bde4-f3e9b1ebf6b9
# ╠═edc6003e-e931-45a0-82c4-4aa394a3f060
# ╟─dd6e4582-d610-4669-8fdf-b458bacb86b3
# ╠═7034614b-eb03-4333-bf69-c3fb64a5c0fc
# ╠═ee01807c-785e-43b8-a8bd-b093aacbf281
# ╠═50751c0d-1b4e-4d56-bb6a-6ce6b2cfef3c
# ╠═c1d64f7b-b7c3-4968-b564-0dba5c119e09
# ╠═ebf4a089-0459-4241-b05d-69e1c514a4a5
# ╠═fb3dd9fe-2f1b-45fd-9d2d-3996859645b3
# ╠═6c106a00-445b-4e0f-9168-66a9aeb35767
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═62407923-e1ba-445c-8726-2d72121a4ed2
# ╠═1443385e-6d9b-4ddc-805c-2c12fb1c9067
# ╠═195582cf-83e3-44c8-9efe-b7d98f092c3c
