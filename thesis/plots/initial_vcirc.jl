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

# ╔═╡ 20c338e3-3d10-41f4-b8ee-d0eda4e755bd
CairoMakie.activate!(type=:png)

# ╔═╡ 50839b67-9514-47b9-add2-0c84d05f12da
module Utils
	include("utils.jl")
end

# ╔═╡ 0dbc4770-fd1d-4a57-9825-3576900dd7a3
function load_vcircs(galaxyname, haloname)
	modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, "stars")

	snap_i = Snapshot(joinpath(modeldir, "iso_initial.hdf5"))
	snap_f = Snapshot(joinpath(modeldir, "iso_final.hdf5"))

	prof_i = MassProfile(snap_i)
	prof_f = MassProfile(snap_f)

	return prof_i, prof_f
end

# ╔═╡ 5a4853a9-1060-4be8-a872-c78297ed586f
function assign_weights!(snap, df)
	snap.weights = df.probability[snap.index]
end

# ╔═╡ b7572e12-85a0-4f80-b7aa-ab2db8a81a13
function get_M_star(galaxyname)
	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))
	M_star = LilGuys.mag_to_L(obs_props["Mv"]) * obs_props["M_L_s"]
	return M_star / M2MSUN
end

# ╔═╡ 890fa950-ca5e-4bec-9ef2-a2e5b3909909
function load_stars(galaxyname, haloname, starsname)
	stars_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, "stars")

	df_probs = LilGuys.read_hdf5_table(joinpath(stars_dir, starsname, "probabilities_stars.hdf5"))

	snap_i = Snapshot(joinpath(stars_dir, "iso_paint.hdf5"))
	assign_weights!(snap_i, df_probs)

	snap_f = Snapshot(joinpath(stars_dir, "iso_final.hdf5"))
	assign_weights!(snap_f, df_probs)

	Mstar = get_M_star(galaxyname)
	@info "stellar mass $Mstar"

	prof_s_i = MassProfile(snap_i, snap_i.weights * Mstar)
	prof_i = MassProfile(snap_i)
	prof_s_f = MassProfile(snap_f, snap_f.weights * Mstar)
	prof_f = MassProfile(snap_f)

	return (;dm_i = prof_i, dm_f = prof_f, stars_i=prof_s_i, stars_f=prof_s_f)
end

# ╔═╡ 8ac764c8-ca9e-400b-bde4-f3e9b1ebf6b9
scl_profs = load_stars("sculptor", "1e7_new_v31_r3.2", "exp2d_rs0.10")

# ╔═╡ 4258c055-e11a-409d-9b53-441c55cd815f
get_M_star("sculptor")

# ╔═╡ edc6003e-e931-45a0-82c4-4aa394a3f060
umi_profs = load_stars("ursa_minor", "1e7_new_v38_r4.0", "exp2d_rs0.10")

# ╔═╡ 1443385e-6d9b-4ddc-805c-2c12fb1c9067
scl_max = LilGuys.fit_v_r_circ_max(radii(scl_profs[1]), middle.(LilGuys.circular_velocity(scl_profs[1])))

# ╔═╡ 195582cf-83e3-44c8-9efe-b7d98f092c3c
umi_max = LilGuys.fit_v_r_circ_max(radii(umi_profs[1]), middle.(LilGuys.circular_velocity(umi_profs[1])))

# ╔═╡ ed95b670-3ba7-438a-8e4d-9fc6362d4d79
function text_along_line!(x, y, x_0; text, h=0.03, kwargs...)
	f = LilGuys.lerp(x, y)
	y_0 = f(x_0)
	dy = Utils.log_derivative(f, x_0, h=h)
	rf = Utils.rotation_factor(Makie.current_axis(), true)
	θ = @lift atan($rf * dy)

	text!(x_0, y_0, text=text, rotation=θ; kwargs...)
end

# ╔═╡ 7034614b-eb03-4333-bf69-c3fb64a5c0fc
function plot_i_f(profs; kwargs...)
	prof_i, prof_f = profs

	x_i, y_i = radii(prof_i), middle.(LilGuys.circular_velocity(prof_i)) * V2KMS
	lines!(x_i, y_i; linestyle=:dot, kwargs...)

	x, y = radii(prof_f), middle.(LilGuys.circular_velocity(prof_f)) * V2KMS

	lines!(x, y; kwargs...)


	text_along_line!(x_i, y_i, 0.05, text="dark matter", color=COLORS[1], align=(:left, :bottom), h=0.3)
	
end


# ╔═╡ ee01807c-785e-43b8-a8bd-b093aacbf281
function plot_i_f_stars(profs; kwargs...)
	x_i, y_i = radii(profs.stars_i), middle.(LilGuys.circular_velocity(profs.stars_f) * V2KMS)
	lines!(x_i, y_i, linestyle=:dot; kwargs...)

	x, y = radii(profs.stars_f), middle.(LilGuys.circular_velocity(profs.stars_f) * V2KMS)
	lines!(x, y; kwargs...)
	#lines!(radii(prof_f), LilGuys.circular_velocity(prof_f) * V2KMS, linestyle=:solid; kwargs...)

	text_along_line!(x_i, y_i, 0.05, text="stars", color=COLORS[2], align=(:left, :bottom))

end


# ╔═╡ 29637bda-3a48-4490-a2df-22ffc7fec59b
function get_obs_props(galaxyname)
	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))
end

# ╔═╡ 430f1b00-9746-4e56-8608-deef86f1d0bb
function load_peris(galaxyname)
	df = read_fits(joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, "EP2020", "orbital_properties.fits"))

	return df.pericentre
end

# ╔═╡ 0813fc5e-236d-4665-b4f2-f4a1790c9aba
import StatsBase: median, quantile

# ╔═╡ fbaaf028-c25c-4e95-900e-c6fbb43e05bb
function ρ_to_v(ρ, r)
	v_scale = sqrt(ρ * 4π/3)
	v = v_scale .* r
end

# ╔═╡ 6c106a00-445b-4e0f-9168-66a9aeb35767
smalllinewidth=theme(:linewidth)[]/2

# ╔═╡ cf449735-0acb-4e83-bb6a-7b0bf5f794c4
function plot_density_line!(ρ, ρ_el, ρ_ep)
	
	r = 10 .^ LinRange(-2, 2, 10)
	v = ρ_to_v(ρ, r)
	lines!(r, v*V2KMS, color=:black, linewidth=smalllinewidth)
	
	v_l = ρ_to_v(ρ + ρ_el, r)
	v_h = ρ_to_v(ρ + ρ_ep, r)
	band!(r, v_l*V2KMS, v_h*V2KMS, alpha=0.5, color=:black)

	text_along_line!(r, v*V2KMS, 2, text=L"\bar\rho_\textrm{MW,\ peri}", align=(:center, :bottom))
end

# ╔═╡ 8e716378-3a9b-4994-87f1-e77a164d9f57
pot = Agama.Potential(file = joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))

# ╔═╡ c47917c7-9eda-489c-bcf2-a602c01dcb24
function get_mean_density(galaxyname)
	peris = load_peris(galaxyname)
	ρs =  Agama.enclosed_mass(pot, peris) ./ (4π/3 * peris .^ 3)

	l, m, h = quantile(ρs, [0.16, 0.5, 0.84])
	return m, l-m, h-m
end

# ╔═╡ ba3943bd-e94c-47fa-9856-827f1f28fc40
get_mean_density("sculptor")

# ╔═╡ 3d0c926b-d2c6-4dcc-842e-f5be4f584fb8
load_peris("sculptor")

# ╔═╡ c1d64f7b-b7c3-4968-b564-0dba5c119e09
α_3d = LilGuys.r_h(LilGuys.Exp2D()) / LilGuys.R_h(LilGuys.Exp2D()) 

# ╔═╡ 4cb47cac-e3c6-41f3-ae8a-81f296414c52
α_3d * 0.24 / LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ d447b0cd-3b58-4c7d-9abe-d1ffc4c6ae47
s = LilGuys.Exp2D(R_s=0.10)

# ╔═╡ b04cbab4-abce-4f3d-8541-ca3f00a6950b
h = LilGuys.NFW(r_circ_max = 3.2, v_circ_max=31/V2KMS)

# ╔═╡ 9b44a224-5fcb-4d28-80ad-83fe38a7ff7b
LilGuys.σv_star_mean(h, s)

# ╔═╡ 8205bd52-8964-4054-9231-67a489b984f5
LilGuys.v_circ(h, LilGuys.r_h(s)) / sqrt(3)

# ╔═╡ ebf4a089-0459-4241-b05d-69e1c514a4a5
function plot_sigma_v!(galaxy, R_s=0.10)
	obs_props = get_obs_props(galaxy) |> LilGuys.collapse_errors

	R_h = [LilGuys.arcmin2kpc(obs_props["R_h"], obs_props["distance"])]
	@info R_h
	σv = [obs_props["sigma_v"]]
	
	errorscatter!(middle.(R_h), middle.(σv), xerror = error_interval.(R_h), yerror=error_interval.(σv), markersize=smalllinewidth, color=COLORS[3])

	shortname = Dict("sculptor"=>"Scl", "ursa_minor"=>"UMi")[galaxy]
	text!(middle.(R_h), middle.(σv), align=(:left, :center), offset=(theme(:fontsize)[]/2, 0), text=L"($\sigma_\textrm{v}$, $R_h$)", color=COLORS[3])
	

	# s = LilGuys.Exp2D(R_s=R_s)
	# @info "R_h = $(LilGuys.r_h(s))"

end

# ╔═╡ cf800bd7-2cd1-4ac1-870b-ba29056e7920
LilGuys.mean_density(NFW(r_circ_max=3.2, v_circ_max = 31/V2KMS), 0.1)

# ╔═╡ fb3dd9fe-2f1b-45fd-9d2d-3996859645b3
function vcirc_axis(gs)
	ax = Axis(gs, 
		xlabel = "radius / kpc",
		ylabel = L"circular velocity / km\,s$^{-1}$",
			  xscale = log10,
			  yscale = log10,
			  xticks = Makie.automatic,
			  yticks = [1; 10:10:40],
			  yminorticks=[1:10; 10:5:20; 20:5:40],
			  limits=(0.03, 100, 1, 40)
	)
end

# ╔═╡ 62407923-e1ba-445c-8726-2d72121a4ed2
function set_limits!(maxes)
	xlims!(0.03, maxes.r_circ_max * 10^1)
	ylims!(maxes.v_circ_max*V2KMS* 10^-1.9, maxes.v_circ_max*V2KMS * 10^0.1)

end

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "initial_velocity" let
	fig = Figure(size=(5*72, 5*72))

	ax = vcirc_axis(fig[1,1])
	ax.title="Sculptor"
	
	plot_i_f(scl_profs, color=COLORS[1], label="Scl DM")
	plot_i_f_stars(scl_profs, 
				   color=COLORS[2], label="Scl stars",  linewidth=smalllinewidth)
	plot_density_line!(get_mean_density("sculptor")...)
	
	plot_sigma_v!("sculptor")
	set_limits!(scl_max)

	ax_umi = vcirc_axis(fig[2,1])
	ax_umi.title = "Ursa Minor"
	
	plot_i_f(umi_profs, color=COLORS[1], label="UMi DM")
	plot_i_f_stars(umi_profs,
				   color=COLORS[2], label="UMi stars", linewidth=smalllinewidth)

	plot_density_line!(get_mean_density("ursa_minor")...)
	
	plot_sigma_v!("ursa_minor")
	set_limits!(umi_max)

	hidexdecorations!(ax, ticks=false, minorticks=false)

	l1 = lines!([NaN], [NaN], linewidth=smalllinewidth, linestyle=:dot, label="initial")
	l2 = lines!([NaN], [NaN], linewidth=smalllinewidth, linestyle=:solid, label="final")

	axislegend(ax_umi, [l1, l2], ["initial conditions", "end of isolation"], position=:rb, backgroundcolor=(:white, 0.8))
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═11cc55a3-d166-4bf3-b147-1c789d690f90
# ╠═b74a4535-1e12-422e-a524-135db4d8a9f7
# ╠═df4f20bf-97d9-408a-859e-e4070edcd0ef
# ╠═20c338e3-3d10-41f4-b8ee-d0eda4e755bd
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═50839b67-9514-47b9-add2-0c84d05f12da
# ╠═0dbc4770-fd1d-4a57-9825-3576900dd7a3
# ╠═890fa950-ca5e-4bec-9ef2-a2e5b3909909
# ╠═5a4853a9-1060-4be8-a872-c78297ed586f
# ╠═b7572e12-85a0-4f80-b7aa-ab2db8a81a13
# ╠═8ac764c8-ca9e-400b-bde4-f3e9b1ebf6b9
# ╠═4258c055-e11a-409d-9b53-441c55cd815f
# ╠═edc6003e-e931-45a0-82c4-4aa394a3f060
# ╠═1443385e-6d9b-4ddc-805c-2c12fb1c9067
# ╠═195582cf-83e3-44c8-9efe-b7d98f092c3c
# ╠═7034614b-eb03-4333-bf69-c3fb64a5c0fc
# ╠═ed95b670-3ba7-438a-8e4d-9fc6362d4d79
# ╠═ee01807c-785e-43b8-a8bd-b093aacbf281
# ╠═29637bda-3a48-4490-a2df-22ffc7fec59b
# ╠═430f1b00-9746-4e56-8608-deef86f1d0bb
# ╠═c47917c7-9eda-489c-bcf2-a602c01dcb24
# ╠═0813fc5e-236d-4665-b4f2-f4a1790c9aba
# ╠═ba3943bd-e94c-47fa-9856-827f1f28fc40
# ╠═fbaaf028-c25c-4e95-900e-c6fbb43e05bb
# ╠═6c106a00-445b-4e0f-9168-66a9aeb35767
# ╠═cf449735-0acb-4e83-bb6a-7b0bf5f794c4
# ╠═8e716378-3a9b-4994-87f1-e77a164d9f57
# ╠═3d0c926b-d2c6-4dcc-842e-f5be4f584fb8
# ╠═c1d64f7b-b7c3-4968-b564-0dba5c119e09
# ╠═4cb47cac-e3c6-41f3-ae8a-81f296414c52
# ╠═d447b0cd-3b58-4c7d-9abe-d1ffc4c6ae47
# ╠═b04cbab4-abce-4f3d-8541-ca3f00a6950b
# ╠═9b44a224-5fcb-4d28-80ad-83fe38a7ff7b
# ╠═8205bd52-8964-4054-9231-67a489b984f5
# ╠═ebf4a089-0459-4241-b05d-69e1c514a4a5
# ╠═cf800bd7-2cd1-4ac1-870b-ba29056e7920
# ╠═fb3dd9fe-2f1b-45fd-9d2d-3996859645b3
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═62407923-e1ba-445c-8726-2d72121a4ed2
