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

# ╔═╡ df4f20bf-97d9-408a-859e-e4070edcd0ef
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 3abb0de9-81b7-4fe8-b457-92ae4b2e78e3
import CSV

# ╔═╡ 7a199993-0622-49cc-867b-8df8504447fe
import DataFrames: DataFrame

# ╔═╡ 11cc55a3-d166-4bf3-b147-1c789d690f90
import Agama

# ╔═╡ b74a4535-1e12-422e-a524-135db4d8a9f7
import TOML

# ╔═╡ 191b6df6-0fa4-4393-b69d-2aefeb3f9373
import StatsBase: quantile

# ╔═╡ d41f4781-5063-48c0-aade-fcd3980e19e7
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ 20c338e3-3d10-41f4-b8ee-d0eda4e755bd
CairoMakie.activate!(type=:png)

# ╔═╡ a4c87a6e-976d-4af6-868e-09cb85e3d424
module Utils
	include("utils.jl")
end

# ╔═╡ 50839b67-9514-47b9-add2-0c84d05f12da
function get_obs_props(galaxyname)
	return TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))
end

# ╔═╡ 5ac62ad0-6413-4d5f-a48f-080590ac862e
function get_M_star(galaxyname)
	obs_props = get_obs_props(galaxyname)

	return LilGuys.mag_to_L(obs_props["Mv"]) * obs_props["M_L_s"] / M2MSUN
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

# ╔═╡ 475816ca-baa6-4b51-a4b1-aeeec13ff105
function calc_σv(velocities, weights)
	σs2 = [LilGuys.std(velocities[i, :], weights)^2 for i in 1:3]

	return sqrt(sum(σs2))
end

# ╔═╡ 5864c349-dab7-46b2-9e62-bfd76d0c3b46
r_i =0.1

# ╔═╡ 48b2d945-e535-4e37-8735-df3c00b154f0
function calc_velocity_profile(snap)
	rs = LilGuys.radii(snap.positions[1:2, :])
	
	bins = LilGuys.bins_both(log10.(rs), snap.weights, bin_width=0.05, num_per_bin=2)

	N_bins = length(bins) - 1

	σv = zeros(N_bins)
	Ns = zeros(N_bins)
	
	for i in 1:N_bins
		filt = bins[i] .<= log10.(rs) .< bins[i+1]
		σ = LilGuys.std(snap.velocities[3, filt], snap.weights[filt])

		Ns[i] = sum(filt)
		σv[i] = σ
	end

	10 .^ midpoints(bins), σv, bins, Ns
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

# ╔═╡ e6796386-4b4d-4607-a53d-ad850b84c395
snap_i = load_snaps("sculptor", "1e6_new_v31_r3.2/orbit_smallperi", "exp2d_rs0.10")[1]

# ╔═╡ 9420a17b-b595-4bec-bef5-d30ee36dc6fc
filt_i = r_i*1.1 .> radii(snap_i.positions[1:2, :]) .> r_i

# ╔═╡ 98bc9898-343a-4328-8b66-2d534f72faa2
LilGuys.std(snap_i.velocities[3, filt_i], snap_i.weights[filt_i])*V2KMS

# ╔═╡ 34f6c3f4-07d2-4d5e-9a91-6102429a619e
LilGuys.std(snap_i.velocities[3, filt_i])*V2KMS

# ╔═╡ 890fa950-ca5e-4bec-9ef2-a2e5b3909909
function load_stars(galaxyname, modelname, starsname)

	snap_i, snap_f = load_snaps(galaxyname, modelname, starsname)
	Mstar = get_M_star(galaxyname)
	@info "stellar mass $Mstar"

	prof_s_i = MassProfile(snap_i, snap_i.weights * Mstar)
	prof_s_vel = calc_velocity_profile(snap_i)
	prof_i = MassProfile(snap_i)
	prof_s_f = MassProfile(snap_f, snap_f.weights * Mstar)
	prof_f = MassProfile(snap_f)
	prof_s_vel_f = calc_velocity_profile(snap_f)


	return (;dm_i = prof_i, dm_f = prof_f, stars_i=prof_s_i, stars_f=prof_s_f, vel_i=prof_s_vel, vel_f=prof_s_vel_f)
end

# ╔═╡ 6e76952a-c9b5-423b-a351-d0cc190270c5
calc_velocity_profile(snap_i)

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

# ╔═╡ 910bf7c0-f53b-49e8-bbc8-4cfc224f76a8
function plot_vel_prof!(profs; kwargs...)
	# prof = profs.vel_i
	# lines!(prof[1], prof[2] * V2KMS; kwargs...)


	prof = profs.vel_f
	lines!(prof[1], prof[2] * V2KMS; color=COLORS[2], kwargs...)

end

# ╔═╡ 50751c0d-1b4e-4d56-bb6a-6ce6b2cfef3c
function plot_prof(prof; text=nothing, color=nothing, kwargs...)
	x = radii(prof)
	y = middle.(LilGuys.circular_velocity(prof) * V2KMS)
	lines!(x, y; color=color, kwargs...)

	if !isnothing(text)
		Utils.text_along_line_log!(x, y, 0.05, text=text, color=color, align=(:left, :bottom), fontsize=smallfontsize)
	end

end

# ╔═╡ 7034614b-eb03-4333-bf69-c3fb64a5c0fc
function plot_i_f(profs; kwargs...)
	# plot_prof(profs.dm_i, text="dark matter", linestyle=:dot; kwargs...)
	plot_prof(profs.dm_f;text="dark matter", kwargs...)
end

# ╔═╡ ee01807c-785e-43b8-a8bd-b093aacbf281
function plot_i_f_stars(profs; kwargs...)
	# plot_prof(profs.stars_i, text="stars", linestyle=:dot; kwargs...)
	plot_prof(profs.stars_f; text="stars", kwargs...)
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

# ╔═╡ 8776586d-5686-4cd9-9c2a-f54b7f16e90e
function sigmav_axis(gs)
	ax = Axis(gs, 
		xlabel = "radius / kpc",
		ylabel = L"$\sigma_\textrm{v, los}$ / km\,s$^{-1}$",
			  xscale = log10,
			  yscale = log10,
			  xticks = [0.1, 1, 10],
			  # yticks = [1; 10:10:40],
			  yminorticks=1:20,
			  limits=(0.03, 3, 5, 15)
	)
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

# ╔═╡ da1a2d1b-d32c-42c6-957d-82e39fcda005
function ρ_to_v(ρ, r)
	v_scale = sqrt(ρ * 4π/3)
	v = v_scale .* r
end

# ╔═╡ cd74f6a3-5d9c-497e-b4e5-6c7693c62bca
function plot_density_line!(ρ, ρ_el, ρ_ep)
	
	r = 10 .^ LinRange(-2, 2, 10)
	v = ρ_to_v(3ρ, r)
	lines!(r, v*V2KMS, color=:black, linewidth=smalllinewidth)
	
	v_l = ρ_to_v(3ρ + 3ρ_el, r)
	v_h = ρ_to_v(3ρ + 3ρ_ep, r)
	band!(r, v_l*V2KMS, v_h*V2KMS, alpha=0.5, color=:black)

	Utils.text_along_line_log!(r, v*V2KMS, 2, text=L"3\bar\rho_\textrm{MW,\ peri}", align=(:center, :bottom), fontsize=smallfontsize)
end

# ╔═╡ 961f2851-8ebb-49ac-9459-5e6f543f80f4
function calc_r_J(halo, ρ_host)
	return LilGuys.find_zero(r -> 3ρ_host - LilGuys.mean_density(halo, r), 1)
end

# ╔═╡ 99675829-e4b2-4124-90c2-595e53290956
function plot_r_J(halo, ρ_host)
	r_J = calc_r_J(halo, ρ_host)
	v0 = LilGuys.v_circ(halo, r_J) * V2KMS

	v0 *= 10^-0.5
	v1 = 10^-0.15 * v0

	ms = 4
	# arrows2d!([r_J], [v1], [0], [v0 - v1], color=:grey,minshaftlength=0, shaftwidth=smalllinewidth, tiplength=ms, tipwidth=ms * 2/sqrt(3))

	vlines!(r_J, color=:grey, linewidth=smalllinewidth,  linestyle=:dot)
	text!(r_J, v1, text=L"R_J", fontsize=smallfontsize, align=(:left, :center), color=:grey,)
end

# ╔═╡ 07a499bb-a861-43fe-b54a-48348141bb0d
calc_r_J(NFW(r_circ_max=3.2, v_circ_max=31/V2KMS), 8.5e-5)

# ╔═╡ da7cc3dd-ee03-42c3-92ee-3e256bd767c3
calc_r_J(NFW(r_circ_max=2.5, v_circ_max=25/V2KMS), 8.5e-5)

# ╔═╡ bbde983a-c98a-4eeb-b04c-41ec02be934d
# get_mean_density("sculptor", "vasiliev24_L3M11", units=Agama.VASILIEV_UNITS)

# ╔═╡ 5ddbbeb3-08d1-4092-9d73-bba36cd18e1b


# ╔═╡ 16ce9a3c-69b3-486f-9c22-a3c483cf9570
df_scl = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/velocities/processed/vz_r_ell_binned.rv_combined_x_wide_2c_psat_0.2.csv"), DataFrame)

# ╔═╡ e7eada5c-6d8a-48b6-a65a-d92126bf37e2
df_umi = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor/velocities/processed/vz_r_ell_binned.rv_combined_x_2c_psat_0.2.csv"), DataFrame)

# ╔═╡ 63f93057-14b5-4df9-8e67-5bed7452bbe4
function plot_σv_obs!(galaxyname)
	if galaxyname == "sculptor"
		df = df_scl
	elseif galaxyname == "ursa_minor"
		df = df_umi
	end

	dist = get_obs_props(galaxyname)["distance"]

	x = LilGuys.arcmin2kpc.(df.x, dist)
	errorscatter!(x, (df.σ), yerror=df.σ_em, color=:black, markersize=6)
end

# ╔═╡ 57519fa1-e8a2-4baa-a4fa-af4482328164
COLORS

# ╔═╡ 9721c592-a0f2-4564-83e4-da2eff4c3c79
let
	fig = Figure()
	ax = Axis(fig[1,1])

	plot_vel_prof!(scl_profs)
	plot_σv_obs!("sculptor")

	xlims!(0, 1.5)

	fig
end

# ╔═╡ 6720afb1-9968-42e3-a22c-c9232309f586
function load_peris(galaxyname, modelname="EP2020"; lmc=false)
	df = read_fits(joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, modelname, "orbital_properties.fits"))

	if lmc
		col = "pericentre_lmc"
	else
		col = "pericentre"
	end
	
	return df[!, col]
end

# ╔═╡ 6b0b3b07-6adf-497b-a513-9ecfa47d97c6
function get_mean_density(galaxyname, modelname="EP2020"; units=Agama.AgamaUnits())
	peris = load_peris(galaxyname, modelname)

	pot = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, modelname, "agama_potential.ini"))
	ρs =  Agama.enclosed_mass(pot, peris, units) ./ (4π/3 * peris .^ 3)

	l, m, h = quantile(ρs, [0.16, 0.5, 0.84])
	return m, l-m, h-m
end

# ╔═╡ cb9f60e7-aa33-4938-a459-eb6145ac555e
get_mean_density("sculptor")

# ╔═╡ 62407923-e1ba-445c-8726-2d72121a4ed2
function set_limits!(maxes)
	xlims!(0.03, maxes.r_circ_max * 10^1)
	ylims!(maxes.v_circ_max*V2KMS* 10^-1.9, maxes.v_circ_max*V2KMS * 10^0.1)

end

# ╔═╡ a0579c53-e244-4e3a-88e1-6569d4c6efff
pot = Agama.Potential(file = joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))


# ╔═╡ b200c0c7-fd55-483e-8a37-6f67f285d68f
pot_lmc = Agama.Potential(file = joinpath(ENV["DWARFS_ROOT"], "agama/potentials/vasiliev24/L3M11/potential_lmc_init.ini"))


# ╔═╡ 192ced23-9b3c-4799-8e8f-d2da7a1a6234
function get_mean_density_lmc(galaxyname)
	peris = load_peris(galaxyname, "vasiliev24_L3M11", lmc=true)
	ρs =  Agama.enclosed_mass(pot_lmc, peris, Agama.VASILIEV_UNITS) ./ (4π/3 * peris .^ 3)

	l, m, h = quantile(ρs, [0.16, 0.5, 0.84])
	return m, l-m, h-m
end

# ╔═╡ 590d652b-6263-486a-9f4d-4f06c1927831
get_mean_density_lmc("sculptor")

# ╔═╡ 1443385e-6d9b-4ddc-805c-2c12fb1c9067
scl_max = LilGuys.fit_v_r_circ_max(radii(scl_profs[1]), middle.(LilGuys.circular_velocity(scl_profs[1])))

# ╔═╡ 195582cf-83e3-44c8-9efe-b7d98f092c3c
umi_max = LilGuys.fit_v_r_circ_max(radii(umi_profs[1]), middle.(LilGuys.circular_velocity(umi_profs[1])))

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "initial_velocity" let
	fig = Figure(size=(3.5, 3.5) .*72)

	ax = vcirc_axis(fig[1,1])
	ax.title="Sculptor"
	
	plot_i_f(scl_profs, color=COLORS[5], label="Scl DM")
	plot_i_f_stars(scl_profs, 
				   color=COLORS[2], label="Scl stars",  linewidth=smalllinewidth)

	plot_r_J(NFW(r_circ_max=3.2, v_circ_max=31/V2KMS), get_mean_density("sculptor")[1])
	
	set_limits!(scl_max)

	ax_scl_sigma = sigmav_axis(fig[2, 1])
	plot_vel_prof!(scl_profs)
	# plot_sigma_v!("sculptor")
	plot_σv_obs!("sculptor")

	
	ax_umi = vcirc_axis(fig[1,2])
	ax_umi.title = "Ursa Minor"
	
	plot_i_f(umi_profs, color=COLORS[5], label="UMi DM")
	plot_i_f_stars(umi_profs,
				   color=COLORS[2], label="UMi stars", linewidth=smalllinewidth)

	plot_r_J(NFW(r_circ_max=4, v_circ_max=38/V2KMS), get_mean_density("ursa_minor")[1])

	set_limits!(umi_max)


	ax_umi_sigma = sigmav_axis(fig[2, 2])
	plot_vel_prof!(umi_profs)
	# plot_sigma_v!("ursa_minor")
	plot_σv_obs!("ursa_minor")

	
	hideydecorations!(ax_umi, ticks=false, minorticks=false)
	hideydecorations!(ax_umi_sigma, ticks=false, minorticks=false)

	# l1 = lines!([NaN], [NaN], linewidth=smalllinewidth, linestyle=:dot, label="initial")
	# l2 = lines!([NaN], [NaN], linewidth=smalllinewidth, linestyle=:solid, label="final")
	# linkyaxes!(ax, ax_umi)
	# colgap!(fig.layout, 0)

	# axislegend(ax, [l1, l2], ["initial conditions", "end of isolation"], position=:lb, backgroundcolor=(:white, 0.8))
	fig
end

# ╔═╡ 7b2c0281-dba8-4e17-9e11-ef8acdc2eba6
@savefig "initial_velocity_nosigma" let
	fig = Figure(size=(3.5, 2) .*72)

	ax = vcirc_axis(fig[1,1])
	ax.title="Sculptor"
	
	plot_i_f(scl_profs, color=COLORS[5], label="Scl DM")
	plot_i_f_stars(scl_profs, 
				   color=COLORS[2], label="Scl stars",  linewidth=smalllinewidth)

	plot_r_J(NFW(r_circ_max=3.2, v_circ_max=31/V2KMS), get_mean_density("sculptor")[1])
	plot_sigma_v!("sculptor")


	set_limits!(scl_max)

	
	ax_umi = vcirc_axis(fig[1,2])
	ax_umi.title = "Ursa Minor"
	
	plot_i_f(umi_profs, color=COLORS[5], label="UMi DM")
	plot_i_f_stars(umi_profs,
				   color=COLORS[2], label="UMi stars", linewidth=smalllinewidth)

	plot_r_J(NFW(r_circ_max=4, v_circ_max=38/V2KMS), get_mean_density("ursa_minor")[1])

	plot_sigma_v!("ursa_minor")
	hideydecorations!(ticks=false, minorticks=false)


	set_limits!(umi_max)


	fig
end

# ╔═╡ ad4bdf3a-4ec5-4113-b2fe-4f3e4247555c
function change_in_vcirc(profs, r)
	prof_i, prof_f = profs.dm_i, profs.dm_f

	M_i = LilGuys.lerp(log10.(prof_i.radii), middle.(prof_i.M_in))(log10(r))
	M_f = LilGuys.lerp(log10.(prof_f.radii), middle.(prof_f.M_in))(log10(r))

	return M_f/M_i
end

# ╔═╡ 7c132ca4-f235-42fd-8495-6d4d12881d07
change_in_vcirc(scl_profs, 0.1 * 3.2 / 4)

# ╔═╡ e084e5e8-e13e-4828-915f-f1c6244e594b
change_in_vcirc(umi_profs, 0.1)

# ╔═╡ 2ea0bba5-5a65-4e95-b5de-1ab6595a94b0
6/4

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═3abb0de9-81b7-4fe8-b457-92ae4b2e78e3
# ╠═7a199993-0622-49cc-867b-8df8504447fe
# ╠═11cc55a3-d166-4bf3-b147-1c789d690f90
# ╠═b74a4535-1e12-422e-a524-135db4d8a9f7
# ╠═df4f20bf-97d9-408a-859e-e4070edcd0ef
# ╠═191b6df6-0fa4-4393-b69d-2aefeb3f9373
# ╠═d41f4781-5063-48c0-aade-fcd3980e19e7
# ╠═20c338e3-3d10-41f4-b8ee-d0eda4e755bd
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═a4c87a6e-976d-4af6-868e-09cb85e3d424
# ╠═50839b67-9514-47b9-add2-0c84d05f12da
# ╠═5ac62ad0-6413-4d5f-a48f-080590ac862e
# ╟─43dad299-d8d2-4146-8320-c90a62f3a3f0
# ╠═0dbc4770-fd1d-4a57-9825-3576900dd7a3
# ╠═475816ca-baa6-4b51-a4b1-aeeec13ff105
# ╠═5864c349-dab7-46b2-9e62-bfd76d0c3b46
# ╠═9420a17b-b595-4bec-bef5-d30ee36dc6fc
# ╠═98bc9898-343a-4328-8b66-2d534f72faa2
# ╠═34f6c3f4-07d2-4d5e-9a91-6102429a619e
# ╠═48b2d945-e535-4e37-8735-df3c00b154f0
# ╠═e6796386-4b4d-4607-a53d-ad850b84c395
# ╠═d66619f4-de41-45f4-9be1-4b938b9a959c
# ╠═890fa950-ca5e-4bec-9ef2-a2e5b3909909
# ╠═6e76952a-c9b5-423b-a351-d0cc190270c5
# ╠═8ac764c8-ca9e-400b-bde4-f3e9b1ebf6b9
# ╠═edc6003e-e931-45a0-82c4-4aa394a3f060
# ╟─dd6e4582-d610-4669-8fdf-b458bacb86b3
# ╠═ec7d88d3-7503-4c18-9088-04319c6f2e98
# ╠═cd74f6a3-5d9c-497e-b4e5-6c7693c62bca
# ╠═7034614b-eb03-4333-bf69-c3fb64a5c0fc
# ╠═ee01807c-785e-43b8-a8bd-b093aacbf281
# ╠═910bf7c0-f53b-49e8-bbc8-4cfc224f76a8
# ╠═50751c0d-1b4e-4d56-bb6a-6ce6b2cfef3c
# ╠═c1d64f7b-b7c3-4968-b564-0dba5c119e09
# ╠═ebf4a089-0459-4241-b05d-69e1c514a4a5
# ╠═fb3dd9fe-2f1b-45fd-9d2d-3996859645b3
# ╠═8776586d-5686-4cd9-9c2a-f54b7f16e90e
# ╠═6c106a00-445b-4e0f-9168-66a9aeb35767
# ╠═da1a2d1b-d32c-42c6-957d-82e39fcda005
# ╠═961f2851-8ebb-49ac-9459-5e6f543f80f4
# ╠═99675829-e4b2-4124-90c2-595e53290956
# ╠═cb9f60e7-aa33-4938-a459-eb6145ac555e
# ╠═590d652b-6263-486a-9f4d-4f06c1927831
# ╠═07a499bb-a861-43fe-b54a-48348141bb0d
# ╠═da7cc3dd-ee03-42c3-92ee-3e256bd767c3
# ╠═bbde983a-c98a-4eeb-b04c-41ec02be934d
# ╠═5ddbbeb3-08d1-4092-9d73-bba36cd18e1b
# ╠═16ce9a3c-69b3-486f-9c22-a3c483cf9570
# ╠═e7eada5c-6d8a-48b6-a65a-d92126bf37e2
# ╠═63f93057-14b5-4df9-8e67-5bed7452bbe4
# ╠═57519fa1-e8a2-4baa-a4fa-af4482328164
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═7b2c0281-dba8-4e17-9e11-ef8acdc2eba6
# ╠═9721c592-a0f2-4564-83e4-da2eff4c3c79
# ╠═6720afb1-9968-42e3-a22c-c9232309f586
# ╠═6b0b3b07-6adf-497b-a513-9ecfa47d97c6
# ╠═192ced23-9b3c-4799-8e8f-d2da7a1a6234
# ╠═62407923-e1ba-445c-8726-2d72121a4ed2
# ╠═a0579c53-e244-4e3a-88e1-6569d4c6efff
# ╠═b200c0c7-fd55-483e-8a37-6f67f285d68f
# ╠═1443385e-6d9b-4ddc-805c-2c12fb1c9067
# ╠═195582cf-83e3-44c8-9efe-b7d98f092c3c
# ╠═ad4bdf3a-4ec5-4113-b2fe-4f3e4247555c
# ╠═7c132ca4-f235-42fd-8495-6d4d12881d07
# ╠═e084e5e8-e13e-4828-915f-f1c6244e594b
# ╠═2ea0bba5-5a65-4e95-b5de-1ab6595a94b0
