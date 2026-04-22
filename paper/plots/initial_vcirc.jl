### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys

	import PyCall
	import PyPlot as plt

	import Agama
	import TOML
	using PyFITS
end

# ╔═╡ 191b6df6-0fa4-4393-b69d-2aefeb3f9373
import StatsBase: quantile

# ╔═╡ 53a732df-3076-465a-97b9-089da190d73e
PyCall.@pyimport arya

# ╔═╡ 282c8c5c-e6f5-4b4e-b4f3-d6df8cb32bf0
arya.init("apj")

# ╔═╡ d41f4781-5063-48c0-aade-fcd3980e19e7
modelnames = TOML.parsefile("model_key.toml")

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

	10 .^ LilGuys.midpoints(bins), σv, bins, Ns
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

# ╔═╡ 8ac764c8-ca9e-400b-bde4-f3e9b1ebf6b9
scl_profs = load_stars(modelnames["scl_smallperi"]...)

# ╔═╡ edc6003e-e931-45a0-82c4-4aa394a3f060
umi_profs = load_stars(modelnames["umi_smallperi"]...)

# ╔═╡ dd6e4582-d610-4669-8fdf-b458bacb86b3
md"""
# Plot utils
"""

# ╔═╡ a0579c53-e244-4e3a-88e1-6569d4c6efff
pot = Agama.Potential(file = joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))


# ╔═╡ c1d64f7b-b7c3-4968-b564-0dba5c119e09
α_3d = LilGuys.r_h(LilGuys.Exp2D()) / LilGuys.R_h(LilGuys.Exp2D()) 

# ╔═╡ 6c106a00-445b-4e0f-9168-66a9aeb35767
smalllinewidth=1

# ╔═╡ ec7d88d3-7503-4c18-9088-04319c6f2e98
smallfontsize=0.8 * 10

# ╔═╡ fb3dd9fe-2f1b-45fd-9d2d-3996859645b3
function vcirc_axis(ax)
	ax.set(
		xlabel = raw"radius $/$ kpc",
		ylabel = raw"circular velocity $/$ km\,s$^{-1}$",
		xscale = "log",
		yscale = "log",
		xticks = [0.1, 1, 10],
		yticks = [1; 10],
		# yminorticks=[1:10; 10:5:20; 20:5:40],
		xlim = (0.03, 100),
		ylim = (1, 40),
	)
end

# ╔═╡ 50751c0d-1b4e-4d56-bb6a-6ce6b2cfef3c
function plot_prof!(prof; text=nothing, color=nothing, kwargs...)
	x = radii(prof)
	y = middle.(LilGuys.circular_velocity(prof) * V2KMS)
	plt.plot(x, y; color=color, kwargs...)

	if !isnothing(text)
		# Utils.text_along_line_log!(x, y, 0.05, text=text, color=color, align=(:left, :bottom), fontsize=smallfontsize)
		idx = argmin(abs.(x .- 0.05))
		plt.annotate(text, (x[idx], y[idx]), color=color, 
					 xytext=(0, smallfontsize/sqrt(2)),
					 textcoords = "offset points",
				ha="left", va="bottom", rotation=35, fontsize=smallfontsize)
	end

end

# ╔═╡ 7034614b-eb03-4333-bf69-c3fb64a5c0fc
function plot_f!(profs; kwargs...)
	plot_prof!(profs.dm_f;text="dark matter", kwargs...)
end

# ╔═╡ ee01807c-785e-43b8-a8bd-b093aacbf281
function plot_f_stars!(profs; kwargs...)
	plot_prof!(profs.stars_f; text="stars", kwargs...)
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
	plt.plot(r, v*V2KMS, color="k", linewidth=smalllinewidth)
	
	v_l = ρ_to_v(3ρ + 3ρ_el, r)
	v_h = ρ_to_v(3ρ + 3ρ_ep, r)
	plt.fill_between(r, v_l*V2KMS, v_h*V2KMS, alpha=0.5, color="k")

	Utils.text_along_line_log!(r, v*V2KMS, 2, text=raw"3\bar\rho_\textrm{MW,\ peri}", align=(:center, :bottom), fontsize=smallfontsize)
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

	plt.axvline(r_J, color=:grey, linewidth=smalllinewidth,  ls=":")
	plt.annotate(raw"$r_J$", (r_J, v1), fontsize=smallfontsize, ha="left",
			 va="center", color=(0.5, 0.5, 0.5),
				   xytext=(smallfontsize/2, 0),
				   textcoords = "offset points",
				  )
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

# ╔═╡ 62407923-e1ba-445c-8726-2d72121a4ed2
function set_limits!(maxes)
	plt.xlim(0.03, maxes.r_circ_max * 10^1)
	plt.ylim(maxes.v_circ_max*V2KMS* 10^-1.9, maxes.v_circ_max*V2KMS * 10^0.1)

end

# ╔═╡ 5ddbbeb3-08d1-4092-9d73-bba36cd18e1b
md"""
# The plot
"""

# ╔═╡ f653954b-26b1-43f9-bed0-6ba9ab7e9456
COLORS = arya.COLORS

# ╔═╡ ebf4a089-0459-4241-b05d-69e1c514a4a5
function plot_sigma_v!(galaxy)
	obs_props = get_obs_props(galaxy) |> LilGuys.collapse_errors

	R_h = [LilGuys.arcmin2kpc(obs_props["R_h"], obs_props["distance"])]
	σv = [obs_props["sigma_v"]]
	
	plt.scatter(middle.(R_h), middle.(σv), 
				  # xerror = error_interval.(R_h), 
				  # yerror=error_interval.(σv), 
				 # fmt="o",
				  # linewidth=smalllinewidth, 
				color=COLORS[3], 
				# markersize=0
			   )
	
	plt.annotate(raw"($\sigma_\textrm{v}$, $R_h$)",
				 (middle(R_h[1]), middle(σv[1])),
				 va="center",
		  xytext=(smallfontsize/2, 0), 
				 textcoords = "offset points",
		  color=COLORS[3], fontsize=smallfontsize,
				  )
end

# ╔═╡ 1443385e-6d9b-4ddc-805c-2c12fb1c9067
scl_max = LilGuys.fit_v_r_circ_max(radii(scl_profs[1]), middle.(LilGuys.circular_velocity(scl_profs[1])))

# ╔═╡ 195582cf-83e3-44c8-9efe-b7d98f092c3c
umi_max = LilGuys.fit_v_r_circ_max(radii(umi_profs[1]), middle.(LilGuys.circular_velocity(umi_profs[1])))

# ╔═╡ 2e729ef4-8a9e-431e-84a0-41c4b766c870
rcParams = plt.PyDict(plt.matplotlib."rcParams")


# ╔═╡ da8eb034-320d-411c-8ded-33275de623f9
rcParams["axes.titlesize"] = 10

# ╔═╡ 7b2c0281-dba8-4e17-9e11-ef8acdc2eba6
let
	fig, axs = plt.subplots(1, 2, sharey=true, gridspec_kw=Dict("wspace" => 0), 
						   figsize = (3.333, 3.333*2/4))

	vcirc_axis(axs[1])
	plt.sca(axs[1])
	axs[1].set_title(raw"\textbf{Sculptor}")
	
	plot_f!(scl_profs, color=arya.COLORS[5], label="Scl DM")
	plot_f_stars!(scl_profs, 
				   color=COLORS[2], label="Scl stars",  linewidth=smalllinewidth)

	plot_r_J(NFW(r_circ_max=3.2, v_circ_max=31/V2KMS), get_mean_density("sculptor")[1])
	plot_sigma_v!("sculptor")


	set_limits!(scl_max)

	
	vcirc_axis(axs[2])
	plt.sca(axs[2])
	axs[2].set_title(raw"\textbf{Ursa Minor}")
	axs[2].set_ylabel("")

	plot_f!(umi_profs, color=COLORS[5], label="UMi DM")
	plot_f_stars!(umi_profs,
				   color=COLORS[2], label="UMi stars", linewidth=smalllinewidth)

	plot_r_J(NFW(r_circ_max=4, v_circ_max=38/V2KMS), get_mean_density("ursa_minor")[1])

	plot_sigma_v!("ursa_minor")


	set_limits!(umi_max)

	plt.savefig("figures/initial_velocity_nosigma.pdf")


	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═191b6df6-0fa4-4393-b69d-2aefeb3f9373
# ╠═53a732df-3076-465a-97b9-089da190d73e
# ╠═282c8c5c-e6f5-4b4e-b4f3-d6df8cb32bf0
# ╠═d41f4781-5063-48c0-aade-fcd3980e19e7
# ╠═a4c87a6e-976d-4af6-868e-09cb85e3d424
# ╠═50839b67-9514-47b9-add2-0c84d05f12da
# ╠═5ac62ad0-6413-4d5f-a48f-080590ac862e
# ╟─43dad299-d8d2-4146-8320-c90a62f3a3f0
# ╠═48b2d945-e535-4e37-8735-df3c00b154f0
# ╠═e6796386-4b4d-4607-a53d-ad850b84c395
# ╠═d66619f4-de41-45f4-9be1-4b938b9a959c
# ╠═890fa950-ca5e-4bec-9ef2-a2e5b3909909
# ╠═8ac764c8-ca9e-400b-bde4-f3e9b1ebf6b9
# ╠═edc6003e-e931-45a0-82c4-4aa394a3f060
# ╟─dd6e4582-d610-4669-8fdf-b458bacb86b3
# ╠═a0579c53-e244-4e3a-88e1-6569d4c6efff
# ╠═c1d64f7b-b7c3-4968-b564-0dba5c119e09
# ╠═6c106a00-445b-4e0f-9168-66a9aeb35767
# ╠═ec7d88d3-7503-4c18-9088-04319c6f2e98
# ╠═fb3dd9fe-2f1b-45fd-9d2d-3996859645b3
# ╠═cd74f6a3-5d9c-497e-b4e5-6c7693c62bca
# ╠═7034614b-eb03-4333-bf69-c3fb64a5c0fc
# ╠═ee01807c-785e-43b8-a8bd-b093aacbf281
# ╠═50751c0d-1b4e-4d56-bb6a-6ce6b2cfef3c
# ╠═ebf4a089-0459-4241-b05d-69e1c514a4a5
# ╠═da1a2d1b-d32c-42c6-957d-82e39fcda005
# ╠═961f2851-8ebb-49ac-9459-5e6f543f80f4
# ╠═99675829-e4b2-4124-90c2-595e53290956
# ╠═6720afb1-9968-42e3-a22c-c9232309f586
# ╠═6b0b3b07-6adf-497b-a513-9ecfa47d97c6
# ╠═62407923-e1ba-445c-8726-2d72121a4ed2
# ╟─5ddbbeb3-08d1-4092-9d73-bba36cd18e1b
# ╠═f653954b-26b1-43f9-bed0-6ba9ab7e9456
# ╠═1443385e-6d9b-4ddc-805c-2c12fb1c9067
# ╠═195582cf-83e3-44c8-9efe-b7d98f092c3c
# ╠═2e729ef4-8a9e-431e-84a0-41c4b766c870
# ╠═da8eb034-320d-411c-8ded-33275de623f9
# ╠═7b2c0281-dba8-4e17-9e11-ef8acdc2eba6
