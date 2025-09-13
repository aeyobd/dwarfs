### A Pluto.jl notebook ###
# v0.20.17

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

# ╔═╡ 38901ab0-64df-4201-9492-0acc305376a1
function get_M_star(galaxyname)
	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))
	M_star = LilGuys.mag_to_L(obs_props["Mv"]) * obs_props["M_L_s"]
	return M_star / M2MSUN
end

# ╔═╡ 7a6d6345-fefb-4017-8230-3ff5bea55c83
function assign_weights!(snap, df)
	snap.weights = df.probability[snap.index]
end

# ╔═╡ 0dbc4770-fd1d-4a57-9825-3576900dd7a3
function load_stars(galaxyname, haloname, starsname)
	stars_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, "stars")

	df_probs = LilGuys.read_hdf5_table(joinpath(stars_dir, starsname, "probabilities_stars.hdf5"))

	snap_i = Snapshot(joinpath(stars_dir, "iso_paint.hdf5"))
	assign_weights!(snap_i, df_probs)

	snap_f = Snapshot(joinpath(stars_dir, "iso_final.hdf5"))
	assign_weights!(snap_f, df_probs)

	Mstar = get_M_star("sculptor")

	prof_s_i = DensityProfile(snap_i, snap_i.weights * Mstar)
	prof_i = DensityProfile(snap_i)
	prof_s_f = DensityProfile(snap_f, snap_f.weights * Mstar)
	prof_f = DensityProfile(snap_f)

	return (;dm_i = prof_i, dm_f = prof_f, stars_i=prof_s_i, stars_f=prof_s_f)
end

# ╔═╡ 8ac764c8-ca9e-400b-bde4-f3e9b1ebf6b9
scl_profs = load_stars("sculptor", "1e7_new_v31_r3.2", "exp2d_rs0.10")

# ╔═╡ ef66f5ee-31c9-461f-a843-5a8506e5fdbe
scl_profs_midres = load_stars("sculptor", "1e6_new_v31_r3.2", "exp2d_rs0.10")

# ╔═╡ 3d78046d-6846-483f-9477-db24747888de
scl_profs_lowres = load_stars("sculptor", "1e5_new_v31_r3.2", "exp2d_rs0.10")

# ╔═╡ edc6003e-e931-45a0-82c4-4aa394a3f060
umi_profs = load_stars("ursa_minor", "1e7_new_v38_r4.0", "exp2d_rs0.10")

# ╔═╡ 7034614b-eb03-4333-bf69-c3fb64a5c0fc
function plot_i_f(prof_i, prof_f; kwargs...)
	lines!(log_radii(prof_i), log_densities(prof_i), linestyle=:dot; kwargs...)
	lines!(log_radii(prof_f), log_densities(prof_f), linestyle=:solid; kwargs...)
end


# ╔═╡ b4e8c800-2c75-4ce5-b2e1-3205256fa4e0
smalllinewidth = theme(:linewidth)[]/2

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "initial_conditions" let
	fig = Figure(size=(5*72, 5*72))

	ax = Axis(fig[1,1], 
		xlabel = "log radius / kpc",
		ylabel = L"log $\rho$ / $10^{10}\,\textrm{M}_\odot\,\textrm{kpc}^{-3}$",
		limits = (-2, 1.5, -9.5, 0.5)
	)
	
	plot_i_f(scl_profs.dm_i, scl_profs.dm_f, color=COLORS[1], label="Scl DM")
	plot_i_f(scl_profs.stars_i, scl_profs.stars_f, color=COLORS[2], label="Scl stars", linewidth=smalllinewidth)
	hidexdecorations!(ticks=false, minorticks=false)

	axislegend(position=:rt, merge=true)

	l1 = lines!([NaN], [NaN], linestyle=:dot, label="initial")
	l2 = lines!([NaN], [NaN], linestyle=:solid, label="end of isolation")
	axislegend(ax, [l1, l2], ["initial profile", "end of isolation"], position=:lb, merge=true)


	ax = Axis(fig[2,1], 
		xlabel = "log radius / kpc",
		ylabel = L"log $\rho$ / $10^{10}\,\textrm{M}_\odot\,\textrm{kpc}^{-3}$",
		limits = (-2, 1.5, -9.5, 0.5)
	)
	

	plot_i_f(umi_profs.dm_i, umi_profs.dm_f, color=COLORS[1], label="UMi DM")
	plot_i_f(umi_profs.stars_i, umi_profs.stars_f, color=COLORS[2], label="UMi stars",  linewidth=smalllinewidth)
	axislegend(position=:rt, merge=true)

	
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═b74a4535-1e12-422e-a524-135db4d8a9f7
# ╠═20c338e3-3d10-41f4-b8ee-d0eda4e755bd
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═0dbc4770-fd1d-4a57-9825-3576900dd7a3
# ╠═38901ab0-64df-4201-9492-0acc305376a1
# ╠═7a6d6345-fefb-4017-8230-3ff5bea55c83
# ╠═8ac764c8-ca9e-400b-bde4-f3e9b1ebf6b9
# ╠═ef66f5ee-31c9-461f-a843-5a8506e5fdbe
# ╠═3d78046d-6846-483f-9477-db24747888de
# ╠═edc6003e-e931-45a0-82c4-4aa394a3f060
# ╠═7034614b-eb03-4333-bf69-c3fb64a5c0fc
# ╠═b4e8c800-2c75-4ce5-b2e1-3205256fa4e0
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
