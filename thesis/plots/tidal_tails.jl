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

# ╔═╡ 7128f9f0-36fd-4ffb-9751-7c66c8cfd8e5
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ cb0f69df-e8c6-4e24-a79c-8e2cccdc3718
obs_labels = Dict(
	:radial_velocity => L"{v}_\textrm{los}\ /\ \textrm{km\,s^{-1}}",
	:radial_velocity_gsr => L"{v}_\textrm{los}'\ /\ \textrm{km\,s^{-1}}",
	:distance => "distance / kpc",
	:pmra => L"{\mu}_{\alpha*}\ /\ \textrm{mas\,yr^{-1}}",
	:pmra_gsr => L"{\mu}_{\alpha*}'\ /\ \textrm{mas\,yr^{-1}}",
	:pmdec => L"{\mu}_{\delta}\ /\ \textrm{mas\,yr^{-1}}",
	:pmdec_gsr => L"{\mu}_{\delta}'\ /\ \textrm{mas\,yr^{-1}}",
	:xi_p => L"\xi'\ /\ \textrm{arcmin}",
	:eta_p => L"\eta'\ /\ \textrm{arcmin}",
)

# ╔═╡ 05830907-1355-49ee-8cd3-df4678f33149
obs_dy = Dict(
	:radial_velocity => 40,
	:radial_velocity_gsr => 40,
	:distance => 20,
	:pmra => 0.2,
	:pmra_gsr =>  0.2,
	:pmdec =>  0.2,
	:pmdec_gsr =>  0.2,
	:eta_p => 5
)

# ╔═╡ 6cbf8dcd-d9e3-49e6-b799-a8bba9bf9211
import StatsBase: sample, weights

# ╔═╡ 0ccda86e-7a70-423c-a873-35b7a347bdfa
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ bcf2d214-84eb-42b4-95b0-da22ed4aecab
function add_xi_eta_p!(df, orbit_props)
	xi_p, eta_p = LilGuys.to_orbit_coords(df.ra, df.dec, orbit_props["ra_f"], orbit_props["dec_f"], orbit_props["theta0"])

	df[!, :xi_p] = xi_p .* 60
	df[!, :eta_p] = eta_p .* 60

	return df
end

# ╔═╡ 2c2b0102-43e6-4314-a6e4-855980f27149
function load_stars(galaxyname, modelname, starsname)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname)
	orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))

	starsfile = joinpath(model_dir, "stars", starsname, "final.fits")
	stars = read_fits(starsfile) 
	add_xi_eta_p!(stars, orbit_props)
	return stars
end

# ╔═╡ e9b30734-0237-4b52-9f3f-add7c5460743
function sample_stars(stars, N=100_000)
	idx = sample(1:size(stars, 1), weights(stars.weights), N)
	return stars[idx, :]
end

# ╔═╡ 904efe1c-6c2f-4e16-8468-ab4c5c79ebe1
function get_R_h(galaxyname)

	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxyname/observed_properties.toml"))
	R_h = obs_props["R_h"]
end

# ╔═╡ 226173f5-93c4-401e-be28-d509880aff00
function get_r_b(galaxyname, modelname, starsname; lmc=false)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$modelname/stars/$starsname/")

	prof_f = SurfaceDensityProfile(model_dir * "final_profile.toml")

	σv = prof_f.annotations["sigma_v"]
	if lmc
		props = TOML.parsefile(model_dir * "../../orbital_properties_lmc.toml")
	else
		props = TOML.parsefile(model_dir * "../../orbital_properties.toml")
	end

	dist_f =  TOML.parsefile(model_dir * "../../orbital_properties.toml")["distance_f"]

	
	dt = props["t_last_peri"]
	r_b = LilGuys.break_radius(σv / V2KMS, dt / T2GYR)

	return LilGuys.kpc2arcmin(r_b, dist_f)	
end

# ╔═╡ f2ec2966-9018-4ad0-8ec3-36e438012ca6
function scatter_coord(gs, stars_mean, ss, sym; xsym = :xi_p, kwargs...)
	x = ss[:, xsym] 
	
	r_max = 10 .* 60
	
	y0 =  stars_mean[sym]
	dy0 = obs_dy[sym]
	y = ss[:, sym]

	limits = (-r_max, r_max,  y0 - dy0, y0 + dy0)
	ax = Axis(gs,
		limits=limits,
		xlabel=obs_labels[xsym],
		ylabel=obs_labels[sym],
		
	)

	scatter!(x, y; rasterize=true, kwargs...)
	scatter!(0, y0, color=:black)
	return ax
end

# ╔═╡ b257827b-ca9c-4abe-bbaa-1097f0034b03
galaxy_labels = Dict(
	"sculptor" => "Sculptor", 
	"ursa_minor" => "Ursa Minor",
)

# ╔═╡ 6be80d25-974f-4bfe-96df-18cb0ce96c5a
function plot_sample(galaxyname, modelname, starsname)
	stars = load_stars(galaxyname, modelname, starsname)
	r_b = get_r_b(galaxyname, modelname, starsname)
	R_h = get_R_h(galaxyname)
	
	ss = sample_stars(stars)
	fig = Figure()

	for i in 1:4
		ysym = [:pmra_gsr, :pmdec_gsr, :distance, :radial_velocity_gsr][i]

		gs = fig[(i-1) ÷ 2 + 1, (i-1)%2 + 1]
		ax = scatter_coord(gs, stars[1, :], ss, ysym, markersize=2, alpha=1, color=COLORS[3])
		if i < 3
			ax.xlabel = ""
		end
		y = ax.finallimits[].origin[2] #- ax.finallimits[].widths[2]/2
		# @info y, r_b
		vlines!([-r_b, r_b], color=:black, linestyle=:dash, linewidth=theme(:linewidth)[]/2)
		if i == 1
			text!(r_b, y, text="break", rotation=π/2, align=(:left, :top), offset=(3, 6))
		end
	end

	Label(fig[0, :], text=galaxy_labels[galaxyname], font=:bold)
	fig

end

# ╔═╡ 8bf905de-5457-4611-8b9d-3ca6039cc582
md"""
# sculptors tidal tails
"""

# ╔═╡ 934b1a11-2646-4860-afa3-593218a9964b
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ 9ae7b57f-7649-443e-a438-4e5bd4b4d048
@savefig "scl_sim_stream" plot_sample(modelnames["scl_smallperi_plummer"]...)

# ╔═╡ ae2c737b-eca9-4ba1-b937-5e3dea3ece65
md"""
# Ursa Minor tidal tails
"""

# ╔═╡ 52c9a0fc-6db3-4e67-8d01-2f0eb0c1663f
@savefig "umi_sim_stream" plot_sample(modelnames["umi_smallperi_plummer"]...)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═cb0f69df-e8c6-4e24-a79c-8e2cccdc3718
# ╠═05830907-1355-49ee-8cd3-df4678f33149
# ╠═6cbf8dcd-d9e3-49e6-b799-a8bba9bf9211
# ╠═0ccda86e-7a70-423c-a873-35b7a347bdfa
# ╠═7128f9f0-36fd-4ffb-9751-7c66c8cfd8e5
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═2c2b0102-43e6-4314-a6e4-855980f27149
# ╠═bcf2d214-84eb-42b4-95b0-da22ed4aecab
# ╠═e9b30734-0237-4b52-9f3f-add7c5460743
# ╠═904efe1c-6c2f-4e16-8468-ab4c5c79ebe1
# ╠═226173f5-93c4-401e-be28-d509880aff00
# ╠═f2ec2966-9018-4ad0-8ec3-36e438012ca6
# ╠═b257827b-ca9c-4abe-bbaa-1097f0034b03
# ╠═6be80d25-974f-4bfe-96df-18cb0ce96c5a
# ╟─8bf905de-5457-4611-8b9d-3ca6039cc582
# ╠═934b1a11-2646-4860-afa3-593218a9964b
# ╠═9ae7b57f-7649-443e-a438-4e5bd4b4d048
# ╟─ae2c737b-eca9-4ba1-b937-5e3dea3ece65
# ╠═52c9a0fc-6db3-4e67-8d01-2f0eb0c1663f
