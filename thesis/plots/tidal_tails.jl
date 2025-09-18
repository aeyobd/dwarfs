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
	:xi_p => L"\xi'\ /\ \textrm{degrees}",
	:eta_p => L"\eta'\ /\ \textrm{degrees}",
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

	df[!, :xi_p] = xi_p
	df[!, :eta_p] = eta_p

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

# ╔═╡ 30ab5003-4218-4a33-a0e6-f97826e56b1d


# ╔═╡ e9b30734-0237-4b52-9f3f-add7c5460743
function sample_stars(stars, N=10_000)
	idx = sample(1:size(stars, 1), weights(stars.weights), 100_000)
	return stars[idx, :]
end

# ╔═╡ 904efe1c-6c2f-4e16-8468-ab4c5c79ebe1


# ╔═╡ f2ec2966-9018-4ad0-8ec3-36e438012ca6
function scatter_coord(gs, stars_mean, ss, sym; xsym = :xi_p, kwargs...)
	x = ss[:, xsym]
	
	r_max = 10
	
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

# ╔═╡ 6be80d25-974f-4bfe-96df-18cb0ce96c5a
function plot_sample(stars)
	ss = sample_stars(stars)
	fig = Figure()

	for i in 1:4
		ysym = [:pmra_gsr, :pmdec_gsr, :distance, :radial_velocity_gsr][i]

		gs = fig[(i-1) ÷ 2 + 1, (i-1)%2 + 1]
		ax = scatter_coord(gs, stars[1, :], ss, ysym, markersize=2, alpha=1, color=COLORS[3])
		if i < 3
			ax.xlabel = ""
		end
	end

	fig

end

# ╔═╡ 8bf905de-5457-4611-8b9d-3ca6039cc582
md"""
# sculptors tidal tails
"""

# ╔═╡ 189db37c-1e59-4e69-9487-8525d45a9c3c
plot_sample(load_stars("sculptor", "1e6_new_v31_r3.2/orbit_mean", "plummer_rs0.20"))

# ╔═╡ f62a2c26-1baa-4a46-98ac-e43f75547a72
plot_sample(load_stars("sculptor", "1e7_new_v31_r3.2/orbit_smallperi", "plummer_rs0.20"))

# ╔═╡ 53c9aaff-51cd-409e-ab1f-c28e940b2ed1
plot_sample(load_stars("sculptor", "1e7_new_v25_r2.5/smallperilmc", "plummer_rs0.20"))

# ╔═╡ 7b274c9f-4cdd-41ef-bdef-8a649331d1d9
plot_sample(load_stars("sculptor", "1e6_new_v31_r3.2/L3M11_9Gyr_smallperi.a4", "plummer_rs0.20"))

# ╔═╡ 5a56e59a-c846-4829-9ddf-496c8a0383fb
plot_sample(load_stars("sculptor", "1e6_new_v31_r3.2/L3M11_9Gyr_smallperi.a4", "exp2d_rs0.10"))

# ╔═╡ ae2c737b-eca9-4ba1-b937-5e3dea3ece65
md"""
# Ursa Minor tidal tails
"""

# ╔═╡ 52c9a0fc-6db3-4e67-8d01-2f0eb0c1663f
@savefig "umi_sim_stream" plot_sample(load_stars("ursa_minor", "1e7_new_v38_r4.0/orbit_smallperi.5", "plummer_rs0.20"))

# ╔═╡ b50e2128-aa69-4ef9-bd7a-a08f412c504e
plot_sample(load_stars("ursa_minor", "1e6_v37_r5.0/orbit_mean.2", "plummer_rs0.20"))

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
# ╠═30ab5003-4218-4a33-a0e6-f97826e56b1d
# ╠═e9b30734-0237-4b52-9f3f-add7c5460743
# ╠═904efe1c-6c2f-4e16-8468-ab4c5c79ebe1
# ╠═f2ec2966-9018-4ad0-8ec3-36e438012ca6
# ╠═6be80d25-974f-4bfe-96df-18cb0ce96c5a
# ╟─8bf905de-5457-4611-8b9d-3ca6039cc582
# ╠═189db37c-1e59-4e69-9487-8525d45a9c3c
# ╠═f62a2c26-1baa-4a46-98ac-e43f75547a72
# ╠═53c9aaff-51cd-409e-ab1f-c28e940b2ed1
# ╠═7b274c9f-4cdd-41ef-bdef-8a649331d1d9
# ╠═5a56e59a-c846-4829-9ddf-496c8a0383fb
# ╟─ae2c737b-eca9-4ba1-b937-5e3dea3ece65
# ╠═52c9a0fc-6db3-4e67-8d01-2f0eb0c1663f
# ╠═b50e2128-aa69-4ef9-bd7a-a08f412c504e
