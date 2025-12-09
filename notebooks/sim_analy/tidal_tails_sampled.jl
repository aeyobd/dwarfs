### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()


	using CairoMakie
	using Arya

end

# ╔═╡ 7dd7aa41-6fa4-49e0-a3d8-43768a4d0446
using PlutoUI

# ╔═╡ 7128f9f0-36fd-4ffb-9751-7c66c8cfd8e5
using PyFITS

# ╔═╡ 69c1ad4b-d6b0-4c0e-8623-9622a3779ec7
import DataFrames

# ╔═╡ 5dff51c5-2c62-4469-8744-ea1a0a7c932b
r_max = 10 .* 60

# ╔═╡ 689209ef-3584-486c-9010-bf9bb59238f6
N = 30_000

# ╔═╡ daf8bd82-89c1-496f-98ed-333b740095d9
function notebook_inputs(; kwargs...)
	return PlutoUI.combine() do Child
		
		user_inputs = [
			md""" $(string(name)): $(
				Child(name, obj)
			)"""
			
			for (name, obj) in kwargs
		]
		
		md"""
		#### Inputs
		$(user_inputs)
		"""
	end
end

# ╔═╡ d593ea8b-47ec-4eca-a21b-ba07a118f468
@bind inputs confirm(notebook_inputs(;
	galaxyname = TextField(default="sculptor"),
	modelname = TextField(default="1e7_new_v31_r3.2/orbit_example", 40),
	starsname = TextField(default="exp2d_rs0.13"),
))

# ╔═╡ 8f0d73ac-58fe-4814-82ce-3295abaaf8bd
galaxyname, modelname, starsname = inputs.galaxyname, inputs.modelname, inputs.starsname

# ╔═╡ 6850284e-a7f1-4836-87bf-726d60e991b4
using LilGuys; FIGDIR = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, starsname, "figures")

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
	:xi_p_deg => L"\xi'\ /\ \textrm{deg}",
	:eta_p_deg => L"\eta'\ /\ \textrm{deg}",
	:xi_deg => L"\xi\ /\ \textrm{deg}",
	:eta_deg => L"\eta\ /\ \textrm{deg}",
	:xi => L"\xi\ /\ \textrm{arcmin}",
	:eta => L"\eta\ /\ \textrm{arcmin}",

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
	:eta_p => 300,
	:eta => 300,
	:eta_p_deg => 5,
	:eta_deg => 5,
)

# ╔═╡ 6cbf8dcd-d9e3-49e6-b799-a8bba9bf9211
import StatsBase: sample, weights

# ╔═╡ 0ccda86e-7a70-423c-a873-35b7a347bdfa
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ b37ee464-cd88-4722-a1c5-90585693476a
Arya.update_figsize!(6, 4)

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

	stars[!, :xi_deg] = stars.xi 
	stars[!, :eta_deg] = stars.eta
	stars.xi .*= 60
	stars.eta .*= 60

	stars[!, :xi_p_deg] = stars.xi_p ./60
	stars[!, :eta_p_deg] = stars.eta_p ./60

	return stars
end

# ╔═╡ e9b30734-0237-4b52-9f3f-add7c5460743
function sample_stars(stars, N=N)
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

# ╔═╡ a1cf1621-2076-4557-96d7-34adf19090d2
r_b = get_r_b(galaxyname, modelname, starsname)

# ╔═╡ 04c37154-fb06-47a8-b2ae-aff5eb7cde26
R_h = get_R_h(galaxyname)


# ╔═╡ 8bf905de-5457-4611-8b9d-3ca6039cc582
md"""
# tidal tails
"""

# ╔═╡ 48c7887a-d8b1-4579-b7c1-fa6ca27c74a1
stars = load_stars(galaxyname, modelname, starsname)

# ╔═╡ 04e2d6ef-3e2d-4cd2-b596-cf0366f71033
stars_sampled = sample_stars(stars)

# ╔═╡ 81638902-86b7-4e3f-b808-d5fd02751e9a
stars_0 = stars[1, :]

# ╔═╡ 6be80d25-974f-4bfe-96df-18cb0ce96c5a
function plot_sample(stars, stars_sampled)

	ss = stars_sampled
	fig = Figure()

	for i in 1:4
		ysym = [:pmra_gsr, :pmdec_gsr, :distance, :radial_velocity_gsr][i]

		gs = fig[(i-1) ÷ 2 + 1, (i-1)%2 + 1]
		ax = scatter_coord(gs, stars_0, ss, ysym, markersize=2, alpha=1, color=COLORS[3])
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

	Label(fig[0, :], text=galaxyname, font=:bold)
	fig

end

# ╔═╡ df04c1ab-5ea5-469b-b96f-e0a0ba542326
function plot_sample_heliocen(stars, ss)

	fig = Figure()

	for i in 1:4
		ysym = [:pmra, :pmdec, :distance, :radial_velocity][i]

		gs = fig[(i-1) ÷ 2 + 1, (i-1)%2 + 1]
		ax = scatter_coord(gs, stars_0, ss, ysym, markersize=2, alpha=1, color=COLORS[3])
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

	Label(fig[0, :], text=galaxyname, font=:bold)
	fig

end

# ╔═╡ 57284b92-0c3c-4e31-82af-65071ea9a3a2
let
	fig = Figure()
	ax = scatter_coord(fig[1, 1], stars_0, stars_sampled, :eta_p, markersize=2, alpha=1, color=COLORS[3])
	
	ax.aspect[] = DataAspect()
	
	ylims!(-300, 300)
	xlims!(-600, 600)

	fig
end

# ╔═╡ f4e74d3d-24b1-417e-82ea-d290df29dbaa
let
	fig = Figure()
	ax = scatter_coord(fig[1, 1], stars_0, stars_sampled, :eta_p_deg, xsym=:xi_p_deg, markersize=2, alpha=1, color=COLORS[3])
	
	ax.aspect[] = DataAspect()
	
	ylims!(-20, 20)
	xlims!(-60, 60)

	fig
end

# ╔═╡ b7758ae1-1f82-4f7f-90b6-f29f6f91c834
let
	fig = Figure()
	ax = scatter_coord(fig[1, 1], stars_0, stars_sampled, :eta_deg, xsym=:xi_deg, markersize=2, alpha=1, color=COLORS[3])
	
	ax.aspect[] = DataAspect()
	
	ylims!(-20, 20)
	xlims!(-20, 20)

	ax.xreversed[] = true

	fig
end

# ╔═╡ c8b2feb7-72ef-4a36-ad34-eba45e3564b5
stars_0[:distance]

# ╔═╡ 9ae7b57f-7649-443e-a438-4e5bd4b4d048
plot_sample(stars, stars_sampled)

# ╔═╡ a9583b46-3314-4aee-a1ea-ed7f433ae18d
DataFrames.median.(eachcol(stars))

# ╔═╡ b48832de-b20f-41ac-a652-145f2672d49a
plot_sample_heliocen(stars, stars_sampled)

# ╔═╡ Cell order:
# ╠═d593ea8b-47ec-4eca-a21b-ba07a118f468
# ╠═69c1ad4b-d6b0-4c0e-8623-9622a3779ec7
# ╠═5dff51c5-2c62-4469-8744-ea1a0a7c932b
# ╠═689209ef-3584-486c-9010-bf9bb59238f6
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═6850284e-a7f1-4836-87bf-726d60e991b4
# ╠═daf8bd82-89c1-496f-98ed-333b740095d9
# ╠═7dd7aa41-6fa4-49e0-a3d8-43768a4d0446
# ╠═8f0d73ac-58fe-4814-82ce-3295abaaf8bd
# ╠═cb0f69df-e8c6-4e24-a79c-8e2cccdc3718
# ╠═05830907-1355-49ee-8cd3-df4678f33149
# ╠═6cbf8dcd-d9e3-49e6-b799-a8bba9bf9211
# ╠═0ccda86e-7a70-423c-a873-35b7a347bdfa
# ╠═7128f9f0-36fd-4ffb-9751-7c66c8cfd8e5
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═b37ee464-cd88-4722-a1c5-90585693476a
# ╠═2c2b0102-43e6-4314-a6e4-855980f27149
# ╠═bcf2d214-84eb-42b4-95b0-da22ed4aecab
# ╠═e9b30734-0237-4b52-9f3f-add7c5460743
# ╠═904efe1c-6c2f-4e16-8468-ab4c5c79ebe1
# ╠═226173f5-93c4-401e-be28-d509880aff00
# ╠═f2ec2966-9018-4ad0-8ec3-36e438012ca6
# ╠═b257827b-ca9c-4abe-bbaa-1097f0034b03
# ╠═a1cf1621-2076-4557-96d7-34adf19090d2
# ╠═04c37154-fb06-47a8-b2ae-aff5eb7cde26
# ╠═6be80d25-974f-4bfe-96df-18cb0ce96c5a
# ╠═df04c1ab-5ea5-469b-b96f-e0a0ba542326
# ╟─8bf905de-5457-4611-8b9d-3ca6039cc582
# ╠═48c7887a-d8b1-4579-b7c1-fa6ca27c74a1
# ╠═04e2d6ef-3e2d-4cd2-b596-cf0366f71033
# ╠═81638902-86b7-4e3f-b808-d5fd02751e9a
# ╠═57284b92-0c3c-4e31-82af-65071ea9a3a2
# ╠═f4e74d3d-24b1-417e-82ea-d290df29dbaa
# ╠═b7758ae1-1f82-4f7f-90b6-f29f6f91c834
# ╠═c8b2feb7-72ef-4a36-ad34-eba45e3564b5
# ╠═9ae7b57f-7649-443e-a438-4e5bd4b4d048
# ╠═a9583b46-3314-4aee-a1ea-ed7f433ae18d
# ╠═b48832de-b20f-41ac-a652-145f2672d49a
