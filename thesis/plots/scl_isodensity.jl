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

# ╔═╡ 9309c10c-6ba3-436c-b975-36d26dafb821
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 2bacd818-4985-4922-85a3-716bdfda5146
import DensityEstimators: histogram2d

# ╔═╡ 3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:svg)

# ╔═╡ b61e0dca-f6a6-498d-9af3-68080ee9eb62
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/observed_properties.toml"))

# ╔═╡ 7254f64c-b216-4720-82c2-85733c2757a7
R_h = obs_props["R_h_inner"]

# ╔═╡ b96ce635-f740-46e7-a0a2-0a4aa6a28e27
gc_obs = LilGuys.transform(Galactocentric, ICRS(obs_props))

# ╔═╡ 65d3653d-5734-40ec-a057-aaa1e509968a
prof_expected = let
	prof = SurfaceDensityProfile(ENV["DWARFS_ROOT"] * "/observations/sculptor/density_profiles/jax_2c_eqw_profile.toml")
		
	prof =  LilGuys.filter_empty_bins(prof)

	prof = LilGuys.scale(prof, 1/R_h, 1)
	prof
end

# ╔═╡ b9e03109-9f49-4897-b26c-e31698a5fe49
function get_r_b(haloname, orbitname, starsname; lmc=false)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/sculptor/$haloname/$orbitname/stars/$starsname/")

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
	return LilGuys.kpc2arcmin(r_b, dist_f)	/ R_h
end

# ╔═╡ 8404ee2c-a99d-4dcc-a468-2629f9a17abf


# ╔═╡ fe3bc6ee-14ed-4006-b3ec-f068d2492da4
function get_r_j(haloname, orbitname; lmc=false)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/sculptor/$haloname/$orbitname/")

	if lmc
		props = TOML.parsefile(model_dir * "jacobi_lmc.toml")
	else
		props = TOML.parsefile(model_dir * "jacobi.toml")
	end
	return props["r_J"] / R_h
end

# ╔═╡ 43c973e6-ad30-4bf8-93b3-4924bc498927
readdir("/cosma/home/durham/dc-boye1/data/dwarfs/analysis/sculptor/1e7_V31_r3.2/orbit_mean/stars/")

# ╔═╡ f875e0e9-2f05-4420-8f40-76284be58e03
function get_normalization(prof)
	return -LilGuys.mean(prof.log_Sigma[1:5])
end

# ╔═╡ 8d90a98b-09df-4834-ab9d-0aed44d6206b
function load_profile(haloname, orbitname, starsname; norm_shift=0)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/sculptor/$haloname/$orbitname/stars/$starsname/")

	prof_i = SurfaceDensityProfile(model_dir * "initial_profile.toml")
	prof_f = SurfaceDensityProfile(model_dir * "final_profile.toml")


	
	prof_i = LilGuys.scale(prof_i, 1/R_h, 1)
	prof_f = LilGuys.scale(prof_f, 1/R_h, 1)

	

	dy = get_normalization(prof_f) + norm_shift .- get_normalization(prof_expected)

	dy = middle(dy)
	prof_i = LilGuys.scale(prof_i, 1, 10^dy)
	prof_f = LilGuys.scale(prof_f, 1, 10^dy)
	
	return prof_i,  prof_f, dy
end

# ╔═╡ 5abb6cec-e947-4fc1-9848-760e50bd5628
function get_stars_final(haloname, orbitname, starsname, filename="final.fits")
	modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis/sculptor/$haloname/$orbitname/stars/$starsname/")

	return read_fits(joinpath(modeldir, filename))
end

# ╔═╡ 05d4f257-9f20-4f5b-93fa-be4742b5060d
s = get_stars_final("1e7_V31_r3.2", "orbit_smallperi", "exp2d_rs0.10",)

# ╔═╡ c59d77a6-d5e1-4d66-87a9-d2f45192c699
bins = 4 * LinRange(-60, 60, 100)

# ╔═╡ 9a0e561a-7ba5-442c-95aa-37ff5a469f4c
logdensityrange = (-6., 4.)

# ╔═╡ 3ed68a46-04aa-4fea-b484-d9f3cc158738
smallfontsize=0.8*theme(:fontsize)[]

# ╔═╡ cb898eeb-0803-42a7-a2c9-d6e2b95f8945
smalllinewidth = theme(:linewidth)[]/2

# ╔═╡ f8c8aff1-bbc2-4be8-9262-6fd25d49414c
import Contour

# ╔═╡ cf872696-9351-4f45-af6d-92ae8f80e3de
import DensityEstimators

# ╔═╡ 829f69d5-9a8f-4464-820e-d2a6d1a83063
function plot_stars_iso(gs, modelname, haloname, starsname; initial=false, norm, r_b=nothing, filename="final.fits")

	stars = get_stars_final(modelname, haloname, starsname, filename)

	@info "loaded stars"
	r_max = 300

	ax = Axis(gs, xlabel=L"\xi \, / \, \textrm{arcmin}", 
			  ylabel=L"\eta \, / \, \textrm{arcmin}",
			  aspect=DataAspect(),
			  limits=(-r_max, r_max, -r_max, r_max),
			  xreversed=true,
			)
	Nbins = 50

	bins = (LinRange(-r_max, r_max, Nbins), LinRange(-r_max, r_max, Nbins))
	
	h = #DensityEstimators.kde2d(stars.xi*60, stars.eta*60, 10, weights=stars.weights * 10^norm, bins=LinRange(-r_max, r_max, 100),)
	
	DensityEstimators.histogram2d(stars.xi*60, stars.eta*60, bins, weights=stars.weights * 10^norm,)
	
	colormax = maximum(h.values)

	@info "calculated kde"

	colorrange = (log10(colormax/1e10), log10(colormax))

	x = midpoints(h.xbins)
	y = midpoints(h.ybins)
	cs = Contour.contours(x, y, log10.(h.values .+ colormax / 1e20), LinRange(colorrange..., 21))

	@info "calculated contours"

	
	for cl in Contour.levels(cs)
		for (i, l) in enumerate(Contour.lines(cl))
			x, y = Contour.coordinates(l)

			lines!(x, y, color=Contour.level(cl), colorrange=colorrange, linestyle=:solid)

			if i == 1
				R_m = LilGuys.mean(radii([x y]'))
				arc!((0,0), R_m, 0, 2π, color=(:black, 0.5), linewidth=smalllinewidth/2)
			end
		end
	end

end

# ╔═╡ 90cd693b-6f1a-469a-b551-f9a2f584257b
let
	fig = Figure()
	cs = plot_stars_iso(fig[1,1], "1e6_new_v31_r3.2", "orbit_smallperi", "exp2d_rs0.13",  norm=0)

	fig

end

# ╔═╡ 3790574e-77a4-43dc-bc06-771726844351
let
	fig = Figure()

	plot_stars_iso(fig[1,1], "1e7_new_v25_r2.5", "smallperilmc", "plummer_rs0.20",  norm=0, filename="pre_lmc.fits")
	


	
	plot_stars_iso(fig[1,2], "1e7_new_v25_r2.5", "smallperilmc", "plummer_rs0.20",  norm=0)

	hideydecorations!(ticks=false, minorticks=false)

	
	fig
end

# ╔═╡ aaad9cde-b588-467d-8d8d-dd1f9ce48b23
let
	fig = Figure()
	plot_stars_iso(fig[1,1], "1e7_new_v25_r2.5", "smallperilmc", "exp2d_rs0.11",  norm=0)

	
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═2bacd818-4985-4922-85a3-716bdfda5146
# ╠═3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
# ╠═9309c10c-6ba3-436c-b975-36d26dafb821
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═b61e0dca-f6a6-498d-9af3-68080ee9eb62
# ╠═7254f64c-b216-4720-82c2-85733c2757a7
# ╠═b96ce635-f740-46e7-a0a2-0a4aa6a28e27
# ╠═65d3653d-5734-40ec-a057-aaa1e509968a
# ╠═8d90a98b-09df-4834-ab9d-0aed44d6206b
# ╠═b9e03109-9f49-4897-b26c-e31698a5fe49
# ╠═8404ee2c-a99d-4dcc-a468-2629f9a17abf
# ╠═fe3bc6ee-14ed-4006-b3ec-f068d2492da4
# ╠═43c973e6-ad30-4bf8-93b3-4924bc498927
# ╠═f875e0e9-2f05-4420-8f40-76284be58e03
# ╠═5abb6cec-e947-4fc1-9848-760e50bd5628
# ╠═05d4f257-9f20-4f5b-93fa-be4742b5060d
# ╠═c59d77a6-d5e1-4d66-87a9-d2f45192c699
# ╠═9a0e561a-7ba5-442c-95aa-37ff5a469f4c
# ╠═3ed68a46-04aa-4fea-b484-d9f3cc158738
# ╠═cb898eeb-0803-42a7-a2c9-d6e2b95f8945
# ╠═f8c8aff1-bbc2-4be8-9262-6fd25d49414c
# ╠═cf872696-9351-4f45-af6d-92ae8f80e3de
# ╠═829f69d5-9a8f-4464-820e-d2a6d1a83063
# ╠═90cd693b-6f1a-469a-b551-f9a2f584257b
# ╠═3790574e-77a4-43dc-bc06-771726844351
# ╠═aaad9cde-b588-467d-8d8d-dd1f9ce48b23
