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

# ╔═╡ e084a02a-f445-422f-b0b3-700a44bf204c
function get_R_h(galaxyname)

	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxyname/observed_properties.toml"))
	R_h = obs_props["R_h"]
end

# ╔═╡ 65d3653d-5734-40ec-a057-aaa1e509968a
function load_expected(galaxyname)
	prof = SurfaceDensityProfile(ENV["DWARFS_ROOT"] * "/observations/$galaxyname/density_profiles/fiducial_profile.toml")
		
	prof =  LilGuys.filter_empty_bins(prof)

	R_h = get_R_h(galaxyname)
	prof = LilGuys.scale(prof, 1/R_h, 1)
	prof
end

# ╔═╡ b9e03109-9f49-4897-b26c-e31698a5fe49
function get_r_b(galaxyname, haloname, orbitname, starsname; lmc=false)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$haloname/$orbitname/stars/$starsname/")

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
	R_h = get_R_h(galaxyname)

	return LilGuys.kpc2arcmin(r_b, dist_f)	/ R_h
end

# ╔═╡ 8404ee2c-a99d-4dcc-a468-2629f9a17abf


# ╔═╡ fe3bc6ee-14ed-4006-b3ec-f068d2492da4
function get_r_j(galaxyname, haloname, orbitname; lmc=false)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$haloname/$orbitname/")

	if lmc
		props = TOML.parsefile(model_dir * "jacobi_lmc.toml")
	else
		props = TOML.parsefile(model_dir * "jacobi.toml")
	end
	R_h = get_R_h(galaxyname)

	return props["r_J"] / R_h
end

# ╔═╡ f875e0e9-2f05-4420-8f40-76284be58e03
function get_normalization(prof)
	return -LilGuys.mean(prof.log_Sigma[1:5])
end

# ╔═╡ 8d90a98b-09df-4834-ab9d-0aed44d6206b
function load_profile(galaxyname, haloname, orbitname, starsname; norm_shift=0)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$haloname/$orbitname/stars/$starsname/")

	prof_i = SurfaceDensityProfile(model_dir * "initial_profile.toml")
	prof_f = SurfaceDensityProfile(model_dir * "final_profile.toml")

	R_h = get_R_h(galaxyname)

	prof_i = LilGuys.scale(prof_i, 1/R_h, 1)
	prof_f = LilGuys.scale(prof_f, 1/R_h, 1)

	prof_expected = load_expected(galaxyname)

	dy = get_normalization(prof_f) + norm_shift .- get_normalization(prof_expected)

	dy = middle(dy)
	prof_i = LilGuys.scale(prof_i, 1, 10^dy)
	prof_f = LilGuys.scale(prof_f, 1, 10^dy)
	
	return prof_i,  prof_f, dy
end

# ╔═╡ 214b6b76-5cc4-44f7-b365-0db9d7a890a8
function compare_profiles(prof_i, prof, r_b; kwargs...)
	fig = Figure()
	compare_profiles(fig[1,1], prof_i, prof, r_b; kwargs...)
	return fig
end

# ╔═╡ 5abb6cec-e947-4fc1-9848-760e50bd5628
function get_stars_final(galaxyname, haloname, orbitname, starsname, filename="final.fits")
	modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$haloname/$orbitname/stars/$starsname/")

	return read_fits(joinpath(modeldir, filename))
end

# ╔═╡ c59d77a6-d5e1-4d66-87a9-d2f45192c699
bins = 4 * LinRange(-60, 60, 100)

# ╔═╡ 9a0e561a-7ba5-442c-95aa-37ff5a469f4c
logdensityrange = (-6., 4.)

# ╔═╡ 89cc7dd4-6643-463d-9fa3-bfc2859476f2
function compare_profiles(gs, prof_i, prof, r_b; galaxyname, t_i, break_height=0, r_j=nothing, plot_final=true) 
	ax = Axis(gs, 
		xlabel=L"log Radius / $R_h$",
		ylabel = "log surface density",
		limits=((-1.3 , 1.5), logdensityrange),
	)
	
	lines!(prof_i.log_R, prof_i.log_Sigma, color=COLORS[3], linestyle=:dot,
			label=L"initial ($t = %$(round(t_i, digits=1))$\,Gyr)")

	if plot_final
		lines!(prof.log_R, prof.log_Sigma, color=COLORS[3], linestyle=:solid,
				label=L"final ($t = 0.0\,$Gyr)")
	end
	prof_expected = load_expected(galaxyname)

	errorscatter!(prof_expected.log_R, prof_expected.log_Sigma,
				  yerror = error_interval.(prof_expected.log_Sigma),
		label="observed",
		color=:black
	)
	
	if !isnothing(break_height)
		annotation!(0, 30, log10(r_b), break_height, 
					color=COLORS[3], 
					linewidth=theme(:linewidth)[]/2,  
					style=Ann.Styles.LineArrow(head = Ann.Arrows.Head()),
				   )
		
		text!(log10(r_b), break_height, 
			  text="break", 
			  rotation=π/2, 
			  offset=(0., 30.), 
			  color=COLORS[3], 
			  align=(:left, :center), 
			  fontsize=0.8 * theme(:fontsize)[]
			)
	end

	if !isnothing(r_j)
		annotation!(0, 30, log10(r_j), break_height,
					linewidth=theme(:linewidth)[]/2,
					style=Ann.Styles.LineArrow(head = Ann.Arrows.Head()),
				   )

		text!(log10(r_j), break_height, 
			  text="Jacobi", 
			  rotation=π/2, 
			  offset=(0., 30.), 
			  align=(:left, :center),
			  fontsize=0.8 * theme(:fontsize)[]
			 )

	end

	axislegend(position=:lb)
end

# ╔═╡ 3ed68a46-04aa-4fea-b484-d9f3cc158738
smallfontsize=0.8*theme(:fontsize)[]

# ╔═╡ cb898eeb-0803-42a7-a2c9-d6e2b95f8945
smalllinewidth = theme(:linewidth)[]/2

# ╔═╡ 6e17334f-86a8-470b-a9aa-30c382b0e5ae
function compare_profiles(halo::String, orbit::String, star; lmc=false, norm_shift=0, r_j=nothing, kwargs...)
	prof_i, prof_f, dy= load_profile(halo, orbit, star; norm_shift=norm_shift)
	r_b = get_r_b(halo, orbit, star, lmc=lmc)
	if !isnothing(r_j)
		r_j = get_r_j(halo, orbit, lmc=lmc)
	end
	compare_profiles(prof_i, prof_f, r_b; galaxyname=galaxyname, r_j=r_j, kwargs...)
end

# ╔═╡ 647b5904-22a5-41d2-98a7-60a75597983d
function get_time_ini(galaxyname, haloname, orbitname)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$haloname/$orbitname")
	
	TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))["t_f_gyr"] 
end


# ╔═╡ 32f4d363-8115-4819-a2ee-4c25d792a14e
function get_recent_orbit(galaxyname, haloname, orbitname)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/$galaxyname/$haloname/$orbitname")
	orbit = LilGuys.Orbit(joinpath(model_dir, "centres.hdf5"))

	idx_f = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))["idx_f"]
	gcs = Galactocentric.(orbit.positions[1, :], orbit.positions[2, :], orbit.positions[3, :], V2KMS*orbit.velocities[1, :], V2KMS*orbit.velocities[2, :], V2KMS*orbit.velocities[3, :])
	icrs = LilGuys.transform.(ICRS, gcs)


	gsr = LilGuys.transform(GSR, icrs[idx_f])
	dx, dy = (gsr.pmra, gsr.pmdec) ./ sqrt(gsr.pmra^2 + gsr.pmdec^2)
	
	R_h = get_R_h(galaxyname)


	xi_past = [-1e3 * dx, -6*R_h * dx]
	eta_past = [-1e3 * dy, -6*R_h * dy]
	xi_future = [6*R_h * dx, 1e3*dx]
	eta_future = [6*R_h * dy,  1e3*dy]
	return xi_past, eta_past, xi_future, eta_future, -T2GYR*(orbit.times[idx_f] - orbit.times[1])
end

# ╔═╡ 7e5b185c-77b4-4066-82f3-3d9f1393db9e
get_recent_orbit("sculptor", "1e6_new_v31_r3.2", "orbit_smallperi",)

# ╔═╡ ece6315f-1dc3-4874-8760-4b32474366f3
colormap = Reverse(:Greys)

# ╔═╡ b15244b0-b395-4cb4-a30e-4ad98cd10d5f
function plot_stars_2d(gs, galaxyname, modelname, haloname, starsname; initial=false, norm, r_b=nothing)
	if initial
		filename = "initial.fits"
	else
		filename = "final.fits"
	end
	stars = get_stars_final(galaxyname, modelname, haloname, starsname, filename)


	ax = Axis(gs, 
			  xlabel = L"\xi \, / \, \textrm{arcmin}", 
			  ylabel = L"\eta \, / \, \textrm{arcmin}",
			  xreversed = true,
			)

	
	h = histogram2d(stars.xi*60, stars.eta*60, bins, weights=stars.weights * 10^norm)

	@info "density maximum: $(maximum(h.values))"
	p = heatmap!(h.xbins, h.ybins, log10.(h.values), colorrange=logdensityrange, colormap=colormap)


	if !isnothing(r_b)
		arc!((0, 0), r_b, 0, 2π, color=COLORS[3], linestyle=:dash, linewidth=smalllinewidth)
		text!(r_b, 0, text="break", color=COLORS[3], fontsize=smallfontsize, offset=(-smallfontsize/2, 0), align=(:center, :bottom), rotation=π/2)

	end

	R_h = get_R_h(galaxyname)
	arc!((0,0), 6R_h, 0, 2π, color=:white, linewidth=theme(:linewidth)[]/2)
	if initial
		text!(6R_h, 0, offset=(smallfontsize/2, 0), text=L"6R_h", align=(:left, :center), color=:white, smallfontsize)
	end
	p
end

# ╔═╡ 24d2dade-4969-4a99-86aa-1e9ee0f64af6
function compare_both(galaxyname, halo::String, orbit::String, star; norm_shift=0, lmc=false, r_j=nothing, title="", kwargs...)
	prof_i, prof_f, norm = load_profile(galaxyname, halo, orbit, star, norm_shift=norm_shift)
	r_b = get_r_b(galaxyname, halo, orbit, star, lmc=lmc)
	if !isnothing(r_j)
		r_j = get_r_j(galaxyname, halo, orbit, lmc=lmc)
	end

	fig = Figure()
	R_h = get_R_h(galaxyname)

	p = plot_stars_2d(fig[1,2], galaxyname, halo, orbit, star, r_b=R_h * r_b, norm=norm)
	hideydecorations!()
	xi, eta, xi_next, eta_next, t_i = get_recent_orbit(galaxyname, halo, orbit) 
	lines!(xi, eta, color=(COLORS[1]), linewidth=smalllinewidth)
	lines!(xi_next, eta_next, color=(COLORS[1]), linestyle=:dot, linewidth=smalllinewidth)

	xlims!(reverse(extrema(bins))...)
	ylims!(extrema(bins)...)
	p = plot_stars_2d(fig[1,1], galaxyname, halo, orbit, star, initial=true,norm=norm)

	Colorbar(fig[1,3], p, label="log surface density", ticks=Makie.automatic)
	compare_profiles(fig[2,1:3], prof_i, prof_f, r_b; galaxyname=galaxyname, r_j=r_j, t_i=t_i, kwargs...)


	Makie.Label(fig[0, :], title)
	rowsize!(fig.layout, 1, Aspect(1, 1.0))
	rowsize!(fig.layout, 2, Aspect(1, 1.5))
	resize_to_layout!(fig)
	fig
end

# ╔═╡ e231dce2-1969-4e86-a519-e9e2eb91d9f6
compare_both("sculptor", "1e6_new_v31_r3.2", "orbit_smallperi", "exp2d_rs0.10", norm_shift=0.0, 	break_height=-0, title="Sculptor low-resolution")

# ╔═╡ fafd07cb-2378-421b-b6e5-8525c71f9c8d
compare_both("sculptor", "1e6_new_v43_r7", "orbit_smallperi.3", "exp2d_rs0.10", norm_shift=-0.7, 	break_height=-0, title="Sculptor heavy")

# ╔═╡ 3317cfec-be70-418c-944a-0a6741f03cf3
@savefig "scl_smallperi_i_f" compare_both("sculptor", "1e7_new_v31_r3.2", "orbit_smallperi", "exp2d_rs0.10",  r_j=true, norm_shift=0.2, 	break_height=-0, title="Sculptor: smallperi-exponential")

# ╔═╡ fc45829c-bf28-4d7b-9bea-e223e4b29e7d
@savefig "scl_plummer_i_f" compare_both("sculptor", "1e7_new_v31_r3.2", "orbit_smallperi", "plummer_rs0.20", norm_shift=0.0,  r_j=true, title="Sculptor: smallperi-Plummer")

# ╔═╡ 5f989e91-8870-462e-a363-767be64dac5c
@savefig "scl_lmc_i_f" compare_both("sculptor", "1e7_new_v25_r2.5", "smallperilmc", "exp2d_rs0.11", norm_shift=0.1, lmc=true, r_j=true, break_height=-6, title="Sculptor: LMC–exponential")

# ╔═╡ e9a0ce9f-27a5-47dc-8c4d-22887f6a7fc1
@savefig "scl_lmc_plummer_i_f" compare_both("sculptor", "1e7_new_v25_r2.5", "smallperilmc", "plummer_rs0.20", norm_shift=0.1, lmc=true, r_j=true, break_height=-6, title="Sculptor: LMC–Plummer")

# ╔═╡ f11acd53-a04d-4e1a-b7f5-d66b6011e3a2
@savefig "scl_impact_i_f" compare_both("sculptor", "1e6_new_v31_r3.2", "L3M11_9Gyr_smallperi.a4", "exp2d_rs0.10",  title="Sculptor: MW impact")

# ╔═╡ f6eaabb0-8508-40fd-9dd1-37b1165e1c6d


# ╔═╡ c7b50c9b-1bfd-4abf-89cd-7c8cd1dc45c0
compare_both("sculptor", "1e6_v38_r7_oblate_0.5", "orbit_smallperi", "exp2d_rs0.13",  title="Sculptor: oblate")

# ╔═╡ 2a3b1df5-7a44-420c-b303-856346df40cf
compare_both("sculptor", "1e6_v43_r3_beta0.2_a4", "orbit_smallperi", "exp2d_rs0.13",  title="Sculptor: anisotropic")

# ╔═╡ 8e4c5116-7de6-4f8d-bc09-ff3aa4ffc50f
compare_both("sculptor", "1e4_exp_M3e-4_r0.1", "orbit_smallperi", "stars", title="Sculptor: dm free", norm_shift=-5)

# ╔═╡ 9c0f4e30-ee69-4b55-a15d-0566c5582200


# ╔═╡ 77c057e4-78da-4201-b9f8-cb3135508a93
md"""
# Ursa Minor
"""

# ╔═╡ 3263677b-e0e3-4d30-8e28-309de99a8603
@savefig "umi_smallperi_i_f" compare_both("ursa_minor", "1e7_new_v38_r4.0", "orbit_smallperi.5", "exp2d_rs0.10",  r_j=true, norm_shift=0.0, 	break_height=-0, title="Ursa Minor: smallperi-exponential")

# ╔═╡ d0ce116f-aa51-43d8-9748-806eb44b8409
@savefig "umi_plummer_i_f" compare_both("ursa_minor", "1e7_new_v38_r4.0", "orbit_smallperi.5", "plummer_rs0.20",  r_j=true, norm_shift=0.0,break_height=-0, title="Ursa Minor: smallperi-Plummer")

# ╔═╡ 5b0df0a1-4c82-45db-a42f-b5ed23aad258
compare_both("ursa_minor", "1e6_new_v38_r4.0", "orbit_smallperi.4", "exp2d_rs0.13", norm_shift=-0.5, 	break_height=-0, title="Ursa Minor")

# ╔═╡ 18f2e3da-b749-4134-9604-4e2318ad82ee
 compare_both("ursa_minor", "1e6_v37_r5.0", "orbit_mean.2", "exp2d_rs0.10",  r_j=true, norm_shift=-0.6, 	break_height=-0)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═2bacd818-4985-4922-85a3-716bdfda5146
# ╠═3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
# ╠═9309c10c-6ba3-436c-b975-36d26dafb821
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═e084a02a-f445-422f-b0b3-700a44bf204c
# ╠═65d3653d-5734-40ec-a057-aaa1e509968a
# ╠═8d90a98b-09df-4834-ab9d-0aed44d6206b
# ╠═b9e03109-9f49-4897-b26c-e31698a5fe49
# ╠═8404ee2c-a99d-4dcc-a468-2629f9a17abf
# ╠═fe3bc6ee-14ed-4006-b3ec-f068d2492da4
# ╠═f875e0e9-2f05-4420-8f40-76284be58e03
# ╠═214b6b76-5cc4-44f7-b365-0db9d7a890a8
# ╠═89cc7dd4-6643-463d-9fa3-bfc2859476f2
# ╠═5abb6cec-e947-4fc1-9848-760e50bd5628
# ╠═c59d77a6-d5e1-4d66-87a9-d2f45192c699
# ╠═9a0e561a-7ba5-442c-95aa-37ff5a469f4c
# ╠═b15244b0-b395-4cb4-a30e-4ad98cd10d5f
# ╠═3ed68a46-04aa-4fea-b484-d9f3cc158738
# ╠═cb898eeb-0803-42a7-a2c9-d6e2b95f8945
# ╠═6e17334f-86a8-470b-a9aa-30c382b0e5ae
# ╠═647b5904-22a5-41d2-98a7-60a75597983d
# ╠═32f4d363-8115-4819-a2ee-4c25d792a14e
# ╠═24d2dade-4969-4a99-86aa-1e9ee0f64af6
# ╠═7e5b185c-77b4-4066-82f3-3d9f1393db9e
# ╠═ece6315f-1dc3-4874-8760-4b32474366f3
# ╠═e231dce2-1969-4e86-a519-e9e2eb91d9f6
# ╠═fafd07cb-2378-421b-b6e5-8525c71f9c8d
# ╠═3317cfec-be70-418c-944a-0a6741f03cf3
# ╠═fc45829c-bf28-4d7b-9bea-e223e4b29e7d
# ╠═5f989e91-8870-462e-a363-767be64dac5c
# ╠═e9a0ce9f-27a5-47dc-8c4d-22887f6a7fc1
# ╠═f11acd53-a04d-4e1a-b7f5-d66b6011e3a2
# ╠═f6eaabb0-8508-40fd-9dd1-37b1165e1c6d
# ╠═c7b50c9b-1bfd-4abf-89cd-7c8cd1dc45c0
# ╠═2a3b1df5-7a44-420c-b303-856346df40cf
# ╠═8e4c5116-7de6-4f8d-bc09-ff3aa4ffc50f
# ╠═9c0f4e30-ee69-4b55-a15d-0566c5582200
# ╟─77c057e4-78da-4201-b9f8-cb3135508a93
# ╠═3263677b-e0e3-4d30-8e28-309de99a8603
# ╠═d0ce116f-aa51-43d8-9748-806eb44b8409
# ╠═5b0df0a1-4c82-45db-a42f-b5ed23aad258
# ╠═18f2e3da-b749-4134-9604-4e2318ad82ee
