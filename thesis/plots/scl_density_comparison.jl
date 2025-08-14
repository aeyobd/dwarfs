### A Pluto.jl notebook ###
# v0.20.13

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

# ╔═╡ 3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:svg)

# ╔═╡ 7254f64c-b216-4720-82c2-85733c2757a7
R_h = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/observed_properties.toml"))["R_h_inner"]

# ╔═╡ b96ce635-f740-46e7-a0a2-0a4aa6a28e27
Mstar = 3e6

# ╔═╡ 65d3653d-5734-40ec-a057-aaa1e509968a
prof_expected = let
	prof = SurfaceDensityProfile(ENV["DWARFS_ROOT"] * "/observations/sculptor/density_profiles/jax_2c_eqw_profile.toml")
		
	prof =  LilGuys.filter_empty_bins(prof)

	prof = LilGuys.scale(prof, 1/R_h, 1)
	prof
end

# ╔═╡ 8d90a98b-09df-4834-ab9d-0aed44d6206b
function load_profile(haloname, orbitname, starsname)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/sculptor/$haloname/$orbitname/stars/$starsname/")

	prof_i = SurfaceDensityProfile(model_dir * "initial_profile.toml")
	prof_f = SurfaceDensityProfile(model_dir * "final_profile.toml")

	prof_i = LilGuys.scale(prof_i, 1/R_h, 1)
	prof_f = LilGuys.scale(prof_f, 1/R_h, 1)
	return prof_i,  prof_f
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

# ╔═╡ 214b6b76-5cc4-44f7-b365-0db9d7a890a8
function compare_profiles(prof_i, prof, r_b; kwargs...)
	fig = Figure()
	compare_profiles(fig[1,1], prof_i, prof, r_b; kwargs...)
	return fig
end

# ╔═╡ 89cc7dd4-6643-463d-9fa3-bfc2859476f2
function compare_profiles(gs, prof_i, prof, r_b; break_height=0, norm_shift = 0, r_j=nothing, plot_final=true) 
	ax = Axis(gs, 
		xlabel=L"log Radius / $R_h$",
		ylabel = "log surface density",
		limits=((-1.3 , 1.5), (-6, 4)),
	)




	dy = get_normalization(prof) + norm_shift .- get_normalization(prof_expected)
	
	@info dy - get_normalization(prof_expected)
	
	lines!(prof_i.log_R, prof_i.log_Sigma .+ dy, color=COLORS[2], linestyle=:dot,
			label="initial")

	if plot_final
		lines!(prof.log_R, prof.log_Sigma .+ dy, color=COLORS[2], linestyle=:solid,
				label="final")
	end
	errorscatter!(prof_expected.log_R, prof_expected.log_Sigma,
				  yerror = error_interval.(prof_expected.log_Sigma),
		label="observed",
		color=:black
	)
	
	if !isnothing(break_height)
		annotation!(0, 30, log10(r_b), break_height, 
					color=COLORS[2], 
					linewidth=theme(:linewidth)[]/2,  
					style=Ann.Styles.LineArrow(head = Ann.Arrows.Head()),
				   )
		
		text!(log10(r_b), break_height, 
			  text="break", 
			  rotation=π/2, 
			  offset=(0., 30.), 
			  color=COLORS[2], 
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
			  text="jacobi", 
			  rotation=π/2, 
			  offset=(0., 30.), 
			  align=(:left, :center),
			  fontsize=0.8 * theme(:fontsize)[]
			 )

	end

	axislegend(position=:lb)
end

# ╔═╡ 5abb6cec-e947-4fc1-9848-760e50bd5628
function get_stars_final(haloname, orbitname, starsname)
	modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis/sculptor/$haloname/$orbitname/stars/$starsname/")

	return read_fits(joinpath(modeldir, "final.fits"))
end

# ╔═╡ 05d4f257-9f20-4f5b-93fa-be4742b5060d
s = get_stars_final("1e7_V31_r3.2", "orbit_smallperi", "exp2d_rs0.10",)

# ╔═╡ c59d77a6-d5e1-4d66-87a9-d2f45192c699
bins = 4 * LinRange(-60, 60, 100)

# ╔═╡ b15244b0-b395-4cb4-a30e-4ad98cd10d5f
function plot_stars_2d(gs, modelname, haloname, starsname; r_b=nothing)
	stars = get_stars_final(modelname, haloname, starsname)


	ax = Axis(gs, xlabel=L"\xi \, / \, \textrm{arcmin}", 
			  ylabel=L"\eta \, / \, \textrm{arcmin}"
			)

	
	hist2d!(stars.xi*60, stars.eta*60, bins=bins, weights=stars.weights * Mstar, colorscale=log10, colorrange=(0.1, nothing))


	if !isnothing(r_b)
		arc!((0, 0), r_b, 0, 2π, color=COLORS[2])
	end
		
	ax
end

# ╔═╡ 9a13a8fc-efa0-4c7d-a232-5ef23f11bb4d
hist2d(s.xi*60, s.eta*60, bins=bins, weights=s.weights, colorscale=log10, colorrange=(1e-10, nothing))

# ╔═╡ 33017dbe-055b-440d-b43d-63a2c5e55441
let 
	fig = Figure()
	plot_stars_2d(fig[1,1], "1e6_new_v31_r3.2", "orbit_smallperi", "exp2d_rs0.13", r_b=60)
	fig
end

# ╔═╡ 6e17334f-86a8-470b-a9aa-30c382b0e5ae
function compare_profiles(halo::String, orbit::String, star; lmc=false, r_j=nothing, kwargs...)
	prof_i, prof_f = load_profile(halo, orbit, star)
	r_b = get_r_b(halo, orbit, star, lmc=lmc)
	if !isnothing(r_j)
		r_j = get_r_j(halo, orbit, lmc=lmc)
	end
	compare_profiles(prof_i, prof_f, r_b; r_j=r_j, kwargs...)
end

# ╔═╡ 24d2dade-4969-4a99-86aa-1e9ee0f64af6
function compare_both(halo::String, orbit::String, star; lmc=false, r_j=nothing, kwargs...)
	prof_i, prof_f = load_profile(halo, orbit, star)
	r_b = get_r_b(halo, orbit, star, lmc=lmc)
	if !isnothing(r_j)
		r_j = get_r_j(halo, orbit, lmc=lmc)
	end

	fig = Figure()

	plot_stars_2d(fig[1,1], halo, orbit, star, r_b=R_h * r_b)
	compare_profiles(fig[1,2], prof_i, prof_f, r_b; r_j=r_j, kwargs...)

	rowsize!(fig.layout, 1, Aspect(1, 1.0))
	resize_to_layout!(fig)
	fig
end

# ╔═╡ e231dce2-1969-4e86-a519-e9e2eb91d9f6
compare_both("1e6_new_v31_r3.2", "orbit_smallperi", "exp2d_rs0.13", norm_shift=0.2, 	break_height=-0)

# ╔═╡ 3317cfec-be70-418c-944a-0a6741f03cf3
@savefig "scl_smallperi_i_f" compare_both("1e7_V31_r3.2", "orbit_smallperi", "exp2d_rs0.10",  r_j=true, norm_shift=0.2, 	break_height=-0)

# ╔═╡ f053b7f9-f934-4cd4-ae46-32a0a373ff7e
# @savefig "scl_mean_i_f" compare_profiles("1e7_V31_r3.2", "orbit_mean", "exp2d_rs0.10", norm_shift=0.15, break_height=-5, r_j=get_r_j("1e7_V31_r3.2", "orbit_mean"))

# ╔═╡ fc45829c-bf28-4d7b-9bea-e223e4b29e7d
@savefig "scl_plummer_i_f" compare_both("1e7_V31_r3.2", "orbit_smallperi", "plummer_rs0.20", norm_shift=0.0,  r_j=true)

# ╔═╡ 5f989e91-8870-462e-a363-767be64dac5c
@savefig "scl_lmc_i_f" compare_both("1e7_V31_r4.2", "vasiliev24_L3M11_2x_smallperilmc", "exp2d_rs0.10", 
	norm_shift=0.1, lmc=true, r_j=true, break_height=-6)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
# ╠═9309c10c-6ba3-436c-b975-36d26dafb821
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═7254f64c-b216-4720-82c2-85733c2757a7
# ╠═b96ce635-f740-46e7-a0a2-0a4aa6a28e27
# ╠═65d3653d-5734-40ec-a057-aaa1e509968a
# ╠═8d90a98b-09df-4834-ab9d-0aed44d6206b
# ╠═b9e03109-9f49-4897-b26c-e31698a5fe49
# ╠═fe3bc6ee-14ed-4006-b3ec-f068d2492da4
# ╠═43c973e6-ad30-4bf8-93b3-4924bc498927
# ╠═f875e0e9-2f05-4420-8f40-76284be58e03
# ╠═214b6b76-5cc4-44f7-b365-0db9d7a890a8
# ╠═89cc7dd4-6643-463d-9fa3-bfc2859476f2
# ╠═5abb6cec-e947-4fc1-9848-760e50bd5628
# ╠═05d4f257-9f20-4f5b-93fa-be4742b5060d
# ╠═c59d77a6-d5e1-4d66-87a9-d2f45192c699
# ╠═b15244b0-b395-4cb4-a30e-4ad98cd10d5f
# ╠═9a13a8fc-efa0-4c7d-a232-5ef23f11bb4d
# ╠═33017dbe-055b-440d-b43d-63a2c5e55441
# ╠═6e17334f-86a8-470b-a9aa-30c382b0e5ae
# ╠═24d2dade-4969-4a99-86aa-1e9ee0f64af6
# ╠═e231dce2-1969-4e86-a519-e9e2eb91d9f6
# ╠═3317cfec-be70-418c-944a-0a6741f03cf3
# ╠═f053b7f9-f934-4cd4-ae46-32a0a373ff7e
# ╠═fc45829c-bf28-4d7b-9bea-e223e4b29e7d
# ╠═5f989e91-8870-462e-a363-767be64dac5c
