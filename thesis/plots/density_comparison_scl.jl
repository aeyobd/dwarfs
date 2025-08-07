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

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:svg)

# ╔═╡ 7254f64c-b216-4720-82c2-85733c2757a7
R_h = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/observed_properties.toml"))["R_h_inner"]

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

# ╔═╡ 991fe214-c10c-44a6-a284-0356ef993412
get_r_j("1e7_V31_r3.2", "orbit_smallperi", )

# ╔═╡ 43c973e6-ad30-4bf8-93b3-4924bc498927
readdir("/cosma/home/durham/dc-boye1/data/dwarfs/analysis/sculptor/1e7_V31_r3.2/orbit_mean/stars/")

# ╔═╡ f875e0e9-2f05-4420-8f40-76284be58e03
function get_normalization(prof)
	return -LilGuys.mean(prof.log_Sigma[1:5])
end

# ╔═╡ 56da34d4-7ccd-41f0-98b4-31bb1afc8a49
theme(:size)

# ╔═╡ d0140b51-02ca-4f25-9ecc-57a20f6a64cd
log10(R_h)

# ╔═╡ 89cc7dd4-6643-463d-9fa3-bfc2859476f2
function compare_profiles(prof_i, prof, r_b; break_height=-6, norm_shift = 0, r_j=nothing, plot_final=true) 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=L"log Radius / $R_h$",
		ylabel = "log surface density",
		limits=((-1.3 , 1.5), (-6, 4)),
	)

	errorscatter!(prof_expected.log_R, prof_expected.log_Sigma,
				  yerror = error_interval.(prof_expected.log_Sigma),
		label="observed",
		color=:black
	)


	dy = get_normalization(prof) + norm_shift .- get_normalization(prof_expected)
	
	@info dy - get_normalization(prof_expected)
	
	lines!(prof_i.log_R, prof_i.log_Sigma .+ dy, color=COLORS[2], linestyle=:dot,
			label="initial")

	if plot_final
		lines!(prof.log_R, prof.log_Sigma .+ dy, color=COLORS[2], linestyle=:solid,
				label="final")
	end
	
	if !isnothing(break_height)
		annotation!(0, 40, log10(r_b), break_height, color=COLORS[2], 
		text="break",  )
	end

	if !isnothing(r_j)
		annotation!(0, 40, log10(r_j), break_height,color=COLORS[1],
		text="jacobi", )
	end

	axislegend(position=:lb, margin=theme(:Legend).padding, patchsize=(36, 24/2))
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

# ╔═╡ 21150e87-f3dd-4965-b4c8-1cca022db114
@savefig "scl_smallperi_i" compare_profiles("1e7_V31_r3.2", "orbit_smallperi", "exp2d_rs0.10", norm_shift=0.15, break_height=nothing, plot_final=false)

# ╔═╡ 7d7ba92c-bb6e-471e-85b5-fbfcbdeee664
@savefig "scl_smallperi_i_f_b" compare_profiles("1e7_V31_r3.2", "orbit_smallperi", "exp2d_rs0.10", norm_shift=0.15, break_height=nothing)

# ╔═╡ 3317cfec-be70-418c-944a-0a6741f03cf3
@savefig "scl_smallperi_i_f" compare_profiles("1e7_V31_r3.2", "orbit_smallperi", "exp2d_rs0.10",  r_j=true, norm_shift=0.2)

# ╔═╡ f053b7f9-f934-4cd4-ae46-32a0a373ff7e
# @savefig "scl_mean_i_f" compare_profiles("1e7_V31_r3.2", "orbit_mean", "exp2d_rs0.10", norm_shift=0.15, break_height=-5, r_j=get_r_j("1e7_V31_r3.2", "orbit_mean"))

# ╔═╡ 4f02015a-dfc0-4c7c-935d-9992b0690576
@savefig "scl_plummer_i" compare_profiles("1e7_V31_r3.2", "orbit_smallperi", "plummer_rs0.20", norm_shift=0.0, plot_final=false, break_height=nothing)

# ╔═╡ fc45829c-bf28-4d7b-9bea-e223e4b29e7d
@savefig "scl_plummer_i_f" compare_profiles("1e7_V31_r3.2", "orbit_smallperi", "plummer_rs0.20", norm_shift=0.0, break_height=1.5, r_j=true)

# ╔═╡ 5f989e91-8870-462e-a363-767be64dac5c
@savefig "scl_lmc_i_f" compare_profiles("1e7_V31_r4.2", "vasiliev24_L3M11_2x_smallperilmc", "exp2d_rs0.10", 
	norm_shift=0.1, break_height=-3.5, lmc=true, r_j=true)

# ╔═╡ fab0c970-7d9a-4ed6-9632-72a7b0d52c5b
LilGuys.Plummer(r_s=0.3)

# ╔═╡ 084e8e1d-daec-4859-bbe2-888d484d6123
@savefig "scl_lmc_i_f_mwb" compare_profiles("1e7_V31_r4.2", "vasiliev24_L3M11_2x_smallperilmc", "exp2d_rs0.10", 
	norm_shift=0.1, break_height=-3.5, lmc=false, r_j=true)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═7254f64c-b216-4720-82c2-85733c2757a7
# ╠═65d3653d-5734-40ec-a057-aaa1e509968a
# ╠═8d90a98b-09df-4834-ab9d-0aed44d6206b
# ╠═b9e03109-9f49-4897-b26c-e31698a5fe49
# ╠═fe3bc6ee-14ed-4006-b3ec-f068d2492da4
# ╠═991fe214-c10c-44a6-a284-0356ef993412
# ╠═43c973e6-ad30-4bf8-93b3-4924bc498927
# ╠═f875e0e9-2f05-4420-8f40-76284be58e03
# ╠═56da34d4-7ccd-41f0-98b4-31bb1afc8a49
# ╠═d0140b51-02ca-4f25-9ecc-57a20f6a64cd
# ╠═89cc7dd4-6643-463d-9fa3-bfc2859476f2
# ╠═6e17334f-86a8-470b-a9aa-30c382b0e5ae
# ╠═21150e87-f3dd-4965-b4c8-1cca022db114
# ╠═7d7ba92c-bb6e-471e-85b5-fbfcbdeee664
# ╠═3317cfec-be70-418c-944a-0a6741f03cf3
# ╠═f053b7f9-f934-4cd4-ae46-32a0a373ff7e
# ╠═4f02015a-dfc0-4c7c-935d-9992b0690576
# ╠═fc45829c-bf28-4d7b-9bea-e223e4b29e7d
# ╠═5f989e91-8870-462e-a363-767be64dac5c
# ╠═fab0c970-7d9a-4ed6-9632-72a7b0d52c5b
# ╠═084e8e1d-daec-4859-bbe2-888d484d6123
