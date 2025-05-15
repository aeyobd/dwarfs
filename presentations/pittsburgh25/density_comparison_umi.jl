### A Pluto.jl notebook ###
# v0.20.8

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
include("./style.jl")

# ╔═╡ b1002c1c-4f9e-4962-b5b4-ecab6ddf098c
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ 65d3653d-5734-40ec-a057-aaa1e509968a
prof_expected = SurfaceDensityProfile(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/density_profiles/jax_2c_eqw_profile.toml") |> LilGuys.filter_empty_bins

# ╔═╡ 8d90a98b-09df-4834-ab9d-0aed44d6206b
function load_profile(haloname, orbitname, starsname)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/ursa_minor/$haloname/$orbitname/stars/$starsname/")

	return SurfaceDensityProfile(model_dir * "initial_profile.toml"), SurfaceDensityProfile(model_dir * "final_profile.toml")
end

# ╔═╡ dd8a0445-8c94-4f43-b22a-1146a7f55d63
function get_r_b(haloname, orbitname, starsname)
	model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/ursa_minor/$haloname/$orbitname/stars/$starsname/")

	prof_f = SurfaceDensityProfile(model_dir * "final_profile.toml")

	σv = prof_f.annotations["sigma_v"]
	props = TOML.parsefile(model_dir * "../../orbital_properties.toml")
	dt = props["t_last_peri"]
	r_b = LilGuys.break_radius(σv / V2KMS, dt / T2GYR)
	return LilGuys.kpc2arcmin(r_b, props["distance_f"])	
end

# ╔═╡ 43c973e6-ad30-4bf8-93b3-4924bc498927
readdir("/cosma/home/durham/dc-boye1/data/dwarfs/analysis/sculptor/1e7_V31_r3.2/orbit_mean/stars/")

# ╔═╡ f875e0e9-2f05-4420-8f40-76284be58e03
function get_normalization(prof)
	return -LilGuys.mean(prof.log_Sigma[1:3])
end

# ╔═╡ 89cc7dd4-6643-463d-9fa3-bfc2859476f2
function compare_profiles(prof_i, prof, r_b) 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=L"\log\ R \ /\ \textrm{arcmin}",
		ylabel = L"\log \Sigma",
		limits=((-0.3, 2.5), (-9, 1))
	)

	errorscatter!(prof_expected.log_R, prof_expected.log_Sigma .+ get_normalization(prof_expected),
				  yerror = error_interval.(prof_expected.log_Sigma),
		label="observed",
		color=:black
	)


	dy = get_normalization(prof)
	
	lines!(prof_i.log_R, prof_i.log_Sigma .+ dy, color=COLORS[2], linestyle=:dot,
			label="initial")

	lines!(prof.log_R, prof.log_Sigma .+ dy, color=COLORS[2], linestyle=:solid,
			label="final")

	vlines!(log10(r_b), color=:grey, label="break radius", linewidth=3, linestyle=:dot)

	axislegend(position=:lb, margin=theme(:Legend).padding, patchsize=(48*1.5, 24))
	fig
end

# ╔═╡ 4851087d-2679-4409-b44f-cb22f778974e
function compare_profiles(halo::String, orbit::String, star)
	prof_i, prof_f = load_profile(halo, orbit, star)
	r_b = get_r_b(halo, orbit, star)
	compare_profiles(prof_i, prof_f, r_b)
end

# ╔═╡ 4e76a22b-631a-4478-b3bd-47be18496b3a
compare_profiles("1e6_v37_r5.0", "orbit_mean.2", "exp2d_rs0.15")

# ╔═╡ 3317cfec-be70-418c-944a-0a6741f03cf3
compare_profiles("1e6_v38_r4.0", "orbit_smallperi.4", "exp2d_rs0.10")

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═b1002c1c-4f9e-4962-b5b4-ecab6ddf098c
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═65d3653d-5734-40ec-a057-aaa1e509968a
# ╠═8d90a98b-09df-4834-ab9d-0aed44d6206b
# ╠═dd8a0445-8c94-4f43-b22a-1146a7f55d63
# ╠═43c973e6-ad30-4bf8-93b3-4924bc498927
# ╠═f875e0e9-2f05-4420-8f40-76284be58e03
# ╠═89cc7dd4-6643-463d-9fa3-bfc2859476f2
# ╠═4851087d-2679-4409-b44f-cb22f778974e
# ╠═4e76a22b-631a-4478-b3bd-47be18496b3a
# ╠═3317cfec-be70-418c-944a-0a6741f03cf3
