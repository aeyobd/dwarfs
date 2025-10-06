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

# ╔═╡ 3d50a77e-7365-41c1-a306-131d9672ca6d
using PyFITS

# ╔═╡ c79a6966-e559-487c-a378-2e258c1fdbfa
using OrderedCollections

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ e0815f0a-0955-4899-9688-c6c70071f510
import HDF5

# ╔═╡ 7845f35d-1112-4c85-bdab-7a6a9f11cb1d
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ 1c5f3f6b-57c2-4f56-90e0-3fc874b86e9f
md"""
# Utilities
"""

# ╔═╡ 45ea911c-b082-483e-b7f5-902088544ede
function starsdir(galaxyname, modelname, starsname)
	return joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, "stars", starsname)
end

# ╔═╡ e08d6905-318e-490c-ad80-a5d9ed8a6a6a
function get_obs_props(galaxyname)
	return TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))
end

# ╔═╡ c2ec04cb-d7c8-4abd-8bae-94b2e5c6af0c
function get_t_f(galaxyname, modelname, _)
	name = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname)
	
    df = TOML.parsefile("$name/orbital_properties.toml")
    idx_f = df["idx_f"]
    times = HDF5.h5open("$name/centres.hdf5") do f
        return f["times"][:]
    end
    
    return times[idx_f]
end

# ╔═╡ a04583f2-fd84-4640-8c48-719c59f717cb
function get_scalars(args...)
	file = joinpath(starsdir(args...), "stellar_profiles_3d_scalars.fits")

	idx_f = TOML.parsefile(joinpath(starsdir(args...), "../../orbital_properties.toml"))["idx_f"]
	df = read_fits(file)

	df.time .-= get_t_f(args...)
	df
end

# ╔═╡ 70b55c06-3e81-4e27-a313-69d7e4028ece
function compare_sigma_v_time(orbits, obs_props)
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = "time / Gyr",
		ylabel = L"velocity disperison / km\,s$^{-1}$"
	)

	hlines!(obs_props["sigma_v"], color=:black)

	hspan!(obs_props["sigma_v"] - obs_props["sigma_v_err"], obs_props["sigma_v"] +obs_props["sigma_v_err"], alpha=0.1, color=:black)

	dfs = OrderedDict(k => get_scalars(v...) for (k, v) in orbits)

	for (label, df) in dfs
		lines!(df.time * T2GYR, df.sigma_v * V2KMS, label=label)
	end


	axislegend(position=:lb)


	ax2 = Axis(fig[2,1], 
		xlabel = "time / Gyr",
		ylabel = L"$r_h$ / kpc"
	)
	
	for (label, df) in dfs
		lines!(df.time * T2GYR, df.r_h, label=label)
	end

	hidexdecorations!(ax, ticks=false, minorticks=false)

	linkxaxes!(ax, ax2)

	R_h = obs_props["R_h"]
	α = LilGuys.r_h(LilGuys.Exp2D()) / LilGuys.R_h(LilGuys.Exp2D())
	r_h_l = LilGuys.arcmin2kpc(R_h, obs_props["distance"] - obs_props["distance_err"]) * α
	r_h_h = LilGuys.arcmin2kpc(R_h, obs_props["distance"] + obs_props["distance_err"]) * α
	
	hspan!(r_h_l, r_h_h, alpha=0.1, color=:black)


	fig
end

# ╔═╡ 82bc10ae-cbf9-4f1f-8167-04eabc136b0b
md"""
# Comparison
"""

# ╔═╡ 55027866-d309-4557-b8a9-9f73a20b7785
scl_orbits = OrderedDict(
	#"mean" => get_scalars("sculptor", "1e7_V31_r3.2/orbit_mean", "exp2d_rs0.10"),
	"smallperi" => get_scalars("sculptor", "1e7_new_v31_r3.2/orbit_smallperi", "exp2d_rs0.10"),
	#"plummer" => get_scalars("sculptor", "1e7_new_v31_r3.2/orbit_smallperi", "plummer_rs0.20"),
	#"lmc" => get_scalars("sculptor", "1e7_V31_r4.2/vasiliev24_L3M11_2x_smallperilmc", "exp2d_rs0.10"),
)
	
	

# ╔═╡ 80d2464d-8a1a-4da5-bb62-4b28e843e7fe
scl_obs = get_obs_props("sculptor")

# ╔═╡ 17eef1aa-9fce-43d8-bd24-bd8613f93ded
@savefig "scl_sigma_v_time" compare_sigma_v_time(OrderedDict(
	"fiducial" => ("sculptor", "1e6_new_v43_r7/orbit_smallperi","exp2d_rs0.10"),
	"cored" => ("sculptor", "1e6_Ms0.54_rs1.08_c1/orbit_smallperi","exp2d_rs0.10"),
	"anisotropic" =>  ("sculptor", "1e6_v43_r5_beta0.2_a4/orbit_smallperi", "exp2d_rs0.10"),
	"oblate" => ("sculptor", "1e6_v48_r7_oblate_0.5/orbit_smallperi", "exp2d_rs0.10"),

), scl_obs)

# ╔═╡ cd7a0e9c-d1ea-43d5-9d52-6b264366d731


# ╔═╡ 0a8968f5-84ff-4bd8-bf92-81e68cb26fb9
@savefig "scl_sigma_v_time_weight" compare_sigma_v_time(OrderedDict(
	"fiducial" => ("sculptor", "1e7_new_v31_r3.2/orbit_smallperi","exp2d_rs0.10"),
	"heavy" => ("sculptor", "1e6_new_v43_r7/orbit_smallperi","exp2d_rs0.10"),
	"light" => ("sculptor", "1e6_new_v25_r2.5/orbit_smallperi","exp2d_rs0.13"),

), scl_obs)

# ╔═╡ bcfd1fd8-f254-4c12-a0aa-434585b473cd
@savefig "scl_sigma_v_time_orbit" compare_sigma_v_time(OrderedDict(
	"fiducial" => ("sculptor", "1e7_new_v31_r3.2/orbit_smallperi","exp2d_rs0.10"),
	# "mean" => ("sculptor", "1e6_new_v31_r3.2/orbit_mean","exp2d_rs0.10"),
	"impact" => ("sculptor", "1e6_new_v31_r3.2/L3M11_9Gyr_smallperi.a4","exp2d_rs0.10"),

), scl_obs)

# ╔═╡ 9cb7699d-e390-4dd5-b051-a5a48d2316bc
umi_orbits = OrderedDict(
	"smallperi" => ("ursa_minor", "1e7_new_v38_r4.0/orbit_smallperi.5", "exp2d_rs0.10"),
)

# ╔═╡ 28c9b087-0101-4710-9940-81138fbf7cef
umi_obs = get_obs_props("ursa_minor")

# ╔═╡ 9d7d5739-4172-4eb8-9a04-dc2a727d9dd4
@savefig "umi_sigma_v_time" compare_sigma_v_time(umi_orbits, umi_obs)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═e0815f0a-0955-4899-9688-c6c70071f510
# ╠═3d50a77e-7365-41c1-a306-131d9672ca6d
# ╠═7845f35d-1112-4c85-bdab-7a6a9f11cb1d
# ╠═c79a6966-e559-487c-a378-2e258c1fdbfa
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═1c5f3f6b-57c2-4f56-90e0-3fc874b86e9f
# ╠═45ea911c-b082-483e-b7f5-902088544ede
# ╠═e08d6905-318e-490c-ad80-a5d9ed8a6a6a
# ╠═a04583f2-fd84-4640-8c48-719c59f717cb
# ╠═70b55c06-3e81-4e27-a313-69d7e4028ece
# ╠═c2ec04cb-d7c8-4abd-8bae-94b2e5c6af0c
# ╟─82bc10ae-cbf9-4f1f-8167-04eabc136b0b
# ╠═55027866-d309-4557-b8a9-9f73a20b7785
# ╠═80d2464d-8a1a-4da5-bb62-4b28e843e7fe
# ╠═17eef1aa-9fce-43d8-bd24-bd8613f93ded
# ╠═cd7a0e9c-d1ea-43d5-9d52-6b264366d731
# ╠═0a8968f5-84ff-4bd8-bf92-81e68cb26fb9
# ╠═bcfd1fd8-f254-4c12-a0aa-434585b473cd
# ╠═9cb7699d-e390-4dd5-b051-a5a48d2316bc
# ╠═28c9b087-0101-4710-9940-81138fbf7cef
# ╠═9d7d5739-4172-4eb8-9a04-dc2a727d9dd4
