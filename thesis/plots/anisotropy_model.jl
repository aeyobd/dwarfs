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

	using PyFITS
end

# ╔═╡ 1482481d-a5f3-48c2-a4d2-1353afe7fd72
using OrderedCollections

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ c1b20bd0-bc18-4c53-b687-5f6b792fdcd0
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ 8081e889-c1d4-4e50-b8e1-91aa8e82adf1
function modeldir(modelname)
	joinpath(ENV["DWARFS_ROOT"], "analysis",modelname,)
end

# ╔═╡ 6161cd69-594b-4eb9-83d9-12a682db7c88


# ╔═╡ f49c9325-4319-4435-ac33-a8980e1c4f7f
function get_profiles(modelname)
	out = LilGuys.Output(modeldir(modelname))
	idx_f = TOML.parsefile(joinpath(modeldir(modelname), "orbital_properties.toml"))["idx_f"]

	r_bins_i = LilGuys.bins_equal_number(radii(out[1]), nothing, num_per_bin=3000)
	prof_i = LilGuys.β_prof(out[1], r_bins=r_bins_i)

	r_bins = LilGuys.bins_equal_number(radii(out[idx_f]), nothing, num_per_bin=3000)

	prof_f = LilGuys.β_prof(out[idx_f], r_bins=r_bins)

	return r_bins_i, prof_i, r_bins, prof_f
end

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
function plot_anisotropies(modelname; title="")
	r_bins_i, prof_i, r_bins, prof_f = get_profiles(modelname)
	
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = L"$r$ / kpc",
		ylabel = L"$\beta$",
		xscale=log10,
		xticks = [0.1, 1, 10],
		title = title
	)


	lines!(midpoints(r_bins_i), prof_i[2], label="initial", linestyle=:dot)
	lines!(midpoints(r_bins), prof_f[2], label="final", linestyle=:solid)

	axislegend(position=:lb)
	ylims!(-1, 1)
	xlims!(1, 100)

	fig
end

# ╔═╡ 415f249c-de75-454f-8ab8-be3bec7e4222
plot_anisotropies("sculptor/1e6_new_v31_r3.2/orbit_smallperi")

# ╔═╡ 767a0c19-c50d-4764-b14a-164eafca8747
@savefig "anisotropy_i_f" plot_anisotropies("sculptor/1e6_v43_r5_beta0.2_a4/orbit_smallperi")

# ╔═╡ 10eb2b80-285e-4572-8cdd-3b09b5316681
h = LilGuys.Snapshot(joinpath(ENV["DWARFS_ROOT"], "analysis/isolation/1e6_beta0.2_a4/fiducial/combined.hdf5/1"))

# ╔═╡ dcda8c6a-66cd-41a5-87db-2ed2a4cc31dc
let
	bins = LilGuys.quantile(radii(h), LinRange(0, 1, 100))
	_, β = LilGuys.β_prof(h, r_bins=bins)

	fig = Figure()
	ax = Axis(fig[1,1], xscale=log10, xticks=[1, 10, 100])
	lines!(midpoints(bins), β)

	fig
end

# ╔═╡ 1bb4bc70-9013-42bd-a985-f032e5bddd42
LilGuys.NFW(r_circ_max=5, v_circ_max=43/V2KMS).r_s*4

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═c1b20bd0-bc18-4c53-b687-5f6b792fdcd0
# ╠═1482481d-a5f3-48c2-a4d2-1353afe7fd72
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═8081e889-c1d4-4e50-b8e1-91aa8e82adf1
# ╠═6161cd69-594b-4eb9-83d9-12a682db7c88
# ╠═f49c9325-4319-4435-ac33-a8980e1c4f7f
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═415f249c-de75-454f-8ab8-be3bec7e4222
# ╠═767a0c19-c50d-4764-b14a-164eafca8747
# ╠═10eb2b80-285e-4572-8cdd-3b09b5316681
# ╠═dcda8c6a-66cd-41a5-87db-2ed2a4cc31dc
# ╠═1bb4bc70-9013-42bd-a985-f032e5bddd42
