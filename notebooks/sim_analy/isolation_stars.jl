### A Pluto.jl notebook ###
# v0.20.8

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

# ╔═╡ 6238c3fa-1974-11ef-1539-2f50fb31fe9a
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya
end

# ╔═╡ b8feefa7-8ec7-4c57-afb7-296aaf6aa24b
using PyFITS

# ╔═╡ 451c25f8-2a47-4503-a5c6-658b008d13c2
using PlutoUI

# ╔═╡ 141d3b97-e344-4760-aa12-a48b1afb125c
begin 
	import DataFrames, CSV
	using HDF5
	import NaNMath as nm
end

# ╔═╡ d9849214-8db0-48d0-bdd1-894e6b0c97f9
using Printf

# ╔═╡ 2845abc7-53cd-4996-a04b-e062ca3e63b7
md"""
Explores the properties and evolution of the stars in an isolation run.
This notebook double checks that there is minimal change to the stellar profile in isolation.
"""

# ╔═╡ 4af1e554-05ed-485f-a15a-c588daf5fd3a
md"""
# Setup
"""

# ╔═╡ 146a5666-cf82-45c7-ba56-bc8014cd41c7
CairoMakie.activate!(type=:png)

# ╔═╡ feb6cc17-a25d-4ba8-a152-78412f422b80
import TOML

# ╔═╡ bb3b3395-a64d-4907-a4bd-1b3f621010b1
import DensityEstimators: histogram

# ╔═╡ dd21570b-9ba5-44b3-a8c1-76251b492eef
import StatsBase: std, weights, mean

# ╔═╡ a3ad8fe3-957e-4e22-809b-9c208400e1fe
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

# ╔═╡ 4340af50-de93-4a8d-b850-6d47d7685995
@bind inputs confirm(notebook_inputs(;
	modelname = TextField(default="ursa_minor/1e6_v37_r5.0"),
	isoname = TextField(default="1e/fiducial"),
	starsname = TextField(default="exp2d_rs0.13"),
))

# ╔═╡ ab861f21-5f5c-4f67-9edb-0a0e79d8ceec
modelname = inputs.modelname

# ╔═╡ 5337775d-8150-4fe9-bdaa-cbac0dcce84b
isoname = inputs.isoname

# ╔═╡ 28ba4f0f-6cc6-44e2-a7bc-4eee460d91b0
starsname = inputs.starsname

# ╔═╡ 5b1dd353-a437-47cd-94be-7da9684581da
modeldir = ENV["DWARFS_ROOT"] * "/analysis/$modelname"

# ╔═╡ 3aac17ed-d5a9-491c-b467-3a32ad8a3c39
isodir = ENV["DWARFS_ROOT"] * "/analysis/isolation/$isoname"

# ╔═╡ 21adbbe7-c8cc-4094-9e75-b68d97fa211a
starsfile = "probabilities_stars.hdf5"

# ╔═╡ 00d66dcc-ba1f-4cd1-bbc0-380fed17aeb4
starsdir = joinpath(modeldir, "stars", starsname)

# ╔═╡ b38a74d3-2a82-412d-981e-8466e403746e
FIGDIR = joinpath(modeldir, "stars", starsname, "figures")

# ╔═╡ b8942c7f-5e20-42ad-bea6-b3c56d0840ba
begin
	
	import LilGuys as lguys
	using LilGuys
	FIGDIR
end

# ╔═╡ bffc3bb2-9be5-40b3-95f9-05b2ed27e680
if isdir(starsdir) && !isdir(FIGDIR)
	mkdir(FIGDIR)
end

# ╔═╡ 7e5edb0d-afd7-4437-86c0-f2532e25971e
md"""
# file loading
"""

# ╔═╡ 83fbe8a1-5d95-458a-87f2-b1d5c66846d2
begin 
	stellar_profiles = lguys.read_ordered_structs(joinpath(modeldir, "stars", starsname * "/stellar_profiles_3d.hdf5"), lguys.DensityProfile)
	
	
	snap_idxs = first.(stellar_profiles)
	stellar_profiles = last.(stellar_profiles)

	snap_idxs, stellar_profiles
end
	

# ╔═╡ 7ebfe65b-e717-4a06-ab73-3f0382af1251
scalars = read_fits(joinpath(modeldir, "stars", starsname * "/stellar_profiles_3d_scalars.fits"))

# ╔═╡ b4e107d1-5b14-4e64-80ce-a02ece3b41b7
prob_df = lguys.read_hdf5_table(joinpath(modeldir, "stars", starsname, starsfile))

# ╔═╡ c76c8acc-ea88-4ce1-81cb-be0b67ef23fd
profile = lguys.load_profile(joinpath(modeldir, "stars", starsname, "profile.toml"))

# ╔═╡ 01742641-9e37-45bc-a0c4-a4a662fde1f9
star_prof_i = stellar_profiles[1]

# ╔═╡ 0297ef30-d37b-48a1-9051-dc7e2a5df801
star_prof_f = stellar_profiles[end]

# ╔═╡ a888cc8f-e27c-4989-a535-6a2862c60c91
probabilities = prob_df.probability

# ╔═╡ e32e2d66-2dad-4f09-84f1-a3081e3891a7
out = lguys.Output(isodir, weights=probabilities)

# ╔═╡ 3150cdfd-7573-4db9-86b7-ef614150a7b9
times = out.times * lguys.T2GYR

# ╔═╡ 2cc047db-ae05-42c5-898f-702ae3b83bd6
idx_i = 20 + 1; idx_f = length(out)

# ╔═╡ ddc895b1-698e-4f89-bf4f-d5ede7ef2c04
out[1].x_cen

# ╔═╡ 052d229a-0362-42eb-b03c-fdab0f0bc6b4
begin
	snap_i = out[idx_i]
	snap_f = out[idx_f]
end

# ╔═╡ 4a32cad7-c18a-4494-aca9-e292d64ef6b7
md"""
# Plots
"""

# ╔═╡ addf19a4-088b-4aff-89a9-df73e8049f2c
let
	bins = LinRange(-1, 1, 100)

	fig = Figure()
	
	ax = Axis(fig[1,1], 
		aspect=DataAspect(),
		xlabel = "x / kpc", 
		ylabel="y/kpc", 
		title="initial"
	)

	Arya.hist2d!(ax, snap_i.positions[1, :], snap_i.positions[2, :],
		weights=snap_i.weights, bins = bins)

	ax2 = Axis(fig[1,2], 
		aspect=DataAspect(),
		xlabel = "x / kpc", 
		ylabel="y/kpc",
		title="final"
	)

	Arya.hist2d!(ax2, snap_f.positions[1, :], snap_f.positions[2, :], 
		weights=snap_f.weights, bins = bins)
	
	hideydecorations!(ax2)
	linkaxes!(ax, ax2)

	@savefig "xy_ini_fin"
	fig
end

# ╔═╡ 6065db1e-1387-4735-ad9d-0f401d9f0f73
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log\, r / \textrm{kpc}", ylabel =  L"\log\, \rho_\star\; [10^{10} M_\odot / \textrm{kpc}^3]", 
		limits=((-2.5, 0.8), (-15, 3))
		)

	lines!(star_prof_i.log_r, log10.(star_prof_i.rho), label="initial")
	lines!(star_prof_f.log_r, log10.(star_prof_f.rho), label="final")

	log_r_pred = LinRange(-2, 2, 1000)
	ρ_s_pred = lguys.density.(profile, 10 .^ log_r_pred)

	lines!(log_r_pred, log10.(ρ_s_pred), label="expected", color="black", linestyle=:dot)

	@savefig "rho_ini_fin"

	axislegend(ax, position=:lb)
	fig
end

# ╔═╡ cea9464d-d620-4fc9-9224-d24e075bf42a
star_prof_i

# ╔═╡ d8f546d3-9e2e-4703-b652-5bea7bbbbd26
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel="log r / kpc containing stellar mass")

	q = star_prof_i.quantiles
	for i in eachindex(q)
		Mq = [log10.(p.r_quantile[i]) for p in stellar_profiles]
		t = times[snap_idxs]
	
		lines!(t, Mq, color=q[i], colorrange=extrema(q), label="$(q[i])")
	end

	Legend(fig[1, 2], ax, "quantile")
	fig
end

# ╔═╡ 193273c9-5b13-4af6-a345-4326cdebcf04
let
	fig = Figure()
	ax = Axis(fig[1,1],
	yscale=log10, xlabel="probability", ylabel="count")
	, yticks=Makie.automatic)
	hist!(probabilities, bins=20, )
	fig
end

# ╔═╡ 7072fd54-28c9-4e0a-b2b1-d2bef00c81bc
md"""
## Line of sight velocity dispersion
Here, we choose to project the snapshot orthographically onto the y-z plane and use the x velocities as the LOS velocities.
"""

# ╔═╡ ca2a6553-71b4-4b2d-a470-65193a42833b
r_max_v_los = 1

# ╔═╡ 93abb048-2a1a-468c-86f5-abba3a0e92e5
v_los = snap_f.velocities[1, :] .* lguys.V2KMS

# ╔═╡ 0f0275b6-473c-4b3e-8e3f-2be7903eec32
r = @. sqrt((snap_f.positions[2, :] .- snap_f.x_cen[2, end])^2 + (snap_f.positions[3, :] .- snap_f.x_cen[3, end])^2)

# ╔═╡ 1e7ac1c1-5d83-4d43-8949-0953a4ef098f
m_los = probabilities[snap_f.index]

# ╔═╡ 5b2f984b-a9cd-4c7f-a901-e2e6df26f5e4
r_filt = (
	(m_los .> 0 * maximum(m_los))
	.& ( r .< r_max_v_los)
)

# ╔═╡ 1e86f26a-df52-4a81-98f6-e12ccea2daee
σ_v = lguys.std(v_los[r_filt], m_los[r_filt])

# ╔═╡ 2ed8758b-688b-42fc-9dbe-9b259adc08dc
μ_v = lguys.mean(v_los[r_filt], m_los[r_filt])

# ╔═╡ 95386397-7303-48a6-9720-c70384b8ec7a
let

	fig, ax = FigAxis(
		xlabel = L"$v_\textrm{los}$ / km\,s$^{-1}$",
		ylabel = "density",
	)
	
	hist!(v_los[r_filt], weights=m_los[r_filt], bins=50, normalization=:pdf)

	x_mod = LinRange(-30, 30, 1000)
	y_mod = lguys.gaussian.(x_mod, μ_v, σ_v)
	lines!(x_mod, y_mod, label=@sprintf("v/kms ~ N(%0.2f,%0.2f)", μ_v, σ_v), color=COLORS[2])

	axislegend()
	fig
end

# ╔═╡ 8bb4aa49-172a-4214-9c13-e644c295e90c
md"""
The below plot shows the velocity dispersion as a function of time.
This should be constant ideally.
"""

# ╔═╡ 2f6da0ff-eef1-407d-a5e4-c3035d049688
let
	fig, ax = FigAxis(
		xlabel = "time / Gyr",
		ylabel = L"\sigma_v / \textrm{km s^{-1}}",
	)
	
	ts = scalars.time
	sigmas = scalars.sigma_v
	scatter!(ts, sigmas * lguys.V2KMS)

	@savefig "sigma_v_t"

	fig
end

# ╔═╡ Cell order:
# ╟─2845abc7-53cd-4996-a04b-e062ca3e63b7
# ╠═4340af50-de93-4a8d-b850-6d47d7685995
# ╟─4af1e554-05ed-485f-a15a-c588daf5fd3a
# ╠═6238c3fa-1974-11ef-1539-2f50fb31fe9a
# ╠═b8feefa7-8ec7-4c57-afb7-296aaf6aa24b
# ╠═b8942c7f-5e20-42ad-bea6-b3c56d0840ba
# ╠═146a5666-cf82-45c7-ba56-bc8014cd41c7
# ╠═451c25f8-2a47-4503-a5c6-658b008d13c2
# ╠═feb6cc17-a25d-4ba8-a152-78412f422b80
# ╠═bb3b3395-a64d-4907-a4bd-1b3f621010b1
# ╠═141d3b97-e344-4760-aa12-a48b1afb125c
# ╠═dd21570b-9ba5-44b3-a8c1-76251b492eef
# ╠═d9849214-8db0-48d0-bdd1-894e6b0c97f9
# ╠═a3ad8fe3-957e-4e22-809b-9c208400e1fe
# ╠═ab861f21-5f5c-4f67-9edb-0a0e79d8ceec
# ╠═5337775d-8150-4fe9-bdaa-cbac0dcce84b
# ╠═28ba4f0f-6cc6-44e2-a7bc-4eee460d91b0
# ╠═5b1dd353-a437-47cd-94be-7da9684581da
# ╠═3aac17ed-d5a9-491c-b467-3a32ad8a3c39
# ╠═21adbbe7-c8cc-4094-9e75-b68d97fa211a
# ╠═00d66dcc-ba1f-4cd1-bbc0-380fed17aeb4
# ╠═b38a74d3-2a82-412d-981e-8466e403746e
# ╠═bffc3bb2-9be5-40b3-95f9-05b2ed27e680
# ╟─7e5edb0d-afd7-4437-86c0-f2532e25971e
# ╠═83fbe8a1-5d95-458a-87f2-b1d5c66846d2
# ╠═7ebfe65b-e717-4a06-ab73-3f0382af1251
# ╠═b4e107d1-5b14-4e64-80ce-a02ece3b41b7
# ╠═e32e2d66-2dad-4f09-84f1-a3081e3891a7
# ╠═c76c8acc-ea88-4ce1-81cb-be0b67ef23fd
# ╠═01742641-9e37-45bc-a0c4-a4a662fde1f9
# ╠═0297ef30-d37b-48a1-9051-dc7e2a5df801
# ╠═3150cdfd-7573-4db9-86b7-ef614150a7b9
# ╠═2cc047db-ae05-42c5-898f-702ae3b83bd6
# ╠═ddc895b1-698e-4f89-bf4f-d5ede7ef2c04
# ╠═052d229a-0362-42eb-b03c-fdab0f0bc6b4
# ╠═a888cc8f-e27c-4989-a535-6a2862c60c91
# ╟─4a32cad7-c18a-4494-aca9-e292d64ef6b7
# ╟─addf19a4-088b-4aff-89a9-df73e8049f2c
# ╠═6065db1e-1387-4735-ad9d-0f401d9f0f73
# ╠═cea9464d-d620-4fc9-9224-d24e075bf42a
# ╠═d8f546d3-9e2e-4703-b652-5bea7bbbbd26
# ╠═193273c9-5b13-4af6-a345-4326cdebcf04
# ╟─7072fd54-28c9-4e0a-b2b1-d2bef00c81bc
# ╠═ca2a6553-71b4-4b2d-a470-65193a42833b
# ╠═93abb048-2a1a-468c-86f5-abba3a0e92e5
# ╠═0f0275b6-473c-4b3e-8e3f-2be7903eec32
# ╠═1e7ac1c1-5d83-4d43-8949-0953a4ef098f
# ╠═5b2f984b-a9cd-4c7f-a901-e2e6df26f5e4
# ╠═1e86f26a-df52-4a81-98f6-e12ccea2daee
# ╠═2ed8758b-688b-42fc-9dbe-9b259adc08dc
# ╟─95386397-7303-48a6-9720-c70384b8ec7a
# ╟─8bb4aa49-172a-4214-9c13-e644c295e90c
# ╟─2f6da0ff-eef1-407d-a5e4-c3035d049688
