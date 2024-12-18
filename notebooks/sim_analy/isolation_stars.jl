### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 6238c3fa-1974-11ef-1539-2f50fb31fe9a
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya

	import LilGuys as lguys
end

# ╔═╡ 141d3b97-e344-4760-aa12-a48b1afb125c
begin 
	import DataFrames, CSV
	using HDF5
	import NaNMath as nm
end

# ╔═╡ 2845abc7-53cd-4996-a04b-e062ca3e63b7
md"""
Explores the properties and evolution of the stars in an isolation run.
"""

# ╔═╡ feb6cc17-a25d-4ba8-a152-78412f422b80
import TOML

# ╔═╡ bb3b3395-a64d-4907-a4bd-1b3f621010b1
import DensityEstimators: histogram

# ╔═╡ 5b1dd353-a437-47cd-94be-7da9684581da
modeldir = "/astro/dboyea/sculptor/isolation/1e6/fiducial/"

# ╔═╡ 28ba4f0f-6cc6-44e2-a7bc-4eee460d91b0
starsname = "stars/exp2d_rs0.10"

# ╔═╡ 21adbbe7-c8cc-4094-9e75-b68d97fa211a
starsfile = "$(starsname)_stars.hdf5"

# ╔═╡ b38a74d3-2a82-412d-981e-8466e403746e
fig_name = joinpath(modeldir, dirname(starsfile), "figures", basename(starsname))

# ╔═╡ f13d237d-ce33-43aa-a1c7-796ac502c9fe
halo_factor = 1

# ╔═╡ 7e5edb0d-afd7-4437-86c0-f2532e25971e
md"""
# file loading
"""

# ╔═╡ 83fbe8a1-5d95-458a-87f2-b1d5c66846d2
begin 
	stellar_profiles, profile_labels = lguys.read_structs_from_hdf5(joinpath(modeldir, starsname * "_stars_profiles.hdf5"), lguys.StellarProfile3D)
	
	snap_idxs = parse.(Int, profile_labels)
	stellar_profiles = stellar_profiles[sortperm(snap_idxs)]
	snap_idxs = sort(snap_idxs)

	snap_idxs, stellar_profiles
end
	

# ╔═╡ b4e107d1-5b14-4e64-80ce-a02ece3b41b7
prob_df = lguys.read_hdf5_table(joinpath(modeldir, starsfile))

# ╔═╡ 312576e2-16da-4285-9c19-a7a8005acf25
paramname = "$(starsname)"

# ╔═╡ c76c8acc-ea88-4ce1-81cb-be0b67ef23fd
profile = lguys.load_profile(joinpath(modeldir, paramname) * ".toml")

# ╔═╡ 7dfc4066-eee5-4118-9d12-48861aa66e03
NamedTuple(d::Dict) = (; zip(Symbol.(keys(d)), values(d))...)

# ╔═╡ 0671aa68-3d02-4b6c-a238-adffede13bd8
ρ_s(r) = lguys.calc_ρ(profile, r)

# ╔═╡ 01742641-9e37-45bc-a0c4-a4a662fde1f9
star_prof_i = stellar_profiles[2]

# ╔═╡ 0297ef30-d37b-48a1-9051-dc7e2a5df801
star_prof_f = stellar_profiles[end]

# ╔═╡ a888cc8f-e27c-4989-a535-6a2862c60c91
probabilities = prob_df.probability

# ╔═╡ e32e2d66-2dad-4f09-84f1-a3081e3891a7
out = lguys.Output("$modeldir/out/", weights=probabilities)

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
	fig = Figure()
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc", title="initial")

	bins = LinRange(-1, 1, 100)
	probs = snap_i.weights
	Arya.hist2d!(ax, snap_i.positions[1, :], snap_i.positions[2, :], weights=probs, bins = bins)

	ax2 = Axis(fig[1,2], aspect=1,
	xlabel = "x / kpc", ylabel="y/kpc",
	title="final")

	probs = snap_f.weights

	Arya.hist2d!(ax2, snap_f.positions[1, :], snap_f.positions[2, :], weights=probs, bins = bins)
	hideydecorations!(ax2)
	linkaxes!(ax, ax2)

	save(fig_name * "_xy_ini_fin.pdf", fig)
	fig
end

# ╔═╡ de8ebcd0-d6e1-4b16-aaec-5bcd47cad1bd
function plot_ρ_s!(snap; bins=200, kwargs...)
	rs = lguys.calc_r(snap)
	ps = probabilities[snap.index]
	h = histogram(log10.(rs), bins, weights=ps)
	ρ = lguys.calc_ρ_from_hist(10 .^ h.bins, h.values)
	x = (lguys.midpoints(h.bins))
	lines!(x, log10.(ρ); kwargs...)
end

# ╔═╡ 98d2168c-f450-41e2-9b9d-2880a662f841
sort(lguys.calc_r(snap_i))

# ╔═╡ 8b69303d-c992-447a-aaf6-af5839173b1a
snap_i.x_cen

# ╔═╡ 7ad553e8-50ac-41c3-b461-3c9ba2cdef17
snap_i.positions

# ╔═╡ b76c91ff-0928-40ad-9263-c455f804b6f5
out.times[end] * lguys.T2GYR

# ╔═╡ 6065db1e-1387-4735-ad9d-0f401d9f0f73
let 
	fig = Figure()

	ax = Axis(fig[1,1], xlabel=L"\log\, r / \textrm{kpc}", ylabel =  L"\log\, \rho_\star\; [10^{10} M_\odot / \textrm{kpc}^3]", 
		limits=((-1.5, 0.8), (-15, 3))
		)

	lines!(star_prof_i.log_r, log10.(star_prof_i.rho))
	lines!(star_prof_f.log_r, log10.(star_prof_f.rho))

	log_r_pred = LinRange(-2, 2, 1000)
	ρ_s_pred = ρ_s.(10 .^ log_r_pred)

	lines!(log_r_pred, log10.(ρ_s_pred), label="expected", color="black", linestyle=:dot)

	save(fig_name * "_rho_ini_fin.pdf", fig)

	axislegend(ax)
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
	
	hist!(probabilities, bins=20, )
	fig
end

# ╔═╡ 93abb048-2a1a-468c-86f5-abba3a0e92e5
v_los = snap_f.velocities[1, :] .* lguys.V2KMS

# ╔═╡ 0f0275b6-473c-4b3e-8e3f-2be7903eec32
r = @. sqrt((snap_f.positions[2, :] .- snap_f.x_cen[2, end])^2 + (snap_f.positions[3, :] .- snap_f.x_cen[3, end])^2)

# ╔═╡ 1e7ac1c1-5d83-4d43-8949-0953a4ef098f
m_los = probabilities[snap_f.index]

# ╔═╡ 5b2f984b-a9cd-4c7f-a901-e2e6df26f5e4
begin
	r_filt = m_los .> 1e-15 * maximum(m_los)
	r_filt .&= r .< 1
end

# ╔═╡ dd21570b-9ba5-44b3-a8c1-76251b492eef
import StatsBase: std, weights, mean

# ╔═╡ 95386397-7303-48a6-9720-c70384b8ec7a
let
	vs = v_los[r_filt]
	w = m_los[r_filt]

	global σ_v = std(vs, weights(w))
	global μ_v = mean(vs, weights(w))

	fig, ax = FigAxis(
		xlabel = L"$v_\textrm{los}$ / km\,s$^{-1}$"
	)
	
	hist!(vs, weights=w, bins=50, normalization=:pdf)

	x_mod = LinRange(-30, 30, 1000)
	y_mod = lguys.gaussian.(x_mod, μ_v, σ_v)
	lines!(x_mod, y_mod, label="N($μ_v, $σ_v)")

	axislegend()
	fig
end

# ╔═╡ 2f6da0ff-eef1-407d-a5e4-c3035d049688
let
	skip = 1

	ts = times[snap_idxs]

	sigmas = [p.sigma_vx for p in stellar_profiles]

	fig, ax = FigAxis(
		xlabel = "time / Gyr",
		ylabel = L"\sigma_v / \textrm{km s^{-1}}",
	)
	
	scatter!(ts, sigmas * lguys.V2KMS)

	save(fig_name * "_sigma_v_t.pdf", fig)


	fig
end

# ╔═╡ Cell order:
# ╠═2845abc7-53cd-4996-a04b-e062ca3e63b7
# ╠═6238c3fa-1974-11ef-1539-2f50fb31fe9a
# ╠═feb6cc17-a25d-4ba8-a152-78412f422b80
# ╠═bb3b3395-a64d-4907-a4bd-1b3f621010b1
# ╠═141d3b97-e344-4760-aa12-a48b1afb125c
# ╠═5b1dd353-a437-47cd-94be-7da9684581da
# ╠═28ba4f0f-6cc6-44e2-a7bc-4eee460d91b0
# ╠═21adbbe7-c8cc-4094-9e75-b68d97fa211a
# ╠═b38a74d3-2a82-412d-981e-8466e403746e
# ╠═f13d237d-ce33-43aa-a1c7-796ac502c9fe
# ╟─7e5edb0d-afd7-4437-86c0-f2532e25971e
# ╠═83fbe8a1-5d95-458a-87f2-b1d5c66846d2
# ╠═b4e107d1-5b14-4e64-80ce-a02ece3b41b7
# ╠═e32e2d66-2dad-4f09-84f1-a3081e3891a7
# ╠═312576e2-16da-4285-9c19-a7a8005acf25
# ╠═c76c8acc-ea88-4ce1-81cb-be0b67ef23fd
# ╠═7dfc4066-eee5-4118-9d12-48861aa66e03
# ╠═0671aa68-3d02-4b6c-a238-adffede13bd8
# ╠═01742641-9e37-45bc-a0c4-a4a662fde1f9
# ╠═0297ef30-d37b-48a1-9051-dc7e2a5df801
# ╠═3150cdfd-7573-4db9-86b7-ef614150a7b9
# ╠═2cc047db-ae05-42c5-898f-702ae3b83bd6
# ╠═ddc895b1-698e-4f89-bf4f-d5ede7ef2c04
# ╠═052d229a-0362-42eb-b03c-fdab0f0bc6b4
# ╠═a888cc8f-e27c-4989-a535-6a2862c60c91
# ╠═4a32cad7-c18a-4494-aca9-e292d64ef6b7
# ╠═addf19a4-088b-4aff-89a9-df73e8049f2c
# ╠═de8ebcd0-d6e1-4b16-aaec-5bcd47cad1bd
# ╠═98d2168c-f450-41e2-9b9d-2880a662f841
# ╠═8b69303d-c992-447a-aaf6-af5839173b1a
# ╠═7ad553e8-50ac-41c3-b461-3c9ba2cdef17
# ╠═b76c91ff-0928-40ad-9263-c455f804b6f5
# ╠═6065db1e-1387-4735-ad9d-0f401d9f0f73
# ╠═cea9464d-d620-4fc9-9224-d24e075bf42a
# ╠═d8f546d3-9e2e-4703-b652-5bea7bbbbd26
# ╠═193273c9-5b13-4af6-a345-4326cdebcf04
# ╠═93abb048-2a1a-468c-86f5-abba3a0e92e5
# ╠═0f0275b6-473c-4b3e-8e3f-2be7903eec32
# ╠═1e7ac1c1-5d83-4d43-8949-0953a4ef098f
# ╠═5b2f984b-a9cd-4c7f-a901-e2e6df26f5e4
# ╠═dd21570b-9ba5-44b3-a8c1-76251b492eef
# ╠═95386397-7303-48a6-9720-c70384b8ec7a
# ╠═2f6da0ff-eef1-407d-a5e4-c3035d049688
