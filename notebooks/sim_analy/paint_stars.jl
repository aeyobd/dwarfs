### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ a893932c-f184-42bc-9a0e-0960f10520aa
begin
	import Pkg
	Pkg.activate()
	
	import LilGuys as lguys
end

# ╔═╡ 641946b3-e6f2-4d6d-8777-7698f353eb3d
begin 
	import QuadGK: quadgk
	using CairoMakie
	using DataFrames, CSV
	using NaNMath; nm = NaNMath
	using Arya
	using HDF5
	import DensityEstimators
end

# ╔═╡ 631a70f3-5284-4c3f-81ef-714455b876ee
using FITSIO

# ╔═╡ 17ffde4b-5796-4915-9741-d594cf0c5ca7
md"""
# Paint Stars

My version of Rapha's codes. The main idea is to use Eddington inversion to calculate the energy distribution function. This allows the relative densities of stars/dark matter to be calculated for each particle, resulting in assigned probabilities for a stellar profile.


The script requires an input toml config file to run. See below for parameter descriptions.
"""

# ╔═╡ b3663c19-c029-4a1e-ab82-a05177e3a5d0
import StatsBase: percentile

# ╔═╡ 530c6c09-4454-4952-8351-dccbb4ed429f
import TOML

# ╔═╡ 93045024-a91d-4b31-9a5a-7c999afdb9ec
md"""
# Inputs
"""

# ╔═╡ da38359e-f138-41ee-9b8f-9cc8a38f77a2
md"""

Parameters
- `profile`: should be a string implementing one of the classes in lguys analytic profiles (Exp2D, Exp3D, ...)
- `num_radial_bins`: int. Number of radial bins to divide particles into
- `num_energy_bins`: number of energy bins at which to sample the distribution function.
- `snapshot_dir` relative directory of all snapshot files from gadget
- `snapshot`: the integer number of the snapshot we would like to load
- `bin_method`: the binning method for radial bins. Should be "equal\_width" or "equal\_number"
"""

# ╔═╡ 48ce69f2-09d5-4166-9890-1ab768f3b59f
# input directory
dir = "/astro/dboyea/sculptor/isolation/1e6/stars"

# ╔═╡ 7809e324-ba5f-4520-b6e4-c7727c227154
paramname = joinpath(dir, "exp2d_rs0.1")

# ╔═╡ 8a8d3180-ab8c-4456-a2cf-6ffb7dc73760
params = TOML.parsefile(paramname * ".toml")

# ╔═╡ 0ede2af5-a572-41c8-b3f0-cb0a24318c5f
profile = lguys.load_profile(paramname * ".toml")

# ╔═╡ 159aa60d-563b-4c85-b75d-502ac944b99c
df_E = lguys.load_fits(paramname * "_energy.fits")

# ╔═╡ f1a7fa1f-bdcd-408c-ab27-b52916b1892f
df_density = lguys.load_fits(paramname * "_density.fits")

# ╔═╡ 1066a445-600d-4508-96a2-aa9b90460097
df_probs = lguys.load_hdf5_table(paramname * "_stars.hdf5")

# ╔═╡ 578c6196-db59-4d5c-96a7-9a8487bbeaae
begin 
	snap = lguys.Snapshot(joinpath(dir, params["snapshot"]))
	snap.weights = df_probs.probability[snap.index]
	snap
end

# ╔═╡ 29930595-5255-4454-8550-22ac6a96f609
r_h = lguys.get_r_h(profile)

# ╔═╡ 4d1991ea-9496-48c7-a400-8fefbecefcd2
md"""
# Plots
"""

# ╔═╡ 5b30475b-b4c4-4c87-817d-0d885546d004
md"""
## Distributions and intermediate quantities
"""

# ╔═╡ 2aa9c02e-621e-4ba2-b69e-a9ed66d834b8
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log radii / kpc", ylabel="counts")
	h = lguys.histogram(log10.(df_probs.radii), log10.(df_density.r), normalization=:none)
	scatter!(h)

	fig
end

# ╔═╡ a5bc5ce3-8e33-4514-bc2d-4b4299f104f9
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log radii", ylabel="PDF")
	stephist!(log10.(df_probs.radii), bins=log10.(df_density.r), normalization=:pdf, label="dark matter")
	stephist!(log10.(df_probs.radii), bins=log10.(df_density.r), weights=df_probs.probability, normalization=:pdf, label="stars (nbody)")
	axislegend()
	fig
end

# ╔═╡ 84fdc265-988c-40db-87e5-44ba55d0e412
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"$\epsilon$ (binding energy)", 
		ylabel="asinh density", 
		)

	scale = 1e-10

	y_trans(x) =  asinh.(x / scale)
	ϵs = df_probs.eps
	ps = df_probs.probability

	
	h1 = lguys.histogram(ϵs, 20, normalization=:pdf)
	h1.values .= y_trans(h1.values)

	h_s = lguys.histogram(ϵs, 20, weights=ps, normalization=:pdf)
	h_s.values .= y_trans(h_s.values)

	
	lines!(df_E.E, df_E.probs, color=Arya.COLORS[3], label="f_s / f_dm")
	
	lines!(h1, label="dark matter")
	lines!(h_s, label="stars (nbody)")
	vlines!([maximum(ϵs)], color="grey", linestyle=:dot, label=L"\epsilon_\textrm{max}")

	axislegend(ax, position=:lb)

	fig
end

# ╔═╡ 8bb8736d-a41b-4dac-a6cd-06d0d4704654
let
	fig = Figure()
	ax = Axis(fig[1,1], limits=(-1.2, 3, -15, 2),
		xlabel="log r", ylabel="log density")
	scatter!(log10.(df_density.r), log10.(df_density.nu_dm), label="DM")
	scatter!(log10.(df_density.r), log10.(df_density.nu_s), label="stars (analytic)")

	axislegend()
	fig
end

# ╔═╡ 75d23b44-71e7-4e28-ad3e-c537f3d4422f
let
	fig = Figure()
	ax = Axis(fig[1,1],xlabel="log ϵ", ylabel="log f", limits=(nothing, (-15, 7)) )
	lines!(log10.(df_E.E), nm.log10.(df_E.f_s_e), label="stars")
	lines!(log10.(df_E.E), nm.log10.(df_E.f_dm_e), label="DM")

	axislegend(ax, position=:lt)
	fig
end

# ╔═╡ 9e2f1606-46aa-4e06-a31f-b03a383cccda
md"""
The calculated energy distribution function for both stars and dark matter (sampled at the specified number of points). Ideally, we would like this curve to be smooth and well-sampled.
"""

# ╔═╡ b625d8a5-7265-4849-9bd6-ca8064d392eb
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"\log\ r / \textrm{kpc}", ylabel=L"\Psi(r)")

	idx = sortperm(df_probs.radii)
	lines!(log10.(df_probs.radii[idx]), df_probs.phi[idx], label="interpolated")

	idx = sortperm(lguys.calc_r(snap))
	lines!(log10.(lguys.calc_r(snap)[idx]), snap.Φs[idx], label="snapshot")
	axislegend()
	fig
end

# ╔═╡ 999df0b7-1ff0-4771-8113-2bfe7a74b646
md"""
# Reconstruction plots
"""

# ╔═╡ ffca9cd9-a2d7-4f52-8467-8f4757ddf445
md"""
For self consistancy, here is the gravitational potential as ccompared to the Gadget4 calculations.
"""

# ╔═╡ 76200404-16aa-4caf-b247-3bc330b82868
function calc_r_h(rs, masses)
	_sort = sortperm(rs)
	M_in = cumsum(masses[_sort])
	M_in ./= M_in[end]
	idx_h = findfirst(M_in .> 0.5)
	return rs[idx_h]
end

# ╔═╡ 65c7cb79-1246-493d-b6d6-4d836746e396
df_density.r

# ╔═╡ 406c942e-d5a5-4ab7-806c-27ce79cefffb
r_e = [df_density.r[1] - df_density.dr[1]; df_density.r .+ df_density.dr]

# ╔═╡ 6fba7fa7-9a50-4379-b376-5c07f3638411
ν_s_nbody = lguys.calc_ρ_hist(df_probs.radii, r_e, weights=df_probs.probability)[2]

# ╔═╡ c5fa22bf-6090-4014-85e6-6681b2bb10a7
df_density.r .- df_density.dr .- r_e[1:end-1]

# ╔═╡ a9335e17-a410-455a-9a9e-d63706a026bd
let
	fig = Figure(size=(700, 500))
	ax = Axis(fig[1,1], ylabel=L"\log \nu", 
		limits=((0, 1), (-15, 3))
		)
	lines!(log10.(df_density.r), nm.log10.(df_density.nu_s), label="stars")
	scatter!(log10.(df_density.r) , nm.log10.(ν_s_nbody), label="nbody", color=COLORS[2])

	ax2 = Axis(fig[2,1], 
		xlabel=L"\log\,r / \textrm{kpc}", ylabel=L"\Delta\log \nu ", 
		limits=((log10(0.5r_h), log10(100r_h)), (-1, 1)))
	
	scatter!(log10.(df_density.r), nm.log10.(ν_s_nbody) .- nm.log10.(df_density.nu_s), label="",
		color=COLORS[2]
	)
	hlines!([0], label="")

	linkxaxes!(ax, ax2, )
	rowsize!(fig.layout, 2, Auto(0.3))
	hidexdecorations!(ax, grid=false)
	fig
end

# ╔═╡ 33a26663-0c08-411b-902b-a509b2afa5ad
let
	fig = Figure()
	Axis(fig[1,1], xlabel="log radii", ylabel="pstar > 0.025")

	hist!(log10.(lguys.calc_r(snap)[snap.weights .> 2e-6]))

	fig
end

# ╔═╡ 77e2c1e3-7756-4ab7-810a-03ccdc635aa1
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="stellar weights", ylabel="frequency", yscale=log10)
	hist!(snap.weights, bins=100)
	fig
end

# ╔═╡ 90856551-fff8-4d15-be66-6353091b5e50
begin 
	r_hist = 5
	N_hist = 50
end

# ╔═╡ a52e5e94-6068-4545-962f-e02a485b62f5
let
	fig = Figure()
	ax = Axis(fig[1,1],
		limits=(-r_hist, r_hist, -r_hist, r_hist), 
		aspect=1,
		title="dark matter",
		xlabel="x / kpc",
		ylabel="y / kpc"
	)
	
	h = Arya.hist2d!(ax, 
		lguys.get_x(snap), lguys.get_y(snap), 
		bins=N_hist,
		colorscale=log10,
		colorrange=(1e-2, nothing)
	)	

	#scatter!(cen.position[1], cen.position[2])

	Colorbar(fig[1, 2],h )

	fig
end

# ╔═╡ cc231e78-cfc0-4876-9f8b-980139f7d27f
let
	fig = Figure()
	filt = snap.weights .> 0
	
	ax = Axis(fig[1,1], 
		limits=(-r_hist, r_hist, -r_hist, r_hist), 
		aspect=1,
		title="stars",
		xlabel="x / kpc",
		ylabel="y / kpc"
	)
	
	h = Arya.hist2d!(ax, 
		lguys.get_x(snap)[filt], lguys.get_y(snap)[filt], 
		weights=snap.weights[filt], bins=N_hist,
		colorscale=log10,
		colorrange=(1e-7, nothing)
	)	

	#scatter!(cen.position[1], cen.position[2])

	Colorbar(fig[1, 2], h)
	fig
end

# ╔═╡ 89b7d969-2294-4d07-a6a3-fcfa92498d48
lguys.arcmin_to_kpc(14, 83.2)

# ╔═╡ 7f7d8cb9-761c-4f30-a336-ab5657144961
let
	r = lguys.calc_r(snap)
	ms = snap.weights

	
	r_h2 = calc_r_h(r, ms)
	println(r_h2)
	
	r_e = 10 .^ DensityEstimators.bins_min_width_equal_number(log10.(r), N_per_bin_min=100, dx_min=0.1)
	r, ν_s_nbody = lguys.calc_ρ_hist(r, r_e, weights=ms)
 
	r = lguys.midpoints(r)
	ν_s = lguys.calc_ρ.(profile, r)
	
	fig = Figure(size=(700, 500))
	ax = Axis(fig[1,1], ylabel=L"\log \nu", 
		limits=(nothing, (-15, 3))
		)
	
	lines!(log10.(r) , nm.log10.(ν_s_nbody), label="nbody")
	lines!(log10.(r), nm.log10.(ν_s), label="stars")
	vlines!(log10(r_h), label="R_h")
	vlines!(log10(r_h2), label="R_h2")

	ax2 = Axis(fig[2,1], 
		xlabel=L"\log r / \textrm{kpc}", ylabel=L"\Delta\log \nu ", 
		limits=((log10(0.01r_h), log10(100r_h)), (-1, 1)))


	
	scatter!(log10.(r), nm.log10.(ν_s_nbody) .- nm.log10.(ν_s), label="")
	hlines!([0], color="black", label="")

	linkxaxes!(ax, ax2, )
	rowsize!(fig.layout, 2, Auto(0.3))
	hidexdecorations!(ax, grid=false)
	fig
end

# ╔═╡ a2f72082-7145-42be-9f40-e00d18deb267


# ╔═╡ bf8305f4-a5b8-4c79-8a01-e2aa18e4a5c5
md"""
## Testing 2D binning
"""

# ╔═╡ 6dd92ee1-374d-47fa-ad61-b54764b23240
let
	x = snap.positions[1, :] 
	y = snap.positions[2, :] 
	R = @. sqrt(x^2 + y^2)

	bins = DensityEstimators.bins_min_width_equal_number(log10.(R), N_per_bin_min=100, dx_min=0.1)

	prof = lguys.calc_properties(R, weights=snap.weights, bins=bins)

	
	fig = Figure()
	ax = Axis(fig[1,1], ylabel=L"\log\ \Sigma", 
		xlabel="log R (2D, xy) / kpc",
		limits=(nothing, (-15, 3))
		)


	profile2 = lguys.Exp2D(R_s = profile.R_s)

	log_Σ(r) = log10(lguys.calc_Σ(profile2, r))

	log_R = LinRange(-2, 2, 1000)
	y = log_Σ.(10 .^ log_R)
	
	errscatter!(prof.log_r, prof.log_Sigma, yerr=prof.log_Sigma_err)
	lines!(log_R, y)

	fig
end

# ╔═╡ 99f274b4-91f3-41d0-b7d3-234badeb43d1
md"""
## and on the sky
"""

# ╔═╡ 4396bced-cae8-4107-ac83-48cc3c4146f2
distance = 83.2

# ╔═╡ 062cea41-3b51-474e-b877-7a4d96813fbc
R_s_arcmin = lguys.kpc_to_arcmin(profile.R_s, distance)

# ╔═╡ 7b7832d1-6739-453e-aaee-fa16a6000d26
function mock_obs(snap_og; distance=86)
	x_sun = [8.122, 0, 0]
	shift_vec = x_sun .+ lguys.rand_unit() * distance
	
	ps = snap.weights
	p_min = 1e-25
	filt = p_min * maximum(ps) .< ps

	ps = ps[filt]
	snap_shift = deepcopy(snap_og[filt])
	snap_shift.positions .+= shift_vec

	snap_shift.weights = ps
	obs = lguys.to_sky(snap_shift)
	return obs
end

# ╔═╡ 87b5b241-db72-45ee-b3a7-a394f99510d9
obs_mock = mock_obs(snap)

# ╔═╡ ff51f97d-2404-49c6-9339-4b201d6a94a9
let
	ra = obs_mock.ra
	dec = obs_mock.dec
	ms = obs_mock.weights

	ra0, dec0 = lguys.calc_centre2D(ra, dec, "mean", ms)
	xi, eta = lguys.to_tangent(ra, dec, ra0, dec0)
	R = @. 60sqrt(xi^2 + eta^2)
	bins = 	 DensityEstimators.bins_min_width_equal_number(log10.(R), N_per_bin_min=30, dx_min=0.1)

	prof = lguys.calc_properties(R, weights=ms, bins=bins)


	fig = Figure()
	ax = Axis(fig[1,1], ylabel=L"\log\ \Sigma", 
		xlabel="R (projected at $distance kpc )/ arcmin",
		limits=((-1, 2.5), (-15, 3))
		)

	errscatter!(prof.log_r, prof.log_Sigma, yerr=prof.log_Sigma_err)

	
	profile2 = lguys.Exp2D(R_s = R_s_arcmin)

	log_Σ(r) = log10(lguys.calc_Σ(profile2, r))

	log_R = LinRange(-2, 2.5, 1000)
	y = log_Σ.(10 .^ log_R)
	
	lines!(log_R, y)
	fig
end

# ╔═╡ 03037c94-7f26-479c-9772-fa2682e3ba37
import StatsBase: mean, std, weights

# ╔═╡ 74845a4f-588b-44d3-a5e1-7155fd6a7b01
let
	x = obs_mock.radial_velocity
	m = obs_mock.weights
	μ = mean(x, weights(m))
	σ = std(x, weights(m))
	println(μ)
	println(σ)

	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)

	stephist!(x, weights=m, normalization=:pdf, bins=50, label="n-body stellar density")
	x_model = LinRange(μ - 5σ, μ + 5σ, 1000)
	y_model = lguys.gaussian.(x_model, μ, σ)

	lines!(x_model, y_model, color=COLORS[2], label="N($(round(μ, digits=2)), $(round(σ, digits=2)))")

	axislegend()
	
	fig
end

# ╔═╡ e4c33671-744b-42bb-87be-2df7add03b96
let
	x = lguys.calc_r(snap.velocities) .^2 * lguys.V2KMS^2
	m = snap.weights
	
	μ = mean(x, weights(m))
	println(sqrt(μ)/ sqrt(3))
end

# ╔═╡ b89fbbb0-b161-4fa7-82ba-c06e3d32f8b8
function save_obs(obs_df, outfile)
	rm(outfile, force=true)
	df = Dict(String(name) => obs_df[:, name] for name in names(obs_df))

	FITS(outfile, "w") do f
		write(f, df)
		println("written to $outfile")

		df
	end
end

# ╔═╡ d29baea8-0e5c-43e0-9aa2-987d6cc08e88
save_obs(obs_mock, params["mock_file"])

# ╔═╡ Cell order:
# ╟─17ffde4b-5796-4915-9741-d594cf0c5ca7
# ╠═a893932c-f184-42bc-9a0e-0960f10520aa
# ╠═641946b3-e6f2-4d6d-8777-7698f353eb3d
# ╠═b3663c19-c029-4a1e-ab82-a05177e3a5d0
# ╠═530c6c09-4454-4952-8351-dccbb4ed429f
# ╠═631a70f3-5284-4c3f-81ef-714455b876ee
# ╟─93045024-a91d-4b31-9a5a-7c999afdb9ec
# ╟─da38359e-f138-41ee-9b8f-9cc8a38f77a2
# ╠═48ce69f2-09d5-4166-9890-1ab768f3b59f
# ╠═7809e324-ba5f-4520-b6e4-c7727c227154
# ╠═8a8d3180-ab8c-4456-a2cf-6ffb7dc73760
# ╠═578c6196-db59-4d5c-96a7-9a8487bbeaae
# ╠═0ede2af5-a572-41c8-b3f0-cb0a24318c5f
# ╠═159aa60d-563b-4c85-b75d-502ac944b99c
# ╠═f1a7fa1f-bdcd-408c-ab27-b52916b1892f
# ╠═1066a445-600d-4508-96a2-aa9b90460097
# ╠═29930595-5255-4454-8550-22ac6a96f609
# ╟─4d1991ea-9496-48c7-a400-8fefbecefcd2
# ╟─5b30475b-b4c4-4c87-817d-0d885546d004
# ╠═2aa9c02e-621e-4ba2-b69e-a9ed66d834b8
# ╠═a5bc5ce3-8e33-4514-bc2d-4b4299f104f9
# ╠═84fdc265-988c-40db-87e5-44ba55d0e412
# ╠═8bb8736d-a41b-4dac-a6cd-06d0d4704654
# ╠═75d23b44-71e7-4e28-ad3e-c537f3d4422f
# ╟─9e2f1606-46aa-4e06-a31f-b03a383cccda
# ╠═b625d8a5-7265-4849-9bd6-ca8064d392eb
# ╟─999df0b7-1ff0-4771-8113-2bfe7a74b646
# ╟─ffca9cd9-a2d7-4f52-8467-8f4757ddf445
# ╠═76200404-16aa-4caf-b247-3bc330b82868
# ╠═6fba7fa7-9a50-4379-b376-5c07f3638411
# ╠═65c7cb79-1246-493d-b6d6-4d836746e396
# ╠═406c942e-d5a5-4ab7-806c-27ce79cefffb
# ╠═c5fa22bf-6090-4014-85e6-6681b2bb10a7
# ╠═a9335e17-a410-455a-9a9e-d63706a026bd
# ╠═33a26663-0c08-411b-902b-a509b2afa5ad
# ╠═77e2c1e3-7756-4ab7-810a-03ccdc635aa1
# ╠═90856551-fff8-4d15-be66-6353091b5e50
# ╠═a52e5e94-6068-4545-962f-e02a485b62f5
# ╠═cc231e78-cfc0-4876-9f8b-980139f7d27f
# ╠═89b7d969-2294-4d07-a6a3-fcfa92498d48
# ╠═7f7d8cb9-761c-4f30-a336-ab5657144961
# ╠═a2f72082-7145-42be-9f40-e00d18deb267
# ╟─bf8305f4-a5b8-4c79-8a01-e2aa18e4a5c5
# ╠═6dd92ee1-374d-47fa-ad61-b54764b23240
# ╟─99f274b4-91f3-41d0-b7d3-234badeb43d1
# ╠═4396bced-cae8-4107-ac83-48cc3c4146f2
# ╠═062cea41-3b51-474e-b877-7a4d96813fbc
# ╠═87b5b241-db72-45ee-b3a7-a394f99510d9
# ╠═7b7832d1-6739-453e-aaee-fa16a6000d26
# ╠═ff51f97d-2404-49c6-9339-4b201d6a94a9
# ╠═03037c94-7f26-479c-9772-fa2682e3ba37
# ╠═74845a4f-588b-44d3-a5e1-7155fd6a7b01
# ╠═e4c33671-744b-42bb-87be-2df7add03b96
# ╠═b89fbbb0-b161-4fa7-82ba-c06e3d32f8b8
# ╠═d29baea8-0e5c-43e0-9aa2-987d6cc08e88
