### A Pluto.jl notebook ###
# v0.20.5

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

# ╔═╡ a893932c-f184-42bc-9a0e-0960f10520aa
begin
	import Pkg
	Pkg.activate()
	
	import LilGuys as lguys
end

# ╔═╡ 1ff14ee1-3f98-4144-9493-7479e56e8790
using PlutoUI

# ╔═╡ 641946b3-e6f2-4d6d-8777-7698f353eb3d
begin 
	using CairoMakie
	using NaNMath; nm = NaNMath
	using Arya
	using HDF5
	import DensityEstimators
end

# ╔═╡ 631a70f3-5284-4c3f-81ef-714455b876ee
using FITSIO

# ╔═╡ 17ffde4b-5796-4915-9741-d594cf0c5ca7
md"""
# Paint Stars Plots

This script takes the outputs from the script `paint_stars.jl` and plots the results for consistency checks.

"""

# ╔═╡ 81dc2694-56ba-496a-85d1-cf2ecbaa2e67
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

# ╔═╡ 420977eb-31df-44f3-a2eb-3209d66b10c2
@bind inputs confirm(notebook_inputs(
	modelname = TextField(60, default="ursa_minor/1e6_v3"),
	starsname = TextField(default="exp2d_rs0.10"),
))

# ╔═╡ 81bf8451-f417-4e37-a234-07bb5317fff1
import DensityEstimators: histogram

# ╔═╡ b3663c19-c029-4a1e-ab82-a05177e3a5d0
import StatsBase: percentile

# ╔═╡ 530c6c09-4454-4952-8351-dccbb4ed429f
import TOML

# ╔═╡ 48ce69f2-09d5-4166-9890-1ab768f3b59f
#dir = "/astro/dboyea/dwarfs/analysis/sculptor/1e7_V31_r3.2/stars/"
dir = joinpath(ENV["DWARFS_ROOT"], "analysis", inputs.modelname, "stars")

# ╔═╡ 939cc89e-7273-4bb5-a13f-241139d922ea
starsname = inputs.starsname

# ╔═╡ 7809e324-ba5f-4520-b6e4-c7727c227154
paramname = joinpath(dir, starsname, "profile.toml")

# ╔═╡ 79f224cf-198c-4527-8423-62d74c753d0a
FIGDIR = joinpath(dir, starsname, "figures")

# ╔═╡ 835cafca-868b-4e94-955f-2eab4e7c2bc4
import LilGuys: @savefig; FIGDIR

# ╔═╡ 66c89b32-e522-48c8-b5aa-a5845715ce28
if !isdir(FIGDIR)
	mkdir(FIGDIR)
end

# ╔═╡ d76e6200-9401-4c2e-bd7c-53e79dd49415
md"""
# File loading
"""

# ╔═╡ 0ede2af5-a572-41c8-b3f0-cb0a24318c5f
profile = lguys.load_profile(paramname)

# ╔═╡ 715f771b-686a-4643-a914-35ba6ca9042d
df_E = lguys.read_hdf5_table(joinpath(dir, starsname, "probabilities_df.hdf5"))

# ╔═╡ 1f10e98c-42ff-4598-9822-3f0af5c065bb
df_E[8000, :]

# ╔═╡ 1066a445-600d-4508-96a2-aa9b90460097
df_probs = lguys.read_hdf5_table(joinpath(dir, starsname, "probabilities_stars.hdf5"))

# ╔═╡ 578c6196-db59-4d5c-96a7-9a8487bbeaae
begin 
	snap = lguys.Snapshot(joinpath(dir, "iso_paint.hdf5"))
	snap.weights = df_probs.probability[snap.index]
	# cen = lguys.Centres.shrinking_spheres(snap.positions)[1]
	# snap.x_cen = cen
	snap.positions
end

# ╔═╡ ed0996a4-8667-49f7-b89f-64ae44f8a3ea
snap.x_cen

# ╔═╡ 0a4521ac-7e35-4976-8781-bdbd4f7242c7
lguys.get_M_tot(profile)

# ╔═╡ 29930595-5255-4454-8550-22ac6a96f609
r_h = lguys.calc_r_h(profile)

# ╔═╡ f5582e2e-6cbf-4b32-9da0-86b4f33c55b6
bins = LinRange(log10(minimum(df_probs.radii)), log10(maximum(df_probs.radii)), 100)

# ╔═╡ 4d1991ea-9496-48c7-a400-8fefbecefcd2
md"""
# Plots
"""

# ╔═╡ 5b30475b-b4c4-4c87-817d-0d885546d004
md"""
## Distributions and intermediate quantities
"""

# ╔═╡ a5bc5ce3-8e33-4514-bc2d-4b4299f104f9
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log radii", ylabel="PDF")
	stephist!(log10.(df_probs.radii), bins=bins, normalization=:pdf, label="dark matter")
	stephist!(log10.(df_probs.radii), bins=bins, weights=df_probs.probability, normalization=:pdf, label="stars (nbody)")
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

	
	bins, h1, _ = lguys.histogram(ϵs, 200, normalization=:pdf)
	h1 = y_trans(h1)

	_, h_s, _ = lguys.histogram(ϵs, 200, weights=ps, normalization=:pdf)
	h_s .= y_trans(h_s)

	
	lines!(df_E.psi, asinh.(df_E.probability), color=Arya.COLORS[3], label="f_s / f_dm")
	
	lines!(midpoints(bins), h1, label="dark matter")
	lines!(midpoints(bins), h_s, label="stars (nbody)")
	vlines!([maximum(ϵs)], color="grey", linestyle=:dot, label=L"\epsilon_\textrm{max}")

	axislegend(ax, position=:lb)

	fig
end

# ╔═╡ 6da679d4-6af6-4f42-b6e9-44ce20faa676
let
	fig = Figure()
	ax = Axis(fig[1,1],xlabel="ϵ", ylabel="f", 
	)

	h = 1e-5
	
	lines!((df_E.psi), asinh.(df_E.f ./ h), label="DM")
	lines!((df_E.psi), asinh.(df_E.f_s ./ h), label="stars")

	axislegend(ax, position=:lt)
	fig
end

# ╔═╡ d8623523-ffcc-4ceb-aafb-fe0cde23f469
md"""
Below shows the functional form of the stellar probabilities as a function of energy. We expect tboth curves to match almost perfectly
"""

# ╔═╡ a462157f-e991-4ee0-9a62-c0bb3c687e87
let
	fig = Figure()
	ax = Axis(fig[1,1],xlabel="ϵ", ylabel="p(e)", 
	)

	h = 0.01
	
	lines!((df_E.psi), asinh.(df_E.probability ./ h), label="interpolated")
	lines!((df_E.psi), asinh.(df_E.f_s ./ df_E.f ./ h), label="ratio")

	axislegend(ax, position=:lt)
	fig
end

# ╔═╡ 49f66f52-5c29-49d1-9e50-2be48f494ceb
maximum(df_E.psi)

# ╔═╡ 061cccde-82c1-4519-aec9-4fc92a67b348
maximum(df_probs.eps)

# ╔═╡ 9e2f1606-46aa-4e06-a31f-b03a383cccda
md"""
The calculated energy distribution function for both stars and dark matter (sampled at the specified number of points). Ideally, we would like this curve to be smooth and well-sampled.
"""

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
	return rs[_sort][idx_h]
end

# ╔═╡ bec1ed0b-7a3d-4523-84ad-3fad4a7296cd
bins

# ╔═╡ 6fba7fa7-9a50-4379-b376-5c07f3638411
ν_s_nbody = lguys.calc_ρ_from_hist(10 .^ bins, histogram(df_probs.radii, 10 .^ bins, weights=df_probs.probability).values)

# ╔═╡ 9e28a6e0-4375-414c-bbe4-287d4a13e22d
md"""
The below plot simply compares the simperical DM and stellar profiles
"""

# ╔═╡ a9335e17-a410-455a-9a9e-d63706a026bd
let
	fig = Figure(size=(700, 500))
	ax = Axis(fig[1,1], ylabel=L"\log \nu", 
		limits=((-1.5, 1), (-15, 3)),
		xlabel = "log r / kpc",
		)

	log_r_m = midpoints(bins)
	lines!(log10.(df_E.radii), nm.log10.(df_E.rho), label="DM profile")
	scatter!((log_r_m) , nm.log10.(ν_s_nbody), label="nbody stellar profile", color=COLORS[2])

	axislegend()
	fig
end

# ╔═╡ b5c587d4-df0b-451b-9e08-bd4a4ead0213
md"""
Below, we show the distribution in log radius / kpc for stars with non-zero stellar weights. Ideally, we are looking for a smooth distribution without huge outliers in pstar
"""

# ╔═╡ 33a26663-0c08-411b-902b-a509b2afa5ad
let
	fig = Figure()

	thresh = 1 / length(snap)

	Axis(fig[1,1], xlabel="log radii", ylabel="pstar > $thresh")

	hist!(log10.(lguys.calc_r(snap)[snap.weights .> thresh]))

	fig
end

# ╔═╡ 2d1782c5-e544-4ee6-a872-c8a180bd0d72
md"""
Distribution of stellar weights. Would like this distribution to be continuous and without outliers towards higher stellar weights. Most particles will have tiny weights.
"""

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

# ╔═╡ 337055c4-55f4-4c02-ad69-473165ccfb64
md"""
Below, we plot the 2D density of DM and stars from the model, just to make sure nothing funny be happening. In both cases, the distribution should be spherical and continuous.
"""

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

# ╔═╡ cdcc6bf6-f15e-45f1-b837-1db4e0fdc654
md"""
This plot is the primary evaluation of how well the reconstructed stellar profile matches our expected version.
We would like the residuals to be << 1 and the profile to qualitatively match well. At large radii, it is acceptable that the residuals diverge because the stellar mass may be entirely negligable here.
"""

# ╔═╡ 7f7d8cb9-761c-4f30-a336-ab5657144961
let
	r = lguys.calc_r(snap)
	ms = snap.weights 

	
	r_h2 = calc_r_h(r, ms)
	println("rh = ", r_h)
	println("rh n-body = ", r_h2)
	
	r_e = 10 .^ DensityEstimators.bins_min_width_equal_number(log10.(r), N_per_bin_min=100, dx_min=0.03)
	ν_s_nbody = lguys.calc_ρ_from_hist(r_e, histogram(r, r_e, weights=ms).values)
 
	r = lguys.midpoints(r_e)
	ν_s = lguys.calc_ρ.(profile, r)
	
	fig = Figure(size=(700, 500))
	ax = Axis(fig[1,1], ylabel=L"\log \nu", 
		limits=(nothing, (-15, 3))
		)
	
	lines!(log10.(r), nm.log10.(ν_s), label="stars", color=:black)
	
	lines!(log10.(r) , nm.log10.(ν_s_nbody), label="nbody", color=COLORS[2])
	vlines!(log10(r_h2), label=L"R_h")

	axislegend(position=:lb)

	ax2 = Axis(fig[2,1], 
		xlabel=L"\log r / \textrm{kpc}", ylabel=L"\Delta\log \nu ", 
		limits=((log10(0.01r_h), log10(100r_h)), (-1, 1)))


	
	scatter!(log10.(r), nm.log10.(ν_s_nbody) .- nm.log10.(ν_s), color=COLORS[2], label="")
	hlines!([0], color="black", label="")

	linkxaxes!(ax, ax2, )
	rowsize!(fig.layout, 2, Auto(0.3))
	hidexdecorations!(ax, grid=false)

	@savefig "density_check"
	fig
end

# ╔═╡ a2f72082-7145-42be-9f40-e00d18deb267
sum(snap.weights .* (lguys.calc_r(snap) .< r_h))

# ╔═╡ bf8305f4-a5b8-4c79-8a01-e2aa18e4a5c5
md"""
## Testing 2D binning
"""

# ╔═╡ 2751e203-8c03-44ec-a072-e3078d9a4f7b
md"""
Similar to the last plot but for 2D binning. Just validating qualitative agreement.
"""

# ╔═╡ 6dd92ee1-374d-47fa-ad61-b54764b23240
let
	x = snap.positions[3, :] 
	y = snap.positions[1, :] 
	R = @. sqrt(x^2 + y^2)

	bins = DensityEstimators.bins_min_width_equal_number(log10.(R), N_per_bin_min=100, dx_min=0.03)

	prof = lguys.StellarProfile(R, weights=snap.weights, bins=bins)

	
	fig = Figure()
	ax = Axis(fig[1,1], ylabel=L"\log\ \Sigma", 
		xlabel="log R (2D, xy) / kpc",
		limits=(nothing, (-15, 3))
		)


	log_Σ(r) = log10(lguys.calc_Σ(profile, r))

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

# ╔═╡ e63b488f-d1d7-488f-a1b8-6c32d2f66829
snap.x_cen, snap.positions

# ╔═╡ 87b5b241-db72-45ee-b3a7-a394f99510d9
obs_mock = lguys.to_gaia(snap, set_to_distance=distance, add_centre=false)

# ╔═╡ 904576b4-2669-45d1-8078-310347f5fec0
sum(.!isnan.(obs_mock.r_ell))

# ╔═╡ 8561415c-7c5d-4c1a-8478-3b237c2cb360
obs_mock

# ╔═╡ ff51f97d-2404-49c6-9339-4b201d6a94a9
let
	ms = obs_mock.weights
	R = obs_mock.r_ell
	
	bins = 	DensityEstimators.bins_min_width_equal_number(log10.(R), N_per_bin_min=5, dx_min=0.03)
	@info "using $(length(bins)) bins"

	prof = lguys.StellarProfile(R, weights=ms, bins=bins, normalization=:central, r_centre=3)


	fig = Figure()
	ax = Axis(fig[1,1], ylabel=L"\log\ \Sigma", 
		xlabel="R (projected at $distance kpc )/ arcmin",
		limits=((-1, 3), (-15, 3))
		)

	errscatter!(prof.log_r, prof.log_Sigma, yerr=prof.log_Sigma_err)

	
	log_Σ(r) = log10(lguys.calc_Σ(profile, lguys.arcmin_to_kpc(r, distance)))

	log_R = LinRange(-2, 3, 1000)
	y = log_Σ.(10 .^ log_R)
	y .-= y[1]
	
	lines!(log_R, y)
	fig
end

# ╔═╡ 03037c94-7f26-479c-9772-fa2682e3ba37
import StatsBase: mean, std, weights

# ╔═╡ ad02a484-11b1-479f-b9e0-0e16e64a5109
md"""
Velocity dispersion!
"""

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

# ╔═╡ Cell order:
# ╟─17ffde4b-5796-4915-9741-d594cf0c5ca7
# ╟─420977eb-31df-44f3-a2eb-3209d66b10c2
# ╠═81dc2694-56ba-496a-85d1-cf2ecbaa2e67
# ╠═1ff14ee1-3f98-4144-9493-7479e56e8790
# ╠═a893932c-f184-42bc-9a0e-0960f10520aa
# ╠═641946b3-e6f2-4d6d-8777-7698f353eb3d
# ╠═81bf8451-f417-4e37-a234-07bb5317fff1
# ╠═b3663c19-c029-4a1e-ab82-a05177e3a5d0
# ╠═530c6c09-4454-4952-8351-dccbb4ed429f
# ╠═631a70f3-5284-4c3f-81ef-714455b876ee
# ╠═835cafca-868b-4e94-955f-2eab4e7c2bc4
# ╠═48ce69f2-09d5-4166-9890-1ab768f3b59f
# ╠═939cc89e-7273-4bb5-a13f-241139d922ea
# ╠═7809e324-ba5f-4520-b6e4-c7727c227154
# ╠═79f224cf-198c-4527-8423-62d74c753d0a
# ╠═66c89b32-e522-48c8-b5aa-a5845715ce28
# ╟─d76e6200-9401-4c2e-bd7c-53e79dd49415
# ╠═0ede2af5-a572-41c8-b3f0-cb0a24318c5f
# ╠═715f771b-686a-4643-a914-35ba6ca9042d
# ╠═1f10e98c-42ff-4598-9822-3f0af5c065bb
# ╠═578c6196-db59-4d5c-96a7-9a8487bbeaae
# ╠═ed0996a4-8667-49f7-b89f-64ae44f8a3ea
# ╠═1066a445-600d-4508-96a2-aa9b90460097
# ╠═0a4521ac-7e35-4976-8781-bdbd4f7242c7
# ╠═29930595-5255-4454-8550-22ac6a96f609
# ╠═f5582e2e-6cbf-4b32-9da0-86b4f33c55b6
# ╟─4d1991ea-9496-48c7-a400-8fefbecefcd2
# ╟─5b30475b-b4c4-4c87-817d-0d885546d004
# ╠═a5bc5ce3-8e33-4514-bc2d-4b4299f104f9
# ╠═84fdc265-988c-40db-87e5-44ba55d0e412
# ╠═6da679d4-6af6-4f42-b6e9-44ce20faa676
# ╠═d8623523-ffcc-4ceb-aafb-fe0cde23f469
# ╠═a462157f-e991-4ee0-9a62-c0bb3c687e87
# ╠═49f66f52-5c29-49d1-9e50-2be48f494ceb
# ╠═061cccde-82c1-4519-aec9-4fc92a67b348
# ╟─9e2f1606-46aa-4e06-a31f-b03a383cccda
# ╟─999df0b7-1ff0-4771-8113-2bfe7a74b646
# ╟─ffca9cd9-a2d7-4f52-8467-8f4757ddf445
# ╠═76200404-16aa-4caf-b247-3bc330b82868
# ╠═bec1ed0b-7a3d-4523-84ad-3fad4a7296cd
# ╠═6fba7fa7-9a50-4379-b376-5c07f3638411
# ╟─9e28a6e0-4375-414c-bbe4-287d4a13e22d
# ╟─a9335e17-a410-455a-9a9e-d63706a026bd
# ╠═b5c587d4-df0b-451b-9e08-bd4a4ead0213
# ╠═33a26663-0c08-411b-902b-a509b2afa5ad
# ╟─2d1782c5-e544-4ee6-a872-c8a180bd0d72
# ╠═77e2c1e3-7756-4ab7-810a-03ccdc635aa1
# ╠═90856551-fff8-4d15-be66-6353091b5e50
# ╠═337055c4-55f4-4c02-ad69-473165ccfb64
# ╠═a52e5e94-6068-4545-962f-e02a485b62f5
# ╠═cc231e78-cfc0-4876-9f8b-980139f7d27f
# ╠═89b7d969-2294-4d07-a6a3-fcfa92498d48
# ╟─cdcc6bf6-f15e-45f1-b837-1db4e0fdc654
# ╠═7f7d8cb9-761c-4f30-a336-ab5657144961
# ╠═a2f72082-7145-42be-9f40-e00d18deb267
# ╟─bf8305f4-a5b8-4c79-8a01-e2aa18e4a5c5
# ╠═2751e203-8c03-44ec-a072-e3078d9a4f7b
# ╠═6dd92ee1-374d-47fa-ad61-b54764b23240
# ╟─99f274b4-91f3-41d0-b7d3-234badeb43d1
# ╠═4396bced-cae8-4107-ac83-48cc3c4146f2
# ╠═e63b488f-d1d7-488f-a1b8-6c32d2f66829
# ╠═87b5b241-db72-45ee-b3a7-a394f99510d9
# ╠═904576b4-2669-45d1-8078-310347f5fec0
# ╠═8561415c-7c5d-4c1a-8478-3b237c2cb360
# ╠═ff51f97d-2404-49c6-9339-4b201d6a94a9
# ╠═03037c94-7f26-479c-9772-fa2682e3ba37
# ╠═ad02a484-11b1-479f-b9e0-0e16e64a5109
# ╠═74845a4f-588b-44d3-a5e1-7155fd6a7b01
