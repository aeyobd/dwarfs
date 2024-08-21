### A Pluto.jl notebook ###
# v0.19.46

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
paramname = "exp2d_rs0.05_ana.toml"

# ╔═╡ 5287d506-2c6e-429c-90b9-9d1574784681
M_s_halo = 0.056

# ╔═╡ 01c0a5fc-f020-4668-a21e-cbfccb9a8826
r_s_halo = 2.76

# ╔═╡ 64578172-da38-49ed-8777-51b538aa9b18
prof_halo = lguys.NFW(M_s=M_s_halo/lguys.A_NFW(1), r_s = r_s_halo)

# ╔═╡ f4a09d4c-e345-4d2c-bf17-c802ea6d139c
lguys.calc_M(prof_halo, r_s_halo)

# ╔═╡ f1a7fa1f-bdcd-408c-ab27-b52916b1892f
overwrite = true

# ╔═╡ 0a2916d6-36f1-4fd7-b0db-d5f68e288283
md"""
# File Loading
"""

# ╔═╡ f1459790-10c5-4538-82d5-23e86dfe3f33
"""
loads in the parameterfile and will 
automatically populate the centres file and output files
"""
function load_params(paramname)
	params = TOML.parsefile(paramname); 

	if "centres_file" ∉ keys(params)
		params["centres_file"] = params["snapshot_dir"] * "/centres.csv"
	end

	if "output_file" ∉ keys(params)
		params["output_file"] = splitext(paramname)[1] * "_stars.hdf5"
	end
	
	if "mock_file" ∉ keys(params)
		params["mock_file"] = splitext(paramname)[1] * "_mock_stars.fits"
	end
	
	params
end

# ╔═╡ 9326c8a6-8b9b-4406-b00f-9febb3dcca46
cd(dir); params=load_params(paramname)

# ╔═╡ 43d45dea-af35-487d-bc37-18bd93b09a1f
"""
Loads in the instance of lguys.Profile from the parameter file.
"""
function load_profile(params)
	profile_class = getproperty(lguys, Symbol(params["profile"]))
	profile = profile_class(;lguys.dict_to_tuple(params["profile_kwargs"])...)

	return profile
end

# ╔═╡ d88cbe6a-b87b-45a2-8a3f-c5ef0b7a8935
profile = load_profile(params)

# ╔═╡ aa69bde5-ab93-4105-9d48-ad0ace88e3f0
r_h = lguys.get_r_h(profile)

# ╔═╡ 1b9f3101-4891-4a8c-8c73-41a222b6c26f
ρ_s(r) = lguys.calc_ρ(profile, r)

# ╔═╡ d77557e5-f7d8-40e9-ae40-a4b6b8df16cd
cen = CSV.read(params["centres_file"], DataFrame)[params["snapshot"] + 1, :]

# ╔═╡ 4cb09115-143d-456f-9c6a-19656f638677
begin 
	snapname = params["snapshot_dir"] * "/snapshot_" * lpad(params["snapshot"], 3, "0") * ".hdf5"
	println("loading snap ", snapname)
	snap_og = lguys.Snapshot(snapname)

	snap_og.x_cen = [cen.x, cen.y, cen.z ]
	snap_og.v_cen = [cen.vx, cen.vy, cen.vz]

	snap_og
end

# ╔═╡ 620a617b-c1d7-4718-a1b8-f4c3c18285d5
md"""
## Preprocessing snapshot
"""

# ╔═╡ 9af14bf8-bd18-48a9-af1e-62f9dab095b7
function calc_phi(snap)
	radii = lguys.calc_r(snap)
	Φs = lguys.calc_radial_discrete_Φ(radii, snap.masses)
	
	return Φs
end

# ╔═╡ 5851e36f-9376-4445-9bac-da3ce9d5ac7e
function calc_eps(snap, Φs)
	ke = lguys.calc_E_spec_kin(snap)
	ϵs = @. -Φs - ke

	return ϵs
end

# ╔═╡ 2c2b5c87-1e5e-4c32-ba64-e6d24323007c
function make_filter(ϵs, radii, params)
	filt = ϵs .> 0
	if "R_t" in keys(params["profile_kwargs"])
		filt .&= radii .< params["profile_kwargs"]["R_t"]
	end
	return filt
end

# ╔═╡ d4a0bd79-5d5e-48e8-bd8c-8865663423e1
function filter_and_sort(snap_og)
	snap = lguys.copy(snap_og)

	# calculate energies
	Φs = calc_phi(snap)

	ϵs = calc_eps(snap, Φs)
	radii = lguys.calc_r(snap)
	filt_snap = make_filter(ϵs, radii, params)

	snap = snap[filt_snap]
end

# ╔═╡ 0f5671e6-deb4-11ee-3178-4d3f920f23a2
begin
	# centre snapshot
	snap = filter_and_sort(snap_og)
	# calculate energies
	Φs = calc_phi(snap)
	ϵs = calc_eps(snap, Φs)
	radii = lguys.calc_r(snap)
end

# ╔═╡ cc562baf-e35b-4ecd-a045-ac5367014805
sum(radii .< r_s_halo) * snap.masses[1]

# ╔═╡ f7f746c9-cd03-4391-988b-dffeb31b2842
println(" $(sum(radii .< r_h)) stars within (3D) half-light radius")

# ╔═╡ 4fe0d79b-92f8-4766-808a-9bfbf4c5ab7f
"""
Given radii, masses, and radii to evaluate at, returns lerp'ed interior mass
"""
function calc_M_in(radii, masses, r)
	Min = cumsum(masses) ./ sum(masses)

	Min_lerp = lguys.lerp(radii, Min)
	M = Min_lerp.(r)
	
	return M
end

# ╔═╡ e935a80a-7f4f-4143-ae6c-bee62b82c30e
md"""
# Core calculation
"""

# ╔═╡ 23158c79-d35c-411a-a35a-950f04214e19
begin 
	epsilon = 1e-10
	r_max = 1e3
	M_s_tot = 4π * quadgk(r-> r^2 * ρ_s(r), epsilon, r_max)[1]
	M_s(r) = 4π * quadgk(r-> r^2 * ρ_s(r) / M_s_tot, epsilon, r)[1]
end

# ╔═╡ 29619cc3-1be3-4b24-92e0-ceccfd4a3f59
"""
given the radii and parameters for stellar profile, returns the
radial bins used in the calculation. 

"""
function make_radius_bins(radii::AbstractVector, params::Dict)
	log_radii = log10.(radii)

	
	lr_min = minimum(log_radii)
	lr_max = maximum(log_radii)
	Nr = params["num_radial_bins"]
	
	r_e = 10 .^ LinRange(lr_min, lr_max, Nr+1)

	return r_e
	
end

# ╔═╡ deb46b0b-3504-4231-9f85-99f27caf924b
r_e = make_radius_bins(radii, params)

# ╔═╡ 36b4adbd-d706-4e72-a922-53080c67946c
r = lguys.midpoint(r_e)

# ╔═╡ f79414b4-620e-440e-a421-7bc13d373546
M = calc_M_in(radii, snap.masses, r)

# ╔═╡ f3e95cfc-0087-4afa-8e88-a4e1628c50a0
"""
Given the stellar mass function M_s, returns total missing stars
"""
function print_missing(r_bins, M_s)
	N_s_out = 1 - M_s(r_e[end])
	N_s_in = M_s(r_e[1])
	println("missing $N_s_out stars outside")
	println("missing $N_s_in stars inside")
end

# ╔═╡ 41a08be6-d822-4a19-afe5-62c7dc9ff118
print_missing(r_e, M_s)

# ╔═╡ dfa675d6-aa32-45c3-a16c-626e16f36083
ψ = -lguys.calc_Φ.(prof_halo, r)

# ╔═╡ ab4d4458-fd1d-417a-8a74-5e06f41af166
ν_dm = lguys.calc_ρ.(prof_halo, r)

# ╔═╡ 20f858d6-a9f5-4880-a431-60b395cc7e50
ν_s = max.(ρ_s.(r) ./ M_s_tot, 0)

# ╔═╡ 1fab14ec-6bfc-4243-bc48-915d1a129925
begin 
	f_dm = lguys.calc_fϵ(ν_dm, ψ, r)
	f_s = lguys.calc_fϵ(ν_s, ψ, r)
end

# ╔═╡ 126c6825-723f-4d13-b5a3-64cba72fc867
md"""
The density distribution of the stars (analytic) and dark matter (calculated) from the snapsho
"""

# ╔═╡ 7481f47a-2b8a-45a3-9b4e-31ea14e9d331
md"""
The potential $\psi = -\Phi$ as a function of log radii (for the spherically calculated & interpolated and actual snapshot from Gadget 4)
"""

# ╔═╡ 1c0899c6-2692-45ab-b9e6-668dc576d679
function make_energy_bins(ψ, params)
	E_max = ψ[1]
	E_min = ψ[end]
	E = LinRange(E_min, E_max, params["num_energy_bins"] + 1)
end

# ╔═╡ 3025546e-3ce1-4a78-824d-24f644238e32
E = make_energy_bins(ψ, params)

# ╔═╡ 500c67b2-8c4a-4d03-b4e4-39beff43a46c
f_dm_e = f_dm.(E)

# ╔═╡ 9e492a55-7b20-4eca-aead-c7afeee63f11
f_s_e = f_s.(E)

# ╔═╡ 8b66d00d-529b-4e8c-9386-b17117996579
md"""
The calculated distribution function as a function of log specific binding energy $\epsilon =  -\Phi - T$
"""

# ╔═╡ 7409a024-cea4-47a5-84d2-846c96d88b7a
begin 
	probs = f_s_e ./ f_dm_e
	probs ./= sum(probs .* lguys.gradient(E)) # pdf, dN/dE
	prob = lguys.lerp(E, probs)
end

# ╔═╡ 3b229c8e-9320-4c07-b948-c34a0c082341
begin 
	ps = prob.(ϵs)
	print(sum(ps .< 0), " negative probabilities")
	ps[ps .< 0] .= 0
	ps[isnan.(ps)] .= 0
	ps ./= sum(ps)
end

# ╔═╡ daf54cb4-06a3-4a8e-a533-354ca8740aec
md"""
A histogram of the assigned (positive) stellar weights. See the number of negative probabilities above
"""

# ╔═╡ b67a17dc-7a59-4b62-9120-4c2ada12298c
md"""
The main result. The reconstructed density profile
"""

# ╔═╡ 06e1b872-ce52-434f-a8d1-3b0a5055eed2
md"""
A histogram of the stellar density
"""

# ╔═╡ ee866164-c6f2-4f70-bde1-360abd5fd80e
md"""
Histogram of same region in dark matter only (over a smaller dynamic range). Dark matter is much more extended (as expected)
"""

# ╔═╡ a3071af2-beff-408d-b109-d4f289f8f7f4
md"""
## Last checks and saving profile
"""

# ╔═╡ 15c911ae-679e-4f0f-aad5-9f3fe51bdeab
p_idx = snap.index

# ╔═╡ 54d6ca0d-a1f8-44b3-895d-a45ad72ff42c
idx_excluded = setdiff(snap_og.index, snap.index)

# ╔═╡ 7bc02274-aa44-4609-b9a6-e409de5172af
begin
	
	idx_all = vcat(p_idx, idx_excluded)
	ps_all = vcat(ps, zeros(length(idx_excluded)))

	_sort = sortperm(idx_all)
	idx_all = idx_all[_sort]
	ps_all = ps_all[_sort]
end

# ╔═╡ f66a4a6d-b31c-481b-bf7a-4801c783ceb4
idx_all == sort(snap_og.index) # check we didn't lose any star particles

# ╔═╡ c18bc8b9-7ea1-47e2-8d7e-971a0943917a
maximum(idx_all) == length(idx_all)

# ╔═╡ 92798e2e-2356-43c2-a4a9-82a70619a5f5
0 == sum(ps_all[idx_excluded]) # should be zero

# ╔═╡ 4671b864-470d-4181-993b-4e64d5687460
sum(ps_all) # should be 1

# ╔═╡ bd8489da-ac53-46c6-979e-06d5dc6e25d1
"""
Writes the stars to an hdf5 file
"""
function write_stars()
	outname = params["output_file"]
	if isfile(outname)
		if overwrite
			rm(outname)
		else
			println("file already exists")
			return
		end
	end


	h5write(outname, "/index", idx_all)
	h5write(outname, "/probabilities", ps_all)
	println("probabilities written to $outname")
end

# ╔═╡ f42cb7e1-64b7-47da-be05-dd50c2471fb3
write_stars()

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
	stephist!(log10.(radii), bins=log10.(r_e), normalization=:pdf, label="dark matter")
	stephist!(log10.(radii), bins=log10.(r_e), weights=ps, normalization=:pdf, label="stars (nbody)")
	axislegend()
	fig
end

# ╔═╡ 84fdc265-988c-40db-87e5-44ba55d0e412
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"$\epsilon$ (binding energy)", 
		ylabel="asinh value", 
		)


	scale = 1e-10

	y_trans(x) =  asinh.(x / scale)
	
	h1 = Arya.histogram(ϵs, normalization=:pdf)
	h1.values .= y_trans(h1.values)

	h_s = Arya.histogram(ϵs, weights=ps, normalization=:pdf)
	h_s.values .= y_trans(h_s.values)

	probs = prob.(E)
	probs = y_trans(probs)
	
	lines!((E), probs, color=Arya.COLORS[3], label="f_s / f_dm")
	
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
	lines!(log10.(r), log10.(ν_dm), label="DM (analytic)")

	y_obs = lguys.calc_ρ_hist(radii, r_e, weights=snap.masses)[2]
	scatter!(log10.(r), log10.(y_obs), label="analy")

	scatter!(log10.(r), log10.(ν_s), label="stars (analytic)")


	
	axislegend()
	fig
end

# ╔═╡ 75d23b44-71e7-4e28-ad3e-c537f3d4422f
let
	fig = Figure()
	ax = Axis(fig[1,1],xlabel="log ϵ", ylabel="log f", limits=(nothing, (-15, 7)) )
	lines!(log10.(E), nm.log10.(f_s_e), label="stars")
	lines!(log10.(E), nm.log10.(f_dm_e), label="DM")

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
	lines!(log10.(r), ψ, label="theoretical")

	skip = 100
	
	scatter!(log10.(radii[1:skip:end]), -snap.Φs[1:skip:end], label="snapshot")
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

# ╔═╡ 6fba7fa7-9a50-4379-b376-5c07f3638411
begin 
	r_nbody, ν_s_nbody = lguys.calc_ρ_hist(radii, 100, weights=ps)
	r_nbody = midpoints(r_nbody)
end

# ╔═╡ a9335e17-a410-455a-9a9e-d63706a026bd
let
	fig = Figure(size=(700, 500))
	ax = Axis(fig[1,1], ylabel=L"\log \nu", 
		limits=((0, 1), (-15, 3))
		)
	lines!(log10.(r), nm.log10.(ν_s), label="stars")
	scatter!(log10.(r_nbody) , nm.log10.(ν_s_nbody), label="nbody", color=COLORS[2])

	ax2 = Axis(fig[2,1], 
		xlabel=L"\log r / \textrm{kpc}", ylabel=L"\Delta\log \nu ", 
		limits=((log10(0.5r_h), log10(100r_h)), (-1, 1)))
	
	scatter!(log10.(r_nbody), nm.log10.(ν_s_nbody) .- nm.log10.(ρ_s.(r_nbody)), label="",
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

	hist!(log10.(radii[ps .> 2e-6]))

	fig
end

# ╔═╡ 77e2c1e3-7756-4ab7-810a-03ccdc635aa1
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="stellar weights", ylabel="frequency", yscale=log10)
	hist!(ps, bins=100)
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

	scatter!(cen.x, cen.y)

	Colorbar(fig[1, 2],h )

	fig
end

# ╔═╡ cc231e78-cfc0-4876-9f8b-980139f7d27f
let
	fig = Figure()
	filt = ps .> 0
	
	ax = Axis(fig[1,1], 
		limits=(-r_hist, r_hist, -r_hist, r_hist), 
		aspect=1,
		title="stars",
		xlabel="x / kpc",
		ylabel="y / kpc"
	)
	
	h = Arya.hist2d!(ax, 
		lguys.get_x(snap)[filt], lguys.get_y(snap)[filt], 
		weights=ps[filt], bins=N_hist,
		colorscale=log10,
		colorrange=(1e-7, nothing)
	)	

	scatter!(cen.x, cen.y)

	Colorbar(fig[1, 2], h)
	fig
end

# ╔═╡ 7f7d8cb9-761c-4f30-a336-ab5657144961
let
	r = lguys.calc_r(snap_og, [cen.x, cen.y, cen.z])
	ms = ps_all[snap_og.index]
	
	r_h2 = calc_r_h(r, ms)
	println(r_h2)
	
	r_e = 10 .^ LinRange(log10.(minimum(r)), log10(maximum(r)), 50)
	r, ν_s_nbody = lguys.calc_ρ_hist(r, r_e, weights=ms)
 
	r = lguys.midpoint(r)
	ν_s = ρ_s.(r)
	
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
	x = snap_og.positions[1, :] .- cen.x
	y = snap_og.positions[2, :] .- cen.y
	R = @. sqrt(x^2 + y^2)

	prof = lguys.calc_properties(R, weights=ps_all[snap_og.index], bins=100)

	
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
distance = 86

# ╔═╡ 062cea41-3b51-474e-b877-7a4d96813fbc
R_s_arcmin = lguys.kpc_to_arcmin(profile.R_s, distance)

# ╔═╡ 7b7832d1-6739-453e-aaee-fa16a6000d26
function mock_obs(snap_og, ps_all; distance=86)
	x_sun = [8.122, 0, 0]
	shift_vec = x_sun .+ lguys.rand_unit() * distance
	
	ps = ps_all[snap_og.index]
	p_min = 1e-25
	filt = p_min * maximum(ps) .< ps

	ps = ps[filt]
	snap_shift = copy(snap_og[filt])
	snap_shift.positions .+= shift_vec

	obs = lguys.to_sky(snap_shift)

	println("number of final stars: ", length(snap_shift))
		
	obs_df = DataFrame(; 
	collect(key => [getproperty(o, key) for o in obs]
			for key in [:ra, :dec, :pm_ra, :pm_dec, :distance, :radial_velocity])...
		
	)

	obs_df[!, "index"] = snap_shift.index
	obs_df[!, "probability"] = ps
	
	rename!(obs_df, "pm_ra"=>"pmra")
	rename!(obs_df, "pm_dec"=>"pmdec")

	return obs_df
end

# ╔═╡ 87b5b241-db72-45ee-b3a7-a394f99510d9
obs_mock = mock_obs(snap_og, ps_all)

# ╔═╡ ff51f97d-2404-49c6-9339-4b201d6a94a9
let
	ra = obs_mock.ra
	dec = obs_mock.dec
	ms = obs_mock.probability

	ra0, dec0 = lguys.calc_centre2D(ra, dec, "mean", ms)
	xi, eta = lguys.to_tangent(ra, dec, ra0, dec0)
	R = @. 60sqrt(xi^2 + eta^2)
	
	prof = lguys.calc_properties(R, weights=ms, bins=50)


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

# ╔═╡ 895dfb3a-c8c4-493f-997e-bc3a41d99686
function make_sample(snap; rel_p_cut=1e-15, r_max=Inf,
	Frame=lguys.HelioRest
)

	ps = probabilities[snap.index]
	
	p_min = rel_p_cut * maximum(ps)
	println("adopting p_min = $p_min")

	filt = ps .> p_min
	snap_stars = snap[filt]
	ps = ps[filt]
	println("sanity check: ", ps == probabilities[snap_stars.index])
	println("number of final stars: ", length(snap_stars))
	
	obs_pred = lguys.to_sky(snap_stars, SkyFrame=Frame)
	
	obs_df = DataFrame(; 
	collect(key => [getproperty(o, key) for o in obs_pred]
			for key in [:ra, :dec, :pm_ra, :pm_dec, :distance, :radial_velocity])...
		
	)

	obs_df[!, "index"] = snap_stars.index
	obs_df[!, "probability"] = ps
	
	rename!(obs_df, "pm_ra"=>"pmra")
	rename!(obs_df, "pm_dec"=>"pmdec")

	return obs_df
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
# ╠═5287d506-2c6e-429c-90b9-9d1574784681
# ╠═01c0a5fc-f020-4668-a21e-cbfccb9a8826
# ╠═64578172-da38-49ed-8777-51b538aa9b18
# ╠═f4a09d4c-e345-4d2c-bf17-c802ea6d139c
# ╠═cc562baf-e35b-4ecd-a045-ac5367014805
# ╠═f1a7fa1f-bdcd-408c-ab27-b52916b1892f
# ╟─0a2916d6-36f1-4fd7-b0db-d5f68e288283
# ╠═f1459790-10c5-4538-82d5-23e86dfe3f33
# ╠═9326c8a6-8b9b-4406-b00f-9febb3dcca46
# ╟─43d45dea-af35-487d-bc37-18bd93b09a1f
# ╠═d88cbe6a-b87b-45a2-8a3f-c5ef0b7a8935
# ╠═aa69bde5-ab93-4105-9d48-ad0ace88e3f0
# ╠═f7f746c9-cd03-4391-988b-dffeb31b2842
# ╠═1b9f3101-4891-4a8c-8c73-41a222b6c26f
# ╠═d77557e5-f7d8-40e9-ae40-a4b6b8df16cd
# ╠═4cb09115-143d-456f-9c6a-19656f638677
# ╟─620a617b-c1d7-4718-a1b8-f4c3c18285d5
# ╠═9af14bf8-bd18-48a9-af1e-62f9dab095b7
# ╠═5851e36f-9376-4445-9bac-da3ce9d5ac7e
# ╠═2c2b5c87-1e5e-4c32-ba64-e6d24323007c
# ╠═d4a0bd79-5d5e-48e8-bd8c-8865663423e1
# ╠═0f5671e6-deb4-11ee-3178-4d3f920f23a2
# ╟─4fe0d79b-92f8-4766-808a-9bfbf4c5ab7f
# ╟─e935a80a-7f4f-4143-ae6c-bee62b82c30e
# ╠═f79414b4-620e-440e-a421-7bc13d373546
# ╠═23158c79-d35c-411a-a35a-950f04214e19
# ╟─29619cc3-1be3-4b24-92e0-ceccfd4a3f59
# ╠═deb46b0b-3504-4231-9f85-99f27caf924b
# ╠═36b4adbd-d706-4e72-a922-53080c67946c
# ╟─f3e95cfc-0087-4afa-8e88-a4e1628c50a0
# ╠═41a08be6-d822-4a19-afe5-62c7dc9ff118
# ╠═dfa675d6-aa32-45c3-a16c-626e16f36083
# ╠═ab4d4458-fd1d-417a-8a74-5e06f41af166
# ╠═20f858d6-a9f5-4880-a431-60b395cc7e50
# ╠═1fab14ec-6bfc-4243-bc48-915d1a129925
# ╟─126c6825-723f-4d13-b5a3-64cba72fc867
# ╟─7481f47a-2b8a-45a3-9b4e-31ea14e9d331
# ╠═1c0899c6-2692-45ab-b9e6-668dc576d679
# ╠═3025546e-3ce1-4a78-824d-24f644238e32
# ╠═500c67b2-8c4a-4d03-b4e4-39beff43a46c
# ╠═9e492a55-7b20-4eca-aead-c7afeee63f11
# ╟─8b66d00d-529b-4e8c-9386-b17117996579
# ╠═7409a024-cea4-47a5-84d2-846c96d88b7a
# ╠═3b229c8e-9320-4c07-b948-c34a0c082341
# ╟─daf54cb4-06a3-4a8e-a533-354ca8740aec
# ╟─b67a17dc-7a59-4b62-9120-4c2ada12298c
# ╟─06e1b872-ce52-434f-a8d1-3b0a5055eed2
# ╟─ee866164-c6f2-4f70-bde1-360abd5fd80e
# ╟─a3071af2-beff-408d-b109-d4f289f8f7f4
# ╠═15c911ae-679e-4f0f-aad5-9f3fe51bdeab
# ╠═54d6ca0d-a1f8-44b3-895d-a45ad72ff42c
# ╠═7bc02274-aa44-4609-b9a6-e409de5172af
# ╠═f66a4a6d-b31c-481b-bf7a-4801c783ceb4
# ╠═c18bc8b9-7ea1-47e2-8d7e-971a0943917a
# ╠═92798e2e-2356-43c2-a4a9-82a70619a5f5
# ╠═4671b864-470d-4181-993b-4e64d5687460
# ╠═bd8489da-ac53-46c6-979e-06d5dc6e25d1
# ╠═f42cb7e1-64b7-47da-be05-dd50c2471fb3
# ╟─4d1991ea-9496-48c7-a400-8fefbecefcd2
# ╟─5b30475b-b4c4-4c87-817d-0d885546d004
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
# ╠═a9335e17-a410-455a-9a9e-d63706a026bd
# ╠═33a26663-0c08-411b-902b-a509b2afa5ad
# ╠═77e2c1e3-7756-4ab7-810a-03ccdc635aa1
# ╠═90856551-fff8-4d15-be66-6353091b5e50
# ╠═a52e5e94-6068-4545-962f-e02a485b62f5
# ╠═cc231e78-cfc0-4876-9f8b-980139f7d27f
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
# ╠═895dfb3a-c8c4-493f-997e-bc3a41d99686
# ╠═b89fbbb0-b161-4fa7-82ba-c06e3d32f8b8
# ╠═d29baea8-0e5c-43e0-9aa2-987d6cc08e88
