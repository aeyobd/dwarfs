### A Pluto.jl notebook ###
# v0.19.42

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
	using Revise
end

# ╔═╡ 17ffde4b-5796-4915-9741-d594cf0c5ca7
md"""
# Daniel's script for painting stars
based on Rapha's codes. 
Eddington inversion

Given the directory, the script changes to there and works there
"""

# ╔═╡ 93045024-a91d-4b31-9a5a-7c999afdb9ec
md"""
# Inputs
"""

# ╔═╡ 530c6c09-4454-4952-8351-dccbb4ed429f
import TOML

# ╔═╡ 48ce69f2-09d5-4166-9890-1ab768f3b59f
dir = "/cosma/home/durham/dc-boye1/data/dwarfs/models/sculptor/isolation/1e6/stars"

# ╔═╡ 7809e324-ba5f-4520-b6e4-c7727c227154
paramname = "exp2d_rs0.16.toml"

# ╔═╡ 9326c8a6-8b9b-4406-b00f-9febb3dcca46
begin 
	cd(dir)
	params = TOML.parsefile(paramname); 

	if "centres_file" ∉ keys(params)
		params["centres_file"] = params["snapshot_dir"] * "/centres.csv"
	end

	if "output_file" ∉ keys(params)
		params["output_file"] = splitext(paramname)[1] * "_stars.hdf5"
	end
	
	params
end

# ╔═╡ 16e0729e-9a75-458e-a56c-73967c819c31
profile_class = getproperty(lguys, Symbol(params["profile"]))

# ╔═╡ 46be1f99-6f64-476e-84cb-11d4b6504a86
NamedTuple(d::Dict) = (; zip(Symbol.(keys(d)), values(d))...)

# ╔═╡ d88cbe6a-b87b-45a2-8a3f-c5ef0b7a8935
profile = profile_class(;NamedTuple(params["profile_kwargs"])...)

# ╔═╡ aa69bde5-ab93-4105-9d48-ad0ace88e3f0
r_h = profile.R_s #lguys.get_r_h(profile)

# ╔═╡ 855fd729-22b6-4845-9d2b-e796d4a15811
begin 
	# parameters 
	ρ_s(r) = lguys.calc_ρ(profile, r)
	Nr = params["num_radial_bins"]
	
	NE = params["num_energy_bins"]


	overwrite = true
end

# ╔═╡ 4cb09115-143d-456f-9c6a-19656f638677
begin 
	snapname = params["snapshot_dir"] * "/snapshot_" * lpad(params["snapshot"], 3, "0") * ".hdf5"
	println(snapname)
	snap = lguys.Snapshot(snapname)
	
end

# ╔═╡ d77557e5-f7d8-40e9-ae40-a4b6b8df16cd
cen = CSV.read(params["centres_file"], DataFrame)[params["snapshot"] + 1, :]

# ╔═╡ cf50c4b4-a634-4310-8aa5-91ab653313a9
cen

# ╔═╡ d43a1be5-7ecd-474b-9580-1191b669ff24
"R_t" in keys(params["profile_kwargs"])

# ╔═╡ 0f5671e6-deb4-11ee-3178-4d3f920f23a2
begin
	# centre snapshot
	snap_i = lguys.copy(snap)
	snap_i.positions .-=  [cen.x, cen.y, cen.z]
	snap_i.velocities .-= [cen.vx, cen.vy, cen.vz]


	# sort by radii
	radii = lguys.calc_r(snap_i.positions)

	snap_i = snap_i[sortperm(radii)]
	radii = sort(radii)
	print(radii[1:10])
	vels = lguys.calc_r(snap_i.velocities)

	# calculate energies
	Φs = lguys.calc_radial_discrete_Φ(radii, snap_i.masses)
	ke = @. 1/2 * vels^2
	ϵs = @. -Φs - ke

	# filter snapshot
	idx = ϵs .> 0
	if "R_t" in keys(params["profile_kwargs"])
		#idx .&= radii .< profile.R_t / 1.3
	end
	idx_excluded = snap_i.index[map(!, idx)]
	
	snap_i = snap_i[idx]
	radii = radii[idx]
	vels = vels[idx]
	
	ke = ke[idx]
	ϵs = ϵs[idx]
	Φs = Φs[idx]
end

# ╔═╡ f7f746c9-cd03-4391-988b-dffeb31b2842
sum(radii .< r_h)

# ╔═╡ eb5ffa80-6959-4a2c-980f-f818f6a03c14
radii

# ╔═╡ 8af6b094-ea45-4f01-abb4-9b4686e81719
snap.positions

# ╔═╡ f79414b4-620e-440e-a421-7bc13d373546
Mins = cumsum(snap_i.masses) ./ sum(snap_i.masses)

# ╔═╡ 45acc05f-85a8-4bbc-bd43-a34583c983b3
Φs

# ╔═╡ 23158c79-d35c-411a-a35a-950f04214e19
begin 
	M_s_tot = 4π * quadgk(r-> r^2 * ρ_s(r), 1e-5, 1.9)[1]
	M_s(r) = 4π * quadgk(r-> r^2 * ρ_s(r) / M_s_tot, 1e-5, r)[1]
end

# ╔═╡ 771413c9-3edf-4687-b4cd-3a60ae0ceda2
M_s_tot

# ╔═╡ b3663c19-c029-4a1e-ab82-a05177e3a5d0
import StatsBase: percentile

# ╔═╡ 36b4adbd-d706-4e72-a922-53080c67946c
begin
	log_radii = log10.(radii)
	r_min = radii[1]
	r_max = 1.8 #radii[end]
	if params["bin_method"] == "equal_width"
		r_e = 10 .^ LinRange(log10.(r_min), log10.(r_max), Nr+1)
	elseif params["bin_method"] == "equal_number"
		r_e = percentile(radii, LinRange(0, 100, Nr+1))
	else
		error("bin method unknown")
	end
	r = lguys.midpoint(r_e)
end

# ╔═╡ e37bf6d7-9445-49bf-8333-f68ad25436b2
r_e

# ╔═╡ 7a39cd4f-9646-4969-9410-b093bca633cb
begin 
	# just a informative note
	N_s_out = 1 - M_s(r_max)
	N_s_in = M_s(r_min)
	println("missing $N_s_out stars outside")
	println("missing $N_s_in stars inside")
end

# ╔═╡ bd1bca1d-0982-47d8-823e-eadc05191b88
begin 
	m_dm = Vector{Float64}(undef, Nr)
	for i in 1:Nr
		filt = r_e[i] .<= radii .< r_e[i+1]
		m_dm[i] = sum(snap_i.masses[filt])
	end
end

# ╔═╡ dfa675d6-aa32-45c3-a16c-626e16f36083
begin 
	ψ = lguys.lerp(radii, -Φs).(r)
	
	ν_dm = lguys.calc_ρ_hist(radii, r_e)[2]
	ν_dm ./= length(snap)

end

# ╔═╡ 830faf5d-da17-4939-ad7c-080503b66990
begin
	M = lguys.lerp(radii, Mins).(r)
	Ms = M_s.(r)
end

# ╔═╡ de614f3d-b3e4-41a0-886b-594b0c93cc5b
r

# ╔═╡ bbd02922-37fc-459d-8702-451127d2b9a2
sum(radii .< 0.1)

# ╔═╡ 20f858d6-a9f5-4880-a431-60b395cc7e50
ν_s = max.(ρ_s.(r) ./ M_s_tot, 0)

# ╔═╡ 1fab14ec-6bfc-4243-bc48-915d1a129925
begin 
	ψ_i = lguys.lerp(radii, -Φs).(r)
	f_dm = lguys.calc_fϵ(ν_dm, ψ_i, r)
	f_s = lguys.calc_fϵ(ν_s, ψ_i, r)
end

# ╔═╡ 126c6825-723f-4d13-b5a3-64cba72fc867
md"""
The density distribution of the stars (analytic) and dark matter (calculated) from the snapsho
"""

# ╔═╡ 8bb8736d-a41b-4dac-a6cd-06d0d4704654
let
	fig = Figure()
	ax = Axis(fig[1,1], limits=(-1.2, 3, -15, 2),
		xlabel="log r", ylabel="log density")
	lines!(log10.(r), log10.(ν_dm), label="DM")
	scatter!(log10.(r), log10.(ν_s), label="stars (analytic)")

	axislegend()
	fig
end

# ╔═╡ c537fdc7-4778-40e6-ab82-7e3151394473
ν_s

# ╔═╡ 46cf718e-79d5-4017-b113-72efa692a014
r

# ╔═╡ b53887a2-0b14-4d95-bdf5-d70960cfe736
pn = lguys.Exp2D(M=1, R_s=1)

# ╔═╡ c2e070b8-6f15-48e1-8ff3-8145f67b2001
lguys.calc_M(pn, 100)

# ╔═╡ 837dce1c-fadb-42df-92f9-49a1bd859d1f
lguys.calc_M_2D(pn, 5)

# ╔═╡ e9129b04-ffca-4b10-8cda-4955841a7655
lguys.calc_Σ_from_ρ(pn, 1)

# ╔═╡ fdf4f368-962a-4c3a-9571-8489740c31b0
lguys.calc_Σ(pn, 1)

# ╔═╡ c9098d39-3573-450b-8c5e-a6d847c85a3e
lguys.calc_ρ_from_Σ(pn, 1)

# ╔═╡ 162e0102-ed49-4648-acf9-c78e52d41ed3
lguys.calc_ρ(pn, 1)

# ╔═╡ 7481f47a-2b8a-45a3-9b4e-31ea14e9d331
md"""
The potential $\psi = -\Phi$ as a function of log radii (for the spherically calculated & interpolated and actual snapshot from Gadget 4)
"""

# ╔═╡ b625d8a5-7265-4849-9bd6-ca8064d392eb
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel=L"\log r / \textrm{kpc}", ylabel=L"\Psi(r)")
	lines!(log10.(r), ψ, label="interpolated")
	lines!(log10.(radii), -snap_i.Φs, label="snapshot")

	fig
end

# ╔═╡ 78ce5a98-fd3f-4e39-981f-2bea58b117bf
begin 
	E_max = ψ[1]
	E_min = ψ[end]
	E = LinRange(E_min, E_max, NE + 1)
	f_dm_e = f_dm.(E)
	f_s_e = f_s.(E)
end

# ╔═╡ 8b66d00d-529b-4e8c-9386-b17117996579
md"""
The calculated distribution function as a function of log specific binding energy $\epsilon =  -\Phi - T$
"""

# ╔═╡ 75d23b44-71e7-4e28-ad3e-c537f3d4422f
let
	fig = Figure()
	ax = Axis(fig[1,1],xlabel="log ϵ", ylabel="log f", limits=(nothing, (-15, 7)) )
	lines!(log10.(E), nm.log10.(f_s_e), label="stars")
	lines!(log10.(E), nm.log10.(f_dm_e), label="DM")

	axislegend(ax, position=:lt)
	fig
end

# ╔═╡ 49c6372d-4ac4-46e2-aace-4cae0564f9b5
f_dm_e

# ╔═╡ 7409a024-cea4-47a5-84d2-846c96d88b7a
begin 
	probs = f_s_e ./ f_dm_e
	probs ./= sum(probs .* lguys.gradient(E)) # pdf, dN/dE
	prob = lguys.lerp(E, probs)
end

# ╔═╡ a389fde4-f05f-48ff-9a8c-8eb0643a849d
ϵs

# ╔═╡ 3b229c8e-9320-4c07-b948-c34a0c082341
begin 
	ps = prob.(ϵs)
	print(sum(ps .< 0), " negative probabilities")
	ps[ps .< 0] .= 0
	ps[isnan.(ps)] .= 0
	ps ./= sum(ps)
end

# ╔═╡ 84fdc265-988c-40db-87e5-44ba55d0e412
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=L"$\epsilon$ (binding energy)", 
		ylabel="relative density", 
		)
	
	lines!((E), 20*prob.(E) ./maximum(prob.(E)), color=Arya.COLORS[3], label="f_s / f_dm")
	stephist!((ϵs), normalization=:pdf, label="dark matter")
	stephist!((ϵs), weights=100ps, label="stars (nbody)")
	vlines!([maximum(ϵs)], color="grey", linestyle=:dot, label=L"\epsilon_\textrm{max}")

	axislegend(ax, position=:lt)

	fig
end

# ╔═╡ daf54cb4-06a3-4a8e-a533-354ca8740aec
md"""
A histogram of the assigned (positive) stellar weights. See the number of negative probabilities above
"""

# ╔═╡ 77e2c1e3-7756-4ab7-810a-03ccdc635aa1
let 
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="stellar weights", ylabel="frequency", yscale=log10)
	hist!(ps, bins=100)
	fig
end

# ╔═╡ 5b9c7242-5d4b-4f5a-83e5-3e4d11017aa5
sum(ps .> 0)

# ╔═╡ b1a34b9b-c6ca-4818-95b2-5b55fb32511e
sum(ps .== 0)

# ╔═╡ a5bc5ce3-8e33-4514-bc2d-4b4299f104f9
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="log radii", ylabel="pdf (dn / dlog r)")
	stephist!(log10.(radii), normalization=:pdf, label="dark matter")
	stephist!(log10.(radii), weights=ps, normalization=:pdf, label="stars")
	axislegend()
	fig
end

# ╔═╡ 33a26663-0c08-411b-902b-a509b2afa5ad
let
	fig = Figure()
	Axis(fig[1,1], xlabel="log radii", ylabel="pstar > 0.025")

	hist!(log10.(radii[ps .> 2e-6]))

	fig
end

# ╔═╡ a73db457-2cc2-4917-bb25-0429f4daecbd
length(radii)

# ╔═╡ 0e7a922f-4e61-4d02-ac28-3915f5d1c9da
length(ps)

# ╔═╡ 7d69e393-c4db-4ff0-ab5a-85ac50c785c2
maximum(ps)

# ╔═╡ 6fba7fa7-9a50-4379-b376-5c07f3638411
begin 
	m_s = Vector{Float64}(undef, Nr)

	for i in 1:Nr
		filt = r_e[i] .<= radii .< r_e[i+1]
		m_s[i] = sum(ps[filt])
	end
	m_s ./= sum(m_s)
	ν_s_nbody = lguys.calc_ρ_hist(radii, r_e, weights=ps)[2]
end

# ╔═╡ b67a17dc-7a59-4b62-9120-4c2ada12298c
md"""
The main result. The reconstructed density profile
"""

# ╔═╡ c01a4635-dd57-4793-8893-1bf5fb6996b3
nm.log10.(ν_s)

# ╔═╡ 131bd303-4f0f-4b24-bd90-5c65cf342f4c
radii

# ╔═╡ 52bf618a-5494-4cfb-9fd0-ee82f9682116
radii[findfirst(x->x>0.5, cumsum(ps))]

# ╔═╡ ec46946a-0bf7-4374-af44-8d8b2ab6a3df
ps

# ╔═╡ 06e1b872-ce52-434f-a8d1-3b0a5055eed2
md"""
A histogram of the stellar density
"""

# ╔═╡ 90856551-fff8-4d15-be66-6353091b5e50
begin 
	r_hist = 2
	N_hist = 100
end

# ╔═╡ 9064b0ca-4d26-4cb8-bbc8-353da44ffc26
profile

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
		lguys.get_x(snap_i)[filt], lguys.get_y(snap_i)[filt], 
		weights=ps[filt], bins=N_hist,
		colorscale=log10,
		colorrange=(1e-7, nothing)
	)	

	Colorbar(fig[1, 2], h)
	fig
end

# ╔═╡ ee866164-c6f2-4f70-bde1-360abd5fd80e
md"""
Histogram of same region in dark matter only (over a smaller dynamic range). Dark matter is much more extended (as expected)
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
		lguys.get_x(snap_i), lguys.get_y(snap_i), 
		bins=N_hist,
		colorscale=log10,
		colorrange=(1e-0, nothing)
	)	

	Colorbar(fig[1, 2],h )

	fig
end

# ╔═╡ a3071af2-beff-408d-b109-d4f289f8f7f4
md"""
## Last checks and saving profile
"""

# ╔═╡ a483545a-42ea-474b-a2ca-6f1bd4c7275b
p_idx = snap_i.index

# ╔═╡ 7bc02274-aa44-4609-b9a6-e409de5172af
begin
	idx_all = vcat(p_idx, idx_excluded)
	ps_all = vcat(ps, zeros(length(idx_excluded)))

	_sort = sortperm(idx_all)
	idx_all = idx_all[_sort]
	ps_all = ps_all[_sort]
end

# ╔═╡ f66a4a6d-b31c-481b-bf7a-4801c783ceb4
idx_all == sort(snap.index) # check we didn't lose any star particles

# ╔═╡ c18bc8b9-7ea1-47e2-8d7e-971a0943917a
maximum(idx_all) == length(idx_all)

# ╔═╡ 92798e2e-2356-43c2-a4a9-82a70619a5f5
sum(ps_all[idx_excluded]) # should be zero

# ╔═╡ 4671b864-470d-4181-993b-4e64d5687460
sum(ps_all) # should be 1

# ╔═╡ bd8489da-ac53-46c6-979e-06d5dc6e25d1
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

# ╔═╡ 793f701e-376d-4ab3-ae74-6abe29e3c3ae
lguys.get_r_h(profile)

# ╔═╡ a9335e17-a410-455a-9a9e-d63706a026bd
let
	fig = Figure(size=(700, 500))
	ax = Axis(fig[1,1], ylabel=L"\log \nu", 
		limits=((0, 1), (-15, 3))
		)
	lines!(log10.(r), nm.log10.(ν_s), label="stars")
	scatter!(log10.(r) , nm.log10.(ν_s_nbody), label="nbody", color=COLORS[2])

	ax2 = Axis(fig[2,1], 
		xlabel=L"\log r / \textrm{kpc}", ylabel=L"\Delta\log \nu ", 
		limits=((log10(0.1r_h), log10(100r_h)), (-1, 1)))
	
	scatter!(log10.(r), nm.log10.(ν_s_nbody) .- nm.log10.(ν_s), label="")
	hlines!([0], color="black", label="")

	linkxaxes!(ax, ax2, )
	rowsize!(fig.layout, 2, Auto(0.3))
	hidexdecorations!(ax, grid=false)
	fig
end

# ╔═╡ 76200404-16aa-4caf-b247-3bc330b82868
function calc_r_h(rs, masses)
	_sort = sortperm(rs)
	M_in = cumsum(masses[_sort])
	M_in ./= M_in[end]
	idx_h = findfirst(M_in .> 0.5)
	return rs[idx_h]
end

# ╔═╡ 7f7d8cb9-761c-4f30-a336-ab5657144961
let
	r = lguys.calc_r(snap, [cen.x, cen.y, cen.z])
	ms = ps_all[snap.index]
	
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
	x = snap.positions[1, :] .- cen.x
	y = snap.positions[2, :] .- cen.y
	R = @. sqrt(x^2 + y^2)

	prof = lguys.calc_properties(R, weights=ps_all[snap.index], bins=50)

	
	fig = Figure()
	ax = Axis(fig[1,1], ylabel=L"\log \Sigma", 
		xlabel="R / kpc",
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

# ╔═╡ 3b585ba2-b9a1-45c4-b9d1-4fbdbd8ccfe6
x_sun = [8.122, 0, 0]

# ╔═╡ 87b5b241-db72-45ee-b3a7-a394f99510d9
shift_vec = x_sun .+ lguys.rand_unit() * distance

# ╔═╡ ff51f97d-2404-49c6-9339-4b201d6a94a9
let
	ms = ps_all[snap.index]
	p_min = 1e-7
	filt = p_min * maximum(ms) .< ms

	ms = ms[filt]
	snap_shift = copy(snap[filt])
	snap_shift.positions .+= shift_vec

	obs = lguys.to_sky(snap_shift)
	ra = [o.ra for o in obs]
	dec = [o.dec for o in obs]

	ra0, dec0 = lguys.calc_centre2D(ra, dec, "mean", ms)
	xi, eta = lguys.to_tangent(ra, dec, ra0, dec0)
	R = @. 60sqrt(xi^2 + eta^2)
	
	prof = lguys.calc_properties(R, weights=ms, bins=50)


	fig = Figure()
	ax = Axis(fig[1,1], ylabel=L"\log \Sigma", 
		xlabel="R / arcmin",
		limits=((-1, 2.5), (-15, 3))
		)

	errscatter!(prof.log_r, prof.log_Sigma, yerr=prof.log_Sigma_err)

	
	profile2 = lguys.Exp2D(R_s = R_s_arcmin)

	log_Σ(r) = log10(lguys.calc_Σ(profile2, r))

	log_R = LinRange(-2, 2, 1000)
	y = log_Σ.(10 .^ log_R)
	
	lines!(log_R, y)
	fig
end

# ╔═╡ Cell order:
# ╟─17ffde4b-5796-4915-9741-d594cf0c5ca7
# ╠═a893932c-f184-42bc-9a0e-0960f10520aa
# ╠═641946b3-e6f2-4d6d-8777-7698f353eb3d
# ╟─93045024-a91d-4b31-9a5a-7c999afdb9ec
# ╠═530c6c09-4454-4952-8351-dccbb4ed429f
# ╠═48ce69f2-09d5-4166-9890-1ab768f3b59f
# ╠═7809e324-ba5f-4520-b6e4-c7727c227154
# ╠═9326c8a6-8b9b-4406-b00f-9febb3dcca46
# ╠═16e0729e-9a75-458e-a56c-73967c819c31
# ╠═d88cbe6a-b87b-45a2-8a3f-c5ef0b7a8935
# ╠═46be1f99-6f64-476e-84cb-11d4b6504a86
# ╠═aa69bde5-ab93-4105-9d48-ad0ace88e3f0
# ╠═f7f746c9-cd03-4391-988b-dffeb31b2842
# ╠═855fd729-22b6-4845-9d2b-e796d4a15811
# ╠═4cb09115-143d-456f-9c6a-19656f638677
# ╠═d77557e5-f7d8-40e9-ae40-a4b6b8df16cd
# ╠═e37bf6d7-9445-49bf-8333-f68ad25436b2
# ╠═cf50c4b4-a634-4310-8aa5-91ab653313a9
# ╠═d43a1be5-7ecd-474b-9580-1191b669ff24
# ╠═0f5671e6-deb4-11ee-3178-4d3f920f23a2
# ╠═eb5ffa80-6959-4a2c-980f-f818f6a03c14
# ╠═8af6b094-ea45-4f01-abb4-9b4686e81719
# ╠═f79414b4-620e-440e-a421-7bc13d373546
# ╠═45acc05f-85a8-4bbc-bd43-a34583c983b3
# ╠═23158c79-d35c-411a-a35a-950f04214e19
# ╠═771413c9-3edf-4687-b4cd-3a60ae0ceda2
# ╠═b3663c19-c029-4a1e-ab82-a05177e3a5d0
# ╠═36b4adbd-d706-4e72-a922-53080c67946c
# ╠═7a39cd4f-9646-4969-9410-b093bca633cb
# ╠═bd1bca1d-0982-47d8-823e-eadc05191b88
# ╠═1fab14ec-6bfc-4243-bc48-915d1a129925
# ╠═dfa675d6-aa32-45c3-a16c-626e16f36083
# ╠═830faf5d-da17-4939-ad7c-080503b66990
# ╠═de614f3d-b3e4-41a0-886b-594b0c93cc5b
# ╠═bbd02922-37fc-459d-8702-451127d2b9a2
# ╠═20f858d6-a9f5-4880-a431-60b395cc7e50
# ╟─126c6825-723f-4d13-b5a3-64cba72fc867
# ╠═8bb8736d-a41b-4dac-a6cd-06d0d4704654
# ╠═c537fdc7-4778-40e6-ab82-7e3151394473
# ╠═46cf718e-79d5-4017-b113-72efa692a014
# ╠═b53887a2-0b14-4d95-bdf5-d70960cfe736
# ╠═c2e070b8-6f15-48e1-8ff3-8145f67b2001
# ╠═837dce1c-fadb-42df-92f9-49a1bd859d1f
# ╠═e9129b04-ffca-4b10-8cda-4955841a7655
# ╠═fdf4f368-962a-4c3a-9571-8489740c31b0
# ╠═c9098d39-3573-450b-8c5e-a6d847c85a3e
# ╠═162e0102-ed49-4648-acf9-c78e52d41ed3
# ╟─7481f47a-2b8a-45a3-9b4e-31ea14e9d331
# ╠═b625d8a5-7265-4849-9bd6-ca8064d392eb
# ╠═78ce5a98-fd3f-4e39-981f-2bea58b117bf
# ╟─8b66d00d-529b-4e8c-9386-b17117996579
# ╠═75d23b44-71e7-4e28-ad3e-c537f3d4422f
# ╠═49c6372d-4ac4-46e2-aace-4cae0564f9b5
# ╠═7409a024-cea4-47a5-84d2-846c96d88b7a
# ╠═84fdc265-988c-40db-87e5-44ba55d0e412
# ╠═a389fde4-f05f-48ff-9a8c-8eb0643a849d
# ╠═3b229c8e-9320-4c07-b948-c34a0c082341
# ╟─daf54cb4-06a3-4a8e-a533-354ca8740aec
# ╠═77e2c1e3-7756-4ab7-810a-03ccdc635aa1
# ╠═5b9c7242-5d4b-4f5a-83e5-3e4d11017aa5
# ╠═b1a34b9b-c6ca-4818-95b2-5b55fb32511e
# ╠═a5bc5ce3-8e33-4514-bc2d-4b4299f104f9
# ╠═33a26663-0c08-411b-902b-a509b2afa5ad
# ╠═a73db457-2cc2-4917-bb25-0429f4daecbd
# ╠═0e7a922f-4e61-4d02-ac28-3915f5d1c9da
# ╠═7d69e393-c4db-4ff0-ab5a-85ac50c785c2
# ╠═6fba7fa7-9a50-4379-b376-5c07f3638411
# ╟─b67a17dc-7a59-4b62-9120-4c2ada12298c
# ╠═c01a4635-dd57-4793-8893-1bf5fb6996b3
# ╠═131bd303-4f0f-4b24-bd90-5c65cf342f4c
# ╠═52bf618a-5494-4cfb-9fd0-ee82f9682116
# ╠═ec46946a-0bf7-4374-af44-8d8b2ab6a3df
# ╟─06e1b872-ce52-434f-a8d1-3b0a5055eed2
# ╠═90856551-fff8-4d15-be66-6353091b5e50
# ╠═9064b0ca-4d26-4cb8-bbc8-353da44ffc26
# ╠═cc231e78-cfc0-4876-9f8b-980139f7d27f
# ╟─ee866164-c6f2-4f70-bde1-360abd5fd80e
# ╠═a52e5e94-6068-4545-962f-e02a485b62f5
# ╟─a3071af2-beff-408d-b109-d4f289f8f7f4
# ╠═a483545a-42ea-474b-a2ca-6f1bd4c7275b
# ╠═7bc02274-aa44-4609-b9a6-e409de5172af
# ╠═f66a4a6d-b31c-481b-bf7a-4801c783ceb4
# ╠═c18bc8b9-7ea1-47e2-8d7e-971a0943917a
# ╠═92798e2e-2356-43c2-a4a9-82a70619a5f5
# ╠═4671b864-470d-4181-993b-4e64d5687460
# ╠═bd8489da-ac53-46c6-979e-06d5dc6e25d1
# ╠═f42cb7e1-64b7-47da-be05-dd50c2471fb3
# ╠═793f701e-376d-4ab3-ae74-6abe29e3c3ae
# ╠═a9335e17-a410-455a-9a9e-d63706a026bd
# ╠═76200404-16aa-4caf-b247-3bc330b82868
# ╠═7f7d8cb9-761c-4f30-a336-ab5657144961
# ╠═a2f72082-7145-42be-9f40-e00d18deb267
# ╟─bf8305f4-a5b8-4c79-8a01-e2aa18e4a5c5
# ╠═6dd92ee1-374d-47fa-ad61-b54764b23240
# ╟─99f274b4-91f3-41d0-b7d3-234badeb43d1
# ╠═4396bced-cae8-4107-ac83-48cc3c4146f2
# ╠═062cea41-3b51-474e-b877-7a4d96813fbc
# ╠═3b585ba2-b9a1-45c4-b9d1-4fbdbd8ccfe6
# ╠═87b5b241-db72-45ee-b3a7-a394f99510d9
# ╠═ff51f97d-2404-49c6-9339-4b201d6a94a9
