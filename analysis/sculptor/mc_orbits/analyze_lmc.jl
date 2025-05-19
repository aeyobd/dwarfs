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

# ╔═╡ e9e2c787-4e0e-4169-a4a3-401fea21baba
begin 
	import Pkg
	Pkg.activate()
	using CairoMakie
	using DataFrames, CSV

	using Printf
	import LilGuys as lguys

	using Arya
end


# ╔═╡ 1ae6d0aa-34f7-4d73-a52f-28ad4e7ac53f
using PlutoUI

# ╔═╡ cdc87dea-cb1a-4f54-9fce-b5008e6decad
using PyFITS

# ╔═╡ e7ad41b9-2845-4a6d-8fc5-7c35544d766c
using Statistics

# ╔═╡ 5acfcb43-9144-4a07-b983-2ee5096f82fc
include("mc_analysis_utils.jl")

# ╔═╡ 7450144e-5464-4036-a215-b6e2cd270405
md"""
This notebook analyzes the result of the MC samples of orbits in the same potential to determine the plausable range of pericentres and apocentres
"""

# ╔═╡ e0d2af42-c9ad-42b8-b223-440d19ac1138
@bind modelname TextField(24, default="example") |> confirm

# ╔═╡ bb90f6a2-cb19-43f1-ad04-a4b258b5812f
using LilGuys; FIGDIR = "./$modelname/figures"

# ╔═╡ b9f469ed-6e4e-41ee-ac75-1b5bfa0a114a
p_value = 0.0014

# ╔═╡ 46348ecb-ee07-4b6a-af03-fc4f2635f57b


# ╔═╡ e4d8e811-6c9a-4500-bc59-c9aedd94bfbb
extract_orbit = true

# ╔═╡ 581df0c4-a509-43d0-beef-c4cc67329d10
md"""
# Setup
"""

# ╔═╡ 9edf5ba2-afc1-422d-98a7-f282869e3253
CairoMakie.activate!(type=:png)

# ╔═╡ 28a113b2-d9a3-4f96-b8d8-e534ff5bd2d0
import TOML

# ╔═╡ 3b83205d-91c1-481e-9305-0d59bc692135
coord_labels = Dict(
	:ra => "ra / degrees",
	:dec => "dec / degrees",
	:pmra => L"$\mu_{\alpha *}$ / mas\,yr$^{-1}$",
	:pmdec => L"$\mu_{\delta}$ / mas\,yr$^{-1}$",
	:radial_velocity => L"$v_\textrm{los}$ / km\,s$^{-1}$",
	:distance => "distance / kpc",
)

# ╔═╡ 8b818798-69fb-481d-ade1-9fd436b1f281
kms_label = L" / km\,s$^{-1}$"

# ╔═╡ 6651d141-f6ca-4e8f-a785-5b14275b367c
T2GYR = lguys.T2GYR

# ╔═╡ 084e5fc4-ccae-41b9-bc35-151d5a67704b
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/sculptor/observed_properties.toml"))

# ╔═╡ fbfece31-31dc-469a-89da-5dbde7e19127
md"""
# Data loading
"""

# ╔═╡ b4f11f71-8857-4a6f-8917-b5e278362b0f
md"""
The code below lets us import just two variables from the sample file. 
These are used to compare the medians of the samples against the input 
means.
"""

# ╔═╡ 302a4ddb-9e4f-434a-9b57-3ee66b6bf528
obs = ICRS(obs_props)

# ╔═╡ f782f563-a771-4652-9406-2bfeb7c624d9
err = get_coord_err(obs_props)

# ╔═╡ 88536e86-cf2a-4dff-ae64-514821957d40
md"""
warning, the simulation output and peris_apos.fits must analyze the same simulation. Checking the individual orbits below should confirm this and the loading code should fail if there are any issues.
"""

# ╔═╡ 26d616da-95ec-4fb9-b9a8-2f095d74c722
"""
	sort_snap(snap)

returns the snapshot sorted by the index
"""
function sort_snap(snap)
	return snap[sortperm(snap.index)]
end

# ╔═╡ 9a22d47b-8474-4596-b418-de33eb07c627
begin 
	out = lguys.Output(modelname);

	df_peris_apos = read_fits("$modelname/peris_apos.fits")
	snap = out[1] |> sort_snap

	@assert all(snap.index  .== df_peris_apos.index) "snapshot and peri apo index must match"
end

# ╔═╡ 73af5376-3b38-4d5f-b6e2-1cafba27dabd
begin
	# loads in trajectory of lmc in Vasiliev 2021
	lmc_file = "$modelname/lmc_traj.csv"
	lmc_traj = CSV.read(lmc_file, DataFrame)
	
	lmc_pos = [lmc_traj.x lmc_traj.y lmc_traj.z]'
	lmc_vel = [lmc_traj.v_x lmc_traj.v_y lmc_traj.v_z]'

end

# ╔═╡ cd2c1a82-8b0b-4abc-8831-f621e5ec483a
df_peris_apos

# ╔═╡ 4481a5e1-d635-4fea-a5c5-c85f0d6df62f
peris = df_peris_apos.pericentre

# ╔═╡ 384be6a6-f9d9-47e0-9792-aef6689dcbdb
apos = df_peris_apos.apocentre

# ╔═╡ a1b2edd5-20e4-456f-a8e5-4634b0a433c2
perilmc = df_peris_apos.peri_lmc

# ╔═╡ b2fffcba-c474-430c-8100-037a0d352d69
Δt_lmc = df_peris_apos.t_last_peri_lmc * T2GYR

# ╔═╡ 850c29d0-f395-448f-87be-d788210594c2
Δt_mw = df_peris_apos.t_last_peri * T2GYR

# ╔═╡ 4df80961-f9c2-49b2-9cef-737670277d84
md"""
## quantile properties & orbit selection
"""

# ╔═╡ 7df6158e-3d20-4995-b506-a621e7af07ca
md"""
We need to select orbits systematically from the probability cloud. Asya's method was to simply take the orbits from the MCMC samples with pericentres closest to specified quantiles. I use this as an initial test, but prefer to take the median values in some percentile rage, such that the orbit tends to be closer to the middle of the cloud.
"""

# ╔═╡ 236c7f58-01da-4e33-91f0-8d578a269b30
quantiles = [0.5, p_value, 1-p_value]

# ╔═╡ e5f728b8-8412-4f57-ad38-a0a35bb08a48
orbit_labels = ["mean", "smallperi", "largeperi"]

# ╔═╡ 413d4e5d-c9cd-4aca-be1e-d132b2bd616d
peri_qs = lguys.quantile(perilmc, quantiles)

# ╔═╡ 17a63cc8-84f4-4248-a7b0-c8378454b1f7
idx = [argmin(abs.(p .- perilmc)) for p in peri_qs]

# ╔═╡ 414c4b71-74f0-4ea6-a472-e3835a0dd9a5
perilmc[idx]

# ╔═╡ bbf3f229-fc3c-46ae-af28-0f8bd81e7d32
perilmc[idx]

# ╔═╡ 3d30c67c-c63a-440b-9a4f-c67a95d3f8bd
md"""
The below functions are for my method of selecting the orbits. The pericentres are filtered to have quantiles between 0.5 and 1.5 times the p_value (between 0.00067 and 0.00202 for current setting) and then the adopted values are those printed out by median_percen. I typically just use the default ra and dec values since the uncertanties are negligable.

Since these values are (likely) not in the random samples, I simply run new orbits for a few selected orbits. These models are named modelnamespecialcases and the associated analyze_with_special_cases.jl makes similar plots to this notebook except shows these special case orbits and their actual trajectories.
"""

# ╔═╡ de2f3380-90df-48f5-ba60-8417e91f4818
function median_residual(observations)
	println("N obs ", length(observations))
	for sym in [:distance, :pmra, :pmdec, :radial_velocity, :ra, :dec]
		md = median(getproperty.(observations, sym))
		@printf "     %-15s  = %7.3f \t \n"  sym md
	end
	println()
	for sym in [:distance, :pmra, :pmdec, :radial_velocity, :ra, :dec]
		md = median(getproperty.(observations, sym))
		res = (md - getproperty(obs, sym) ) / err[string(sym)]
		@printf "Δ ln %-15s  = %6.2f \t \n"  sym res
	end
end

# ╔═╡ 1acef60e-60d6-47ba-85fd-f9780934788b
md"""
# plots
"""

# ╔═╡ 0f3bf339-ea93-44f8-99fa-ded41a3a9980
md"""
## Histograms
"""

# ╔═╡ ca1c236e-795a-408b-845b-9c13bc838619
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "pericentre / kpc",
		ylabel = "count"
	)

	bins, counts, err = lguys.histogram(peris)
	scatter!(lguys.midpoints(bins), counts)

	fig
end

# ╔═╡ 7b7b8357-ebf1-4257-8a38-272b245842cc
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "perilmc / kpc",
		ylabel = "count"
	)

	bins, counts, err = lguys.histogram(perilmc)
	scatter!(lguys.midpoints(bins), counts)

	fig
end

# ╔═╡ c43412d1-3a9b-4cd0-901b-05048c0b07ec
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "delta t MW / Gyr",
		ylabel = "count"
	)

	hist!(Δt_mw)

	fig
end

# ╔═╡ 4cbbebec-f231-4af4-a10f-e01d090ea2ca
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "delta t LMC / Gyr",
		ylabel = "count"
	)

	hist!(Δt_lmc)

	fig
end

# ╔═╡ 46b4242b-8af7-4233-8ecf-d86740b4c884
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "apocentre / kpc",
		ylabel = "count"
	)

	bins, counts, err = lguys.histogram(apos)
	scatter!(lguys.midpoints(bins), counts)

	fig
end

# ╔═╡ 471a5501-bc7b-4114-be1a-9088b4e674cd
hist(lguys.radii(snap.positions),
	axis=(; xlabel="initial galactocentric radius / kpc")
)

# ╔═╡ 68b0383a-3d5a-4b94-967c-f0e31e8a0ce1
hist(lguys.radii(snap.velocities) * lguys.V2KMS,
	axis=(; xlabel="initial galactocentric velocity / km/s")
)

# ╔═╡ aa8ed2ce-472c-4b91-b371-219ce5b63774


# ╔═╡ 5f11f6ab-c9ab-4600-acca-a0bb84d81a12
begin
	points = [lguys.Galactocentric(
		snap.positions[:, i]*lguys.R2KPC, 
		-snap.velocities[:, i]*lguys.V2KMS)
		for i in 1:length(snap)]
	
	observations = lguys.transform.(lguys.ICRS, points)
end


# ╔═╡ 13f03c68-0edb-418b-b1f7-7a7f24659be2
let 
	peri1 = lguys.quantile(perilmc, 0.5p_value)
	peri2 = lguys.quantile(perilmc, 1.5p_value)

	idx_s = @. peri1 < perilmc < peri2

	println("peri lmc ~ ", median(perilmc[idx_s]))
	median_residual(observations[idx_s])
end

# ╔═╡ 34104429-05e0-40a6-83e5-078dbe346504
let
	peri2 = lguys.quantile(perilmc, 1 - 0.5p_value)
	peri1 = lguys.quantile(perilmc, 1 - 1.5p_value)

	idx_s = @. peri1 < perilmc < peri2
	println("peri lmc ~ ", median(perilmc[idx_s]))

	median_residual(observations[idx_s])
end

# ╔═╡ 4ee33ce2-c00a-4fcf-b7fc-b78c1677c9e4
let 
	peri1 = lguys.quantile(peris, 0.5p_value)
	peri2 = lguys.quantile(peris, 1.5p_value)

	idx_s = @. peri1 < peris < peri2

	println("peri MW ~ ", median(peris[idx_s]))
	median_residual(observations[idx_s])
end

# ╔═╡ e5825c4a-b446-44a3-8fd5-d94664965aca
println("sample median"); median_residual(observations)

# ╔═╡ 92aac8e8-d599-4a1e-965e-e653bc54c509
dists = getproperty.(observations, :distance)

# ╔═╡ ede3836c-740d-4ac7-bbc7-3165981a1878
normal_dist(x, μ, σ) = 1/√(2π) * 1/σ * exp(-(x-μ)^2/2σ^2)

# ╔═╡ 44660b2f-6220-473b-bb2f-07e23b176491
columns = [:ra, :dec, :distance, :pmra, :pmdec, :radial_velocity]

# ╔═╡ d3063e30-2cb3-4f1b-8546-0d5e81d90d9f
let
	
	for sym in columns[1:2]
	    x = [getproperty(o, sym) for o in observations]
	    y = peris

		
	    p = scatter(x, y, alpha=0.1,
			axis=(; xlabel=coord_labels[sym], ylabel="pericentre / kpc")
		)
	    @info p 
	end

end

# ╔═╡ 99cbc08e-5b0e-4748-9727-667e0eca2f74
let
	fig = Figure(size=(600, 600))
	plot_kwargs = Dict(
		:color => :black,
		:alpha => 0.1,
		:markersize => 1,
	)

	
	orbit_points_kwargs = [Dict(
		:color => COLORS[i],
		:alpha => 1,
		:markersize => 10,
		:label => orbit_labels[i]
	) for i in eachindex(idx)
		]
	
	ax_kwargs = Dict(
		:xgridvisible => false,
		:ygridvisible => false,
		:ylabel => "perilmc / kpc",
	)

	ax_idx = Dict(
		:pmra => [1, 1],
		:pmdec => [1, 2],
		:distance => [2, 1],
		:radial_velocity => [2, 2],
	)

	for sym in [:pmra, :pmdec, :distance, :radial_velocity]
		ax = Axis(fig[ax_idx[sym]...],
			xlabel=coord_labels[sym];
			ax_kwargs...
		)
		x = [getproperty(o, sym) for o in observations]
		scatter!(x, perilmc; plot_kwargs...)
		for i in eachindex(idx)
			scatter!(x[idx[i]], perilmc[idx[i]]; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)
	@savefig "perilmc_corr"

	fig
end

# ╔═╡ 456acc60-e7e1-428e-bd5d-bc72a18a146f
let
	fig = Figure(size=(600, 600))
	plot_kwargs = Dict(
		:color => :black,
		:alpha => 0.1,
		:markersize => 1,
	)

	
	orbit_points_kwargs = [Dict(
		:color => COLORS[i],
		:alpha => 1,
		:markersize => 10,
		:label => orbit_labels[i]
	) for i in eachindex(idx)
		]
	
	ax_kwargs = Dict(
		:xgridvisible => false,
		:ygridvisible => false,
		:ylabel => "t last LMC peri / Gyr",
	)

	ax_idx = Dict(
		:pmra => [1, 1],
		:pmdec => [1, 2],
		:distance => [2, 1],
		:radial_velocity => [2, 2],
	)

	for sym in [:pmra, :pmdec, :distance, :radial_velocity]
		ax = Axis(fig[ax_idx[sym]...],
			xlabel=coord_labels[sym];
			ax_kwargs...
		)
		x = [getproperty(o, sym) for o in observations]
		scatter!(x, Δt_lmc; plot_kwargs...)
		for i in eachindex(idx)
			scatter!(x[idx[i]], Δt_lmc[idx[i]]; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)
	@savefig "t_last_perilmc"

	fig
end

# ╔═╡ c48b4e73-480e-4a50-b5fc-db5f6c5b040e
let
	fig = Figure(size=(600, 600))
	plot_kwargs = Dict(
		:color => :black,
		:alpha => 0.1,
		:markersize => 1,
	)

	
	orbit_points_kwargs = [Dict(
		:color => COLORS[i],
		:alpha => 1,
		:markersize => 10,
		:label => orbit_labels[i]
	) for i in eachindex(idx)
		]
	
	ax_kwargs = Dict(
		:xgridvisible => false,
		:ygridvisible => false,
		:ylabel => "pericentre / kpc",
	)

	ax_idx = Dict(
		:pmra => [1, 1],
		:pmdec => [1, 2],
		:distance => [2, 1],
		:radial_velocity => [2, 2],
	)

	for sym in [:pmra, :pmdec, :distance, :radial_velocity]
		ax = Axis(fig[ax_idx[sym]...],
			xlabel=coord_labels[sym];
			ax_kwargs...
		)
		x = [getproperty(o, sym) for o in observations]
		scatter!(x, peris; plot_kwargs...)
		for i in eachindex(idx)
			scatter!(x[idx[i]], peris[idx[i]]; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)


	@savefig "peri_mc_orbits_corr"
	fig
end

# ╔═╡ e9311003-b6af-4e09-907f-c3388d59287d
let
	fig = Figure(size=(600, 600))
	plot_kwargs = Dict(
		:color => :black,
		:alpha => 0.1,
		:markersize => 1,
	)

	
	orbit_points_kwargs = [Dict(
		:color => COLORS[i],
		:alpha => 1,
		:markersize => 10,
		:label => orbit_labels[i]
	) for i in eachindex(idx)
		]
	
	ax_kwargs = Dict(
		:xgridvisible => false,
		:ygridvisible => false,
		:ylabel => "t last peri",
	)

	ax_idx = Dict(
		:pmra => [1, 1],
		:pmdec => [1, 2],
		:distance => [2, 1],
		:radial_velocity => [2, 2],
	)

	for sym in [:pmra, :pmdec, :distance, :radial_velocity]
		ax = Axis(fig[ax_idx[sym]...],
			xlabel=coord_labels[sym];
			ax_kwargs...
		)
		x = [getproperty(o, sym) for o in observations]
		scatter!(x, Δt_mw; plot_kwargs...)
		for i in eachindex(idx)
			scatter!(x[idx[i]], Δt_mw[idx[i]]; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)

	fig
end

# ╔═╡ 43d43f63-4c13-4b23-950e-ada59aa86bc9
let
	
	for sym in [:ra, :dec]
	    x = [getproperty(o, sym) for o in observations]
	    y = perilmc

		
	    p = scatter(x, y, alpha=0.1, 
			axis=(; xlabel=String(sym), ylabel="apocentre / kpc")
		)
	    @info p 
	end

end

# ╔═╡ 16f4ac20-d8cf-4218-8c01-c15e04e567fb
md"""
# The selected orbits
"""

# ╔═╡ d31f91e8-6db6-4771-9544-8e54a816ecc1
if extract_orbit
	positions = [lguys.extract_vector(out, :positions, i) for i in idx]
	velocities = [lguys.extract_vector(out, :velocities, i) for i in idx]
	accelerations = [lguys.extract_vector(out, :accelerations, i) for i in idx]
end

# ╔═╡ d3de3a33-8f9a-44a5-9851-1992c88f04ae
if extract_orbit
	Φs_ext = [lguys.extract(out, :potential_ext, i) for i in idx]
	Φs = [lguys.extract(out, :potential, i) for i in idx]
end

# ╔═╡ 5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
if extract_orbit
	rs = lguys.radii.(positions)
	vs = lguys.radii.(velocities)
	accs = lguys.radii.(accelerations)
	rs_lmc = lguys.radii.(positions, [lmc_pos])
end

# ╔═╡ e5d40e2f-ac47-4827-853d-2f94bc39a624
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="Scl—MW distance / kpc",
		limits=(nothing, nothing, 0, nothing)
	)

	for i in eachindex(idx)
		lines!(-out.times * lguys.T2GYR, rs[i], label=orbit_labels[i])
	
		#hlines!([peris[idx[i]], apos[idx[i]]], linestyle=:dot)
	end

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 06791b4d-8e39-4057-8b8b-a530ab96d30d
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="radius / kpc",
		limits=(nothing, nothing, 0, 100)
	)

	for i in eachindex(idx)
		lines!(-out.times * lguys.T2GYR, rs_lmc[i], label=orbit_labels[i])
	
		#hlines!([peris[idx[i]], apos[idx[i]]], linestyle=:dot)
	end

	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
lguys.plot_xyz(positions..., labels=orbit_labels)

# ╔═╡ aec104ef-a4ca-4d69-aa72-1657c83c037d
lguys.plot_xyz([p .- lmc_pos for p in positions]..., labels=orbit_labels)

# ╔═╡ 34efe422-e672-4d17-a558-ce32fb704a8e
lguys.plot_xyz(velocities..., units=" / km / s")

# ╔═╡ 756f6d16-bc47-47d8-a090-32b4d7b3f3a0
times = out.times * lguys.T2GYR

# ╔═╡ 52e8c357-1269-4165-a6d7-d7e3c129f123
idx_perilmc = [argmin(r) for r in rs_lmc]

# ╔═╡ c82d023f-8123-4b20-ace1-cfb48dd94a5e
idx_peri = [argmin(r) for r in rs]

# ╔═╡ 385eb53c-011b-431e-993b-ca12fece56e2
let
	fig = Figure()
	ax_acc = Axis(fig[1, 1], 
		yticksmirrored = false,
		xlabel="time / Gyr", 
		ylabel="Scl acceleration"
	)
	ax_dist = Axis(fig[1, 1], yaxisposition=:right,
		yticksmirrored = false,
		ylabel="distance",
		
	)
	
	hidexdecorations!(ax_dist)
	linkxaxes!(ax_acc, ax_dist)

	i = 1
	lines!(ax_acc, out.times * lguys.T2GYR, accs[1],
		color="black"
	)

	lines!(ax_dist, out.times * lguys.T2GYR, rs[1],
	)
	vlines!(ax_dist, times[idx_peri[i]], alpha=0.3)

	lines!(ax_dist, out.times * lguys.T2GYR, rs_lmc[1],
	)

	vlines!(ax_dist, times[idx_perilmc[i]], alpha=0.3)

	
	lguys.hide_grid!(ax_acc)
	lguys.hide_grid!(ax_dist)
	
	fig
end

# ╔═╡ 09bbae0d-ca3e-426d-b77c-69dd68ca42cc
begin
	M_b = 115
	r_b = 20
	
	a_exp(r) =  M_b / (r+r_b)^2
	phi_exp(r) = -M_b / (r+r_b)
end

# ╔═╡ 0fc3644a-edd6-4bf2-8fb5-8ad91fd779d8
md"""
## Additional info
"""

# ╔═╡ de1e5245-0946-47cd-8e2c-ba1914cfeb74
begin 
	# orbit info
	for i in 1:length(idx)
		t =1
		@printf "orbit: \t\t %i\n" i
		@printf "index: \t\t %i\n" idx[i]
		
		@printf "pericentre:\t %0.1f\n" peris[idx[i]]
		@printf "apocentre: \t %0.1f\n" apos[idx[i]]

		@printf "time of first apocentre: %f \n" out.times[end] - out.times[t]
		@printf "radius of first apocentre: %f\n" rs[i][t]
		@printf "intial position: [%f, %f, %f]\n" positions[i][:, t]...
		@printf "intial velocity: [%f, %f, %f]\n" -lguys.V2KMS* velocities[i][:, t]...
		@printf "final position: [%f, %f, %f]\n" positions[i][:, 1]...
		@printf "final velocity: [%f, %f, %f]\n" -lguys.V2KMS * velocities[i][:, 1]...

		o = observations[idx[i]]
		@printf "ra: \t %f\n" o.ra
		@printf "dec: \t %f\n" o.dec
		@printf "pmra: \t %f\n" o.pmra
		@printf "pmdec: \t %f\n" o.pmdec
		@printf "dist: \t %f\n" o.distance
		@printf "rv: \t %f\n" o.radial_velocity

		println()
	end
end

# ╔═╡ Cell order:
# ╟─7450144e-5464-4036-a215-b6e2cd270405
# ╠═e0d2af42-c9ad-42b8-b223-440d19ac1138
# ╠═b9f469ed-6e4e-41ee-ac75-1b5bfa0a114a
# ╠═414c4b71-74f0-4ea6-a472-e3835a0dd9a5
# ╠═46348ecb-ee07-4b6a-af03-fc4f2635f57b
# ╠═e4d8e811-6c9a-4500-bc59-c9aedd94bfbb
# ╟─581df0c4-a509-43d0-beef-c4cc67329d10
# ╠═9edf5ba2-afc1-422d-98a7-f282869e3253
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═1ae6d0aa-34f7-4d73-a52f-28ad4e7ac53f
# ╠═28a113b2-d9a3-4f96-b8d8-e534ff5bd2d0
# ╠═cdc87dea-cb1a-4f54-9fce-b5008e6decad
# ╠═bb90f6a2-cb19-43f1-ad04-a4b258b5812f
# ╠═e7ad41b9-2845-4a6d-8fc5-7c35544d766c
# ╠═3b83205d-91c1-481e-9305-0d59bc692135
# ╠═8b818798-69fb-481d-ade1-9fd436b1f281
# ╠═6651d141-f6ca-4e8f-a785-5b14275b367c
# ╠═084e5fc4-ccae-41b9-bc35-151d5a67704b
# ╠═5acfcb43-9144-4a07-b983-2ee5096f82fc
# ╟─fbfece31-31dc-469a-89da-5dbde7e19127
# ╟─b4f11f71-8857-4a6f-8917-b5e278362b0f
# ╠═302a4ddb-9e4f-434a-9b57-3ee66b6bf528
# ╠═f782f563-a771-4652-9406-2bfeb7c624d9
# ╟─88536e86-cf2a-4dff-ae64-514821957d40
# ╠═26d616da-95ec-4fb9-b9a8-2f095d74c722
# ╠═9a22d47b-8474-4596-b418-de33eb07c627
# ╠═73af5376-3b38-4d5f-b6e2-1cafba27dabd
# ╠═cd2c1a82-8b0b-4abc-8831-f621e5ec483a
# ╠═4481a5e1-d635-4fea-a5c5-c85f0d6df62f
# ╠═384be6a6-f9d9-47e0-9792-aef6689dcbdb
# ╠═a1b2edd5-20e4-456f-a8e5-4634b0a433c2
# ╠═b2fffcba-c474-430c-8100-037a0d352d69
# ╠═850c29d0-f395-448f-87be-d788210594c2
# ╟─4df80961-f9c2-49b2-9cef-737670277d84
# ╟─7df6158e-3d20-4995-b506-a621e7af07ca
# ╠═236c7f58-01da-4e33-91f0-8d578a269b30
# ╠═e5f728b8-8412-4f57-ad38-a0a35bb08a48
# ╠═413d4e5d-c9cd-4aca-be1e-d132b2bd616d
# ╠═17a63cc8-84f4-4248-a7b0-c8378454b1f7
# ╠═bbf3f229-fc3c-46ae-af28-0f8bd81e7d32
# ╟─3d30c67c-c63a-440b-9a4f-c67a95d3f8bd
# ╠═de2f3380-90df-48f5-ba60-8417e91f4818
# ╠═13f03c68-0edb-418b-b1f7-7a7f24659be2
# ╠═34104429-05e0-40a6-83e5-078dbe346504
# ╠═4ee33ce2-c00a-4fcf-b7fc-b78c1677c9e4
# ╠═e5825c4a-b446-44a3-8fd5-d94664965aca
# ╟─1acef60e-60d6-47ba-85fd-f9780934788b
# ╟─0f3bf339-ea93-44f8-99fa-ded41a3a9980
# ╠═ca1c236e-795a-408b-845b-9c13bc838619
# ╠═7b7b8357-ebf1-4257-8a38-272b245842cc
# ╠═c43412d1-3a9b-4cd0-901b-05048c0b07ec
# ╠═4cbbebec-f231-4af4-a10f-e01d090ea2ca
# ╠═46b4242b-8af7-4233-8ecf-d86740b4c884
# ╠═92aac8e8-d599-4a1e-965e-e653bc54c509
# ╠═471a5501-bc7b-4114-be1a-9088b4e674cd
# ╠═68b0383a-3d5a-4b94-967c-f0e31e8a0ce1
# ╠═aa8ed2ce-472c-4b91-b371-219ce5b63774
# ╠═5f11f6ab-c9ab-4600-acca-a0bb84d81a12
# ╠═ede3836c-740d-4ac7-bbc7-3165981a1878
# ╠═44660b2f-6220-473b-bb2f-07e23b176491
# ╠═d3063e30-2cb3-4f1b-8546-0d5e81d90d9f
# ╟─99cbc08e-5b0e-4748-9727-667e0eca2f74
# ╟─456acc60-e7e1-428e-bd5d-bc72a18a146f
# ╟─c48b4e73-480e-4a50-b5fc-db5f6c5b040e
# ╟─e9311003-b6af-4e09-907f-c3388d59287d
# ╠═43d43f63-4c13-4b23-950e-ada59aa86bc9
# ╟─16f4ac20-d8cf-4218-8c01-c15e04e567fb
# ╠═d31f91e8-6db6-4771-9544-8e54a816ecc1
# ╠═d3de3a33-8f9a-44a5-9851-1992c88f04ae
# ╠═5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
# ╟─e5d40e2f-ac47-4827-853d-2f94bc39a624
# ╠═06791b4d-8e39-4057-8b8b-a530ab96d30d
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
# ╠═aec104ef-a4ca-4d69-aa72-1657c83c037d
# ╠═34efe422-e672-4d17-a558-ce32fb704a8e
# ╠═756f6d16-bc47-47d8-a090-32b4d7b3f3a0
# ╠═385eb53c-011b-431e-993b-ca12fece56e2
# ╠═52e8c357-1269-4165-a6d7-d7e3c129f123
# ╠═c82d023f-8123-4b20-ace1-cfb48dd94a5e
# ╠═09bbae0d-ca3e-426d-b77c-69dd68ca42cc
# ╟─0fc3644a-edd6-4bf2-8fb5-8ad91fd779d8
# ╠═de1e5245-0946-47cd-8e2c-ba1914cfeb74
