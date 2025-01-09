### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ e9e2c787-4e0e-4169-a4a3-401fea21baba
begin 
	import Pkg
	Pkg.activate()
	using CairoMakie
	using DataFrames, CSV

	using Printf
	using Revise # so we can overwrite file if anything happens
	import LilGuys as lguys

	using Arya
end


# ╔═╡ d975d00c-fd69-4dd0-90d4-c4cbe73d9754
using Statistics, Distributions

# ╔═╡ 7450144e-5464-4036-a215-b6e2cd270405
md"""
This notebook analyzes the result of the MC samples of orbits in the same potential to determine the plausable range of pericentres and apocentres.
In general, this notebook is meant to plot and calculate main properties.
If you would like to investigate a special case further, it is likely best to make a new notebook in the target analysis directory.
Additionally, see analyze_lmc.jl in this directory for a version which also plots additional plots for an MW-LMC potential.
"""

# ╔═╡ 2b9d49c6-74cc-4cce-b29e-04e94776863f
md"""
The most important variable is to set the modelname to the appropriate directory.
"""

# ╔═╡ 6ca3fe17-3f13-43fe-967b-881078135ead
modelname = "vasiliev24_L3M11_plmc_err"

# ╔═╡ 46348ecb-ee07-4b6a-af03-fc4f2635f57b
fig_dir = "./$modelname/figures"

# ╔═╡ b9f469ed-6e4e-41ee-ac75-1b5bfa0a114a
p_value = 0.001349898031630093 # 3sigma

# ╔═╡ 7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
md"""
## Setup
"""

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

# ╔═╡ 8f70add4-effe-437d-a10a-4e15228f9fec
md"""
# Data loading
"""

# ╔═╡ b15fb3ac-2219-419a-854a-31a783acf891
md"""
The code below lets us import the obs variable from the sample.jl file in the simulation directory.
This is primarily used for tests later.
"""

# ╔═╡ 9583b7d0-0d86-4346-998b-000ea68e94b6
module SampleSetup
	import ..modelname 
	include(joinpath(modelname, "simulation/sample.jl"))
end

# ╔═╡ fa790a4d-e74f-479b-8ff6-aa2f23cb573d
obs = SampleSetup.obs

# ╔═╡ eff58c52-a32b-4faa-9b98-c8234d9b21fc
err = SampleSetup.err

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

	df_peris_apos = lguys.read_fits("$modelname/peris_apos.fits")
	snap = out[1] |> sort_snap

	@assert all(snap.index  .== df_peris_apos.index) "snapshot and peri apo index must match"
end

# ╔═╡ da6e5566-f2df-4feb-9188-53eca9a1a0d5
df_peris_apos

# ╔═╡ 384be6a6-f9d9-47e0-9792-aef6689dcbdb
apos = df_peris_apos.apocentre

# ╔═╡ 4481a5e1-d635-4fea-a5c5-c85f0d6df62f
peris = df_peris_apos.pericentre

# ╔═╡ 392315ee-72d2-4a14-9afc-5fd6424b3e83
md"""
## quantile properties & orbit selection
"""

# ╔═╡ 950e0210-e4fe-4bad-b82f-69247fd0edd8
md"""
We need to select orbits systematically from the probability cloud. 
Asya's method was to simply take the orbits from the MCMC samples with pericentres
closest to specified quantiles. 
I use this as an initial test, but prefer to take the 
median values in some percentile rage, such that the orbit tends to be closer to the middle of the cloud.
"""

# ╔═╡ 4a4b8b73-92c3-4e1f-93d8-e44369b8f148
quantiles = [0.5, p_value, 1-p_value]

# ╔═╡ e5f728b8-8412-4f57-ad38-a0a35bb08a48
orbit_labels = ["mean", "smallperi", "largeperi"]

# ╔═╡ 413d4e5d-c9cd-4aca-be1e-d132b2bd616d
peri_qs = lguys.quantile(peris, quantiles)

# ╔═╡ 17a63cc8-84f4-4248-a7b0-c8378454b1f7
idx = [argmin(abs.(p .- peris)) for p in peri_qs]

# ╔═╡ bbf3f229-fc3c-46ae-af28-0f8bd81e7d32
peris[idx]

# ╔═╡ d54733e9-7552-4d9c-b8e8-670433469385
md"""
The below functions are for my method of selecting the orbits.
The pericentres are filtered to have quantiles between 0.5 and 1.5 times the p\_value
(between $(round(0.5p_value, digits=5)) and $(round(1.5p_value, digits=5)) for current setting)
and then the adopted values are those printed out by median_percen.
I typically just use the default ra and dec values since the uncertanties are negligable.

Since these values are (likely) not in the random samples, I simply run new orbits for a few selected orbits. These models are named modelname_special_cases and the associated `analyze_with_special_cases.jl` makes similar plots to this notebook except shows these special case orbits and their actual trajectories.
"""

# ╔═╡ de2f3380-90df-48f5-ba60-8417e91f4818
function median_residual(observations)
	for sym in [:ra, :dec, :distance, :pmra, :pmdec, :radial_velocity]
		md = median(getproperty.(observations, sym))
		res = (md - getproperty(obs, sym) ) / getproperty(err, sym)
		@printf "Δ ln %-15s  = %6.2f \t \n"  sym res
	end
end

# ╔═╡ c8aec4f8-975f-4bbc-b874-bf0172d35868
function median_percen(observations)
	for sym in [:ra, :dec, :distance, :pmra, :pmdec, :radial_velocity]
		md = median(getproperty.(observations, sym))
		err = std(getproperty.(observations, sym)) / sqrt(length(observations))
		@printf "     %-15s  = %8.3f ± %8.3f \n"  sym md err
	end
end

# ╔═╡ 1acef60e-60d6-47ba-85fd-f9780934788b
md"""
# plots
"""

# ╔═╡ 50baf5a6-fb5b-494e-95f3-53414a9f1cc0
md"""
## Histograms
I have some histograms of orbital properties below, to get a sense of the overal distribution.
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

# ╔═╡ 5a903509-e2cc-4bac-9c59-d1689ccc408e
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel = "time of last pericentre / Gyr",
		ylabel = "count"
	)

	bins, counts, err = lguys.histogram(df_peris_apos.t_last_peri .* lguys.T2GYR)
	scatter!(lguys.midpoints(bins), counts)

	fig
end

# ╔═╡ 471a5501-bc7b-4114-be1a-9088b4e674cd
hist(lguys.calc_r(snap.positions),
	axis=(; xlabel="initial galactocentric radius / kpc",
	ylabel="count")
)

# ╔═╡ 68b0383a-3d5a-4b94-967c-f0e31e8a0ce1
hist(lguys.calc_r(snap.velocities) * lguys.V2KMS,
	axis=(; xlabel="initial galactocentric velocity / km/s",
		ylabel="count"
	)
)

# ╔═╡ 5f11f6ab-c9ab-4600-acca-a0bb84d81a12
begin
	# calculates initial conditions as ICRS coordinates
	
	points = [lguys.Galactocentric(
		snap.positions[:, i]*lguys.R2KPC, 
		-snap.velocities[:, i]*lguys.V2KMS)
		for i in 1:length(snap)]
	
	observations = lguys.transform.(lguys.ICRS, points)
end

# ╔═╡ 4ee33ce2-c00a-4fcf-b7fc-b78c1677c9e4
let 
	idx_s =  peris .< quantile(peris, 2*p_value)

	median_residual(observations[idx_s])
	median_percen(observations[idx_s])
end

# ╔═╡ 34104429-05e0-40a6-83e5-078dbe346504
let
	idx_s =  peris .> quantile(peris, 1-2*p_value)

	median_residual(observations[idx_s])
	median_percen(observations[idx_s])

end

# ╔═╡ e5825c4a-b446-44a3-8fd5-d94664965aca
median_residual(observations)

# ╔═╡ ef57611c-2986-4b03-aa5a-ab45003edd72
median_percen(observations)

# ╔═╡ 92aac8e8-d599-4a1e-965e-e653bc54c509
dists = getproperty.(observations, :distance)

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


	save(joinpath(fig_dir, "peri_mc_orbits_corr.pdf"), fig)
	fig
end

# ╔═╡ 43d43f63-4c13-4b23-950e-ada59aa86bc9
let
	
	for sym in columns
	    x = [getproperty(o, sym) for o in observations]
	    y = peris

		
	    p = scatter(x, y, alpha=0.1, 
			axis=(; xlabel=String(sym), ylabel="apocentre / kpc")
		)
	    @info p 
	end

end

# ╔═╡ 8b6f95f7-284f-4133-b6f1-a22dd9c405f0
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
		:ylabel => "time of last pericentre / Gyr",
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
		scatter!(x, df_peris_apos.t_last_peri * lguys.T2GYR .+ 0.001 .* randn(length(snap)); plot_kwargs...)
		for i in eachindex(idx)
			scatter!(x[idx[i]], df_peris_apos.t_last_peri[idx[i]]* lguys.T2GYR; orbit_points_kwargs[i]...)
		end
	end


	linkyaxes!(fig.content...)


	save(joinpath(fig_dir, "t_last_peri.pdf"), fig)
	fig
end

# ╔═╡ 69e77193-29cc-4304-98a1-44828eaedf9f
md"""
# Validation
"""

# ╔═╡ 89b81ed0-82a1-4a81-a7fd-b6be0644a79d
md"""
The histograms below are for the initial conditions.
Ideally, these should look like gaussians described by the `obs` and `err` variables. If not, something went fairly awry.
"""

# ╔═╡ ede3836c-740d-4ac7-bbc7-3165981a1878
normal_dist(x, μ, σ) = 1/√(2π) * 1/σ * exp(-(x-μ)^2/2σ^2)

# ╔═╡ 6b95d3b2-38db-4376-83b5-8c6e6f1fdfa2
let
	fig = Figure(size=(600, 600))

	ax_kwargs = Dict(
		:xgridvisible => false,
		:ygridvisible => false,
		:ylabel => "density",
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
	    x = getproperty.(observations, sym)
		
		    stephist!(x, bins=50, normalization=:pdf, label="MC samples", color=:black)
	
		
		x_mod = LinRange(minimum(x), maximum(x), 1000)

		mu_exp = getproperty(obs, sym)
		err_exp = getproperty(err, sym)
		y_mod = normal_dist.(x_mod, mu_exp, err_exp)
		lines!(x_mod, y_mod, label="expected gaussian")
			
		axislegend(labelsize=10, padding=(6, 6, 6, 6), patchlabelgap=1, patchsize=(6, 6))
	end

	

	save(joinpath(fig_dir, "peri_mc_orbits_corr.pdf"), fig)
	fig
end

# ╔═╡ ac81acd8-4a78-4230-bc70-3b78a861b618
let

	for sym in [:ra, :dec]
		
		fig = Figure()
		ax = Axis(fig[1,1], 
			xlabel=String(sym),
			ylabel="density",
			#limits=((μ - 5σ, μ + 5σ), nothing),
		)
		
	    x = getproperty.(observations, sym)
		
	    stephist!(x, bins=50, normalization=:pdf, label="MC samples", color=:black)


		μ = getproperty(obs, sym)
		σ = getproperty(err, sym)
		x_mod = LinRange(μ - 3σ, μ + 3σ, 10_000)
		y_mod = normal_dist.(x_mod, μ, σ) 
		lines!(x_mod, y_mod, label="expected")
			
		axislegend()

		@info fig
	end

end

# ╔═╡ 16f4ac20-d8cf-4218-8c01-c15e04e567fb
md"""
# The example orbits
"""

# ╔═╡ d31f91e8-6db6-4771-9544-8e54a816ecc1
begin
	
	positions = [lguys.extract_vector(out, :positions, i) for i in idx]
	velocities = [lguys.extract_vector(out, :velocities, i) for i in idx]
	accelerations = [lguys.extract_vector(out, :accelerations, i) for i in idx]

	Φs_ext = [lguys.extract(out, :Φs_ext, i) for i in idx]
	Φs = [lguys.extract(out, :Φs, i) for i in idx]


end

# ╔═╡ 5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
begin
	rs = lguys.calc_r.(positions)
	vs = lguys.calc_r.(velocities)
	accs = lguys.calc_r.(accelerations)
end

# ╔═╡ 14c36202-66ca-46b3-b282-3895b72311fe
md"""
The plots below are designed to show the special orbits in a variety of frames.
"""

# ╔═╡ e5d40e2f-ac47-4827-853d-2f94bc39a624
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel="Scl–MW distance / kpc"
	)

	for i in eachindex(idx)
		lines!(-out.times * lguys.T2GYR, rs[i], label=orbit_labels[i])
	
		hlines!([peris[idx[i]], apos[idx[i]]], linestyle=:dot)
		scatter!(-df_peris_apos.t_last_peri[idx[i]] .* lguys.T2GYR, peris[idx[i]], color=COLORS[i])
		scatter!(-df_peris_apos.t_last_apo[idx[i]] .* lguys.T2GYR, apos[idx[i]], color=COLORS[i], marker=:utriangle)

	end

	scatter!([NaN], [NaN], color=:black, label="pericentre")
	scatter!([NaN], [NaN], color=:black, label="apocentre", marker=:utriangle)
	hlines!([NaN], linestyle=:dot, color=:black, label="global peri/apocentre")
	
	Legend(fig[1, 2], ax)
	fig
end

# ╔═╡ ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
positions

# ╔═╡ 130fca42-cee8-4d88-a764-cdded04a636e
lguys.plot_xyz(positions..., labels=orbit_labels)

# ╔═╡ 57a8d1c8-3940-4430-8b46-375fb2bf1695
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "R / kpc",
		ylabel = "z / kpc",
		aspect = DataAspect(),
		limits = (0, nothing, nothing, nothing)
	)
	for i in eachindex(positions)
		x = positions[i][1, :]
		y = positions[i][2, :]
		z = positions[i][3, :]
		R = @. sqrt(x^2 + y^2)
	
		plot!(R, z, 
			
		)
	end

	fig
end

# ╔═╡ 34efe422-e672-4d17-a558-ce32fb704a8e
lguys.plot_xyz(velocities..., units=" / km / s", labels=orbit_labels)

# ╔═╡ ad078920-225d-436e-835b-d87a9db53c49
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "r / kpc",
		ylabel = L"$v$ %$kms_label"
	)

	for i in eachindex(idx)
		scatter!(rs[i], vs[i], color=out.times)
	end

	fig
end

# ╔═╡ 0fb21216-3faf-4a40-9a8e-67ee5f31f933
md"""
The next two plots compare the acceleration and potential along each orbit as compared to a simple static NFW halo (similar to EP20 potential).
"""

# ╔═╡ 09bbae0d-ca3e-426d-b77c-69dd68ca42cc
begin
	M_b = 115
	r_b = 20
	
	a_exp(r) =  M_b / (r+r_b)^2
	phi_exp(r) = -M_b / (r+r_b)
end

# ╔═╡ fa4e7992-1ac6-4d71-a923-8b3cf81d0030
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "time / Gyr",
		ylabel = "relative residual acceleration",
	)
	
	for i in 1:length(idx)
		res = (accs[i] .- a_exp.(rs[i])) ./ a_exp.(rs[i])
		scatter!(-out.times * lguys.T2GYR, res, 
			label = orbit_labels[i]
		)
	end
	fig
end 

# ╔═╡ 35f0ea14-a945-4745-910c-365b730676c5
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "time / Gyr",
		ylabel = "relative residual  potential"
		
	)
	
	for i in 1:length(idx)
		res =  (Φs_ext[i] .- phi_exp.(rs[i])) ./ phi_exp.(rs[i])
		scatter!(-out.times * lguys.T2GYR, res,
			label = orbit_labels[i]
		)
	end
	axislegend()
	
	fig
end

# ╔═╡ 0c90dc59-6641-4094-b283-fcef68271019
md"""
The two plots below show the 3D cloud of initial positions and velocities for visualization purposes.
"""

# ╔═╡ f6b27164-ee7c-428b-aefb-75e89d178f3e
let
	fig = Figure()
	ax = Axis3(fig[1, 1], 
		xlabel = L"$x$ / kpc",
		ylabel = L"$y$ / kpc",
		zlabel = L"$z$ / kpc",
	)
	
	scatter!(snap.positions,
)
	fig
end

# ╔═╡ 5fdd8307-d528-4cd7-a5e4-1f15aba75cd5
let
	fig = Figure()
	ax = Axis3(fig[1, 1], 
		xlabel = L"$v_x$ %$kms_label",
		ylabel = L"$v_y$ %$kms_label",
		zlabel = L"$v_z$ %$kms_label",
	)
	
	scatter!(-snap.velocities .* lguys.V2KMS,
)
	fig
end

# ╔═╡ c4a1a691-51bb-4c3e-89ab-398841b1d155
md"""
# Orbit Info
"""

# ╔═╡ 519a88f0-8e2d-4c09-83e0-3cc2ee147e35
function get_initial_t(j)
	i = idx[j]
	N = length(rs[j])
	t_ini = -1
	for t in N-2:-1:1
		if diff(rs[j])[t] >= 0 && diff(rs[j])[t+1] <= 0
			t_ini = t
			break
		end
	end
	if t_ini == -1
		t_ini = length(rs[j])
	end
	return t_ini
end
	

# ╔═╡ de1e5245-0946-47cd-8e2c-ba1914cfeb74
begin 
	# orbit info
	for i in 1:length(idx)
		t = get_initial_t(i)
		@printf "orbit: \t\t %i\n" i
		@printf "index: \t\t %i\n" idx[i]
		
		@printf "pericentre:\t %0.1f\n" peris[idx[i]]
		@printf "apocentre: \t %0.1f\n" apos[idx[i]]

		@printf "time of first apocentre: %f \n" out.times[end] - out.times[t]
		@printf "radius of first apocentre: %f\n" rs[i][t]
		@printf "intial position: [%f, %f, %f]\n" positions[i][:, t]...
		@printf "intial velocity: [%f, %f, %f]\n" -1* velocities[i][:, t]...
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
# ╟─2b9d49c6-74cc-4cce-b29e-04e94776863f
# ╠═6ca3fe17-3f13-43fe-967b-881078135ead
# ╠═46348ecb-ee07-4b6a-af03-fc4f2635f57b
# ╠═b9f469ed-6e4e-41ee-ac75-1b5bfa0a114a
# ╟─7edf0c89-cc4e-4dc2-b339-b95ad173d7e7
# ╠═e9e2c787-4e0e-4169-a4a3-401fea21baba
# ╠═d975d00c-fd69-4dd0-90d4-c4cbe73d9754
# ╠═3b83205d-91c1-481e-9305-0d59bc692135
# ╠═8b818798-69fb-481d-ade1-9fd436b1f281
# ╟─8f70add4-effe-437d-a10a-4e15228f9fec
# ╠═b15fb3ac-2219-419a-854a-31a783acf891
# ╠═9583b7d0-0d86-4346-998b-000ea68e94b6
# ╠═fa790a4d-e74f-479b-8ff6-aa2f23cb573d
# ╠═eff58c52-a32b-4faa-9b98-c8234d9b21fc
# ╟─88536e86-cf2a-4dff-ae64-514821957d40
# ╟─26d616da-95ec-4fb9-b9a8-2f095d74c722
# ╠═9a22d47b-8474-4596-b418-de33eb07c627
# ╠═da6e5566-f2df-4feb-9188-53eca9a1a0d5
# ╠═384be6a6-f9d9-47e0-9792-aef6689dcbdb
# ╠═4481a5e1-d635-4fea-a5c5-c85f0d6df62f
# ╠═392315ee-72d2-4a14-9afc-5fd6424b3e83
# ╟─950e0210-e4fe-4bad-b82f-69247fd0edd8
# ╠═4a4b8b73-92c3-4e1f-93d8-e44369b8f148
# ╠═e5f728b8-8412-4f57-ad38-a0a35bb08a48
# ╠═413d4e5d-c9cd-4aca-be1e-d132b2bd616d
# ╠═17a63cc8-84f4-4248-a7b0-c8378454b1f7
# ╠═bbf3f229-fc3c-46ae-af28-0f8bd81e7d32
# ╟─d54733e9-7552-4d9c-b8e8-670433469385
# ╠═de2f3380-90df-48f5-ba60-8417e91f4818
# ╠═c8aec4f8-975f-4bbc-b874-bf0172d35868
# ╠═4ee33ce2-c00a-4fcf-b7fc-b78c1677c9e4
# ╠═34104429-05e0-40a6-83e5-078dbe346504
# ╠═e5825c4a-b446-44a3-8fd5-d94664965aca
# ╠═ef57611c-2986-4b03-aa5a-ab45003edd72
# ╟─1acef60e-60d6-47ba-85fd-f9780934788b
# ╟─50baf5a6-fb5b-494e-95f3-53414a9f1cc0
# ╠═ca1c236e-795a-408b-845b-9c13bc838619
# ╠═46b4242b-8af7-4233-8ecf-d86740b4c884
# ╠═5a903509-e2cc-4bac-9c59-d1689ccc408e
# ╠═92aac8e8-d599-4a1e-965e-e653bc54c509
# ╠═471a5501-bc7b-4114-be1a-9088b4e674cd
# ╠═68b0383a-3d5a-4b94-967c-f0e31e8a0ce1
# ╠═5f11f6ab-c9ab-4600-acca-a0bb84d81a12
# ╠═44660b2f-6220-473b-bb2f-07e23b176491
# ╠═d3063e30-2cb3-4f1b-8546-0d5e81d90d9f
# ╟─c48b4e73-480e-4a50-b5fc-db5f6c5b040e
# ╠═43d43f63-4c13-4b23-950e-ada59aa86bc9
# ╠═8b6f95f7-284f-4133-b6f1-a22dd9c405f0
# ╟─69e77193-29cc-4304-98a1-44828eaedf9f
# ╟─89b81ed0-82a1-4a81-a7fd-b6be0644a79d
# ╠═ede3836c-740d-4ac7-bbc7-3165981a1878
# ╟─6b95d3b2-38db-4376-83b5-8c6e6f1fdfa2
# ╟─ac81acd8-4a78-4230-bc70-3b78a861b618
# ╟─16f4ac20-d8cf-4218-8c01-c15e04e567fb
# ╠═d31f91e8-6db6-4771-9544-8e54a816ecc1
# ╠═5be3fdaa-5c87-4fef-b0eb-06dfa780cb11
# ╟─14c36202-66ca-46b3-b282-3895b72311fe
# ╟─e5d40e2f-ac47-4827-853d-2f94bc39a624
# ╠═ee01b25e-c32e-4f6e-96d6-cb9c6f3ea95c
# ╠═130fca42-cee8-4d88-a764-cdded04a636e
# ╟─57a8d1c8-3940-4430-8b46-375fb2bf1695
# ╠═34efe422-e672-4d17-a558-ce32fb704a8e
# ╠═ad078920-225d-436e-835b-d87a9db53c49
# ╟─0fb21216-3faf-4a40-9a8e-67ee5f31f933
# ╠═09bbae0d-ca3e-426d-b77c-69dd68ca42cc
# ╠═fa4e7992-1ac6-4d71-a923-8b3cf81d0030
# ╟─35f0ea14-a945-4745-910c-365b730676c5
# ╟─0c90dc59-6641-4094-b283-fcef68271019
# ╟─f6b27164-ee7c-428b-aefb-75e89d178f3e
# ╟─5fdd8307-d528-4cd7-a5e4-1f15aba75cd5
# ╟─c4a1a691-51bb-4c3e-89ab-398841b1d155
# ╠═519a88f0-8e2d-4c09-83e0-3cc2ee147e35
# ╠═de1e5245-0946-47cd-8e2c-ba1914cfeb74
