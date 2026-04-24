### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 1fee9c80-2157-11f1-9d09-f949bf60fb6c
begin
	import Pkg; Pkg.activate()

	using LilGuys
	import TOML
	using CairoMakie
	using Arya
end

# ╔═╡ f373ca79-fcba-45fe-851b-5166af77caf8
using Printf

# ╔═╡ 0a983127-d215-4257-a0f3-f01bed02fd8f
using OrderedCollections

# ╔═╡ 85e0ecf3-a3e1-4c66-9bd7-85b72091d30d
include("paper_style.jl")

# ╔═╡ 88b25ecb-99a4-4176-b1b7-c4cec818cc16
import Agama

# ╔═╡ 4f8ee02f-10f4-43c3-aa0d-f279cec5615d
simsdir = ENV["DWARFS_SIMS"]

# ╔═╡ e60a0096-0df9-4cf3-b559-b905324a9c04
anadir = joinpath(ENV["DWARFS_ROOT"], "analysis")

# ╔═╡ bb062773-2ae7-43ae-b522-e11ea0a46d7a
modelnames = [
	"ursa_minor/1e5_new_v38_r4.0/orbit_smallperi.1",
	"ursa_minor/1e5_new_v38_r4.0/orbit_smallperi.2",
	"ursa_minor/1e5_new_v38_r4.0/orbit_smallperi.3",
	"ursa_minor/1e6_new_v38_r4.0/orbit_smallperi.3",
	"ursa_minor/1e6_new_v38_r4.0/orbit_smallperi.4",
	"ursa_minor/1e6_new_v38_r4.0/orbit_smallperi.5",
	"ursa_minor/1e7_new_v38_r4.0/orbit_smallperi.5",
]

# ╔═╡ e7f58e8b-6af1-467f-98d6-c934e2ae257f
model_labels = ["1", "2", "3", "3-mr", "4-mr", "5-mr", "5-hr"]

# ╔═╡ c8739cdd-aac9-4aa5-a8b8-c620ec41302b
model_offsets = 1:length(modelnames)

# ╔═╡ 5450d97c-9d37-48fd-8bed-e5980c3159ef
resolutions = [1, 1, 1, 2, 2, 2, 4]

# ╔═╡ da848ea2-beb0-4c60-8230-ca4ba604827e
md"""
# Utilities
"""

# ╔═╡ 62f4b647-5351-4c8a-8fa1-191313662485
function get_orbit_props(modelname)
	return TOML.parsefile(joinpath(simsdir, modelname, "orbit.toml"))
end

# ╔═╡ 617140da-db2e-418c-b128-d5896b009866
function get_initial_posvel(modelname)
	props = get_orbit_props(modelname)

	return Galactocentric(props["x_i"], props["y_i"], props["z_i"], 
						 props["v_x_i"], props["v_y_i"], props["v_z_i"])
end

# ╔═╡ 71cd9e1c-f498-48dc-8ee1-a6a7a4c4d4a1
function get_target_coord(modelname)
	return LilGuys.transform(Galactocentric, ICRS(get_orbit_props(modelname)))
end

# ╔═╡ 1c669e57-c252-4ec3-85a3-1d1c91071f00
function get_chi2(modelname)
	props = TOML.parsefile(joinpath(anadir, modelname, "orbital_properties.toml"))

	println(keys(props))
	return props["chi2_best"], props["chi2_best_interp"]
end

# ╔═╡ 994ef5f0-1712-4e79-adac-4caf63a5a31f
get_chi2(modelnames[1])

# ╔═╡ 0d88f27b-6c0b-4277-8aff-bd5681a0691a
function get_idx_f(modelname)
	idx_f = TOML.parsefile(joinpath(anadir, modelname, "orbital_properties.toml"))["idx_f"]
	if modelname == "ursa_minor/1e5_new_v38_r4.0/orbit_smallperi.3"
		idx_f += 1
	end

	idx_f
end


# ╔═╡ 3339aef6-09e5-4a81-81f7-ff7de3fedf7c
function get_final_posvel(modelname)
	idx_f = get_idx_f(modelname)
	orbit = Orbit(joinpath(anadir, modelname, "centres.hdf5"))

	return Galactocentric(orbit.positions[:, idx_f], orbit.velocities[:, idx_f] * V2KMS)
end

# ╔═╡ bab03f4a-e1eb-4f47-9030-df6ec6779fda
function get_t_f(modelname)
	idx_f = get_idx_f(modelname)
	orbit = Orbit(joinpath(anadir, modelname, "centres.hdf5"))

	return orbit.times[idx_f]
end


# ╔═╡ 4399c8f2-e1c6-4f4c-aa6d-605c3417f4c8
function get_action_angle_change(modelname)
	props = get_orbit_props(modelname)

	if "dJ_R" ∉ keys(props)
		return [NaN, NaN, NaN], [NaN, NaN, NaN]
	end
	actions = [props["dJ_R"], props["dJ_z"], props["dJ_phi"]]
	angles = [props["dtheta_R"], props["dtheta_z"], props["dtheta_phi"]]

	return actions, angles
end

# ╔═╡ 7f6d22fc-7234-48ee-8836-be4bfbaf7269
function get_action_angle_ini(modelname)
	props = get_orbit_props(modelname)

	actions = [props["J_R"], props["J_z"], props["J_phi"]]
	angles = [props["theta_R"], props["theta_z"], props["theta_phi"]]

	return actions, angles
end

# ╔═╡ 2e1c5867-6bd6-49df-8ed1-85d034a1f373
FIGDIR = "figures"

# ╔═╡ b4ffad22-96d0-414c-8629-d57c770f7903
md"""
# Setup
"""

# ╔═╡ 3b5c2a7d-3f01-4ec2-a47c-e7b8ab957902
pot = Agama.Potential(file=joinpath(simsdir, modelnames[1], "agama_potential.ini"))

# ╔═╡ 021af6f9-70a1-4197-b00f-d4bb5ac461f2
action_finder = Agama.ActionFinder(pot)

# ╔═╡ 26c75b36-bf51-4b5f-ab09-5a94ad00e410
function to_action_angles(coord)
	actions, angles, frequencies =  Agama.actions_angles(action_finder, LilGuys.position(coord), LilGuys.velocity(coord) / V2KMS)

	return actions * V2KMS, angles
end

# ╔═╡ 768868f9-4eb3-4332-bb90-963d09576da9
initial_coords = get_initial_posvel.(modelnames)

# ╔═╡ 168dbab2-6202-439d-b846-fc9ba06ab352
final_coords = get_final_posvel.(modelnames)

# ╔═╡ 78e0ecb9-6b37-4237-a25c-1add3045199c
initial_actions_angles = to_action_angles.(initial_coords)

# ╔═╡ a770f174-3a13-4116-94d4-eddb6391f50e
initial_actions_2 = get_action_angle_ini.(modelnames[2:end])

# ╔═╡ 0ec1e13f-2914-4587-9d91-b34bf17df5a5
d_action_angles = [(aa2[1] .- aa1[1], aa2[2] .- aa1[2]) for (aa1, aa2) in zip(initial_actions_angles[1:end-1], initial_actions_angles[2:end])]

# ╔═╡ e9f2d926-0e9e-4fa6-9ff7-7b754002917e
[73.6872
432.405
-143.106] / V2KMS

# ╔═╡ e0b0b01f-3bec-4328-a919-0cc006116021
d_action_angles_2 = get_action_angle_change.(modelnames[2:end])

# ╔═╡ f319be76-fd4e-4fec-8cbe-b364add47a14
md"""
# Mean action averaging
"""

# ╔═╡ d8fa1e8c-b4a1-480e-a842-d7ae80e4bd31
function get_actions(af, pos, vel, pos_err, vel_err; units=units, N=1000)

	N = size(pos, 2)
	actions = Matrix{Float64}(undef, 3, N)
	angles = Matrix{Float64}(undef, 3, N)

	actions_err = Matrix{Float64}(undef, 3, N)
	angles_err = Matrix{Float64}(undef, 3, N)

	for i in 1:size(pos, 2)
		pos_samples = pos[:, i] .+ pos_err[i] * randn(3, N)
		vel_samples = vel[:, i] .+ vel_err[i] * randn(3, N)
		act, ang, freq_py = Agama.actions_angles(af, pos_samples, vel_samples, units)
		
		act_m = LilGuys.mean.(eachrow(act))
		act_err = LilGuys.std.(eachrow(act))
		ang_m = LilGuys.mean.(eachrow(ang))
		ang_err = LilGuys.std.(eachrow(ang))
		actions[:, i] = act_m
		actions_err[:, i] = act_err
		angles[:, i] = ang_m
		angles_err[:, i] = ang_err

	end

	return actions, angles, actions_err, angles_err
end

# ╔═╡ 4390a009-1238-4515-9faf-61923c827798
function get_actions(pos, vel)

	N = size(pos, 2)
	actions = Matrix{Float64}(undef, 3, N)
	angles = Matrix{Float64}(undef, 3, N)

	actions, angles = Agama.actions_angles(action_finder, pos, vel)[1:2]
		


	return actions, angles
end

# ╔═╡ 384e7dea-dd8f-4696-9bde-172dd957332d
function get_all_actions(modelname)
	local x_cen, v_cen, times
	idx_f = get_idx_f(modelname)

	orbit = Orbit(joinpath(anadir, modelname, "centres.hdf5"))
	window = 10

	act_nbody, ang_nbody = get_actions(orbit.positions, orbit.velocities)

	return act_nbody[:, 1:idx_f], ang_nbody[:, 1:idx_f], orbit.times[1:idx_f]
end

# ╔═╡ a51612d1-8d83-4c27-9275-35f29fd7f0c7
function get_final_actions(modelname; window=10)
	act_nbody, ang_nbody, _ = get_all_actions(modelname)
	acts = [LilGuys.mean((act_nbody)[i, end-window:end]) for i in 1:3]

	angs = [ang_nbody[i, end] for i in 1:3]

	return V2KMS*acts, angs
end

# ╔═╡ 25d6000f-b6ce-455c-bbd3-79e5e4558d60
md"""
# Plotting
"""

# ╔═╡ 90288fdc-ad5c-4384-ae4c-5accf20cb12b
target_coord = get_target_coord(modelnames[5])

# ╔═╡ 0030b7ce-a59f-479b-b109-d6e6da3c34c1
final_action_angles = get_final_actions.(modelnames)

# ╔═╡ ff9e7006-bf70-4b83-a34c-f8f8bc2d2e0e
final_action_angles[2]

# ╔═╡ e750f858-9bec-4b1c-b6d8-0cbf9884d6de
target_actions = to_action_angles(target_coord)

# ╔═╡ fe735338-9d73-4737-9df3-99182ba71dcc
Niter = length(modelnames)

# ╔═╡ 940e61d2-5ffe-4611-8693-920459ec1c45
all_actions_angles = get_all_actions.(modelnames)

# ╔═╡ b9cc963c-9c1a-478c-af52-54584201ec85
all_actions_angles[1]

# ╔═╡ 59e455f2-4bab-433f-9470-1f62456ce9e5
frequencies = Agama.actions_angles(action_finder, LilGuys.position(target_coord), LilGuys.velocity(target_coord) / V2KMS)[3]

# ╔═╡ 5c141c94-71c8-448b-84c3-2fb971c66e05
9 * frequencies

# ╔═╡ dfba0c75-3c10-42de-b372-3d8cc634c4f7
lw = theme(:linewidth)[]

# ╔═╡ a8f098d0-2f48-4476-a7c6-bc28d6d39e8a
smallfontsize = theme(:fontsize)[]*0.8

# ╔═╡ 06990077-03b7-47f6-a263-9050f2046f46
slw = lw / 2

# ╔═╡ 94d41453-329e-479b-ae34-bf9394cd8166
target_icrs = LilGuys.transform(ICRS, target_coord)

# ╔═╡ c027bb8a-7b57-420a-922f-f2c11e165080
final_icrs = LilGuys.transform.(ICRS, final_coords)

# ╔═╡ 9df6995b-68ce-4621-9f4b-749eb2382bfb
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor", "observed_properties.toml"))

# ╔═╡ 6fa9de41-6efc-49c0-acd9-46829d4841ba
r_h = obs_props["r_h"] / 60

# ╔═╡ f3ca8188-d9a2-4e71-9e07-14aa8b153457
model_colors = COLORS[[1, 2, 3, 3, 4, 5, 5]]

# ╔═╡ c436fbb4-8d37-4e6d-9a8a-392926bc4c90
let
	fig = Figure(
		size = (7, 4) .* 72
	)
	xspan = 0.8
	colors = model_colors
	linewidths = 1.5*lw ./ resolutions

	for j in 1:3
		coord = ["r", "z", "\\phi"][j]
		ax = Axis(fig[j,1],
				  xticks = model_offsets,
				  xtickformat = l -> model_labels,
				  xminorticksvisible=false,
				  ylabel = L"$J_%$(coord)$ / kpc km s$^{-1}$",
				 )
		
		hlines!(ax, target_actions[1][j], color=:black, label="target value", linewidth=lw/2)
		ylims!(target_actions[1][j]-600, target_actions[1][j]+600)

		if j == 1
			x0 = 5
			y0 = target_actions[1][1] - 450
			arrows2d!([x0], [y0], [0.8], [0], shaftwidth=slw, tiplength=3slw, tipwidth=3slw)
			text!(x0, y0, text="time", offset=(0, smallfontsize/4), fontsize=smallfontsize)
			
		end

		ax_2 = Axis(fig[j,2], 
				  xticks = model_offsets,
				  xtickformat = l -> model_labels,
					xminorticksvisible=false,
				  ylabel = L"($\theta_%$(coord) - \Omega_%$coord t$ )/ rad",
				   )
		hlines!(ax_2, target_actions[2][j], color=:black, label="target value", linewidth=lw/2)
		ylims!(target_actions[2][j]-1.5, target_actions[2][j]+1.5)

		if j == 3
			ax.xlabel = "simulation"
			ax_2.xlabel="simulation"
		end
		

		for i in eachindex(model_labels)
			t0 = model_offsets[i]
			acts, angs, ts = all_actions_angles[i]
	
			x = @. xspan * (ts - ts[1]) / (ts[end] - ts[1]) + t0

			lines!(ax, x, V2KMS * acts[j, :], label="model", color=colors[i], linewidth=linewidths[i])
			if i < length(model_offsets)
				arrows2d!(ax, [model_offsets[i] + xspan], [final_action_angles[i][1][j]], 
						  [0.0], [d_action_angles[i][1][j]], 
						  shaftwidth=1, tiplength=3, tipwidth=3, minshaftlength=0, color=:black)
			end



			lines!(ax_2, x, mod.(angs[j, :] .- frequencies[j] .* (ts .- ts[end]), 2π), color=colors[i], linewidth=linewidths[i])
			if i <  length(model_offsets)
				arrows2d!(ax_2, [model_offsets[i] + xspan], [final_action_angles[i][2][j]], 
						  [0.0], [d_action_angles[i][2][j]], shaftwidth=1, tiplength=3, tipwidth=3, minshaftlength=0, color=:black)
			end

		end
	end

	@savefig "action_corrections"

	fig

end

# ╔═╡ 24985ba2-6c30-4168-a43f-1560f45c09b6


# ╔═╡ 36f4647f-b471-4406-8b7d-e234ecfb59a3
let
	fig = Figure(size=(3.3, 2.7) .* 72)

	local ax
	for j in 1:2
		for i in 1:3
			ax = Axis(fig[j,i],
					 xminorticksvisible=false, limits=((0.5, 5.5), nothing))
			hlines!(target_actions[j][i], color=:black, label="target value")

			y_f = [aa[j][i] for aa in final_action_angles]
			y_i = [aa[j][i] for aa in initial_actions_angles]
			scatter!(model_offsets, y_i, label="initial")
			scatter!(model_offsets, y_f, label="final")
			y = [aa[j][i] for aa in initial_actions_angles[1:4]]
			dy = [aa[j][i] for aa in d_action_angles_2]
			if j == 1
				dy *= V2KMS
			end
			arrows2d!(model_offsets[1:Niter-1], y_f[1:Niter-1], zeros(), dy, shaftwidth=1, tiplength=3, tipwidth=3, color=COLORS[3], minshaftlength=0)
			
			if j == 1
				ax.title = [L"R", L"z", L"\theta"][i]
				ax.xticklabelsvisible = false
				ylims!(target_actions[j][i]-2.5*V2KMS, target_actions[j][i]+2.5*V2KMS)

			elseif j == 2
				ax.xlabel = "simulation"
				ylims!(target_actions[j][i]-1.5, target_actions[j][i]+1.5)

			end
			if (i == 1) && (j == 1)
				ax.ylabel = L"action ($\textrm{kpc}\ \textrm{km}\ \textrm{s}^{-1}$)"
			elseif (i==1) && (j==2)
				ax.ylabel = "action angle"
			end
	
		end
	end

	Legend(fig[3, 2], ax, tellwidth=false, tellheight=true, nbanks=3 )

	fig
end

# ╔═╡ 53c6a042-baae-46fd-958d-66cf377df43f
let
	fig = Figure(size=(3.0, 6) .* 72)
	ms = 3 ./ sqrt.(resolutions) * theme(:markersize)[] / 2

	ax_radec = Axis(fig[1,1],
					autolimitaspect = 1/cosd(target_icrs.dec),
					xlabel="RA / degrees",
					ylabel="Dec / degrees",
				   )
	errorscatter!([target_icrs.ra], [target_icrs.dec], 
				  xerror = [r_h / cosd(target_icrs.dec)], yerror = [r_h],
				  color=:black)

	colors = model_colors
	x = [c.ra for c in final_icrs]
	y = [c.dec for c in final_icrs]
	scatter!(x, y, color=colors, markersize=ms)
	text!(x, y, text=model_labels, color=colors, align=(:left, :center), offset=(5, 0))



	ax_pm = Axis(fig[2,1], autolimitaspect=1,
				xlabel=L"$\mu_{\alpha*}$ / mas yr$^{-1}$",
				ylabel=L"$\mu_\delta$ / mas yr$^{-1}$",)

	errorscatter!([target_icrs.pmra], [target_icrs.pmdec], 
				  xerror = [obs_props["pmra_err"]], yerror = [obs_props["pmdec_err"]],
				  color=:black)

	x = [c.pmra for c in final_icrs]
	y = [c.pmdec for c in final_icrs]
	
	scatter!(x, y, color=colors, markersize=ms)
	text!(x, y, text=model_labels, color=colors, align=(:left, :center), offset=(5, 0))


	ax_rv = Axis(fig[3,1], 
				xlabel="distance / kpc",
				ylabel=L"$\textrm{v}_\textrm{los}$ / km s$^{-1}$",)

	errorscatter!([target_icrs.distance], [target_icrs.radial_velocity], 
				  xerror = [obs_props["distance_err"]], yerror = [obs_props["radial_velocity_err"]],
				  color=:black)

	x = [c.distance for c in final_icrs]
	y = [c.radial_velocity for c in final_icrs]
	
	scatter!(x, y, color=colors, markersize=ms)
	text!(x, y, text=model_labels, color=colors, align=(:left, :center), offset=(5, 0))

	@savefig "observed_action_corrections"
	fig
end

# ╔═╡ 4b129935-7707-4e36-86f0-fe35bc537cdb
md"""
# Action table
"""

# ╔═╡ ac79f899-142e-4ad6-8518-80ead6aa0772
function get_all_props(modelname)
	df = OrderedDict()
	initial_coords = get_initial_posvel(modelname)
	final_coords = get_final_posvel(modelname)
	final_icrs = LilGuys.transform(ICRS, final_coords)

	final_actions, final_angles = get_final_actions(modelname)
	initial_actions, initial_angles = to_action_angles(initial_coords)
	d_actions, d_angles = get_action_angle_change(modelname)

	for sym in [:ra, :dec, :distance, :pmra, :pmdec, :radial_velocity]
		df[sym] = getproperty(final_icrs, sym)
	end

	df[:chi2], df[:chi2_interp] = get_chi2(modelname)

	for i in (1:3)
		coord = ["x", "y", "z"][i]
		df[Symbol("$(coord)_i")] = LilGuys.position(initial_coords)[i]
	end

	for i in (1:3)
		coord = ["x", "y", "z"][i]
		df[Symbol("v_$(coord)_i")] = LilGuys.velocity(initial_coords)[i]
	end


	# for i in (1:3)
	# 	coord = ["x", "y", "z"][i]
	# 	df[Symbol("$(coord)_f")] = LilGuys.position(final_coords)[i]
	# end

	# for i in (1:3)
	# 	coord = ["x", "y", "z"][i]
	# 	df[Symbol("v_$(coord)_f")] = LilGuys.velocity(final_coords)[i]
	# end

	# for i in (1:3)
	# 	coord = ["R", "z", "phi"][i]
	# 	df[Symbol("J_$(coord)_i")] = initial_actions[i]
	# end

	# for i in (1:3)
	# 	coord = ["R", "z", "phi"][i]
	# 	df[Symbol("theta_$(coord)_i")] = initial_angles[i]
	# end


	# for i in (1:3)
	# 	coord = ["R", "z", "phi"][i]
	# 	df[Symbol("J_$(coord)_f")] = final_actions[i]
	# end

	# for i in (1:3)
	# 	coord = ["R", "z", "phi"][i]
	# 	df[Symbol("theta_$(coord)_f")] = final_angles[i]
	# end


	for i in 1:3
		coord = ["R", "z", "phi"][i]
		df[Symbol("dJ_$(coord)_i")] = d_actions[i] * V2KMS
	end

	for i in 1:3
		coord = ["R", "z", "phi"][i]
		df[Symbol("dtheta_$(coord)_i")] = d_angles[i]
	end

	df

end

# ╔═╡ 1bd7c8cb-538c-4a6d-a30f-d3393933a9d9
get_all_props(modelnames[5])

# ╔═╡ 4b9e72f2-6f1c-4adc-ab86-08163d4fd93c
modelnames[5]

# ╔═╡ f73f7c85-2a13-4e03-9439-08056c9a8b06
df_out = get_all_props.(modelnames)

# ╔═╡ 7604f735-4f00-4ab7-9210-29d229933b03
df_target = let 
	df = OrderedDict()
	target_icrs = LilGuys.transform(ICRS, target_coord)
	for sym in [:ra, :dec, :distance, :pmra, :pmdec, :radial_velocity]
		df[sym] = getproperty(target_icrs, sym)
	end

	for i in (1:3)
		coord = ["R", "z", "phi"][i]
		df[Symbol("J_$(coord)_i")] = NaN
		df[Symbol("dJ_$(coord)_i")] = NaN

		df[Symbol("theta_$(coord)_i")] = NaN
		df[Symbol("dtheta_$(coord)_i")] = NaN
	end


	for i in (1:3)
		coord = ["R", "z", "phi"][i]
		df[Symbol("J_$(coord)_f")] = target_actions[1][i]
	end

	for i in (1:3)
		coord = ["R", "z", "phi"][i]
		df[Symbol("theta_$(coord)_f")] = target_actions[2][i]
	end

	df
end

# ╔═╡ 7fd40021-8689-464a-9c85-c47105e00a84
keys(df_out[1])

# ╔═╡ 212dab6b-2afd-4b6b-b75e-2a42e7b273f7
function print_row(name, args...)
	@printf "%32s" name
	for arg in args
		s = @sprintf "\$%0.2f\$" arg

		@printf "  & %12s" s
	end

end

# ╔═╡ ec2e14ed-c42e-4869-a82c-aa968f417d27
nice_names = OrderedDict(
	# :ra => raw"$\alpha\ / \ ^\circ$",
	# :dec => raw"$\delta\ / \ ^\circ$",
	# :distance => raw"distance / $\kpc$",
	# :pmra => raw"$\mu_{\alpha*}\ / \ \masyr$",
	# :pmdec => raw"$\mu_{\delta} \ / \ \masyr$",
	# :radial_velocity => raw"$\V_textrm{los} \ / \ \kms$",
	:dJ_R_i => raw"$\Delta J_r\ / \ \kpc\ \kms$",
	:dJ_z_i => raw"$\Delta J_z\ / \ \kpc\ \kms$",
	:dJ_phi_i => raw"$\Delta J_\phi\ / \ \kpc\ \kms$",
	:dtheta_R_i => raw"$\Delta \theta_r$",
	:dtheta_z_i => raw"$\Delta \theta_z$",
	:dtheta_phi_i => raw"$\Delta \theta_\phi$",
	:x_i => raw"$\Delta x_i\ / \ \kpc$",
	:y_i => raw"$\Delta y_i\ / \ \kpc$",
	:z_i => raw"$\Delta z_i\ / \ \kpc$",
	:v_x_i => raw"$\Delta \V_{x, i} \ / \ \kms$",
	:v_y_i => raw"$\Delta \V_{y, i} \ / \ \kms$",
	:v_z_i => raw"$\Delta \V_{z, i} \ / \ \kms$",
	:chi2 => raw"$\tilde\chi^2$",
	:chi2_interp => raw"$\tilde\chi^2_\textrm{interp}$",

)

# ╔═╡ 2a57ffbb-f521-4039-b73e-ff9716b6d1cc
for (name, label) in nice_names
	print_row(label, [df[name] for df in df_out]...,)
	println("  \\\\")
end

# ╔═╡ Cell order:
# ╠═1fee9c80-2157-11f1-9d09-f949bf60fb6c
# ╠═88b25ecb-99a4-4176-b1b7-c4cec818cc16
# ╠═2e1c5867-6bd6-49df-8ed1-85d034a1f373
# ╠═85e0ecf3-a3e1-4c66-9bd7-85b72091d30d
# ╠═4f8ee02f-10f4-43c3-aa0d-f279cec5615d
# ╠═e60a0096-0df9-4cf3-b559-b905324a9c04
# ╠═bb062773-2ae7-43ae-b522-e11ea0a46d7a
# ╠═e7f58e8b-6af1-467f-98d6-c934e2ae257f
# ╠═c8739cdd-aac9-4aa5-a8b8-c620ec41302b
# ╠═5450d97c-9d37-48fd-8bed-e5980c3159ef
# ╟─da848ea2-beb0-4c60-8230-ca4ba604827e
# ╠═62f4b647-5351-4c8a-8fa1-191313662485
# ╠═617140da-db2e-418c-b128-d5896b009866
# ╠═71cd9e1c-f498-48dc-8ee1-a6a7a4c4d4a1
# ╠═3339aef6-09e5-4a81-81f7-ff7de3fedf7c
# ╠═1c669e57-c252-4ec3-85a3-1d1c91071f00
# ╠═994ef5f0-1712-4e79-adac-4caf63a5a31f
# ╠═0d88f27b-6c0b-4277-8aff-bd5681a0691a
# ╠═bab03f4a-e1eb-4f47-9030-df6ec6779fda
# ╠═26c75b36-bf51-4b5f-ab09-5a94ad00e410
# ╠═4399c8f2-e1c6-4f4c-aa6d-605c3417f4c8
# ╠═7f6d22fc-7234-48ee-8836-be4bfbaf7269
# ╟─b4ffad22-96d0-414c-8629-d57c770f7903
# ╠═3b5c2a7d-3f01-4ec2-a47c-e7b8ab957902
# ╠═021af6f9-70a1-4197-b00f-d4bb5ac461f2
# ╠═768868f9-4eb3-4332-bb90-963d09576da9
# ╠═168dbab2-6202-439d-b846-fc9ba06ab352
# ╠═78e0ecb9-6b37-4237-a25c-1add3045199c
# ╠═a770f174-3a13-4116-94d4-eddb6391f50e
# ╠═0ec1e13f-2914-4587-9d91-b34bf17df5a5
# ╠═e9f2d926-0e9e-4fa6-9ff7-7b754002917e
# ╠═e0b0b01f-3bec-4328-a919-0cc006116021
# ╠═f319be76-fd4e-4fec-8cbe-b364add47a14
# ╠═d8fa1e8c-b4a1-480e-a842-d7ae80e4bd31
# ╠═384e7dea-dd8f-4696-9bde-172dd957332d
# ╠═a51612d1-8d83-4c27-9275-35f29fd7f0c7
# ╠═4390a009-1238-4515-9faf-61923c827798
# ╠═25d6000f-b6ce-455c-bbd3-79e5e4558d60
# ╠═90288fdc-ad5c-4384-ae4c-5accf20cb12b
# ╠═0030b7ce-a59f-479b-b109-d6e6da3c34c1
# ╠═ff9e7006-bf70-4b83-a34c-f8f8bc2d2e0e
# ╠═e750f858-9bec-4b1c-b6d8-0cbf9884d6de
# ╠═fe735338-9d73-4737-9df3-99182ba71dcc
# ╠═940e61d2-5ffe-4611-8693-920459ec1c45
# ╠═b9cc963c-9c1a-478c-af52-54584201ec85
# ╠═59e455f2-4bab-433f-9470-1f62456ce9e5
# ╠═5c141c94-71c8-448b-84c3-2fb971c66e05
# ╠═dfba0c75-3c10-42de-b372-3d8cc634c4f7
# ╠═a8f098d0-2f48-4476-a7c6-bc28d6d39e8a
# ╠═06990077-03b7-47f6-a263-9050f2046f46
# ╠═c436fbb4-8d37-4e6d-9a8a-392926bc4c90
# ╠═94d41453-329e-479b-ae34-bf9394cd8166
# ╠═c027bb8a-7b57-420a-922f-f2c11e165080
# ╠═9df6995b-68ce-4621-9f4b-749eb2382bfb
# ╠═6fa9de41-6efc-49c0-acd9-46829d4841ba
# ╠═f3ca8188-d9a2-4e71-9e07-14aa8b153457
# ╠═24985ba2-6c30-4168-a43f-1560f45c09b6
# ╠═36f4647f-b471-4406-8b7d-e234ecfb59a3
# ╠═53c6a042-baae-46fd-958d-66cf377df43f
# ╟─4b129935-7707-4e36-86f0-fe35bc537cdb
# ╠═f373ca79-fcba-45fe-851b-5166af77caf8
# ╠═0a983127-d215-4257-a0f3-f01bed02fd8f
# ╠═ac79f899-142e-4ad6-8518-80ead6aa0772
# ╠═1bd7c8cb-538c-4a6d-a30f-d3393933a9d9
# ╠═4b9e72f2-6f1c-4adc-ab86-08163d4fd93c
# ╠═f73f7c85-2a13-4e03-9439-08056c9a8b06
# ╠═7604f735-4f00-4ab7-9210-29d229933b03
# ╠═7fd40021-8689-464a-9c85-c47105e00a84
# ╠═212dab6b-2afd-4b6b-b75e-2a42e7b273f7
# ╠═ec2e14ed-c42e-4869-a82c-aa968f417d27
# ╠═2a57ffbb-f521-4039-b73e-ff9716b6d1cc
