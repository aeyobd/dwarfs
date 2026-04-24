### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 156c47a4-39b1-11f1-861a-df30869de806
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using CairoMakie, Arya
end

# ╔═╡ fdc2fb9c-b900-4ea0-94b0-7dcd9192ad51
using PyFITS

# ╔═╡ 9c3b167a-a3fd-4738-b303-743587f90755
using OrderedCollections

# ╔═╡ cd0bd3a3-63dd-4cb7-bda0-05dd615a48f7
using Printf

# ╔═╡ 9f80a0a7-d83b-4a27-8038-dc375a82dc15
CairoMakie.activate!(type=:png)

# ╔═╡ fddd99d4-fd2f-4f90-bf04-8f8d20c39c48
dwarfs_dir = ENV["DWARFS_ROOT"]

# ╔═╡ cc264eb9-65df-4fc0-9be1-26142ecd2690
module Rapha
	include(joinpath(ENV["DWARFS_ROOT"], "utils/rapha_utils.jl"))
end

# ╔═╡ fcf74f25-3456-4c57-b10f-81885e05ae9b
import TOML

# ╔═╡ 67580761-e143-4824-b829-64d34496675a
import Agama

# ╔═╡ 769ac4b9-cde0-4b18-b4e2-95f790c7fa11
md"""
# Model loading
"""

# ╔═╡ e8da258e-dddb-4e74-99c5-b465f90687ff


# ╔═╡ 929b599c-38a3-48ad-9a24-41e85e35ff3e
md"""
# Tidal prediction
"""

# ╔═╡ 6702b75a-dfba-4894-8729-10af81b0034c
function make_track(df)
	rmax_i=df["r_max_i"] 
	vmax_i= df["v_max_i"]
	peri = df["peri"]
	apo = df["apo"]
	V0 = df["V0"]

	f_ecc = Rapha.ecc_factor(peri, apo)
	period = df["T_orb"]
	t_rel_max = 10 / T2GYR / (f_ecc * period)

	t_rel = LinRange(0, t_rel_max, 1000) |> collect

	
	rvmax = Rapha.rapha_final_halo.(rmax_i, vmax_i, peri, apo, t_rel .* f_ecc, V0=V0)


	T_peri = 2π * peri / V0

	
	rmax = first.(rvmax)
	vmax = last.(rvmax)

	t_max = @. rmax * 2π / vmax

	LilGuys.DataFrame(
		r_max = rmax,
		v_max = vmax,
		t_max = t_max,
		t_rel = t_rel,
		time = period * t_rel .* f_ecc,
		M_max = vmax .^2 .* rmax, 
	)
end

# ╔═╡ df0e52c1-a8df-4321-ac8d-7ae178b6a7b9
function make_orbit_track(df)
	rmax_i=df["r_max_i"] 
	vmax_i= df["v_max_i"]
	peri = df["peri"]
	apo = df["apo"]
	V0 = df["V0"]

	f_ecc = Rapha.ecc_factor(peri, apo)
	period = df["T_orb"]
	
	t_rel_max = 10 / T2GYR / (f_ecc * period)


	n_orbits = 5
	
	f_ecc = Rapha.ecc_factor(peri, apo)

	t_orb = 0:round(Int, n_orbits)

	
	rvmax = Rapha.rapha_final_halo.(rmax_i, vmax_i, peri, apo, t_orb, V0=V0)


	T_peri = 2π * peri / V0

	
	rmax = first.(rvmax)
	vmax = last.(rvmax)

	t_max = @. rmax * 2π / vmax

	LilGuys.DataFrame(
		r_max = rmax,
		v_max = vmax,
		t_max = t_max,
		time = period * t_orb,
		M_max = vmax .^2 .* rmax, 
	)
end

# ╔═╡ ae868b06-0294-4560-8881-0e7e5dabef42
make_tracks(props) = OrderedDict(
	label => make_orbit_track(prop) for (label, prop) in props
)

# ╔═╡ e37235f9-3e65-4c55-933a-0d214df1a976
function make_orbit(pot, pos_i, vel_i)
	orbit = LilGuys.agama_orbit(pot, Galactocentric(pos_i, vel_i * V2KMS), timerange=(0, 10 / T2GYR))
	return orbit.pericenter, orbit.apocenter
end

# ╔═╡ 7d071e1f-9d49-48a2-b568-5772b662b8ad
function peri_apo_to_E_L(pot, peri, apo)
	Φ(r) = Agama.potential(pot, [0, 0, r])
	
	L = sqrt(
		(Φ(peri) - Φ(apo)) / 
			(1/2 * (1/apo^2 - 1/peri^2))
	)

	E = Φ(peri) + L^2 / 2peri^2

	@assert Φ(apo) + L^2 / 2apo^2 ≈ E
	return E, L
end

# ╔═╡ a9bbb390-fbc9-475b-bbc1-1d557ad9c232
function orbital_period(pot, peri, apo)
	Φ(r) = Agama.potential(pot, [0, 0, r])
	
	E, L = peri_apo_to_E_L(pot, peri, apo)
	return 2*LilGuys.integrate(
		r -> 1 / sqrt(2*(E - Φ(r)) - L^2 / r^2),
		peri * (1 + 1e-12), apo*(1 - 1e-12)
	)
end

# ╔═╡ fbb3d0d2-3367-4bba-b07f-66e93ceb2f51
function load_model_props(modelname)
	modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", modelname)
	
	orbit_props = TOML.parsefile(joinpath(modeldir, "orbital_properties.toml"))
	orbit_ideal = TOML.parsefile(joinpath(modeldir, "simulation/orbit.toml"))
	

	pot = Agama.Potential(file=joinpath(modeldir, "simulation/agama_potential.ini"))

	peri, apo = make_orbit(pot, orbit_ideal["position_i"], orbit_ideal["velocity_i"])
	T_orb = orbital_period(pot, peri, apo)
	V0 = Agama.circular_velocity(pot, peri)

	model_track = read_fits(joinpath(modeldir, "profiles_scalars.fits"))


	
	@info "peri", peri, orbit_props["pericentre"]
	@info "apo", apo, orbit_props["apocentre"]
	@info "period", T_orb, orbit_props["period"]
	

	return OrderedDict(
		"peri" => peri,
		"apo" => apo,
		"f_ecc" => Rapha.ecc_factor(peri, apo),
		"V0" => V0,
		"T_orb_obs" => orbit_props["period"] / T2GYR,
		"T_orb" => T_orb,
		"T_peri" => 2π * peri / V0,
		"r_max_i" => model_track.r_circ_max[1],
		"v_max_i" => model_track.v_circ_max[1],
		"track_obs" => model_track,
	)
end

# ╔═╡ a1645a65-b25c-456e-ac88-e19681aac396
props_peri = OrderedDict(
	"5" => load_model_props("isothermal/1e5_v30_r3.0/orbit_5_150"),
	"10" => load_model_props("isothermal/1e5_v30_r3.0/orbit_10_150"),
	"15" => load_model_props("isothermal/1e5_v30_r3.0/orbit_15_150"),
	# "17" => load_model_props("isothermal/1e5_v30_r3.0/orbit_17_143"),
	"30" => load_model_props("isothermal/1e5_v30_r3.0/orbit_30_150"),
	"30x300" => load_model_props("isothermal/1e5_v15_r3.0/orbit_30_300"),
)

# ╔═╡ cc460907-c5b3-44c2-a4ca-051520f80ab0
props_peri["15"]["T_orb"]

# ╔═╡ 4ebac403-2b1d-4c5c-bb10-bbfbe497f527
tracks_peri = make_tracks(props_peri)

# ╔═╡ 67c4ae51-4e87-4584-b2f5-293f9abcf02a
props_halo = OrderedDict(
	"v30_r3" => load_model_props("isothermal/1e5_v30_r3.0/orbit_15_150"),
	"v10_r1" => load_model_props("isothermal/1e5_v10_r1.0/orbit_15_150"),
	"v40_r3" => load_model_props("isothermal/1e5_v40_r3.0/orbit_15_150"),
	"v20_r3" => load_model_props("isothermal/1e5_v20_r3.0/orbit_15_150"),
)

# ╔═╡ 381ec7c0-9cc3-44b8-86bb-018b5b3ccdc4
tracks_halo = make_tracks(props_halo)

# ╔═╡ f2c04684-d489-46ac-a2e7-deafcffecfd5
md"""
# Plots
"""

# ╔═╡ 0289dae4-9664-43a0-873d-249967171b27
function plot_vmax_time(models, tracks)
	fig = Figure()
	
	ax = Axis(fig[1,1],
			 xlabel = "orbit time / Torb fecc",
			 ylabel = "log vmax")

	for (label, props) in models
		df = tracks[label]
		model_track = props["track_obs"]

		t_scale = props["T_orb"] * props["f_ecc"]
		scatter!(model_track.time / t_scale, log10.(model_track.v_circ_max .* V2KMS), markersize=3, label=label)

	
		scatterlines!(df.time / t_scale, log10.(df.v_max .* V2KMS ), )
	end

	Legend(fig[1, 2], ax)
	fig

end

# ╔═╡ 6dda351d-abc1-47f1-a541-7da7650a8f59
function plot_vmax_logtime(models, tracks; kwargs...)
	fig = Figure()
	
	ax = Axis(fig[1,1],
			 xlabel = "log orbit time / Torb fecc",
			 ylabel = "log vmax"; 
			 kwargs...)

	for (label, props) in models
		df = tracks[label]
		model_track = props["track_obs"]

		t_scale = props["T_orb"] * props["f_ecc"]
		scatter!(log10.(model_track.time / t_scale), log10.(model_track.v_circ_max .* V2KMS), markersize=3, label=label)

	
		scatterlines!(log10.(df.time / t_scale), log10.(df.v_max .* V2KMS ), )
	end

	Legend(fig[1, 2], ax)
	fig

end

# ╔═╡ 5c408855-843e-4e15-864f-4498c4b07916
function plot_vmax_rmax_t(models, tracks)
	fig = Figure(size=(5, 3) .* 72)
	
	ax = Axis(fig[1,1],
			 xlabel = L"$\log\, {r}_\textrm{max}\ / \ r_\textrm{max, 0}$",
			 ylabel = L"\log\, \textrm{v}_\textrm{max}\ / \ \textrm{v}_\textrm{max, 0}")

	for (label, props) in models
		df = tracks[label]
		model_track = props["track_obs"]
	
		scatter!(
			log10.(model_track.r_circ_max / props["r_max_i"]),
			log10.(model_track.v_circ_max / props["v_max_i"]), 
			markersize=3,
			label = label
		)

	
		scatterlines!(
			log10.(df.r_max / props["r_max_i"]), 
			log10.(df.v_max / props["v_max_i"])
			)
	end

	axislegend(position=:rb)


	ax = Axis(fig[1, 2],
			 xlabel = L"orbit time / $T_\textrm{orb}\,f_\textrm{ecc}$",
			 ylabel = L"$\log\, \textrm{v}_\textrm{max} / $km\,s$^{-1}$")

	for (label, props) in models
		df = tracks[label]
		model_track = props["track_obs"]

		t_scale = props["T_orb"] * props["f_ecc"]
		scatter!(model_track.time / t_scale, log10.(model_track.v_circ_max .* V2KMS), markersize=3, label=label)

	
		scatterlines!(df.time / t_scale, log10.(df.v_max .* V2KMS ), )
	end

	fig

end

# ╔═╡ 0fc9451d-95f6-4ffb-a051-bd629429a96f
function plot_tmax_time(models, tracks)
	fig = Figure()
	
	ax = Axis(fig[1,1],
			 xlabel = "orbit time / Torb fecc",
			 ylabel = "log tmax")

	for (label, props) in models
		df = tracks[label]
		model_track = props["track_obs"]


		t_scale = props["T_orb"] * props["f_ecc"]

		scatter!(
			model_track.time / t_scale, 
			log10.(model_track.r_circ_max * 2π ./ model_track.v_circ_max / props["T_peri"]),
			markersize=3
		)

	
		scatterlines!(df.time / t_scale,
					  log10.(df.t_max  / props["T_peri"]))

		@info "tmax/tperi" df.t_max[1]  / props["T_peri"]
	end

	fig

end

# ╔═╡ 3b0a420d-1588-4a61-8fb4-15cabe941041
md"""
# Effect of pericentre
"""

# ╔═╡ 3317f14f-ca02-4d93-884b-b1ef7ef9ff60
plot_vmax_rmax_t(props_peri, tracks_peri)

# ╔═╡ dcd67e2b-4132-4288-8f48-d1dd495c2a78
let 
	f = plot_vmax_time(props_peri, tracks_peri)

	xlims!(0, 1)
	ylims!(0.5, 1.5)
	f
end

# ╔═╡ 53bb8bc7-2bec-4c29-bca9-8166fd8602da
let
	fig = plot_tmax_time(props_peri, tracks_peri)
	ylims!(-0.5, 0.7)
	xlims!(0, 1.0)
	fig

end

# ╔═╡ 05bb861f-475f-4c71-a8c2-7a79949342dc
plot_vmax_logtime(props_peri, tracks_peri, limits=(-1.5, 0, 0.5, 1.5))
	

# ╔═╡ c7f1efcb-c805-41f3-accf-270e07698e0d
md"""
## Halo ICs
"""

# ╔═╡ 74a7e8e6-a7fd-48ae-9e51-e67ecdab4d67
plot_vmax_rmax_t(props_halo, tracks_halo)

# ╔═╡ da032027-1fd9-48c8-bebe-e5ea340b84e3
plot_vmax_time(props_halo, tracks_halo)

# ╔═╡ 66862b22-f84b-4f25-81ca-a84903ca3c65
let 
	f = plot_tmax_time(props_halo, tracks_halo)
	xlims!(0, 1)
	ylims!(-0.3, 0.4)
	f
end

# ╔═╡ cbce03bc-e2e6-47df-ab69-48b1edec66a7
props_peri["15"]["T_peri"]

# ╔═╡ 2064778d-ddf9-4597-8b7e-e614e45f2390
function print_table_t_tmx_peri(props)
	df = props["track_obs"]
	t_max_rel = df.r_circ_max * 2π ./ df.v_circ_max ./ props["T_peri"]
	time_rel = df.time / props["T_orb"]
	@printf "%12s,%12s\n" "t_torb" "tmx_tperi"
	for i in eachindex(time_rel)
		@printf "%12.8f,%12.8f\n" time_rel[i] t_max_rel[i]
	end
end


# ╔═╡ 87bdc286-c194-4d0d-948d-c2a8a4a8ecec
let
	track = tracks_peri["15"]
	props = props_peri["15"]
		
	t_max_rel = track.t_max ./ props["T_peri"]
	time_rel = track.time / props["T_orb"]
	@printf "%12s,%12s\n" "t_torb" "tmx_tperi"
	for i in eachindex(time_rel)
		@printf "%12.8f,%12.8f\n" time_rel[i] t_max_rel[i]
	end

end

# ╔═╡ cd5c9410-a4e6-4b4b-9bdc-6df3430c17c9
print_table_t_tmx_peri(props_peri["15"])

# ╔═╡ Cell order:
# ╠═156c47a4-39b1-11f1-861a-df30869de806
# ╠═9f80a0a7-d83b-4a27-8038-dc375a82dc15
# ╠═fddd99d4-fd2f-4f90-bf04-8f8d20c39c48
# ╠═cc264eb9-65df-4fc0-9be1-26142ecd2690
# ╠═fdc2fb9c-b900-4ea0-94b0-7dcd9192ad51
# ╠═fcf74f25-3456-4c57-b10f-81885e05ae9b
# ╠═67580761-e143-4824-b829-64d34496675a
# ╠═9c3b167a-a3fd-4738-b303-743587f90755
# ╟─769ac4b9-cde0-4b18-b4e2-95f790c7fa11
# ╠═fbb3d0d2-3367-4bba-b07f-66e93ceb2f51
# ╠═ae868b06-0294-4560-8881-0e7e5dabef42
# ╠═cc460907-c5b3-44c2-a4ca-051520f80ab0
# ╠═a1645a65-b25c-456e-ac88-e19681aac396
# ╠═67c4ae51-4e87-4584-b2f5-293f9abcf02a
# ╠═381ec7c0-9cc3-44b8-86bb-018b5b3ccdc4
# ╠═4ebac403-2b1d-4c5c-bb10-bbfbe497f527
# ╠═e8da258e-dddb-4e74-99c5-b465f90687ff
# ╟─929b599c-38a3-48ad-9a24-41e85e35ff3e
# ╠═6702b75a-dfba-4894-8729-10af81b0034c
# ╠═df0e52c1-a8df-4321-ac8d-7ae178b6a7b9
# ╠═e37235f9-3e65-4c55-933a-0d214df1a976
# ╠═a9bbb390-fbc9-475b-bbc1-1d557ad9c232
# ╠═7d071e1f-9d49-48a2-b568-5772b662b8ad
# ╟─f2c04684-d489-46ac-a2e7-deafcffecfd5
# ╠═0289dae4-9664-43a0-873d-249967171b27
# ╠═6dda351d-abc1-47f1-a541-7da7650a8f59
# ╠═5c408855-843e-4e15-864f-4498c4b07916
# ╠═0fc9451d-95f6-4ffb-a051-bd629429a96f
# ╠═3b0a420d-1588-4a61-8fb4-15cabe941041
# ╠═3317f14f-ca02-4d93-884b-b1ef7ef9ff60
# ╠═dcd67e2b-4132-4288-8f48-d1dd495c2a78
# ╠═53bb8bc7-2bec-4c29-bca9-8166fd8602da
# ╠═05bb861f-475f-4c71-a8c2-7a79949342dc
# ╠═c7f1efcb-c805-41f3-accf-270e07698e0d
# ╠═74a7e8e6-a7fd-48ae-9e51-e67ecdab4d67
# ╠═da032027-1fd9-48c8-bebe-e5ea340b84e3
# ╠═66862b22-f84b-4f25-81ca-a84903ca3c65
# ╠═cd0bd3a3-63dd-4cb7-bda0-05dd615a48f7
# ╠═cbce03bc-e2e6-47df-ab69-48b1edec66a7
# ╠═2064778d-ddf9-4597-8b7e-e614e45f2390
# ╠═87bdc286-c194-4d0d-948d-c2a8a4a8ecec
# ╠═cd5c9410-a4e6-4b4b-9bdc-6df3430c17c9
