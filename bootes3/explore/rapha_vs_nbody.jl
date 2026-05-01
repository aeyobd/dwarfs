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

# ╔═╡ 3ee45a93-2b84-47f1-8502-7fb58ea2a8be
modelname = let
	# "bootes3/1e6_v30_r3.0/5_peri_18kpc"
	# "sculptor/1e7_new_v31_r3.2/orbit_smallperi"
	# "ursa_minor/1e7_new_v38_r4.0/orbit_smallperi.5"
	# "ursa_minor/1e6_new_v38_r4.0/orbit_smallperi.5"
	"isothermal/1e5_v30_r3.0/orbit_15_150"
end

# ╔═╡ ff87ec65-2f65-48ad-9a51-32731d9cbeb7
starsname = "exp2d_rs0.20"

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

# ╔═╡ d843b373-3afc-4aae-8358-1c718a43d885
modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", modelname)

# ╔═╡ a829e274-961f-456d-89e6-fc5f1014af1e
pot = Agama.Potential(file=joinpath(modeldir, "simulation/agama_potential.ini"))

# ╔═╡ 01b9bd98-ff58-4658-99f8-c49efe41337f
orbit_props = TOML.parsefile(joinpath(modeldir, "orbital_properties.toml"))

# ╔═╡ 9fade8be-9249-4eb0-b12f-b0d65ff9d1f2
orbit_ideal = TOML.parsefile(joinpath(modeldir, "simulation/orbit.toml"))

# ╔═╡ 4dc7f54b-af2f-4a29-a85d-ca4c96d4baec
peri = orbit_props["pericentre"]

# ╔═╡ c6b167c0-b784-4368-a356-c64501fc26be
apo = let
	a = orbit_props["apocentre"][1]
	if isnan(a)
		a = orbit_ideal["apocentre"]
	end

	a
end

# ╔═╡ fd18ca6c-c652-4e7d-94fe-9239fe9f1ae2
period = orbit_props["period"] / T2GYR

# ╔═╡ 4ca685d1-e62c-49c3-b316-46ee661d3086
orbit_props["t_f_gyr"] / orbit_props["period"], length(orbit_props["apocentre"])

# ╔═╡ f42a6a93-bf30-4f82-bf37-2adbb20c0ad3
V0 = Agama.circular_velocity(pot, peri)

# ╔═╡ c37ea901-cd5a-4800-99da-4c0bf6e9cbe8
model_track = read_fits(joinpath(modeldir, "profiles_scalars.fits"))

# ╔═╡ 913c7814-8b18-433c-8db9-b10f452874cc
stars_track = read_fits(joinpath(modeldir, "stars",starsname,  "stellar_profiles_3d_scalars.fits"))

# ╔═╡ b4cf0a16-0d0f-4988-9314-9e0e02dc3dbc
rmax_i = model_track.r_circ_max[1]

# ╔═╡ 33868ab2-ea24-42d0-936f-7e452898e106
vmax_i = model_track.v_circ_max[1]

# ╔═╡ 929b599c-38a3-48ad-9a24-41e85e35ff3e
md"""
# Tidal prediction
"""

# ╔═╡ 6b617252-7929-4421-bcad-d05c64bfc388
T_peri = 2π * peri / V0

# ╔═╡ 23ef0e79-15cd-4c23-80c5-8a10cd2341a2
Φ(r) = Agama.potential(pot, [0, 0, r])

# ╔═╡ 7d071e1f-9d49-48a2-b568-5772b662b8ad
function peri_apo_to_E_L(peri, apo)
	L = sqrt(
		(Φ(peri) - Φ(apo)) / 
			(1/2 * (1/apo^2 - 1/peri^2))
	)

	E = Φ(peri) + L^2 / 2peri^2

	@info Φ(apo) + L^2 / 2apo^2, E
	@assert Φ(apo) + L^2 / 2apo^2 ≈ E
	return E, L
end

# ╔═╡ a9bbb390-fbc9-475b-bbc1-1d557ad9c232
function orbital_period(peri, apo)
	E, L = peri_apo_to_E_L(peri, apo)
	return 2*LilGuys.integrate(
		r -> 1 / sqrt(2*(E - Φ(r)) - L^2 / r^2),
		peri * (1 + 1e-12), apo*(1 - 1e-12)
	)
end

# ╔═╡ 6702b75a-dfba-4894-8729-10af81b0034c
function make_track(; rmax_i=rmax_i, vmax_i=vmax_i, peri=peri, apo=apo, t_rel_max=nothing, V0=V0)

	f_ecc = Rapha.ecc_factor(peri, apo)
	period = orbital_period(peri, apo)
	if isnothing(t_rel_max)
		t_rel_max = 10 / T2GYR / (f_ecc * period)
	end

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
function make_orbit_track(; period=period,  rmax_i=rmax_i, vmax_i=vmax_i, peri=peri, apo=apo, t_rel_max=20, V0=V0, n_orbits=10 / T2GYR / period)

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
		time = orbital_period(peri, apo) * t_orb,
		M_max = vmax .^2 .* rmax, 
	)
end

# ╔═╡ 24b98381-2d50-405e-97ca-d041d12625db
orbital_period(peri, apo)

# ╔═╡ c76e56e8-1425-4aaf-bb10-fafd70a148b4
T_orb = period

# ╔═╡ 1733ed18-3466-458f-bfd2-3cf34ec69540
E_i, L_i = peri_apo_to_E_L(peri, apo)

# ╔═╡ f397349c-b23e-42d0-a7d1-933f25f3e4c6
orbit = LilGuys.agama_orbit(pot, Galactocentric([0, 0, apo], [0, L_i/apo * V2KMS, 0]), timerange=(0, 10*T_orb))

# ╔═╡ 520ae3b5-5a01-425e-887b-3b210954b123
let
	fig = Figure()
	ax = Axis(fig[1,1])

	lines!(LilGuys.times(orbit) / T_orb, radii(orbit))

	fig
end

# ╔═╡ 99df26cf-1707-4c9a-815f-d367652fa746
E0 = Φ(apo) + 1/2 * (L_i/apo)^2

# ╔═╡ f2c04684-d489-46ac-a2e7-deafcffecfd5
md"""
# Reproduces Fig. 4
"""

# ╔═╡ 51dfba63-93ae-47bb-b63b-7b521f388b13
df_orbits = make_orbit_track()

# ╔═╡ 221ce0a0-498c-477f-8fe7-66950b5aca4f
hcat(log10.(df_orbits.r_max), log10.(df_orbits.v_max .* V2KMS), log10.(df_orbits.t_max * T2GYR))

# ╔═╡ b59b921f-248a-43a2-b2df-b1a516287277
2π * df_orbits.r_max ./ df_orbits.v_max ./ df_orbits.t_max

# ╔═╡ 2a9ee177-9dc6-4432-a5c4-b882efb9ca06
T_peri

# ╔═╡ 0d8a0334-7fdf-42de-b87b-f8c49ecf3bfa
T_orb

# ╔═╡ e8f798f5-9cfa-494d-a047-45084779f65b
df_orbits.t_max

# ╔═╡ 4e314385-3b54-4185-a7ee-d4f2ff133cd1
df_orbits.t_max ./ T_peri

# ╔═╡ 35303b9c-3217-4022-946f-3238aad1aa54
println("python tipy_script.py $rmax_i $(vmax_i * V2KMS) $peri $apo $(V0 * V2KMS)")

# ╔═╡ 4df8f988-d160-4727-94f3-e7be08333648
let
	# this code reproduces the output from tipy
	df_orbits = make_orbit_track(n_orbits = 15)
	time = df_orbits.time * T2GYR
	r_max = df_orbits.r_max
	v_max = df_orbits.v_max * V2KMS
	M_max = df_orbits.M_max * M2MSUN
	t_max = df_orbits.t_max * T2GYR

	println("time tmax rmax vmax M_max")
	for i in eachindex(time)
		println("$i $(t_max[i]) $(r_max[i]) $(v_max[i]) $(M_max[i])")
	end
end

# ╔═╡ 7d9b6f14-c2e3-4d73-869d-43f94f6cfe94
let
	fig = Figure(size=(6, 5) .* 72)
	ax = Axis(fig[1,1],
			 xlabel = "orbit time / Gyr",
			 ylabel = "log tmax / Gyr")

	df = make_track()
	df_orbits = make_orbit_track()

	filt = df.v_max .> 0


	scatter!(model_track.time .* T2GYR, log10.(2π * model_track.r_circ_max ./ model_track.v_circ_max .* T2GYR), color=COLORS[3], markersize=3)
	

	lines!(df.time[filt] * T2GYR, log10.(df.t_max[filt] * T2GYR ), color=:black)
	scatter!(df_orbits.time .* T2GYR, log10.(df_orbits.t_max .* T2GYR), color=:black)


	ax_r = Axis(fig[1,2],
			 xlabel = "log rmax",
			 ylabel = "log tmax")



	scatter!(log10.(model_track.r_circ_max), log10.(2π * model_track.r_circ_max ./ model_track.v_circ_max  .* T2GYR), color=COLORS[3], markersize=3)
	
	lines!(log10.(df.r_max[filt]), log10.(df.t_max[filt] * T2GYR ), color=:black)
	scatter!(log10.(df_orbits.r_max), log10.(df_orbits.t_max * T2GYR ), color=:black)

	linkyaxes!(ax, ax_r)



	ax_tv = Axis(fig[2,1],
			 xlabel = "orbit time / Gyr",
			 ylabel = "log vmax")

	scatter!(model_track.time .* T2GYR, log10.(model_track.v_circ_max .* V2KMS), color=COLORS[3], markersize=3)

	
	lines!(df.time[filt] * T2GYR, log10.(df.v_max[filt] .* V2KMS ), color=:black)

	
	scatter!(df_orbits.time .* T2GYR, log10.(df_orbits.v_max .* V2KMS), color=:black)



	ax_rv = Axis(fig[2,2],
			 xlabel = "log rmax",
			 ylabel = "log vmax")

	scatter!(log10.(model_track.r_circ_max), log10.(model_track.v_circ_max .* V2KMS), color=COLORS[3], markersize=3)


	lines!(log10.(df.r_max[filt]), log10.(df.v_max[filt] * V2KMS ), color=:black)
	scatter!(log10.(df_orbits.r_max), log10.(df_orbits.v_max * V2KMS ), color=:black)


	linkyaxes!(ax_tv, ax_rv)

	fig
end

# ╔═╡ 517b7035-2d8e-43c4-a03a-9fccbd7f2d84
df = let
	df = make_track()
	df[df.v_max .> model_track.v_circ_max[end], :]
end

# ╔═╡ dd1f6c90-f63c-4ac0-bea6-eecb20a8f1ff
@assert df.t_max ≈ df.r_max ./ df.v_max * 2π

# ╔═╡ ba0ff3f2-e60e-4709-a0e6-7b93e66c00a7
md"""
# Stellar evolution
"""

# ╔═╡ 2923913e-f8fa-4516-ae00-7255bcc3c6af
halo_i = NFW(r_circ_max=rmax_i, v_circ_max=vmax_i)

# ╔═╡ 0599fecd-485b-4d5e-acfd-69f7e2423d73
prof_stars = LilGuys.load_profile(joinpath(modeldir, "../stars", starsname, "profile.toml"))

# ╔═╡ a6be20d6-814e-43c6-aac8-be98bbc99553
stars_props_i = Rapha.find_initial_df(halo_i, prof_stars, N_r=100)

# ╔═╡ f4ea1949-50df-4386-b4cf-5de425061355
ρ_f = Rapha.final_ρ_star(stars_props_i, 0.3)

# ╔═╡ 59f0da50-29df-4b06-928e-a9b2e56ef095
idx_star = 1:10:size(df, 1)

# ╔═╡ 6f2723e7-0c00-4b16-ad9e-40265c66f142
# ρ_stars = [Rapha.final_ρ_star(stars_props_i, row.M_max / df.M_max[1]) for row in eachrow(df[idx_star, :])]

# ╔═╡ deef3fe4-b04e-4c62-a1c5-6d3eae82943c
function calc_r_h(halo, E_s, α_star=3, β_star = 3; r_range=(1e-4, 1e4), N_r=100)
	r_dm = logrange(r_range[1], r_range[2], N_r) |> collect
	Ψ_i(r) = -LilGuys.potential(halo_i, r)
    ϵ_i = Ψ_i.(r_dm)

    ϵ_max_i = Ψ_i(r_range[1] * 0.5) # 
    x_i = 1 .- ϵ_i / ϵ_max_i


    dNde_star_i = Rapha.dN_dE_poly.(x_i, E_s, α_star, β_star)
    g_i = Rapha.density_of_states.(Ψ_i, ϵ_i)
	
	f_star_i = dNde_star_i ./ g_i

	ρ = Rapha.make_stellar_density(LilGuys.lerp(ϵ_i, f_star_i), Ψ_i, r_dm)
	return Rapha.calc_r_h(ρ.interp, R_end=1e4)
end


# ╔═╡ 48fb401b-e6eb-44fb-9e86-1641aeab5e7c
r_h = [Rapha.calc_r_h(ρ_f.interp, R_end=1e4) for ρ_f in ρ_stars]


# ╔═╡ 99dfc3b6-8bed-4308-9813-3906bb323219
M = [Rapha.calc_mass(ρ_f.interp, 1e4) for ρ_f in ρ_stars]

# ╔═╡ 6593bb93-ef90-4b97-9dc7-147b303c63a8
function final_mass(M_max_rel)
	x_f, dNde_star_f = Rapha.final_stellar_energies(stars_props_i.x_i, stars_props_i.dNde_star_i, M_max_rel)
	return sum(LilGuys.gradient(x_f) .* dNde_star_f)
end

# ╔═╡ 1b1c9405-411c-43ae-b522-ffc94ee3b402
L2 = final_mass.(df.M_max / df.M_max[1])[idx_star]

# ╔═╡ 2c4dc40f-f7be-4f6f-87db-a7e0708bd75d
let
	fig = Figure(size=(5, 3) .* 72)
	ax = Axis(fig[1,1],
			 xlabel = "time",
			 ylabel = L"\log r_h")

	scatterlines!((stars_track.time * T2GYR), log10.(stars_track.r_h ), color=:black)

	lines!((df.time[idx_star]*T2GYR), log10.(r_h))

	xlims!(0, 12)
	ylims!(-1.5, 0.5)


	ax = Axis(fig[1,2],
			 xlabel = "log r_h",
			 ylabel = "log stellar mass")

	scatterlines!(log10.(stars_track.r_h), log10.(stars_track.bound_mass ./ stars_track.bound_mass[1]), color=:black)

	lines!(log10.(r_h), log10.(L2 ./ L2[1]))



	fig
	
end

# ╔═╡ 837c1520-99f0-464a-aa6a-a67a8886a4ee
let
	fig = Figure()
	ax = Axis(fig[1,1])

	x = LinRange(-1.5, 1, 1000)

	lines!(x, log10.(LilGuys.v_circ.(halo_i, rmax_i .* 10 .^ x) ./ vmax_i))

	x, y = LilGuys.EN21_tidal_track(1, 1)
	lines!(log10.(x), log10.(y))

	scatterlines!(log10.(stars_track.r_h ./ rmax_i), log10.(stars_track.sigma_v ./ vmax_i .* sqrt(3)), color=:black)

	x = []
	y = []
	log_M = []
	for i in eachindex(ρ_stars)
		M_rel = df.M_max[idx_star[i]] / df.M_max[1]
		h = Rapha.final_halo(rmax_i, vmax_i, M_rel)
		
		prof = ρ_stars[i]
		r = r_h[i]

		push!(x, log10(r / rmax_i))
		push!(y, log10(LilGuys.v_circ(h, r) ./ vmax_i))
		push!(log_M, log10(M_rel))
	end

	lines!(x, y, color=log_M)

	fig
end

# ╔═╡ Cell order:
# ╠═3ee45a93-2b84-47f1-8502-7fb58ea2a8be
# ╠═ff87ec65-2f65-48ad-9a51-32731d9cbeb7
# ╠═156c47a4-39b1-11f1-861a-df30869de806
# ╠═9f80a0a7-d83b-4a27-8038-dc375a82dc15
# ╠═fddd99d4-fd2f-4f90-bf04-8f8d20c39c48
# ╠═cc264eb9-65df-4fc0-9be1-26142ecd2690
# ╠═fdc2fb9c-b900-4ea0-94b0-7dcd9192ad51
# ╠═fcf74f25-3456-4c57-b10f-81885e05ae9b
# ╠═67580761-e143-4824-b829-64d34496675a
# ╟─769ac4b9-cde0-4b18-b4e2-95f790c7fa11
# ╠═d843b373-3afc-4aae-8358-1c718a43d885
# ╠═a829e274-961f-456d-89e6-fc5f1014af1e
# ╠═01b9bd98-ff58-4658-99f8-c49efe41337f
# ╠═9fade8be-9249-4eb0-b12f-b0d65ff9d1f2
# ╠═4dc7f54b-af2f-4a29-a85d-ca4c96d4baec
# ╠═c6b167c0-b784-4368-a356-c64501fc26be
# ╠═fd18ca6c-c652-4e7d-94fe-9239fe9f1ae2
# ╠═4ca685d1-e62c-49c3-b316-46ee661d3086
# ╠═f42a6a93-bf30-4f82-bf37-2adbb20c0ad3
# ╠═c37ea901-cd5a-4800-99da-4c0bf6e9cbe8
# ╠═913c7814-8b18-433c-8db9-b10f452874cc
# ╠═b4cf0a16-0d0f-4988-9314-9e0e02dc3dbc
# ╠═33868ab2-ea24-42d0-936f-7e452898e106
# ╟─929b599c-38a3-48ad-9a24-41e85e35ff3e
# ╠═6702b75a-dfba-4894-8729-10af81b0034c
# ╠═df0e52c1-a8df-4321-ac8d-7ae178b6a7b9
# ╠═6b617252-7929-4421-bcad-d05c64bfc388
# ╠═23ef0e79-15cd-4c23-80c5-8a10cd2341a2
# ╠═a9bbb390-fbc9-475b-bbc1-1d557ad9c232
# ╠═7d071e1f-9d49-48a2-b568-5772b662b8ad
# ╠═24b98381-2d50-405e-97ca-d041d12625db
# ╠═c76e56e8-1425-4aaf-bb10-fafd70a148b4
# ╠═f397349c-b23e-42d0-a7d1-933f25f3e4c6
# ╠═99df26cf-1707-4c9a-815f-d367652fa746
# ╠═520ae3b5-5a01-425e-887b-3b210954b123
# ╠═1733ed18-3466-458f-bfd2-3cf34ec69540
# ╠═f2c04684-d489-46ac-a2e7-deafcffecfd5
# ╠═51dfba63-93ae-47bb-b63b-7b521f388b13
# ╠═221ce0a0-498c-477f-8fe7-66950b5aca4f
# ╠═b59b921f-248a-43a2-b2df-b1a516287277
# ╠═2a9ee177-9dc6-4432-a5c4-b882efb9ca06
# ╠═0d8a0334-7fdf-42de-b87b-f8c49ecf3bfa
# ╠═e8f798f5-9cfa-494d-a047-45084779f65b
# ╠═4e314385-3b54-4185-a7ee-d4f2ff133cd1
# ╠═35303b9c-3217-4022-946f-3238aad1aa54
# ╠═4df8f988-d160-4727-94f3-e7be08333648
# ╠═7d9b6f14-c2e3-4d73-869d-43f94f6cfe94
# ╠═517b7035-2d8e-43c4-a03a-9fccbd7f2d84
# ╠═dd1f6c90-f63c-4ac0-bea6-eecb20a8f1ff
# ╟─ba0ff3f2-e60e-4709-a0e6-7b93e66c00a7
# ╠═2923913e-f8fa-4516-ae00-7255bcc3c6af
# ╠═0599fecd-485b-4d5e-acfd-69f7e2423d73
# ╠═a6be20d6-814e-43c6-aac8-be98bbc99553
# ╠═f4ea1949-50df-4386-b4cf-5de425061355
# ╠═59f0da50-29df-4b06-928e-a9b2e56ef095
# ╠═6f2723e7-0c00-4b16-ad9e-40265c66f142
# ╠═deef3fe4-b04e-4c62-a1c5-6d3eae82943c
# ╠═48fb401b-e6eb-44fb-9e86-1641aeab5e7c
# ╠═99dfc3b6-8bed-4308-9813-3906bb323219
# ╠═2c4dc40f-f7be-4f6f-87db-a7e0708bd75d
# ╠═6593bb93-ef90-4b97-9dc7-147b303c63a8
# ╠═1b1c9405-411c-43ae-b522-ffc94ee3b402
# ╠═837c1520-99f0-464a-aa6a-a67a8886a4ee
