### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ bb4f68d6-a3ad-11ef-3e0b-53acf2fe8870
begin 
	using Pkg; Pkg.activate()

	using LilGuys
	using Arya
	using CairoMakie
end

# ╔═╡ b0c7516c-c5ec-4e29-a5ff-ecedccc0f50d
using CSV, DataFrames

# ╔═╡ 33b80e4c-09d8-414e-9464-c3fe891e7081
out = Output(".")

# ╔═╡ 987eceae-b2ff-4b0f-af5c-36f4a6b88a7d
out_nbody = Output("../all_galaxies_nbody")

# ╔═╡ 805a81e4-d76c-4306-891e-9d6cc105cb25
out_nolmc = Output("../all_galaxies_nolmc")

# ╔═╡ 1a05f83d-da8f-46f5-9c0d-2b2cf7f7e203
function split_dim(M; dim=2)
	N = size(M, dim)

	return [M[:, i, :] for i in 1:N]
end

# ╔═╡ 065439a0-2e10-463f-9a84-7a5ac1595ab0
positions = LilGuys.extract_vector(out, :positions) |> split_dim

# ╔═╡ 6f22e772-0996-417d-b477-d475caef1ba6
positions_nbody = LilGuys.extract_vector(out_nbody, :positions) |> split_dim

# ╔═╡ c5c269fe-1c5d-45ee-9e3b-ae674bc3721c
positions_nolmc = LilGuys.extract_vector(out_nolmc, :positions) |> split_dim

# ╔═╡ 51ca8ac6-1bbe-46db-85d5-7d0ef718cc70
velocities = LilGuys.extract_vector(out, :velocities) |> split_dim

# ╔═╡ afe4fc6a-e53a-4b05-bf02-de1171ff43a5
snap_i = out[1][sortperm(out[1].index)]

# ╔═╡ 02a57d8f-297c-4304-a6e1-88dfbe60aa8a
masses = LilGuys.extract(out_nbody[1], :masses)

# ╔═╡ 5535a429-502b-4609-acb8-f59dd1fcc57e
LilGuys.to_gaia(snap_i, add_centre=false, invert_velocity=true)

# ╔═╡ 6f6cdd69-e700-49ca-910c-38d95199ad8e
times = out.times

# ╔═╡ 5a27d8b0-c0b0-45e5-b2f8-cf9267627629
orbit_lmc = CSV.read("orbit_lmc.csv", DataFrame)

# ╔═╡ 76991988-2dc4-4eeb-9518-1f9ee4d9ddc5
pos_lmc = hcat(orbit_lmc.x, orbit_lmc.y, orbit_lmc.z)'

# ╔═╡ 6c1aa92c-7ccc-480a-a3b3-0b9418d032da
obs_props = CSV.read(ENV["DWARFS_ROOT"] * "/observations/pace_complete_nolmc.csv", DataFrame)


# ╔═╡ 744f3812-1548-4f79-b704-fb037157f0a0
galaxies = obs_props.name

# ╔═╡ 79f72640-4565-4028-873a-4da50f6147f5
obs_props

# ╔═╡ d2ad4ca8-0fac-406d-8f02-323134bd2772
LP = LilGuys.Plots

# ╔═╡ 6bd28b7c-d170-4914-b662-989a33a7a361
md"""
# ORBITSS
"""

# ╔═╡ 2b03c722-b5e7-4948-8704-6c4749c5a59f
r_max = 300

# ╔═╡ 9ad636bd-e51e-4c68-a885-1b4aff7b8b83
LP.plot_xyz(positions..., limits=tuple(fill((-r_max, r_max), 3)...))

# ╔═╡ aa463ef7-9348-499b-979c-93360234bb32
LP.plot_xyz(positions_nbody..., limits=tuple(fill((-r_max, r_max), 3)...))

# ╔═╡ 4d14a245-92a9-4b17-939b-6362ac33d635
LP.plot_xyz(positions_nolmc..., limits=tuple(fill((-r_max, r_max), 3)...))

# ╔═╡ bc657b37-ef38-4d02-99f8-e79df5afa1fa
r_nbody = [calc_r(pos, pos1) for (pos, pos1) in zip(positions, positions_nbody)]

# ╔═╡ 135a45a1-ab83-4916-bc39-478328627ba2
idx_worst = sortperm([r[end] for r in r_nbody], rev=true)

# ╔═╡ 7e6f9db8-cdf2-4f7d-af73-eee908fb68ca
let
	fig = Figure()
	ax = Axis(fig[1,1])

	for i in idx_worst[10:25]
		lines!(-times * T2GYR, r_nbody[i], label=galaxies[i])

	end
	
	Legend(fig[1, 2], ax)

	fig
end

# ╔═╡ 4e6ac68a-55f1-4d67-966f-d90664e4c5e9
r_diff_lmc= [calc_r(pos, pos1) for (pos, pos1) in zip(positions, positions_nolmc)]

# ╔═╡ 7a9d0ea8-f09b-4e95-a826-90752467b832
idx_worst_lmc = sortperm([r[end] for r in r_diff_lmc], rev=true)

# ╔═╡ 03d9395c-6c0d-4c2f-8113-5f7165858daf
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "diff in LMC & no LMC orbit" )

	for i in idx_worst_lmc[1:10]
		lines!(-times * T2GYR, r_diff_lmc[i], label=galaxies[i])

	end
	
	Legend(fig[1, 2], ax)

	fig
end

# ╔═╡ d16b5151-9a17-4b7b-ada9-98235cfac327
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "diff in LMC & no LMC orbit" )

	for i in idx_worst_lmc[50:end]
		lines!(-times * T2GYR, r_diff_lmc[i], label=galaxies[i])

	end
	
	Legend(fig[1, 2], ax)

	fig
end

# ╔═╡ 5cd7fd4e-43c0-40c5-a7d5-463346554044
md"""
# Orbits of select galaxies
"""

# ╔═╡ d4afbe65-d27f-4c7b-8b76-ee850e37d7da
function extract_galaxy(name)
	if name ∉ galaxies
		error("galaxy not known $name")
	end
	
	i = argmax(galaxies .== name)

	return (;
		galaxy = name,
		positions = positions[i],
		pos_nolmc = positions_nolmc[i],
		pos_lmc = positions[i] .- pos_lmc,
		positions_nbody = positions_nbody[i],
		velocities = velocities[i],
		mass = masses[i],
	)
end

# ╔═╡ 4698ceae-a636-45b3-a271-cc476309d4b6


# ╔═╡ 2bbefc84-e251-4c53-97ce-e354c44c9b54
function plot_radii(gal_orbits)

	r_mw = calc_r(gal_orbits.positions)
	r_nbody = calc_r(gal_orbits.positions_nbody)
	r_lmc = calc_r(gal_orbits.pos_lmc)
	r_nolmc = calc_r(gal_orbits.pos_nolmc)

	fig = Figure(
	)
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel = "distance / kpc",
		title=gal_orbits.galaxy,
		limits=(nothing, nothing, 0, nothing)
	)

	lines!(-T2GYR * times, r_mw, label="MW")
	lines!(-T2GYR * times, r_nbody, label="MW (+nbody)")
	lines!(-T2GYR * times, r_nolmc, label="MW (no LMC)")
	lines!(-T2GYR * times, r_lmc, label="LMC")

	axislegend()
	fig
end

# ╔═╡ 55186229-9200-4915-818c-2413e0a06cc7
function plot_orbit(name)
	gal = extract_galaxy(name)

	fig = LP.plot_xyz(gal.positions, gal.positions_nbody, gal.pos_nolmc, pos_lmc, labels=["galaxy", "+nbody", "-LMC", "lmc orbit"])
	Label(fig.layout[0, :], name)
	@info fig

	fig = LP.plot_xyz(gal.pos_lmc, gal.positions_nbody .- pos_lmc, labels=["galaxy", "+nbody"])
	Label(fig.layout[0, :], "$name in LMC frame")

	resize_to_layout!(fig)
	@info fig

	@info plot_radii(gal)

end

# ╔═╡ 7cb9f053-b1a9-4f18-82e0-9f780ef26b4f
pos_lmc

# ╔═╡ 19771df1-aac3-4809-ba0d-180cef942b39
plot_orbit("Sculptor")

# ╔═╡ 44729954-4083-42e5-9281-1281edc3f221
plot_orbit("SMC")

# ╔═╡ 3fcc543b-e5c9-4c56-afa0-916cac51b40c
plot_orbit("Hydrus I")

# ╔═╡ d4555cbd-526b-4d12-8ef7-4ad1d0aa40d3
plot_orbit("Horologium I")

# ╔═╡ 93e105b7-c788-442b-8db6-7d804e8960a5
plot_orbit("Carina III")

# ╔═╡ 68acaf76-7b30-444f-9150-a570e4c844cf
plot_orbit("Leo IV")

# ╔═╡ 0d4d84d5-e5b9-4be9-86fb-68c9abfe36dc
plot_orbit("Segue 1")

# ╔═╡ b95aa267-5d24-4176-a296-c92b10cb0eb7
plot_orbit("Sagittarius")

# ╔═╡ f30fd4ea-2856-48f1-a355-e13612b003c9
galaxies

# ╔═╡ f3e2b76b-87a4-4dc2-a4b8-e2cb5b590ce9
md"""
### Segue I
"""

# ╔═╡ d0f62081-2fb1-48d0-af65-9896cf296617
segue1 = extract_galaxy("Segue 1")

# ╔═╡ f5e6245e-fb40-4546-8f4d-01321e752210
LP.plot_xyz(segue1.positions .- segue1.pos_nolmc)

# ╔═╡ 7f864e00-03d2-4ce8-a1e5-f88aaad136f5
LP.plot_xyz(segue1.positions, segue1.pos_nolmc)

# ╔═╡ 31c88fd8-a5e7-4dac-817a-97bf4b50bfee
md"""
# Sculptor
"""

# ╔═╡ 621ff1be-26a8-41a1-9ec9-9311d5ca7832
scl = extract_galaxy("Sculptor")

# ╔═╡ 458ada32-0ca3-4bac-848e-a79daa2871f3
r_scl = [minimum(calc_r(pos, scl.positions_nbody)) for pos in positions_nbody]

# ╔═╡ 6b72ea04-fdba-4230-9202-efeb7b4f1f32
sort_scl = sortperm(r_scl)

# ╔═╡ 5e220294-4ff5-4dcf-bf4e-5f6a40d22d68
galaxies[sort_scl]

# ╔═╡ 9f4c8ba4-c9bc-4497-af8a-f8e8d37b0ebd
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "distance from Scl")

	for i in 2:10 # skip sculptor
		gal = extract_galaxy(galaxies[sort_scl][i])
		@info "galaxy $(gal.galaxy), mass $(gal.mass)"
		r = calc_r(scl.positions_nbody, gal.positions_nbody)

		lines!(-T2GYR * times, r, label = gal.galaxy)
	end

	Legend(fig[1,2], ax)
	fig
end
		
	

# ╔═╡ Cell order:
# ╠═bb4f68d6-a3ad-11ef-3e0b-53acf2fe8870
# ╠═b0c7516c-c5ec-4e29-a5ff-ecedccc0f50d
# ╠═33b80e4c-09d8-414e-9464-c3fe891e7081
# ╠═987eceae-b2ff-4b0f-af5c-36f4a6b88a7d
# ╠═805a81e4-d76c-4306-891e-9d6cc105cb25
# ╠═1a05f83d-da8f-46f5-9c0d-2b2cf7f7e203
# ╠═065439a0-2e10-463f-9a84-7a5ac1595ab0
# ╠═6f22e772-0996-417d-b477-d475caef1ba6
# ╠═c5c269fe-1c5d-45ee-9e3b-ae674bc3721c
# ╠═51ca8ac6-1bbe-46db-85d5-7d0ef718cc70
# ╠═afe4fc6a-e53a-4b05-bf02-de1171ff43a5
# ╠═02a57d8f-297c-4304-a6e1-88dfbe60aa8a
# ╠═5535a429-502b-4609-acb8-f59dd1fcc57e
# ╠═6f6cdd69-e700-49ca-910c-38d95199ad8e
# ╠═5a27d8b0-c0b0-45e5-b2f8-cf9267627629
# ╠═76991988-2dc4-4eeb-9518-1f9ee4d9ddc5
# ╠═744f3812-1548-4f79-b704-fb037157f0a0
# ╠═6c1aa92c-7ccc-480a-a3b3-0b9418d032da
# ╠═79f72640-4565-4028-873a-4da50f6147f5
# ╠═d2ad4ca8-0fac-406d-8f02-323134bd2772
# ╠═6bd28b7c-d170-4914-b662-989a33a7a361
# ╠═2b03c722-b5e7-4948-8704-6c4749c5a59f
# ╠═9ad636bd-e51e-4c68-a885-1b4aff7b8b83
# ╠═aa463ef7-9348-499b-979c-93360234bb32
# ╠═4d14a245-92a9-4b17-939b-6362ac33d635
# ╠═bc657b37-ef38-4d02-99f8-e79df5afa1fa
# ╠═135a45a1-ab83-4916-bc39-478328627ba2
# ╠═7e6f9db8-cdf2-4f7d-af73-eee908fb68ca
# ╠═4e6ac68a-55f1-4d67-966f-d90664e4c5e9
# ╠═7a9d0ea8-f09b-4e95-a826-90752467b832
# ╠═03d9395c-6c0d-4c2f-8113-5f7165858daf
# ╠═d16b5151-9a17-4b7b-ada9-98235cfac327
# ╠═5cd7fd4e-43c0-40c5-a7d5-463346554044
# ╠═d4afbe65-d27f-4c7b-8b76-ee850e37d7da
# ╠═4698ceae-a636-45b3-a271-cc476309d4b6
# ╠═2bbefc84-e251-4c53-97ce-e354c44c9b54
# ╠═55186229-9200-4915-818c-2413e0a06cc7
# ╠═7cb9f053-b1a9-4f18-82e0-9f780ef26b4f
# ╠═19771df1-aac3-4809-ba0d-180cef942b39
# ╠═44729954-4083-42e5-9281-1281edc3f221
# ╠═3fcc543b-e5c9-4c56-afa0-916cac51b40c
# ╠═d4555cbd-526b-4d12-8ef7-4ad1d0aa40d3
# ╠═93e105b7-c788-442b-8db6-7d804e8960a5
# ╠═68acaf76-7b30-444f-9150-a570e4c844cf
# ╠═0d4d84d5-e5b9-4be9-86fb-68c9abfe36dc
# ╠═b95aa267-5d24-4176-a296-c92b10cb0eb7
# ╠═f30fd4ea-2856-48f1-a355-e13612b003c9
# ╠═f3e2b76b-87a4-4dc2-a4b8-e2cb5b590ce9
# ╠═d0f62081-2fb1-48d0-af65-9896cf296617
# ╠═f5e6245e-fb40-4546-8f4d-01321e752210
# ╠═7f864e00-03d2-4ce8-a1e5-f88aaad136f5
# ╠═31c88fd8-a5e7-4dac-817a-97bf4b50bfee
# ╠═621ff1be-26a8-41a1-9ec9-9311d5ca7832
# ╠═458ada32-0ca3-4bac-848e-a79daa2871f3
# ╠═6b72ea04-fdba-4230-9202-efeb7b4f1f32
# ╠═5e220294-4ff5-4dcf-bf4e-5f6a40d22d68
# ╠═9f4c8ba4-c9bc-4497-af8a-f8e8d37b0ebd
