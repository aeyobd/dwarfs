### A Pluto.jl notebook ###
# v0.20.5

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

# ╔═╡ ec1e618d-d760-4649-b241-443983fcf327
md"""
# Introduction
"""

# ╔═╡ 0c7f39a3-bca6-4881-9a87-7dcecbabb23d
md"""
Main questions:
How much does the LMC affect galaxy orbits? Do LMC satellites stay bound to the LMC under our potential assumptions? Do other galaxy-galaxy interactions affect the orbits of any dwarf galaxies?

In this notebook, I analyze a few very simplistic models of galaxy orbits under the LMC potential. I use the Pace LG satellite catalogue for local group properties and select only MW dwarfs which have RV measurements. I then integrate these galaxies in the potential from Vasiliev+21 and analyze the resulting orbits here.


"""

# ╔═╡ dc02d2bc-edc2-48ae-ae80-f743d3815c94
md"""
I have three different models
- fiducial: This is the one I read in as simply `out` and store the positions in `positions`. This model uses the time-dependent Vasiliev+21 LMC potential and integrates each galaxy as a massless point particle
- `nbody`: this model's variables are always suffixed with a `_nbody`. In this model, I took the stellar masses from the Pace catalogue and used the mean stellar mass-v circ max relation from fahatti+18 and the mass-concentration relation from ludlow+16 to create NFW profiles for each galaxy. I then calculated the total mass within 10kpc assuming an exponential truncation radius of 10x r_s. I then used these masses in Gadget and assumed a softening length of 1kpc to order of magnitude approximate the galaxies. 
- `nolmc`: This one uses the initial MW profile from Vasiliev+21 and assumes the potential remains static and does not include a LMC component.

"""

# ╔═╡ 75a0e369-4a4a-48a0-8504-5a9744dc292a
CairoMakie.activate!(type=:png)

# ╔═╡ 17bc4283-3673-4ba6-8d8a-b4928d2745fe
md"""
## Data Loading
"""

# ╔═╡ 33b80e4c-09d8-414e-9464-c3fe891e7081
out = Output(".")

# ╔═╡ 987eceae-b2ff-4b0f-af5c-36f4a6b88a7d
out_nbody = Output("../vasiliev+21_nbody")

# ╔═╡ cb5066cd-cd18-4e1a-a54a-882f1ce26321
out_ep20 = Output("../EP2020")

# ╔═╡ 805a81e4-d76c-4306-891e-9d6cc105cb25
out_nolmc = Output("../vasiliev+21_nolmc")

# ╔═╡ 326a759d-11e0-4311-a575-ca8a96d3e058
md"""
## Analysis
"""

# ╔═╡ 015c5e87-b5f6-4766-a8ff-bd26dd5f6046
peris_ep20 = LilGuys.peris_apos(out_ep20)

# ╔═╡ 1a05f83d-da8f-46f5-9c0d-2b2cf7f7e203
function split_dim(M; dim=2)
	N = size(M, dim)

	return [M[:, i, :] for i in 1:N]
end

# ╔═╡ 065439a0-2e10-463f-9a84-7a5ac1595ab0
positions = LilGuys.extract_vector(out, :positions) |> split_dim

# ╔═╡ 70e32c80-e4cf-4f11-a236-70e9502019c5
positions_ep20 = LilGuys.extract_vector(out_ep20, :positions) |> split_dim

# ╔═╡ 6f22e772-0996-417d-b477-d475caef1ba6
positions_nbody = LilGuys.extract_vector(out_nbody, :positions) |> split_dim

# ╔═╡ c5c269fe-1c5d-45ee-9e3b-ae674bc3721c
positions_nolmc = LilGuys.extract_vector(out_nolmc, :positions) |> split_dim

# ╔═╡ 1eec667a-edb6-4aa1-873a-8852a5e0a823
sort(out[1].index) == out[7].index[sortperm(out[7].index)]

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
obs_props = CSV.read(ENV["DWARFS_ROOT"] * "/observations/observed_properties_complete.csv", DataFrame)

# ╔═╡ cab8b442-2b8e-468f-b20d-80a02ac54cb5
obs_props[obs_props.galaxyname .== "bootes1", :]

# ╔═╡ 227e725e-6f1e-437f-b626-51aa77a73e57
obs_props[obs_props.galaxyname .== "sculptor", :]

# ╔═╡ 744f3812-1548-4f79-b704-fb037157f0a0
galaxies = obs_props.galaxyname

# ╔═╡ 929d138d-3336-453c-95c0-0d8a9196f725
peris_ep20[galaxies .== ["sculptor"], :].t_last_peri * T2GYR

# ╔═╡ 6bd28b7c-d170-4914-b662-989a33a7a361
md"""
# ORBITS
"""

# ╔═╡ 2b03c722-b5e7-4948-8704-6c4749c5a59f
r_max = 300

# ╔═╡ 579b8946-44b0-4af9-b8b3-6ed4176851b3
md"""
The plots below simply show the orbits of all of the galaxies for each of the models. This is more illustrative and shows the wiide range of orbit. For the two models with the LMC, you can see a group of galaxies moving in a helix-like pattern around the LMC orbit.
"""

# ╔═╡ 9ad636bd-e51e-4c68-a885-1b4aff7b8b83
LilGuys.plot_xyz(positions..., limits=tuple(fill((-r_max, r_max), 3)...))

# ╔═╡ aa463ef7-9348-499b-979c-93360234bb32
LilGuys.plot_xyz(positions_nbody..., limits=tuple(fill((-r_max, r_max), 3)...))

# ╔═╡ 4d14a245-92a9-4b17-939b-6362ac33d635
LilGuys.plot_xyz(positions_nolmc..., limits=tuple(fill((-r_max, r_max), 3)...))

# ╔═╡ bc657b37-ef38-4d02-99f8-e79df5afa1fa
r_nbody = [LilGuys.radii(pos, pos1) for (pos, pos1) in zip(positions, positions_nbody)]

# ╔═╡ 135a45a1-ab83-4916-bc39-478328627ba2
idx_worst = sortperm([r[end] for r in r_nbody], rev=true)

# ╔═╡ deaa101a-3f04-4776-a9d4-099ccb59dcd8
md"""
The below plots explore the differences in orbits in each model by simply calculating the radius between the orbit positions at each time
"""

# ╔═╡ 7e6f9db8-cdf2-4f7d-af73-eee908fb68ca
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time / Gyr",
		ylabel = "distance nbody & LMC only orbits / kpc"
	)

	for i in idx_worst[1:10]
		lines!(-times * T2GYR, r_nbody[i], label=galaxies[i])

	end
	
	Legend(fig[1, 2], ax)

	fig
end

# ╔═╡ 4e6ac68a-55f1-4d67-966f-d90664e4c5e9
r_diff_lmc= [radii(pos, pos1) for (pos, pos1) in zip(positions, positions_nolmc)]

# ╔═╡ 7a9d0ea8-f09b-4e95-a826-90752467b832
idx_worst_lmc = sortperm([r[end] for r in r_diff_lmc], rev=true)

# ╔═╡ 03d9395c-6c0d-4c2f-8113-5f7165858daf
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "distance in LMC & no LMC orbits / kpc" )

	for i in idx_worst_lmc[1:10]
		lines!(-times * T2GYR, r_diff_lmc[i], label=galaxies[i])

	end
	
	Legend(fig[1, 2], ax)

	fig
end

# ╔═╡ d16b5151-9a17-4b7b-ada9-98235cfac327
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "distance in LMC & no LMC orbits / kpc" )

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
		pos_ep20 = positions_ep20[i],
		pos_lmc = positions[i] .- pos_lmc,
		positions_nbody = positions_nbody[i],
		velocities = velocities[i],
		mass = masses[i],
	)
end

# ╔═╡ 4698ceae-a636-45b3-a271-cc476309d4b6
radii

# ╔═╡ 2bbefc84-e251-4c53-97ce-e354c44c9b54
function plot_radii(gal_orbits)

	r_mw = radii(gal_orbits.positions)
	r_nbody = radii(gal_orbits.positions_nbody)
	r_lmc = radii(gal_orbits.pos_lmc)
	r_nolmc = radii(gal_orbits.pos_nolmc)
	r_lmc_nbody = radii(gal_orbits.positions_nbody, pos_lmc)

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
	lines!(-T2GYR * times, r_lmc_nbody, label="LMC (+nbody)")

	axislegend()
	fig
end

# ╔═╡ 55186229-9200-4915-818c-2413e0a06cc7
function plot_orbit(name)
	gal = extract_galaxy(name)

	fig = LilGuys.plot_xyz(gal.positions, gal.positions_nbody, gal.pos_nolmc, pos_lmc, labels=["galaxy", "+nbody", "-LMC", "lmc orbit"])
	Label(fig.layout[0, :], name)
	@info fig

	fig = LilGuys.plot_xyz(gal.pos_lmc, gal.positions_nbody .- pos_lmc, labels=["galaxy", "+nbody"])
	Label(fig.layout[0, :], "$name in LMC frame")

	resize_to_layout!(fig)
	@info fig

	@info plot_radii(gal)

end

# ╔═╡ 7cb9f053-b1a9-4f18-82e0-9f780ef26b4f
pos_lmc

# ╔═╡ 19771df1-aac3-4809-ba0d-180cef942b39
plot_orbit("sculptor")

# ╔═╡ 44729954-4083-42e5-9281-1281edc3f221
plot_orbit("smc")

# ╔═╡ 3fcc543b-e5c9-4c56-afa0-916cac51b40c
plot_orbit("hydrus1")

# ╔═╡ d4555cbd-526b-4d12-8ef7-4ad1d0aa40d3
plot_orbit("horologium1")

# ╔═╡ 93e105b7-c788-442b-8db6-7d804e8960a5
plot_orbit("carina3")

# ╔═╡ 68acaf76-7b30-444f-9150-a570e4c844cf
plot_orbit("leo4")

# ╔═╡ 0d4d84d5-e5b9-4be9-86fb-68c9abfe36dc
plot_orbit("segue1")

# ╔═╡ b95aa267-5d24-4176-a296-c92b10cb0eb7
plot_orbit("sagittarius2")

# ╔═╡ 5c588311-809b-4ba7-ac59-19b19b6e81a8
plot_orbit("ursa_minor")

# ╔═╡ 595c0582-76f9-413e-91a2-ded9794c31f0
plot_orbit("draco2")

# ╔═╡ c5c6c637-dbe1-4fe4-aabe-4785a19e8c45
plot_orbit("bootes1")

# ╔═╡ 7ca7ff6e-1300-4846-9a3a-9eda31e0c452
plot_orbit("bootes3")

# ╔═╡ d252cffe-176a-4d21-81ba-b27b07611a9d
plot_orbit("grus2")

# ╔═╡ a91cd2ea-62a0-41bb-bbb6-ab4f2f2224fb
plot_orbit("segue1")

# ╔═╡ 472a4d93-ebe2-4dd0-93ad-d980c826965e
plot_orbit("tucana2")

# ╔═╡ 38b2b26f-9a6f-45eb-90a8-293cffe6185a
plot_orbit("tucana3")

# ╔═╡ f30fd4ea-2856-48f1-a355-e13612b003c9
galaxies

# ╔═╡ 52ffacd6-d123-4aa2-9d6b-c186d065ed1d
1

# ╔═╡ f3e2b76b-87a4-4dc2-a4b8-e2cb5b590ce9
md"""
### Segue I
"""

# ╔═╡ d0f62081-2fb1-48d0-af65-9896cf296617
segue1 = extract_galaxy("segue1")

# ╔═╡ f5e6245e-fb40-4546-8f4d-01321e752210
LilGuys.plot_xyz(segue1.positions .- segue1.pos_nolmc)

# ╔═╡ 7f864e00-03d2-4ce8-a1e5-f88aaad136f5
LilGuys.plot_xyz(segue1.positions, segue1.pos_nolmc)

# ╔═╡ 31c88fd8-a5e7-4dac-817a-97bf4b50bfee
md"""
# Sculptor
"""

# ╔═╡ cc30ae64-56dc-4a87-8f22-d078a6cd30c7


# ╔═╡ 621ff1be-26a8-41a1-9ec9-9311d5ca7832
scl = extract_galaxy("sculptor")

# ╔═╡ 1f4a0dfc-8bc0-4484-a11f-ad19b3b3313a
smc = extract_galaxy("smc")

# ╔═╡ 458ada32-0ca3-4bac-848e-a79daa2871f3
r_scl = [minimum(radii(pos, scl.positions_nbody)) for pos in positions_nbody]

# ╔═╡ 6b72ea04-fdba-4230-9202-efeb7b4f1f32
sort_scl = sortperm(r_scl)

# ╔═╡ 5e220294-4ff5-4dcf-bf4e-5f6a40d22d68
galaxies[sort_scl]

# ╔═╡ c4f2d244-f9cc-4c4d-905f-fda73f4cb9f3
md"""
Which satellites become the closest to Sculptor?
"""

# ╔═╡ 9f4c8ba4-c9bc-4497-af8a-f8e8d37b0ebd
let
	fig = Figure(
		size=(400, 300)
	)
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "distance from Scl")

	for i in 2:10 # skip sculptor
		gal = extract_galaxy(galaxies[sort_scl][i])
		@info "galaxy $(gal.galaxy), mass $(gal.mass)"
		r = radii(scl.positions_nbody, gal.positions_nbody)

		lines!(-T2GYR * times, log10.(r), label = gal.galaxy)
	end

	Legend(fig[1,2], ax)
	fig
end
		
	

# ╔═╡ 125009ef-4f3c-4869-a892-306091769a1f
(LilGuys.TruncNFW(v_circ_max=20/V2KMS, r_circ_max=2.2, trunc=10), 10)

# ╔═╡ b2da19e3-e570-40e6-90fe-21c1afaae45a
filt_scl_peri = 20:40

# ╔═╡ a42e1019-ce0d-4c9e-a765-a72e4faaee1d
LilGuys.plot_xyz(scl.pos_lmc[:, filt_scl_peri], smc.pos_lmc[:, filt_scl_peri], labels=["Scl", "SMC"])

# ╔═╡ f21009ee-1385-4592-8b27-900d7dd9a07d
LilGuys.plot_xyz(scl.pos_ep20)

# ╔═╡ 2f695966-ebdb-4d32-b2df-b3b68fc3d2af
LilGuys.plot_xyz(scl.positions .- smc.positions, scl.positions .- pos_lmc, labels=["Scl", "SMC"])

# ╔═╡ 91635892-504e-4337-a9c9-d56466c9f10c
LilGuys.plot_xyz(scl.positions[:, filt_scl_peri], pos_lmc[:, filt_scl_peri], smc.positions[:, filt_scl_peri], labels=["Scl", "LMC", "SMC"], limits=((-20, 20), (-40, 0), (-80, -40)))

# ╔═╡ ed6a1c88-32ef-4932-b1ef-f026e8cc8ae2
LilGuys.plot_xyz(scl.positions[:, filt_scl_peri], pos_lmc[:, filt_scl_peri], smc.positions[:, filt_scl_peri],  limits=((-20, 20), (-40, 0), (-80, -40)), times=times[filt_scl_peri])

# ╔═╡ 3df80581-dbd3-4cc0-bb06-e468e5677075
LilGuys.plot_xyz(scl.pos_lmc[:, 20:40] .- smc.pos_lmc[:, 20:40])

# ╔═╡ c99520f6-42d7-4b9f-8041-60bb1abecb81
r_scl_smc = radii(scl.positions .- smc.positions)

# ╔═╡ 086237e5-9e88-4647-8804-0aac36291a04
times[20:30] * T2GYR

# ╔═╡ d0948ce0-79ee-4ed5-b346-4875b0345ead
argmin(r_scl_smc)

# ╔═╡ 6133036d-788e-4d82-9cd4-efe5dd68a2aa
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "distance from Scl",
		limits=(-0.5, 0, 0, 100)
	)

	t = -times * T2GYR

	lines!(t, radii(scl.positions .- pos_lmc), label="LMC")
	lines!(t, r_scl_smc, label="SMC")
	lines!(t, radii(scl.positions_nbody .- smc.positions_nbody), linestyle=:dot, label="SMC + nbody")
	lines!(t, radii(scl.positions_nbody .- pos_lmc), linestyle=:dot, label="LMC + nbody")

	axislegend()
	fig
end


# ╔═╡ bb9cdb67-01ff-4527-b0d4-3add90e3385b
md"""
### U Mi
"""

# ╔═╡ 55522f17-0ded-4df2-bec2-2ffe04ff3f86
umi = extract_galaxy("ursa_minor")

# ╔═╡ f473336d-2ee6-488e-b9b4-960a1133cba0
LilGuys.plot_xyz(umi.positions .- umi.pos_nolmc)

# ╔═╡ 552ed892-fd37-42ca-a22f-942a48b6ca15
LilGuys.plot_xyz(umi.positions, umi.pos_nolmc)

# ╔═╡ 8d9666be-f03b-4e90-bf6c-ecea5ea1fb9e
md"""
## Boo 1
"""

# ╔═╡ d7bdd638-7e2e-4383-8b0a-243666fcfa69
boo1 = extract_galaxy("bootes1")

# ╔═╡ 364ff46e-1a0c-41f3-a38d-a0a08b424511
r_boo1 = [minimum(radii(pos, boo1.positions_nbody)) for pos in positions_nbody]

# ╔═╡ a205beda-024b-4a5d-8b10-47ec03520b9e
sort_boo1 = sortperm(r_boo1)

# ╔═╡ c7d35c3b-53f0-48e8-8f30-da15a97e0d63
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "distance from Boo I",
	)

	for i in 2:10 # skip sculptor
		gal = extract_galaxy(galaxies[sort_boo1][i])
		@info "galaxy $(gal.galaxy), mass $(gal.mass)"
		r = radii(boo1.positions_nbody, gal.positions_nbody)

		lines!(-T2GYR * times, log10.(r), label = "$(gal.galaxy), $(gal.mass)")
	end

	Legend(fig[1,2], ax)
	fig
end
		

# ╔═╡ 6709941b-63b4-4ea4-9255-1dfda12e9428
boo3 = extract_galaxy("bootes3")

# ╔═╡ 156f8f25-edce-4334-95b7-1b12c7b282b7
r_boo3 = [minimum(radii(pos, boo3.positions_nbody)) for pos in positions_nbody]

# ╔═╡ 5bc24806-3f5c-4fb2-8f20-6cbcc5ac4d25
sort_boo3 = sortperm(r_boo3)

# ╔═╡ 6fd2ceb8-cd02-4fdd-a9d4-ae63cad2518d
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "distance from Boo 3",
	)

	for i in 2:10 # skip sculptor
		gal = extract_galaxy(galaxies[sort_boo3][i])
		@info "galaxy $(gal.galaxy), mass $(gal.mass)"
		r = radii(boo3.positions_nbody, gal.positions_nbody)

		lines!(-T2GYR * times, log10.(r), label = "$(gal.galaxy), $(round(gal.mass, digits=2))")
	end

	Legend(fig[1,2], ax)
	fig
end

# ╔═╡ f576d1c4-3e35-4b65-9d1c-927febfb2dfc
md"""
# Adding orbits to properties
"""

# ╔═╡ c15caa3a-b8f1-40a8-8145-1ad1bc725124
let
	global new_properties = copy(obs_props)

	peris = Float64[]
	apos = Float64[]
	
	perilmc = Float64[]
	apolmc = Float64[]

	peris_ep = Float64[]
	apos_ep = Float64[]

	for (i, gal) in enumerate(galaxies)
		orbit = extract_galaxy(gal)
		peri, apo = extrema(radii(orbit.positions))
		push!(peris, peri)
		push!(apos, apo)

		peri, apo = extrema(radii(orbit.pos_ep20))
		push!(peris_ep, peri)
		push!(apos_ep, apo)

		@info (peri, peris_ep20.pericentre[i])
		
		peri, apo = extrema(radii(orbit.positions, pos_lmc))
		push!(perilmc, peri)
		push!(apolmc, apo)

	end

	new_properties[!, :peri] = peris_ep20.pericentre
	new_properties[!, :t_last_peri] = peris_ep20.t_last_peri
	new_properties[!, :t_last_apo] = peris_ep20.t_last_apo
	new_properties[!, :apos] = peris_ep20.apocentre
	new_properties[!, :peri_v21] = perilmc
	new_properties[!, :apo_v21] = apolmc

	new_properties[!, :perilmc_v21] = perilmc
	new_properties[!, :apolmc_v21] = apolmc
	
	new_properties
end

# ╔═╡ 0dca3434-ecdb-4da3-a0dc-8b8e4725c35d
CSV.write("properties_w_orbits.csv", new_properties)

# ╔═╡ Cell order:
# ╟─ec1e618d-d760-4649-b241-443983fcf327
# ╟─0c7f39a3-bca6-4881-9a87-7dcecbabb23d
# ╟─dc02d2bc-edc2-48ae-ae80-f743d3815c94
# ╠═bb4f68d6-a3ad-11ef-3e0b-53acf2fe8870
# ╠═b0c7516c-c5ec-4e29-a5ff-ecedccc0f50d
# ╠═75a0e369-4a4a-48a0-8504-5a9744dc292a
# ╟─17bc4283-3673-4ba6-8d8a-b4928d2745fe
# ╠═33b80e4c-09d8-414e-9464-c3fe891e7081
# ╠═987eceae-b2ff-4b0f-af5c-36f4a6b88a7d
# ╠═cb5066cd-cd18-4e1a-a54a-882f1ce26321
# ╠═805a81e4-d76c-4306-891e-9d6cc105cb25
# ╠═326a759d-11e0-4311-a575-ca8a96d3e058
# ╠═015c5e87-b5f6-4766-a8ff-bd26dd5f6046
# ╠═929d138d-3336-453c-95c0-0d8a9196f725
# ╠═1a05f83d-da8f-46f5-9c0d-2b2cf7f7e203
# ╠═065439a0-2e10-463f-9a84-7a5ac1595ab0
# ╠═70e32c80-e4cf-4f11-a236-70e9502019c5
# ╠═6f22e772-0996-417d-b477-d475caef1ba6
# ╠═c5c269fe-1c5d-45ee-9e3b-ae674bc3721c
# ╠═1eec667a-edb6-4aa1-873a-8852a5e0a823
# ╠═51ca8ac6-1bbe-46db-85d5-7d0ef718cc70
# ╠═afe4fc6a-e53a-4b05-bf02-de1171ff43a5
# ╠═02a57d8f-297c-4304-a6e1-88dfbe60aa8a
# ╠═5535a429-502b-4609-acb8-f59dd1fcc57e
# ╠═6f6cdd69-e700-49ca-910c-38d95199ad8e
# ╠═5a27d8b0-c0b0-45e5-b2f8-cf9267627629
# ╠═76991988-2dc4-4eeb-9518-1f9ee4d9ddc5
# ╠═cab8b442-2b8e-468f-b20d-80a02ac54cb5
# ╠═227e725e-6f1e-437f-b626-51aa77a73e57
# ╠═744f3812-1548-4f79-b704-fb037157f0a0
# ╠═6c1aa92c-7ccc-480a-a3b3-0b9418d032da
# ╟─6bd28b7c-d170-4914-b662-989a33a7a361
# ╠═2b03c722-b5e7-4948-8704-6c4749c5a59f
# ╟─579b8946-44b0-4af9-b8b3-6ed4176851b3
# ╠═9ad636bd-e51e-4c68-a885-1b4aff7b8b83
# ╠═aa463ef7-9348-499b-979c-93360234bb32
# ╠═4d14a245-92a9-4b17-939b-6362ac33d635
# ╠═bc657b37-ef38-4d02-99f8-e79df5afa1fa
# ╠═135a45a1-ab83-4916-bc39-478328627ba2
# ╟─deaa101a-3f04-4776-a9d4-099ccb59dcd8
# ╟─7e6f9db8-cdf2-4f7d-af73-eee908fb68ca
# ╠═4e6ac68a-55f1-4d67-966f-d90664e4c5e9
# ╠═7a9d0ea8-f09b-4e95-a826-90752467b832
# ╟─03d9395c-6c0d-4c2f-8113-5f7165858daf
# ╟─d16b5151-9a17-4b7b-ada9-98235cfac327
# ╟─5cd7fd4e-43c0-40c5-a7d5-463346554044
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
# ╠═5c588311-809b-4ba7-ac59-19b19b6e81a8
# ╠═595c0582-76f9-413e-91a2-ded9794c31f0
# ╠═c5c6c637-dbe1-4fe4-aabe-4785a19e8c45
# ╠═7ca7ff6e-1300-4846-9a3a-9eda31e0c452
# ╠═d252cffe-176a-4d21-81ba-b27b07611a9d
# ╠═a91cd2ea-62a0-41bb-bbb6-ab4f2f2224fb
# ╠═472a4d93-ebe2-4dd0-93ad-d980c826965e
# ╠═38b2b26f-9a6f-45eb-90a8-293cffe6185a
# ╠═f30fd4ea-2856-48f1-a355-e13612b003c9
# ╠═52ffacd6-d123-4aa2-9d6b-c186d065ed1d
# ╠═f3e2b76b-87a4-4dc2-a4b8-e2cb5b590ce9
# ╠═d0f62081-2fb1-48d0-af65-9896cf296617
# ╠═f5e6245e-fb40-4546-8f4d-01321e752210
# ╠═7f864e00-03d2-4ce8-a1e5-f88aaad136f5
# ╠═31c88fd8-a5e7-4dac-817a-97bf4b50bfee
# ╠═cc30ae64-56dc-4a87-8f22-d078a6cd30c7
# ╠═621ff1be-26a8-41a1-9ec9-9311d5ca7832
# ╠═1f4a0dfc-8bc0-4484-a11f-ad19b3b3313a
# ╠═458ada32-0ca3-4bac-848e-a79daa2871f3
# ╠═6b72ea04-fdba-4230-9202-efeb7b4f1f32
# ╠═5e220294-4ff5-4dcf-bf4e-5f6a40d22d68
# ╟─c4f2d244-f9cc-4c4d-905f-fda73f4cb9f3
# ╠═9f4c8ba4-c9bc-4497-af8a-f8e8d37b0ebd
# ╠═125009ef-4f3c-4869-a892-306091769a1f
# ╠═b2da19e3-e570-40e6-90fe-21c1afaae45a
# ╠═a42e1019-ce0d-4c9e-a765-a72e4faaee1d
# ╠═f21009ee-1385-4592-8b27-900d7dd9a07d
# ╠═2f695966-ebdb-4d32-b2df-b3b68fc3d2af
# ╠═91635892-504e-4337-a9c9-d56466c9f10c
# ╠═ed6a1c88-32ef-4932-b1ef-f026e8cc8ae2
# ╠═3df80581-dbd3-4cc0-bb06-e468e5677075
# ╠═c99520f6-42d7-4b9f-8041-60bb1abecb81
# ╠═086237e5-9e88-4647-8804-0aac36291a04
# ╠═d0948ce0-79ee-4ed5-b346-4875b0345ead
# ╠═6133036d-788e-4d82-9cd4-efe5dd68a2aa
# ╠═bb9cdb67-01ff-4527-b0d4-3add90e3385b
# ╠═55522f17-0ded-4df2-bec2-2ffe04ff3f86
# ╠═f473336d-2ee6-488e-b9b4-960a1133cba0
# ╠═552ed892-fd37-42ca-a22f-942a48b6ca15
# ╠═8d9666be-f03b-4e90-bf6c-ecea5ea1fb9e
# ╠═d7bdd638-7e2e-4383-8b0a-243666fcfa69
# ╠═364ff46e-1a0c-41f3-a38d-a0a08b424511
# ╠═a205beda-024b-4a5d-8b10-47ec03520b9e
# ╠═c7d35c3b-53f0-48e8-8f30-da15a97e0d63
# ╠═6709941b-63b4-4ea4-9255-1dfda12e9428
# ╠═156f8f25-edce-4334-95b7-1b12c7b282b7
# ╠═5bc24806-3f5c-4fb2-8f20-6cbcc5ac4d25
# ╠═6fd2ceb8-cd02-4fdd-a9d4-ae63cad2518d
# ╠═f576d1c4-3e35-4b65-9d1c-927febfb2dfc
# ╠═c15caa3a-b8f1-40a8-8145-1ad1bc725124
# ╠═0dca3434-ecdb-4da3-a0dc-8b8e4725c35d
