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

# ╔═╡ 061b1886-1878-11ef-3806-b91643300982
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	using LilGuys

	using Arya
end

# ╔═╡ cd43d649-1d52-46a7-a621-68c8122fce6d
using PyFITS

# ╔═╡ 6b3df4f3-de5b-40b1-8aec-1cfce9aa843a
using PlutoUI

# ╔═╡ 2c702eb7-ebb6-44c9-8e01-ca52d011c014
using HDF5

# ╔═╡ 9e2cbc6b-58ec-45d0-9afa-568a7bc8a33e
using Printf

# ╔═╡ 8b41af50-9ae0-475b-bacc-3799e2949b30
md"""
Analyzes the orbit of a n-body halo in a gravitational potential.
Requires the centres to be calculated prior.

This notebook only uses the centres of the orbit, but then calculates useful quantities such as the pericentre, the orbital period, the time of each peri and apocentre and outputs this into the model analysis directory.

"""

# ╔═╡ 27577252-53dc-415c-b9a1-82155ef9e4ca
md"""
# Setup
"""

# ╔═╡ 882d4fc5-07ae-4b06-8da5-67f0894595db
import LinearAlgebra: dot

# ╔═╡ ab57edae-2292-4aef-9c1f-53802dbc0600
import TOML

# ╔═╡ 7bf1a839-b24e-4478-bdd3-200989779f68
save = Makie.save

# ╔═╡ 7834415a-b479-495b-998a-1d12f42f0dc6
Slider = PlutoUI.Slider

# ╔═╡ 56d15948-4fd8-4541-9925-99837d9584f5
CairoMakie.activate!(type=:png)

# ╔═╡ 643cd0bf-77b3-4201-9ff7-09dd5aee277c
md"""
# inputs
"""

# ╔═╡ 0e2e7f93-09d1-4902-ad24-223df50d37cb
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

# ╔═╡ d94663e8-b30e-4712-8a3e-6ef7f954f141
@bind inputs confirm(notebook_inputs(;
	galaxyname = TextField(default="ursa_minor"),
	haloname = TextField(default="1e6_v37_r5.0"),
	orbitname = TextField(default="orbit_"),
	t_min = NumberField(-10:0.1:10),
))

# ╔═╡ 3953b76f-a726-4211-a6b8-5cf38149dcdf
galaxyname = inputs.galaxyname

# ╔═╡ 24ed3bc1-37c6-4601-b941-0780c53a9630
haloname = inputs.haloname

# ╔═╡ 079afa70-c7dc-4347-9ae3-8459ba2fa941
orbitname =  inputs.orbitname

# ╔═╡ 94344455-d1d2-4ef9-af11-2d79ee4729ee
t_min = inputs.t_min

# ╔═╡ dd56b7ec-be11-447f-acc1-12750d82879b
md"""
##### the below should hopefully be always the same
"""

# ╔═╡ bc6deac8-b70a-483b-9fd7-1413c6f17aa7
mc_name = ""

# ╔═╡ 69d83e00-7eb6-4271-838f-80e4d1654dac
modelname = "$galaxyname/$haloname/$orbitname"

# ╔═╡ ac2c7484-9acd-4fda-9699-fdf17da507c2
parentdir = ENV["DWARFS_ROOT"]

# ╔═╡ d142b7bd-3002-4331-a725-577873c42f28
properties_file = joinpath(parentdir, "analysis", modelname, "simulation/orbit.toml")

# ╔═╡ 61d788ec-3518-4e3a-8eef-59c86ae5fc1a
obs_file =  "$parentdir/observations/$galaxyname/observed_properties.toml"

# ╔═╡ 0dd476fd-be53-4e9b-a686-a4462485c64c
orbit_file = joinpath(parentdir, "analysis", modelname, "simulation/orbit.csv")

# ╔═╡ 2bc762ad-e590-443e-b3c2-91dc42a8a4d9
outfile = joinpath(parentdir, "analysis", modelname, "orbital_properties.toml")

# ╔═╡ bb9ef388-bb5a-45a3-836e-4c06dbe0ab65
centresfile = joinpath(parentdir, "analysis", modelname, "centres.hdf5")

# ╔═╡ 9c427388-c657-4bb7-bc0a-b4de3597c645
skyorbit_outfile = joinpath(parentdir, "analysis", modelname, "skyorbit.fits")

# ╔═╡ 30969f77-667e-4ae4-9897-82c1c1182652
md"""
# File loading
"""

# ╔═╡ 96a57df5-a7b7-447a-a4a6-2b05e391a5c6
begin 
	obs_today = TOML.parsefile(properties_file)
	if isdefined(@__MODULE__, :ic)
		println("updating")
		for (k, v) in ic
			obs_today[k] = ic[k]
		end
	end

	rh = TOML.parsefile(obs_file)["r_h"]

	obs_today["ra_err"] = rh / 60 
	obs_today["dec_err"] = rh / 60
	obs_today
end

# ╔═╡ b250bf10-c228-4b14-938a-35561ae871d7
h5open(centresfile, "r") do  f
	global x_cen, v_cen, t
	x_cen = read(f["positions"])
	v_cen = read(f["velocities"])
	t = read(f["times"])
end

# ╔═╡ bb0cb8c2-2cbd-4205-a39e-4c4c0ff98b8a
begin 
	orbit_expected = CSV.read(orbit_file, DataFrame)
	x_cen_exp = transpose(hcat(orbit_expected.x, orbit_expected.y, orbit_expected.z))
	v_cen_exp = transpose(hcat(orbit_expected.v_x, orbit_expected.v_y, orbit_expected.v_z))
	orbit_expected.t .-= orbit_expected.t[1]
end

# ╔═╡ 08c3df42-738b-47c4-aa6b-fc39a9cfc02f
md"""
# plots
"""

# ╔═╡ a1c992c6-ad12-4968-b105-adfa1f327e76
let
	fig = LilGuys.plot_xyz(x_cen, x_cen_exp, labels=["n body", "point particle"])
	@savefig "centre_xyz"
	fig
end

# ╔═╡ 5255c605-56ea-4eb3-bd20-5134e3a96705
LilGuys.plot_xyz(v_cen, v_cen_exp, units=" / km/ s")

# ╔═╡ 15293cb8-61d3-478d-a2ae-5a5b2006db44
T2GYR = LilGuys.T2GYR

# ╔═╡ f88b909f-c3dc-41e0-bdb1-25e229964d27
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "radius / kpc")
	r = radii(x_cen)
	lines!(t * T2GYR, r, label="n-body")
	lines!(T2GYR*(orbit_expected.t), radii(x_cen_exp),
		label="point particle"
	)

	axislegend(ax)

	@savefig "centre_r_t"
	fig
end

# ╔═╡ 7d29a3bd-dc83-4eb3-ae65-fce5270ed8d5
md"""
# Sky Properties
"""

# ╔═╡ f134b3ce-53f0-47e9-84e9-1e73064d5191
snap_cen = Snapshot(x_cen, v_cen, zeros(size(x_cen, 2)))

# ╔═╡ 64e558da-2928-4815-ad5a-7528516311f9
obs_c_gr = LilGuys.to_gaia(snap_cen, SkyFrame=LilGuys.GSR, add_centre=false)

# ╔═╡ defc4184-2613-4480-9adf-fa135f168382
obs_c = LilGuys.to_gaia(snap_cen, add_centre=false)

# ╔═╡ 71ef4c9b-1284-4451-bcdd-eff79d334539


# ╔═╡ 5ec062c3-3815-4cf7-b45a-f97332d1b800
snap_cen.masses

# ╔═╡ a179323f-4878-4021-b8d4-69ca733658cb
function calc_χ2s(obs_c, obs_today)
	χ2 = zeros(size(obs_c, 1))
	for name in ["ra", "dec", "pmra", "pmdec", "distance", "radial_velocity"]
		μ = obs_today[name]
		x = obs_c[:, name]
		σ = obs_today[name * "_err"]
		χ2 .+= @. (x -μ)^2/(σ)^2
	end
	return χ2 ./ 6
end

# ╔═╡ ecf7c820-81a4-4cb7-a794-b7835c77811e
χ2 = calc_χ2s(obs_c, obs_today)

# ╔═╡ 319b905c-2d08-4a95-9d95-9cd26e2f5b1f
times = t * T2GYR

# ╔═╡ 9d60d54b-70e8-4b3c-a7c7-7caaa2f94a1c
t_end = times[times .> t_min][argmin(χ2[times .> t_min])]

# ╔═╡ 7646ea5b-f1b1-4934-be68-330139f7f838
length(χ2)

# ╔═╡ d9df3376-6ca1-4701-afb5-2df994bb3442
idx_f = searchsortedfirst(times, t_end)

# ╔═╡ aa2c3a93-19a3-43d8-82de-ae6ed8c4b9f7
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="R / kpc", ylabel="z / kpc",
		aspect=DataAspect(),
		xgridvisible=false,
		ygridvisible=false,
	)
	x = x_cen[1, :]
	y = x_cen[2, :]
	z = x_cen[3, :]
	R = @. sqrt(x^2 + y^2)
	lines!(R, z, color=times)

	scatter!(R[idx_f], z[idx_f], color=COLORS[2])
	scatter!(R[1], z[1], color = COLORS[3], marker=:rtriangle, )

	@savefig "R_z_centre"
	fig
end

# ╔═╡ cdde517a-1b3e-4d96-9156-4a8f72b795e9
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="time",
		ylabel="chi2 property fit",
		yscale=log10,
			  yticks = Makie.automatic,
	)
	
	lines!(t * T2GYR, χ2,)
	scatter!(t[idx_f] * T2GYR, χ2[idx_f])

	fig
end

# ╔═╡ cb601d4f-0fd6-4d1d-8bf0-910471c729c6
χ2[idx_f]

# ╔═╡ 9530e936-1225-4cfc-aa9a-bf7644d612f5
r = radii(x_cen)

# ╔═╡ 6d0612d6-8609-4832-80ae-4e8e78c557cc
minimum(r)

# ╔═╡ 7a30bd90-946e-418c-8339-be64c37cda76
vr = [dot(x_cen[:, i], v_cen[:, i]) / r[i] for i in 1:size(x_cen, 2)]

# ╔═╡ 0c69519c-9650-46b9-89f9-cc37227f5b1a
v_cen

# ╔═╡ d95c457b-c9be-4570-bc90-b4bbb7de56e2
plot(t, vr)

# ╔═╡ 88a9bb2d-9f86-4a96-a3b8-7a057480c7c3
"""

Given the velocities, finds when the velocity changes sign, 
i.e. a local extrema in r
"""
function find_all_peris(r)
	is_local_max(i) = (r[i] >= r[i-1]) && (r[i] >= r[i+1])
	is_local_min(i) = (r[i] <= r[i-1]) && (r[i] <= r[i+1])
	
	peris = []
	apos = []
	
	for i in 2:length(r)-1
		if is_local_max(i)
			push!(apos, i)
		elseif is_local_min(i)
			push!(peris, i)
		end
	end
	return peris, apos
end

# ╔═╡ 8d1508af-1715-4ef6-aab9-e95a02265913
idx_peris, idx_apos = find_all_peris(r)

# ╔═╡ ddd74bbb-df27-4150-b497-b1df5243f518
r[idx_apos]

# ╔═╡ efcbae60-cf7c-4e74-aae4-39d19b74b6fa
idx_peri = maximum(idx_peris[idx_peris .< idx_f])

# ╔═╡ d81c1455-728a-4023-ad65-e3cce37a69f9
r[idx_peris[end]], minimum(r)

# ╔═╡ a5ce5442-73ca-4aaf-915a-72fe9936e791
d_idx = 20

# ╔═╡ 7e14d1c7-b37f-4c0c-8e2a-1ac7cce27aaa
begin
	peri_filt = idx_f-d_idx:idx_f
	t_last_peri_arg = argmin(r[peri_filt])
	t_last_peri = t[peri_filt[t_last_peri_arg]] * T2GYR
	delta_t_peri = t[idx_f] * T2GYR - t_last_peri
end

# ╔═╡ f64dcd49-b0ca-4319-b615-5520b23d7818
lcm

# ╔═╡ 04d29fcb-70a0-414b-a487-7a18c44b9d58
let
	fig = Figure(size=(400, 150))
	ax = Axis(fig[1,1], xlabel="time / Gyr", ylabel = "radius / kpc")
	r = radii(x_cen)
	lines!(t * T2GYR, r)
	scatter!(t[idx_f] * T2GYR, r[idx_f], 
		label="adpoted end", marker=:rect
	)
	scatter!(t[idx_f] * T2GYR, radii(x_cen_exp)[end], 
		marker=:+, markersize=10, label="expected"
	)
	
	scatter!(t[idx_peris] * T2GYR, r[idx_peris], 
		label="pericentrs"
	)

	scatter!(t[idx_apos] * T2GYR, r[idx_apos], 
		label="apocentres"
	)
	
	scatter!(t[idx_peri] * T2GYR, r[idx_peri], 
		label=" last pericentre"
	)
	
	Legend(fig[1, 2], ax)

	@savefig "r_t_orbit"

	fig
end

# ╔═╡ af8a50bd-e761-4439-9fc9-80048c264d5b
begin 
	if idx_peri > 0
		t_peri = T2GYR * t[idx_peri]
		r_peri = r[idx_peri]

	else 
		t_peri = NaN
		r_peri = NaN
	end

end

# ╔═╡ 73bb2d61-37f3-4782-ae89-d36d1ff8f9ff
t_f = t[idx_f] * T2GYR

# ╔═╡ 14eebce8-04f7-493b-824a-7808c7fa35dd
md"""
# validating today
"""

# ╔═╡ 54cf5233-a955-4831-86ad-23b72f15789d
for property in ["ra", "dec", "pmra", "pmdec", "distance", "radial_velocity"]
	dy = (obs_c[idx_f, property] - obs_c[idx_f-1, property])/2
	
	@printf "%20s\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n" property obs_today[property] obs_today[property * "_err"]  obs_c[idx_f, property]  dy
end

# ╔═╡ cdabdc7d-76a1-45f5-b83a-2454576d3964
let

	fig = Figure(size=(4*72, 6*72))
	for i in 1:3
		x, y = [("ra", "dec"), ("pmra", "pmdec"), ("distance", "radial_velocity")][i]
		ax = Axis(fig[i,1],
			xlabel=x,
			ylabel=y
		)
	
		idx = idx_f-3:idx_f
	
		scatterlines!(obs_c[idx, x], obs_c[idx, y], color=log10.(χ2[idx]))
		
		errorscatter!([obs_today[x]], [obs_today[y]],
			xerror=[obs_today[x * "_err"]], yerror=[obs_today[y * "_err"]]
		)
	
	end

	@savefig "skyorbit_agreement_today"
	fig
end

# ╔═╡ 891f31ed-5565-4f43-ab82-b8bc3a76c1cf
let

	fig = Figure(size=(3*72, 6*72))
	for i in 1:3
		x, y = [("ra", "dec"), ("pmra", "pmdec"), ("distance", "radial_velocity")][i]
		ax = Axis(fig[i,1],
			xlabel=x,
			ylabel=y
		)
	
		idx = idx_f-3:idx_f
	
		scatterlines!(obs_c[idx, x], obs_c[idx, y], color=times[idx])
		
		errorscatter!([obs_today[x]], [obs_today[y]],
			xerror=[obs_today[x * "_err"]], yerror=[obs_today[y * "_err"]]
		)
	
	end

	@savefig "skyorbit_agreement_today"
	fig
end

# ╔═╡ ca193997-61c3-4302-a200-4e7d4e777521
Arya.UNITS_PER_INCH

# ╔═╡ 293090a7-8ee0-442e-b70f-e6b7750ab319
lerps = Dict(c => LilGuys.lerp(times, obs_c[:, c]) for c in ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity"])

# ╔═╡ 052fb31a-788a-485f-b8d2-1cf49f6ffb4b


# ╔═╡ e96a9789-22b8-433a-8a3c-182d7ec4e82a
let
	global obs_c_new, t_best, idx_best

	
	ts = LinRange(times[idx_f-1], times[idx_f], 1000)

	obs_c_new = DataFrame()

	for (k, v) in lerps
		obs_c_new[!, k] = v.(ts)
	end
	obs_c_new[!, :times] = ts

	obs_c_new

	χ2s_new = calc_χ2s(obs_c_new, obs_today)

	idx_best = argmin(χ2s_new)
	t_best = ts[idx_best]
	println(minimum(χ2s_new))

	lines(ts, log10.(χ2s_new))

end

# ╔═╡ 225a8adb-82ee-4479-806e-15796e2b08e2
let

	fig = Figure(size=(3*72, 8*72))
	for i in 1:6
		x = "times"
		y = ["ra", "dec", "pmra", "pmdec", "distance", "radial_velocity"][i]

		ax = Axis(fig[i,1],
			xlabel=x,
			ylabel=y
		)
	
		idx = idx_f-3:idx_f
	
		scatterlines!(times[idx], obs_c[idx, y], color=log10.(χ2[idx]))
		
		errorscatter!([times[idx_f]], [obs_today[y]], yerror=[obs_today[y * "_err"]]
		)
	
		errorscatter!([obs_c_new[idx_best, :times]], [obs_today[y]], yerror=[obs_today[y * "_err"]]
		)


		dy = obs_c_new[idx_best, y] - obs_today[y]
		χ2i = (dy) / obs_today[y*"_err"]
		
		text!(0.1, 0.5, space=:relative, text="dy = $(round(dy, digits=2)) \n σ=$(round(χ2i, digits=2))")
		if i < 6
			hidexdecorations!(ax, grid=false, ticks=false)
		end

		
	
	end

	@savefig "skyorbit_agreement_byvar"
	fig
end

# ╔═╡ 5858e6d0-a501-43ee-9669-6a225a1df1c9
Makie.spaces()

# ╔═╡ a48b78df-588d-4f94-9252-badd27179deb
obs_c

# ╔═╡ 3448ffc5-41e6-4208-b11f-2c00168bf50a
md"""
# Stream Coordinate frame
"""

# ╔═╡ 66435478-0619-47aa-a659-c06089951f72
ra0 = obs_c_gr.ra[idx_f]

# ╔═╡ 78271e36-12b7-4edc-bfb0-20ecc597ab20
dec0 = obs_c_gr.dec[idx_f]

# ╔═╡ 661ca87c-c8da-49b1-b8a3-72c81050590b
θ0 = atand(obs_c_gr.pmra[idx_f], obs_c_gr.pmdec[idx_f])

# ╔═╡ 77aa1d73-0c90-4f6c-9383-99e9e9f0379a
sind(θ0), cosd(θ0)

# ╔═╡ 2318b289-d6c4-44bf-b2ff-36f448faf97b
obs_c_gr.pmra[idx_f], obs_c_gr.pmdec[idx_f]

# ╔═╡ afdb058a-ebbd-4f07-b0d8-a85bb1070737
90 - atand(0, 1)

# ╔═╡ 31a11704-1dad-4007-b704-9312b81a5bad
begin
	dr = 0
	dθ = 180
	ra1 = ra0 + sind(θ0 + dθ) * dr / cosd(dec0)
	dec1 = dec0 + cosd(θ0 + dθ) * dr

	ra1 = round(ra1, digits=3)
	dec1 = round(dec1, digits=3)
end

# ╔═╡ bc1a33f4-cf2e-44a6-a20b-1344f47e75c6
md"""
Blue point below is at 


( $ra1, $dec1 )

"""

# ╔═╡ 2d1017a9-03a0-4aa2-829f-72174eaaa363
θ0

# ╔═╡ 63d2d908-2f12-483b-bbad-833b2aecc4e3
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="ra",
		ylabel="dec",
		xreversed=true
	)
	
	idx = idx_f-10:idx_f

	h = scatter!(obs_c[idx, :ra], obs_c[idx, :dec], color=t[idx] * T2GYR)

	arrows!([ra0], [dec0], [sind(θ0)], [cosd(θ0)])
	scatter!(ra1, dec1)

	Colorbar(fig[1, 2], h, label="time / Gyr")

	fig
end

# ╔═╡ 5676b72f-a981-4efb-8328-c25d3c5d6fb0
periods = [diff(times[idx_peris]); diff(times[idx_apos])]

# ╔═╡ 592a18e9-ee9a-4638-9454-f0bda3a0a3f2
period = LilGuys.mean(periods)

# ╔═╡ 179c3c32-1368-4a58-b4b8-26d9d3f19f8c
md"""
## Identification of interesting locations to search....
- for this project, want to find stars between 2 and 10 degrees away along stream path.
"""

# ╔═╡ 5704daca-a8c4-4292-a6c0-ea294f4373fd
md"""
## Saving
"""

# ╔═╡ 76d5bed6-f1ba-4a2d-8425-ceb40d18abdc
let
	
	orbital_properties = Dict(
		"pericentre" => r_peri,
		"apocentre" => r[idx_apos[idx_apos .< idx_f]],
		"period" => period,
		"idx_f" => idx_f,
		"distance_f" => obs_c.distance[idx_f],
		"idx_peri" => idx_peri,
		"t_last_peri" => t_f - t_peri,
		# these three are for the final stream coordinate frame
		"ra0" => ra0,
		"dec0" => dec0,
		"theta0" => θ0,
		"idx_peris" => idx_peris, 
		"idx_apos" => idx_apos,
		"t_peris" => t[idx_peris], 
		"t_apos" => t[idx_apos]
	)


	open(outfile, "w") do f
		TOML.print(f, orbital_properties)
	end

	println("saved properties to $outfile")
	orbital_properties
end

# ╔═╡ fcf93f45-f4a1-4bec-bf4f-b4e515bf5d67
write_fits(skyorbit_outfile, obs_c, overwrite=true)

# ╔═╡ Cell order:
# ╟─8b41af50-9ae0-475b-bacc-3799e2949b30
# ╠═d94663e8-b30e-4712-8a3e-6ef7f954f141
# ╟─27577252-53dc-415c-b9a1-82155ef9e4ca
# ╠═061b1886-1878-11ef-3806-b91643300982
# ╠═cd43d649-1d52-46a7-a621-68c8122fce6d
# ╠═882d4fc5-07ae-4b06-8da5-67f0894595db
# ╠═ab57edae-2292-4aef-9c1f-53802dbc0600
# ╠═7bf1a839-b24e-4478-bdd3-200989779f68
# ╠═6b3df4f3-de5b-40b1-8aec-1cfce9aa843a
# ╠═7834415a-b479-495b-998a-1d12f42f0dc6
# ╠═56d15948-4fd8-4541-9925-99837d9584f5
# ╟─643cd0bf-77b3-4201-9ff7-09dd5aee277c
# ╠═0e2e7f93-09d1-4902-ad24-223df50d37cb
# ╠═3953b76f-a726-4211-a6b8-5cf38149dcdf
# ╠═24ed3bc1-37c6-4601-b941-0780c53a9630
# ╠═079afa70-c7dc-4347-9ae3-8459ba2fa941
# ╠═94344455-d1d2-4ef9-af11-2d79ee4729ee
# ╟─dd56b7ec-be11-447f-acc1-12750d82879b
# ╠═bc6deac8-b70a-483b-9fd7-1413c6f17aa7
# ╠═69d83e00-7eb6-4271-838f-80e4d1654dac
# ╠═d142b7bd-3002-4331-a725-577873c42f28
# ╠═61d788ec-3518-4e3a-8eef-59c86ae5fc1a
# ╠═ac2c7484-9acd-4fda-9699-fdf17da507c2
# ╠═0dd476fd-be53-4e9b-a686-a4462485c64c
# ╠═2bc762ad-e590-443e-b3c2-91dc42a8a4d9
# ╠═bb9ef388-bb5a-45a3-836e-4c06dbe0ab65
# ╠═9c427388-c657-4bb7-bc0a-b4de3597c645
# ╟─30969f77-667e-4ae4-9897-82c1c1182652
# ╠═96a57df5-a7b7-447a-a4a6-2b05e391a5c6
# ╠═2c702eb7-ebb6-44c9-8e01-ca52d011c014
# ╠═b250bf10-c228-4b14-938a-35561ae871d7
# ╠═bb0cb8c2-2cbd-4205-a39e-4c4c0ff98b8a
# ╟─08c3df42-738b-47c4-aa6b-fc39a9cfc02f
# ╠═a1c992c6-ad12-4968-b105-adfa1f327e76
# ╠═5255c605-56ea-4eb3-bd20-5134e3a96705
# ╠═aa2c3a93-19a3-43d8-82de-ae6ed8c4b9f7
# ╠═15293cb8-61d3-478d-a2ae-5a5b2006db44
# ╠═ddd74bbb-df27-4150-b497-b1df5243f518
# ╠═f88b909f-c3dc-41e0-bdb1-25e229964d27
# ╠═6d0612d6-8609-4832-80ae-4e8e78c557cc
# ╟─7d29a3bd-dc83-4eb3-ae65-fce5270ed8d5
# ╠═f134b3ce-53f0-47e9-84e9-1e73064d5191
# ╠═64e558da-2928-4815-ad5a-7528516311f9
# ╠═defc4184-2613-4480-9adf-fa135f168382
# ╠═71ef4c9b-1284-4451-bcdd-eff79d334539
# ╠═5ec062c3-3815-4cf7-b45a-f97332d1b800
# ╠═a179323f-4878-4021-b8d4-69ca733658cb
# ╠═ecf7c820-81a4-4cb7-a794-b7835c77811e
# ╠═cdde517a-1b3e-4d96-9156-4a8f72b795e9
# ╠═cb601d4f-0fd6-4d1d-8bf0-910471c729c6
# ╠═319b905c-2d08-4a95-9d95-9cd26e2f5b1f
# ╠═9d60d54b-70e8-4b3c-a7c7-7caaa2f94a1c
# ╠═7646ea5b-f1b1-4934-be68-330139f7f838
# ╠═d9df3376-6ca1-4701-afb5-2df994bb3442
# ╠═9530e936-1225-4cfc-aa9a-bf7644d612f5
# ╠═7a30bd90-946e-418c-8339-be64c37cda76
# ╠═0c69519c-9650-46b9-89f9-cc37227f5b1a
# ╠═d95c457b-c9be-4570-bc90-b4bbb7de56e2
# ╠═88a9bb2d-9f86-4a96-a3b8-7a057480c7c3
# ╠═8d1508af-1715-4ef6-aab9-e95a02265913
# ╠═efcbae60-cf7c-4e74-aae4-39d19b74b6fa
# ╠═d81c1455-728a-4023-ad65-e3cce37a69f9
# ╠═a5ce5442-73ca-4aaf-915a-72fe9936e791
# ╠═7e14d1c7-b37f-4c0c-8e2a-1ac7cce27aaa
# ╠═f64dcd49-b0ca-4319-b615-5520b23d7818
# ╠═04d29fcb-70a0-414b-a487-7a18c44b9d58
# ╠═af8a50bd-e761-4439-9fc9-80048c264d5b
# ╠═73bb2d61-37f3-4782-ae89-d36d1ff8f9ff
# ╟─14eebce8-04f7-493b-824a-7808c7fa35dd
# ╠═9e2cbc6b-58ec-45d0-9afa-568a7bc8a33e
# ╠═54cf5233-a955-4831-86ad-23b72f15789d
# ╠═cdabdc7d-76a1-45f5-b83a-2454576d3964
# ╠═891f31ed-5565-4f43-ab82-b8bc3a76c1cf
# ╠═ca193997-61c3-4302-a200-4e7d4e777521
# ╠═225a8adb-82ee-4479-806e-15796e2b08e2
# ╠═293090a7-8ee0-442e-b70f-e6b7750ab319
# ╠═052fb31a-788a-485f-b8d2-1cf49f6ffb4b
# ╠═e96a9789-22b8-433a-8a3c-182d7ec4e82a
# ╠═5858e6d0-a501-43ee-9669-6a225a1df1c9
# ╠═a48b78df-588d-4f94-9252-badd27179deb
# ╟─3448ffc5-41e6-4208-b11f-2c00168bf50a
# ╠═66435478-0619-47aa-a659-c06089951f72
# ╠═78271e36-12b7-4edc-bfb0-20ecc597ab20
# ╠═661ca87c-c8da-49b1-b8a3-72c81050590b
# ╠═77aa1d73-0c90-4f6c-9383-99e9e9f0379a
# ╠═2318b289-d6c4-44bf-b2ff-36f448faf97b
# ╠═afdb058a-ebbd-4f07-b0d8-a85bb1070737
# ╠═31a11704-1dad-4007-b704-9312b81a5bad
# ╠═bc1a33f4-cf2e-44a6-a20b-1344f47e75c6
# ╠═2d1017a9-03a0-4aa2-829f-72174eaaa363
# ╠═63d2d908-2f12-483b-bbad-833b2aecc4e3
# ╠═5676b72f-a981-4efb-8328-c25d3c5d6fb0
# ╠═592a18e9-ee9a-4638-9454-f0bda3a0a3f2
# ╟─179c3c32-1368-4a58-b4b8-26d9d3f19f8c
# ╟─5704daca-a8c4-4292-a6c0-ea294f4373fd
# ╠═76d5bed6-f1ba-4a2d-8425-ceb40d18abdc
# ╠═fcf93f45-f4a1-4bec-bf4f-b4e515bf5d67
