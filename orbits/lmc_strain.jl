### A Pluto.jl notebook ###
# v0.20.13

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

# ╔═╡ a9f69692-73c1-11f0-0f32-07a13b194774
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ ec032107-0cd7-4083-89ab-c5df1669cc08
using PlutoUI

# ╔═╡ 5b22ede7-b654-49d6-aefb-d42e532ad618
using Agama

# ╔═╡ 1b94b9cc-9392-4951-a02e-c2939882e5e3
using CSV, DataFrames

# ╔═╡ 9b5f7216-1aec-49a5-9a6e-506ab99b1455
@bind modelname TextField(24, default="vasiliev24_L3M11") |> confirm

# ╔═╡ 48507e9a-d114-4319-9aa7-a5505295f879
@bind galaxyname TextField(18, default="sculptor") |> confirm

# ╔═╡ c8d4b259-94b2-4709-88d6-946ec9be94ec
Nmax=50

# ╔═╡ 8cb946ac-e82d-4e0e-924a-cf7a5eb4dd2b
module OrbitUtils
	include("orbit_utils.jl")
end

# ╔═╡ 15d4d472-ca0c-4f34-bb0d-fcbc2130e526
CairoMakie.activate!(type=:png)

# ╔═╡ 5f143abb-7316-4901-9b3f-f6bb782d92ab
md"""
# Data Loading
"""

# ╔═╡ 1894ed6a-81d0-42ad-95bb-748cbabb6c7d
modeldir = "$galaxyname/$modelname"

# ╔═╡ 5234002f-c690-463b-a4dc-578cbd88cbfa
pot = OrbitUtils.get_potential(modeldir)

# ╔═╡ 00ad8f61-de82-45c8-8ff6-5b46943cd5e4
pot_lmc = Agama.Potential(file=joinpath(galaxyname, modelname, "potential_lmc.ini"    ))

# ╔═╡ 698f4940-0e45-4282-9ebd-52089b462ea4
pot_lmc_static = Agama.Potential(file=joinpath(galaxyname, modelname, "potential_lmc_static.ini"    ))

# ╔═╡ 91ffe8b5-48ba-4390-9a5c-30d34f4c2af2
units = OrbitUtils.get_units(modeldir)

# ╔═╡ c6a4738c-d515-4557-9cbf-93f3a3d6881c
idx, orbits = let
	structs = LilGuys.read_ordered_structs(joinpath(galaxyname, modelname, "orbits.hdf5"), LilGuys.Orbit)

	filt = 1:min(Nmax, length(structs))
	first.(structs)[filt], last.(structs)[filt]
end

# ╔═╡ 5a01d531-a7d7-4098-b965-857f5c8b5ff2
lmc_orbit = LilGuys.resample(OrbitUtils.get_lmc_orbit(modeldir), orbits[1].times)

# ╔═╡ cc04ec75-2dd5-458b-86a3-28e2fb49737a
md"""
# calculations
"""

# ╔═╡ ce2b607c-d5ee-435b-a851-fb0945da0867
pos = LilGuys.positions.(orbits)

# ╔═╡ 296f7eb8-f0d4-4a8c-9935-5cac11cf23e7
pos_lmc = LilGuys.positions.(orbits .- [lmc_orbit])

# ╔═╡ 0465b522-615d-4a42-8711-a032424785ad
T_max = OrbitUtils.scalar_tidal_forces.(pot, orbits, agama_units=units)

# ╔═╡ 5fe81374-16a9-4524-b0c8-544f2b36dbd8
T_lmc_max =  OrbitUtils.scalar_tidal_forces.(pot_lmc, orbits, agama_units=units)

# ╔═╡ d53303dd-7a5e-4876-a11e-66da0ac609aa
times = orbits[1].times

# ╔═╡ 4a729d8a-43fe-492c-bd69-2315038bcea5
T = [Agama.stress(pot, x, units, t=times) for x in pos]

# ╔═╡ 67f2069c-fc5b-45e9-a0f8-b5966c128695
T_lmc = [Agama.stress(pot_lmc, x, units, t=times) for x in pos]

# ╔═╡ efc9b19d-1cbb-4d83-af47-c93b9846cc69
T_lmc_static = [Agama.stress(pot_lmc_static, x, units, t=times) for x in pos]

# ╔═╡ 67dab4cf-8d6b-4dcf-9c07-a8f7542232a4
boundmass = CSV.read("sculptor/vasiliev24_L3M11/boundmass.txt", DataFrame, skipto=2, header=["time", "boundmass"], delim=" ", ignorerepeated=true)

# ╔═╡ 746beaa4-5b3b-4b03-83f4-66489be2ef93
lines(boundmass.time, log10.(boundmass.boundmass), axis=(;xlabel="time / Gyr", ylabel = "log10 boundmass / solar masses"))

# ╔═╡ 2514a322-da6e-4f67-9f25-cef3d357dc85
bm = LilGuys.lerp(boundmass.time ./ T2GYR, boundmass.boundmass)

# ╔═╡ 0dbc84f7-b3f5-4a7c-87bc-988ad886f0fd
md"""
# Plots
"""

# ╔═╡ 1f314edf-3cc0-4b7a-90f7-372c14fe5ef2
T_lmc_norm = [t ./ bm.(times)' .* 1e11 for t in T_lmc]

# ╔═╡ a14ca353-f71a-412a-83b1-a0b2caba796c
tides_lmc = OrbitUtils.scalar_tidal_forces.(pot_lmc, pos, agama_units=units)                 

# ╔═╡ 36dae006-f1ab-4c22-9f01-104c1a423457
if modelname == "vasiliev24_L3M11"                                                
	modelname_no = "vasiliev24_M11"                                               
end

# ╔═╡ 24d0c057-14fc-4080-a8b6-34be134e5435
pot_no = Agama.Potential(file=joinpath(galaxyname, modelname_no, "agama_potential.ini"))

# ╔═╡ 282ea824-06b3-4fe2-83c6-175f3fb02952
T_no = [Agama.stress(pot_no, x, units, t=times) for x in pos]

# ╔═╡ 68e074d5-8f14-46f5-8a03-f7a7575a8a2a
let
	fig = Figure(size=(6*72, 5*72))

	# diagonal terms (row 1)
	ax_x = Axis(fig[1,1], title=L"x ($i=1$)", ylabel = L"strain $T_{i, i}$")

	for i in eachindex(orbits)
		lines!(times * T2GYR, T[i][1, :], color=COLORS[1], alpha=0.1)
		lines!(times * T2GYR, T_lmc[i][1, :], color=COLORS[2], alpha=0.1)
		lines!(times * T2GYR, T_no[i][1, :], color=COLORS[3], alpha=0.1)
		lines!(times * T2GYR, T_lmc_static[i][1, :], color=COLORS[4], alpha=0.1)
	end


	ax_y = Axis(fig[1,2], title=L"y ($i=2$)")
	for i in eachindex(orbits)
		lines!(times * T2GYR, T[i][2, :], color=COLORS[1], alpha=0.1)
		lines!(times * T2GYR, T_lmc[i][2, :], color=COLORS[2], alpha=0.1)
		lines!(times * T2GYR, T_no[i][2, :], color=COLORS[3], alpha=0.1)
		lines!(times * T2GYR, T_lmc_static[i][2, :], color=COLORS[4], alpha=0.1)

	end

	ax_z = Axis(fig[1,3], title=L"z ($i=3$)")
	for i in eachindex(orbits)
		lines!(times * T2GYR, T[i][3, :], color=COLORS[1], alpha=0.1)
		lines!(times * T2GYR, T_lmc[i][3, :], color=COLORS[2], alpha=0.1)
		lines!(times * T2GYR, T_no[i][3, :], color=COLORS[3], alpha=0.1)
		lines!(times * T2GYR, T_lmc_static[i][3, :], color=COLORS[4], alpha=0.1)


	end


	# cross terms (row 2)

	ax_xy = Axis(fig[2,1], ylabel=L"strain $T_{i, i+1}$")
	for i in eachindex(orbits)
		lines!(times * T2GYR, T[i][4, :], color=COLORS[1], alpha=0.1)
		lines!(times * T2GYR, T_lmc[i][4, :], color=COLORS[2], alpha=0.1)
		lines!(times * T2GYR, T_no[i][4, :], color=COLORS[3], alpha=0.1)
	end

	# create legend
	lines!([NaN], [NaN], color=COLORS[1], label="total")
	lines!([NaN], [NaN], color=COLORS[2], label="lmc only")
	lines!([NaN], [NaN], color=COLORS[3], label="mw only")
	lines!([NaN], [NaN], color=COLORS[4], label="static lmc")
	axislegend(position=:lb)
	

	ax_yz = Axis(fig[2,2])
	for i in eachindex(orbits)
		lines!(times * T2GYR, T[i][5, :], color=COLORS[1], alpha=0.1)
		lines!(times * T2GYR, T_lmc[i][5, :], color=COLORS[2], alpha=0.1)
		lines!(times * T2GYR, T_no[i][5, :], color=COLORS[3], alpha=0.1)
	end

	ax_zx = Axis(fig[2,3])
	for i in eachindex(orbits)
		lines!(times * T2GYR, T[i][6, :], color=COLORS[1], alpha=0.1)
		lines!(times * T2GYR, T_lmc[i][6, :], color=COLORS[2], alpha=0.1)		
		lines!(times * T2GYR, T_no[i][6, :], color=COLORS[3], alpha=0.1)
	end


	# positions (row 3)

	ax_xx = Axis(fig[3,1], ylabel="(relative) position")

	for i in eachindex(orbits)
		lines!(times * T2GYR, pos_lmc[i][1, :], color=COLORS[2], alpha=0.1)
		lines!(times * T2GYR, pos[i][1, :], color=COLORS[3], alpha=0.1)
	end

	

	ax_yy = Axis(fig[3,2], xlabel="time / Gyr")
	for i in eachindex(orbits)
		lines!(times * T2GYR, pos_lmc[i][2, :], color=COLORS[2], alpha=0.1)
		lines!(times * T2GYR, pos[i][2, :], color=COLORS[3], alpha=0.1)
	end

	ax_zz = Axis(fig[3,3])
	for i in eachindex(orbits)
		lines!(times * T2GYR, pos_lmc[i][3, :], color=COLORS[2], alpha=0.1)
		lines!(times * T2GYR, pos[i][3, :], color=COLORS[3], alpha=0.1)
	end


	# formatting
	linkxaxes!(ax_x, ax_y, ax_z, ax_xy, ax_yz, ax_zx, ax_xx, ax_yy, ax_zz)
	linkyaxes!(ax_x, ax_y, ax_z, ax_xy, ax_yz, ax_zx)
	linkyaxes!(ax_xx, ax_yy, ax_zz)

	for ax in [ax_y, ax_z, ax_yz, ax_zx, ax_yy, ax_zz]
		hideydecorations!(ax, ticks=false)
	end

	for ax in [ax_x, ax_y, ax_z, ax_xy, ax_yz, ax_zx]
		hidexdecorations!(ax, ticks=false)
	end

	
	xlims!(-1, 0)
	
	fig
end

# ╔═╡ d8348db6-638c-4897-b227-f6bc7247432a
T_no_max =  OrbitUtils.scalar_tidal_forces.(pot_no, orbits, agama_units=units)

# ╔═╡ 6bb5a6ca-d1d3-4f57-be67-ecd39861adb0
_, orbits_no = let
	structs = LilGuys.read_ordered_structs(joinpath(galaxyname, modelname_no, "orbits.hdf5"), LilGuys.Orbit)

	filt = 1:min(Nmax, length(structs))
	first.(structs)[filt], last.(structs)[filt]
end

# ╔═╡ 159bb1d3-f44b-463d-89a4-59b036ecd7b9
T_orbit_no_max =  OrbitUtils.scalar_tidal_forces.(pot_no, orbits_no, agama_units=units)

# ╔═╡ 59baf6d9-3acb-43ea-9b44-05743932efc5
let
	fig = Figure(size=(6*72, 3*72))
	ax_xy = Axis(fig[1,1], limits=(-300, 300, -300, 300), aspect=DataAspect())
	ax_yz = Axis(fig[1,2], limits=(-300, 300, -300, 300), aspect=DataAspect())
	ax_xz = Axis(fig[1,3], limits=(-300, 300, -300, 300), aspect=DataAspect())


	i = 1
	x_lmc = Observable(lmc_orbit.positions[1, i])
	y_lmc = Observable(lmc_orbit.positions[2, i]) 
	z_lmc = Observable(lmc_orbit.positions[3, i])
	
	scatter!(ax_xy, x_lmc, y_lmc)
	scatter!(ax_yz, y_lmc, z_lmc)
	scatter!(ax_xz, x_lmc, z_lmc)

	x = Observable([p[1, i] for p in pos])
	y = Observable([p[2, i] for p in pos])
	z = Observable([p[3, i] for p in pos])
	
	scatter!(ax_xy, x, y, alpha=0.1)
	scatter!(ax_yz, y, z, alpha=0.1)
	scatter!(ax_xz, x, z, alpha=0.1)
	
	scatter!(ax_xy, 0, 0, color=:black)
	scatter!(ax_yz, 0, 0, color=:black)
	scatter!(ax_xz, 0, 0, color=:black)
	t = Observable("0.0 Gyr")

	text!(0.1, 0.9, text=t, space=:relative, align=(:left, :top))


	CairoMakie.Makie.Record(fig, eachindex(lmc_orbit.times)[end:-30:1]) do i
		#l[1][] = Point2f.([p[2, i] for p in pos], [p[3, i] for p in pos])
		x_lmc[] = lmc_orbit.positions[1, i]
		y_lmc[] = lmc_orbit.positions[2, i]
		z_lmc[] = lmc_orbit.positions[3, i]
		
		x[] = [p[1, i] for p in pos]
		y[] = [p[2, i] for p in pos]
		z[] = [p[3, i] for p in pos]
		t[] = "$(round(lmc_orbit.times[i] * T2GYR, digits=2)), Gyr"

	end
end

# ╔═╡ 4fede773-ce74-48a0-94d8-498e3afb765d
let
	fig = Figure()
	ax_x = Axis(fig[1,1], xlabel = "time / Gyr", ylabel = "tidal strain")

	for i in eachindex(orbits)
		lines!(times * T2GYR, T_max[i], color=COLORS[1], alpha=0.1)
		lines!(times * T2GYR, T_lmc_max[i], color=COLORS[2], alpha=0.1)
		lines!(times * T2GYR, T_no_max[i], color=COLORS[3], alpha=0.1)
		lines!(times * T2GYR, T_orbit_no_max[i], color=COLORS[4], alpha=0.1)
	end

	fig


	
	lines!([NaN], [NaN], color=COLORS[1], label="total")
	lines!([NaN], [NaN], color=COLORS[2], label="lmc only")
	lines!([NaN], [NaN], color=COLORS[3], label="mw only")
	lines!([NaN], [NaN], color=COLORS[4], label="no LMC")
	axislegend(position=:lt)


	fig
end

# ╔═╡ Cell order:
# ╠═9b5f7216-1aec-49a5-9a6e-506ab99b1455
# ╠═48507e9a-d114-4319-9aa7-a5505295f879
# ╠═c8d4b259-94b2-4709-88d6-946ec9be94ec
# ╠═a9f69692-73c1-11f0-0f32-07a13b194774
# ╠═ec032107-0cd7-4083-89ab-c5df1669cc08
# ╠═5b22ede7-b654-49d6-aefb-d42e532ad618
# ╠═8cb946ac-e82d-4e0e-924a-cf7a5eb4dd2b
# ╠═15d4d472-ca0c-4f34-bb0d-fcbc2130e526
# ╟─5f143abb-7316-4901-9b3f-f6bb782d92ab
# ╠═1894ed6a-81d0-42ad-95bb-748cbabb6c7d
# ╠═5234002f-c690-463b-a4dc-578cbd88cbfa
# ╠═00ad8f61-de82-45c8-8ff6-5b46943cd5e4
# ╠═698f4940-0e45-4282-9ebd-52089b462ea4
# ╠═24d0c057-14fc-4080-a8b6-34be134e5435
# ╠═91ffe8b5-48ba-4390-9a5c-30d34f4c2af2
# ╠═c6a4738c-d515-4557-9cbf-93f3a3d6881c
# ╠═6bb5a6ca-d1d3-4f57-be67-ecd39861adb0
# ╠═5a01d531-a7d7-4098-b965-857f5c8b5ff2
# ╠═1b94b9cc-9392-4951-a02e-c2939882e5e3
# ╠═cc04ec75-2dd5-458b-86a3-28e2fb49737a
# ╠═ce2b607c-d5ee-435b-a851-fb0945da0867
# ╠═296f7eb8-f0d4-4a8c-9935-5cac11cf23e7
# ╠═4a729d8a-43fe-492c-bd69-2315038bcea5
# ╠═0465b522-615d-4a42-8711-a032424785ad
# ╠═67f2069c-fc5b-45e9-a0f8-b5966c128695
# ╠═5fe81374-16a9-4524-b0c8-544f2b36dbd8
# ╠═efc9b19d-1cbb-4d83-af47-c93b9846cc69
# ╠═282ea824-06b3-4fe2-83c6-175f3fb02952
# ╠═d8348db6-638c-4897-b227-f6bc7247432a
# ╠═159bb1d3-f44b-463d-89a4-59b036ecd7b9
# ╠═d53303dd-7a5e-4876-a11e-66da0ac609aa
# ╠═67dab4cf-8d6b-4dcf-9c07-a8f7542232a4
# ╠═746beaa4-5b3b-4b03-83f4-66489be2ef93
# ╠═2514a322-da6e-4f67-9f25-cef3d357dc85
# ╟─0dbc84f7-b3f5-4a7c-87bc-988ad886f0fd
# ╠═1f314edf-3cc0-4b7a-90f7-372c14fe5ef2
# ╠═68e074d5-8f14-46f5-8a03-f7a7575a8a2a
# ╠═a14ca353-f71a-412a-83b1-a0b2caba796c
# ╠═36dae006-f1ab-4c22-9f01-104c1a423457
# ╠═59baf6d9-3acb-43ea-9b44-05743932efc5
# ╠═4fede773-ce74-48a0-94d8-498e3afb765d
