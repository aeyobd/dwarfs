### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ ef5f301c-23b7-11f1-befa-6f29a56e80d0
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya
	using CairoMakie
end

# ╔═╡ b24c494a-f632-49e6-ac64-58ff9002db5c
using PyFITS

# ╔═╡ bef4b575-9fb6-4007-9be2-88af0e10ddb6
using OrderedCollections

# ╔═╡ 0b36910e-a8ac-46c1-8a43-09c54aa3f885
CairoMakie.activate!(type=:png)

# ╔═╡ c6134a0e-12c9-40e6-8d88-0c6d019bf359


# ╔═╡ c46210a2-a330-463d-aecd-6f3df17c8785
function load_samples(orbitname)
	read_fits(joinpath(orbitname, "orbital_properties.fits"))
end

# ╔═╡ 42d0dbbb-2cfe-4bfe-a9f3-e13451fc0774
orbit_props = OrderedDict(
	"EP2020" => load_samples("EP2020"),
	"M11" => load_samples("vasiliev24_M11"),
	"L3M11" => load_samples("vasiliev24_L3M11"),
	"L2M10" => load_samples("vasiliev24_L2M10"),
)

# ╔═╡ e58bd8fe-1528-4f03-abe4-b06b4b1b73a8
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="pericentre / kpc")

	for (label, props) in orbit_props
		stephist!(props.pericentre, label=label,  bins=100, normalization=:pdf)
	end

	axislegend()
	fig
end

# ╔═╡ a1edc0f2-a30e-44ac-8ca1-0aaa0ba6668b
orbit_props["L2M10"].period

# ╔═╡ efc8c6c2-193d-469f-a397-61385277ea11
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="period / Gyr")

	for (label, props) in orbit_props
		stephist!(-T2GYR * props.period[.!ismissing.(props.period)], label=label,  bins=LinRange(0.5, 4, 50), normalization=:pdf)
	end

	axislegend(position=:rt)
	fig
end

# ╔═╡ 089546cf-f50d-44f8-9e17-25f67377b9fa
let
	fig = Figure(size=(5, 2) .* 72)
	ax = Axis(fig[1,1], xlabel="pericentre / kpc")

	hist!(orbit_props["EP2020"].pericentre)


	ax = Axis(fig[1,2], xlabel="time last peri / Gyr")

	hist!(orbit_props["EP2020"].time_last_peri * T2GYR)

	fig
end

# ╔═╡ b6971b90-09bb-4b65-bea5-ef407675d105
import TOML

# ╔═╡ 64280527-3c5a-4381-bdbb-313a9ca2c2a5
halo = NFW(v_circ_max = 22 / V2KMS, r_circ_max = 3.9)

# ╔═╡ c9943e95-2401-42e7-b35c-59ff4da72bb4
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/bootes3/observed_properties.toml"))

# ╔═╡ 84512d5b-f065-40d0-8618-0d2021ea56f6
Nsamples = length(orbit_props["EP2020"].time_last_peri)

# ╔═╡ c759f16a-09a9-4b18-b05b-b8f9992b0316
σv_samples = obs_props["sigma_v"] .+ randn(Nsamples) *  obs_props["sigma_v_ep"]

# ╔═╡ 294dc46c-df81-4215-82ed-668d6ad00bd1
r_break = LilGuys.break_radius.(orbit_props["EP2020"].time_last_peri, -σv_samples ./ V2KMS
								)

# ╔═╡ 0bdffe08-2514-4de6-a034-dec09b0e35eb
orbit_props["EP2020"].distance

# ╔═╡ 9bf97e31-de49-4cbc-9663-6f6ff12e9c8e
r_break_arcmin = LilGuys.kpc2arcmin.(r_break, orbit_props["EP2020"].distance)

# ╔═╡ e8356b0f-00f7-4384-8c50-c7f3fbc1f49b
hist(r_break)

# ╔═╡ 4edce737-e441-4e58-a6f9-8019e87a6bfd


# ╔═╡ cd2883db-7e60-4b1c-965e-b2e0ef76d586
function calc_ρ_mean(prof::LilGuys.SphericalProfile, r)
    y = LilGuys.mass.(prof, r)
    rho_mean = @. y / (4π/3*r^3)
end

# ╔═╡ e7a8024f-e1f0-4e41-b02e-d9acd3e73c4e
import Agama

# ╔═╡ 3025afd1-90e2-4030-b475-027541ac86ce
function calc_ρ_mean(prof::Agama.Potential, r; )
    M = Agama.enclosed_mass(prof, r)
    

    y = @. M / (4π/3 * r^3)
end

# ╔═╡ db4968d1-8763-4c5b-8660-d32a88394e7a
modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", "bootes3", "EP2020")

# ╔═╡ 871fc06b-ae34-4e1a-92f2-39d917e5e05e
mw_halo = Agama.Potential(file=joinpath(modeldir, "agama_potential.ini"))

# ╔═╡ 4061977b-3a03-448c-9804-d5acfebb8b0a
function r_J(peri)
	ρ_host = calc_ρ_mean(mw_halo, [peri, peri])[1]

	r_J = LilGuys.find_zero(r -> calc_ρ_mean(halo, r) - 3*ρ_host, 0.001)
end

# ╔═╡ b73d75cb-b503-4dfe-9cd2-c2d92f80e18d
r_jacobii = r_J.(orbit_props["EP2020"].pericentre[1:1000])

# ╔═╡ e64241f7-fd43-4514-9fdd-f0b993144f25
r_J_arcmin = LilGuys.kpc2arcmin.(r_jacobii, orbit_props["EP2020"].distance[1:1000])

# ╔═╡ 3527b54b-c30d-461c-8ba7-07b38bcd0aac
r_J(1)

# ╔═╡ 57d96424-54c9-4cd2-8ca9-8bc782713e82
let
	fig = Figure(size=(5, 2) .* 72)
	ax = Axis(fig[1,1], xlabel="break radius / arcmin")

	hist!(r_break_arcmin)

	ax = Axis(fig[1,2], xlabel="jacobi radius / arcmin")
	
	hist!(r_J_arcmin)

	fig
end

# ╔═╡ 6d4116f5-7114-4c2a-bad3-08e6f7298c33
let
	fig = Figure()
	ax = Axis(fig[1,1])

	for (label, props) in orbit_props
		scatter!(props.pericentre, props.time_last_peri .* T2GYR, label=label => (;markersize=5), markersize=1, )
	end
	axislegend()

	fig
end

# ╔═╡ f6208f80-16e0-4923-b18c-a7e09e16d5dd
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="pericentre / kpc", ylabel="period  / Gyr")

	for (label, props) in orbit_props
		scatter!(props.pericentre, -props.period .* T2GYR, label=label => (;markersize=5), markersize=1, )
	end
	axislegend(position=:lt)

	fig
end

# ╔═╡ 13403957-8511-4009-a648-3113251e8bb3
Nmax = 100

# ╔═╡ 1963059f-1f91-4761-9c26-295d744f6cdc
function read_orbits(modeldir)
	structs = LilGuys.read_ordered_structs(joinpath(modeldir, "orbits.hdf5"), LilGuys.Orbit)

	filt = 1:min(Nmax, length(structs))
	
	last.(structs)[filt]
end

# ╔═╡ c65ae398-47ec-45dd-83cd-663b952cc056
orbits = OrderedDict(
	"EP2020" => read_orbits("EP2020"),
	"M11" => read_orbits("vasiliev24_M11"),
	"L3M11" => read_orbits("vasiliev24_L3M11"),
	"L2M10" => read_orbits("vasiliev24_L2M10"),
)

# ╔═╡ 0b5e25fc-780f-4046-b6ba-957d3bfbada9
let
	fig = Figure()


	for (i, (label, orbit)) in enumerate(orbits)
		ax_rt = Axis(fig[i,1],
			xlabel = "time / Gyr",
			ylabel = L"$r_\textrm{Galcen}$ / kpc",
			limits=(-10, 0.1, 0, 150),
					   )	
		color = COLORS[i]

		for o in orbit
			lines!(o.times * T2GYR, radii(o), label=label, color=color, alpha=0.1, linewidth=1)
		end

		if i < length(orbits)
			hidexdecorations!(ticks=false, minorticks=false)
		end
	end

	fig
end

# ╔═╡ Cell order:
# ╠═ef5f301c-23b7-11f1-befa-6f29a56e80d0
# ╠═b24c494a-f632-49e6-ac64-58ff9002db5c
# ╠═0b36910e-a8ac-46c1-8a43-09c54aa3f885
# ╠═c6134a0e-12c9-40e6-8d88-0c6d019bf359
# ╠═c46210a2-a330-463d-aecd-6f3df17c8785
# ╠═bef4b575-9fb6-4007-9be2-88af0e10ddb6
# ╠═42d0dbbb-2cfe-4bfe-a9f3-e13451fc0774
# ╠═e58bd8fe-1528-4f03-abe4-b06b4b1b73a8
# ╠═a1edc0f2-a30e-44ac-8ca1-0aaa0ba6668b
# ╠═efc8c6c2-193d-469f-a397-61385277ea11
# ╠═089546cf-f50d-44f8-9e17-25f67377b9fa
# ╠═b6971b90-09bb-4b65-bea5-ef407675d105
# ╠═64280527-3c5a-4381-bdbb-313a9ca2c2a5
# ╠═c9943e95-2401-42e7-b35c-59ff4da72bb4
# ╠═c759f16a-09a9-4b18-b05b-b8f9992b0316
# ╠═84512d5b-f065-40d0-8618-0d2021ea56f6
# ╠═294dc46c-df81-4215-82ed-668d6ad00bd1
# ╠═0bdffe08-2514-4de6-a034-dec09b0e35eb
# ╠═9bf97e31-de49-4cbc-9663-6f6ff12e9c8e
# ╠═e64241f7-fd43-4514-9fdd-f0b993144f25
# ╠═e8356b0f-00f7-4384-8c50-c7f3fbc1f49b
# ╠═4edce737-e441-4e58-a6f9-8019e87a6bfd
# ╠═b73d75cb-b503-4dfe-9cd2-c2d92f80e18d
# ╠═3527b54b-c30d-461c-8ba7-07b38bcd0aac
# ╠═4061977b-3a03-448c-9804-d5acfebb8b0a
# ╠═cd2883db-7e60-4b1c-965e-b2e0ef76d586
# ╠═3025afd1-90e2-4030-b475-027541ac86ce
# ╠═e7a8024f-e1f0-4e41-b02e-d9acd3e73c4e
# ╠═db4968d1-8763-4c5b-8660-d32a88394e7a
# ╠═871fc06b-ae34-4e1a-92f2-39d917e5e05e
# ╠═57d96424-54c9-4cd2-8ca9-8bc782713e82
# ╠═6d4116f5-7114-4c2a-bad3-08e6f7298c33
# ╠═f6208f80-16e0-4923-b18c-a7e09e16d5dd
# ╠═13403957-8511-4009-a648-3113251e8bb3
# ╠═1963059f-1f91-4761-9c26-295d744f6cdc
# ╠═c65ae398-47ec-45dd-83cd-663b952cc056
# ╠═0b5e25fc-780f-4046-b6ba-957d3bfbada9
