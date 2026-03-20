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
# ╠═6d4116f5-7114-4c2a-bad3-08e6f7298c33
# ╠═f6208f80-16e0-4923-b18c-a7e09e16d5dd
# ╠═13403957-8511-4009-a648-3113251e8bb3
# ╠═1963059f-1f91-4761-9c26-295d744f6cdc
# ╠═c65ae398-47ec-45dd-83cd-663b952cc056
# ╠═0b5e25fc-780f-4046-b6ba-957d3bfbada9
