### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 53f06974-1fcc-4c90-86a0-5dd0cec4e4b8
begin
	import Pkg; Pkg.activate()

	using Arya
	using GLMakie

	using DataFrames
	using Measurements
	import LilGuys as lguys
end

# ╔═╡ 3d32bef6-de64-4d1e-93a8-46c921c86011
using StatsBase: mean, std

# ╔═╡ 72eb810c-114e-11ef-30ea-fd46e923ea49
md"""
# Observations of Sculptor

A quick collection of past literature observations and measurments of the Sculptor DSph
"""

# ╔═╡ 498440b6-5572-4ddf-b4b4-0a96107b1417
begin
	value(x::Real) = x
	value(x::Measurement) = x.val

	err(x::Real) = 0
	err(x::Measurement) = x.err
end

# ╔═╡ 010b8c68-8ae1-4e18-b001-32c24c2e43e4
function obs_row(study; ra=NaN, dec=NaN, pm_ra=NaN, pm_dec=NaN, distance=NaN, radial_velocity=NaN, fe_h=NaN)
	return (study=study, ra=ra, dec=dec, pm_ra=pm_ra, pm_dec=pm_dec, distance=distance, radial_velocity=radial_velocity, fe_h=fe_h)
end

# ╔═╡ 722bd144-d047-4abf-b82f-f733134d3eb7
dm_to_d(dm) = 10 * 10^(dm / 5) / 1e3

# ╔═╡ 26ae0d94-698c-4e9d-bac6-3e91d0d197ab
md"""
# References

"""

# ╔═╡ 351731df-8d84-41ba-83e7-793898b9a148
md"""
## Compilations
- McChonnachie 2012: For sculptor distance: Pietrzy´nski 2008, Walker + 2009 for RV. 
"""

# ╔═╡ 0e1d8bc0-103a-4623-8546-aa7a32ee4504
md"""
## Kirby 2013
Metallicity, SFR ??
"""

# ╔═╡ cb3ef5d8-8231-42f1-b54f-42a6b08694c1
md"""
## Tolstoy+2023
VLT/FLAMES spectoscopic survey with Gaia DR3 proper motions (and parallax cuts)

fe_h=-1.82 ± 0.45,
"""

# ╔═╡ 1ecf4404-caa5-418d-9170-2d539df69275
tolstoy23 = obs_row("T+23",
	radial_velocity = 111.2 ± 0.25,
	pm_ra = 0.097 ± 0.006,
	pm_dec = -0.148 ± 0.004 # is proper?
)

# ╔═╡ 8f9fe40c-dd18-437f-ae9e-f0a50fba8155
md"""
## Pietrzy´nski + 2008
RR Lyrae stars in J and K to determine distance. Get a distance modulus of 19.67 pm 0.02 pm 0.12 (sys).
Also mention Rizzi's 2002 thesis (not easy to find) which measures similar DM from TRGB (19.64pm0.08) and HB (19.66 pm 0.15)
"""

# ╔═╡ 750dd472-05e4-4c4b-b347-4f0581ec6f55
pietryznski2008 = obs_row("P+08",
	distance=dm_to_d(19.67 ± 0.12)
)

# ╔═╡ 0ea649d4-c260-450f-a51c-7598da5bfe2e
md"""
## Walker + 2009
https://ui.adsabs.harvard.edu/abs/2009AJ....137.3109W/abstract

Use an algorithm with expectation maximum to determine memberships with RVs and metallicities from Mike. Also record a velocity dispersion of 9.2 \pm 1.1 km / s
"""

# ╔═╡ 6eed53c4-7798-4ced-ae58-59679f4cb380
walker2009 =  obs_row("W+09",
	radial_velocity=111.4 ± 0.1
)

# ╔═╡ b8055846-efa5-4bcc-a82e-e3b4cd40788d
md"""
## Martínez-Vázquez + 2015 
RR Lyrae distance modulus
"""

# ╔═╡ f7132d82-eb69-4e56-b61b-dd60e6cf8fdb
mv15 = obs_row("MV+15", 
	distance= dm_to_d(19.62±0.04)
)

# ╔═╡ 8a1832ae-6335-44c8-bc56-d8a05e3e301d
md"""
## Battaglia + 2022
references for basic properties: Battaglia+2008, Martinez-Vazquez 2015, Muñoz + 2018

With orbits in the perturbed (LMC) potential, find peries between44.3 to 51.1, and apos between 274.3 and 552.4m. The LMC is likely very important for the evolution of sculptor
"""

# ╔═╡ 30f6a515-76dc-47fc-abc0-a90e28f64d9f
battaglia2022 = obs_row("B+22",
	pm_ra=0.099±0.002,
	pm_dec = -0.159 ± 0.002
)

# ╔═╡ 7a5e54c2-6186-4975-9cd4-b971a982a8cf
md"""
## Muñoz + 2018

"""

# ╔═╡ 28fe217f-9bfb-450a-9e58-f552adb580d7
md"""
## McChonnachie & Venn 2020 a
"""

# ╔═╡ 5540a9a2-54b6-42cf-ad68-470cdc9fc667
md"""
# Comparisons
"""

# ╔═╡ 777ec196-e193-456d-8d44-cba200a366dd
obs = DataFrame([pietryznski2008, walker2009, battaglia2022, mv15, tolstoy23])

# ╔═╡ e7a27898-0707-43f2-86e6-58864ccfc419
@recipe(ErrScatter) do scene
    Attributes(
        color = theme(scene, :markercolor),
		marker = :circle
    )
end

# ╔═╡ ed73b1da-ed6a-43c9-b608-591451114642
function Makie.plot!(sc::ErrScatter)
	x = sc[1]
	y = sc[2]
	xerr = sc[3]
	yerr = sc[4]

	println(x)
	println(y)
	println(xerr)
	println(yerr)
	errorbars!(sc, x, y, yerr, color=sc.color)
	errorbars!(sc, x, y, xerr, direction=:x, color=sc.color)
	scatter!(sc, x, y, color=sc.color, marker=sc.marker)

	sc
end

# ╔═╡ d03034b9-da5d-47c1-9d3e-157ba9c695cf
function Makie.convert_arguments(::Type{<: ErrScatter}, x::AbstractArray, y::AbstractArray)
	return (value.(x), value.(y), err.(x), err.(y))
end

# ╔═╡ a1f7d351-4f25-475b-939f-8e03ba5a10a0
md"""
# Positions
"""

# ╔═╡ f41f6639-1523-47c7-96f8-82f1fdafb7a1
md"""
### Distances
- Tully et al. (2013): Updates earlier Tully work: For sculptor, distance = 80 kpc +- 0.08  Tip of RGB measurement
- 
- McChonnachie & Venn 2012

"""

# ╔═╡ 51d76b15-28d1-4871-a7dc-2a8a43802bff
let
	filt = value.(obs.distance) .!== NaN
	df = obs[filt, :]
	N = size(df, 1)
	x = collect(1:N)
	
	fig = Figure()

	ax = Axis(fig[1,1],
		xticks =(x, df.study), 
		xminorticksvisible=false,
		xticklabelrotation=-0π/6,
	ylabel="distance / kpc")

	tight_xticklabel_spacing!(ax)

	println(err.(x))
	errscatter!(x, df.distance)

	
	fig
end

# ╔═╡ 4314e298-f4d3-41d5-88e3-d9e60ba966b1
md"""

### Proper motions

"""

# ╔═╡ 5f28c6a9-0efd-4afd-9881-a3eaefc98f35
let
	fig = Figure()
 
	ax = Axis(fig[1,1], 
		xlabel=L"\mu_{\alpha*}", ylabel=L"\mu_\delta"
	)

	i = 0
	for row in eachrow(obs)
		if value.(row.pm_ra) !== NaN
			i += 1
			errscatter!(ax, [row.pm_ra], [row.pm_dec], color=Arya.COLORS[i], label=row.study)
		end
	end

	axislegend(ax)
	fig
end

# ╔═╡ c3f68bcf-7795-4bf6-95c6-8228d83ff79e
md"""
### Radial velocities

"""

# ╔═╡ 21bc8939-1e8e-420b-ae41-a465203aa2e3
let
	filt = value.(obs.radial_velocity) .!== NaN
	df = obs[filt, :]
	N = size(df, 1)
	x = collect(1:N)
	
	fig = Figure()

	ax = Axis(fig[1,1],
		xticks =(x, df.study), 
		xminorticksvisible=false,
	ylabel="radial velocity / km / s")

	println(err.(x))
	errscatter!(x, df.radial_velocity)

	fig
end

# ╔═╡ 2d5a9bcb-972a-4f00-86ab-b92b69921f11
md"""

### Density profiles

"""

# ╔═╡ 6421acc3-2689-40bc-8be5-3c24468910b6
md"""


### Masses & Magnitudes
"""

# ╔═╡ Cell order:
# ╟─72eb810c-114e-11ef-30ea-fd46e923ea49
# ╠═53f06974-1fcc-4c90-86a0-5dd0cec4e4b8
# ╠═498440b6-5572-4ddf-b4b4-0a96107b1417
# ╠═3d32bef6-de64-4d1e-93a8-46c921c86011
# ╠═010b8c68-8ae1-4e18-b001-32c24c2e43e4
# ╠═722bd144-d047-4abf-b82f-f733134d3eb7
# ╟─26ae0d94-698c-4e9d-bac6-3e91d0d197ab
# ╠═351731df-8d84-41ba-83e7-793898b9a148
# ╠═0e1d8bc0-103a-4623-8546-aa7a32ee4504
# ╟─cb3ef5d8-8231-42f1-b54f-42a6b08694c1
# ╠═1ecf4404-caa5-418d-9170-2d539df69275
# ╟─8f9fe40c-dd18-437f-ae9e-f0a50fba8155
# ╠═750dd472-05e4-4c4b-b347-4f0581ec6f55
# ╟─0ea649d4-c260-450f-a51c-7598da5bfe2e
# ╠═6eed53c4-7798-4ced-ae58-59679f4cb380
# ╟─b8055846-efa5-4bcc-a82e-e3b4cd40788d
# ╠═f7132d82-eb69-4e56-b61b-dd60e6cf8fdb
# ╟─8a1832ae-6335-44c8-bc56-d8a05e3e301d
# ╠═30f6a515-76dc-47fc-abc0-a90e28f64d9f
# ╠═7a5e54c2-6186-4975-9cd4-b971a982a8cf
# ╠═28fe217f-9bfb-450a-9e58-f552adb580d7
# ╟─5540a9a2-54b6-42cf-ad68-470cdc9fc667
# ╠═777ec196-e193-456d-8d44-cba200a366dd
# ╠═e7a27898-0707-43f2-86e6-58864ccfc419
# ╠═ed73b1da-ed6a-43c9-b608-591451114642
# ╠═d03034b9-da5d-47c1-9d3e-157ba9c695cf
# ╟─a1f7d351-4f25-475b-939f-8e03ba5a10a0
# ╟─f41f6639-1523-47c7-96f8-82f1fdafb7a1
# ╠═51d76b15-28d1-4871-a7dc-2a8a43802bff
# ╟─4314e298-f4d3-41d5-88e3-d9e60ba966b1
# ╠═5f28c6a9-0efd-4afd-9881-a3eaefc98f35
# ╟─c3f68bcf-7795-4bf6-95c6-8228d83ff79e
# ╠═21bc8939-1e8e-420b-ae41-a465203aa2e3
# ╟─2d5a9bcb-972a-4f00-86ab-b92b69921f11
# ╟─6421acc3-2689-40bc-8be5-3c24468910b6
