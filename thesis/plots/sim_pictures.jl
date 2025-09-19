### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ b0a0dddc-fb5a-4bb2-b048-a54859d0b703
using DataFrames, CSV

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 913e0316-a04b-4270-ba31-0ba0f7fdd705
galaxyname = "ursa_minor"

# ╔═╡ bf49209c-fbfc-4439-a7d8-cfad5ceba8cc
import TOML

# ╔═╡ 4e5ddbe9-e90a-42cf-a0f1-fabb3d3f1675
xy_iso = CSV.read("resources/EP2020_iso_xy.csv", DataFrame, tasks=1)

# ╔═╡ 53cdcc20-96d4-4bd5-8028-df9a38af71ae
xz_iso = CSV.read("resources/EP2020_iso_xz.csv", DataFrame, tasks=1)

# ╔═╡ 8d7ba7af-e064-45a5-ba6a-4b12deafa693


# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ b9600d93-5946-4380-a2ae-9b5f673bbaf5
modelname = if galaxyname == "sculptor"
	"sculptor/1e7_new_v31_r3.2/orbit_smallperi"
elseif galaxyname == "ursa_minor"
	"ursa_minor/1e7_new_v38_r4.0/orbit_smallperi.5"
end

# ╔═╡ 0f71807d-d698-4164-9f30-49af8dd8ba55
point_orbit = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "simulations", modelname, "orbit.toml"))

# ╔═╡ b593188f-1618-44d4-b38c-d5be448e927c
pos_final = LilGuys.position(LilGuys.transform(Galactocentric, ICRS(point_orbit)))

# ╔═╡ 32db23d9-7959-41ac-aff4-b63df5e4b94a
figname = Dict(
	"sculptor" => "scl",
	"ursa_minor" => "umi"
)[galaxyname]

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
out = Output(joinpath(ENV["DWARFS_ROOT"], "analysis", modelname));

# ╔═╡ a3be2d61-98eb-4037-afb4-4155ba24cc21
orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", modelname, "orbital_properties.toml"))

# ╔═╡ d50c120c-d72a-4fdc-b18a-9ec185c25c04
orbit = Orbit(joinpath(ENV["DWARFS_ROOT"], "analysis", modelname, "centres.hdf5"))

# ╔═╡ 5a40b893-021b-46e5-a115-0284e13ae7ae
bins = LinRange(-150, 150, 512)

# ╔═╡ f9976cbe-ee5c-4c4b-95d7-94302bfbf7aa
function get_xy(snap)
	return snap.positions[2, :] .- 0snap.x_cen[2], snap.positions[3, :] .- 0snap.x_cen[3]
end

# ╔═╡ dac0e03d-7fee-455d-9fd6-1e2a03ce65da
function max_density(out)
	return maximum(Arya.histogram2d(get_xy(out[1])..., bins).values)
end

# ╔═╡ a9e6cdcc-d31b-404a-bb32-4ceb83b90efa
colorrange=colorrange=(1e-5, 1.) .* max_density(out)

# ╔═╡ 666eff06-c57a-4219-ac5e-e68e6b860882
function plot_xy_density!(snap)
	Arya.hist2d!(get_xy(snap)...,  bins = bins, colorscale=log10, colorrange=colorrange)
end


# ╔═╡ 4ae004da-7e17-4f7c-a11d-6ad97dcbe6dd
orbit_props["idx_apos"], orbit_props["idx_peris"]

# ╔═╡ d57c2310-706f-4a0f-af6e-f4dae0cd8376
orbit_props["idx_f"]

# ╔═╡ 962d5024-6afa-452c-bb04-85046971ee2a
idxs = if galaxyname == "sculptor"
	[1, 87, 171, 254, 338, orbit_props["idx_f"]]
elseif galaxyname == "ursa_minor"
	[1, 28, 53, 103, 152, 206]
end

# ╔═╡ 971bb9b8-6807-40f7-ad3e-c9e62ab4f3f1
function plot_today!()
	scatter!(pos_final[2], pos_final[3], color=COLORS[3], markersize=theme(:markersize)[]/2)
end

# ╔═╡ 800caee7-deee-4884-92f8-54ce6bb439f2
idx_apos = [1; orbit_props["idx_apos"]]

# ╔═╡ 2e0655e6-1351-4495-b56c-bc70bc478d01
@assert idx_apos |> issorted

# ╔═╡ d5591e1b-ec3d-46fc-856e-ce8d1cf6bb03
function plot_orbit_trace!(orbit, idx; d_idx=5)
	if idx-d_idx < idx_apos[1]
		idx_last = 1
	else
		idx_last = filter(x->x .<= idx - d_idx, idx_apos)[end]
	end
	if 	idx + d_idx > idx_apos[end]
		idx_next = length(orbit)
	else
		idx_next = filter(x->x .>= idx .+ d_idx, idx_apos)[1]
	end

	x = orbit.positions[2, idx_last:idx]
	y = orbit.positions[3, idx_last:idx]
	lines!(x, y, linewidth=theme(:linewidth)[]/2, color=(:white, 0.3))

	x = orbit.positions[2, idx:idx_next]
	y = orbit.positions[3, idx:idx_next]
	lines!(x, y, linestyle=:dot, linewidth=theme(:linewidth)[]/2, color=(:white, 0.3))

end

# ╔═╡ c689f2b1-eded-4284-9b40-67c72016de8c
function plot_scalebar!(scale_length=50)
	length_relative = scale_length / (bins[end] - bins[1])

	x0, y0 = 0.05, 0.05
	lines!([x0, x0 + length_relative], [y0, y0], color=:white, space=:relative, linewidth=theme(:linewidth)[] / 2)
	text!(x0, y0, text="$scale_length kpc", color=:white, space=:relative, fontsize=0.8 * theme(:fontsize)[], )
end

# ╔═╡ 5e48acfc-085f-43c8-a4ae-d0aaabcae4ee
let
	fig = Figure()

	for (i, idx) in enumerate(idxs)
		ii, jj = (i-1) ÷ 3 , (i-1) % 3 
		axis = Axis(fig[ii, jj])
		plot_xy_density!(out[idx])
		poly!(xz_iso.x, xz_iso.z, color=COLORS[9])
		if i == 1
			plot_scalebar!()
		end

		plot_today!()
		plot_orbit_trace!(orbit, idx)

		hidexdecorations!()
		hideydecorations!()
		t = (out.times[idx] - out.times[orbit_props["idx_f"]]) * T2GYR
		text!(0.05, 0.95, text= "$(round(t, digits=1)) Gyr", space=:relative, color=:white, align=(:left, :top))
	end

	rowsize!(fig.layout, 1, Aspect(1, 1.0))
	rowsize!(fig.layout, 0, Aspect(1, 1.0))

	rowgap!(fig.layout, 6)
	colgap!(fig.layout, 6)

	@savefig "$(figname)_sim_images"
	fig
end


# ╔═╡ Cell order:
# ╠═913e0316-a04b-4270-ba31-0ba0f7fdd705
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═b0a0dddc-fb5a-4bb2-b048-a54859d0b703
# ╠═bf49209c-fbfc-4439-a7d8-cfad5ceba8cc
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═4e5ddbe9-e90a-42cf-a0f1-fabb3d3f1675
# ╠═53cdcc20-96d4-4bd5-8028-df9a38af71ae
# ╠═0f71807d-d698-4164-9f30-49af8dd8ba55
# ╠═8d7ba7af-e064-45a5-ba6a-4b12deafa693
# ╠═b593188f-1618-44d4-b38c-d5be448e927c
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═b9600d93-5946-4380-a2ae-9b5f673bbaf5
# ╠═32db23d9-7959-41ac-aff4-b63df5e4b94a
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═a3be2d61-98eb-4037-afb4-4155ba24cc21
# ╠═d50c120c-d72a-4fdc-b18a-9ec185c25c04
# ╠═5a40b893-021b-46e5-a115-0284e13ae7ae
# ╠═f9976cbe-ee5c-4c4b-95d7-94302bfbf7aa
# ╠═dac0e03d-7fee-455d-9fd6-1e2a03ce65da
# ╠═666eff06-c57a-4219-ac5e-e68e6b860882
# ╠═a9e6cdcc-d31b-404a-bb32-4ceb83b90efa
# ╠═4ae004da-7e17-4f7c-a11d-6ad97dcbe6dd
# ╠═d57c2310-706f-4a0f-af6e-f4dae0cd8376
# ╠═962d5024-6afa-452c-bb04-85046971ee2a
# ╠═971bb9b8-6807-40f7-ad3e-c9e62ab4f3f1
# ╠═800caee7-deee-4884-92f8-54ce6bb439f2
# ╠═2e0655e6-1351-4495-b56c-bc70bc478d01
# ╠═d5591e1b-ec3d-46fc-856e-ce8d1cf6bb03
# ╠═c689f2b1-eded-4284-9b40-67c72016de8c
# ╠═5e48acfc-085f-43c8-a4ae-d0aaabcae4ee
