### A Pluto.jl notebook ###
# v0.20.13

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

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 913e0316-a04b-4270-ba31-0ba0f7fdd705
galaxyname = "ursa_minor"

# ╔═╡ bf49209c-fbfc-4439-a7d8-cfad5ceba8cc
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ b9600d93-5946-4380-a2ae-9b5f673bbaf5
modelname = if galaxyname == "sculptor"
	"sculptor/1e6_V31_r3.2/orbit_smallperi"
elseif galaxyname == "ursa_minor"
	"ursa_minor/1e7_new_v38_r4.0/orbit_smallperi.5"
end

# ╔═╡ 32db23d9-7959-41ac-aff4-b63df5e4b94a
figname = Dict(
	"sculptor" => "scl",
	"ursa_minor" => "umi"
)[galaxyname]

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
out = Output(joinpath(ENV["DWARFS_ROOT"], "analysis", modelname))

# ╔═╡ a3be2d61-98eb-4037-afb4-4155ba24cc21
orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", modelname, "orbital_properties.toml"))

# ╔═╡ 5a40b893-021b-46e5-a115-0284e13ae7ae
bins = LinRange(-150, 150, 100)

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

# ╔═╡ 962d5024-6afa-452c-bb04-85046971ee2a
idxs = [1, 40, 78, 115, 152, orbit_props["idx_f"]]

# ╔═╡ 5e48acfc-085f-43c8-a4ae-d0aaabcae4ee
let
	fig = Figure()

	for (i, idx) in enumerate(idxs)
		ii, jj = (i-1) ÷ 3 , (i-1) % 3 
		axis = Axis(fig[ii, jj])
		plot_xy_density!(out[idx])

		hidexdecorations!()
		hideydecorations!()
		t = (out.times[idx] - out.times[orbit_props["idx_f"]]) * T2GYR
		text!(0.05, 0.95, text= "$(round(t, digits=1)) Gyr", space=:relative, color=:white, align=(:left, :top))
	end

	rowsize!(fig.layout, 1, Aspect(1, 1.0))
	rowsize!(fig.layout, 0, Aspect(1, 1.0))

	rowgap!(fig.layout, 12)
	colgap!(fig.layout, 12)

	@savefig "$(figname)_sim_images"
	fig
end


# ╔═╡ Cell order:
# ╠═913e0316-a04b-4270-ba31-0ba0f7fdd705
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═bf49209c-fbfc-4439-a7d8-cfad5ceba8cc
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═b9600d93-5946-4380-a2ae-9b5f673bbaf5
# ╠═32db23d9-7959-41ac-aff4-b63df5e4b94a
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═a3be2d61-98eb-4037-afb4-4155ba24cc21
# ╠═5a40b893-021b-46e5-a115-0284e13ae7ae
# ╠═f9976cbe-ee5c-4c4b-95d7-94302bfbf7aa
# ╠═dac0e03d-7fee-455d-9fd6-1e2a03ce65da
# ╠═666eff06-c57a-4219-ac5e-e68e6b860882
# ╠═a9e6cdcc-d31b-404a-bb32-4ceb83b90efa
# ╠═4ae004da-7e17-4f7c-a11d-6ad97dcbe6dd
# ╠═962d5024-6afa-452c-bb04-85046971ee2a
# ╠═5e48acfc-085f-43c8-a4ae-d0aaabcae4ee
