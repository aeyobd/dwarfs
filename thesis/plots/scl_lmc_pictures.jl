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

# ╔═╡ b0a0dddc-fb5a-4bb2-b048-a54859d0b703
using DataFrames, CSV

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 913e0316-a04b-4270-ba31-0ba0f7fdd705
galaxyname = "sculptor"

# ╔═╡ bf49209c-fbfc-4439-a7d8-cfad5ceba8cc
import TOML

# ╔═╡ 9c7a671e-917e-473a-855e-26291e216676
import Agama

# ╔═╡ 4e5ddbe9-e90a-42cf-a0f1-fabb3d3f1675
xy_iso = CSV.read("resources/L3M11_iso_xy.csv", DataFrame)

# ╔═╡ 53cdcc20-96d4-4bd5-8028-df9a38af71ae
xz_iso = CSV.read("resources/L3M11_iso_xz.csv", DataFrame)

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ b9600d93-5946-4380-a2ae-9b5f673bbaf5
modelname = if galaxyname == "sculptor"
	"sculptor/1e6_new_v31_r4.2/smallperilmc"
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

# ╔═╡ 185cfe93-b946-4e72-be6f-853205c49d22
module AgamaProjection 
	include(joinpath(ENV["HOME"], "LilGuys.jl/scripts", "project_potential.jl"))

end

# ╔═╡ 5f036156-64e5-4127-a80d-614f8f526559
module Animation 
	include(joinpath(ENV["HOME"], "LilGuys.jl/scripts", "animate_dm.jl"))

end

# ╔═╡ d7771353-51a7-4461-a8d8-17b18623f46c
ENV["LGUYS_SCRIPTS"]

# ╔═╡ 147e195d-99f1-4c33-a65d-2ccfd7fa9bca
t_scale = Agama.time_scale(Agama.VASILIEV_UNITS)

# ╔═╡ 4078aa8f-4f2f-4f09-9706-4600838d9960
projected_potential(pot, idx) = AgamaProjection.project_agama_potential(pot._py, (bins, bins), time=out.times[idx] / t_scale)

# ╔═╡ 336bc60f-9cdc-4ad8-ab87-255654350c61
transparant_cmap = to_colormap([(COLORS[3], ceil(i)) for i in LinRange(0, 1, 1000)])

# ╔═╡ 653f3273-6c1f-49f9-a48b-bd395318a551
transparant_scl = to_colormap([(COLORS[2], i) for i in LinRange(0, 1, 1000)])

# ╔═╡ df6c58bd-385a-4f04-a652-c5de934916d8
scl_norm_max = maximum(Arya.histogram2d(get_xy(out[1])..., bins).values)

# ╔═╡ a9e6cdcc-d31b-404a-bb32-4ceb83b90efa
colorrange=colorrange=(1e-5, 1.) .* max_density(out)

# ╔═╡ 666eff06-c57a-4219-ac5e-e68e6b860882
function plot_xy_density!(snap)
	Arya.hist2d!(get_xy(snap)...,  bins = bins, colorscale=log10, colorrange=colorrange)
end


# ╔═╡ 962d5024-6afa-452c-bb04-85046971ee2a
pot_lmc = Agama.Potential(file = joinpath(ENV["DWARFS_ROOT"], "agama/potentials", "vasiliev24/L3M11/potential.ini"))

# ╔═╡ 8f864640-4880-4d40-b5dd-9441d110f0ca
norm_max = maximum(projected_potential(pot_lmc, 130))

# ╔═╡ 781bd5dc-ae4b-4a11-87e4-01a7a8f5e804
function plot_xy_lmc!(pot, idx)
	h = log10.(projected_potential(pot, idx)) .- log10.( norm_max)

	heatmap!(bins, bins, h, colormap=transparant_cmap, colorrange=(-0.5, 0.0))
end

# ╔═╡ 82ae6fc1-c4f5-4ab8-9272-3722d432c46b
let
	fig = Figure()
	ax = Axis(fig[1,1])

	plot_xy_lmc!(pot_lmc, 200)
	fig
end

# ╔═╡ 4ccdd469-57fc-4896-a634-76c3a0816455
traj_lmc = CSV.read(joinpath(ENV["DWARFS_ROOT"], "agama/potentials", "vasiliev24/L3M11/trajlmc.txt"), DataFrame, delim=" ", ignorerepeated=true, header=["t", "x", "y", "z", "v_x", "v_y", "v_z"])

# ╔═╡ 5992a8eb-69c0-4f25-bd37-80cebc4353bf
y_lmc = LilGuys.lerp(traj_lmc.t .* t_scale, traj_lmc.y)

# ╔═╡ 6f6a5974-a82b-4c9e-a38e-9d0c61327748
z_lmc = LilGuys.lerp(traj_lmc.t .* t_scale, traj_lmc.z)

# ╔═╡ 463da824-603a-49c2-81ee-431e18e50337
length(out)

# ╔═╡ 50d85a70-5a31-4d0a-b46e-0d2c7671f8d2
idxs = [1, 106, 150, 180, 200, 213]

# ╔═╡ 7197e7f9-0d05-45ad-9d52-849f566bba04
function plot_scalebar!(scale_length=50)
	length_relative = scale_length / (bins[end] - bins[1])

	x0, y0 = 0.05, 0.05
	lines!([x0, x0 + length_relative], [y0, y0], color=:white, space=:relative, linewidth=theme(:linewidth)[] / 2)
	text!(x0, y0, text="$scale_length kpc", color=:white, space=:relative, fontsize=0.8 * theme(:fontsize)[], )
end

# ╔═╡ c5731120-b0e7-4cd5-aef1-4f3e9c86d50b
colors = Animation.parse_colors.(["green", "purple"])  .|> Makie.to_color

# ╔═╡ 7ee77dcd-ba58-4909-878b-a61de0820bac
Γ = 0.5

# ╔═╡ 5e48acfc-085f-43c8-a4ae-d0aaabcae4ee
let
	fig = Figure()

	for (i, idx) in enumerate(idxs)
		ii, jj = (i-1) ÷ 3 , (i-1) % 3 
		axis = Axis(fig[ii, jj], backgroundcolor=:black, limits=(extrema(bins), extrema(bins)))

		
		x, y = get_xy(out[idx])
		h = Arya.histogram2d(x, y, bins)
		h2 = projected_potential(pot_lmc, idx)

		v1 = h.values ./ scl_norm_max
		@info extrema(v1)
		v2 =  h2 ./ norm_max
		@info extrema(v2)
		img = Animation.combine_densities([v1,v2], colors, dm_power=-2)
		
		image!(extrema(bins), extrema(bins), img)
		
		hidexdecorations!()
		hideydecorations!()
		
		t = (out.times[idx] - out.times[orbit_props["idx_f"]]) * T2GYR
		text!(0.05, 0.95, text= "$(round(t, digits=2)) Gyr", space=:relative, color=:white, align=(:left, :top))
	end

	rowsize!(fig.layout, 1, Aspect(1, 1.0))
	rowsize!(fig.layout, 0, Aspect(1, 1.0))

	rowgap!(fig.layout, 12)
	colgap!(fig.layout, 12)

	@savefig "$(figname)_sim_images_composite"
	fig
end


# ╔═╡ 8e391b46-5bb6-4288-8750-82fe2ae8d0f7
let
	fig = Figure()

	for (i, idx) in enumerate(idxs)
		ii, jj = (i-1) ÷ 3 , (i-1) % 3 
		axis = Axis(fig[ii, jj], backgroundcolor=:black, limits=(extrema(bins), extrema(bins)))
		plot_xy_density!(out[idx])
		poly!(xz_iso.x, xz_iso.z, color=COLORS[9])
		#plot_xy_lmc!(pot_lmc, idx)
		scatter!(y_lmc(out.times[idx]), z_lmc(out.times[idx]), color=COLORS[3], markersize=1.5*theme(:markersize)[])

		if i == 1
			plot_scalebar!()
		end

		if i == 3
			text!(y_lmc(out.times[idx]), z_lmc(out.times[idx]), text="LMC", align=(:center, 1.2), color=COLORS[3])
		end
		
		hidexdecorations!()
		hideydecorations!()
		t = (out.times[idx] - out.times[orbit_props["idx_f"]]) * T2GYR
		text!(0.05, 0.95, text= "$(round(t, digits=2)) Gyr", space=:relative, color=:white, align=(:left, :top))
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
# ╠═b0a0dddc-fb5a-4bb2-b048-a54859d0b703
# ╠═bf49209c-fbfc-4439-a7d8-cfad5ceba8cc
# ╠═9c7a671e-917e-473a-855e-26291e216676
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═4e5ddbe9-e90a-42cf-a0f1-fabb3d3f1675
# ╠═53cdcc20-96d4-4bd5-8028-df9a38af71ae
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═b9600d93-5946-4380-a2ae-9b5f673bbaf5
# ╠═32db23d9-7959-41ac-aff4-b63df5e4b94a
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═a3be2d61-98eb-4037-afb4-4155ba24cc21
# ╠═5a40b893-021b-46e5-a115-0284e13ae7ae
# ╠═f9976cbe-ee5c-4c4b-95d7-94302bfbf7aa
# ╠═dac0e03d-7fee-455d-9fd6-1e2a03ce65da
# ╠═666eff06-c57a-4219-ac5e-e68e6b860882
# ╠═185cfe93-b946-4e72-be6f-853205c49d22
# ╠═5f036156-64e5-4127-a80d-614f8f526559
# ╠═d7771353-51a7-4461-a8d8-17b18623f46c
# ╠═147e195d-99f1-4c33-a65d-2ccfd7fa9bca
# ╠═4078aa8f-4f2f-4f09-9706-4600838d9960
# ╠═336bc60f-9cdc-4ad8-ab87-255654350c61
# ╠═653f3273-6c1f-49f9-a48b-bd395318a551
# ╠═781bd5dc-ae4b-4a11-87e4-01a7a8f5e804
# ╠═8f864640-4880-4d40-b5dd-9441d110f0ca
# ╠═df6c58bd-385a-4f04-a652-c5de934916d8
# ╠═82ae6fc1-c4f5-4ab8-9272-3722d432c46b
# ╠═a9e6cdcc-d31b-404a-bb32-4ceb83b90efa
# ╠═962d5024-6afa-452c-bb04-85046971ee2a
# ╠═4ccdd469-57fc-4896-a634-76c3a0816455
# ╠═5992a8eb-69c0-4f25-bd37-80cebc4353bf
# ╠═6f6a5974-a82b-4c9e-a38e-9d0c61327748
# ╠═463da824-603a-49c2-81ee-431e18e50337
# ╠═50d85a70-5a31-4d0a-b46e-0d2c7671f8d2
# ╠═7197e7f9-0d05-45ad-9d52-849f566bba04
# ╠═c5731120-b0e7-4cd5-aef1-4f3e9c86d50b
# ╠═7ee77dcd-ba58-4909-878b-a61de0820bac
# ╠═5e48acfc-085f-43c8-a4ae-d0aaabcae4ee
# ╠═8e391b46-5bb6-4288-8750-82fe2ae8d0f7
