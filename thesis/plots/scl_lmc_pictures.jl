### A Pluto.jl notebook ###
# v0.20.19

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

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ f419aa46-1f1c-49d4-ae19-95ca7fa90217
modelnames = TOML.parsefile("model_key.toml")

# ╔═╡ 60f13e7b-240b-4f0b-8ffa-515c99a8632e
md"""
# Setup
"""

# ╔═╡ 4e5ddbe9-e90a-42cf-a0f1-fabb3d3f1675
xy_iso = CSV.read("resources/L3M11_iso_xy.csv", DataFrame)

# ╔═╡ 53cdcc20-96d4-4bd5-8028-df9a38af71ae
xz_iso = CSV.read("resources/L3M11_iso_xz.csv", DataFrame)

# ╔═╡ b9600d93-5946-4380-a2ae-9b5f673bbaf5
modelname = joinpath(modelnames["scl_lmc"][1:2]...)

# ╔═╡ 714e6067-3ffa-49a3-93a3-a45eb9d24c97
orbit = Orbit(joinpath(ENV["DWARFS_ROOT"], "analysis", modelname, "centres.hdf5"))

# ╔═╡ a3be2d61-98eb-4037-afb4-4155ba24cc21
orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "analysis", modelname, "orbital_properties.toml"))

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
out = Output(joinpath(ENV["DWARFS_ROOT"], "analysis", modelname));

# ╔═╡ 4ccdd469-57fc-4896-a634-76c3a0816455
traj_lmc = CSV.read(joinpath(ENV["DWARFS_ROOT"], "agama/potentials", "vasiliev24/L3M11/trajlmc.txt"), DataFrame, delim=" ", ignorerepeated=true, header=["t", "x", "y", "z", "v_x", "v_y", "v_z"])

# ╔═╡ 93b3145f-8ac5-41a2-bf4e-fccd5a280076
orbit_lmc = let 
	o = Orbit(positions=[traj_lmc.x traj_lmc.y traj_lmc.z]', velocities=[traj_lmc.v_x traj_lmc.v_y traj_lmc.v_z]' / V2KMS, times=traj_lmc.t / T2GYR)

	LilGuys.resample(o, orbit.times)
end

# ╔═╡ b0ac7331-a39f-48d1-8813-fcd5059bf9e4
md"""
# Plot utils
"""

# ╔═╡ f9976cbe-ee5c-4c4b-95d7-94302bfbf7aa
function get_xy(snap)
	return snap.positions[2, :] .- 0snap.x_cen[2], snap.positions[3, :] .- 0snap.x_cen[3]
end

# ╔═╡ 5a40b893-021b-46e5-a115-0284e13ae7ae
bins = LinRange(-150, 150, 512)

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


# ╔═╡ 147e195d-99f1-4c33-a65d-2ccfd7fa9bca
t_scale = Agama.time_scale(Agama.VASILIEV_UNITS)

# ╔═╡ 5992a8eb-69c0-4f25-bd37-80cebc4353bf
y_lmc = LilGuys.lerp(traj_lmc.t .* t_scale, traj_lmc.y)

# ╔═╡ 6f6a5974-a82b-4c9e-a38e-9d0c61327748
z_lmc = LilGuys.lerp(traj_lmc.t .* t_scale, traj_lmc.z)

# ╔═╡ 50d85a70-5a31-4d0a-b46e-0d2c7671f8d2
idxs = [1, 212, 300, 360, 400, length(out)]

# ╔═╡ 7197e7f9-0d05-45ad-9d52-849f566bba04
function plot_scalebar!(scale_length=50)
	length_relative = scale_length / (bins[end] - bins[1])

	x0, y0 = 0.05, 0.05
	lines!([x0, x0 + length_relative], [y0, y0], color=:white, space=:relative, linewidth=theme(:linewidth)[] / 2)
	text!(x0, y0, text="$scale_length kpc", color=:white, space=:relative, fontsize=0.8 * theme(:fontsize)[], )
end

# ╔═╡ bd81c08d-0e9d-4a70-98b4-6381e60fbd1e
function plot_orbit_trace!(orbit, idx; d_idx=5, color=(:white, 0.3), future=false)

	x = orbit.positions[2, 1:idx]
	y = orbit.positions[3, 1:idx]
	lines!(x, y, linewidth=theme(:linewidth)[]/2, color=color)

	if future
		x = orbit.positions[2, idx:end]
		y = orbit.positions[3, idx:end]
		lines!(x, y, linestyle=:dot, linewidth=theme(:linewidth)[]/2, color=color)
	end
end

# ╔═╡ 199559d8-0d33-4716-9cab-1e43a8a42a75
md"""
# Plot
"""

# ╔═╡ 32db23d9-7959-41ac-aff4-b63df5e4b94a
figname = Dict(
	"sculptor" => "scl",
	"ursa_minor" => "umi"
)[galaxyname]

# ╔═╡ 8e391b46-5bb6-4288-8750-82fe2ae8d0f7
let
	fig = Figure()

	for (i, idx) in enumerate(idxs)
		ii, jj = (i-1) ÷ 3 , (i-1) % 3 
		axis = Axis(fig[ii, jj], backgroundcolor=:black, limits=(extrema(bins), extrema(bins)))
		plot_xy_density!(out[idx])
		poly!(xz_iso.x, xz_iso.z, color=COLORS[9])
		#plot_xy_lmc!(pot_lmc, idx)
		scatter!(y_lmc(out.times[idx]), z_lmc(out.times[idx]), color=COLORS[3], markersize=1*theme(:markersize)[])

		if i == 1
			plot_scalebar!()
		end

		if i == 3
			text!(y_lmc(out.times[idx]), z_lmc(out.times[idx]), text="LMC", align=(:center, 1.2), color=COLORS[3])
		end
		plot_orbit_trace!(orbit, idx, future=true)

		plot_orbit_trace!(orbit_lmc, idx, color=COLORS[3])
		
		hidexdecorations!()
		hideydecorations!()
		t = (out.times[idx] - out.times[orbit_props["idx_f"]]) * T2GYR
		text!(0.05, 0.95, text= "$(round(t, digits=2)) Gyr", space=:relative, color=:white, align=(:left, :top))
	end

	rowsize!(fig.layout, 1, Aspect(1, 1.0))
	rowsize!(fig.layout, 0, Aspect(1, 1.0))

	rowgap!(fig.layout, 6)
	colgap!(fig.layout, 6)

	resize_to_layout!()
	@savefig "$(figname)_lmc_sim_images"
	fig
end


# ╔═╡ Cell order:
# ╠═913e0316-a04b-4270-ba31-0ba0f7fdd705
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═b0a0dddc-fb5a-4bb2-b048-a54859d0b703
# ╠═bf49209c-fbfc-4439-a7d8-cfad5ceba8cc
# ╠═9c7a671e-917e-473a-855e-26291e216676
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═f419aa46-1f1c-49d4-ae19-95ca7fa90217
# ╟─60f13e7b-240b-4f0b-8ffa-515c99a8632e
# ╠═4e5ddbe9-e90a-42cf-a0f1-fabb3d3f1675
# ╠═53cdcc20-96d4-4bd5-8028-df9a38af71ae
# ╠═b9600d93-5946-4380-a2ae-9b5f673bbaf5
# ╠═714e6067-3ffa-49a3-93a3-a45eb9d24c97
# ╠═a3be2d61-98eb-4037-afb4-4155ba24cc21
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═4ccdd469-57fc-4896-a634-76c3a0816455
# ╠═93b3145f-8ac5-41a2-bf4e-fccd5a280076
# ╟─b0ac7331-a39f-48d1-8813-fcd5059bf9e4
# ╠═f9976cbe-ee5c-4c4b-95d7-94302bfbf7aa
# ╠═5a40b893-021b-46e5-a115-0284e13ae7ae
# ╠═dac0e03d-7fee-455d-9fd6-1e2a03ce65da
# ╠═666eff06-c57a-4219-ac5e-e68e6b860882
# ╠═a9e6cdcc-d31b-404a-bb32-4ceb83b90efa
# ╠═147e195d-99f1-4c33-a65d-2ccfd7fa9bca
# ╠═5992a8eb-69c0-4f25-bd37-80cebc4353bf
# ╠═6f6a5974-a82b-4c9e-a38e-9d0c61327748
# ╠═50d85a70-5a31-4d0a-b46e-0d2c7671f8d2
# ╠═7197e7f9-0d05-45ad-9d52-849f566bba04
# ╠═bd81c08d-0e9d-4a70-98b4-6381e60fbd1e
# ╟─199559d8-0d33-4716-9cab-1e43a8a42a75
# ╠═8e391b46-5bb6-4288-8750-82fe2ae8d0f7
# ╠═32db23d9-7959-41ac-aff4-b63df5e4b94a
