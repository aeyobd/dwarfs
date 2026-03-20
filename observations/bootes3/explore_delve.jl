### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 3b4af864-23b1-11f1-8d07-0f1be0ec4672
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using CairoMakie; CairoMakie.activate!(type=:png)
	using Arya
	using PyFITS
end

# ╔═╡ 6d238910-967a-43ce-b35b-fb6c25ea8846
using CSV, DataFrames

# ╔═╡ 8874ace4-f625-4487-82a8-eeb08ef0b3b7
import TOML

# ╔═╡ 151f3201-1331-46de-8282-10b10939cc55
module GaiaFilters
	include("../../utils/gaia_filters.jl")
end

# ╔═╡ 6a443c21-7365-405f-91ae-07ae9e979433
obs_props = TOML.parsefile("observed_properties.toml")

# ╔═╡ 2b9c4e36-23f6-4a1b-9621-cd56507f8c20
function get_isochrone()
    iso_columns = string.(split("Zini     MH   logAge Mini        int_IMF         Mass   logL    logTe  logg  label   McoreTP C_O  period0  period1  period2  period3  period4  pmode  Mloss  tau1m   X   Y   Xc  Xn  Xo  Cexcess  Z 	mbolmag  umag    gmag    rmag    imag    zmag    Ymag", 
	r"\s+"))

    all_isochrones = CSV.read(joinpath("isochrone.iso"), DataFrame,
					 comment="#", ignorerepeated=true, delim = ' ', header=iso_columns)

	
	filt = all_isochrones.label .< 5 # only keep through RGB

	all_isochrones[filt, :]
end

# ╔═╡ 207db1e9-836b-4056-9d1a-a7807d3d8757
iso = get_isochrone()

# ╔═╡ a3c6ca4c-b44d-4c14-aaaf-bbeb921c2a65
stars = let
	df = read_fits("data/delve_dr2_good.fits")
	df[!, :bp_rp] = df.gmag .- df.rmag
	df[!, :G] = df.gmag

	df[!, :xi], df[!, :eta] = 60 .* LilGuys.to_tangent(df.ra, df.dec, obs_props["ra"], obs_props["dec"])
	df
end

# ╔═╡ 608a9262-518a-454f-af4b-5ca209d8130e
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 yreversed = true,
			 xlabel = "g - r",
			 ylabel = "g",
			 limits=(-0.5, 1.0, 14, 24))

	scatter!(stars.gmag .- stars.rmag, stars.gmag, color=:black, alpha=0.3, markersize=0.3)

	lines!(iso.gmag .- iso.rmag, iso.gmag .+ obs_props["distance_modulus"])
	fig
end

# ╔═╡ ec16621c-bb5e-441b-a81f-41ac32b3a911


# ╔═╡ 93cae4a2-84c2-42d1-8737-2cf92a35ae3e
err_model(x, params) = @. (params[1] + x*params[2] + x^2 * params[3])

# ╔═╡ 30d23c64-01a8-459d-8dbe-45771d2a9a85
gr_err = @. sqrt(stars.magerr_psf_g^2 + stars.magerr_psf_r^2)

# ╔═╡ fb4e16bd-533f-417f-b6f8-6f7afc0e4abf
popt, covt = LilGuys.curve_fit(err_model, stars.gmag, log10.(gr_err), [0., 1., 1.])

# ╔═╡ d14c2a16-9a21-4a62-a5c8-ab6cecd30d77
color_err(gmag) = 10 ^ err_model(gmag, popt)

# ╔═╡ 68de652c-93e2-4ce1-aa45-82a3b0590714
let
	f = scatter(stars.gmag, gr_err)

	x = LinRange(14, 24, 1000)
	y = err_model.(x, [popt])
	lines!(x, color_err.(x), color=COLORS[2])
	f
end

# ╔═╡ 7038f7f4-208a-48dd-a62a-36163fffc4a9
function plot_tangent(stars)
	f = scatter(stars.xi, stars.eta, markersize=1, alpha=0.3, color=:black, 
	   axis=(;
			xreversed = true,
			xlabel = "xi / arcmin", 
			ylabel = "eta / arcmin",
			aspect = DataAspect()))


	arc!((0, 0), 60, 0, 2π)
	f

end

# ╔═╡ cd65a0cd-5bac-458c-afd5-0d48553ada61
plot_tangent(stars)

# ╔═╡ acea97d7-34c8-4607-8a59-757ddc3c83b6


# ╔═╡ 2fba1a8e-8052-48b0-9344-4b48bc6ca8a2
function get_iso_selection(dm=obs_props["distance_modulus"])
	x, y = iso.gmag .- iso.rmag, iso.gmag .+ dm

	filt = (y .< 25) .& (iso.label .< 4)
	x = x[filt]
	y = y[filt]

	return x, y
end

# ╔═╡ f3688825-0bef-408d-b6e8-022a357bb80e
function make_iso_selection(; dm=obs_props["distance_modulus"], color_width=0.2)
	x_i, y_i = get_iso_selection(dm)

	return [x_i .- color_width .- color_err.(y_i); reverse(x_i .+ color_width .+ color_err.(y_i))], [y_i; reverse(y_i)]
end

# ╔═╡ cb47e65e-917a-4e7c-9f56-917a990f5212
cmd_cut_x, cmd_cut_y = make_iso_selection()

# ╔═╡ df252a18-5c5f-4575-8fa9-31da97061bed
cmd_cut = vcat(([x, y] for (x, y) in zip(cmd_cut_x, cmd_cut_y))...)

# ╔═╡ 990199df-c227-4fa5-8355-cee4484534dc
print("cmd_cut = $cmd_cut")

# ╔═╡ 9b5a64e6-3936-4cb0-907f-c4668e521ea1
let
	x, y = get_iso_selection()
	f = lines(x, y)
	lines!(x .- 0.2, y)
	lines!(x .+ 0.2, y)
	lines!(x, y .-0.2)
	lines!(x, y .+0.2)

	x, y = make_iso_selection()

	poly!(x, y, alpha=0.2)
	f
end

# ╔═╡ e439ba74-178a-401d-b5f7-6017f78aac81
function filter_iso(stars, dm=obs_props["distance_modulus"], color_width=0.2)
	cmd_cut_x, cmd_cut_y = make_iso_selection(dm=dm, color_width=color_width)

	cmd_cut = vcat(([x, y] for (x, y) in zip(cmd_cut_x, cmd_cut_y))...)

	return stars[GaiaFilters.cmd_filter(stars, cmd_cut), :]
end

# ╔═╡ 0e912255-caa2-4749-8e58-0ef4fde040c9
let
	df = filter_iso(stars)
	plot_tangent(df)
end

# ╔═╡ ccb91c11-53bb-4339-a492-8b1f9c25abc4
let
	df = filter_iso(stars, obs_props["distance_modulus"] + 0.5)
	plot_tangent(df)
end

# ╔═╡ cad05ae7-2562-42b9-8fb7-49ed557671b8
let
	df = filter_iso(stars, obs_props["distance_modulus"] - 0.5)
	plot_tangent(df)
end

# ╔═╡ da7a9d5c-82d8-4e70-aafe-608507e98cbf
let
	df = filter_iso(stars, obs_props["distance_modulus"], 0.1)
	plot_tangent(df)
end

# ╔═╡ 02a7bce0-68b1-4bf7-9bdc-a59d3c027d6f
let
	df = filter_iso(stars, obs_props["distance_modulus"], 0.05)
	plot_tangent(df)
end

# ╔═╡ Cell order:
# ╠═3b4af864-23b1-11f1-8d07-0f1be0ec4672
# ╠═6d238910-967a-43ce-b35b-fb6c25ea8846
# ╠═8874ace4-f625-4487-82a8-eeb08ef0b3b7
# ╠═151f3201-1331-46de-8282-10b10939cc55
# ╠═6a443c21-7365-405f-91ae-07ae9e979433
# ╠═2b9c4e36-23f6-4a1b-9621-cd56507f8c20
# ╠═207db1e9-836b-4056-9d1a-a7807d3d8757
# ╠═a3c6ca4c-b44d-4c14-aaaf-bbeb921c2a65
# ╠═608a9262-518a-454f-af4b-5ca209d8130e
# ╠═ec16621c-bb5e-441b-a81f-41ac32b3a911
# ╠═93cae4a2-84c2-42d1-8737-2cf92a35ae3e
# ╠═30d23c64-01a8-459d-8dbe-45771d2a9a85
# ╠═fb4e16bd-533f-417f-b6f8-6f7afc0e4abf
# ╠═d14c2a16-9a21-4a62-a5c8-ab6cecd30d77
# ╠═68de652c-93e2-4ce1-aa45-82a3b0590714
# ╠═cd65a0cd-5bac-458c-afd5-0d48553ada61
# ╠═7038f7f4-208a-48dd-a62a-36163fffc4a9
# ╠═acea97d7-34c8-4607-8a59-757ddc3c83b6
# ╠═2fba1a8e-8052-48b0-9344-4b48bc6ca8a2
# ╠═f3688825-0bef-408d-b6e8-022a357bb80e
# ╠═cb47e65e-917a-4e7c-9f56-917a990f5212
# ╠═df252a18-5c5f-4575-8fa9-31da97061bed
# ╠═990199df-c227-4fa5-8355-cee4484534dc
# ╠═9b5a64e6-3936-4cb0-907f-c4668e521ea1
# ╠═e439ba74-178a-401d-b5f7-6017f78aac81
# ╠═0e912255-caa2-4749-8e58-0ef4fde040c9
# ╠═ccb91c11-53bb-4339-a492-8b1f9c25abc4
# ╠═cad05ae7-2562-42b9-8fb7-49ed557671b8
# ╠═da7a9d5c-82d8-4e70-aafe-608507e98cbf
# ╠═02a7bce0-68b1-4bf7-9bdc-a59d3c027d6f
