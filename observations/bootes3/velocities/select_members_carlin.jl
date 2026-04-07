### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	using Arya
	using CairoMakie

	using LilGuys; FIGDIR = "figures"

	using OrderedCollections
	using PyFITS
	import TOML
	import CSV
	using DataFrames

	import StatsBase: quantile, mean, std, median, sem
end

# ╔═╡ 3114c0ba-3332-4da2-aba5-f6e9108e6215
md"""
This notebook loads in a file, selects the members according to J+24 probabilities, and finally corrects the velocities and saves the resulting, cleaned file.
"""

# ╔═╡ 3ed8c28f-5908-42dc-a56b-24a9b2685a07
md"""
## Inputs
"""

# ╔═╡ 50488b8f-6886-4191-8778-af66929f1445
begin 
	rv_file = "rv_deimos.fits"
	j24_sample = "2c"
end

# ╔═╡ 680e7f76-cb4d-40d6-9a9f-d4672427a633
md"""
## derived
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "processed"

# ╔═╡ 86fe351f-ef12-474a-85cc-c10c22a65e77
outname = "carlin18_x_" * j24_sample

# ╔═╡ 7330c75e-1bf9-476a-8274-ebc86d555e6f
md"""
# RV sample models
"""

# ╔═╡ 6bec7416-40c8-4e2b-9d3d-14aa19e5642d
md"""
This notebook takes the dataset created from velocity_xmatch.jl and analyzes it using MCMC to estimate the velocity dispersion and search for any possible velocity gradients.
"""

# ╔═╡ b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
CairoMakie.activate!(type=:png)

# ╔═╡ 2ab0018f-1628-4f5e-b7da-370eb20c00d0
module RVUtils
	include("../../rv_utils.jl")
end

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ feca883d-cbb9-435a-966d-89fefb69ca49
j24 = RVUtils.read_gaia_stars("../data/j24_2c.fits", obs_properties,)

# ╔═╡ c28071fe-6077-43ad-b930-604483d5eb28
md"""
## Membership
"""

# ╔═╡ 55ce0f69-8a96-4bbb-a59f-ee6503624ea6
md"""
# Numbers
"""

# ╔═╡ 6152be53-df8f-47e7-9427-0108ced29a78
to_float(x) = @. parse(Float64, replace(x, ["−" => "-"]))

# ╔═╡ 517e8f01-fc64-4512-9bf1-1723fdaaf7d0
rv_all = let
	df = CSV.read("../data/carlin_sand2018.txt", DataFrame, ignorerepeated=true, delim=" ")

	filt_xmatch_carlin, idx_xmatch_carlin = RVUtils.xmatch(df, j24, 1)

	df[!, :source_id] = j24.source_id[idx_xmatch_carlin]

	df[!, :fe_h] = to_float.(df.fe_h)
	df[!, :pmra] = to_float.(df.pmra)
	df[!, :pmdec] = to_float.(df.pmdec)
	rename!(df, Dict(:pmra=>:pmra_carlin, :pmdec=>:pmdec_carlin))
	df

end

# ╔═╡ 295b0abf-b562-4d80-a9e6-8a946d04e5d7
boo3_carlin_memb = rv_all[rv_all.member .== Ref("Y"), :]

# ╔═╡ b20720ac-2787-4cf7-a44b-0cb5293a00b9
@assert :L_PM_SAT ∉ names(rv_all)

# ╔═╡ 4dac920b-8252-48a7-86f5-b9f96de6aaa0
rv_meas = RVUtils.xmatch_and_clean(rv_all, j24, obs_properties)

# ╔═╡ 733fe42e-b7a5-4285-8c73-9a41e4488d40
memb_filt =  rv_all.member .== Ref("Y")

# ╔═╡ cc6c65db-ef57-4745-8ada-e11427274a77
memb_stars = rv_meas[memb_filt, :]

# ╔═╡ 31a6c2e4-538c-4adc-bbda-5043680b17f7
extrema(memb_stars.RV)

# ╔═╡ 3377f632-713d-4fec-84a9-b0211b02cb43
median_err = median(memb_stars.RV_err)

# ╔═╡ a1938588-ca40-4844-ab82-88c4254c435b
number_memb = length(memb_stars.RV)

# ╔═╡ 064dff05-192d-41bb-993f-5fb5abc36ecd
number_qual = length(rv_all.RV)

# ╔═╡ a162219f-df5a-41d8-bf54-927a355f6431
write_fits(joinpath(data_dir, "$outname.fits"), memb_stars, overwrite=true)

# ╔═╡ b8ee2ea3-98aa-44c4-a309-1a356feb0686
sample_info = OrderedDict(
	"median_err" => median_err,
	"number_qual" => number_qual,
	"number_memb" => number_memb,
)

# ╔═╡ b1230b9d-e3ea-4330-b7f5-e708f08db51c
open(joinpath(data_dir, "$outname.info.toml"), "w") do f
	TOML.print(f, sample_info)
end

# ╔═╡ cb3bc2ab-8ed7-493d-9868-f793fe24bc42
md"""
# Plots
"""

# ╔═╡ ea3d420f-00f8-4ca2-a49d-e26b48e50afd
nonmemb_stars = rv_meas[.!memb_filt, :]

# ╔═╡ 081c30c7-28e9-4155-80a5-d35317ba926a
hist(nonmemb_stars.RV)

# ╔═╡ 68edc01b-496e-466e-9980-83a586b0bb82
hist(memb_stars.RV)

# ╔═╡ fb52ac04-1483-471f-a164-9bbe15464378
scatter(nonmemb_stars.xi, nonmemb_stars.eta)

# ╔═╡ 6493a62d-cdb7-4831-9da6-25035e3cb7c5
scatter(memb_stars.xi, memb_stars.eta, markersize=2, alpha=0.4)

# ╔═╡ 9e23b685-8a84-421e-8539-e5c1ae87d53b
let
	fig = Figure(
		size=(5*72, 3*72)
	)
	
	ax = Axis(fig[1,1],
		xlabel = "R / arcmin",
		ylabel = L"RV / km s$^{-1}$",
		#limits=(nothing, (60, 150))
	)

	scatter!(rv_meas.R_ell[.!memb_filt], rv_meas.RV[.!memb_filt])

	scatter!(memb_stars.R_ell, memb_stars.RV)



	fig

end

# ╔═╡ 831750b9-ef44-47a9-bced-5bccba841a11
let
	fig = Figure()
	ax = Axis(fig[1,1],
			  xlabel = "RV deimos",
			  ylabel="RV carlin",
			  aspect=DataAspect()
			  
			 )

	x = rv_all[idx_xmatch_carlin[filt_xmatch_carlin], "RV"]
	y = boo3_carlin[filt_xmatch_carlin, "RV"]
	x_err = boo3_deimos[idx_xmatch_carlin[filt_xmatch_carlin], "RV_err"]
	y_err = boo3_carlin[filt_xmatch_carlin, "RV_err"]
	errorscatter!(x, y, xerror=x_err, yerror=y_err)

	lines!([140, 220], [140, 220] .+ 11)
	fig
end

# ╔═╡ Cell order:
# ╟─3114c0ba-3332-4da2-aba5-f6e9108e6215
# ╟─3ed8c28f-5908-42dc-a56b-24a9b2685a07
# ╠═50488b8f-6886-4191-8778-af66929f1445
# ╟─680e7f76-cb4d-40d6-9a9f-d4672427a633
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═86fe351f-ef12-474a-85cc-c10c22a65e77
# ╟─7330c75e-1bf9-476a-8274-ebc86d555e6f
# ╟─6bec7416-40c8-4e2b-9d3d-14aa19e5642d
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
# ╠═2ab0018f-1628-4f5e-b7da-370eb20c00d0
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╠═517e8f01-fc64-4512-9bf1-1723fdaaf7d0
# ╠═feca883d-cbb9-435a-966d-89fefb69ca49
# ╠═295b0abf-b562-4d80-a9e6-8a946d04e5d7
# ╠═b20720ac-2787-4cf7-a44b-0cb5293a00b9
# ╠═4dac920b-8252-48a7-86f5-b9f96de6aaa0
# ╟─c28071fe-6077-43ad-b930-604483d5eb28
# ╠═733fe42e-b7a5-4285-8c73-9a41e4488d40
# ╠═cc6c65db-ef57-4745-8ada-e11427274a77
# ╟─55ce0f69-8a96-4bbb-a59f-ee6503624ea6
# ╠═31a6c2e4-538c-4adc-bbda-5043680b17f7
# ╠═3377f632-713d-4fec-84a9-b0211b02cb43
# ╠═064dff05-192d-41bb-993f-5fb5abc36ecd
# ╠═a1938588-ca40-4844-ab82-88c4254c435b
# ╠═6152be53-df8f-47e7-9427-0108ced29a78
# ╠═a162219f-df5a-41d8-bf54-927a355f6431
# ╠═b8ee2ea3-98aa-44c4-a309-1a356feb0686
# ╠═b1230b9d-e3ea-4330-b7f5-e708f08db51c
# ╟─cb3bc2ab-8ed7-493d-9868-f793fe24bc42
# ╠═ea3d420f-00f8-4ca2-a49d-e26b48e50afd
# ╠═081c30c7-28e9-4155-80a5-d35317ba926a
# ╠═68edc01b-496e-466e-9980-83a586b0bb82
# ╠═fb52ac04-1483-471f-a164-9bbe15464378
# ╠═6493a62d-cdb7-4831-9da6-25035e3cb7c5
# ╠═9e23b685-8a84-421e-8539-e5c1ae87d53b
# ╠═831750b9-ef44-47a9-bced-5bccba841a11
