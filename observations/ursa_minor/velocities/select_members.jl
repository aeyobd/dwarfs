### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	using Arya
	using CairoMakie

	import LilGuys as lguys
end

# ╔═╡ 93838644-cad6-4df3-b554-208b7afeb3b8
using PyFITS

# ╔═╡ bd6dfd17-02ee-4855-be37-fecfdab6776f
using LilGuys; FIGDIR = "figures"

# ╔═╡ 72f1febc-c6ea-449a-8cec-cd0e49c4e20c
using DataFrames

# ╔═╡ d7d439be-b77b-4a3d-9fb6-7dd7583dc52e
using OrderedCollections

# ╔═╡ 035cdedb-da23-4cfd-aa19-a3aff089d3ac
using Measurements

# ╔═╡ 3ed8c28f-5908-42dc-a56b-24a9b2685a07
md"""
## Inputs
"""

# ╔═╡ 50488b8f-6886-4191-8778-af66929f1445
begin 
	rv_file = "rv_apogee.fits"
	j24_sample = "2c"
	psat_min = 0.2
end

# ╔═╡ 680e7f76-cb4d-40d6-9a9f-d4672427a633
md"""
## derived
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "processed"

# ╔═╡ 86fe351f-ef12-474a-85cc-c10c22a65e77
outname = splitext(basename(rv_file))[1] * "_x_" * j24_sample * "_psat_$psat_min"

# ╔═╡ 7330c75e-1bf9-476a-8274-ebc86d555e6f
md"""
# RV sample models
"""

# ╔═╡ 6bec7416-40c8-4e2b-9d3d-14aa19e5642d
md"""
This notebook takes the dataset created from velocity_xmatch.jl and analyzes it using MCMC to estimate the velocity dispersion and search for any possible velocity gradients.
"""

# ╔═╡ 34e43f4a-bcff-41cb-92c4-0c8d600fd053
import CSV

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, sem

# ╔═╡ d2888213-61e3-4a6f-872b-48a075640ef5
import TOML

# ╔═╡ b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
CairoMakie.activate!(type=:png)

# ╔═╡ 2ab0018f-1628-4f5e-b7da-370eb20c00d0
module RVUtils
	include("../../rv_utils.jl")
end

# ╔═╡ 5ec475a1-14bb-40f6-856a-69fa9efe087a
⊕ = RVUtils.:⊕

# ╔═╡ 7e086680-983e-429b-a917-5330f88c7c55
module GaiaFilters
	include("../../../utils/gaia_filters.jl")
end

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
obs_properties = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/ursa_minor/observed_properties.toml")

# ╔═╡ c4b08cab-0039-4782-8080-b9ef3a6f6d98
filt_params = GaiaFilters.GaiaFilterParams(obs_properties, filename="../data/jensen+24_$j24_sample.fits")

# ╔═╡ feca883d-cbb9-435a-966d-89fefb69ca49
j24 = GaiaFilters.read_gaia_stars(filt_params)

# ╔═╡ 1bc6b7f5-4884-479d-b4be-54f28c2e0a8a
f_sat = TOML.parsefile("../j+24_fsat.toml")[j24_sample]["f_sat"]

# ╔═╡ 019d1bfd-b6d0-4178-b644-cc6ca45b66ae
mean(j24.PSAT[j24.F_BEST .=== 1.])

# ╔═╡ 66c35421-f6d3-4b2d-86e4-319f5476b222
σv = obs_properties["sigma_v"]

# ╔═╡ 2d110151-f7f3-4b09-8684-a2fc4327814b
Δv_gsr = RVUtils.rv_gsr_shift(obs_properties["ra"], obs_properties["dec"])

# ╔═╡ 84509e42-8484-410a-8a76-38473b9f4b71
rv0 = obs_properties["radial_velocity"] .- Δv_gsr

# ╔═╡ 66682e56-4ad9-4823-99da-fc599882eb41
rv_all = read_fits(joinpath(data_dir, rv_file))

# ╔═╡ b20720ac-2787-4cf7-a44b-0cb5293a00b9
@assert :L_PM_SAT ∉ names(rv_all)

# ╔═╡ 097a102b-d6a6-444d-9761-ecb06d64d07f
icrs0 = RVUtils.icrs(obs_properties)

# ╔═╡ fc3ec4d1-5812-4177-90b8-29bbba72d720
gsr0 = lguys.transform(GSR, icrs0)

# ╔═╡ 4dac920b-8252-48a7-86f5-b9f96de6aaa0
rv_meas = let
	rv_meas = rv_all[.!ismissing.(rv_all.source_id), :]
	disallowmissing!(rv_meas, :source_id)
	leftjoin!(rv_meas, j24, on=:source_id, makeunique=true)
	rv_meas = rv_meas[rv_meas.F_BEST .== 1.0, :]

	RVUtils.add_gsr!(rv_meas, 
					 distance=obs_properties["distance"],
					pmra=obs_properties["pmra"], pmdec=obs_properties["pmdec"])

	Δrv = RVUtils.rv_correction(rv_meas.ra, rv_meas.dec, gsr0)

	rv_meas[!, :delta_rv] = Measurements.value.(Δrv)
	rv_meas[!, :delta_rv_err] = Measurements.uncertainty.(Δrv)

	
	rv_meas[:, :vz] .= rv_meas.radial_velocity_gsr .+ rv_meas.delta_rv
	rv_meas[:, :vz_err] .= rv_meas.RV_err

	RVUtils.add_PSAT_RV!(rv_meas; sigma_v=σv, radial_velocity_gsr=rv0, f_sat=f_sat)
	
	rv_meas
end

# ╔═╡ b76a094e-11a4-4d43-869b-27c7f0c2eaee
hist(rv_meas.delta_rv_err ./ rv_meas.RV_err)

# ╔═╡ c28071fe-6077-43ad-b930-604483d5eb28
md"""
## Membership
"""

# ╔═╡ 733fe42e-b7a5-4285-8c73-9a41e4488d40
memb_filt = (rv_meas.PSAT_RV .> psat_min) .&
	rv_meas.F_scatter .&
	(rv_meas.F_BEST .== 1)

# ╔═╡ cc6c65db-ef57-4745-8ada-e11427274a77
memb_stars = rv_meas[memb_filt, :]

# ╔═╡ 2c69332e-9341-4805-9066-10df38fa32fe
"$(outname)_memb.fits"

# ╔═╡ f9f5759a-9a4a-4765-81c8-d513c8a0a181
write_fits(joinpath(data_dir, "$(outname).fits"), memb_stars, overwrite=true)

# ╔═╡ 92f4c8b3-442d-4a56-b05f-1fc547231508
"RV_gmos" ∈ names(memb_stars) && memb_stars[.!ismissing.(memb_stars.RV_gmos), [:source_id, :ra, :dec, :PSAT]]

# ╔═╡ ea3d420f-00f8-4ca2-a49d-e26b48e50afd
nonmemb_stars = rv_meas[.!memb_filt, :]

# ╔═╡ a162219f-df5a-41d8-bf54-927a355f6431
write_fits(joinpath(data_dir, "$(outname)_nonmemb.fits"), nonmemb_stars, overwrite=true)

# ╔═╡ 226fe323-0e2b-4adb-ae53-5b50feafc05a


# ╔═╡ da83eff8-1459-428d-94b5-d0c82489c6c4
memb_filt_bin = (rv_meas.PSAT_RV .> psat_min) .&
	 (.!rv_meas.F_scatter) .&
	(rv_meas.F_BEST .== 1)

# ╔═╡ ee2f1f6f-5b31-4ff9-9fb9-f0f4768ca56f
memb_bin = rv_meas[memb_filt_bin, :]

# ╔═╡ 7f62cdb7-3003-46e5-8cec-e4c75b98d577
write_fits(joinpath(data_dir, "$(outname)_bin.fits"), memb_bin, overwrite=true)

# ╔═╡ 55ce0f69-8a96-4bbb-a59f-ee6503624ea6
md"""
# Numbers
"""

# ╔═╡ 31a6c2e4-538c-4adc-bbda-5043680b17f7
extrema(memb_stars.RV)

# ╔═╡ 3377f632-713d-4fec-84a9-b0211b02cb43
median_err = median(memb_stars.RV_err)

# ╔═╡ 064dff05-192d-41bb-993f-5fb5abc36ecd
number_qual = sum(rv_meas.F_scatter .& (rv_meas.F_BEST .== 1.0))

# ╔═╡ a1938588-ca40-4844-ab82-88c4254c435b
number_memb = length(memb_stars.RV)

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

# ╔═╡ 081c30c7-28e9-4155-80a5-d35317ba926a
hist(nonmemb_stars.RV)

# ╔═╡ 68edc01b-496e-466e-9980-83a586b0bb82
hist(memb_stars.RV)

# ╔═╡ fb52ac04-1483-471f-a164-9bbe15464378
scatter(nonmemb_stars.xi, nonmemb_stars.eta)

# ╔═╡ 6493a62d-cdb7-4831-9da6-25035e3cb7c5
scatter(memb_stars.xi, memb_stars.eta, markersize=2, alpha=0.4)

# ╔═╡ Cell order:
# ╟─3ed8c28f-5908-42dc-a56b-24a9b2685a07
# ╠═50488b8f-6886-4191-8778-af66929f1445
# ╟─680e7f76-cb4d-40d6-9a9f-d4672427a633
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═86fe351f-ef12-474a-85cc-c10c22a65e77
# ╟─7330c75e-1bf9-476a-8274-ebc86d555e6f
# ╟─6bec7416-40c8-4e2b-9d3d-14aa19e5642d
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═34e43f4a-bcff-41cb-92c4-0c8d600fd053
# ╠═93838644-cad6-4df3-b554-208b7afeb3b8
# ╠═bd6dfd17-02ee-4855-be37-fecfdab6776f
# ╠═72f1febc-c6ea-449a-8cec-cd0e49c4e20c
# ╠═d7d439be-b77b-4a3d-9fb6-7dd7583dc52e
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═035cdedb-da23-4cfd-aa19-a3aff089d3ac
# ╠═d2888213-61e3-4a6f-872b-48a075640ef5
# ╠═b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
# ╠═2ab0018f-1628-4f5e-b7da-370eb20c00d0
# ╠═5ec475a1-14bb-40f6-856a-69fa9efe087a
# ╠═7e086680-983e-429b-a917-5330f88c7c55
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╠═c4b08cab-0039-4782-8080-b9ef3a6f6d98
# ╠═feca883d-cbb9-435a-966d-89fefb69ca49
# ╠═1bc6b7f5-4884-479d-b4be-54f28c2e0a8a
# ╠═019d1bfd-b6d0-4178-b644-cc6ca45b66ae
# ╠═66c35421-f6d3-4b2d-86e4-319f5476b222
# ╠═2d110151-f7f3-4b09-8684-a2fc4327814b
# ╠═84509e42-8484-410a-8a76-38473b9f4b71
# ╠═66682e56-4ad9-4823-99da-fc599882eb41
# ╠═b20720ac-2787-4cf7-a44b-0cb5293a00b9
# ╠═fc3ec4d1-5812-4177-90b8-29bbba72d720
# ╠═097a102b-d6a6-444d-9761-ecb06d64d07f
# ╠═4dac920b-8252-48a7-86f5-b9f96de6aaa0
# ╠═b76a094e-11a4-4d43-869b-27c7f0c2eaee
# ╟─c28071fe-6077-43ad-b930-604483d5eb28
# ╠═733fe42e-b7a5-4285-8c73-9a41e4488d40
# ╠═cc6c65db-ef57-4745-8ada-e11427274a77
# ╠═2c69332e-9341-4805-9066-10df38fa32fe
# ╠═f9f5759a-9a4a-4765-81c8-d513c8a0a181
# ╠═92f4c8b3-442d-4a56-b05f-1fc547231508
# ╠═ea3d420f-00f8-4ca2-a49d-e26b48e50afd
# ╠═a162219f-df5a-41d8-bf54-927a355f6431
# ╠═226fe323-0e2b-4adb-ae53-5b50feafc05a
# ╠═7f62cdb7-3003-46e5-8cec-e4c75b98d577
# ╠═da83eff8-1459-428d-94b5-d0c82489c6c4
# ╠═ee2f1f6f-5b31-4ff9-9fb9-f0f4768ca56f
# ╟─55ce0f69-8a96-4bbb-a59f-ee6503624ea6
# ╠═31a6c2e4-538c-4adc-bbda-5043680b17f7
# ╠═3377f632-713d-4fec-84a9-b0211b02cb43
# ╠═064dff05-192d-41bb-993f-5fb5abc36ecd
# ╠═a1938588-ca40-4844-ab82-88c4254c435b
# ╠═b8ee2ea3-98aa-44c4-a309-1a356feb0686
# ╠═b1230b9d-e3ea-4330-b7f5-e708f08db51c
# ╟─cb3bc2ab-8ed7-493d-9868-f793fe24bc42
# ╠═081c30c7-28e9-4155-80a5-d35317ba926a
# ╠═68edc01b-496e-466e-9980-83a586b0bb82
# ╠═fb52ac04-1483-471f-a164-9bbe15464378
# ╠═6493a62d-cdb7-4831-9da6-25035e3cb7c5
