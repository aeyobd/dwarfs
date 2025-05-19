### A Pluto.jl notebook ###
# v0.20.5

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

# ╔═╡ 72f1febc-c6ea-449a-8cec-cd0e49c4e20c
using DataFrames

# ╔═╡ bd6dfd17-02ee-4855-be37-fecfdab6776f
using LilGuys; FIGDIR = "figures"

# ╔═╡ 6b527f92-22b6-4a79-8a4b-9add60243ddc
using Measurements

# ╔═╡ 91c33dbb-efb9-47a9-8a1d-6b00d50939d9
f_sat = 0.25

# ╔═╡ 6bec7416-40c8-4e2b-9d3d-14aa19e5642d
md"""
This notebook takes the dataset created from velocity_xmatch.jl and analyzes it using MCMC to estimate the velocity dispersion and search for any possible velocity gradients.
"""

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, sem

# ╔═╡ c3298b45-7c5d-4937-8e4c-c87de36a1354
import DensityEstimators: histogram, bins_equal_number

# ╔═╡ d2888213-61e3-4a6f-872b-48a075640ef5
import TOML

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

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "processed"

# ╔═╡ 66c35421-f6d3-4b2d-86e4-319f5476b222
# ╠═╡ disabled = true
#=╠═╡
σv = obs_properties["sigma_v"]
  ╠═╡ =#

# ╔═╡ 5ec475a1-14bb-40f6-856a-69fa9efe087a
⊕ = RVUtils.:⊕

# ╔═╡ 84509e42-8484-410a-8a76-38473b9f4b71
# ╠═╡ disabled = true
#=╠═╡
rv0 = obs_properties["radial_velocity"] .- Δv_gsr
  ╠═╡ =#

# ╔═╡ 9e2420ea-8d47-4eab-a4bd-0caeb09d9ebb


# ╔═╡ 3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
obs_properties = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml")

# ╔═╡ 2d110151-f7f3-4b09-8684-a2fc4327814b
Δv_gsr = RVUtils.rv_gsr_shift(obs_properties["ra"], obs_properties["dec"])

# ╔═╡ 8ef1fbe4-34c1-414b-92ae-66c1974f58f5
σv = obs_properties["sigma_v"]

# ╔═╡ 062994bc-fb0c-4f08-b7f7-7cc6714bad1e
md"""
## Gradient
"""

# ╔═╡ 88f2918e-e126-420a-96a2-5746a8010f73
icrs0 = lguys.ICRS(obs_properties)

# ╔═╡ 4a473039-79f0-4d77-aa0c-681e2fba4f4c
gsr0 = lguys.transform(lguys.GSR, icrs0)

# ╔═╡ 3c651d98-bbfc-49cb-a2ff-b95a62684310
icrs1 = lguys.ICRS(ra=icrs0.ra + 2, dec=icrs0.dec + 2, distance=icrs0.distance, 
				  pmra=icrs0.pmra, pmdec=icrs0.pmdec, radial_velocity=icrs0.radial_velocity)

# ╔═╡ adb55fa9-8917-470f-bdac-e385c0d704e5
gsr1 = lguys.transform(lguys.GSR, icrs1)

# ╔═╡ c3daa23e-8782-46b5-8bcb-787cdda117f0
gsr1.pmra, gsr1.pmdec

# ╔═╡ 05c9a232-41d7-4ae7-b080-2ab7ddb831ec
gsr0.pmra, gsr0.pmdec

# ╔═╡ ed35eb68-74f7-4009-9b68-dfca2ea547af
pm_gsr_induced = lguys.transform(lguys.GSR, lguys.ICRS(ra=icrs0.ra, dec=icrs0.dec, distance=icrs0.distance, pmra=0, pmdec=0, radial_velocity=0))

# ╔═╡ 10911e3c-b9d8-4165-92d3-c1e6ae44b341
Δv = RVUtils.rv_gsr_shift(icrs0.ra, icrs0.dec)

# ╔═╡ 80d489b4-b269-4042-8060-4582a6589458
rv0 = obs_properties["radial_velocity"] - Δv

# ╔═╡ 4dac920b-8252-48a7-86f5-b9f96de6aaa0
rv_meas = let
	rv_meas = read_fits(joinpath(data_dir, "rv_combined.fits"))


	obs = [lguys.ICRS(ra=r.ra, dec=r.dec, 
		distance=obs_properties["distance"], 
						pmra=r.pmra, pmdec=r.pmdec, radial_velocity=r.RV)
	for r in eachrow(rv_meas)]
		
	obs_gsr = lguys.to_frame(lguys.transform.(lguys.GSR, obs))


	rv_meas[!, :pmra_gsr] = obs_gsr.pmra
	rv_meas[!, :pmdec_gsr] = obs_gsr.pmdec
	rv_meas[!, :radial_velocity_gsr] = obs_gsr.radial_velocity

	rv_meas[:, :vz] .= rv_meas.radial_velocity_gsr .+ rv_meas.delta_rv
	rv_meas[:, :vz_err] .= rv_meas.RV_err .⊕ rv_meas.delta_rv_err
	

	rv_meas[:, :L_RV_SAT] = RVUtils.L_RV_SAT.(rv_meas.vz, rv_meas.RV_err, rv0, σv)
	rv_meas[:, :L_RV_BKD] = RVUtils.L_RV_BKD.(rv_meas.vz, rv_meas.ra, rv_meas.dec)

	rv_meas[:, :PSAT_RV] = RVUtils.PSAT_RV(rv_meas, f_sat)
	
	rv_meas
end

# ╔═╡ ac3d3fc2-73c6-428c-a9b6-a983d6aa6fd4
let
	fig = Figure()

	ax=Axis(fig[1,1];
		 aspect=DataAspect(),
		 xlabel = L"$\xi$ / arcmin",
		 ylabel = L"$\eta$ / arcmin",
		xreversed=true
		)
	
	p = scatter!(rv_meas.xi, rv_meas.eta, color=rv_meas.delta_rv,
			colormap=:bluesreds, 
		colorrange=(-4, 4)
	  )

	arrows!([0], [0], [20gsr0.pmra], [20gsr0.pmdec])

	Colorbar(fig[1,2], p, label="projection correction (km/s)")

	fig
end

# ╔═╡ 3702b5c0-d276-41dd-ab0c-33aff9ee3b12
let
	fig = Figure()

	ax=Axis(fig[1,1];
		 aspect=DataAspect(),
		 xlabel = L"$\xi$ / arcmin",
		 ylabel = L"$\eta$ / arcmin",
		xreversed=true
		)
	
	p = scatter!(rv_meas.xi, rv_meas.eta, color=-RVUtils.rv_gsr_shift.(rv_meas.ra, rv_meas.dec) .+ Δv,
			colormap=:bluesreds, 
		colorrange=(-4, 4)
	  )

	arrows!([0], [0], [-20pm_gsr_induced.pmra], [-20pm_gsr_induced.pmdec])
	
	Colorbar(fig[1,2], p, label="gsr correction (km/s)")

	fig
end

# ╔═╡ b44863f5-a051-4ee1-a6c5-7d2efc63552c
Δvs_gsr = RVUtils.rv_gsr_shift.(rv_meas.ra, rv_meas.dec)

# ╔═╡ a4584ae3-e295-49b9-8beb-8e644a724836
let
	fig = Figure()

	ax=Axis(fig[1,1];
		 aspect=DataAspect(),
		 xlabel = L"$\xi$ / arcmin",
		 ylabel = L"$\eta$ / arcmin",
		xreversed=true
		)
	
	p =	scatter!(rv_meas.xi, rv_meas.eta, color=rv_meas.delta_rv .- (Δvs_gsr .- Δv), colormap=:bluesreds, colorrange=(-1, 1))
	
	arrows!([0], [0], [20icrs0.pmra], [20icrs0.pmdec])


	Colorbar(fig[1,2], p, label="both correction gradient")
	fig
end

# ╔═╡ 95ec0615-d1a3-4021-9012-ab2d76b9e406
md"""
# Proper motion method
"""

# ╔═╡ 97c2ff3f-1d1b-4c59-be9c-1e793fcc23e6
σ_pm = lguys.kms2pm(σv, obs_properties["distance"])

# ╔═╡ 409fd1f7-b2a1-43a7-b1d3-98910e1993ff


# ╔═╡ 9586a846-df6a-4903-9392-b35e84890bed
pmra_0 = gsr0.pmra ± (σ_pm ⊕ obs_properties["pmdec_err"])

# ╔═╡ 2eca9a91-6c23-4e41-b238-200ae4a7e9e4
pmdec_0 = gsr0.pmdec ± (σ_pm ⊕ obs_properties["pmdec_err"])

# ╔═╡ 0a3a44e6-f465-4fd4-9f35-2bc351e2a301
distance = obs_properties["distance"] ± obs_properties["distance_err"]

# ╔═╡ ca67b4fa-26bb-4f0d-ad3a-ea6413940fb6
pmras = map(eachrow(rv_meas)) do row

	w1 = 1 / row.pmra_error^2
	w2 = 1/obs_properties["pmra_err"]^2 * row.PSAT_RV
	ws = [w1, w2]
	xs = [row.pmra, obs_properties["pmra"]]

	lguys.mean(xs, ws) ± sqrt(1/sum(ws))
end

# ╔═╡ e6cfe55e-0dbd-485b-b106-eaf3215a82c5
pmdecs = map(eachrow(rv_meas)) do row

	w1 = 1 / row.pmdec_error^2
	w2 = 1/obs_properties["pmdec_err"]^2 * row.PSAT_RV
	ws = [w1, w2]
	xs = [row.pmdec, obs_properties["pmdec"]]

	@info ws
	lguys.mean(xs, ws) ± sqrt(1/sum(ws))
end

# ╔═╡ 898b6bf5-eca2-4806-b6b7-cdb39d3d5b17
icrss = [lguys.ICRS(ra=rv_meas.ra[i], dec=rv_meas.dec[i], distance=distance, 
				 pmra=pmras[i], pmdec=pmdecs[i], radial_velocity=rv_meas.RV[i])
for i in eachindex(pmras)]

# ╔═╡ aa0df3da-8b28-49a5-b5a9-865fde369f30
gsrs = lguys.transform.(lguys.GSR, icrss)

# ╔═╡ a7eb2112-8294-46cd-89d4-f425d193dc42
scatter(rv_meas.pmra, Measurements.value.(pmras) .- rv_meas.pmra)

# ╔═╡ cb975d25-a552-4569-b250-bf0a5f46a778
scatter(rv_meas.pmdec, Measurements.value.(pmdecs) .- rv_meas.pmdec)

# ╔═╡ a6c627c5-2f4d-4af5-9b0b-2b3255613d1f
pmras_gsr = [g.pmra for g in gsrs]

# ╔═╡ bbdc7335-3038-47ac-b3bf-46069898a962
pmdecs_gsr = [g.pmdec for g in gsrs]

# ╔═╡ c404f5b0-7d0e-4d3b-8f0f-600c5ae7ef7e
Δrv_solar = RVUtils.rv_correction.(rv_meas.ra, rv_meas.dec, icrs0.ra, icrs0.dec,
						   pmras, pmdecs, distance)

# ╔═╡ e486c6a8-a52f-4411-a93b-818463b821db
Δrv_gsr = RVUtils.rv_correction.(rv_meas.ra, rv_meas.dec, icrs0.ra, icrs0.dec,
						   pmras_gsr, pmdecs_gsr, distance)

# ╔═╡ 0c6406b2-8c24-47c7-9c4d-c6ab5016e88b
Δrv_gsr_fixed = RVUtils.rv_correction.(rv_meas.ra, rv_meas.dec, icrs0.ra, icrs0.dec,
						   gsr0.pmra, gsr0.pmdec, Measurements.value.(distance))

# ╔═╡ 4cc199ca-c668-400e-97fa-88aa19184735
scatter(Measurements.value.(pmdecs_gsr .- gsr0.pmra))

# ╔═╡ f4e4ae6f-1794-4e57-af20-3aaf69fa1e9e
lguys.pm2kms(icrs0.pmra ⊕ icrs0.pmdec,  distance) / (180/π)

# ╔═╡ 2d7bc925-c1b0-418d-954f-ca5c12443aad
import LinearAlgebra: normalize, ⋅

# ╔═╡ 388ec2ef-f622-43be-9a06-6f17a239c76f
pm_hat = normalize([icrs0.pmra, icrs0.pmdec])

# ╔═╡ aea39cf7-b3de-49c0-9dc6-b6d2343a4b6a
pm_axis = sum(pm_hat' .* hcat(rv_meas.xi, rv_meas.eta), dims=2)

# ╔═╡ ddbf285e-efca-4ef6-9e4a-33612b68e7e3
scatter(vec(pm_axis) / 60, rv_meas.delta_rv .- (Δvs_gsr .- Δv))

# ╔═╡ e24110f5-7a0e-42bf-9955-51aaf9bcbe66
pm_hat

# ╔═╡ 3129a943-0a91-4bb0-8f19-9f964544ae23
scatter(rv_meas.xi, Measurements.value.(Δrv_gsr_fixed) .- Measurements.value.(Δrv_gsr), color=:transparent, strokewidth=1, alpha=0.1)

# ╔═╡ 26870502-1c37-4f18-b990-5acee5cbd646
hist(Measurements.value.(Δrv_gsr_fixed .- Δrv_gsr)[rv_meas.PSAT_RV .> 0.2])

# ╔═╡ dd5b64cf-0720-45ba-b241-adb1baa19324
let
	fig = Figure()

	ax=Axis(fig[1,1];
		 aspect=DataAspect(),
		 xlabel = L"$\xi$ / arcmin",
		 ylabel = L"$\eta$ / arcmin",
		xreversed=true
		)
	
	p =	scatter!(rv_meas.xi, rv_meas.eta, color=Measurements.value.(Δrv_solar), colormap=:bluesreds, colorrange=(-1, 1), markersize=3)
	
	arrows!([0], [0], [20icrs0.pmra], [20icrs0.pmdec])


	Colorbar(fig[1,2], p, label="naieve")
	fig
end

# ╔═╡ 4599e728-52f1-41a4-b4d5-c0fa908d465e
Δrv_tot = Measurements.value.(Δrv_gsr) .- (Δvs_gsr .- Δv)

# ╔═╡ d7c3fea5-4504-4320-8dd4-d9bedf300fde
scatter(vec(pm_axis) / 60, Measurements.value.(Δrv_solar .- Δrv_tot), alpha=0.1, markersize=1)

# ╔═╡ 6a4018e6-1c62-4281-ac67-5ddb4e61a4d5
let
	fig = Figure()

	ax=Axis(fig[1,1];
		 aspect=DataAspect(),
		 xlabel = L"$\xi$ / arcmin",
		 ylabel = L"$\eta$ / arcmin",
		xreversed=true
		)
	
	p =	scatter!(rv_meas.xi, rv_meas.eta, color=Measurements.value.(Δrv_gsr) .- (Δvs_gsr .- Δv), colormap=:bluesreds, colorrange=(-1, 1), markersize=3)
	
	arrows!([0], [0], [20icrs0.pmra], [20icrs0.pmdec])


	Colorbar(fig[1,2], p, label="both correction gradient")
	fig
end

# ╔═╡ ebe93c7c-1a3e-46df-ba41-9786ae2d7535
let
	fig = Figure()

	ax=Axis(fig[1,1];
		 aspect=DataAspect(),
		 xlabel = L"$\xi$ / arcmin",
		 ylabel = L"$\eta$ / arcmin",
		xreversed=true
		)
	
	p =	scatter!(rv_meas.xi, rv_meas.eta, color=Measurements.value.(Δrv_gsr), colormap=:bluesreds, colorrange=(-1, 1), markersize=3)
	
	arrows!([0], [0], [20icrs0.pmra], [20icrs0.pmdec])


	Colorbar(fig[1,2], p, label="both correction gradient")
	fig
end

# ╔═╡ Cell order:
# ╠═91c33dbb-efb9-47a9-8a1d-6b00d50939d9
# ╠═6bec7416-40c8-4e2b-9d3d-14aa19e5642d
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═93838644-cad6-4df3-b554-208b7afeb3b8
# ╠═72f1febc-c6ea-449a-8cec-cd0e49c4e20c
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═c3298b45-7c5d-4937-8e4c-c87de36a1354
# ╠═d2888213-61e3-4a6f-872b-48a075640ef5
# ╠═b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
# ╠═bd6dfd17-02ee-4855-be37-fecfdab6776f
# ╠═2ab0018f-1628-4f5e-b7da-370eb20c00d0
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═66c35421-f6d3-4b2d-86e4-319f5476b222
# ╠═2d110151-f7f3-4b09-8684-a2fc4327814b
# ╠═5ec475a1-14bb-40f6-856a-69fa9efe087a
# ╠═84509e42-8484-410a-8a76-38473b9f4b71
# ╠═9e2420ea-8d47-4eab-a4bd-0caeb09d9ebb
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╠═8ef1fbe4-34c1-414b-92ae-66c1974f58f5
# ╠═80d489b4-b269-4042-8060-4582a6589458
# ╠═4dac920b-8252-48a7-86f5-b9f96de6aaa0
# ╠═062994bc-fb0c-4f08-b7f7-7cc6714bad1e
# ╠═4a473039-79f0-4d77-aa0c-681e2fba4f4c
# ╠═88f2918e-e126-420a-96a2-5746a8010f73
# ╠═3c651d98-bbfc-49cb-a2ff-b95a62684310
# ╠═adb55fa9-8917-470f-bdac-e385c0d704e5
# ╠═c3daa23e-8782-46b5-8bcb-787cdda117f0
# ╠═05c9a232-41d7-4ae7-b080-2ab7ddb831ec
# ╠═ed35eb68-74f7-4009-9b68-dfca2ea547af
# ╠═ac3d3fc2-73c6-428c-a9b6-a983d6aa6fd4
# ╠═10911e3c-b9d8-4165-92d3-c1e6ae44b341
# ╠═3702b5c0-d276-41dd-ab0c-33aff9ee3b12
# ╠═b44863f5-a051-4ee1-a6c5-7d2efc63552c
# ╠═a4584ae3-e295-49b9-8beb-8e644a724836
# ╠═95ec0615-d1a3-4021-9012-ab2d76b9e406
# ╠═6b527f92-22b6-4a79-8a4b-9add60243ddc
# ╠═97c2ff3f-1d1b-4c59-be9c-1e793fcc23e6
# ╠═409fd1f7-b2a1-43a7-b1d3-98910e1993ff
# ╠═9586a846-df6a-4903-9392-b35e84890bed
# ╠═2eca9a91-6c23-4e41-b238-200ae4a7e9e4
# ╠═0a3a44e6-f465-4fd4-9f35-2bc351e2a301
# ╠═ca67b4fa-26bb-4f0d-ad3a-ea6413940fb6
# ╠═e6cfe55e-0dbd-485b-b106-eaf3215a82c5
# ╠═898b6bf5-eca2-4806-b6b7-cdb39d3d5b17
# ╠═aa0df3da-8b28-49a5-b5a9-865fde369f30
# ╠═a7eb2112-8294-46cd-89d4-f425d193dc42
# ╠═cb975d25-a552-4569-b250-bf0a5f46a778
# ╠═a6c627c5-2f4d-4af5-9b0b-2b3255613d1f
# ╠═bbdc7335-3038-47ac-b3bf-46069898a962
# ╠═c404f5b0-7d0e-4d3b-8f0f-600c5ae7ef7e
# ╠═e486c6a8-a52f-4411-a93b-818463b821db
# ╠═0c6406b2-8c24-47c7-9c4d-c6ab5016e88b
# ╠═4cc199ca-c668-400e-97fa-88aa19184735
# ╠═f4e4ae6f-1794-4e57-af20-3aaf69fa1e9e
# ╠═2d7bc925-c1b0-418d-954f-ca5c12443aad
# ╠═388ec2ef-f622-43be-9a06-6f17a239c76f
# ╠═aea39cf7-b3de-49c0-9dc6-b6d2343a4b6a
# ╠═ddbf285e-efca-4ef6-9e4a-33612b68e7e3
# ╠═e24110f5-7a0e-42bf-9955-51aaf9bcbe66
# ╠═3129a943-0a91-4bb0-8f19-9f964544ae23
# ╠═26870502-1c37-4f18-b990-5acee5cbd646
# ╠═dd5b64cf-0720-45ba-b241-adb1baa19324
# ╠═4599e728-52f1-41a4-b4d5-c0fa908d465e
# ╠═d7c3fea5-4504-4320-8dd4-d9bedf300fde
# ╠═6a4018e6-1c62-4281-ac67-5ddb4e61a4d5
# ╠═ebe93c7c-1a3e-46df-ba41-9786ae2d7535
