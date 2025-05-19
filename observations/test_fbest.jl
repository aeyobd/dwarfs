### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 944a654e-d955-4ea8-94ab-8aa5b73f9ea6
import Pkg; Pkg.activate()

# ╔═╡ 3d9e96f5-7d9b-4f07-8ca8-6da1bf731e70
using Arya

# ╔═╡ 8cff9c1e-1a19-11f0-2dd6-7bddd47fe438
using PyFITS

# ╔═╡ 3b0e285e-cefd-42ff-8c0e-dbe71ef2b720
using CairoMakie

# ╔═╡ 6a1415e7-2cdc-4fe9-88cb-43563f2141f1
using PythonCall

# ╔═╡ 25728329-5b72-4815-a37f-d24873fd408a
using DataFrames

# ╔═╡ 73596566-f0bf-4966-9e4b-a8d451ce307d
md"""
This notebook sanity checks all of the ingredients in the J+24 model.
"""

# ╔═╡ 061b5be8-fc49-4fc4-9c5d-239079740900
import TOML

# ╔═╡ 31e6bedf-b347-4c83-89f6-c5c1adc7b481
import LilGuys as lguys

# ╔═╡ 0653d025-c0e1-4b98-9701-f3db2f6ebe8d
galaxy = "sculptor"

# ╔═╡ c3933f11-ec53-49b0-a7c6-fc34ed8ef5c6
df = read_fits("$galaxy.fits")

# ╔═╡ 9beeec68-ca34-49a8-8f5b-6be279bacd13
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxy, "observed_properties.toml"))

# ╔═╡ 8515cbfb-1385-40f6-a53a-36183be38375
md"""
# Total
"""

# ╔═╡ ead04430-df91-42f2-9a66-9e47563411db
md"""
Combination of all filters, positive satellite CMD likelihood, and within 2 arcmin.
"""

# ╔═╡ 4b12c940-5202-4903-9b5b-7598b11e4726
md"""
Also includes nonzero satellite likelihood terms and requires distance to be less than field apeture.
"""

# ╔═╡ f35de603-cf15-4a1b-9c00-f4c70308af60
filt_missed = df.F_BEST .!= (df.F_NOTANAN .* df.F_CPAR .* df.F_FLUXEXCESS .* df.F_INGRID .* df.F_ASTROMETRIC .* (df.L_CMD_SAT .> 0) .* (df.L_S_SAT .> 0))

# ╔═╡ 2d407500-5ed6-4f14-9a19-db68616202b8
df_missed = df[filt_missed, :]

# ╔═╡ dd6a2a61-910f-4a14-afeb-ea6016125478
scatter(df_missed.xi, df_missed.eta)

# ╔═╡ a7e52514-c2f7-4969-ac5c-df9c60b1d4bc
md"""
Missed stars are on boundary of selection (rounding error?)
"""

# ╔═╡ ba34ed2a-2a1c-494a-ba1a-8ce45b1d6a04
md"""
## NaN filter
"""

# ╔═╡ 7ffbdc53-5c3d-438c-8071-59e801590b6b
md"""
Remove nans in magnitudes and astrometry. Works exactly
"""

# ╔═╡ 1b1098fd-d7cc-437f-91d2-62c90fcd97e0
sum(ismissing.(df.phot_g_mean_mag))

# ╔═╡ 69e3c284-cc10-4941-a6bf-655f3cbe79b0
sum(ismissing.(df.phot_bp_mean_mag))

# ╔═╡ 9666b674-88c2-4466-b9ce-f19682368362
sum(ismissing.(df.phot_rp_mean_mag))

# ╔═╡ 0898688e-d627-4174-a01c-192d0e7b2158
sum(ismissing.(df.pmra))

# ╔═╡ ecbb1ee7-8ef6-4fc7-a612-c09d29a2aaed
sum(ismissing.(df.parallax))

# ╔═╡ 3eca57d2-3049-4f36-ad15-3219981629dd
missed_nan = (df.F_NOTANAN .== 1) .!= .!(
	ismissing.(df.pmra) .| ismissing.(df.pmdec)
	.| ismissing.(df.phot_g_mean_mag) .| ismissing.(df.bp_rp)
)


# ╔═╡ 02836247-daef-4cd1-9fe9-e963a840132a
sum(missed_nan)

# ╔═╡ 2c344123-e54c-4a94-aa58-7eb9e305084e
md"""
## CPAR filter
"""

# ╔═╡ eccc30af-0447-4456-9026-16a3aa13f964
md"""
Missed parallax values are along boundary (rounding?)
"""

# ╔═╡ 4b086ff9-262e-400c-886c-37ea3cb2e2e8
π_zeropoint = -0.017

# ╔═╡ 1f6d3e45-b7e0-4b6f-915a-27fbb4a30ff5
md"""
# Grid filter
"""

# ╔═╡ 2f3fbd8e-96c4-40e1-91c0-d9e3bd9f2e1a
flag_pm = (abs.(df.pmra) .< 10) .& (abs.(df.pmdec) .< 10) .&
	.!ismissing.(df.pmra) .& .!ismissing.(df.pmdec)

# ╔═╡ 67ab8518-3b43-4d28-840d-690febff5959
minimum(skipmissing(df.phot_g_mean_mag[.!ismissing.(df.PSAT) .& (df.PSAT .> 0.2)]))

# ╔═╡ ccc18cbd-43b3-4d4e-a308-727561f81c30
md"""
Use padova v1.1 with M/H =-1.70 at 12Gyr.
"""

# ╔═╡ bdc91b1f-b8a9-4150-b86a-00de4c25d519
dm = Dict(
	"sculptor" => 19.68
)[galaxy]

# ╔═╡ 6c512f14-93ac-468a-be8c-9479fcb3eca9
TRGB_min = Dict(
	"sculptor" => -3.064
)[galaxy] + dm

# ╔═╡ 2be03779-adac-4d46-af22-2f54f1fec593
dm_err = Dict(
	"sculptor" => 0.17
)[galaxy]

# ╔═╡ c049e7c6-e6e8-4c91-b563-e34545b6a71f
parallax = 1 / lguys.dm2dist(dm)

# ╔═╡ 2da0de80-0409-4f5d-93d5-9757d64c233c
π_galaxy = parallax

# ╔═╡ cfe395e1-ed32-4435-bedf-7d25a443762a
parallax_err = parallax .*  dm_err .* log(10)

# ╔═╡ 874cf235-809f-43e4-9d3c-b7b9311d1d54
π_galaxy_err = parallax_err

# ╔═╡ 1cef7e19-53ab-4c8e-9613-113e86925b10
n_sigma_pi = (df.parallax .- π_galaxy .- π_zeropoint) ./ sqrt.(π_galaxy_err .^2 .+ df.parallax_error .^ 2)

# ╔═╡ 05c56974-6a0c-4d90-a903-7156b7fae3fe
missed_parallax = (df.F_CPAR .== 1) .!= (
	abs.(n_sigma_pi) .< 3
) .| (ismissing.(df.parallax))

# ╔═╡ a58ba633-8a6c-4ff0-b041-4b5a5c8d84c6
sum(missed_parallax)

# ╔═╡ 29d5536d-1b82-46d0-9b01-75f12579e196
df.F_CPAR[missed_parallax]

# ╔═╡ dd28242b-b8fa-4bae-8fc9-2d62ed5ab6cd
scatter(df.parallax[missed_parallax], df.parallax_error[missed_parallax])

# ╔═╡ 87578292-10b6-4022-b398-7adfeb22003f
df[missed_parallax, :].parallax_over_error

# ╔═╡ 666f9963-d75b-4543-84cc-120c6efab7de
TRGB_min .- 5*dm_err  

# ╔═╡ 2b1fd97b-2d31-4fd1-83a2-5eb26482544f
scatter((df.bp_rp )[flag_pm], (df.phot_g_mean_mag)[flag_pm], color=df.F_INGRID[flag_pm], markersize=1, axis=(; limits=(0.3, 2.6, 15.5, 22)))

# ╔═╡ ecf7e193-cacc-46f9-83b6-85a61333ec8c
TRGB_min

# ╔═╡ 63fde46a-2fcb-4e14-8a6a-0d61e4339679
md"""
# Dustmap
"""

# ╔═╡ 55bd8a2a-73c6-41fa-ba3a-857cf9c7327c
dustmaps = pyimport("dustmaps.sfd")

# ╔═╡ 34a183a1-e2d9-487f-83e4-cf539f102aac
dustmap = dustmaps.SFDQuery()

# ╔═╡ 478e5197-8969-4f59-a961-022de16c5e0a
SkyCoord = pyimport("astropy.coordinates").SkyCoord

# ╔═╡ 49168719-8a5b-4e14-8566-dff652abedcd
icrss = SkyCoord(df.ra, df.dec, unit="degree", frame="icrs")

# ╔═╡ eab76e65-3c36-42ae-8d16-71783f217f70
ebv = pyconvert(Vector, dustmap(icrss))

# ╔═╡ e47d4fbd-e036-4838-abc0-400a118faa61
A0=3.1*ebv

# ╔═╡ da7e92d8-4182-4656-ab45-5efda19ea0ff
Ag, Ab, Ar = let
	cg=[0.9761, -0.1704, 0.0086, 0.0011, -0.0438, 0.0013, 0.0099]
	cb=[1.1517, -0.0871, -0.0333, 0.0173, -0.0230, 0.0006, 0.0043]
	cr=[0.6104, -0.0170, -0.0026, -0.0017, -0.0078, 0.00005, 0.0006]
	terms = hcat(ones(length(df.xi)), df.bp_rp, (df.bp_rp) .^2.,(df.bp_rp) .^3., A0, A0 .^2., (df.bp_rp) .*A0)

	kg = terms * (cg)
	kb =  terms * cb
	kr = terms * cr

	A0 .*kg, A0 .*kb, A0 .*kr
end

# ╔═╡ 46118ebf-6537-4611-8703-cde3aeacc8ed
flag_cmd = (TRGB_min .- 5*dm_err  .< df.phot_g_mean_mag .- Ag) .& (df.phot_g_mean_mag .- Ag .< 22 ) .&
	(-0.5   .< df.bp_rp .- Ab .+ Ar .< 2.5  ) .& 
	.!ismissing(df.phot_g_mean_mag) .& .!ismissing.(df.bp_rp)

# ╔═╡ dd5d1b18-d37d-4c30-acbf-eae5078014c9
missed_grid = (flag_cmd .& flag_pm) .!= (df.F_INGRID .== 1.0)

# ╔═╡ fcf14f68-fbeb-4ca4-86d6-18fef427cda8
sum(skipmissing(missed_grid))

# ╔═╡ a3ad3081-0534-4aee-ad63-66ce6c635ad0
df[missed_grid, [:pmra, :pmdec, :phot_g_mean_mag, :bp_rp, :F_INGRID]]

# ╔═╡ b042228a-474a-4bab-aed7-fc34f4072cb9
scatter(df.bp_rp[missed_grid], df.phot_g_mean_mag[missed_grid], color=df.F_INGRID[missed_grid])

# ╔═╡ 156642e2-84e2-4c54-8fb5-46d071a89bf5
scatter(df.xi[missed_grid .& flag_cmd], df.eta[missed_grid .& flag_cmd] )

# ╔═╡ f295bee9-a1ec-4442-aa18-1da153a2241f
df[missed_grid .& flag_cmd, :].L_PM_SAT |> unique

# ╔═╡ 64443ab9-a442-409d-bbb0-dace080f2c42
scatter((df.bp_rp .- Ab .+ Ar)[flag_pm], (df.phot_g_mean_mag .- Ag)[flag_pm], color=df.F_INGRID[flag_pm], markersize=1, axis=(; limits=(0.3, 2.6, 15.5, 22)))

# ╔═╡ ab4e0800-f8bb-41fb-bd3e-4566ffecdac6
scatter(df.bp_rp[missed_grid .& .!flag_cmd], (df.phot_g_mean_mag .- Ag)[missed_grid .& .!flag_cmd], )

# ╔═╡ 787241a9-861b-4b71-ae22-fb8d637b96c0
hist(filter(x->!ismissing(x), Ag))

# ╔═╡ 7d1eb898-e5d8-4b6e-90fb-75878866adeb
ext = (Ab .- Ar)

# ╔═╡ 80491f2c-55bc-49a6-a312-883192e0303e
k_filt = .!ismissing.(Ab .- Ar)

# ╔═╡ 1bb6e702-d186-4231-bcdb-84d3d4c30246
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			  xlabel = L"\xi", 
			  ylabel = L"\eta",
		aspect=DataAspect()
			 )
	
	
	p = scatter!(df.xi[k_filt], df.eta[k_filt], color=Float64.((Ab .- Ar)[k_filt]))

	Colorbar(fig[1,2], p, label="ABP-ARP")

	Makie.save("dustmap_$galaxy.pdf", fig)
	fig
end

# ╔═╡ 04b6aa5c-5a60-4de6-a59b-c584e8eed6ba
BP = df.phot_bp_mean_mag .- Ab

# ╔═╡ a13e4a12-35a8-44da-802c-8259518be133
RP = df.phot_rp_mean_mag .- Ar

# ╔═╡ 51b5cdc1-0a14-407b-bfaa-ac912e6f0010
GP = df.phot_g_mean_mag .- Ag

# ╔═╡ 0aeff51d-ce79-40d3-ad17-c1e5afc18d60
md"""
# Flux excess
"""

# ╔═╡ b91e8566-ce3a-458c-9a8e-0ef440ef177d
md"""
Works excatly. Reillo 2021.
"""

# ╔═╡ a04d229c-5001-44ff-991e-99804fbd3a06
flag_excess = @. (abs.(Cc) < 3 *( 0.0059898 + 8.817481e-12 * df.phot_g_mean_mag .^ 7.618399)) .& (df.phot_g_mean_mag .> 4)

# ╔═╡ b8b32430-7ced-4349-9093-f08545d87043
missed_flux = (flag_excess .& .!ismissing.(flag_excess)) .!= (df.F_FLUXEXCESS .== 1)

# ╔═╡ 77d5b1ca-4f92-4191-9142-f66450f14eb5
@assert sum(missed_flux)==0 "flux excess filter wonk"

# ╔═╡ 6c4b0d0d-dab0-48c6-9c94-5b566fde6cd4
function flux_correction(x)
	if x === missing
		return NaN
	end
	if x < 0.5
		return 1.154360 + 0.033772*x + 0.032277*x^2
	elseif x < 4.0
		return 1.162004 + 0.011464*x + 0.049255*x^2 + -0.005879x^3
	else
		return 1.057572 + 0.140537*x
	end
end

# ╔═╡ c4f01609-86dd-4920-bb5b-5e2e8650d474
Cc = df.phot_bp_rp_excess_factor .- flux_correction.(df.bp_rp)

# ╔═╡ f118cf45-ed74-46d0-97bd-6b7118a1a633


# ╔═╡ 27cced0c-9f2b-443b-ab5d-e43bf349bce2
md"""
# RUWE filter
"""

# ╔═╡ 13cc34a7-8025-4acb-83b2-f38799a05a2c
md"""
Works exactly
"""

# ╔═╡ 6dddbddc-d080-46ef-91a4-3670eadaa8f8
flag_ruwe = df.ruwe .< 1.3

# ╔═╡ d2e08976-ec4f-4c7f-9d4a-9b260f05bc49
sum(skipmissing(flag_ruwe .!= (df.F_ASTROMETRIC .== 1)))

# ╔═╡ 927e1156-009a-4923-bf61-4d5f4d7cf942
md"""
# Isochrone
"""

# ╔═╡ 345a9ac3-95cf-4c3d-ab0c-822fa55cbd8e
import CSV

# ╔═╡ 1fef960a-84f2-49d4-8273-4961512a1be6
header = string.(split("Zini     MH   logAge Mini        int_IMF         Mass   logL    logTe  logg  label  mbolmag  Gmag    G_BPmag  G_RPmag", r"\s+"))

# ╔═╡ 29921cb9-5197-4104-95a7-9da92be8e3cc
md"""
Isochrone should lay along main track (it does!)
"""

# ╔═╡ 92c6bb19-4b54-467f-8b58-c2ea2119a551
iso = CSV.read("padova_mh-1.7.txt", DataFrame, ignorerepeated=true, delim=" ", comment="#", header=header)

# ╔═╡ 5b8cee3e-047e-4aaf-ab74-af97d8e0c84d
let
	fig = Figure()
	 ax = Axis(fig[1,1],
			   yreversed=true,
			   limits=(-0.5, 2.5, 15, 22)
			  )

	scatter!(df.bp_rp .- Ab .+ Ar, df.phot_g_mean_mag .- Ag, color=df.L_CMD_SAT)


	lines!(iso.G_BPmag .- iso.G_RPmag, iso.Gmag .+ dm)

	fig
end

# ╔═╡ Cell order:
# ╠═73596566-f0bf-4966-9e4b-a8d451ce307d
# ╠═944a654e-d955-4ea8-94ab-8aa5b73f9ea6
# ╠═3d9e96f5-7d9b-4f07-8ca8-6da1bf731e70
# ╠═8cff9c1e-1a19-11f0-2dd6-7bddd47fe438
# ╠═3b0e285e-cefd-42ff-8c0e-dbe71ef2b720
# ╠═061b5be8-fc49-4fc4-9c5d-239079740900
# ╠═31e6bedf-b347-4c83-89f6-c5c1adc7b481
# ╠═0653d025-c0e1-4b98-9701-f3db2f6ebe8d
# ╠═c3933f11-ec53-49b0-a7c6-fc34ed8ef5c6
# ╠═9beeec68-ca34-49a8-8f5b-6be279bacd13
# ╟─8515cbfb-1385-40f6-a53a-36183be38375
# ╟─ead04430-df91-42f2-9a66-9e47563411db
# ╠═4b12c940-5202-4903-9b5b-7598b11e4726
# ╠═f35de603-cf15-4a1b-9c00-f4c70308af60
# ╠═2d407500-5ed6-4f14-9a19-db68616202b8
# ╠═dd6a2a61-910f-4a14-afeb-ea6016125478
# ╟─a7e52514-c2f7-4969-ac5c-df9c60b1d4bc
# ╟─ba34ed2a-2a1c-494a-ba1a-8ce45b1d6a04
# ╠═7ffbdc53-5c3d-438c-8071-59e801590b6b
# ╠═1b1098fd-d7cc-437f-91d2-62c90fcd97e0
# ╠═69e3c284-cc10-4941-a6bf-655f3cbe79b0
# ╠═9666b674-88c2-4466-b9ce-f19682368362
# ╠═0898688e-d627-4174-a01c-192d0e7b2158
# ╠═ecbb1ee7-8ef6-4fc7-a612-c09d29a2aaed
# ╠═3eca57d2-3049-4f36-ad15-3219981629dd
# ╠═02836247-daef-4cd1-9fe9-e963a840132a
# ╟─2c344123-e54c-4a94-aa58-7eb9e305084e
# ╟─eccc30af-0447-4456-9026-16a3aa13f964
# ╠═4b086ff9-262e-400c-886c-37ea3cb2e2e8
# ╠═2da0de80-0409-4f5d-93d5-9757d64c233c
# ╠═874cf235-809f-43e4-9d3c-b7b9311d1d54
# ╠═1cef7e19-53ab-4c8e-9613-113e86925b10
# ╠═05c56974-6a0c-4d90-a903-7156b7fae3fe
# ╠═a58ba633-8a6c-4ff0-b041-4b5a5c8d84c6
# ╠═29d5536d-1b82-46d0-9b01-75f12579e196
# ╠═dd28242b-b8fa-4bae-8fc9-2d62ed5ab6cd
# ╠═87578292-10b6-4022-b398-7adfeb22003f
# ╟─1f6d3e45-b7e0-4b6f-915a-27fbb4a30ff5
# ╠═2f3fbd8e-96c4-40e1-91c0-d9e3bd9f2e1a
# ╠═67ab8518-3b43-4d28-840d-690febff5959
# ╠═ccc18cbd-43b3-4d4e-a308-727561f81c30
# ╠═6c512f14-93ac-468a-be8c-9479fcb3eca9
# ╠═bdc91b1f-b8a9-4150-b86a-00de4c25d519
# ╠═2be03779-adac-4d46-af22-2f54f1fec593
# ╠═c049e7c6-e6e8-4c91-b563-e34545b6a71f
# ╠═cfe395e1-ed32-4435-bedf-7d25a443762a
# ╠═666f9963-d75b-4543-84cc-120c6efab7de
# ╠═46118ebf-6537-4611-8703-cde3aeacc8ed
# ╠═dd5d1b18-d37d-4c30-acbf-eae5078014c9
# ╠═fcf14f68-fbeb-4ca4-86d6-18fef427cda8
# ╠═a3ad3081-0534-4aee-ad63-66ce6c635ad0
# ╠═b042228a-474a-4bab-aed7-fc34f4072cb9
# ╠═64443ab9-a442-409d-bbb0-dace080f2c42
# ╠═2b1fd97b-2d31-4fd1-83a2-5eb26482544f
# ╠═ecf7e193-cacc-46f9-83b6-85a61333ec8c
# ╠═ab4e0800-f8bb-41fb-bd3e-4566ffecdac6
# ╠═156642e2-84e2-4c54-8fb5-46d071a89bf5
# ╠═f295bee9-a1ec-4442-aa18-1da153a2241f
# ╠═63fde46a-2fcb-4e14-8a6a-0d61e4339679
# ╠═787241a9-861b-4b71-ae22-fb8d637b96c0
# ╠═6a1415e7-2cdc-4fe9-88cb-43563f2141f1
# ╠═55bd8a2a-73c6-41fa-ba3a-857cf9c7327c
# ╠═34a183a1-e2d9-487f-83e4-cf539f102aac
# ╠═478e5197-8969-4f59-a961-022de16c5e0a
# ╠═49168719-8a5b-4e14-8566-dff652abedcd
# ╠═eab76e65-3c36-42ae-8d16-71783f217f70
# ╠═e47d4fbd-e036-4838-abc0-400a118faa61
# ╠═da7e92d8-4182-4656-ab45-5efda19ea0ff
# ╠═7d1eb898-e5d8-4b6e-90fb-75878866adeb
# ╟─1bb6e702-d186-4231-bcdb-84d3d4c30246
# ╠═80491f2c-55bc-49a6-a312-883192e0303e
# ╠═04b6aa5c-5a60-4de6-a59b-c584e8eed6ba
# ╠═a13e4a12-35a8-44da-802c-8259518be133
# ╠═51b5cdc1-0a14-407b-bfaa-ac912e6f0010
# ╟─0aeff51d-ce79-40d3-ad17-c1e5afc18d60
# ╠═b91e8566-ce3a-458c-9a8e-0ef440ef177d
# ╠═b8b32430-7ced-4349-9093-f08545d87043
# ╠═77d5b1ca-4f92-4191-9142-f66450f14eb5
# ╠═a04d229c-5001-44ff-991e-99804fbd3a06
# ╠═6c4b0d0d-dab0-48c6-9c94-5b566fde6cd4
# ╠═c4f01609-86dd-4920-bb5b-5e2e8650d474
# ╠═f118cf45-ed74-46d0-97bd-6b7118a1a633
# ╟─27cced0c-9f2b-443b-ab5d-e43bf349bce2
# ╠═13cc34a7-8025-4acb-83b2-f38799a05a2c
# ╠═6dddbddc-d080-46ef-91a4-3670eadaa8f8
# ╠═d2e08976-ec4f-4c7f-9d4a-9b260f05bc49
# ╠═927e1156-009a-4923-bf61-4d5f4d7cf942
# ╠═345a9ac3-95cf-4c3d-ab0c-822fa55cbd8e
# ╠═25728329-5b72-4815-a37f-d24873fd408a
# ╠═1fef960a-84f2-49d4-8273-4961512a1be6
# ╟─29921cb9-5197-4104-95a7-9da92be8e3cc
# ╠═92c6bb19-4b54-467f-8b58-c2ea2119a551
# ╠═5b8cee3e-047e-4aaf-ab74-af97d8e0c84d
