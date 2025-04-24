### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# ╔═╡ 944a654e-d955-4ea8-94ab-8aa5b73f9ea6
import Pkg; Pkg.activate()

# ╔═╡ 8cff9c1e-1a19-11f0-2dd6-7bddd47fe438
using PyFITS

# ╔═╡ 3b0e285e-cefd-42ff-8c0e-dbe71ef2b720
using CairoMakie

# ╔═╡ c8f9d3d9-8d91-4e63-8ba1-dc2139f9dca4
using DustExtinction

# ╔═╡ 061b5be8-fc49-4fc4-9c5d-239079740900
import TOML

# ╔═╡ 0653d025-c0e1-4b98-9701-f3db2f6ebe8d
galaxy = "ursa_minor"

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

# ╔═╡ 4b086ff9-262e-400c-886c-37ea3cb2e2e8
π_zeropoint = -0.029

# ╔═╡ 2da0de80-0409-4f5d-93d5-9757d64c233c
π_galaxy = 1 / 80

# ╔═╡ 874cf235-809f-43e4-9d3c-b7b9311d1d54
π_galaxy_err = obs_props["distance_err"] ./ obs_props["distance"] .^ 2

# ╔═╡ 1cef7e19-53ab-4c8e-9613-113e86925b10
n_sigma_pi = (df.parallax .- π_galaxy .- π_zeropoint) ./ sqrt.(π_galaxy_err .^2 .+ df.parallax_error .^ 2)

# ╔═╡ 05c56974-6a0c-4d90-a903-7156b7fae3fe
missed_parallax = (df.F_CPAR .== 1) .!= (
	abs.(n_sigma_pi) .< 3
) .| (ismissing.(df.parallax))

# ╔═╡ a58ba633-8a6c-4ff0-b041-4b5a5c8d84c6
sum(missed_parallax)

# ╔═╡ dd28242b-b8fa-4bae-8fc9-2d62ed5ab6cd
scatter(df.parallax[missed_parallax], df.parallax_error[missed_parallax])

# ╔═╡ 87578292-10b6-4022-b398-7adfeb22003f
df[missed_parallax, :].parallax_over_error

# ╔═╡ 1f6d3e45-b7e0-4b6f-915a-27fbb4a30ff5
md"""
# Grid filter
"""

# ╔═╡ 2f3fbd8e-96c4-40e1-91c0-d9e3bd9f2e1a
flag_pm = (abs.(df.pmra) .< 10) .& (abs.(df.pmdec) .< 10) .&
	.!ismissing.(df.pmra) .& .!ismissing.(df.pmdec)

# ╔═╡ 6c512f14-93ac-468a-be8c-9479fcb3eca9
TRGB_min = minimum(df[(df.PSAT .> 0.2) .& (df.F_BEST .== 1), :phot_g_mean_mag]) - 0.2

# ╔═╡ 46118ebf-6537-4611-8703-cde3aeacc8ed
flag_cmd = (TRGB_min .- df.dG  .< df.phot_g_mean_mag .< 22 .+ df.dG) .&
	(-0.5 .- df.dBP/2  .< df.bp_rp .< 2.5  .+ df.dBP/2 ) .& 
	.!ismissing(df.phot_g_mean_mag) .& .!ismissing.(df.bp_rp)

# ╔═╡ dd5d1b18-d37d-4c30-acbf-eae5078014c9
missed_grid = (flag_cmd .& flag_pm) .!= (df.F_INGRID .== 1.0)

# ╔═╡ fcf14f68-fbeb-4ca4-86d6-18fef427cda8
sum(skipmissing(missed_grid))

# ╔═╡ a3ad3081-0534-4aee-ad63-66ce6c635ad0
df[missed_grid, [:pmra, :pmdec, :phot_g_mean_mag, :bp_rp, :F_INGRID]]

# ╔═╡ b042228a-474a-4bab-aed7-fc34f4072cb9
scatter(df.bp_rp[missed_grid], df.phot_g_mean_mag[missed_grid], )

# ╔═╡ 63fde46a-2fcb-4e14-8a6a-0d61e4339679
md"""
# Excess filter
"""

# ╔═╡ a0fa0916-57a9-49dd-b173-3de5166d8701
dustmap = SFD98Map()

# ╔═╡ 0121ad8c-6ef5-4b59-98e6-3d77ee92e0e7
import SkyCoords

# ╔═╡ 49168719-8a5b-4e14-8566-dff652abedcd
icrss = @. SkyCoords.ICRSCoords.(deg2rad(df.ra), deg2rad(df.dec))

# ╔═╡ 6e6eae76-6faf-4997-b151-20f05d40ca3e
galcs = SkyCoords.convert.(SkyCoords.GalCoords, icrss)

# ╔═╡ 52c460f7-482d-4c44-9b68-f61a31ba8c23
icrsss2 = SkyCoords.convert.(SkyCoords.ICRSCoords, galcs)

# ╔═╡ 10232bcb-d8f7-41fc-992d-193259ab9790
gall = [g.l for g in galcs]

# ╔═╡ 99d77ae0-fc36-450b-b1e6-a3450b69e199
galb = [g.b for g in galcs]

# ╔═╡ eab76e65-3c36-42ae-8d16-71783f217f70
ebv = @. dustmap(gall, galb)

# ╔═╡ e47d4fbd-e036-4838-abc0-400a118faa61
A0=3.1*ebv

# ╔═╡ 3d8a2406-b632-407c-9171-51660b1086d9
import LinearAlgebra: ⋅

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

# ╔═╡ 307c9958-863e-48a4-aba7-95addab99371
ebv ./ 3.1

# ╔═╡ ae7f0bf3-7659-411d-b03c-51b3eedcbf0b
[1 2 3
 4 5 6] * [1,2,3]

# ╔═╡ ac5c8687-9e1d-4d4f-92fd-f0b371ac1a6d
df.dec

# ╔═╡ 7d1eb898-e5d8-4b6e-90fb-75878866adeb
ext = (Ab .- Ar)

# ╔═╡ a931f2a3-ca28-4de0-80cd-56d93ab62cc7
df[:, [:ra, :dec, :phot_bp_mean_mag, :phot_g_mean_mag, :phot_rp_mean_mag]]

# ╔═╡ 80491f2c-55bc-49a6-a312-883192e0303e
k_filt = .!ismissing.(Ab .- Ar)

# ╔═╡ 1bb6e702-d186-4231-bcdb-84d3d4c30246
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			  xlabel = "xi", 
			  ylabel = "eta",
		aspect=DataAspect()
			 )
	
	
	p = scatter!(df.xi[k_filt], df.eta[k_filt], color=Float64.((Ab .- Ar)[k_filt]))

	Colorbar(fig[1,2], p, label="ABP-ARP")

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

# ╔═╡ a04d229c-5001-44ff-991e-99804fbd3a06
flag_excess = @. (abs.(Cc) < 3 *( 0.0059898 + 8.817481e-12 * df.phot_g_mean_mag .^ 7.618399)) .& (df.phot_g_mean_mag .> 4)

# ╔═╡ b8b32430-7ced-4349-9093-f08545d87043
missed_flux = (flag_excess .& .!ismissing.(flag_excess)) .!= (df.F_FLUXEXCESS .== 1)

# ╔═╡ 77d5b1ca-4f92-4191-9142-f66450f14eb5
sum(missed_flux)

# ╔═╡ ec19db9e-7479-4a86-83bb-fdb9fb1794e0
scatter(((BP .- RP)) .^ 2, Cc , color=df.F_FLUXEXCESS, 
	   axis = (;limits=(nothing, nothing,-0.5, 5)), markersize=1)

# ╔═╡ 27cced0c-9f2b-443b-ab5d-e43bf349bce2
md"""
# RUWE filter
"""

# ╔═╡ 6dddbddc-d080-46ef-91a4-3670eadaa8f8
flag_ruwe = df.ruwe .< 1.3

# ╔═╡ d2e08976-ec4f-4c7f-9d4a-9b260f05bc49
sum(skipmissing(flag_ruwe .!= (df.F_ASTROMETRIC .== 1)))

# ╔═╡ 66181e0b-fa47-4063-8116-f1dad80a1f59
df[filt_missed, :]

# ╔═╡ Cell order:
# ╠═944a654e-d955-4ea8-94ab-8aa5b73f9ea6
# ╠═8cff9c1e-1a19-11f0-2dd6-7bddd47fe438
# ╠═3b0e285e-cefd-42ff-8c0e-dbe71ef2b720
# ╠═061b5be8-fc49-4fc4-9c5d-239079740900
# ╠═0653d025-c0e1-4b98-9701-f3db2f6ebe8d
# ╠═c3933f11-ec53-49b0-a7c6-fc34ed8ef5c6
# ╠═9beeec68-ca34-49a8-8f5b-6be279bacd13
# ╟─8515cbfb-1385-40f6-a53a-36183be38375
# ╟─ead04430-df91-42f2-9a66-9e47563411db
# ╠═f35de603-cf15-4a1b-9c00-f4c70308af60
# ╠═2d407500-5ed6-4f14-9a19-db68616202b8
# ╠═dd6a2a61-910f-4a14-afeb-ea6016125478
# ╟─a7e52514-c2f7-4969-ac5c-df9c60b1d4bc
# ╟─ba34ed2a-2a1c-494a-ba1a-8ce45b1d6a04
# ╠═3eca57d2-3049-4f36-ad15-3219981629dd
# ╠═02836247-daef-4cd1-9fe9-e963a840132a
# ╟─2c344123-e54c-4a94-aa58-7eb9e305084e
# ╠═4b086ff9-262e-400c-886c-37ea3cb2e2e8
# ╠═2da0de80-0409-4f5d-93d5-9757d64c233c
# ╠═874cf235-809f-43e4-9d3c-b7b9311d1d54
# ╠═1cef7e19-53ab-4c8e-9613-113e86925b10
# ╠═05c56974-6a0c-4d90-a903-7156b7fae3fe
# ╠═a58ba633-8a6c-4ff0-b041-4b5a5c8d84c6
# ╠═dd28242b-b8fa-4bae-8fc9-2d62ed5ab6cd
# ╠═87578292-10b6-4022-b398-7adfeb22003f
# ╠═1f6d3e45-b7e0-4b6f-915a-27fbb4a30ff5
# ╠═2f3fbd8e-96c4-40e1-91c0-d9e3bd9f2e1a
# ╠═6c512f14-93ac-468a-be8c-9479fcb3eca9
# ╠═46118ebf-6537-4611-8703-cde3aeacc8ed
# ╠═dd5d1b18-d37d-4c30-acbf-eae5078014c9
# ╠═fcf14f68-fbeb-4ca4-86d6-18fef427cda8
# ╠═a3ad3081-0534-4aee-ad63-66ce6c635ad0
# ╠═b042228a-474a-4bab-aed7-fc34f4072cb9
# ╟─63fde46a-2fcb-4e14-8a6a-0d61e4339679
# ╠═c8f9d3d9-8d91-4e63-8ba1-dc2139f9dca4
# ╠═a0fa0916-57a9-49dd-b173-3de5166d8701
# ╠═0121ad8c-6ef5-4b59-98e6-3d77ee92e0e7
# ╠═49168719-8a5b-4e14-8566-dff652abedcd
# ╠═6e6eae76-6faf-4997-b151-20f05d40ca3e
# ╠═52c460f7-482d-4c44-9b68-f61a31ba8c23
# ╠═10232bcb-d8f7-41fc-992d-193259ab9790
# ╠═99d77ae0-fc36-450b-b1e6-a3450b69e199
# ╠═eab76e65-3c36-42ae-8d16-71783f217f70
# ╠═e47d4fbd-e036-4838-abc0-400a118faa61
# ╠═3d8a2406-b632-407c-9171-51660b1086d9
# ╠═da7e92d8-4182-4656-ab45-5efda19ea0ff
# ╠═307c9958-863e-48a4-aba7-95addab99371
# ╠═ae7f0bf3-7659-411d-b03c-51b3eedcbf0b
# ╠═ac5c8687-9e1d-4d4f-92fd-f0b371ac1a6d
# ╠═7d1eb898-e5d8-4b6e-90fb-75878866adeb
# ╠═1bb6e702-d186-4231-bcdb-84d3d4c30246
# ╠═a931f2a3-ca28-4de0-80cd-56d93ab62cc7
# ╠═80491f2c-55bc-49a6-a312-883192e0303e
# ╠═04b6aa5c-5a60-4de6-a59b-c584e8eed6ba
# ╠═a13e4a12-35a8-44da-802c-8259518be133
# ╠═51b5cdc1-0a14-407b-bfaa-ac912e6f0010
# ╠═0aeff51d-ce79-40d3-ad17-c1e5afc18d60
# ╠═b8b32430-7ced-4349-9093-f08545d87043
# ╠═77d5b1ca-4f92-4191-9142-f66450f14eb5
# ╠═a04d229c-5001-44ff-991e-99804fbd3a06
# ╠═6c4b0d0d-dab0-48c6-9c94-5b566fde6cd4
# ╠═c4f01609-86dd-4920-bb5b-5e2e8650d474
# ╠═ec19db9e-7479-4a86-83bb-fdb9fb1794e0
# ╟─27cced0c-9f2b-443b-ab5d-e43bf349bce2
# ╠═6dddbddc-d080-46ef-91a4-3670eadaa8f8
# ╠═d2e08976-ec4f-4c7f-9d4a-9b260f05bc49
# ╠═66181e0b-fa47-4063-8116-f1dad80a1f59
