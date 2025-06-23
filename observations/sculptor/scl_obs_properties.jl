### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 53f06974-1fcc-4c90-86a0-5dd0cec4e4b8
begin
	using Arya
	using CairoMakie

	using DataFrames
	using Measurements
end

# ╔═╡ de046588-5a2b-41d4-a8cc-90dac36318e6
using Unitful, UnitfulAstro

# ╔═╡ 3d32bef6-de64-4d1e-93a8-46c921c86011
using StatsBase: mean, std

# ╔═╡ 72eb810c-114e-11ef-30ea-fd46e923ea49
md"""
# Observations of Sculptor

A quick collection of past literature observations and measurments of the Sculptor DSph
"""

# ╔═╡ 722bd144-d047-4abf-b82f-f733134d3eb7
dm_to_d(dm) = 10 * 10^(dm / 5) / 1e3

# ╔═╡ 23022242-98d6-4509-bb60-b992c56a3bb0
md"""
# index
"""

# ╔═╡ ff2a1d65-c80f-4eb4-939f-f52e0e29db56


# ╔═╡ 26ae0d94-698c-4e9d-bac6-3e91d0d197ab
md"""
# References

"""

# ╔═╡ fef51bf1-0022-478c-ac62-0ffd127c1bd0
md"""
## Me
"""

# ╔═╡ d13c99eb-2cb2-4093-85bb-19a3dd18675d
myobs = Dict(
	:study => "this",
	:ra => 15.0209 ± 0.017, # gaia systematics
	:dec => -33.7127 ± 0.017,
	:pm_ra => 0.111299 ± 0.0078,
	:pm_dec => -0.148202 ± 0.00586,
	:σ_v => (9.61 ± 0.16)*u"km/s",
	:radial_velocity => (111.03 ± 0.23)*u"km/s",
)

# ╔═╡ 39880bc6-3e37-46df-b0d8-742afe684075
md"""
## Hodge ....
"""

# ╔═╡ 3e9a4212-8735-4483-b1fd-97e5d5196546
md"""
# 1982 ESO/Uppsala survey of the ESO(B) atlas

"""

# ╔═╡ 2392b8d8-a763-4efd-b0df-735466731001
uppsala = Dict(
	:study => "ESO82",
	:ra => 015.03899 ± 4e-5,
	:dec => -33.70903 ± 3e-5
)

# ╔═╡ e5d18175-f4c8-4da5-9246-fd58bcbd344e
md"""
## Webbink 1985
## Caldwell 1992
"""

# ╔═╡ 28941a38-1583-417d-aecf-5dfc0c4fccd8
md"""
## Pryor 1992
"""

# ╔═╡ f3525bab-9241-4484-9a0e-b1db0d986c07
md"""
## Armandroff and Da Costa 1986
Radial velocities on 16 stars

"""

# ╔═╡ ae53cae7-d547-4d1d-9101-b064ec6616cf
md"""
## Kunkel & Demers 1977
"""

# ╔═╡ 5a27c2b9-add0-40dc-93c6-2530a60f927d
ad1986 = Dict(
	:study => "A&DC86",
	:σ_v => (6.3 ± 1.2)*u"km/s",
	:radial_velocity => (107.4 ± 2)*u"km/s",
)

# ╔═╡ 346f542f-d7dc-490f-aadb-32ba109e277d
md"""
## Mateo 1993
"""

# ╔═╡ 3e476519-94ae-4463-94a2-1d26e0a7969b
md"""
## Irwin & Hatzidimitriou 1995
photometric structure of scolptor (and other dsph) by number density
"""

# ╔═╡ c522a3da-c394-4ab3-8a13-dcd1810a818f
iw1995 = Dict(
	:study => "I&W95",
	:ell => 0.32 ± 0.03,
	:PA => 99 ± 1, 
	:r_h => 6.8, # double check this one
	:L => (1.4 ± 0.6) * 1e6u"Lsun"
)

# ╔═╡ e8775a3e-afab-48c7-a781-1fcc877ddcc1
md"""
## battaglia + 2008
VLT/FLAMES of 470 candidates. Notices a velocity gradient along the major axis. Mass measurements from kinimatics.

- NFW: $c = 20$, $M(<1.8{\rm kpc}) = (2.2 \pm 1) \times 10^9 M_\odot$
- M/L < 1.8 kpc = $150 \pm 30$
"""

# ╔═╡ 36d90fee-8368-465f-baba-e714bf141efc
battaglia2008 = Dict(
	:study => "battaglia+08",
	:σ_v => (10.1 ± 0.3)*u"km/s",
	:radial_velocity => (110.5 ± 0.5)*u"km/s",
)

# ╔═╡ 8f9fe40c-dd18-437f-ae9e-f0a50fba8155
md"""
## Pietrzy´nski + 2008
RR Lyrae stars in J and K to determine distance. Get a distance modulus of 19.67 pm 0.02 pm 0.12 (sys).
Also mention Rizzi's 2002 thesis (not easy to find) which measures similar DM from TRGB (19.64pm0.08) and HB (19.66 pm 0.15)
"""

# ╔═╡ 750dd472-05e4-4c4b-b347-4f0581ec6f55
pietryznski2008 = Dict(
	:study => "P+08",
	:distance => dm_to_d(19.67 ± 0.12)
)

# ╔═╡ 7fa4817f-246d-4956-aee2-42ea577412f2
md"""
# Kirby + 2009
DEIMOS of 388 stars?
"""

# ╔═╡ 47e8e035-f0a0-4f68-9296-12b7939d83b2
kirby2009 = Dict(
	:study => "K+09",
	:radial_velocity => (111.6 ± 0.5) * u"km/s",
	:σ_r => (8 ± 0.7)
)

# ╔═╡ d90a79bc-0f53-4e70-aba5-0253624f6049
md"""
# Łokas + 2009
Dynamical analysis of Walker+2009 data.

Adopt distances from Mateo 98 and shape fits from Irwin & Hatzidimitriou (1995). 
Fit radial velocity dispersion profile and anisotropy using Jeans equations (?). 
"""

# ╔═╡ 1833bf84-f05d-461a-a84e-2f28efaffd17
lokas2009 = Dict(
	:study => "Ł+09",
	:β => −0.09 ± 0.2,
	:M_tot => (3.1 ± 0.2) * 1e7u"Msun"
)

# ╔═╡ 0ea649d4-c260-450f-a51c-7598da5bfe2e
md"""
## Walker + 2009

Actually two papers. Walker, Mateo & Olszewski presents observations of the stars in several dsph. 

Walker, Mateo, Olszewski, Sen & Woodroofe describe and apply an expectation-value maximum likelyhood to filter out contaminates and derive the radial velocities and velocity dispersions.
"""

# ╔═╡ 6eed53c4-7798-4ced-ae58-59679f4cb380
walker2009 =  Dict(
	:study => "W+09",
	:radial_velocity => (111.4 ± 0.1) * u"km/s",
	:σ_v => (9.2 ± 1.1) * u"km / s"
)

# ╔═╡ e71ff05c-bbad-4e1c-995c-1196bb5fb462
md"""
## Tully et al. (2013)
Large sample of galaxy distances
Updates earlier Tully work: For sculptor, distance = 80 kpc +- 0.08  Tip of RGB measurement + misc?
"""

# ╔═╡ d570bb54-e2d6-4eda-ae2e-f1c97cf5f401
80 * 0.08

# ╔═╡ 96ba64bc-90d2-4beb-8bde-f5e661051644
tully2013 = Dict(
	:study => "Tully+13",
	:distance => (80 ± 6.4)
)

# ╔═╡ 351731df-8d84-41ba-83e7-793898b9a148
md"""
## McChonnachie 2012
For sculptor distance: Pietrzy´nski 2008, Walker + 2009 for RV. 
"""

# ╔═╡ 0e1d8bc0-103a-4623-8546-aa7a32ee4504
md"""
## Kirby 2013
Metallicity, SFR ??
"""

# ╔═╡ b8055846-efa5-4bcc-a82e-e3b4cd40788d
md"""
## Martínez-Vázquez + 2015 
RR Lyrae distance modulus
"""

# ╔═╡ f7132d82-eb69-4e56-b61b-dd60e6cf8fdb
mv2015 = Dict(
	:study => "MV+15", 
	:distance => dm_to_d(19.62±0.04)
)

# ╔═╡ 7eae9400-72bf-4c96-95b7-3c335dddd40a
md"""
# Sohn + 2017
HST proper motions
"""

# ╔═╡ d7277867-73d6-486a-a622-d5a70cc89e85
sohn2017 = Dict(
	:study => "Sohn+17",
	:pm_ra => -0.0296 ± 0.0209,
	:pm_dec => -0.1358 ± 0.0214,
)

# ╔═╡ 73ed55db-165e-4fc4-9271-2ca4027a6268
md"""
## Strigari, Frenk, White 2017
"""

# ╔═╡ 7a0b8c50-2881-4c1f-8de1-10c865849980
md"""
## Garofalo + 2018
RRL distance estimate
- μ = 19.60 ± 0.02 (statistical) ± 0.04 (photometric) mag (with σsys = 0.09 mag)
"""

# ╔═╡ e2a4f7d9-dc24-4929-afb9-c1e9e2ead9e2
sqrt(0.02^2 + 0.04^2 + 0.09^2)

# ╔═╡ 98fe3c7a-2ee8-4600-a8d4-eacb4bd78054
garofalo2018 = Dict(
	:study => "G+18",
	:distance => dm_to_d(19.60 ± 0.1),
)

# ╔═╡ e078ea6f-9fd3-4bf5-a3ae-665a8d64e046
md"""
## Massari + 2018
Also proper motion analysis of HST + Gaia. Find that for an orbit  r = 73+8 kpc and r = 222+170 kpc. 

σR = 11.5 ± 4.3 km s−1 and σT = 8.5 ± 3.2 km s−1
- β~0.86+0.12-0.83

- absolute proper motions: (μα cos(δ), μδ )=(−0.20±0.14, −0.33±0.11) masyr−1. References several older studies...

"""

# ╔═╡ 5b435aac-7fcd-403f-8388-5f635bbadd9a
massari18 = Dict(
	:study => "Massari+18",
	:pm_ra => 0.1615 ± 0.14,
	:pm_dec => -0.805 ± 0.11
)

# ╔═╡ 7a5e54c2-6186-4975-9cd4-b971a982a8cf
md"""
## Muñoz + 2018

Uses megacam to do photometric survey on dwarfs and globular clusters. Solve for ra, dec, r_h, PA, ecc, and sigma0, and background for each satalite using several different profiles. Also derive total magnitudes.

Absolute magnitudes are derived from star counts + a theoretical luminosity function.
"""

# ╔═╡ ac92bf08-1a86-4964-951a-97c8f738ed1f
muñoz2018 = Dict(
	:study => "Muñoz+18",
	:ra => 15.0183 ± 0.30 / 3600 * 360 / 24, # seconds to arcseconds
	:dec => -33.7186 ± 2.6 / 3600,
	:L => 10 .^ (6.262 ± 0.056) * u"Lsun", # derived from Mv=−10.82±0.14
	:PA => 92 ± 1, # for exponential, includes several others
	:ell => 0.36 ± 0.01,
	:r_h => 12.43 ± 0.18
)

# ╔═╡ 38e9075e-f827-4a0e-9e4b-06dd5b21f70e
2.6 / 3600 / 360

# ╔═╡ 7fdfd324-950e-4c81-9a43-c35e748594ac
0.30 / 3600 / 24 * cosd(-33.7182) # corrected ra error

# ╔═╡ 202abcd8-1b8a-4f18-aca5-7b92fca82a56
md"""
## Strigari, Frenk, White 2018
Gaia proper motion analysis + HST
 Measure tangential and radial dispersion..., in agreement.
In general in agreement with Walker & Penarrubia 2011
"""

# ╔═╡ 28fe217f-9bfb-450a-9e58-f552adb580d7
md"""
## McChonnachie & Venn 2020 (&2020a)

The B versiun is just an update from DR2 to Gaia EDR3.

Other properties besides proper motions and profiles are simply updated from McChonnachie 2012. I use the Gaia systematics from Lindegren+2018 (DR2) or Lindegren+2021 (EDR3).
"""

# ╔═╡ ab329391-9ed0-44f8-bf2b-520c87dab362
mv2020 = Dict(
	:study => "MV2020", 
	:pm_ra => 0.081±0.027,
	:pm_dec => -0.136±0.027,
	:ra => 15.0392,
	:dec => -33.7092
)

# ╔═╡ cde8a009-958e-4511-be4f-b903913e6b49
mv2020a = Dict(
	:study => "MV2020a", 
	:pm_ra => 0.099±0.017,
	:pm_dec => -0.160±0.017
)

# ╔═╡ 6a7761ba-cfe2-40e7-8d2d-78f7a47ff164
md"""
## Pace & Li
"""

# ╔═╡ 8a1832ae-6335-44c8-bc56-d8a05e3e301d
md"""
## Battaglia + 2022
references for basic properties: Battaglia+2008, Martinez-Vazquez 2015, Muñoz + 2018

With orbits in the perturbed (LMC) potential, find peries between44.3 to 51.1, and apos between 274.3 and 552.4m. The LMC is likely very important for the evolution of sculptor
"""

# ╔═╡ 30f6a515-76dc-47fc-abc0-a90e28f64d9f
battaglia2022 = Dict(
	:study => "B+22",
	:ra => 15.01830,
	:dec => -33.71860,
	:pm_ra => 0.099±0.017, # gaia systematic err
	:pm_dec => -0.159 ± 0.017
)

# ╔═╡ dcf1884d-9f8a-4ff9-9437-acce8e447709
md"""
## Pace + 2022

7k members

"""

# ╔═╡ 74531f13-33b5-4821-8991-13c57a44d956
pace2022 = Dict(
	:study => "Pace+22",
	:pm_ra => 0.100 ± 0.017,
	:pm_dec => -0.158 ± 0.017 # gaia sys err
)

# ╔═╡ 687c621a-1103-4fbb-bb88-17162c6c30c9
md"""
# Tran + 2022
Distance estimates from RRL, TRGB, and HB
- RRL = 19.60 pm 0.01 pm 0.05
- TRGB = 19.59 pm 0.07 pm 0.05
- HB = 19.54 pm 0.03 pm 0.09
"""

# ╔═╡ 0bf2d300-82ca-4ff4-a57e-4f50fb319605
sqrt(0.01^2 + 0.05^2)

# ╔═╡ 93207859-349a-426f-b77e-c310e1fee59c
tran2022 = Dict(
	:study => "tran+22",
	:distance => dm_to_d(19.60 ± 0.051)
)

# ╔═╡ c36ea4b5-323a-451f-a620-961ab8456e11
md"""
## del Pino + 2023

Rotation study
"""

# ╔═╡ cb3ef5d8-8231-42f1-b54f-42a6b08694c1
md"""
## Tolstoy+2023
VLT/FLAMES spectoscopic survey with Gaia DR3 proper motions (and parallax cuts)
Contains a large sample of stars (1701 incdividual stars). Builds on Battaglia 2022.
Claim that there is a gradient (decreasing) of velocity dispersion and radial velocity (increasing) in the outer regions. Do not see a rotation signature.

fe_h=-1.82 ± 0.45,

"""

# ╔═╡ 1ecf4404-caa5-418d-9170-2d539df69275
tolstoy2023 = Dict(
	:study => "T+23",
	:radial_velocity => (111.2 ± 0.25) * u"km/s",
	:pm_ra => 0.097 ± 0.017, # gaia systematic err
	:pm_dec => -0.148 ± 0.017 # is proper?
)

# ╔═╡ 777ec196-e193-456d-8d44-cba200a366dd
obs = [
	myobs,
	uppsala,
	ad1986,
	iw1995,
	battaglia2008,
	pietryznski2008, 
	kirby2009,
	walker2009, 
	lokas2009,
	tully2013,
	mv2015, 
	sohn2017,
	muñoz2018,
	massari18, #not comparable to other proper motion estimates
	garofalo2018,
	mv2020, 
	mv2020a, 
	pace2022, 
	battaglia2022,
	tran2022,
	tolstoy2023,

]

# ╔═╡ 5540a9a2-54b6-42cf-ad68-470cdc9fc667
md"""
# Comparisons
"""

# ╔═╡ 4a0bd7c4-9b37-4625-a31c-2f4f85087eed
function get_properties(obs, key)
	filt = haskey.(obs, key)
	values = [v[key] for v in obs[filt]]
	studies = string.([v[:study] for v in obs[filt]])
	return studies, values
end

# ╔═╡ 5928daad-d96d-4228-8ec3-0315c1c3cf2d
Arya.value(a::Measurement) = a.val

# ╔═╡ 684841fb-4da5-4e9a-89e9-01b425feae5e
Arya.err(a::Measurement) = a.err

# ╔═╡ 30dd27e2-fe7b-4695-99d3-0ad0b769628f
function compare_measurements(key, label; units=1, kwargs...)
	fig = Figure()

	x, y = get_properties(obs, key)
	N = length(x)
	y = y / units
	println(y)

	xt = collect(1:N)
	
	ax = Axis(fig[1,1],
		xticks =(xt, x), 
		xminorticksvisible=false,
		xticklabelrotation=-0π/6,
		ylabel=label, 
		kwargs...
	)

	tight_xticklabel_spacing!(ax)

	errscatter!(xt, Arya.value.(y), yerr=Arya.err.(y))

	
	fig

end

# ╔═╡ bc848e9f-db31-457c-80d5-65cf5397e906
let
	fig = Figure()
 
	ax = Axis(fig[1,1], 
		xlabel="RA / degrees", ylabel="Dec / degrees",
		aspect=1/cosd(-33.72)
	)

	study, ra = get_properties(obs, :ra)
	study, dec = get_properties(obs, :dec)
	
	N = length(ra)
	for i in 1:N
		x = [ra[i]]
		y = [dec[i]]
		errscatter!(ax, x, y, 
			yerr=Arya.err.(y), xerr=Arya.err.(x),
			color=Arya.COLORS[i], label=study[i])
	end

	axislegend(ax, position=:lt)
	fig
end

# ╔═╡ 5f28c6a9-0efd-4afd-9881-a3eaefc98f35
let
	fig = Figure()
 
	ax = Axis(fig[1,1], 
		xlabel=L"\mu_{\alpha*}\;/\;\textrm{mas\,yr^{-1}}", ylabel=L"\mu_\delta\;/\;\textrm{mas\,yr^{-1}}"
	)

	study, pmra = get_properties(obs, :pm_ra)
	study, pmdec = get_properties(obs, :pm_dec)

	N = length(pmra)
	for i in 1:N
		x = [pmra[i]]
		y = [pmdec[i]]
		errscatter!(ax, x, y, 
			yerr=Arya.err.(y), xerr=Arya.err.(x),
			color=Arya.COLORS[i % 7 + 1], label=study[i])
	end

	axislegend(ax)
	fig
end

# ╔═╡ 4ba652e8-2d24-4859-9305-06d9db40c40c
compare_measurements(:distance, "heliocentric distance / kpc")

# ╔═╡ 60639986-1ffc-400f-9fcd-e44ad33cb033
compare_measurements(:radial_velocity, "radial velocity", units=u"km/s")

# ╔═╡ a3c6615c-fb23-4692-a04b-f4e327faf9cb
compare_measurements(:σ_v, "velocity dispersion", units=u"km/s")

# ╔═╡ efac8358-96c8-4204-85e4-00ad25fae1a4
compare_measurements(:L, "luminosity", units=u"Lsun")

# ╔═╡ 23329577-3769-4f63-af49-c4d95ee1f476
compare_measurements(:ell, "ellipticity")

# ╔═╡ 3baa83fc-932d-4695-bd93-de162d3ef5df
compare_measurements(:PA, "PA / degrees")

# ╔═╡ 0457b20e-4b0d-4ec4-99fe-aac076361bf4
compare_measurements(:r_h, "rh")

# ╔═╡ Cell order:
# ╟─72eb810c-114e-11ef-30ea-fd46e923ea49
# ╠═53f06974-1fcc-4c90-86a0-5dd0cec4e4b8
# ╠═de046588-5a2b-41d4-a8cc-90dac36318e6
# ╠═3d32bef6-de64-4d1e-93a8-46c921c86011
# ╠═722bd144-d047-4abf-b82f-f733134d3eb7
# ╟─23022242-98d6-4509-bb60-b992c56a3bb0
# ╠═777ec196-e193-456d-8d44-cba200a366dd
# ╠═ff2a1d65-c80f-4eb4-939f-f52e0e29db56
# ╟─26ae0d94-698c-4e9d-bac6-3e91d0d197ab
# ╠═fef51bf1-0022-478c-ac62-0ffd127c1bd0
# ╠═d13c99eb-2cb2-4093-85bb-19a3dd18675d
# ╟─39880bc6-3e37-46df-b0d8-742afe684075
# ╠═3e9a4212-8735-4483-b1fd-97e5d5196546
# ╠═2392b8d8-a763-4efd-b0df-735466731001
# ╟─e5d18175-f4c8-4da5-9246-fd58bcbd344e
# ╟─28941a38-1583-417d-aecf-5dfc0c4fccd8
# ╟─f3525bab-9241-4484-9a0e-b1db0d986c07
# ╠═ae53cae7-d547-4d1d-9101-b064ec6616cf
# ╠═5a27c2b9-add0-40dc-93c6-2530a60f927d
# ╠═346f542f-d7dc-490f-aadb-32ba109e277d
# ╠═3e476519-94ae-4463-94a2-1d26e0a7969b
# ╠═c522a3da-c394-4ab3-8a13-dcd1810a818f
# ╟─e8775a3e-afab-48c7-a781-1fcc877ddcc1
# ╠═36d90fee-8368-465f-baba-e714bf141efc
# ╟─8f9fe40c-dd18-437f-ae9e-f0a50fba8155
# ╠═750dd472-05e4-4c4b-b347-4f0581ec6f55
# ╟─7fa4817f-246d-4956-aee2-42ea577412f2
# ╠═47e8e035-f0a0-4f68-9296-12b7939d83b2
# ╠═d90a79bc-0f53-4e70-aba5-0253624f6049
# ╠═1833bf84-f05d-461a-a84e-2f28efaffd17
# ╠═0ea649d4-c260-450f-a51c-7598da5bfe2e
# ╠═6eed53c4-7798-4ced-ae58-59679f4cb380
# ╟─e71ff05c-bbad-4e1c-995c-1196bb5fb462
# ╠═d570bb54-e2d6-4eda-ae2e-f1c97cf5f401
# ╠═96ba64bc-90d2-4beb-8bde-f5e661051644
# ╠═351731df-8d84-41ba-83e7-793898b9a148
# ╟─0e1d8bc0-103a-4623-8546-aa7a32ee4504
# ╟─b8055846-efa5-4bcc-a82e-e3b4cd40788d
# ╠═f7132d82-eb69-4e56-b61b-dd60e6cf8fdb
# ╠═7eae9400-72bf-4c96-95b7-3c335dddd40a
# ╠═d7277867-73d6-486a-a622-d5a70cc89e85
# ╠═73ed55db-165e-4fc4-9271-2ca4027a6268
# ╠═7a0b8c50-2881-4c1f-8de1-10c865849980
# ╠═e2a4f7d9-dc24-4929-afb9-c1e9e2ead9e2
# ╠═98fe3c7a-2ee8-4600-a8d4-eacb4bd78054
# ╟─e078ea6f-9fd3-4bf5-a3ae-665a8d64e046
# ╠═5b435aac-7fcd-403f-8388-5f635bbadd9a
# ╟─7a5e54c2-6186-4975-9cd4-b971a982a8cf
# ╠═ac92bf08-1a86-4964-951a-97c8f738ed1f
# ╠═38e9075e-f827-4a0e-9e4b-06dd5b21f70e
# ╠═7fdfd324-950e-4c81-9a43-c35e748594ac
# ╟─202abcd8-1b8a-4f18-aca5-7b92fca82a56
# ╠═28fe217f-9bfb-450a-9e58-f552adb580d7
# ╠═ab329391-9ed0-44f8-bf2b-520c87dab362
# ╠═cde8a009-958e-4511-be4f-b903913e6b49
# ╟─6a7761ba-cfe2-40e7-8d2d-78f7a47ff164
# ╟─8a1832ae-6335-44c8-bc56-d8a05e3e301d
# ╠═30f6a515-76dc-47fc-abc0-a90e28f64d9f
# ╟─dcf1884d-9f8a-4ff9-9437-acce8e447709
# ╠═74531f13-33b5-4821-8991-13c57a44d956
# ╠═687c621a-1103-4fbb-bb88-17162c6c30c9
# ╠═0bf2d300-82ca-4ff4-a57e-4f50fb319605
# ╠═93207859-349a-426f-b77e-c310e1fee59c
# ╟─c36ea4b5-323a-451f-a620-961ab8456e11
# ╟─cb3ef5d8-8231-42f1-b54f-42a6b08694c1
# ╠═1ecf4404-caa5-418d-9170-2d539df69275
# ╟─5540a9a2-54b6-42cf-ad68-470cdc9fc667
# ╠═4a0bd7c4-9b37-4625-a31c-2f4f85087eed
# ╠═30dd27e2-fe7b-4695-99d3-0ad0b769628f
# ╠═5928daad-d96d-4228-8ec3-0315c1c3cf2d
# ╠═684841fb-4da5-4e9a-89e9-01b425feae5e
# ╟─bc848e9f-db31-457c-80d5-65cf5397e906
# ╟─5f28c6a9-0efd-4afd-9881-a3eaefc98f35
# ╠═4ba652e8-2d24-4859-9305-06d9db40c40c
# ╠═60639986-1ffc-400f-9fcd-e44ad33cb033
# ╠═a3c6615c-fb23-4692-a04b-f4e327faf9cb
# ╠═efac8358-96c8-4204-85e4-00ad25fae1a4
# ╠═23329577-3769-4f63-af49-c4d95ee1f476
# ╠═3baa83fc-932d-4695-bd93-de162d3ef5df
# ╠═0457b20e-4b0d-4ec4-99fe-aac076361bf4
