### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	
	using CSV, DataFrames
	using Arya
	using CairoMakie
	#using FITSIO
	using Measurements

	import LilGuys as lguys
end

# ╔═╡ 202e0f8b-b417-4597-a737-7c60c0575fd3
using FITSIO

# ╔═╡ d4123ffd-cb32-486e-a93d-f48d7112a831
include("filter_utils.jl")

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, kurtosis, sem

# ╔═╡ dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
not = !

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "../../data/"

# ╔═╡ da5a3b57-72bc-46e1-b1a0-6c02eb101626
function sigma_clip(x, nσ=5)

	dN = 1
	filt = not.(isnan.(x))
	while dN > 0
		println("iteration, $dN")
		σ = std(x[filt])
		μ = median(x[filt])
		N = sum(filt)
		filt .&= x .> μ - nσ * σ
		filt .&= x .< μ + nσ * σ
		dN = N - sum(filt)
	end
	
	return filt
end

# ╔═╡ 3d8b3c55-153b-4a4a-9289-78f34df07abc
"""
Cross matches 2 dataframes given angular seperation in arcmin
"""
function xmatch(df1::DataFrame, df2::DataFrame, max_sep=2)
	max_sep = max_sep / 3600
	dists = lguys.angular_distance.(df1.ra, df1.dec, df2.ra', df2.dec')

	idxs = [i[2] for i in dropdims(argmin(dists, dims=2), dims=2)]

	filt = dropdims(minimum(dists, dims=2), dims=2) .< max_sep
	return filt, idxs
end

# ╔═╡ 248c2d5f-85cc-44be-bc63-43d4af470182
begin 
	params = read_file("sculptor/fiducial.toml")
	params["filename"] = "$data_dir/Sculptor.GAIASOURCE.RUWE.VELS.PROB.fits"
	params["PA"] = 92
	params["ellipticity"] = 0.36
	params["PSAT_min"] = nothing
	params = DensityParams(params)
end

# ╔═╡ dfd2fbee-9993-4553-b693-6fb71a7b11a2
FITS

# ╔═╡ 7a50a176-96a5-4098-88d6-0fa2874d0f90
begin #TODO: implement J24 crossmatching
	j24 = load_stars(params.filename, params) # does not filter
end

# ╔═╡ 6f2359d8-332b-11ef-0db9-f1f06474c561
md"""
From Fed

Here attached some tables that I used for my paper on Sculptor. The literature compilations are from the DART survey and APOGEE DR17. There are some member stars in common between the two compilations. Therefore, I have divided into 4 files, I guess that depending on your needs, they might give you the same results.  I think all of them has been crossmatched with the Jaclyn’s member list and probabilities. I selected stars with a probability of being members of > 40 percent.

- The `Sculptor_DARTS_BEST_psat40.csv` contains all the DART members. 
- The table `Sculptor_DARTS_BEST_psat40_notAPOGEE` is the same but without the stars that are also in APOGEE. 

- `sculptor_apogeeDR17_xmatch.csv` has all the APOGEE stars,
- while DARt stars are removed in `sculptor_apogeeDR17_xmatch_notDART.csv`. 

The total members from DART and APOGEE (removing duplicates) should be 617 stars, either from `Sculptor_DARTS_BEST_psat40` + `sculptor_apogeeDR17_xmatch_notDART`  or Sculptor_DARTS_BEST_psat40_notAPOGEE + `sculptor_apogeeDR17_xmatch`.

Then, you should also include the stars observed with GMOS in the paper, they are in the table `targets_stellar_met.csv`. Since GMOS has relatively large systematics uncertainties on the RV, you should sum in quadrature to the RV_err the value 13.3 km/s (as in the paper).

Let me know if you have questions on the tables or things in Sculptor.
"""

# ╔═╡ 7db47590-b82b-4822-8a01-eaff96c9389f
md"""
## Dart
"""

# ╔═╡ 97277013-a8fc-4f70-b016-ef7f9025ce97
md"""
DART sample given from Federico. 

The sample is used in Battaglia + 2008 (Kinematic status of Scl). I cannot find the whole dataset publically available. 

Measurements taken using VLT/FLAMES. 

The paper reports $v_{\rm hel, sys} = 110.6\pm0.5$ and a velocity dispersion of $\sigma=10.1 \pm 0.3$
"""

# ╔═╡ 6964faa7-17d2-412c-b4a2-5f981f8c4b54
md"""
we reproduce the mean, but slightly lower velocity dispersion, maybe due to different selection criteria?
"""

# ╔═╡ 3dfa6b39-49a1-413c-a0a2-950b3bb22a0b
md"""
should be included in Tolstoy+23
"""

# ╔═╡ 23402e30-2064-494c-bd4a-8243cb474b61
begin 
	dart = CSV.read("$data_dir/Sculptor_DARTS_BEST_psat40.csv", DataFrame)
	rename!(dart, 
		"vel"=>"RV",
		"evel"=>"RV_err",
		"feh" => "Fe_H",		
		# "feh_lo" => "Fe_H_lo",		
		# "feh_hi" => "Fe_H_hi",		
	)

	dart = dart[:, [:ra, :dec, :source_id, :RV, :RV_err, ]]
end

# ╔═╡ 89729357-a2e7-4684-bde4-5b1684080f27
dart

# ╔═╡ 96eda577-555e-4ec4-bd92-3700f9b668e5
mean(dart.RV) ± sem(dart.RV)

# ╔═╡ 61a678df-7bb3-4daf-a9eb-1eae3ec1fbaf
std(dart.RV)

# ╔═╡ b540138f-eb4c-40be-ad0b-62bc0daebe56
dart.RV

# ╔═╡ baacb491-9dff-4ff4-b62c-5842d79794da
let
	fig, ax = FigAxis()
	hist!(dart.RV)
	fig
end

# ╔═╡ e731c5ea-ef35-4973-af22-09d8ca6e3624
begin
	fig, ax = FigAxis()

	p = scatter!(dart.ra, dart.dec, color=dart.RV,
		colormap=:bluesreds
	)
	
	Colorbar(fig[1, 2], p)
	fig
end

# ╔═╡ c470c2b9-093d-42ab-b96d-9dac231ccabc
md"""
## Apogee
"""

# ╔═╡ 45d6b4aa-ca44-4a71-afb8-ba6b2e674c7a
md"""
APOGEE DR 17 sample from federico's paper
"""

# ╔═╡ bb7f6769-ec92-460d-8423-449029175f79
begin 
	apogee_all = CSV.read("$data_dir/sculptor_apogeeDR17_xmatch.csv", DataFrame)

	rename!(apogee_all, 
		"Elliptical half-light radii"=>"r_h",
		"VHELIO_AVG"=>"RV",
		"VERR"=>"RV_err",
		"GAIAEDR3_PMRA"=>"pmra",
		"GAIAEDR3_PMDEC"=>"pmdec",
		"GAIAEDR3_SOURCE_ID" => "source_id",
		"RA (deg)"=>"ra",
		"Dec (deg)"=>"dec"
	)

	_apogee_all_filt = not.(ismissing.(apogee_all.RV))

	apogee_all = DataFrame(apogee_all[_apogee_all_filt, :])

	apogee_all[!, :RV] = float.(apogee_all.RV)

	apogee_all = apogee_all[:, [:source_id, :ra, :dec, :RV, :RV_err]]
end

# ╔═╡ fcc95615-f24a-40fa-a7d1-9e474df9f798
filt_apogee = sigma_clip(apogee_all.RV, 3)

# ╔═╡ 4d63430e-f59c-4c68-97e1-7eaa5679e55f
apogee = apogee_all[filt_apogee, :]

# ╔═╡ 1bc2541c-69ac-40aa-94f5-743707ca5a66
hist(apogee.RV)

# ╔═╡ 27063de7-01d4-48a8-a06f-cc24aec662d2
md"""
### Check apogee-dart xmatch against fed
"""

# ╔═╡ 48c62811-136f-4962-a42c-b1dd1fc74f8c
begin 
	apogee_notdart = CSV.read("$data_dir/sculptor_apogeeDR17_xmatch_notDART.csv", DataFrame)

	rename!(apogee_notdart, 
		"Elliptical half-light radii"=>"r_h",
		"VHELIO_AVG"=>"RV",
		"VERR"=>"RV_err",
		"GAIAEDR3_PMRA"=>"pmra",
		"GAIAEDR3_PMDEC"=>"pmdec",
		"RA (deg)"=>"ra",
		"Dec (deg)"=>"dec",
		"GAIAEDR3_SOURCE_ID" => "source_id"
	)

	_apogee_filt = not.(ismissing.(apogee_notdart.RV))

	apogee_notdart = DataFrame(apogee_notdart[_apogee_filt, :])

	nothing
end

# ╔═╡ 2e7ce524-573e-45c9-b0f4-ce9fea68e026
a2 = apogee_all[not.(xmatch(apogee_all, dart)[1]), :];

# ╔═╡ 0cf2d6b0-d1b9-42fb-8eca-0bc7bccfea70
lguys.angular_distance.([14.479505360741172], [-33.99473291943769
], [14.479505360741172], [-33.99473291943769])

# ╔═╡ 4dec130d-d325-4cc8-9be5-8771366878b5
14.479505360741172 == 14.479505360741172

# ╔═╡ 83cb3df5-813b-4fb6-b22a-08742105b919
apogee_all[5:5, :].ra == [14.479505360741172]

# ╔═╡ 854dd3cf-7c8d-4a9c-bb86-cd5bfb775c4b
dart[435:435, :].dec == apogee_all[5:5, :].dec

# ╔═╡ 5de5206d-c746-40fd-9089-5ca0e1b14dbc
println(dart[435, :dec])

# ╔═╡ f8775eb1-2cb9-4d81-8c02-39c43ceb9b45
sort(a2.source_id) == sort(apogee_notdart.source_id)

# ╔═╡ b6de4afb-9362-4618-a114-df460031e4f9
md"""
# GMOS alla Federico
"""

# ╔═╡ b7345279-4f80-47ad-a726-537571849eae
let
	global gmos = CSV.read("$data_dir/targets_stellar_met.csv", DataFrame)
	rename!(gmos,
		"RA"=>"ra",
		"Dec"=>"dec"
	)

	gmos_rv_sys_err = 13.3
	gmos.RV_err .= @. sqrt(gmos.RV_err^2 + gmos_rv_sys_err^2)

	filt, idx = xmatch(gmos, j24)
	println(filt)
	gmos[!, :source_id] = j24.source_id[idx]

	gmos = gmos
end

# ╔═╡ c31f9e07-c520-443b-94dc-787519021d01
md"""
## Tolstoy et al. 2023
"""

# ╔═╡ 15f2a8e2-90df-48a9-a7bf-e86955f566ce

begin 
	tolstoy23_all = lguys.load_fits("$data_dir/tolstoy_2023.fits")

	rename!(tolstoy23_all, 
		:Vlos => :RV,
		:e_Vlos => :RV_err,
		:RAJ2000 => :ra,
		:DEJ2000 => :dec,
		:GaiaDR3 => :source_id
	)


	tolstoy23_all = tolstoy23_all
end

# ╔═╡ 1864bc12-cbc2-4fd4-84a8-6bb300f17c1b
filt_tolstoy = sigma_clip(tolstoy23_all.RV, 3)

# ╔═╡ b36c2a37-2359-4a1a-98fc-3a1cd17fd790
tolstoy23 = tolstoy23_all[filt_tolstoy, :]

# ╔═╡ 1d01f9a2-2b77-4ad0-9893-b31562de7924
filt_tolstoy .&= tolstoy23_all.Mem .== "m"

# ╔═╡ 180fac98-678c-4a14-966c-385387c60ac3
md"""
## Walker et al. (2009)

Data downloaded from CDS associated with observation paper.
Should have 1365 scl members.
Weighted mean of measurements for repeated stars.

Note that we combine two dataframes, one with individual measurements and a second using weighted sums. However, only the second one contains the stellar coordinates, so we crossmatch them.
"""

# ╔═╡ 77830e77-50a4-48b7-b070-8fbd7508c173
let 
	global walker09_single
	walker09_single = lguys.load_fits("$data_dir/walker_2009.fits") 
	
	walker09_single[!, "Galaxy"] = [s[1:3] for s in walker09_single.Target];

end

# ╔═╡ dee71790-ffeb-477d-adbe-112731dfe134
let 
	global walker09_averaged = lguys.load_fits("$data_dir/walker+09.fits")

	rename!(walker09_averaged, 
		"RAJ2000"=>"ra",
		"DEJ2000" => "dec",
		)
	walker09_averaged[!, "Galaxy"] = [s[1:3] for s in walker09_averaged.Target]

end

# ╔═╡ 56fd0ddd-e28b-4cb8-8928-111664a6ec92
sum(walker09_single.Galaxy .== "Scl")

# ╔═╡ 9842ab04-903c-4750-a1c2-20383286ef6d
sum(walker09_averaged.Galaxy .== "Scl")

# ╔═╡ 89206511-c278-45f4-b5ad-11f63811200b
walker09_averaged.__HV_

# ╔═╡ 398ac7fd-8136-48fd-8f0f-dc5b8d2b004e
walker09_averaged.Target

# ╔═╡ 08755a1d-27d2-4144-9060-58cdd56a25ed
walker09_single.Target

# ╔═╡ 91ca7d42-4537-4777-b6d0-ebfd964e0171
unique(walker09_single.Target[walker09_single.Galaxy .== "Scl"])

# ╔═╡ 6c4b17df-74ba-46ed-814b-eb286372e824
begin
	walker09_all = leftjoin(walker09_single, walker09_averaged, on="Target", makeunique=true)


	_walker_filt = walker09_all.Galaxy .== "Scl"
	_walker_filt .&= not.(isnan.(walker09_all.Mmb))

	walker09_all = walker09_all[_walker_filt, :]

	filt, idxs = xmatch(walker09_all, j24)
	walker09_all[!, "source_id"] = j24.source_id[idxs]
	walker09_all[not.(filt), "source_id"].= -1

	println("failed to match ", sum(not.(filt)))
	
	rename!(walker09_all, "__HV_" => "RV", "e__HV_" => "RV_err")

	walker09_all[isnan.(walker09_all.RV_err), :RV_err] .= walker09_all.e_HV[isnan.(walker09_all.RV_err)]

	walker09_all
end

# ╔═╡ 5dd59d8b-d3f1-448a-a63c-8dca9e27c18e
sum(walker09_all.Mmb) # about 1365 :)

# ╔═╡ 86ffdea6-5f02-446c-9dc8-b6c18aa835fb
length(unique(walker09_all.Target)) == size(walker09_all, 1)

# ╔═╡ 3026841e-3479-4dcb-ae0d-462166100c2b
hist(walker09_all.Mmb)

# ╔═╡ 5d5fbb43-ea1d-4b1c-a567-7143be1a9a5a
filt_walker = walker09_all.Mmb .> 0.5

# ╔═╡ e1c2e02e-05de-4faf-af2f-f93759ddadfe
filt_walker_2 = sigma_clip(walker09_all[filt_walker, :].RV, 3)

# ╔═╡ a5a95eba-d282-4881-a84e-a25d4c83f114
walker09 = walker09_all[filt_walker, :]

# ╔═╡ 8c59d607-d484-497c-9ee7-751a4edd4992
hist(walker09.RV)

# ╔═╡ fe6fa910-f41e-4657-836b-7eda2f0cddb2
function add_xmatch!(df, new, suffix)
	leftjoin!(df, rename(n->"$(n)_$suffix", new), on="source_id"=>"source_id_$suffix")
end

# ╔═╡ e472cbb6-258e-4303-85e2-56f26358c97b
let
	global all_stars 

	all_stars = copy(j24)
	# add_xmatch!(all_stars, dart, "dart")
	add_xmatch!(all_stars, apogee_all, "apogee")
	add_xmatch!(all_stars, walker09_all, "w09")
	add_xmatch!(all_stars, tolstoy23_all, "t23")
	add_xmatch!(all_stars, gmos, "gmos")


	rename!(all_stars,
		:dr2_radial_velocity => "RV_gaia",
		:dr2_radial_velocity_error => "RV_err_gaia",
	)

end

# ╔═╡ 0d2dbc73-1ded-46e3-b142-9bc7b777728d
all_stars.RV_gmos[not.(ismissing.(all_stars.RV_gmos))]

# ╔═╡ bd51fe42-7e39-4cd8-8065-58ab9814f966
sum(not.(ismissing.(all_stars.RV_apogee))) == length(apogee_all.RV)

# ╔═╡ ab49efb3-2ab7-47d0-a3a5-a342c789aa9b
sum(not.(ismissing.(all_stars.RV_w09))) , length(walker09_all.RV)

# ╔═╡ de762a39-b430-4452-ba87-8b8cf1ad9852
sum(not.(ismissing.(all_stars.RV_t23))) == length(tolstoy23_all.RV)

# ╔═╡ 89552b84-d12e-4d88-a58e-8b89ad4b2569
md"""
# Validation for xmatch
"""

# ╔═╡ e6f2de3b-ce32-4d61-851f-4e42fcce95c0
function plot_xmatch_radec(suffix)

	ra2 = all_stars[:, "ra_$suffix"]
	filt = not.(ismissing.(ra2))
	ra2 = ra2[filt]

	ra1 = all_stars.ra[filt]
	dec1 = all_stars.dec[filt]
	dec2 = all_stars[filt, "dec_$suffix"]

	fig, ax = FigAxis(
		xlabel="ra / degrees",
		ylabel="dec / degrees",
		aspect = secd(-33.7)
	)

	scatter!(ra1, dec1)
	scatter!(float.(ra2), float.(dec2), markersize=5)

	fig
end

# ╔═╡ 5b2fceff-9c3e-472d-9310-31e920137e41
plot_xmatch_radec("apogee")

# ╔═╡ 4b305b83-1a3b-48a6-b19f-6f3ebed0768f
plot_xmatch_radec("t23")

# ╔═╡ c218cfa8-2f55-4957-bcdd-8b3970fe639a
plot_xmatch_radec("w09")

# ╔═╡ f7213298-3dcb-49ac-a9f0-a129a03423aa
plot_xmatch_radec("gmos")

# ╔═╡ f7ec8bba-9f45-435b-b67c-33182e992dfd
md"""
# Cross study RV
"""

# ╔═╡ 93d185f2-1eaa-4d35-87dd-b84f385483de
function filt_missing(col, verbose=false; low=-Inf, high=Inf)
	filt =  not.(ismissing.(col)) .& not.(isnan.(col))
	filt1 = high .> col .> low
	if verbose
		println("excluding ", sum(not.(filt1)[filt]), " outliers")
	end
	return filt .& filt1
end

# ╔═╡ 7ebe7832-b9f0-4d67-81e4-27c45b3aa760
let
	filt, idx = xmatch(dart, tolstoy23)
	println("matched ", sum(filt), "out of", length(filt))

	df_matched = copy(dart[filt, :])
	df_matched[:, :RV_tolstoy] = tolstoy23[idx[filt], :RV]
	df_matched[:, :RV_err_tolstoy] = tolstoy23[idx[filt], :RV_err]
	
	rv1  = df_matched.RV
	rv1_err  = df_matched.RV_err
	
	rv2  = df_matched.RV_tolstoy
	rv2_err  = df_matched.RV_err_tolstoy

	filt = filt_missing(rv1, true; low=75, high=150) .& filt_missing(rv2, true;  low=75, high=150)


	fig, ax = FigAxis(
		xlabel = "RV DART",
		ylabel = "RV Tolstoy",
		aspect=DataAspect()
	)

	errscatter!(rv1[filt], rv2[filt], xerr=rv1_err[filt], yerr=rv2_err[filt])

	lines!([80, 150], [80, 150], color=:black)
	
	fig
end

# ╔═╡ 36d2d86f-6a75-46f4-b48f-36137e95e90d
filt_missing(all_stars.RV_gaia)

# ╔═╡ 9213cc36-74d1-452f-bd9a-eb5c5cab1f87
function compare_rv(study1, study2)
	rv1  = all_stars[:, "RV_$study1"]
	rv1_err  = all_stars[:, "RV_err_$study1"]
	
	rv2  = all_stars[:, "RV_$study2"]
	rv2_err  = all_stars[:, "RV_err_$study2"]

	filt = filt_missing(rv1, true; low=75, high=150) .& filt_missing(rv2, true;  low=75, high=150)

	println("matched ", sum(filt), " stars")


	if sum(filt) == 0
		println("nothing to plot")
		return
	end
	
	println("plotting ", sum(filt), " stars")


	fig, ax = FigAxis(
		xlabel = "RV $study1",
		ylabel = "RV $study2",
		aspect=DataAspect()
	)

	errscatter!(rv1[filt], rv2[filt], xerr=rv1_err[filt], yerr=rv2_err[filt])

	lines!([90, 150], [90, 150], color=:black)
	return fig
end

# ╔═╡ f9e1bcbd-2def-4e64-b69e-d1cfcf7ffedd
compare_rv("t23", "w09")

# ╔═╡ f61debe4-8b23-4415-b77c-43c4465ccfcb
compare_rv("t23", "apogee")

# ╔═╡ 102c73ef-2c95-4784-80df-ed0312511c00
compare_rv("apogee", "w09")

# ╔═╡ d3333b48-aa4e-42c1-9e0a-bbff98e3647d
all_studies = ["gmos", "apogee", "w09", "t23"]

# ╔═╡ 8bc140fa-6f4e-4ec5-bc95-989fc0ea48d1
function safe_weighted_mean(values, errors)
	filt = filt_missing(values)
	filt .&= filt_missing(errors)
	if sum(filt) == 0
		return missing
	end

	x = float.(values[filt]) .± float.(errors[filt])

	return weightedmean(x), std(values[filt])
end

# ╔═╡ 88ed48b4-baa9-4ac0-86e1-8348edcd59b4
begin 
	rvs = [all_stars[:, "RV_$study"] for study in all_studies]
	rv_errs = [all_stars[:, "RV_err_$study"] for study in all_studies]

	rvs = hcat(rvs...)
	rv_errs = hcat(rv_errs...)
end

# ╔═╡ 222bb254-8b65-44d3-b3d2-b08fbcbe9950
all_studies

# ╔═╡ e3f05ee2-cc5f-437e-801d-3c7d842af709
all_stars[:, "RV_w09"]

# ╔═╡ a13a33f1-65c4-4217-8636-639b1e14d109
sum(filt_missing(all_stars[:, "RV_err_w09"]))

# ╔═╡ ee3c22db-6b6b-4c30-8d1f-86b103c018fc
sum(filt_missing(walker09_all.RV_err))

# ╔═╡ 6bc02c4f-31a2-4e4e-8612-e66f8cc9c93e
rvs

# ╔═╡ 11fcf4f8-fcd5-4426-a4e9-b046138bde1b
rv_errs

# ╔═╡ 877706f1-08a0-496a-890d-622e3a2fd9ec
RV_a_std = [safe_weighted_mean(rvs[i, :], rv_errs[i, :]) for i in 1:size(rvs, 1)]

# ╔═╡ 73bd1553-e2ae-4bfb-aac1-0880346f5054
begin
	filt_is_meas = not.(ismissing.(RV_a_std))
	rv_meas = all_stars[filt_is_meas, :]
	RV = Measurements.value.(first.(RV_a_std[filt_is_meas]))
	RV_err = Measurements.uncertainty.(first.(RV_a_std[filt_is_meas]))
	RV_std = (last.(RV_a_std[filt_is_meas]))

	rv_meas[!, :RV] = RV
	rv_meas[!, :RV_err] = RV_err
	rv_meas[!, :RV_std] = RV_std

end

# ╔═╡ 01b3135f-7a72-4669-b586-4bc5894464ad
sum(filt_is_meas)

# ╔═╡ d11edca7-b9ae-4269-9e1b-661d59bd965e
all_stars[not.(ismissing.(all_stars.RV_gmos)), :].source_id

# ╔═╡ c724e720-18ca-4905-a4b6-39fc47abe39d
j24[j24.source_id .== 5006419626331394048, :]

# ╔═╡ 733fe42e-b7a5-4285-8c73-9a41e4488d40
begin 
	memb_filt = rv_meas.PSAT .> 0.2
	memb_filt .&= 150 .> rv_meas.RV .> 60

	quality_score = rv_meas.RV_std ./ rv_meas.RV_err
	quality_score[isnan.(rv_meas.RV_std)] .= 0

	memb_filt .&= quality_score .< 5
end

# ╔═╡ 74b10a3e-1342-454f-8eed-77b371f81edf
lines(Arya.histogram(quality_score, normalization=:none))

# ╔═╡ d8800a31-1ed3-422f-ac51-90f18cf61c29
sum(quality_score .< Inf)

# ╔═╡ 015f2345-74d4-4296-b0dc-35d6e13ad8cc
sum(quality_score .< 5)

# ╔═╡ cc6c65db-ef57-4745-8ada-e11427274a77
memb_stars = rv_meas[memb_filt, :]

# ╔═╡ 7f4a5254-ed6f-4faa-a71e-4b4986a99d45
hist(memb_stars.RV)

# ╔═╡ a1938588-ca40-4844-ab82-88c4254c435b
length(memb_stars.RV)

# ╔═╡ 3655b6a0-ac9a-4a49-86fd-6a6484409819
md"""
## Check RV means make sense
"""

# ╔═╡ 7f13339e-6a0c-4944-8a36-5ab136fd8415
function compare_rv_mean(study1, rv_meas=rv_meas)
	rv1  = rv_meas[:, "RV_$study1"]
	rv1_err  = rv_meas[:, "RV_err_$study1"]
	
	rv2  = rv_meas[:, "RV"]
	rv2_err  = rv_meas[:, "RV_err"]

	filt = filt_missing(rv1, true; low=75, high=150) .& filt_missing(rv2, true;  low=75, high=150)

	println("matched ", sum(filt), " stars")


	if sum(filt) == 0
		println("nothing to plot")
		return
	end
	
	println("plotting ", sum(filt), " stars")


	fig, ax = FigAxis(
		xlabel = "RV mean",
		ylabel = "RV $study1 - RV mean",
		aspect=DataAspect()
	)

	x = float.(rv2[filt])
	y = float.(rv1[filt] .- rv2[filt])
	xerr = float.(rv2_err[filt])
	yerr = float.(rv1_err[filt])
	errscatter!(x, y, xerr=xerr, yerr=yerr, alpha=0.1)
	
	hlines!(0, color=:black)
	return fig
end

# ╔═╡ aaaf5ba2-c9ed-41ec-a22a-d78ed96fd84e
compare_rv_mean("w09", memb_stars)

# ╔═╡ a2370516-e7ec-4502-ae4b-b111bcf68d36
compare_rv_mean("apogee", memb_stars)

# ╔═╡ 78e7aff0-3658-4d64-b117-58189a85307a
compare_rv_mean("t23", memb_stars)

# ╔═╡ f4f9dd06-1a1a-458b-be75-05d52623580c
compare_rv_mean("gmos", )

# ╔═╡ 74152829-27ef-4d8d-8b32-ed30a18f30e4
rv_meas.RV_gmos[not.(ismissing.(rv_meas.RV_gmos))]

# ╔═╡ abbd2a53-e077-4af7-a168-b571e1a906b8
xlabel = L"radial velocity / km s$^{-1}$"

# ╔═╡ 5f71b8e9-5540-4431-8028-4ce14c8d7856
let 
	fig, ax = FigAxis(xlabel=xlabel)
	bins = Arya.make_bins(apogee_all.RV, Arya.calc_limits(apogee_all.RV), bandwidth=1)

	h = Arya.histogram(apogee.RV, bins, normalization=:none)
	scatter!(h)
	
	h = Arya.histogram(apogee_all.RV, bins, normalization=:none)

	scatter!(h, markersize=5)


	fig
end

# ╔═╡ 9d1495e8-5f8a-4892-b5b5-b22f3eb6ab7c
let 
	fig, ax = FigAxis(xlabel=xlabel)
	bins = Arya.make_bins(walker09_all.RV, Arya.calc_limits(walker09_all.RV), bandwidth=3)


	
	h = Arya.histogram(walker09_all.RV, bins, normalization=:none)

	lines!(h, label="all")
	
	h = Arya.histogram(walker09.RV, bins, normalization=:none)
	lines!(h, label="selected")

	axislegend()

	fig
end

# ╔═╡ 95ef12d2-174a-47c3-a106-9b4005b2a80d
md"""
# Full fits to data
"""

# ╔═╡ 9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
function fit_rv_sigma(rv, rv_err; N=3_000, p=0.16, burn=0.2)
	μ = median(rv)
	μ_p = NaN
	μ_err = NaN

	σ = std(rv)
	σ_p = NaN
	σ_err = NaN

	return μ, σ, μ_err, σ_err
end

# ╔═╡ 9380d9d1-58b2-432d-8528-d247cf5724e9
function plot_sample_normal_fit(sample, props; kwargs...)
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density";
		kwargs...
	)
	h = Arya.histogram(Float64.(sample.RV), normalization=:pdf)
	
	errscatter!(midpoints(h.bins), h.values, yerr=h.err, color=:black)

	μ, σ, _, _ = props
	x_model = LinRange(80, 140, 1000)
	y_model = lguys.gaussian.(x_model, μ, σ)
	lines!(x_model, y_model)

	fig
end

# ╔═╡ 764b5306-20f9-4810-8188-1bdf9482260f
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = Arya.histogram(Float64.(memb_stars.RV), 30, normalization=:pdf)
	
	errscatter!(midpoints(h.bins), h.values, yerr=h.err, color=COLORS[6])

	fig
end

# ╔═╡ 61e15c47-c454-48db-95f9-02abe052676e
mean(memb_stars.RV)

# ╔═╡ d5938fc3-9c8a-4e33-8401-500b4201df14
std(memb_stars.RV) / sqrt(length(memb_stars.RV))

# ╔═╡ 6514815e-e1f3-44a1-8c0b-a90f13e0077e
std(memb_stars.RV)

# ╔═╡ 5688e9c9-515c-4973-820f-1215031253f2
import StatsBase: var

# ╔═╡ be10c541-246f-4fcf-a5e5-416ecb69d7e7
function sesd(x)
	k = kurtosis(x)
	n = length(x)
	s = var(x)
	sev = sqrt(1/n * (k .- (n-3)/(n-1)) .* s .^2)

	return 1 / 2sqrt(s) .* sev
end

# ╔═╡ 30cf37dc-a042-4ed1-88d4-21dd15539236
sesd(memb_stars.RV)

# ╔═╡ 5bee3041-fe7f-4f67-a73d-fb60ed006959
md"""
## For each dataset
"""

# ╔═╡ ebdf8cfc-ebff-4331-bab2-998a734f2cd1
datasets = Dict(
	:apogee => apogee,
	:dart => dart,
	:walker => walker09,
	:tolstoy => tolstoy23
)

# ╔═╡ ec4bc80f-22e6-49f9-a589-5c5bc9e50a8b
props = Dict(name => fit_rv_sigma(float.(data.RV), float.(data.RV_err)) for (name, data) in datasets)

# ╔═╡ eaf47845-95dd-4f7d-bf76-3d9b1154711a
for key in keys(props)
	println(key)
	println(props[key])
	@info plot_sample_normal_fit(datasets[key], props[key], title=string(key))
end

# ╔═╡ b5faa873-d5c9-447b-b90e-f693db26e6c2
md"""
### Comparing sample properties
"""

# ╔═╡ 46434fa6-09df-4b02-9270-cbdcc9648e38
let
	ks = keys(props)
	y = [props[k][1] for k in ks]
	yerr = [props[k][3] for k in ks]

	x = collect(1:length(ks))
	
	fig, ax = FigAxis(
		xticks=(x, string.(ks)),
		xminorticksvisible=false,
		ylabel=xlabel
	)


	errscatter!(x, y, yerr=yerr)
	fig
end

# ╔═╡ dc4c0453-fed0-4733-9f59-0b2df506b45c
let
	ks = keys(props)
	y = [props[k][2] for k in ks]
	yerr = [props[k][4] for k in ks]

	x = collect(1:length(ks))
	
	fig, ax = FigAxis(
		xticks=(x, string.(ks)),
		xminorticksvisible=false,
		ylabel=L"$\sigma_{v}$ / km s$^{-1}$"
	)


	errscatter!(x, y, yerr=yerr)
	fig
end

# ╔═╡ 7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
md"""
## Binned properties
"""

# ╔═╡ 30e3dc5b-3ce6-4dd7-9c2a-c82774909a8c
sum(not.(ismissing.(memb_stars.RV_gmos)))

# ╔═╡ 33f54afc-cdb9-4eb8-887f-5a43281b837c
let
	fig = Figure(figsize=(900, 400))
	ax = Axis(fig[1,1],
		xlabel = "R / arcmin",
		ylabel = L"RV / km s$^{-1}$",
		limits=(nothing, (60, 150))
	)

	#scatter!(memb_stars.r_ell, memb_stars.RV_dart, label="DART", markersize=5)
	scatter!(memb_stars.r_ell, memb_stars.RV_apogee, label="APOGEE", markersize=5)
	scatter!(rv_meas.r_ell, rv_meas.RV_gmos, label="GMOS", markersize=10)
	scatter!(memb_stars.r_ell, memb_stars.RV_w09, label="Walker+09", markersize=5)
	scatter!(memb_stars.r_ell, memb_stars.RV_t23, label="tolstoy + 23", markersize=3)

	Legend(fig[1,2], ax)


	resize_to_layout!(fig)
	fig
end

# ╔═╡ 725b9c7d-93f9-48c3-b97b-8a1147a16f78
begin 
	out_df = copy(rv_meas)
	for col in names(out_df)
		val = out_df[:, col]
		if eltype(val) <: Union{Missing, <:AbstractFloat}
			out_df[:, col] = replace(out_df[:, col], missing=>NaN)
		elseif eltype(val) <: Union{Missing, String}
			out_df[:, col] = replace(out_df[:, col], missing=>"")
		elseif eltype(val) <: Union{Missing, Integer}
			out_df[:, col] = replace(out_df[:, col], missing=>0)
		end	end

	select!(out_df, Not([:o_Target_w09, :Nspec_t23, ]))
	disallowmissing!(out_df)
end

# ╔═╡ 62301344-869a-4ec0-8299-29f0ff5d1c15
sum(sum(eachcol(ismissing.(out_df))))

# ╔═╡ 70290654-394f-4679-aaab-38345874e2e3
sum(sum(eachcol(ismissing.(rv_meas))))

# ╔═╡ f1912692-d890-4f97-ba98-7b226f29e9c8
sum(sum(eachcol(ismissing.(memb_stars))))

# ╔═╡ 4fdccf0a-adfa-4f05-bb48-df8a1efba940
out_df.:q__Fe_H__t23

# ╔═╡ 1a0eac4e-4100-4da6-820b-19c34e283118
lguys.write_fits(joinpath(data_dir, "sculptor_all_rv.fits"), out_df)

# ╔═╡ a5a5a6a9-1b02-4f30-b6e3-426d1c996eef
# if ispath(joinpath(data_dir, "sculptor_all_rv.fits"))
# 	rm(joinpath(data_dir, "sculptor_all_rv.fits"))
# end

# ╔═╡ Cell order:
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═202e0f8b-b417-4597-a737-7c60c0575fd3
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═da5a3b57-72bc-46e1-b1a0-6c02eb101626
# ╠═89729357-a2e7-4684-bde4-5b1684080f27
# ╠═3d8b3c55-153b-4a4a-9289-78f34df07abc
# ╠═d4123ffd-cb32-486e-a93d-f48d7112a831
# ╠═248c2d5f-85cc-44be-bc63-43d4af470182
# ╠═dfd2fbee-9993-4553-b693-6fb71a7b11a2
# ╠═7a50a176-96a5-4098-88d6-0fa2874d0f90
# ╟─6f2359d8-332b-11ef-0db9-f1f06474c561
# ╟─7db47590-b82b-4822-8a01-eaff96c9389f
# ╟─97277013-a8fc-4f70-b016-ef7f9025ce97
# ╟─6964faa7-17d2-412c-b4a2-5f981f8c4b54
# ╟─3dfa6b39-49a1-413c-a0a2-950b3bb22a0b
# ╠═7ebe7832-b9f0-4d67-81e4-27c45b3aa760
# ╠═96eda577-555e-4ec4-bd92-3700f9b668e5
# ╠═61a678df-7bb3-4daf-a9eb-1eae3ec1fbaf
# ╠═23402e30-2064-494c-bd4a-8243cb474b61
# ╠═b540138f-eb4c-40be-ad0b-62bc0daebe56
# ╠═baacb491-9dff-4ff4-b62c-5842d79794da
# ╠═e731c5ea-ef35-4973-af22-09d8ca6e3624
# ╟─c470c2b9-093d-42ab-b96d-9dac231ccabc
# ╟─45d6b4aa-ca44-4a71-afb8-ba6b2e674c7a
# ╠═bb7f6769-ec92-460d-8423-449029175f79
# ╠═5f71b8e9-5540-4431-8028-4ce14c8d7856
# ╠═fcc95615-f24a-40fa-a7d1-9e474df9f798
# ╠═4d63430e-f59c-4c68-97e1-7eaa5679e55f
# ╠═1bc2541c-69ac-40aa-94f5-743707ca5a66
# ╟─27063de7-01d4-48a8-a06f-cc24aec662d2
# ╠═48c62811-136f-4962-a42c-b1dd1fc74f8c
# ╠═2e7ce524-573e-45c9-b0f4-ce9fea68e026
# ╠═0cf2d6b0-d1b9-42fb-8eca-0bc7bccfea70
# ╠═4dec130d-d325-4cc8-9be5-8771366878b5
# ╠═83cb3df5-813b-4fb6-b22a-08742105b919
# ╠═854dd3cf-7c8d-4a9c-bb86-cd5bfb775c4b
# ╠═5de5206d-c746-40fd-9089-5ca0e1b14dbc
# ╠═f8775eb1-2cb9-4d81-8c02-39c43ceb9b45
# ╟─b6de4afb-9362-4618-a114-df460031e4f9
# ╠═b7345279-4f80-47ad-a726-537571849eae
# ╟─c31f9e07-c520-443b-94dc-787519021d01
# ╠═15f2a8e2-90df-48a9-a7bf-e86955f566ce
# ╠═b36c2a37-2359-4a1a-98fc-3a1cd17fd790
# ╠═1864bc12-cbc2-4fd4-84a8-6bb300f17c1b
# ╠═1d01f9a2-2b77-4ad0-9893-b31562de7924
# ╠═180fac98-678c-4a14-966c-385387c60ac3
# ╠═5dd59d8b-d3f1-448a-a63c-8dca9e27c18e
# ╠═77830e77-50a4-48b7-b070-8fbd7508c173
# ╠═dee71790-ffeb-477d-adbe-112731dfe134
# ╠═56fd0ddd-e28b-4cb8-8928-111664a6ec92
# ╠═9842ab04-903c-4750-a1c2-20383286ef6d
# ╠═89206511-c278-45f4-b5ad-11f63811200b
# ╠═398ac7fd-8136-48fd-8f0f-dc5b8d2b004e
# ╠═08755a1d-27d2-4144-9060-58cdd56a25ed
# ╠═91ca7d42-4537-4777-b6d0-ebfd964e0171
# ╠═6c4b17df-74ba-46ed-814b-eb286372e824
# ╠═86ffdea6-5f02-446c-9dc8-b6c18aa835fb
# ╠═3026841e-3479-4dcb-ae0d-462166100c2b
# ╠═9d1495e8-5f8a-4892-b5b5-b22f3eb6ab7c
# ╠═e1c2e02e-05de-4faf-af2f-f93759ddadfe
# ╠═5d5fbb43-ea1d-4b1c-a567-7143be1a9a5a
# ╠═a5a95eba-d282-4881-a84e-a25d4c83f114
# ╠═8c59d607-d484-497c-9ee7-751a4edd4992
# ╠═fe6fa910-f41e-4657-836b-7eda2f0cddb2
# ╠═e472cbb6-258e-4303-85e2-56f26358c97b
# ╠═0d2dbc73-1ded-46e3-b142-9bc7b777728d
# ╠═bd51fe42-7e39-4cd8-8065-58ab9814f966
# ╠═ab49efb3-2ab7-47d0-a3a5-a342c789aa9b
# ╠═de762a39-b430-4452-ba87-8b8cf1ad9852
# ╟─89552b84-d12e-4d88-a58e-8b89ad4b2569
# ╠═e6f2de3b-ce32-4d61-851f-4e42fcce95c0
# ╠═5b2fceff-9c3e-472d-9310-31e920137e41
# ╠═4b305b83-1a3b-48a6-b19f-6f3ebed0768f
# ╠═c218cfa8-2f55-4957-bcdd-8b3970fe639a
# ╠═f7213298-3dcb-49ac-a9f0-a129a03423aa
# ╠═f7ec8bba-9f45-435b-b67c-33182e992dfd
# ╠═93d185f2-1eaa-4d35-87dd-b84f385483de
# ╠═36d2d86f-6a75-46f4-b48f-36137e95e90d
# ╠═9213cc36-74d1-452f-bd9a-eb5c5cab1f87
# ╠═f9e1bcbd-2def-4e64-b69e-d1cfcf7ffedd
# ╠═f61debe4-8b23-4415-b77c-43c4465ccfcb
# ╠═102c73ef-2c95-4784-80df-ed0312511c00
# ╠═d3333b48-aa4e-42c1-9e0a-bbff98e3647d
# ╠═8bc140fa-6f4e-4ec5-bc95-989fc0ea48d1
# ╠═88ed48b4-baa9-4ac0-86e1-8348edcd59b4
# ╠═222bb254-8b65-44d3-b3d2-b08fbcbe9950
# ╠═e3f05ee2-cc5f-437e-801d-3c7d842af709
# ╠═a13a33f1-65c4-4217-8636-639b1e14d109
# ╠═ee3c22db-6b6b-4c30-8d1f-86b103c018fc
# ╠═6bc02c4f-31a2-4e4e-8612-e66f8cc9c93e
# ╠═11fcf4f8-fcd5-4426-a4e9-b046138bde1b
# ╠═877706f1-08a0-496a-890d-622e3a2fd9ec
# ╠═73bd1553-e2ae-4bfb-aac1-0880346f5054
# ╠═01b3135f-7a72-4669-b586-4bc5894464ad
# ╠═d11edca7-b9ae-4269-9e1b-661d59bd965e
# ╠═c724e720-18ca-4905-a4b6-39fc47abe39d
# ╠═733fe42e-b7a5-4285-8c73-9a41e4488d40
# ╠═74b10a3e-1342-454f-8eed-77b371f81edf
# ╠═d8800a31-1ed3-422f-ac51-90f18cf61c29
# ╠═015f2345-74d4-4296-b0dc-35d6e13ad8cc
# ╠═cc6c65db-ef57-4745-8ada-e11427274a77
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╠═a1938588-ca40-4844-ab82-88c4254c435b
# ╟─3655b6a0-ac9a-4a49-86fd-6a6484409819
# ╠═7f13339e-6a0c-4944-8a36-5ab136fd8415
# ╠═aaaf5ba2-c9ed-41ec-a22a-d78ed96fd84e
# ╠═a2370516-e7ec-4502-ae4b-b111bcf68d36
# ╠═78e7aff0-3658-4d64-b117-58189a85307a
# ╠═f4f9dd06-1a1a-458b-be75-05d52623580c
# ╠═74152829-27ef-4d8d-8b32-ed30a18f30e4
# ╠═abbd2a53-e077-4af7-a168-b571e1a906b8
# ╟─95ef12d2-174a-47c3-a106-9b4005b2a80d
# ╠═9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
# ╠═9380d9d1-58b2-432d-8528-d247cf5724e9
# ╠═764b5306-20f9-4810-8188-1bdf9482260f
# ╠═61e15c47-c454-48db-95f9-02abe052676e
# ╠═d5938fc3-9c8a-4e33-8401-500b4201df14
# ╠═6514815e-e1f3-44a1-8c0b-a90f13e0077e
# ╠═5688e9c9-515c-4973-820f-1215031253f2
# ╠═30cf37dc-a042-4ed1-88d4-21dd15539236
# ╠═be10c541-246f-4fcf-a5e5-416ecb69d7e7
# ╟─5bee3041-fe7f-4f67-a73d-fb60ed006959
# ╠═ebdf8cfc-ebff-4331-bab2-998a734f2cd1
# ╠═ec4bc80f-22e6-49f9-a589-5c5bc9e50a8b
# ╠═eaf47845-95dd-4f7d-bf76-3d9b1154711a
# ╟─b5faa873-d5c9-447b-b90e-f693db26e6c2
# ╠═46434fa6-09df-4b02-9270-cbdcc9648e38
# ╠═dc4c0453-fed0-4733-9f59-0b2df506b45c
# ╟─7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═30e3dc5b-3ce6-4dd7-9c2a-c82774909a8c
# ╠═33f54afc-cdb9-4eb8-887f-5a43281b837c
# ╠═725b9c7d-93f9-48c3-b97b-8a1147a16f78
# ╠═62301344-869a-4ec0-8299-29f0ff5d1c15
# ╠═70290654-394f-4679-aaab-38345874e2e3
# ╠═f1912692-d890-4f97-ba98-7b226f29e9c8
# ╠═4fdccf0a-adfa-4f05-bb48-df8a1efba940
# ╠═1a0eac4e-4100-4da6-820b-19c34e283118
# ╠═a5a5a6a9-1b02-4f30-b6e3-426d1c996eef
