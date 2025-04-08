### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	
	using CSV, DataFrames
	using PythonCall
	using Arya
	using CairoMakie

	import LilGuys as lguys
end

# ╔═╡ 7ed5bcf5-dfc3-4e79-a608-d503124a1e96
using LilGuys

# ╔═╡ d4123ffd-cb32-486e-a93d-f48d7112a831
include("../../utils/gaia_filters.jl")

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Radial Velocity CrossMatch

This notebook loads in radial velocity data from several sources (including Tolstoy+23, APOGEE, DART, and Fedrico++'s GMOS observations) and consolidates this data into a single combined radial velocity catalogue cross-matched with J+24's membership calculations.

With this catalogue, detailed radial velocity analysis can then be calculated.
"""

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, kurtosis, sem

# ╔═╡ 257caa92-ccce-44ff-88ee-6a9c550ae5d9
CairoMakie.activate!(type=:png)

# ╔═╡ dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
not = !

# ╔═╡ ef1bb6f5-00b8-405b-8702-244eca618644
import DensityEstimators: histogram, calc_limits, make_bins

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "data/"

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

# ╔═╡ 77e7884c-0360-4b7f-b9ad-e81be2516552
obs_properties = TOML.parsefile("observed_properties.toml")

# ╔═╡ 248c2d5f-85cc-44be-bc63-43d4af470182
params = GaiaFilterParams(obs_properties, filename="data/jensen+24_wide.fits")

# ╔═╡ 7a50a176-96a5-4098-88d6-0fa2874d0f90
j24 = read_gaia_stars(params)

# ╔═╡ dd1eee07-bb37-4de3-8781-47edd794c1f3
md"""
## Federico's samples
"""

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
APOGEE DR 17 sample from federico's paper is apogee_raw_2. I have also xmatched the field against apogee which I use for self-consistency.
"""

# ╔═╡ 0e14808d-df92-419b-b272-f3a08f4b86b1
apogee_raw_2 = CSV.read("$data_dir/sculptor_apogeeDR17_xmatch.csv", DataFrame)

# ╔═╡ d4d0a488-1c0e-4bc6-9a88-94b27d84e8ce
apogee_raw = LilGuys.read_fits("$data_dir/apogee_xmatch.fits")

# ╔═╡ bb7f6769-ec92-460d-8423-449029175f79
begin 
	apogee_all = copy(apogee_raw)
	rename!(apogee_all, 
		"FE_H" => "fe_h",
		"FE_H_ERR" => "fe_h_err"
	)

	_apogee_all_filt = not.(ismissing.(apogee_all.RV))

	apogee_all = DataFrame(apogee_all[_apogee_all_filt, :])

	apogee_all[!, :RV] = float.(apogee_all.RV)

	apogee_all = apogee_all[:, [:source_id, :ra, :dec, :RV, :RV_err, :RV_sigma, :RV_count, :fe_h, :fe_h_err]]
	apogee_all[ismissing.(apogee_all.fe_h), :fe_h] .= NaN
	apogee_all[ismissing.(apogee_all.fe_h_err), :fe_h_err] .= NaN

	apogee_all
end

# ╔═╡ 344f3c29-0873-4180-9025-51aeaeb2c681
CairoMakie.activate!(type="svg", pt_per_unit=1)

# ╔═╡ bf088da8-a7c5-46e7-8e5f-ef5b91c10fb8
filt_apogee = sigma_clip(apogee_all.RV, 3)

# ╔═╡ 4d63430e-f59c-4c68-97e1-7eaa5679e55f
apogee = apogee_all[filt_apogee, :]

# ╔═╡ 1bc2541c-69ac-40aa-94f5-743707ca5a66
hist(apogee.RV)

# ╔═╡ 63792b12-dde2-4361-88b2-1f5a789631ff
hist(filter(isfinite, apogee.fe_h))

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

# ╔═╡ f8775eb1-2cb9-4d81-8c02-39c43ceb9b45
sort(a2.source_id) == sort(apogee_notdart.source_id)

# ╔═╡ 6009cc0e-6936-4589-85aa-dd2864948b7c
let
	fig = Figure()
    ax = Axis(fig[1,1])

	scatter!(apogee.ra, apogee.dec, markersize=5)
	scatter!(apogee_raw_2[:, "RA (deg)"], apogee_raw_2[:, "Dec (deg)"], markersize=3)

	fig
end

# ╔═╡ b6de4afb-9362-4618-a114-df460031e4f9
md"""
## GMOS alla Federico
"""

# ╔═╡ 8e0095e7-2198-439d-bd19-976bdb1d3766
gmos_raw = CSV.read("$data_dir/sestito+23_gmos.csv", DataFrame)

# ╔═╡ b7345279-4f80-47ad-a726-537571849eae
begin
	gmos = copy(gmos_raw) 
	rename!(gmos,
		"RA"=>"ra",
		"Dec"=>"dec",
		"met_mean" => "fe_h",
		"met_err" => "fe_h_err"
	)

	let
		gmos_rv_sys_err = 13.3
		gmos.RV_err .= @. sqrt(gmos.RV_err^2 + gmos_rv_sys_err^2)
	
		filt, idx = xmatch(gmos, j24)
		println(filt)
		gmos[!, :source_id] = j24.source_id[idx]

	end

	gmos[!, :RV_count] .= 1.
	gmos[!, :RV_sigma] .= 0.
	
	gmos
end

# ╔═╡ 1151914d-b20f-4e46-8b25-e0d995e8b866
hist(filter(isfinite, gmos.fe_h))

# ╔═╡ c31f9e07-c520-443b-94dc-787519021d01
md"""
## Tolstoy et al. 2023
"""

# ╔═╡ 15f2a8e2-90df-48a9-a7bf-e86955f566ce

begin 
	tolstoy23_all = lguys.read_fits("$data_dir/tolstoy+23.fits")

	rename!(tolstoy23_all, 
		:Vlos => :RV,
		:e_Vlos => :RV_err,
		:s_Vlos => :RV_sigma,
		:Nspec => :RV_count,
		:RAJ2000 => :ra,
		:DEJ2000 => :dec,
		:GaiaDR3 => :source_id,
		:__Fe_H_ => :fe_h,
		:e__Fe_H_ => :fe_h_err,
	)

	tolstoy23_all.RV_err = @. max((1.448 - 0.021*tolstoy23_all.S_N)*tolstoy23_all.RV_err, 0.65)
	
	tolstoy23_all = tolstoy23_all
end

# ╔═╡ 624af663-758f-4b68-8c71-5e5e324ad3b1
md"""
Velocity correction from @arroyo-polonio 2024

``
	(∆v)′_{\rm mlm} = \max((1.448 - 0.021\,SNR) (∆v)_{\rm mlm} , 0.65{\rm km\,s}^{−1}) ,
``
"""

# ╔═╡ 1864bc12-cbc2-4fd4-84a8-6bb300f17c1b
filt_tolstoy = sigma_clip(tolstoy23_all.RV, 3)

# ╔═╡ b36c2a37-2359-4a1a-98fc-3a1cd17fd790
tolstoy23 = tolstoy23_all[filt_tolstoy, :]

# ╔═╡ 1d01f9a2-2b77-4ad0-9893-b31562de7924
filt_tolstoy .&= tolstoy23_all.Mem .== "m"

# ╔═╡ f064e833-5a6b-44be-9a5b-75599b39d150
hist(filter(isfinite, tolstoy23_all.fe_h))

# ╔═╡ fd85bfca-1ed2-4324-8933-5070c5c426a0
hist(filter(isfinite, tolstoy23_all.RV_count))

# ╔═╡ 31c61c48-86d7-4db9-a479-eb25c058c281
hist(filter(isfinite, tolstoy23_all.RV_err))

# ╔═╡ 9e488fb1-d88b-48b3-b0d5-8d306cfd9f83
hist(filter(x->x>0, tolstoy23_all.RV_sigma ./ tolstoy23_all.RV_err))

# ╔═╡ 180fac98-678c-4a14-966c-385387c60ac3
md"""
## Walker et al. (2009)

Data downloaded from CDS associated with observation paper.
Should have 1365 scl members.
Weighted mean of measurements for repeated stars.

Note that we combine two dataframes, one with individual measurements and a second using weighted sums. However, only the second one contains the stellar coordinates, so we crossmatch them.
"""

# ╔═╡ 77830e77-50a4-48b7-b070-8fbd7508c173
begin 
	walker09_single = lguys.read_fits("$data_dir/walker+09_observations.fits") 
	
	walker09_single[:, "Galaxy"] = [s[1:3] for s in walker09_single.Target];

	walker09_single
end

# ╔═╡ dee71790-ffeb-477d-adbe-112731dfe134
begin 
	walker09_averaged = lguys.read_fits("$data_dir/walker+09_summary.fits")

	rename!(walker09_averaged, 
		"RAJ2000"=>"ra",
		"DEJ2000" => "dec",
		)
	walker09_averaged[:, "Galaxy"] = [s[1:3] for s in walker09_averaged.Target]

	walker09_averaged
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
	
	rename!(walker09_all, "__HV_" => "RV", "e__HV_" => "RV_err",
		   )

	walker09_all[isnan.(walker09_all.RV_err), :RV_err] .= walker09_all.e_HV[isnan.(walker09_all.RV_err)]

	walker09_all[:, :RV_count] .= 0.
	walker09_all[:, :RV_sigma] .= 0.

	for (i, target) in enumerate(walker09_all.Target)
		RVs = walker09_single[walker09_single.Target .== target, :HV]
		walker09_all[i, :RV_count] = length(RVs)
		walker09_all[i, :RV_sigma] = std(RVs)
	end

	walker09_all
end

# ╔═╡ 5dd59d8b-d3f1-448a-a63c-8dca9e27c18e
sum(walker09_all.Mmb) # about 1365 :)

# ╔═╡ baa5285a-edb5-4e98-bbe4-d40a5fbe8ab9
apogee.RV_sigma

# ╔═╡ 1a6173db-6d5d-49c3-a2eb-9a77e591e747
sum(walker09_all.source_id .== -1)

# ╔═╡ 86ffdea6-5f02-446c-9dc8-b6c18aa835fb
length(unique(walker09_all.Target)) == size(walker09_all, 1)

# ╔═╡ 3026841e-3479-4dcb-ae0d-462166100c2b
hist(walker09_all.Mmb)

# ╔═╡ cf2c1c32-1fbc-4395-88cb-e1d9e410ef7d
hist(walker09_all.RV_count)

# ╔═╡ 10478b34-2af4-4833-9dd4-03192b4f1576
hist(filter(isfinite, walker09_all.RV_sigma))

# ╔═╡ a678426b-f895-46f0-92fd-be4c35e68698
hist(filter(isfinite, walker09_all.RV_err))

# ╔═╡ 6bdcb307-ee27-4a40-acbb-86787a4b5567
hist(filter(isfinite, walker09_all.RV_sigma ./ walker09_all.RV_err))

# ╔═╡ 5d5fbb43-ea1d-4b1c-a567-7143be1a9a5a
filt_walker = walker09_all.Mmb .>= 0.5

# ╔═╡ e1c2e02e-05de-4faf-af2f-f93759ddadfe
filt_walker_2 = sigma_clip(walker09_all[filt_walker, :].RV, 3)

# ╔═╡ a5a95eba-d282-4881-a84e-a25d4c83f114
walker09 = walker09_all[filt_walker, :]

# ╔═╡ 8c59d607-d484-497c-9ee7-751a4edd4992
hist(walker09.RV)

# ╔═╡ bbf49122-11b1-4272-a660-0437c6aa2b3f
md"""
# Combined stars
"""

# ╔═╡ fe6fa910-f41e-4657-836b-7eda2f0cddb2
function add_xmatch!(df, new, suffix)
	leftjoin!(df, rename(n->"$(n)_$suffix", new), on="source_id"=>"source_id_$suffix")
end

# ╔═╡ e472cbb6-258e-4303-85e2-56f26358c97b
all_stars = let

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
	all_stars

end

# ╔═╡ 0d2dbc73-1ded-46e3-b142-9bc7b777728d
all_stars.RV_gmos[not.(ismissing.(all_stars.RV_gmos))]

# ╔═╡ bd51fe42-7e39-4cd8-8065-58ab9814f966
@assert sum(not.(ismissing.(all_stars.RV_apogee))) == length(apogee_all.RV)

# ╔═╡ 8f243944-3d0b-48ce-9f90-8d448c089239
md"""
Some of the stars in Walker are xmatched to the same star, but error of order 5 stars
"""

# ╔═╡ ab49efb3-2ab7-47d0-a3a5-a342c789aa9b
@assert sum(not.(ismissing.(all_stars.RV_w09))) == sum(walker09_all.source_id .> 0)

# ╔═╡ de762a39-b430-4452-ba87-8b8cf1ad9852
@assert sum(not.(ismissing.(all_stars.RV_t23))) == length(tolstoy23_all.RV)

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

	scatter!(ra1, dec1, markersize=4)
	scatter!(float.(ra2), float.(dec2), markersize=2, marker=:circle, color=COLORS[2], alpha=1)

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

	w = @. 1 / (rv1_err^2 + rv2_err^2)
	μ = LilGuys.mean(rv1 .- rv2, w)
	δμ = sqrt(1/sum(w))
	@info "mean = $μ ± $δμ"
	errorscatter!(rv1[filt], rv2[filt], xerror=rv1_err[filt], yerror=rv2_err[filt])

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

	errorscatter!(rv1[filt], rv2[filt], xerror=rv1_err[filt], yerror=rv2_err[filt])

	w = @. 1/(rv1_err[filt]^2 + rv2_err[filt]^2)
	μ = LilGuys.mean(rv1[filt] .- rv2[filt], w)
	δμ = sqrt(1/sum(w))
	sigma = LilGuys.std(rv1[filt] .- rv2[filt] .- μ, w)
	@info "mean difference = $μ ± $δμ"
	@info "std = $sigma"
	@info "mean uncertainty = $(mean(1 ./ sqrt.(w)))"
	lines!([90, 150], [90, 150], color=:black)
	return fig
end

# ╔═╡ f9e1bcbd-2def-4e64-b69e-d1cfcf7ffedd
compare_rv("t23", "w09")

# ╔═╡ aa8ccf4d-fde4-4a73-a93c-8c11ff747da7
median(apogee.RV_count)

# ╔═╡ aa9b5b6f-f187-498e-b86d-f3b8a16a73ba
median(apogee.RV_sigma)

# ╔═╡ 30461ef6-1778-4612-b55d-8d16f0c6db12
mean(tolstoy23.RV_count), median(tolstoy23.RV_err), mean(tolstoy23.RV_sigma)

# ╔═╡ ca3c2080-af1b-4218-a29e-26eed51f0ef8
mean(walker09.RV_count),  median(walker09.RV_err), median(walker09.RV_sigma)

# ╔═╡ f61debe4-8b23-4415-b77c-43c4465ccfcb
compare_rv("t23", "apogee")

# ╔═╡ 102c73ef-2c95-4784-80df-ed0312511c00
compare_rv("apogee", "w09")

# ╔═╡ d3333b48-aa4e-42c1-9e0a-bbff98e3647d
all_studies = ["gmos", "apogee", "w09", "t23"]

# ╔═╡ 8bc140fa-6f4e-4ec5-bc95-989fc0ea48d1
function safe_weighted_mean(values, errors, sigmas, counts)
	filt = filt_missing(values)
	filt .&= filt_missing(errors)
	if sum(filt) == 0
		return missing
	end

	w = disallowmissing(1 ./ errors[filt] .^ 2)
	x = disallowmissing(values[filt])
	n = disallowmissing(counts[filt])
	s = disallowmissing(sigmas[filt])

	m = LilGuys.mean(x, w)
	ss_individual = @. s^2 * (n-1)
	ss_individual[isnan.(ss_individual)] .= 0
	ss_between = sum((x .- m) .^ 2 .* n)

	std_err = sqrt(1/sum(w))
	
	counts_tot = sum(n)
	sigma = sqrt((ss_between + sum(ss_individual) )/ (counts_tot - 1))
	return m, std_err, sigma, counts_tot 
end

# ╔═╡ d3ef6a98-8770-4359-9e37-94f08d2d9aa5
safe_weighted_mean([1.0, 2.0], [0.2, 0.1], [0., 0.,], [1., 1.])

# ╔═╡ 78d39e37-44b5-4171-ab45-0672f607b517
sqrt(1/(1/0.2^2 + 1/0.1^2))

# ╔═╡ 88ed48b4-baa9-4ac0-86e1-8348edcd59b4
begin 
	rvs = [all_stars[:, "RV_$study"] for study in all_studies]
	rv_errs = [all_stars[:, "RV_err_$study"] for study in all_studies]
	rv_sigmas = [all_stars[:, "RV_sigma_$study"] for study in all_studies]
	rv_counts = [all_stars[:, "RV_count_$study"] for study in all_studies]

	rvs = hcat(rvs...)
	rv_errs = hcat(rv_errs...)
	rv_counts = hcat(rv_counts...)
	rv_sigmas = hcat(rv_sigmas...)
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

# ╔═╡ 1798bdc1-1770-4dd1-a7ed-5b105509a31d


# ╔═╡ 877706f1-08a0-496a-890d-622e3a2fd9ec
RV_a_std = [safe_weighted_mean(rvs[i, :], rv_errs[i, :], rv_sigmas[i, :], rv_counts[i, :]) for i in 1:size(rvs, 1)]

# ╔═╡ 73bd1553-e2ae-4bfb-aac1-0880346f5054
begin
	filt_is_meas = not.(ismissing.(RV_a_std))
	rv_meas = copy(all_stars[filt_is_meas, :])
	RV = first.(RV_a_std[filt_is_meas])
	RV_err = [x[2] for x in RV_a_std[filt_is_meas]]
	RV_sigma = [x[3] for x in RV_a_std[filt_is_meas]]
	RV_count = [x[4] for x in RV_a_std[filt_is_meas]]

	rv_meas[!, :RV] = RV
	rv_meas[!, :RV_err] = RV_err
	rv_meas[!, :RV_sigma] = RV_sigma
	rv_meas[!, :RV_count] = RV_count

end

# ╔═╡ 8c1457d6-32d4-4727-b45e-1a46244708a6
RV_err

# ╔═╡ 63e28d21-c1ab-413e-ac7e-e6776bcf6ca5
rvs[filt_is_meas, :]

# ╔═╡ dbb51473-6506-4119-a11d-50df2ceb6dcf
rv_errs[filt_is_meas, :]

# ╔═╡ d913bf9b-c03b-4bce-a39d-7d23eca9b355
rv_sigmas[filt_is_meas, :]

# ╔═╡ 24b2455e-1c9a-4b9a-9d88-19b98d4b2165
RV_sigma

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

	quality_score = rv_meas.RV_sigma ./ rv_meas.RV_err
	quality_score[isnan.(rv_meas.RV_sigma)] .= 0

	memb_filt .&= quality_score .< 5
end

# ╔═╡ 74b10a3e-1342-454f-8eed-77b371f81edf
lines(histogram(filter(isfinite, quality_score), normalization=:none))

# ╔═╡ d8800a31-1ed3-422f-ac51-90f18cf61c29
sum(quality_score .< Inf)

# ╔═╡ 015f2345-74d4-4296-b0dc-35d6e13ad8cc
sum(quality_score .< 5)

# ╔═╡ 08e46c2e-d7c6-4bdc-9243-e18e6b9e9dc2
sum(quality_score .> 5)

# ╔═╡ cc6c65db-ef57-4745-8ada-e11427274a77
memb_stars = rv_meas[memb_filt, :]

# ╔═╡ 7f4a5254-ed6f-4faa-a71e-4b4986a99d45
hist(memb_stars.RV)

# ╔═╡ a1938588-ca40-4844-ab82-88c4254c435b
length(memb_stars.RV)

# ╔═╡ ab27964f-b070-4328-81c2-6ecf4b8cec1e
md"""
# Corrected coordinate frame
"""

# ╔═╡ 7dab615b-e0a1-431e-88c0-e1e9d9c29568
ra0, dec0 = obs_properties["ra"], obs_properties["dec"]

# ╔═╡ 57a77527-a14c-47a2-8949-1cbdafc49994
import Measurements: ±

# ╔═╡ 96eda577-555e-4ec4-bd92-3700f9b668e5
mean(dart.RV) ± sem(dart.RV)

# ╔═╡ 23d4b7c4-61f9-4748-b3c4-80eaad551914
distance = obs_properties["distance"] ± float.(obs_properties["distance_err"])

# ╔═╡ 5930b133-3e3c-43be-9ed6-f125ea4418d4
import Measurements

# ╔═╡ c4859522-5f50-4042-aa24-6bf9b214cb22
pmra = map(eachrow(rv_meas)) do row
	n_sigma = abs(row.pmra - obs_properties["pmra"]) / row.pmra_error
	if n_sigma < 3
		return row.pmra ± row.pmra_error
	else
		w1 = 1/row.pmra^2
		w2 = 1/obs_properties["pmra_err"]^2
		m = (w1*row.pmra + w2 * obs_properties["pmra"] ) / (w1+w2)
		err = 1/sqrt(w1 + w2)
		return m ± err
	end
end

# ╔═╡ ad7ae1a8-c5a4-4d79-bd86-35431fa71856
pmdec = map(eachrow(rv_meas)) do row
	n_sigma = abs(row.pmdec - obs_properties["pmdec"]) / row.pmdec_error
	if n_sigma < 3
		return row.pmdec ± row.pmdec_error
	else
		w1 = 1/row.pmdec^2
		w2 = 1/obs_properties["pmdec_err"]^2
		m = (w1*row.pmdec + w2 * obs_properties["pmdec"] ) / (w1+w2)
		err = 1/sqrt(w1 + w2)
		return m ± err
	end
end

# ╔═╡ d2c28cfd-652d-41fd-8ae2-4b5454761377
sind(1) *9

# ╔═╡ c0accd82-17a7-437c-9ba8-4db437071a5b
icrs = [
	lguys.ICRS(
		ra=rv_meas.ra[i], dec=rv_meas.dec[i], 
		distance=distance, pmra=pmra[i],
		pmdec=pmdec[i], 
		radial_velocity = rv_meas.RV[i] ± rv_meas.RV_err[i]
	) for i in 1:size(rv_meas, 1)]

# ╔═╡ 490b14a7-719d-4123-9173-e620288ff146
hist(filter(isfinite, rv_meas.pmra_error), bins=10 .^ LinRange(-2, 0, 30), axis=(;xscale=log10, xticks=[0.01, 0.03, 0.1, 0.3]))

# ╔═╡ 9ae17758-a548-48da-9bc2-a1b4f05e1818
hist(filter(isfinite, (rv_meas.pmra .- obs_properties["pmra"]) ./ (rv_meas.pmra_error .+ 0.017)), bins=100, axis=(; yscale=log10, yticks=Makie.automatic))

# ╔═╡ e06b8dc3-5e14-4c8a-80a2-e187de4fe0aa
hist(filter(isfinite, (rv_meas.pmdec .- obs_properties["pmdec"]) ./ (rv_meas.pmdec_error .+ 0.017)), bins=100, axis=(yscale=log10,yticks=Makie.automatic))

# ╔═╡ 383a9a9c-59da-4d3d-9cb9-253d9c685057
md"""
number of candidate members which significantly differ from systematic PM
"""

# ╔═╡ 7bf999a8-86b7-4e9f-b9c2-aa9de943ccad
sum((abs.(rv_meas.pmra .- obs_properties["pmra"]) ./ (rv_meas.pmra_error) .> 3)[rv_meas.PSAT .> 0.01])

# ╔═╡ 4180e4bd-944c-4817-a36d-35bc98a19bed
sum(rv_meas.pmra_error .< 0.02)

# ╔═╡ d0c9ad40-4bdc-4269-a1ab-3e501c84ce13
sum(rv_meas.pmdec_error .< 0.02)

# ╔═╡ 3f51257c-0cb6-4825-a41f-21fc6757060e
LilGuys.kms2pm(9, obs_properties["distance"])

# ╔═╡ d0842a73-dd1b-452f-97f9-50a23668fe14
methods(lguys.ICRS)

# ╔═╡ 11461563-dbcd-48fc-a97a-1a0391538462
gsr = lguys.transform.(lguys.GSR{Measurement{Float64}}, icrs)

# ╔═╡ 718ec8bc-cae6-4905-a336-b04964699b61
vra = [o.pmra for o in gsr] .* distance

# ╔═╡ 379bff0e-e65c-4abd-bce1-d6b776656bc8
vdec = [o.pmdec for o in gsr] .* distance

# ╔═╡ 6c1aad21-80cc-4dfb-80d5-4c1de4ecfe23
ϕ_pm = lguys.angular_distance.(rv_meas.ra, rv_meas.dec, ra0, dec0)

# ╔═╡ 65e5ecd2-f184-40cf-850d-326683560784
maximum(ϕ_pm)

# ╔═╡ 596bb8fc-b5f7-413e-bdef-8bd1968ea7cb
θ_pm = @. atand(rv_meas.xi, rv_meas.eta)

# ╔═╡ 9492f559-2a37-4d86-8247-2d217a387d8c
v_tan = @. sqrt(vra^2 + vdec^2)

# ╔═╡ 6f7a98dd-77cf-42a0-888d-8e34e2a0490f
rv = [o.radial_velocity for o in gsr]

# ╔═╡ cc68349e-da2f-4a73-92f5-6a021237accc
vz = @. rv * cosd(ϕ_pm) + sind(ϕ_pm) * (vdec * cosd(θ_pm) + vra * sind(θ_pm))

# ╔═╡ ae1f5ac1-265a-4f81-abe8-38d882b494bc
vx = @. rv * sind(ϕ_pm) * sind(θ_pm) +   vra * cosd(ϕ_pm)

# ╔═╡ f61c0c4c-446f-466f-8f64-95250a9177f1
vy = @. rv * sind(ϕ_pm) * cosd(θ_pm) +  vdec * cosd(ϕ_pm)

# ╔═╡ 64ea3168-59a4-4db2-8fec-f7facf4f0df6
vtot1 = @. sqrt(vx^2 + vy^2 + vz^2)

# ╔═╡ 0ef130e2-14b0-4fc6-9b0c-76ab18f3775d
vtot2 = @. sqrt(rv^2 + vra^2 + vdec^2)

# ╔═╡ 58e5855a-9279-4f1c-bc63-52d471e63389
vtot1 ./ vtot2

# ╔═╡ 13f3ac83-2fb3-4537-a3e6-bdb83d8ba72a
let
	fig = Figure()
	ax = Axis(fig[1, 1], 
		xlabel=L"$v_\textrm{los}$ / km\,s$^{-1}$",
		ylabel = L"$v_x - v_\textrm{los}$ / km\,s$^{-1}$",
		xgridvisible=false,
		ygridvisible=false,
	)

	errorscatter!(Measurements.value.(rv), Measurements.value.(vz) .- Measurements.value.(rv), 
		#xerr=Measurements.uncertainty.(rv), yerr=Measurements.uncertainty.(vz),
		alpha=0.03, 
		color=:black
	)
	fig
end

# ╔═╡ 0eecd302-185b-4467-8614-d25fa26a9b5d
let
	fig = Figure()
	ax = Axis(fig[1, 1], 
		xlabel="angular distance (degrees)",
		ylabel = L"$v_x - v_\textrm{los}$ / km\,s$^{-1}$",
		xgridvisible=false,
		ygridvisible=false,
	)

	errorscatter!(ϕ_pm, Measurements.value.(vz) .- Measurements.value.(rv), 
		#xerr=Measurements.uncertainty.(rv), yerr=Measurements.uncertainty.(vz),
		alpha=0.03, 
		color=:black
	)
	fig
end

# ╔═╡ e0111a22-8dd9-40d2-81e9-fcffe04adf4e
let
	fig = Figure()
	ax = Axis(fig[1, 1], 
		xlabel=L"$v_\textrm{los}$ / km\,s$^{-1}$",
		ylabel = L"$\delta v_x / \delta v_\textrm{los}$ ",
		xgridvisible=false,
		ygridvisible=false,
	)

	errorscatter!(Measurements.value.(rv), Measurements.uncertainty.(vz) 
 ./ Measurements.uncertainty.(rv), 
		alpha=0.03, 
		color=:black
	)
	fig
end

# ╔═╡ 79ffc43f-be70-40c7-9419-8666d33b5947
3*9 * sind(2)

# ╔═╡ aa1edde8-2b95-4e9b-b844-ba3150e81a0a
hist(rv_meas.RV_count, axis=(;yscale=log10, yticks=Makie.automatic))

# ╔═╡ 9cd82f65-ad8e-4e32-a9f9-c8096be2becd
mean(rv_meas.RV_count)

# ╔═╡ a7d7edca-e7c9-4187-95ab-0d3df2f8a75b
hist(filter(x->!ismissing(x), memb_stars.RV_err ./ memb_stars.RV_err_t23), bins=60)

# ╔═╡ 3b80b8ee-22b8-4fc7-b2f9-67274b494c03
  !ismissing(missing) && missing

# ╔═╡ 4adf3ac6-b11c-4391-8202-e5a67765e192
sqrt(1 / (1/0.65^2 + 1/1.0^2))

# ╔═╡ b1741955-2839-4afd-b381-997951030fd3
memb_stars[ .!ismissing.(memb_stars.RV_err_t23) .&& (memb_stars.RV_err .> memb_stars.RV_err_t23), [:RV, :RV_err, :RV_t23, :RV_err_t23, :RV_w09, :RV_err_w09, :RV_apogee, :RV_err_apogee]]

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
	errorscatter!(x, y, xerror=xerr, yerror=yerr, alpha=0.1)

	w = 1 ./ xerr .^2
	μ = LilGuys.mean(y, w)
	δμ = sqrt( 1 / sum(w))
	s = LilGuys.std(y, w)
	@info "μ = $μ ± $δμ"
	@info "σ = $s"
	hlines!(0, color=:black)
	return fig
end

# ╔═╡ a2370516-e7ec-4502-ae4b-b111bcf68d36
compare_rv_mean("apogee", memb_stars)

# ╔═╡ aaaf5ba2-c9ed-41ec-a22a-d78ed96fd84e
compare_rv_mean("w09", memb_stars)

# ╔═╡ 78e7aff0-3658-4d64-b117-58189a85307a
compare_rv_mean("t23", memb_stars)

# ╔═╡ 92bfc52c-05f5-45e3-ae26-90c579ecfe2c
hist(filter(not ∘ ismissing, memb_stars.RV_err_t23 ./ memb_stars.RV_err))

# ╔═╡ f4f9dd06-1a1a-458b-be75-05d52623580c
compare_rv_mean("gmos", )

# ╔═╡ 74152829-27ef-4d8d-8b32-ed30a18f30e4
rv_meas.RV_gmos[not.(ismissing.(rv_meas.RV_gmos))]

# ╔═╡ abbd2a53-e077-4af7-a168-b571e1a906b8
xlabel = L"radial velocity / km s$^{-1}$"

# ╔═╡ 5f71b8e9-5540-4431-8028-4ce14c8d7856
let 
	fig, ax = FigAxis(xlabel=xlabel)
	bins = make_bins(apogee_all.RV, calc_limits(apogee_all.RV), bandwidth=1)

	h = histogram(apogee.RV, bins, normalization=:none)
	scatter!(h)
	
	h = histogram(apogee_all.RV, bins, normalization=:none)

	scatter!(h, markersize=5)


	fig
end

# ╔═╡ 9d1495e8-5f8a-4892-b5b5-b22f3eb6ab7c
let 
	fig, ax = FigAxis(xlabel=xlabel)
	bins = make_bins(walker09_all.RV, calc_limits(walker09_all.RV), bandwidth=3)


	
	h = histogram(walker09_all.RV, bins, normalization=:none)

	lines!(h, label="all")
	
	h = histogram(walker09.RV, bins, normalization=:none)
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
	h = histogram(Float64.(sample.RV), normalization=:pdf)
	
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=:black)

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
	h = histogram(Float64.(memb_stars.RV), 30, normalization=:pdf)
	
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

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

# ╔═╡ fcf8f848-f8ce-41ab-9aa5-0507e90b28c7


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


	errorscatter!(x, y, yerror=yerr)
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


	errorscatter!(x, y, yerror=yerr)
	fig
end

# ╔═╡ ce3067d3-55e6-43a1-9b7b-4cf53c09ee88
md"""
# J24 purity checks
"""

# ╔═╡ 01e0cfcd-92af-4c04-a7d4-9c9d8fd5a9f1
rv_mean = 111.03; sigma_los = 9.6 # from mcmc modeling in next section

# ╔═╡ 03462abb-32d6-433c-a9c4-65ef2dd58965
function is_rv_member(rv, err)
	return abs(rv - rv_mean) / (err ⊕ sigma_los ⊕ 0.23) < 3
end

# ╔═╡ 82d4e049-7f47-4f8b-a3a0-7191ef82912b
filt_rv = is_rv_member.(rv_meas.RV, rv_meas.RV_err) 

# ╔═╡ e343dcd0-078f-476d-abe0-fd84001b72ae
mean(filt_rv[.!isnan.(rv_meas.RV_err)])

# ╔═╡ cc4a1b32-ab48-4381-bb53-4234f80638ad
rv_meas.RV_err .⊕ sigma_los

# ╔═╡ 0105a7ec-392e-426b-9bf9-a4b6144fe829
rv_meas.RV .- rv_mean

# ╔═╡ 7b22fa97-c311-42ae-ad5a-41d053c5334e
rv_meas.RV

# ╔═╡ e899faa9-580b-4aad-902e-382008048908
function purity_given_p(psatlow, psathigh)
	filt = psathigh .> rv_meas.PSAT .>= psatlow
	filt .&= .!isnan.(rv_meas.PSAT)
	filt .&= .!ismissing.(filt_rv)
	return sum(filt_rv[filt]) / sum(filt)
end

# ╔═╡ c0a4e207-9778-4d43-8197-5fc8887b2b82
function number_satisfying(psatlow, psathigh)
	filt = psathigh .> rv_meas.PSAT .>= psatlow
	filt .&= .!isnan.(rv_meas.PSAT)
	filt .&= .!ismissing.(filt_rv)

	return sum(filt)
end

# ╔═╡ 69161135-1f19-4d9a-8ff6-ae63a79e3eb5
purity_given_p(0.2, 1)

# ╔═╡ ab12510f-9393-4acd-8ed0-dcffa06c65e6
purity_given_p(0.0, 1.1)

# ╔═╡ 678a13d6-c892-43fe-a806-f3534661f785
let
	fig, ax = FigAxis(
		xlabel = "P SAT",
		ylabel = "purity"
	)

	ps = 0:0.1:1

	y = purity_given_p.(ps[1:end-1], ps[2:end])

	lines!(midpoints(ps), y)
	fig
end

# ╔═╡ 19f2cf18-7c90-4c66-bc81-94aea2ac8a3a
let
	fig, ax = FigAxis(
		xlabel = "P SAT",
		ylabel = "number in bin",
		yscale=log10,
		yticks = Makie.automatic,
	)

	ps = 0:0.1:1

	y = number_satisfying.(ps[1:end-1], ps[2:end])

	lines!(midpoints(ps), y)
	fig
end

# ╔═╡ 5bef0a45-2a4c-49bc-889f-b8bea315bebe
tanh(3)

# ╔═╡ 7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
md"""
## Binned properties
"""

# ╔═╡ 30e3dc5b-3ce6-4dd7-9c2a-c82774909a8c
sum(not.(ismissing.(memb_stars.RV_gmos)))

# ╔═╡ 33f54afc-cdb9-4eb8-887f-5a43281b837c
let
	fig = Figure(
		size=(5*72, 3*72)
	)
	ax = Axis(fig[1,1],
		xlabel = L"$R_\textrm{ell}$ / arcmin",
		ylabel = L"RV / km s$^{-1}$",
		#limits=(nothing, (60, 150))
	)

	#scatter!(memb_stars.r_ell, memb_stars.RV_dart, label="DART", markersize=5)
	scatter!(memb_stars.R_ell, memb_stars.RV_apogee, label="APOGEE", markersize=3)
	scatter!(memb_stars.R_ell, memb_stars.RV_gmos, label="Sestito+23", markersize=5)
	scatter!(memb_stars.R_ell, memb_stars.RV_w09, label="Walker+09", markersize=3)
	scatter!(memb_stars.R_ell, memb_stars.RV_t23, label="Tolstoy+23", markersize=3)

	Legend(fig[1,2], ax)
	#axislegend(ax, position=:rb)


	resize_to_layout!(fig)
	fig
end

# ╔═╡ dcb1bcf1-41e3-4389-b879-be6a53b5e670
mean(.!ismissing.(memb_stars.RV_gmos))

# ╔═╡ a70a2546-3aaa-4d78-a060-2a9a39940dc8
mean(.!ismissing.(memb_stars.RV_apogee))

# ╔═╡ 382e7503-0e31-49d5-b593-a6a1a0f9c5e0
mean(.!ismissing.(memb_stars.RV_t23))

# ╔═╡ 0b7107f8-38d2-4b0c-aad1-e07d3927cf87
mean(.!ismissing.(memb_stars.RV_w09))

# ╔═╡ e32f0d5a-9a6a-407c-a466-586c5d63fda8
mean(ismissing.(memb_stars.RV_t23) .& .!ismissing.(memb_stars.RV_w09))

# ╔═╡ 500815e1-f9e1-4401-ba2d-72f326cfa783
let
	fig = Figure(
		size=(5*72, 3*72)
	)
	
	ax = Axis(fig[1,1],
		xlabel = "R / arcmin",
		ylabel = L"RV / km s$^{-1}$",
		#limits=(nothing, (60, 150))
	)

	#scatter!(memb_stars.r_ell, memb_stars.RV_dart, label="DART", markersize=5)
	scatter!(rv_meas.R_ell, rv_meas.RV_apogee, label="APOGEE", markersize=3)
	scatter!(rv_meas.R_ell, rv_meas.RV_gmos, label="GMOS", markersize=3)
	scatter!(rv_meas.R_ell, rv_meas.RV_w09, label="Walker+09", markersize=3)
	scatter!(rv_meas.R_ell, rv_meas.RV_t23, label="tolstoy + 23", markersize=3)

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

	out_df[!, :vx] = Measurements.value.(vx)
	out_df[!, :vy] = Measurements.value.(vy)
	out_df[!, :vz] = Measurements.value.(vz)
	out_df[!, :vx_err] = Measurements.uncertainty.(vx)
	out_df[!, :vy_err] = Measurements.uncertainty.(vy)
	out_df[!, :vz_err] = Measurements.uncertainty.(vz)

	
	
	select!(out_df, Not([:o_Target_w09, ]))
	disallowmissing!(out_df)
end

# ╔═╡ 62301344-869a-4ec0-8299-29f0ff5d1c15
sum(sum(eachcol(ismissing.(out_df))))

# ╔═╡ 1deb1520-194a-40b5-8968-bf73b756ba3d
rv_meas[ .! ismissing.(rv_meas.RV_gmos), :PSAT]

# ╔═╡ 70290654-394f-4679-aaab-38345874e2e3
sum(sum(eachcol(ismissing.(rv_meas))))

# ╔═╡ f1912692-d890-4f97-ba98-7b226f29e9c8
sum(sum(eachcol(ismissing.(memb_stars))))

# ╔═╡ 4fdccf0a-adfa-4f05-bb48-df8a1efba940
out_df.:q__Fe_H__t23

# ╔═╡ 615f05fb-4907-40bc-9251-3065f565929b
let
	filename = "processed/sculptor_all_rv.fits"
	if isfile(filename)
		rm(filename); 
	end
	lguys.write_fits(filename, out_df)
	@info "wrote data"
end

# ╔═╡ Cell order:
# ╟─811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═257caa92-ccce-44ff-88ee-6a9c550ae5d9
# ╠═dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
# ╠═ef1bb6f5-00b8-405b-8702-244eca618644
# ╠═7ed5bcf5-dfc3-4e79-a608-d503124a1e96
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═da5a3b57-72bc-46e1-b1a0-6c02eb101626
# ╠═89729357-a2e7-4684-bde4-5b1684080f27
# ╠═3d8b3c55-153b-4a4a-9289-78f34df07abc
# ╠═d4123ffd-cb32-486e-a93d-f48d7112a831
# ╠═77e7884c-0360-4b7f-b9ad-e81be2516552
# ╠═248c2d5f-85cc-44be-bc63-43d4af470182
# ╠═7a50a176-96a5-4098-88d6-0fa2874d0f90
# ╟─dd1eee07-bb37-4de3-8781-47edd794c1f3
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
# ╠═0e14808d-df92-419b-b272-f3a08f4b86b1
# ╠═d4d0a488-1c0e-4bc6-9a88-94b27d84e8ce
# ╠═bb7f6769-ec92-460d-8423-449029175f79
# ╠═344f3c29-0873-4180-9025-51aeaeb2c681
# ╠═5f71b8e9-5540-4431-8028-4ce14c8d7856
# ╠═4d63430e-f59c-4c68-97e1-7eaa5679e55f
# ╠═bf088da8-a7c5-46e7-8e5f-ef5b91c10fb8
# ╠═1bc2541c-69ac-40aa-94f5-743707ca5a66
# ╠═63792b12-dde2-4361-88b2-1f5a789631ff
# ╟─27063de7-01d4-48a8-a06f-cc24aec662d2
# ╠═48c62811-136f-4962-a42c-b1dd1fc74f8c
# ╠═2e7ce524-573e-45c9-b0f4-ce9fea68e026
# ╠═f8775eb1-2cb9-4d81-8c02-39c43ceb9b45
# ╠═6009cc0e-6936-4589-85aa-dd2864948b7c
# ╠═b6de4afb-9362-4618-a114-df460031e4f9
# ╠═8e0095e7-2198-439d-bd19-976bdb1d3766
# ╠═b7345279-4f80-47ad-a726-537571849eae
# ╠═1151914d-b20f-4e46-8b25-e0d995e8b866
# ╟─c31f9e07-c520-443b-94dc-787519021d01
# ╠═15f2a8e2-90df-48a9-a7bf-e86955f566ce
# ╠═624af663-758f-4b68-8c71-5e5e324ad3b1
# ╠═b36c2a37-2359-4a1a-98fc-3a1cd17fd790
# ╠═1864bc12-cbc2-4fd4-84a8-6bb300f17c1b
# ╠═1d01f9a2-2b77-4ad0-9893-b31562de7924
# ╠═f064e833-5a6b-44be-9a5b-75599b39d150
# ╠═fd85bfca-1ed2-4324-8933-5070c5c426a0
# ╠═31c61c48-86d7-4db9-a479-eb25c058c281
# ╠═9e488fb1-d88b-48b3-b0d5-8d306cfd9f83
# ╟─180fac98-678c-4a14-966c-385387c60ac3
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
# ╠═baa5285a-edb5-4e98-bbe4-d40a5fbe8ab9
# ╠═1a6173db-6d5d-49c3-a2eb-9a77e591e747
# ╠═86ffdea6-5f02-446c-9dc8-b6c18aa835fb
# ╠═3026841e-3479-4dcb-ae0d-462166100c2b
# ╠═cf2c1c32-1fbc-4395-88cb-e1d9e410ef7d
# ╠═10478b34-2af4-4833-9dd4-03192b4f1576
# ╠═a678426b-f895-46f0-92fd-be4c35e68698
# ╠═6bdcb307-ee27-4a40-acbb-86787a4b5567
# ╠═9d1495e8-5f8a-4892-b5b5-b22f3eb6ab7c
# ╠═e1c2e02e-05de-4faf-af2f-f93759ddadfe
# ╠═5d5fbb43-ea1d-4b1c-a567-7143be1a9a5a
# ╠═a5a95eba-d282-4881-a84e-a25d4c83f114
# ╠═8c59d607-d484-497c-9ee7-751a4edd4992
# ╟─bbf49122-11b1-4272-a660-0437c6aa2b3f
# ╠═fe6fa910-f41e-4657-836b-7eda2f0cddb2
# ╠═e472cbb6-258e-4303-85e2-56f26358c97b
# ╠═0d2dbc73-1ded-46e3-b142-9bc7b777728d
# ╠═bd51fe42-7e39-4cd8-8065-58ab9814f966
# ╠═8f243944-3d0b-48ce-9f90-8d448c089239
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
# ╠═aa8ccf4d-fde4-4a73-a93c-8c11ff747da7
# ╠═aa9b5b6f-f187-498e-b86d-f3b8a16a73ba
# ╠═30461ef6-1778-4612-b55d-8d16f0c6db12
# ╠═ca3c2080-af1b-4218-a29e-26eed51f0ef8
# ╠═f61debe4-8b23-4415-b77c-43c4465ccfcb
# ╠═102c73ef-2c95-4784-80df-ed0312511c00
# ╠═d3333b48-aa4e-42c1-9e0a-bbff98e3647d
# ╠═8bc140fa-6f4e-4ec5-bc95-989fc0ea48d1
# ╠═d3ef6a98-8770-4359-9e37-94f08d2d9aa5
# ╠═78d39e37-44b5-4171-ab45-0672f607b517
# ╠═88ed48b4-baa9-4ac0-86e1-8348edcd59b4
# ╠═222bb254-8b65-44d3-b3d2-b08fbcbe9950
# ╠═e3f05ee2-cc5f-437e-801d-3c7d842af709
# ╠═a13a33f1-65c4-4217-8636-639b1e14d109
# ╠═ee3c22db-6b6b-4c30-8d1f-86b103c018fc
# ╠═6bc02c4f-31a2-4e4e-8612-e66f8cc9c93e
# ╠═11fcf4f8-fcd5-4426-a4e9-b046138bde1b
# ╠═1798bdc1-1770-4dd1-a7ed-5b105509a31d
# ╠═877706f1-08a0-496a-890d-622e3a2fd9ec
# ╠═8c1457d6-32d4-4727-b45e-1a46244708a6
# ╠═63e28d21-c1ab-413e-ac7e-e6776bcf6ca5
# ╠═dbb51473-6506-4119-a11d-50df2ceb6dcf
# ╠═d913bf9b-c03b-4bce-a39d-7d23eca9b355
# ╠═24b2455e-1c9a-4b9a-9d88-19b98d4b2165
# ╠═73bd1553-e2ae-4bfb-aac1-0880346f5054
# ╠═01b3135f-7a72-4669-b586-4bc5894464ad
# ╠═d11edca7-b9ae-4269-9e1b-661d59bd965e
# ╠═c724e720-18ca-4905-a4b6-39fc47abe39d
# ╠═733fe42e-b7a5-4285-8c73-9a41e4488d40
# ╠═74b10a3e-1342-454f-8eed-77b371f81edf
# ╠═d8800a31-1ed3-422f-ac51-90f18cf61c29
# ╠═015f2345-74d4-4296-b0dc-35d6e13ad8cc
# ╠═08e46c2e-d7c6-4bdc-9243-e18e6b9e9dc2
# ╠═cc6c65db-ef57-4745-8ada-e11427274a77
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╠═a1938588-ca40-4844-ab82-88c4254c435b
# ╟─ab27964f-b070-4328-81c2-6ecf4b8cec1e
# ╠═23d4b7c4-61f9-4748-b3c4-80eaad551914
# ╠═7dab615b-e0a1-431e-88c0-e1e9d9c29568
# ╠═57a77527-a14c-47a2-8949-1cbdafc49994
# ╠═5930b133-3e3c-43be-9ed6-f125ea4418d4
# ╠═c4859522-5f50-4042-aa24-6bf9b214cb22
# ╠═ad7ae1a8-c5a4-4d79-bd86-35431fa71856
# ╠═d2c28cfd-652d-41fd-8ae2-4b5454761377
# ╠═c0accd82-17a7-437c-9ba8-4db437071a5b
# ╠═490b14a7-719d-4123-9173-e620288ff146
# ╠═9ae17758-a548-48da-9bc2-a1b4f05e1818
# ╠═e06b8dc3-5e14-4c8a-80a2-e187de4fe0aa
# ╠═383a9a9c-59da-4d3d-9cb9-253d9c685057
# ╠═7bf999a8-86b7-4e9f-b9c2-aa9de943ccad
# ╠═4180e4bd-944c-4817-a36d-35bc98a19bed
# ╠═d0c9ad40-4bdc-4269-a1ab-3e501c84ce13
# ╠═3f51257c-0cb6-4825-a41f-21fc6757060e
# ╠═d0842a73-dd1b-452f-97f9-50a23668fe14
# ╠═11461563-dbcd-48fc-a97a-1a0391538462
# ╠═718ec8bc-cae6-4905-a336-b04964699b61
# ╠═379bff0e-e65c-4abd-bce1-d6b776656bc8
# ╠═6c1aad21-80cc-4dfb-80d5-4c1de4ecfe23
# ╠═65e5ecd2-f184-40cf-850d-326683560784
# ╠═596bb8fc-b5f7-413e-bdef-8bd1968ea7cb
# ╠═9492f559-2a37-4d86-8247-2d217a387d8c
# ╠═cc68349e-da2f-4a73-92f5-6a021237accc
# ╠═ae1f5ac1-265a-4f81-abe8-38d882b494bc
# ╠═f61c0c4c-446f-466f-8f64-95250a9177f1
# ╠═64ea3168-59a4-4db2-8fec-f7facf4f0df6
# ╠═0ef130e2-14b0-4fc6-9b0c-76ab18f3775d
# ╠═58e5855a-9279-4f1c-bc63-52d471e63389
# ╠═6f7a98dd-77cf-42a0-888d-8e34e2a0490f
# ╠═13f3ac83-2fb3-4537-a3e6-bdb83d8ba72a
# ╠═0eecd302-185b-4467-8614-d25fa26a9b5d
# ╠═e0111a22-8dd9-40d2-81e9-fcffe04adf4e
# ╠═79ffc43f-be70-40c7-9419-8666d33b5947
# ╠═aa1edde8-2b95-4e9b-b844-ba3150e81a0a
# ╠═9cd82f65-ad8e-4e32-a9f9-c8096be2becd
# ╠═a7d7edca-e7c9-4187-95ab-0d3df2f8a75b
# ╠═3b80b8ee-22b8-4fc7-b2f9-67274b494c03
# ╠═4adf3ac6-b11c-4391-8202-e5a67765e192
# ╠═b1741955-2839-4afd-b381-997951030fd3
# ╟─3655b6a0-ac9a-4a49-86fd-6a6484409819
# ╠═7f13339e-6a0c-4944-8a36-5ab136fd8415
# ╠═a2370516-e7ec-4502-ae4b-b111bcf68d36
# ╠═aaaf5ba2-c9ed-41ec-a22a-d78ed96fd84e
# ╠═78e7aff0-3658-4d64-b117-58189a85307a
# ╠═92bfc52c-05f5-45e3-ae26-90c579ecfe2c
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
# ╠═fcf8f848-f8ce-41ab-9aa5-0507e90b28c7
# ╠═eaf47845-95dd-4f7d-bf76-3d9b1154711a
# ╟─b5faa873-d5c9-447b-b90e-f693db26e6c2
# ╠═46434fa6-09df-4b02-9270-cbdcc9648e38
# ╠═dc4c0453-fed0-4733-9f59-0b2df506b45c
# ╠═ce3067d3-55e6-43a1-9b7b-4cf53c09ee88
# ╠═01e0cfcd-92af-4c04-a7d4-9c9d8fd5a9f1
# ╠═03462abb-32d6-433c-a9c4-65ef2dd58965
# ╠═82d4e049-7f47-4f8b-a3a0-7191ef82912b
# ╠═e343dcd0-078f-476d-abe0-fd84001b72ae
# ╠═cc4a1b32-ab48-4381-bb53-4234f80638ad
# ╠═0105a7ec-392e-426b-9bf9-a4b6144fe829
# ╠═7b22fa97-c311-42ae-ad5a-41d053c5334e
# ╠═e899faa9-580b-4aad-902e-382008048908
# ╠═c0a4e207-9778-4d43-8197-5fc8887b2b82
# ╠═69161135-1f19-4d9a-8ff6-ae63a79e3eb5
# ╠═ab12510f-9393-4acd-8ed0-dcffa06c65e6
# ╠═678a13d6-c892-43fe-a806-f3534661f785
# ╠═19f2cf18-7c90-4c66-bc81-94aea2ac8a3a
# ╠═5bef0a45-2a4c-49bc-889f-b8bea315bebe
# ╟─7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═30e3dc5b-3ce6-4dd7-9c2a-c82774909a8c
# ╠═33f54afc-cdb9-4eb8-887f-5a43281b837c
# ╠═dcb1bcf1-41e3-4389-b879-be6a53b5e670
# ╠═a70a2546-3aaa-4d78-a060-2a9a39940dc8
# ╠═382e7503-0e31-49d5-b593-a6a1a0f9c5e0
# ╠═0b7107f8-38d2-4b0c-aad1-e07d3927cf87
# ╠═e32f0d5a-9a6a-407c-a466-586c5d63fda8
# ╠═500815e1-f9e1-4401-ba2d-72f326cfa783
# ╠═725b9c7d-93f9-48c3-b97b-8a1147a16f78
# ╠═62301344-869a-4ec0-8299-29f0ff5d1c15
# ╠═1deb1520-194a-40b5-8968-bf73b756ba3d
# ╠═70290654-394f-4679-aaab-38345874e2e3
# ╠═f1912692-d890-4f97-ba98-7b226f29e9c8
# ╠═4fdccf0a-adfa-4f05-bb48-df8a1efba940
# ╠═615f05fb-4907-40bc-9251-3065f565929b
