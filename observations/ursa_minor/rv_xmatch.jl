### A Pluto.jl notebook ###
# v0.20.5

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

	import PythonCall
	import LilGuys as lguys
end

# ╔═╡ d4123ffd-cb32-486e-a93d-f48d7112a831
include(ENV["DWARFS_ROOT"] * "/utils/gaia_filters.jl")

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Radial Velocity CrossMatch

This notebook loads in radial velocity data from several sources (including Tolstoy+23, APOGEE, DART, and Fedrico++'s GMOS observations) and consolidates this data into a single combined radial velocity catalogue cross-matched with J+24's membership calculations.

With this catalogue, detailed radial velocity analysis can then be calculated.
"""

# ╔═╡ 12e0e7ad-925b-47dc-a1d8-c56e65ffaa4c


# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, kurtosis, sem

# ╔═╡ dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
not = !

# ╔═╡ ef1bb6f5-00b8-405b-8702-244eca618644
import DensityEstimators: histogram, calc_limits, make_bins

# ╔═╡ 923674f9-a5fa-454d-a5e2-a4d30d047a11
import TOML

# ╔═╡ 68c13192-58c0-4e77-9f86-567d55a036b5
rv_low = -280

# ╔═╡ 17629b76-2cb8-4a56-9687-3a630950872b
rv_high = -210

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

# ╔═╡ 7a50a176-96a5-4098-88d6-0fa2874d0f90
j24 = lguys.read_fits("data/jensen+24_2c.fits")

# ╔═╡ c470c2b9-093d-42ab-b96d-9dac231ccabc
md"""
## Apogee
"""

# ╔═╡ bb7f6769-ec92-460d-8423-449029175f79
apogee_all = lguys.read_fits("$data_dir/apogee_xmatch.fits")

# ╔═╡ 4e8b7258-76fa-404f-91ca-8cfda241a332
sum(apogee_all.RV_flag)

# ╔═╡ 4ebb0411-1264-4297-be0c-cacfdb2fd59d
filt_apogee_good = apogee_all.RV_sigma ./ apogee_all.RV_err ./ sqrt.(apogee_all.RV_count) .< 5

# ╔═╡ 14f866c4-8c65-461f-9834-329b281a3c0a
sum(.!filt_apogee_good)

# ╔═╡ c941ef83-42b3-4bdf-b1db-6f0b9b3f1cd0
apogee_good = apogee_all[filt_apogee_good, :]

# ╔═╡ bf088da8-a7c5-46e7-8e5f-ef5b91c10fb8
filt_apogee = sigma_clip(apogee_all.RV, 3)

# ╔═╡ 4d63430e-f59c-4c68-97e1-7eaa5679e55f
apogee = apogee_all[filt_apogee, :]

# ╔═╡ 1bc2541c-69ac-40aa-94f5-743707ca5a66
hist(apogee.RV)

# ╔═╡ 27063de7-01d4-48a8-a06f-cc24aec662d2
md"""
### Check apogee-dart xmatch against fed
"""

# ╔═╡ b6de4afb-9362-4618-a114-df460031e4f9
md"""
## Graces alla Federico
"""

# ╔═╡ b7345279-4f80-47ad-a726-537571849eae
let
	
	global graces = CSV.read("$data_dir/sestito+23_graces.csv", DataFrame)


	filt, idx = xmatch(graces, j24)
	println(filt)
	graces[!, :source_id] = j24.source_id[idx]

end

# ╔═╡ c31f9e07-c520-443b-94dc-787519021d01
md"""
## Spencer et al. 2018
"""

# ╔═╡ 15f2a8e2-90df-48a9-a7bf-e86955f566ce

begin 
	spencer18_all = lguys.read_fits("$data_dir/umi_rv_spencer+18.fits")
	rename!(spencer18_all, 
		:e_RV => :RV_err,
		:RAJ2000 => :ra,
		:DEJ2000 => :dec,
	)
	
	filt, idx = xmatch(spencer18_all, j24, 2)
	println("matched / missed\t ", sum(filt), "\t", sum(.!filt))
	
	spencer18_all[!, :source_id] = j24.source_id[idx]

	
end

# ╔═╡ c50824f9-652f-47e5-9a05-fa052c36a247
spencer18_all

# ╔═╡ 7922960e-a02a-4c2b-b7f0-0aca1d514626
unique(spencer18_all.source_id)

# ╔═╡ c405ede3-258c-4b3d-8c4e-e4bab68f3bf0
hist(spencer18_all.RV)

# ╔═╡ b48ea14d-48d7-4953-8d6a-b84bce1a4788
std(spencer18_all.RV)

# ╔═╡ 180fac98-678c-4a14-966c-385387c60ac3
md"""
## Pace et al. 2020
"""

# ╔═╡ 19d8aa00-14d1-4270-8e8c-9b8907779e3b
begin
	pace20 = lguys.read_fits("$data_dir/pace_xmatch.fits")

	rename!(pace20, Dict(
		"vLOSu" => "RV",
		"e_vLOSu" => "RV_err",
	))

	pace20.ra = pace20.RAJ2000
	pace20.dec = pace20.DEJ2000
	
	pace20.source_id[ismissing.(pace20.source_id)] .= -pace20.ID[ismissing.(pace20.source_id)]
	pace20
end

# ╔═╡ d415ba9e-75c0-4192-8761-dbda70e43a73
filt_p20, idx_p20 = xmatch(j24, pace20)

# ╔═╡ 1e6ff00b-eee8-4e03-8b54-e27d6cbb6f16
round(Int, pace20.source_id[5] )

# ╔═╡ 3056ca82-9ffe-4065-b207-ad958aaab08b
pace20.source_id[1] 

# ╔═╡ 0313553f-7a05-492b-99b9-9648e98baad4
pace20[pace20.ID .== 721, :].ra

# ╔═╡ dd136efd-ed53-4472-ac73-2cd97a4cdad2
length(unique(pace20.GaiaDR2))

# ╔═╡ 6f1ff560-c458-48cd-8d7f-9f3dd7631b08
sum(pace20.GaiaDR2 .> 0)

# ╔═╡ 26092ed5-6283-4bf9-8ee1-60a8cf15e1ef
length(unique(pace20.source_id))

# ╔═╡ 1c6074eb-bafa-439c-aa55-ee7a5c3a58d5
length(unique(pace20.ID))

# ╔═╡ cdf27cbc-fd37-412f-a926-e1f86f836bb8
duplicate_ids = pace20[nonunique(pace20, ["source_id"]), :source_id]

# ╔═╡ 36cdb5ff-65a3-4c5c-b274-4d7da53d2ab4
hist(pace20.RV)

# ╔═╡ 1ca9b968-87d0-46cc-a1fe-634717ac7446
md"""
# Xmatched
"""

# ╔═╡ fe6fa910-f41e-4657-836b-7eda2f0cddb2
function add_xmatch(df, new, suffix; cols=["ra", "dec", "RV", "RV_err"])
	copynew = copy(new[:, ["source_id"; cols]])
	
	rename!(copynew, Dict(
		col => "$(col)_$suffix" for col in cols
	)
	)
	
	result = rightjoin(df, copynew, on="source_id")

	result[!, "study"] .= suffix

	@assert size(new, 1) == size(result, 1)

	return result
end

# ╔═╡ 842e2c79-5eee-4bf3-b86d-eda1204336d7
unique(apogee_all.source_id)

# ╔═╡ bc4fd1bb-fb30-4581-8952-891293476345
unique(pace20.source_id) == pace20.source_id

# ╔═╡ c235cd75-cada-4d8f-a380-43af619eaa24
sum(.!isnan.(j24.dr2_radial_velocity_error))

# ╔═╡ a85d8d3f-197a-4f7f-97cd-5e2523a885d7
let
	global rv_gaia =  j24[.!isnan.(j24.dr2_radial_velocity_error), :]
	rv_gaia[:, :RV] = rv_gaia.dr2_radial_velocity
	rv_gaia[:, :RV_err] = rv_gaia.dr2_radial_velocity_error

end

# ╔═╡ b94b7591-018b-47d7-8af8-be0c4b3ea49f


# ╔═╡ e472cbb6-258e-4303-85e2-56f26358c97b
let
	global all_stars 

	all_stars = copy(j24[1:0, :])
	append!(all_stars, add_xmatch(j24, apogee_good, "apogee"), cols = :union)
	append!(all_stars, add_xmatch(j24, pace20, "p20", cols=["ra", "dec", "RV", "RV_err", "ID","PMR"]), cols = :union)
	append!(all_stars, add_xmatch(j24, spencer18_all, "s18"), cols = :union)
	append!(all_stars, add_xmatch(j24, rv_gaia, "gaia"), cols=:union)
	append!(all_stars, add_xmatch(j24, graces, "s24"), cols=:union)


	filt_pace = (.!ismissing.(all_stars.ID_p20)) .&& ismissing.(all_stars.ra)
	all_stars[filt_pace, :ra] = all_stars.ra_p20[filt_pace]
	all_stars[filt_pace, :dec] = all_stars.dec_p20[filt_pace]
	all_stars[filt_pace, :pmra] .= NaN
	all_stars[filt_pace, :pmra_error] .= NaN
	all_stars[filt_pace, :pmdec] .= NaN
	all_stars[filt_pace, :pmdec_error] .= NaN
	all_stars[filt_pace, :PSAT] = all_stars.PMR_p20[filt_pace]

	all_stars[!, :star_id] .= ""

	# pace 2020 includes stars not in gaia
	in_gaia = abs.(all_stars.source_id) .> 10000

	all_stars[in_gaia, :star_id] = ["GaiaDR3 $i" for i in all_stars.source_id[in_gaia]]

	
	all_stars[filt_pace, :star_id] = ["pace $i" for i in all_stars.ID_p20[filt_pace]]

	all_stars[:, "RV"] = [row["RV_$(row.study)"] for row in eachrow(all_stars)]
	all_stars[:, "RV_err"] = [row["RV_err_$(row.study)"] for row in eachrow(all_stars)]

	add_xi_eta!(all_stars, obs_properties["ra"], obs_properties["dec"])

	all_stars[:, "R_ell"] = lguys.calc_R_ell(all_stars.xi, all_stars.eta, obs_properties["ellipticity"], obs_properties["position_angle"])
	all_stars
end

# ╔═╡ 51468818-3e4e-4445-b23f-4df757c80478
unique(all_stars.study)

# ╔═╡ ac200cf0-714c-4b66-855a-e204735c7c04
sum(ismissing.(all_stars.study))

# ╔═╡ 0c8bceff-8ec3-471e-ac1b-e0d6caecfa66
all_stars[ismissing.(all_stars.ra) .&& ((abs.(all_stars.source_id) .> 10000)), :].source_id

# ╔═╡ 28bf5ad2-26f0-4a25-890f-2f8c80f9e6e6
all_stars[(ismissing.(all_stars.ra)), :ID_p20]

# ╔═╡ eec5c292-b7d0-49cf-83dc-283c400b9c05
"RV" ∈ names(rv_gaia)

# ╔═╡ 3e37a9dd-fa57-45f8-999e-e499b760a420
j24[.!isnan.(j24.dr2_radial_velocity), :]

# ╔═╡ 14d30389-c25c-426c-99f5-68fe98aabc7f
all_stars[ismissing.(all_stars.RV), :].dec_p20

# ╔═╡ 85511c49-479f-44a0-b88d-b33ffef99e0a
size(all_stars, 1) - length(unique(all_stars.star_id)) 

# ╔═╡ 29e0bd8d-d66f-44d5-81af-b0f54ec060e6
sum(all_stars.star_id .== "")

# ╔═╡ 51b3371a-87ef-4c10-a312-d5fa00ff7844
pace20

# ╔═╡ 2eba0d13-c688-4a70-95e4-59bc0d21020b
sum(ismissing.(all_stars.ra))

# ╔═╡ bd51fe42-7e39-4cd8-8065-58ab9814f966
sum(not.(ismissing.(all_stars.RV_apogee))) == length(apogee_all.RV)

# ╔═╡ 89552b84-d12e-4d88-a58e-8b89ad4b2569
md"""
# Validation for xmatch
"""

# ╔═╡ ff236dd1-3097-44f6-8e12-00803a8bee42
md"""
Each of these plots should exactly line up within the xmatch tolerance.
"""

# ╔═╡ e6f2de3b-ce32-4d61-851f-4e42fcce95c0
function plot_xmatch_radec(suffix)

	ra2 = all_stars[:, "ra_$suffix"]
	filt = not.(ismissing.(ra2)) .& (all_stars.source_id .> 0)
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

# ╔═╡ 7a77bf54-fff8-4906-b738-0d3c91bdf6aa
sum(ismissing.(all_stars.source_id))

# ╔═╡ 5b2fceff-9c3e-472d-9310-31e920137e41
plot_xmatch_radec("apogee")

# ╔═╡ 4b305b83-1a3b-48a6-b19f-6f3ebed0768f
plot_xmatch_radec("p20")

# ╔═╡ 6ed9fabe-1bdd-4ecf-90b4-084d5296150e
plot_xmatch_radec("s18")

# ╔═╡ f7213298-3dcb-49ac-a9f0-a129a03423aa
plot_xmatch_radec("s24")

# ╔═╡ a4fd8597-8706-45aa-a8c2-f72d90cfbbfd
pace20.ra[1]

# ╔═╡ f7ec8bba-9f45-435b-b67c-33182e992dfd
md"""
# Cross study RV
"""

# ╔═╡ 93d185f2-1eaa-4d35-87dd-b84f385483de


# ╔═╡ 66e5693f-90cc-4b49-8a0c-8d92ea9e8cdd
function get_matching_stars(study1, study2)
	rv_col_1 = "RV_$study1"
	rv_col_2 = "RV_$study2"
	
	df1 = all_stars[.!ismissing.(all_stars[:, rv_col_1]), ["ra", "dec", "star_id", "RV_$study1", "RV_err_$(study1)"]]
	df2 = all_stars[.!ismissing.(all_stars[:, rv_col_2]), ["ra", "dec", "star_id", "RV_$study2", "RV_err_$(study2)"]]

	df = innerjoin(df1, df2, on="star_id", makeunique=true)
end

# ╔═╡ ce62265b-e411-4e0a-942c-281f82c03c6e
df_p20_apogee = get_matching_stars("p20", "apogee")

# ╔═╡ ee4d33fd-1bd8-41e8-8dc3-17092e4c267d
all_stars[all_stars.star_id .== "GaiaDR3 1645376716690536192", [:ra, :dec, :RV, :RV_err, :study,]]

# ╔═╡ f3ab19a6-95c6-4756-b22c-c7abfcf330c1
println(all_stars[all_stars.star_id .== "GaiaDR3 1645376716690536192", :ra_p20])

# ╔═╡ a49845c5-0b66-4316-bf25-febb57b77a51
println(all_stars[all_stars.star_id .== "GaiaDR3 1645376716690536192", :dec])

# ╔═╡ 02921108-b708-4a90-b7f9-6903632da601
println(all_stars[all_stars.star_id .== "GaiaDR3 1645376716690536192", :dec_p20])

# ╔═╡ 9213cc36-74d1-452f-bd9a-eb5c5cab1f87
function compare_rv(study1, study2; limits=(nothing, nothing, nothing, nothing))
	rv_col_1 = "RV_$study1"
	rv_col_2 = "RV_$study2"
	
	df1 = all_stars[.!ismissing.(all_stars[:, rv_col_1]), ["ra", "dec", "star_id", "RV_$study1", "RV_err_$(study1)"]]
	df2 = all_stars[.!ismissing.(all_stars[:, rv_col_2]), ["ra", "dec", "star_id", "RV_$study2", "RV_err_$(study2)"]]

	df = innerjoin(df1, df2, on="star_id", makeunique=true)

	if size(df, 1) == 0
		println("nothing to plot")
		return
	end
	
	println("plotting ", size(df,1), " stars")


	fig, ax = FigAxis(
		xlabel = "RV $study1",
		ylabel = "RV $study2",
		aspect=DataAspect(),
		limits=limits
	)

	errorscatter!(df[:, rv_col_1], df[:, rv_col_2], xerror=df[:, "RV_err_$(study1)"], yerror=df[:, "RV_err_$(study2)"])

	ll = [min(minimum(df[:, rv_col_1]), minimum(df[:, rv_col_2])),
		max(maximum(df[:, rv_col_1]), maximum(df[:, rv_col_2])),
		]
	lines!(ll, ll, color=:black)
	return fig
end

# ╔═╡ db06c81e-1456-4e12-a95a-c87bb23eaec3
Base.summarysize(all_stars)

# ╔═╡ f9e1bcbd-2def-4e64-b69e-d1cfcf7ffedd
compare_rv("s18", "p20")

# ╔═╡ f61debe4-8b23-4415-b77c-43c4465ccfcb
compare_rv("p20", "apogee")

# ╔═╡ 3ad3bd41-56f4-4fa7-8345-c275303c637d
compare_rv("p20", "apogee", limits=(-280, -200, -280, -200))

# ╔═╡ 102c73ef-2c95-4784-80df-ed0312511c00
compare_rv("apogee", "s18")

# ╔═╡ 3a9841cc-4e19-4962-8b5f-dd25b223189e
compare_rv("s18", "apogee", limits=(-290, -200, -290, -200))

# ╔═╡ cd256d12-5c9e-436b-9473-ed3e7e713225
compare_rv("apogee", "gaia")

# ╔═╡ dfb4db7f-1d27-4e0b-9eeb-0afc2ca3d92a
compare_rv("p20", "gaia")

# ╔═╡ 578b9efe-9cd5-4b55-a321-4eea37bc43d1
md"""
# Weighted Mean
"""

# ╔═╡ f4505ad1-6026-4f07-b849-6f55401f2d80
"""
	filt_missing(vector; verbose, low, high)

Returns a filter excluding NaN and missing values. Additionally, 
can truncate values past low and high (defaulting to Infs).
"""
function filt_missing(col; verbose=false, low=-Inf, high=Inf)
	filt =  not.(ismissing.(col)) .& not.(isnan.(col))
	filt1 = high .> col .> low
	if verbose
		println("excluding ", sum(not.(filt1)[filt]), " outliers")
	end
	return filt .& filt1
end

# ╔═╡ 36d2d86f-6a75-46f4-b48f-36137e95e90d
filt_missing(all_stars.RV_gaia)

# ╔═╡ 8bc140fa-6f4e-4ec5-bc95-989fc0ea48d1
"""
	safe_weighted_mean(values, errors)

Given a set of values with standard errors, computes
the weighted mean between non-missing (and non-nan) values,
returning the inverse varience weighted mean and the standard deviation.
"""
function safe_weighted_mean(values, errors)
	filt = filt_missing(values)
	filt .&= filt_missing(errors)
	if sum(filt) == 0
		return missing
	end

	x = float.(values[filt]) .± float.(errors[filt])

	return weightedmean(x), std(values[filt])
end

# ╔═╡ db9a60bf-a7db-410d-ada2-0b6ed011e2c4


# ╔═╡ 0d0c854c-b794-469d-8dba-2bf83a9f6d2f
length(unique(all_stars.star_id))

# ╔═╡ 36f74c41-4190-4705-a656-61a1223ffed8
names(all_stars)

# ╔═╡ 07070037-e64a-46eb-a23a-4cc7e35b4551
unique(all_stars.study)

# ╔═╡ 7b9a5d6d-48b0-4426-8145-19d5a414af8d
let
		filt = .!isnan.(all_stars.RV_err)
	filt[ismissing.(all_stars.RV)] .= false
	filt = disallowmissing(filt)
	all_stars[filt, :]
end

# ╔═╡ d867234f-dfc2-436c-b83d-02e87089897a
sum(ismissing.(all_stars.RV))

# ╔═╡ 29f47769-8d2c-4e9d-b5eb-79b168580a04
sum(isnan.(all_stars.RV[.!ismissing.(all_stars.RV)]))

# ╔═╡ 3655088f-8ae3-4edc-899d-466d909e5e6c
groupby(all_stars, "star_id")[1]

# ╔═╡ d3333b48-aa4e-42c1-9e0a-bbff98e3647d
all_studies = ["gaia", "apogee", "p20", "s18"]

# ╔═╡ 77130b64-383e-4b3c-bb55-76f925db9986
function average_by_star(all_stars; id_col="star_id", study_rv_names=all_studies)
	averaged_stars = all_stars[1:0, :]

	filt = .!isnan.(all_stars.RV_err)
	filt[ismissing.(all_stars.RV)] .= false
	filt = disallowmissing(filt)

	println("excludes $(sum(.!filt)) missings")
	for group in groupby(all_stars[filt, :], id_col)
		if all(isnan.(group.RV))
			println("skiped nan")
			continue
		end

		μ, σ = safe_weighted_mean(group.RV, group.RV_err)

		newrow = copy(group[1:1, :])
		newrow[:, "RV"] .= μ.val
		newrow[:, "RV_err"] .= μ.err
		newrow[:, "RV_std"] .= σ
		newrow[:, "RV_count"] .= size(group, 1)

		for study in study_rv_names
			ave = safe_weighted_mean(group[:, "RV_$study"], group[:, "RV_err_$(study)"])
			if ismissing(ave)
				μ = NaN ± NaN
				σ = NaN
			else
				μ, σ = ave
			end
			newrow[:, "RV_$study"] .= μ.val
			newrow[:, "RV_err_$(study)"] .= μ.err
			
			newrow[:, "RV_std_$(study)"] .= σ
			newrow[:, "RV_count_$(study)"] .= sum(.!isnan.(group[:, "RV_$study"]))
		end


		append!(averaged_stars, newrow, cols=:union)
	end

	averaged_stars
end


# ╔═╡ 782691d3-dd01-4a34-9db2-a18b0b6566e1
averaged_stars = average_by_star(all_stars)

# ╔═╡ 2365da83-7be8-44ec-a399-d520d610014b
hist(averaged_stars.RV)

# ╔═╡ 828cf7c4-d91e-4d02-8649-dd0194177f67
hist(asinh.(averaged_stars.RV_std[averaged_stars.RV_count .> 1] ./ averaged_stars.RV_err[averaged_stars.RV_count .> 1]))

# ╔═╡ 2078b2ab-e973-40fc-b19c-83ae860b626d
hist(averaged_stars.RV_count, axis=(;yscale=log10))

# ╔═╡ 1f40dc2a-0a93-4b41-a550-1eb5b77fc037
averaged_stars[argmax(averaged_stars.RV_count), :].star_id

# ╔═╡ 0008053d-77ec-4b74-a701-1e584254ed11
sum(averaged_stars.RV_count .> 1)

# ╔═╡ 88ed48b4-baa9-4ac0-86e1-8348edcd59b4
begin 
	rvs = [all_stars[:, "RV_$study"] for study in all_studies]
	rv_errs = [all_stars[:, "RV_err_$study"] for study in all_studies]

	rvs = hcat(rvs...)
	rv_errs = hcat(rv_errs...)
end

# ╔═╡ 222bb254-8b65-44d3-b3d2-b08fbcbe9950
all_studies

# ╔═╡ 6bc02c4f-31a2-4e4e-8612-e66f8cc9c93e
rvs

# ╔═╡ 11fcf4f8-fcd5-4426-a4e9-b046138bde1b
rv_errs

# ╔═╡ 877706f1-08a0-496a-890d-622e3a2fd9ec
RV_a_std = [safe_weighted_mean(rvs[i, :], rv_errs[i, :]) for i in 1:size(rvs, 1)]

# ╔═╡ 733fe42e-b7a5-4285-8c73-9a41e4488d40
begin 
	memb_filt = averaged_stars.PSAT .> 0.2
	memb_filt .&= -276 .< averaged_stars.RV .< -216

	quality_score = averaged_stars.RV_std ./ averaged_stars.RV_err
	quality_score[isnan.(averaged_stars.RV_std)] .= 0.0

	memb_filt .&= quality_score .< 5
end

# ╔═╡ 74b10a3e-1342-454f-8eed-77b371f81edf
hist(asinh.(quality_score), axis=(; yscale=log10))

# ╔═╡ d8800a31-1ed3-422f-ac51-90f18cf61c29
sum(quality_score .< Inf)

# ╔═╡ 015f2345-74d4-4296-b0dc-35d6e13ad8cc
sum(quality_score .< 5)

# ╔═╡ cc6c65db-ef57-4745-8ada-e11427274a77
memb_stars = averaged_stars[memb_filt, :]

# ╔═╡ 7f4a5254-ed6f-4faa-a71e-4b4986a99d45
hist(memb_stars.RV)

# ╔═╡ a1938588-ca40-4844-ab82-88c4254c435b
length(memb_stars.RV)

# ╔═╡ 1823f690-81c9-471b-bdec-c9fe261544c1
rv_test = [1.2, 1.5, 0.3, 0.9]

# ╔═╡ 8a5615a0-c34c-4de5-8a62-543cdb613988
rv_test_err = [0.01, 0.03, 0.025, 0.05]

# ╔═╡ 8fb2fe7f-9f74-4cf6-9cef-4bac3ebc9ae4
w_test = 1 ./ rv_test_err .^ 2

# ╔═╡ a7071b79-f36d-453a-9f9e-2a16832f0785
lguys.std(rv_test, w_test)

# ╔═╡ c96b5112-7f07-4e7c-b42b-18b26055561a
lguys.mean(rv_test, w_test)

# ╔═╡ 2be7ce8f-f44b-4e57-b8a6-53b623376c4a
1 / sqrt(sum(w_test)) # mean uncertainty

# ╔═╡ 4ed42682-b81c-4635-b5a9-7a73857cc927
safe_weighted_mean(rv_test, rv_test_err)

# ╔═╡ ab27964f-b070-4328-81c2-6ecf4b8cec1e
md"""
# Corrected coordinate frame
"""

# ╔═╡ 23d4b7c4-61f9-4748-b3c4-80eaad551914
distance = obs_properties["distance"] ± obs_properties["distance_err"]

# ╔═╡ 7dab615b-e0a1-431e-88c0-e1e9d9c29568
ra0, dec0 = obs_properties["ra"], obs_properties["dec"]

# ╔═╡ c0accd82-17a7-437c-9ba8-4db437071a5b
icrs = [
	lguys.ICRS(
		ra=row.ra ± 0, dec=row.dec ± 0, 
		distance=82.3 ± 2, pmra=row.pmra ± row.pmra_error,
		pmdec=row.pmdec ± row.pmdec_error, 
		radial_velocity = row.RV ± row.RV_err
	) for row in eachrow(averaged_stars)]

# ╔═╡ d0842a73-dd1b-452f-97f9-50a23668fe14
methods(lguys.ICRS)

# ╔═╡ 11461563-dbcd-48fc-a97a-1a0391538462
gsr = lguys.transform.(lguys.GSR{Measurement{Float64}}, icrs)

# ╔═╡ 718ec8bc-cae6-4905-a336-b04964699b61
vra = [o.pmra for o in gsr] .* distance

# ╔═╡ 379bff0e-e65c-4abd-bce1-d6b776656bc8
vdec = [o.pmdec for o in gsr] .* distance

# ╔═╡ 9991361e-0c93-41d0-810c-b501f7ed1201
rv_meas = averaged_stars

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

# ╔═╡ 1971afe2-6fd2-4b3b-a228-4264ab34fc51
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = L"\xi",
		ylabel = L"\eta",
		aspect=DataAspect()
	)


	filt = memb_filt
	dv = Measurements.value.(vz) .- Measurements.value.(rv)
	p =scatter!(averaged_stars.xi[filt], averaged_stars.eta[filt], color=dv[filt], colorrange=0.3 .* (-1, 1), markersize= 5, colormap=:RdBu)

	Colorbar(fig[1,2], p, label=L"dV / km s$^{-1}$")
	fig
end

# ╔═╡ 0eecd302-185b-4467-8614-d25fa26a9b5d
let
	fig = Figure()
	ax = Axis(fig[1, 1], 
		xlabel=L"angular distance (degrees)",
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
		ylabel = L"$v_x - v_\textrm{los}$ / km\,s$^{-1}$",
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

# ╔═╡ 94f371f5-bb3a-4400-b7bc-6cc76cb43927
begin
	averaged_stars_w_vels = copy(averaged_stars)
	averaged_stars_w_vels[!, :VX] = Measurements.value.(vx)
	averaged_stars_w_vels[!, :VY] = Measurements.value.(vy)
	averaged_stars_w_vels[!, :VZ] = Measurements.value.(vz)
	averaged_stars_w_vels[!, :VX_err] = Measurements.uncertainty.(vx)
	averaged_stars_w_vels[!, :VY_err] = Measurements.uncertainty.(vy)
	averaged_stars_w_vels[!, :VZ_err] = Measurements.uncertainty.(vz)

	averaged_stars_w_vels[!, :RV_gsr] = Measurements.value.(rv)
	averaged_stars_w_vels[!, :pmra_gsr] = [Measurements.value.(o.pmra) for o in gsr]
	averaged_stars_w_vels[!, :pmdec_gsr] = [Measurements.value.(o.pmdec) for o in gsr]
end

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

	filt = filt_missing(rv1; verbose=true, low=rv_low, high=rv_high) .& filt_missing(rv2;  verbose=true, low=rv_low, high=rv_high)

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
	
	hlines!(0, color=:black)
	return fig
end

# ╔═╡ aaaf5ba2-c9ed-41ec-a22a-d78ed96fd84e
compare_rv_mean("p20", memb_stars)

# ╔═╡ a2370516-e7ec-4502-ae4b-b111bcf68d36
compare_rv_mean("apogee", memb_stars)

# ╔═╡ 78e7aff0-3658-4d64-b117-58189a85307a
compare_rv_mean("s18", memb_stars)

# ╔═╡ f4f9dd06-1a1a-458b-be75-05d52623580c
compare_rv_mean("s24")

# ╔═╡ 74152829-27ef-4d8d-8b32-ed30a18f30e4
rv_meas.RV[not.(ismissing.(rv_meas.RV_s24))]

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

# ╔═╡ 95ef12d2-174a-47c3-a106-9b4005b2a80d
md"""
# Full fits to data
"""

# ╔═╡ 9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
function fit_rv_sigma(rv, rv_err; N=3_000, p=0.16, burn=0.2)
	μ = median(rv)
	μ_p = NaN
	μ_err = std(rv) / sqrt(length(rv))

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
	x_model = LinRange(rv_low, rv_high, 1000)
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

# ╔═╡ 5bee3041-fe7f-4f67-a73d-fb60ed006959
md"""
## For each dataset
"""

# ╔═╡ 2de62c9d-f76b-4e3b-9ad2-69d81842e507
begin 
	pace20_filt = pace20.PdSph .> 0.2
end

# ╔═╡ 88e5fea7-d42c-4429-873f-93ddad293b19
let 
	global pace20_cleaned = average_by_star(pace20[pace20_filt, :], id_col="ID", study_rv_names=[])

	filt_good = pace20_cleaned.RV_std .< 5pace20_cleaned.RV_err  .* sqrt.(pace20_cleaned.RV_count)

	println("num bad: ", sum(.!filt_good))

	pace20_cleaned = pace20_cleaned[filt_good, :]
end

# ╔═╡ 721080b2-e57b-4353-84a3-41f1fb740d92
pace20_cleaned

# ╔═╡ 6942bdca-1300-46be-8aa2-a3725d7fdfba
apogee_memb_filt = (apogee.PSAT .> 0.2 ) .&& (apogee.RV .> rv_low) .&& (apogee.RV .< rv_high) .&& (apogee.RV_sigma ./ apogee.RV_err .< 5)

# ╔═╡ 133a324f-9a29-4f4e-a8b0-8a712e4f2493
apogee_cleaned = average_by_star(apogee[apogee_memb_filt, :], id_col="source_id", study_rv_names=[])

# ╔═╡ 904ff35b-dbfa-4ab8-ae9a-59f87abf6c9d
let 
	global spencer_cleaned = average_by_star(spencer18_all, id_col="ID", study_rv_names=[])

	filt_good = spencer_cleaned.RV_std ./ spencer_cleaned.RV_err .< 5 .* sqrt.(spencer_cleaned.RV_count)

	println(sum(.!filt_good))
	spencer_cleaned = spencer_cleaned[filt_good, :]

end
	

# ╔═╡ 1901e3d8-dcb5-4bb2-8f13-7472b1c5b6e4
all_studies

# ╔═╡ ebdf8cfc-ebff-4331-bab2-998a734f2cd1
datasets = Dict(
	:apogee => apogee_cleaned,
	:p20 => pace20_cleaned,
	:s18 => spencer_cleaned,
	:sestito => graces,
)

# ╔═╡ 781a4656-76e1-4c6e-9057-611b90bb5649
gaia_memb_filt = (rv_gaia.PSAT .> 0.2 ) .&& (rv_gaia.RV .> rv_low) .&& (rv_gaia.RV .< rv_high)

# ╔═╡ 2bae2760-76f1-4930-85cd-8001a8f6de7f
sum(gaia_memb_filt)

# ╔═╡ 2a1bde56-9cb3-42b8-909b-a7078da574e1
apogee_memb_filt

# ╔═╡ ec4bc80f-22e6-49f9-a589-5c5bc9e50a8b
props = Dict(name => fit_rv_sigma(float.(data.RV), float.(data.RV_err)) for (name, data) in datasets)

# ╔═╡ eaf47845-95dd-4f7d-bf76-3d9b1154711a
for key in keys(props)
	println(key)
	println(props[key])
	fig =  plot_sample_normal_fit(datasets[key], props[key], title=string(key))
	if key !== :sestito
		rv = averaged_stars[:, "RV_$key"]
		filt = filt_missing(rv) .& memb_filt
		rv_err = averaged_stars[filt, "RV_err_$key"]
		rv = rv[filt]

		println(sum(filt))
		h = histogram(Float64.(rv), normalization=:pdf)
		
		errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=:red)
		println(fit_rv_sigma(rv, rv_err))
		println()
	end

	@info fig
end

# ╔═╡ b5faa873-d5c9-447b-b90e-f693db26e6c2
md"""
### Comparing sample properties
"""

# ╔═╡ 8e7c650c-ebc5-4775-a725-0810b735cc36
0.3 ⊕ 0.5

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
	ks = keys(props) |> collect
	ks = ks[ks .!= :sestito]
	y = [props[k][2] for k in ks]
	yerr = [props[k][4] for k in ks]

	x = collect(1:length(ks))
	
	fig, ax = FigAxis(
		xticks=(x, string.(ks)),
		xminorticksvisible=false,
		ylabel=L"$\sigma_{v}$ / km s$^{-1}$",
	)


	errorscatter!(x, y, yerror=yerr)
	fig
end

# ╔═╡ ce3067d3-55e6-43a1-9b7b-4cf53c09ee88
md"""
# J24 purity checks
"""

# ╔═╡ 01e0cfcd-92af-4c04-a7d4-9c9d8fd5a9f1
rv_mean = -246.9; sigma_los = 8 # from mcmc modeling in next section

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

# ╔═╡ 7110364b-67ea-4e2a-b59f-a3dc9b292757
sum(isnan.(averaged_stars_w_vels.RV_gsr))

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

# ╔═╡ 43b5cf29-6645-4758-b73c-cfe77367da07
logit(x) = log(x / (1-x))

# ╔═╡ 678a13d6-c892-43fe-a806-f3534661f785
let
	fig, ax = FigAxis(
		xlabel = "logit P SAT",
		ylabel = "purity",
		limits=(-10, 15, nothing, nothing)
	)

	ps = quantile(rv_meas.PSAT[rv_meas.PSAT .>= 0], LinRange(0, 1, 30))

	y = purity_given_p.(ps[1:end-1], ps[2:end])

	scatterlines!(logit.(midpoints(ps)), y)
	vlines!(logit(0.2), linestyle=:dot, color=:black)
	fig
end

# ╔═╡ b95daa04-4f93-4920-ae96-2368b6fd82bf
logit(1-0.001)

# ╔═╡ 5bef0a45-2a4c-49bc-889f-b8bea315bebe
tanh(3)

# ╔═╡ 7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
md"""
## Binned properties
"""

# ╔═╡ 33f54afc-cdb9-4eb8-887f-5a43281b837c
let
	fig = Figure(figsize=(900, 400))
	ax = Axis(fig[1,1],
		xlabel = "R / arcmin",
		ylabel = L"RV / km s$^{-1}$",
		#limits=(nothing, (60, 150))
	)

	scatter!(memb_stars.r_ell, memb_stars.RV_apogee, label="APOGEE", markersize=5)
	scatter!(memb_stars.r_ell, memb_stars.RV_s24, label="sestito+23", markersize=10)
	scatter!(memb_stars.r_ell, memb_stars.RV_p20, label="pace+20", markersize=5)
	scatter!(memb_stars.r_ell, memb_stars.RV_s18, label="spencer + 18", markersize=3)

	Legend(fig[1,2], ax)


	resize_to_layout!(fig)
	fig
end

# ╔═╡ 500815e1-f9e1-4401-ba2d-72f326cfa783
let
	fig = Figure(figsize=(900, 400))
	ax = Axis(fig[1,1],
		xlabel = "R / arcmin",
		ylabel = L"RV / km s$^{-1}$",
		#limits=(nothing, (rv_low, rv_high))
	)

	rv_meas = averaged_stars
	#scatter!(memb_stars.r_ell, memb_stars.RV_dart, label="DART", markersize=5)
	scatter!(rv_meas.r_ell, rv_meas.RV_apogee, label="APOGEE", markersize=5)
	#scatter!(rv_meas.r_ell, rv_meas.RV_gmos, label="GMOS", markersize=10)
	scatter!(rv_meas.r_ell, rv_meas.RV_p20, label="Walker+09", markersize=5)
	scatter!(rv_meas.r_ell, rv_meas.RV_s18, label="tolstoy + 23", markersize=3)

	Legend(fig[1,2], ax)


	resize_to_layout!(fig)
	fig
end

# ╔═╡ 725b9c7d-93f9-48c3-b97b-8a1147a16f78


# ╔═╡ f6a9738c-94cc-449d-92ed-805c5579d579
function replace_missings(df)
	out_df = copy(df)
	for col in names(out_df)
		val = out_df[:, col]
		if eltype(val) <: Union{Missing, <:AbstractFloat}
			out_df[:, col] = replace(out_df[:, col], missing=>NaN)
		elseif eltype(val) <: Union{Missing, String}
			out_df[:, col] = replace(out_df[:, col], missing=>"")
		elseif eltype(val) <: Union{Missing, Integer}
			out_df[:, col] = replace(out_df[:, col], missing=>0)
		end
		
	end
	disallowmissing!(out_df)

	return out_df
end

# ╔═╡ fd3996c2-c6e9-4245-8fee-95097ef8f6bf
begin 
	df_averaged = replace_missings(averaged_stars_w_vels)
	df_averaged.RV = convert.(Float64, df_averaged.RV)
	df_averaged.RV_err = convert.(Float64, df_averaged.RV_err)

	select!(df_averaged, Not(["RV_$study" for study in all_studies]))
	select!(df_averaged, Not(["RV_err_$study" for study in all_studies]))
	df_averaged
end

# ╔═╡ e5d9eb9e-33c2-49dc-acf1-05a9e54c53b7
begin 
	df_all = replace_missings(all_stars)
	df_all.RV = convert.(Float64, df_all.RV)
	df_all.RV_err = convert.(Float64, df_all.RV_err)

	df_all
end

# ╔═╡ 79344fee-f3ca-41fc-b650-872d56e8cd0f
for col in names(df_averaged)
	println(col, "\t\t", typeof(df_averaged[:, col]))
end

# ╔═╡ a44cb867-b1c9-4380-b871-74a6b3eb8804
df_all.RV

# ╔═╡ 1a0eac4e-4100-4da6-820b-19c34e283118
lguys.write_fits(joinpath("processed", "umi_all_rv.fits"), df_all, overwrite=true)

# ╔═╡ 46b539eb-c195-4ca4-922d-63298afaebbc
lguys.write_fits(joinpath("processed", "umi_averaged_rv.fits"), df_averaged, overwrite=true)

# ╔═╡ Cell order:
# ╟─811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═12e0e7ad-925b-47dc-a1d8-c56e65ffaa4c
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
# ╠═ef1bb6f5-00b8-405b-8702-244eca618644
# ╠═923674f9-a5fa-454d-a5e2-a4d30d047a11
# ╠═68c13192-58c0-4e77-9f86-567d55a036b5
# ╠═17629b76-2cb8-4a56-9687-3a630950872b
# ╠═d4123ffd-cb32-486e-a93d-f48d7112a831
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═da5a3b57-72bc-46e1-b1a0-6c02eb101626
# ╠═3d8b3c55-153b-4a4a-9289-78f34df07abc
# ╠═d415ba9e-75c0-4192-8761-dbda70e43a73
# ╠═77e7884c-0360-4b7f-b9ad-e81be2516552
# ╠═7a50a176-96a5-4098-88d6-0fa2874d0f90
# ╟─c470c2b9-093d-42ab-b96d-9dac231ccabc
# ╠═bb7f6769-ec92-460d-8423-449029175f79
# ╠═4e8b7258-76fa-404f-91ca-8cfda241a332
# ╠═4ebb0411-1264-4297-be0c-cacfdb2fd59d
# ╠═14f866c4-8c65-461f-9834-329b281a3c0a
# ╠═c941ef83-42b3-4bdf-b1db-6f0b9b3f1cd0
# ╠═5f71b8e9-5540-4431-8028-4ce14c8d7856
# ╠═4d63430e-f59c-4c68-97e1-7eaa5679e55f
# ╠═bf088da8-a7c5-46e7-8e5f-ef5b91c10fb8
# ╠═1bc2541c-69ac-40aa-94f5-743707ca5a66
# ╟─27063de7-01d4-48a8-a06f-cc24aec662d2
# ╠═b6de4afb-9362-4618-a114-df460031e4f9
# ╠═b7345279-4f80-47ad-a726-537571849eae
# ╟─c31f9e07-c520-443b-94dc-787519021d01
# ╠═15f2a8e2-90df-48a9-a7bf-e86955f566ce
# ╠═c50824f9-652f-47e5-9a05-fa052c36a247
# ╠═7922960e-a02a-4c2b-b7f0-0aca1d514626
# ╠═c405ede3-258c-4b3d-8c4e-e4bab68f3bf0
# ╠═b48ea14d-48d7-4953-8d6a-b84bce1a4788
# ╠═180fac98-678c-4a14-966c-385387c60ac3
# ╠═19d8aa00-14d1-4270-8e8c-9b8907779e3b
# ╠═1e6ff00b-eee8-4e03-8b54-e27d6cbb6f16
# ╠═3056ca82-9ffe-4065-b207-ad958aaab08b
# ╠═0313553f-7a05-492b-99b9-9648e98baad4
# ╠═dd136efd-ed53-4472-ac73-2cd97a4cdad2
# ╠═6f1ff560-c458-48cd-8d7f-9f3dd7631b08
# ╠═26092ed5-6283-4bf9-8ee1-60a8cf15e1ef
# ╠═1c6074eb-bafa-439c-aa55-ee7a5c3a58d5
# ╠═cdf27cbc-fd37-412f-a926-e1f86f836bb8
# ╠═721080b2-e57b-4353-84a3-41f1fb740d92
# ╠═36cdb5ff-65a3-4c5c-b274-4d7da53d2ab4
# ╠═1ca9b968-87d0-46cc-a1fe-634717ac7446
# ╠═fe6fa910-f41e-4657-836b-7eda2f0cddb2
# ╠═842e2c79-5eee-4bf3-b86d-eda1204336d7
# ╠═51468818-3e4e-4445-b23f-4df757c80478
# ╠═ac200cf0-714c-4b66-855a-e204735c7c04
# ╠═bc4fd1bb-fb30-4581-8952-891293476345
# ╠═c235cd75-cada-4d8f-a380-43af619eaa24
# ╠═a85d8d3f-197a-4f7f-97cd-5e2523a885d7
# ╠═b94b7591-018b-47d7-8af8-be0c4b3ea49f
# ╠═e472cbb6-258e-4303-85e2-56f26358c97b
# ╠═0c8bceff-8ec3-471e-ac1b-e0d6caecfa66
# ╠═28bf5ad2-26f0-4a25-890f-2f8c80f9e6e6
# ╠═eec5c292-b7d0-49cf-83dc-283c400b9c05
# ╠═3e37a9dd-fa57-45f8-999e-e499b760a420
# ╠═14d30389-c25c-426c-99f5-68fe98aabc7f
# ╠═85511c49-479f-44a0-b88d-b33ffef99e0a
# ╠═29e0bd8d-d66f-44d5-81af-b0f54ec060e6
# ╠═51b3371a-87ef-4c10-a312-d5fa00ff7844
# ╠═2eba0d13-c688-4a70-95e4-59bc0d21020b
# ╠═bd51fe42-7e39-4cd8-8065-58ab9814f966
# ╟─89552b84-d12e-4d88-a58e-8b89ad4b2569
# ╠═ff236dd1-3097-44f6-8e12-00803a8bee42
# ╠═e6f2de3b-ce32-4d61-851f-4e42fcce95c0
# ╠═7a77bf54-fff8-4906-b738-0d3c91bdf6aa
# ╠═5b2fceff-9c3e-472d-9310-31e920137e41
# ╠═4b305b83-1a3b-48a6-b19f-6f3ebed0768f
# ╠═6ed9fabe-1bdd-4ecf-90b4-084d5296150e
# ╠═f7213298-3dcb-49ac-a9f0-a129a03423aa
# ╠═a4fd8597-8706-45aa-a8c2-f72d90cfbbfd
# ╠═f7ec8bba-9f45-435b-b67c-33182e992dfd
# ╠═93d185f2-1eaa-4d35-87dd-b84f385483de
# ╠═36d2d86f-6a75-46f4-b48f-36137e95e90d
# ╠═66e5693f-90cc-4b49-8a0c-8d92ea9e8cdd
# ╠═ce62265b-e411-4e0a-942c-281f82c03c6e
# ╠═ee4d33fd-1bd8-41e8-8dc3-17092e4c267d
# ╠═f3ab19a6-95c6-4756-b22c-c7abfcf330c1
# ╠═a49845c5-0b66-4316-bf25-febb57b77a51
# ╠═02921108-b708-4a90-b7f9-6903632da601
# ╠═9213cc36-74d1-452f-bd9a-eb5c5cab1f87
# ╠═db06c81e-1456-4e12-a95a-c87bb23eaec3
# ╠═f9e1bcbd-2def-4e64-b69e-d1cfcf7ffedd
# ╠═f61debe4-8b23-4415-b77c-43c4465ccfcb
# ╠═3ad3bd41-56f4-4fa7-8345-c275303c637d
# ╠═102c73ef-2c95-4784-80df-ed0312511c00
# ╠═3a9841cc-4e19-4962-8b5f-dd25b223189e
# ╠═cd256d12-5c9e-436b-9473-ed3e7e713225
# ╠═dfb4db7f-1d27-4e0b-9eeb-0afc2ca3d92a
# ╠═578b9efe-9cd5-4b55-a321-4eea37bc43d1
# ╠═f4505ad1-6026-4f07-b849-6f55401f2d80
# ╠═8bc140fa-6f4e-4ec5-bc95-989fc0ea48d1
# ╠═77130b64-383e-4b3c-bb55-76f925db9986
# ╠═db9a60bf-a7db-410d-ada2-0b6ed011e2c4
# ╠═782691d3-dd01-4a34-9db2-a18b0b6566e1
# ╠═0d0c854c-b794-469d-8dba-2bf83a9f6d2f
# ╠═36f74c41-4190-4705-a656-61a1223ffed8
# ╠═2365da83-7be8-44ec-a399-d520d610014b
# ╠═828cf7c4-d91e-4d02-8649-dd0194177f67
# ╠═2078b2ab-e973-40fc-b19c-83ae860b626d
# ╠═1f40dc2a-0a93-4b41-a550-1eb5b77fc037
# ╠═07070037-e64a-46eb-a23a-4cc7e35b4551
# ╠═0008053d-77ec-4b74-a701-1e584254ed11
# ╠═7b9a5d6d-48b0-4426-8145-19d5a414af8d
# ╠═d867234f-dfc2-436c-b83d-02e87089897a
# ╠═29f47769-8d2c-4e9d-b5eb-79b168580a04
# ╠═3655088f-8ae3-4edc-899d-466d909e5e6c
# ╠═d3333b48-aa4e-42c1-9e0a-bbff98e3647d
# ╠═88ed48b4-baa9-4ac0-86e1-8348edcd59b4
# ╠═222bb254-8b65-44d3-b3d2-b08fbcbe9950
# ╠═6bc02c4f-31a2-4e4e-8612-e66f8cc9c93e
# ╠═11fcf4f8-fcd5-4426-a4e9-b046138bde1b
# ╠═877706f1-08a0-496a-890d-622e3a2fd9ec
# ╠═733fe42e-b7a5-4285-8c73-9a41e4488d40
# ╠═74b10a3e-1342-454f-8eed-77b371f81edf
# ╠═d8800a31-1ed3-422f-ac51-90f18cf61c29
# ╠═015f2345-74d4-4296-b0dc-35d6e13ad8cc
# ╠═cc6c65db-ef57-4745-8ada-e11427274a77
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╠═a1938588-ca40-4844-ab82-88c4254c435b
# ╠═1823f690-81c9-471b-bdec-c9fe261544c1
# ╠═8a5615a0-c34c-4de5-8a62-543cdb613988
# ╠═8fb2fe7f-9f74-4cf6-9cef-4bac3ebc9ae4
# ╠═a7071b79-f36d-453a-9f9e-2a16832f0785
# ╠═c96b5112-7f07-4e7c-b42b-18b26055561a
# ╠═2be7ce8f-f44b-4e57-b8a6-53b623376c4a
# ╠═4ed42682-b81c-4635-b5a9-7a73857cc927
# ╟─ab27964f-b070-4328-81c2-6ecf4b8cec1e
# ╠═23d4b7c4-61f9-4748-b3c4-80eaad551914
# ╠═7dab615b-e0a1-431e-88c0-e1e9d9c29568
# ╠═c0accd82-17a7-437c-9ba8-4db437071a5b
# ╠═d0842a73-dd1b-452f-97f9-50a23668fe14
# ╠═11461563-dbcd-48fc-a97a-1a0391538462
# ╠═718ec8bc-cae6-4905-a336-b04964699b61
# ╠═379bff0e-e65c-4abd-bce1-d6b776656bc8
# ╠═9991361e-0c93-41d0-810c-b501f7ed1201
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
# ╠═1971afe2-6fd2-4b3b-a228-4264ab34fc51
# ╠═0eecd302-185b-4467-8614-d25fa26a9b5d
# ╠═e0111a22-8dd9-40d2-81e9-fcffe04adf4e
# ╠═79ffc43f-be70-40c7-9419-8666d33b5947
# ╠═94f371f5-bb3a-4400-b7bc-6cc76cb43927
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
# ╠═be10c541-246f-4fcf-a5e5-416ecb69d7e7
# ╟─5bee3041-fe7f-4f67-a73d-fb60ed006959
# ╠═2de62c9d-f76b-4e3b-9ad2-69d81842e507
# ╠═88e5fea7-d42c-4429-873f-93ddad293b19
# ╠═6942bdca-1300-46be-8aa2-a3725d7fdfba
# ╠═133a324f-9a29-4f4e-a8b0-8a712e4f2493
# ╠═904ff35b-dbfa-4ab8-ae9a-59f87abf6c9d
# ╠═1901e3d8-dcb5-4bb2-8f13-7472b1c5b6e4
# ╠═ebdf8cfc-ebff-4331-bab2-998a734f2cd1
# ╠═781a4656-76e1-4c6e-9057-611b90bb5649
# ╠═2bae2760-76f1-4930-85cd-8001a8f6de7f
# ╠═2a1bde56-9cb3-42b8-909b-a7078da574e1
# ╠═ec4bc80f-22e6-49f9-a589-5c5bc9e50a8b
# ╠═eaf47845-95dd-4f7d-bf76-3d9b1154711a
# ╟─b5faa873-d5c9-447b-b90e-f693db26e6c2
# ╠═8e7c650c-ebc5-4775-a725-0810b735cc36
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
# ╠═7110364b-67ea-4e2a-b59f-a3dc9b292757
# ╠═e899faa9-580b-4aad-902e-382008048908
# ╠═c0a4e207-9778-4d43-8197-5fc8887b2b82
# ╠═69161135-1f19-4d9a-8ff6-ae63a79e3eb5
# ╠═ab12510f-9393-4acd-8ed0-dcffa06c65e6
# ╠═678a13d6-c892-43fe-a806-f3534661f785
# ╠═43b5cf29-6645-4758-b73c-cfe77367da07
# ╠═b95daa04-4f93-4920-ae96-2368b6fd82bf
# ╠═5bef0a45-2a4c-49bc-889f-b8bea315bebe
# ╟─7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═33f54afc-cdb9-4eb8-887f-5a43281b837c
# ╠═500815e1-f9e1-4401-ba2d-72f326cfa783
# ╠═725b9c7d-93f9-48c3-b97b-8a1147a16f78
# ╠═f6a9738c-94cc-449d-92ed-805c5579d579
# ╠═fd3996c2-c6e9-4245-8fee-95097ef8f6bf
# ╠═e5d9eb9e-33c2-49dc-acf1-05a9e54c53b7
# ╠═79344fee-f3ca-41fc-b650-872d56e8cd0f
# ╠═a44cb867-b1c9-4380-b871-74a6b3eb8804
# ╠═1a0eac4e-4100-4da6-820b-19c34e283118
# ╠═46b539eb-c195-4ca4-922d-63298afaebbc
