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

	using PyFITS
	import LilGuys as lguys
end

# ╔═╡ 7ed5bcf5-dfc3-4e79-a608-d503124a1e96
using LilGuys

# ╔═╡ 02c55b13-5681-475b-bead-9e0e0b9e9656
using Measurements

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Radial Velocity CrossMatch

This notebook loads in radial velocity data preprocessed and matched to gaia from different sources and combines these into a signle catalogue.
"""

# ╔═╡ 9485165e-0cb3-452d-a0ec-953b790c9d7f
FIGDIR = "figures"

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, kurtosis, sem

# ╔═╡ 257caa92-ccce-44ff-88ee-6a9c550ae5d9
CairoMakie.activate!(type=:png)

# ╔═╡ ef1bb6f5-00b8-405b-8702-244eca618644
import DensityEstimators: histogram, calc_limits, make_bins

# ╔═╡ 58caad92-a887-4527-ac84-be02fba5b23c
import TOML

# ╔═╡ 4a0d2b79-be8e-460e-9cd2-f873d1af237a
module GaiaFilters
	include(joinpath(ENV["DWARFS_ROOT"], "utils", "gaia_filters.jl"))
end

# ╔═╡ 9a59505b-039b-41e1-8907-a8107ce68177
module RVUtils
	include("../../rv_utils.jl")
end

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "../data/"

# ╔═╡ 300aa043-3750-4d01-ae48-1bf765cd92b1
md"""
we can load  stars without fbest, but this causes havoc later on since the probabilities are not defined.
"""

# ╔═╡ 7a50a176-96a5-4098-88d6-0fa2874d0f90
j24 = read_fits("../processed/best_sample.fits")

# ╔═╡ 77e7884c-0360-4b7f-b9ad-e81be2516552
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 6964faa7-17d2-412c-b4a2-5f981f8c4b54
md"""
we reproduce the mean, but slightly lower velocity dispersion, maybe due to different selection criteria?
"""

# ╔═╡ d4d0a488-1c0e-4bc6-9a88-94b27d84e8ce
apogee_all = read_fits("processed/rv_apogee.fits")

# ╔═╡ 9ed91d41-e650-4d3c-9e74-100ff3d57d82
walker09_all = let
	df = read_fits("processed/rv_walker+09.fits")

	df = df[.!ismissing.(df.source_id), :]
end

# ╔═╡ 15f2a8e2-90df-48a9-a7bf-e86955f566ce
tolstoy23_all = read_fits("processed/rv_tolstoy+23.fits")

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
	
		filt, idx = RVUtils.xmatch(gmos, j24)
		println(filt)
		gmos[!, :source_id] = j24.source_id[idx]

	end

	gmos[!, :RV_count] .= 1
	gmos[!, :RV_sigma] .= 0.
	
	gmos
end

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
all_stars.RV_gmos[.!ismissing.(all_stars.RV_gmos)]

# ╔═╡ ff644334-bf90-49e4-a3da-fc4ecc840929
length(unique(apogee_all.source_id))

# ╔═╡ 8f243944-3d0b-48ce-9f90-8d448c089239
md"""
double check we have all stars
"""

# ╔═╡ bd51fe42-7e39-4cd8-8065-58ab9814f966
 sum(.!ismissing.(all_stars.RV_apogee)), length(apogee_all.RV)

# ╔═╡ 0c6c6d84-fb14-4eaf-858e-22fb90ab4d8d
@assert length(unique(walker09_all.source_id)) == length(walker09_all.source_id)

# ╔═╡ ab49efb3-2ab7-47d0-a3a5-a342c789aa9b
length(walker09_all.source_id),  sum(.!ismissing.(all_stars.RV_w09))

# ╔═╡ de762a39-b430-4452-ba87-8b8cf1ad9852
length(tolstoy23_all.RV), sum(.!ismissing.(all_stars.RV_t23))

# ╔═╡ 04ff1abd-c584-41de-8c83-7503482c3731
md"""
Gaia stars are too bright to be members :(
"""

# ╔═╡ 0f50ae12-e207-4748-9c45-780c9215be73
sum(skipmissing(.!ismissing.(all_stars.RV_err_gaia) .* (all_stars.PSAT .> 001)))

# ╔═╡ 89552b84-d12e-4d88-a58e-8b89ad4b2569
md"""
# Validation for xmatch
"""

# ╔═╡ e6f2de3b-ce32-4d61-851f-4e42fcce95c0
function plot_xmatch_radec(suffix)

	ra2 = all_stars[:, "ra_$suffix"]
	filt = .!ismissing.(ra2)
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
	filt =  @. !ismissing(col) & !isnan(col)
	filt1 = high .> col .> low
	if verbose
		println("excluding ", sum(.!(filt1)[filt]), " outliers")
	end
	return filt .& filt1
end

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

	text!(0.1, 0.8, space=:relative, text="μ = $(μ ± δμ)\nσ = $(round(sigma, digits=2))", fontsize=10)

	
	lines!([75, 150], [75, 150], color=:black)
	return fig
end

# ╔═╡ f9e1bcbd-2def-4e64-b69e-d1cfcf7ffedd
@savefig "t23_w09_xmatch" compare_rv("t23", "w09")

# ╔═╡ 609f55cd-aa49-4a5e-b871-413ec7ef0990
import StatsBase as sb

# ╔═╡ f61debe4-8b23-4415-b77c-43c4465ccfcb
@savefig "t23_apogee_xmatch" compare_rv("t23", "apogee")

# ╔═╡ 102c73ef-2c95-4784-80df-ed0312511c00
@savefig "apogee_w09_xmatch" compare_rv("apogee", "w09")

# ╔═╡ a488eeb2-bf24-408e-a2ba-7df1b305a2b8
md"""
# Combining observations
"""

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

	x_mean = LilGuys.mean(x, w)
	std_err = sqrt(1/sum(w))
	sigma = LilGuys.std(x, w)
	counts_tot = sum(n)
	
	return x_mean, std_err, sigma, counts_tot 
end

# ╔═╡ aee8f66c-ff93-496e-a21d-48faefd75a43
function sigma_interstudy(values, errors, sigmas, counts)
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
	return sigma 
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

# ╔═╡ 877706f1-08a0-496a-890d-622e3a2fd9ec
RV_a_std = [safe_weighted_mean(rvs[i, :], rv_errs[i, :], rv_sigmas[i, :], rv_counts[i, :]) for i in 1:size(rvs, 1)]

# ╔═╡ 73bd1553-e2ae-4bfb-aac1-0880346f5054
begin
	filt_is_meas = .!ismissing.(RV_a_std)
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

# ╔═╡ 0cf94e27-b200-4516-9495-dcac7f10d9a0
sum(.!ismissing.(rv_meas.RV_apogee) .& ismissing.(rv_meas.RV_t23) .& ismissing.(rv_meas.RV_w09)) # most apogee stars are covered by other surveys.

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

# ╔═╡ 8c1457d6-32d4-4727-b45e-1a46244708a6
RV_err

# ╔═╡ 63e28d21-c1ab-413e-ac7e-e6776bcf6ca5
rvs[filt_is_meas, :]

# ╔═╡ dbb51473-6506-4119-a11d-50df2ceb6dcf
rv_errs[filt_is_meas, :]

# ╔═╡ d913bf9b-c03b-4bce-a39d-7d23eca9b355
rv_sigmas[filt_is_meas, :]

# ╔═╡ e0905cda-92ab-4073-9632-ce6857ecb999
size(rv_sigmas)

# ╔═╡ 24b2455e-1c9a-4b9a-9d88-19b98d4b2165
RV_sigma

# ╔═╡ 70e75283-b6f2-4c37-95ea-1c0b2e3f2110
N_studies = sum.(.!ismissing.(eachrow(rvs)))[filt_is_meas]

# ╔═╡ 01b3135f-7a72-4669-b586-4bc5894464ad
sum(filt_is_meas)

# ╔═╡ d11edca7-b9ae-4269-9e1b-661d59bd965e
all_stars[.!ismissing.(all_stars.RV_gmos), :].source_id

# ╔═╡ 22564a47-9b03-4778-b30c-d092581ec107
md"""
# Quality flag
"""

# ╔═╡ 4cda9764-f208-4532-b51f-5deb62992467
md"""
We want to remove stars with statistically large standard deviations. So if the standard deviation of a given study is larger than 5 times the statistical uncertainty, or our interstudy std is > 5 times sqrt(n) times the standard error, we have a problem
"""

# ╔═╡ 9aafd6b9-6ec8-4b6c-908e-1faac4615e0b
sigma_sigma = ifelse.(isfinite.(rv_meas.RV_sigma), 
	rv_meas.RV_sigma ./ (rv_meas.RV_err .* sqrt.(N_studies)),
	0.				
	)

# ╔═╡ c099b81d-f987-47f7-85db-bdf883ac62f0
rv_meas.RV_sigma

# ╔═╡ d665895c-9a97-4c36-9ea9-c6e30d0a17ad
mean(sigma_sigma .> 5)

# ╔═╡ 15b1baf7-d9bc-4d3e-a9fa-9d8b8a4dbc6e
F_qual_inter = sigma_sigma .< 5

# ╔═╡ c833ebb2-9f5e-4fb4-b812-dbc7fe57d8b5
function filter_qual_study(rv_meas, study)
	map(eachrow(rv_meas)) do row
		if ismissing(row["RV_sigma_$study"])
			return true
		elseif isnan(row["RV_sigma_$study"])
			return true
		else
			row["RV_sigma_$study"] .< 5 * row["RV_err_$study"] * sqrt(row["RV_count_$study"])
		end

	end
end

# ╔═╡ dcec54f8-27b4-405e-a386-6e07470dc073
F_qual_study = filter_qual_study(rv_meas, "w09") .& filter_qual_study(rv_meas, "apogee") .& filter_qual_study(rv_meas, "t23")

# ╔═╡ 73a41773-39bc-4bd5-9d99-39804000631a
F_qual = F_qual_study .& F_qual_inter

# ╔═╡ 5c86f548-4838-41e1-b9af-fdd93d900940
sum(.!F_qual)

# ╔═╡ 870571f3-f201-47fd-a215-995a416bc223
sum(.!F_qual_study)

# ╔═╡ 06b5c8d8-e531-4f00-a3cf-8d6202bb71f2
let
	fig = Figure()
	ax = Axis(fig[1,1],
		yscale=log10, yticks=Makie.automatic,
		xlabel = "log sigma / std",
		ylabel = "count",
			  limits=(-3, nothing, nothing, nothing)
			 )

	
	hist!(filter(isfinite, log10.(sigma_sigma)), bins=120)
	vlines!(log10(5), color=COLORS[2])

	fig

end

# ╔═╡ db3177e7-7132-4f62-846a-f4416a804009
md"""
# write data
"""

# ╔═╡ 7faa2813-e502-4187-855a-047a2f5dd48d
σ_pm = lguys.kms2pm(obs_properties["sigma_v"], obs_properties["distance"])

# ╔═╡ 21daa18e-88c3-43ee-8dfa-ea8d116e5ff4
sum(.!ismissing.(rv_meas.RV_gmos))

# ╔═╡ 7567ad1b-3f31-49cc-8efa-ad5d9e0c3bf2
icrs1 = RVUtils.icrs(obs_properties)

# ╔═╡ 6a406ddc-3fe9-40ff-adcf-9e2e429016a7
md"""
The velocity correction is magnitude 3 kms with uncertainty < 5x measurement uncertainties 
"""

# ╔═╡ 520f0ae4-eae4-4ab1-86e6-5758a2d7ebad


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
		#aspect=DataAspect()
	)

	x = float.(rv2[filt])
	y = float.(rv1[filt] .- rv2[filt])
	xerr = float.(rv2_err[filt])
	yerr = float.(rv1_err[filt])
	errorscatter!(x, y, xerror=xerr, yerror=yerr, color=(COLORS[1], 0.2), size=1)

	w = 1 ./ xerr .^2
	μ = LilGuys.mean(y, w)
	δμ = sqrt( 1 / sum(w))
	s = LilGuys.std(y, w)
	text!(0.1, 0.1, space=:relative, text="μ = $(μ ± δμ)\nσ = $(round(s, digits=2))", fontsize=10)

	hlines!(0, color=:black)
	return fig
end

# ╔═╡ f4f9dd06-1a1a-458b-be75-05d52623580c
compare_rv_mean("gmos", )

# ╔═╡ 74152829-27ef-4d8d-8b32-ed30a18f30e4
rv_meas.RV_gmos[.!ismissing.(rv_meas.RV_gmos)]

# ╔═╡ 5688e9c9-515c-4973-820f-1215031253f2
import StatsBase: var

# ╔═╡ ce3067d3-55e6-43a1-9b7b-4cf53c09ee88
md"""
# J24 purity checks
"""

# ╔═╡ 01e0cfcd-92af-4c04-a7d4-9c9d8fd5a9f1
rv_mean = obs_properties["radial_velocity"]; sigma_los = obs_properties["sigma_v"] # from mcmc modeling in next section

# ╔═╡ f15ad626-bb0d-4484-8004-d6f81ccf0825
⊕ = GaiaFilters.:⊕

# ╔═╡ a6bc7223-3267-4760-9b6a-d886ac6f4544
obs_properties["pmra_err"] ⊕ σ_pm

# ╔═╡ 0e098399-85ae-4cdf-9d82-654a9fd1cd35
icrs0 = lguys.ICRS(
	ra=obs_properties["ra"],
	dec=obs_properties["dec"],
	distance = obs_properties["distance"] ± obs_properties["distance_err"],
	pmra = obs_properties["pmra"] ± (obs_properties["pmra_err"] ⊕ σ_pm),
	pmdec = obs_properties["pmdec"] ± (obs_properties["pmdec_err"] ⊕ σ_pm),
	radial_velocity = obs_properties["radial_velocity"] ± obs_properties["sigma_v"]
			   )

# ╔═╡ 2771f91a-a132-4fb1-ac93-36525690a3eb
icrs1.radial_velocity, icrs0.radial_velocity

# ╔═╡ 42fa6c2c-c9cd-4b0d-b1be-701a2db318ee
gsr0 = lguys.transform(lguys.GSR, icrs0)

# ╔═╡ d7c091aa-a0d5-4c78-bd50-5e3c521fa1ed
sind(1) * lguys.pm2kms(gsr0.pmdec, gsr0.distance)

# ╔═╡ 2b6a78d1-e69e-4661-85de-1f73c2ab1dcb
gsr0.pmdec

# ╔═╡ 2a07cc9f-e5a8-4baf-af1c-e69f1def7bc0
gsr0.pmdec

# ╔═╡ 6ea1eff9-14d4-459c-b086-199e0aaf84f0
Δrv = RVUtils.rv_correction(rv_meas.ra, rv_meas.dec, Measurements.value(icrs0.ra), Measurements.value(icrs0.dec), 
						   gsr0.pmra, gsr0.pmdec, gsr0.distance)

# ╔═╡ 2ecbc309-2b17-4e75-a09b-6483b8b0a711
hist(Measurements.uncertainty.(Δrv) ./ rv_meas.RV_err )

# ╔═╡ 615f05fb-4907-40bc-9251-3065f565929b
let
	filename = "processed/rv_combined.fits"
	if isfile(filename)
		rm(filename); 
	end
	rv_meas[!, :F_RV] = F_qual
	rv_meas[!, :delta_rv] = Measurements.value.(Δrv)
	rv_meas[!, :delta_rv_err] = Measurements.uncertainty.(Δrv)
	write_fits(filename, rv_meas)
	@info "wrote data"
end

# ╔═╡ 4d3cd232-81c0-4c32-bfc9-7ba8098d2e4d
icrs0.ra, icrs0.dec, icrs0.pmra, icrs0.pmdec

# ╔═╡ 18adfc17-ee1f-4cfb-b194-cd8494850200
δ_gsr_0 = gsr0.radial_velocity .- icrs0.radial_velocity

# ╔═╡ ca3f437e-6c55-462f-b85d-68b19a241e8a
sort(rv_meas.xi .⊕ rv_meas.eta) ./ 60

# ╔═╡ 03462abb-32d6-433c-a9c4-65ef2dd58965
function is_rv_member(rv, err)
	return abs(rv - rv_mean) / (err ⊕ sigma_los) < 3
end

# ╔═╡ 82d4e049-7f47-4f8b-a3a0-7191ef82912b
filt_rv = is_rv_member.(rv_meas.RV, rv_meas.RV_err) 

# ╔═╡ 3920e27f-48fd-41a5-bf53-1f80edf3d56d
memb_stars = rv_meas[rv_meas.PSAT .> 0.2 .&& F_qual .&& filt_rv .|| .!ismissing.(rv_meas.RV_gmos), :]

# ╔═╡ a2370516-e7ec-4502-ae4b-b111bcf68d36
compare_rv_mean("apogee", memb_stars)

# ╔═╡ aaaf5ba2-c9ed-41ec-a22a-d78ed96fd84e
compare_rv_mean("w09", memb_stars)

# ╔═╡ 78e7aff0-3658-4d64-b117-58189a85307a
compare_rv_mean("t23", memb_stars)

# ╔═╡ 92bfc52c-05f5-45e3-ae26-90c579ecfe2c
hist(filter((!) ∘ ismissing, memb_stars.RV_err_t23 ./ memb_stars.RV_err))

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

# ╔═╡ e343dcd0-078f-476d-abe0-fd84001b72ae
mean(filt_rv[.!isnan.(rv_meas.RV_err)])

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
		ylabel = "purity"
	)

	ps = 0:0.1:1

	y = purity_given_p.(ps[1:end-1], ps[2:end])

	lines!(midpoints(ps), y)

	hidexdecorations!(ax, ticks=false, minorticks=false)

	ax2 = Axis(fig[2, 1], ylabel = "# in bin", xlabel="PSAT",
			  yscale=log10, yticks=Makie.automatic,
			  )
	ylims!(1, nothing)

	y = number_satisfying.(ps[1:end-1], ps[2:end])

	lines!(midpoints(ps), y)

	rowgap!(fig.layout, 0)
	rowsize!(fig.layout, 2, Relative(1/3))
	
	@savefig "purity_vs_psat"
end

# ╔═╡ 7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
md"""
## Quick properties
"""

# ╔═╡ 30e3dc5b-3ce6-4dd7-9c2a-c82774909a8c
sum(.!ismissing.(memb_stars.RV_gmos))

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

	@savefig "rv_scatter_alstudy"
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

# ╔═╡ 62301344-869a-4ec0-8299-29f0ff5d1c15
sum(sum(eachcol(ismissing.(rv_meas))))

# ╔═╡ 1deb1520-194a-40b5-8968-bf73b756ba3d
rv_meas[ .! ismissing.(rv_meas.RV_gmos), :PSAT]

# ╔═╡ Cell order:
# ╟─811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═9485165e-0cb3-452d-a0ec-953b790c9d7f
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═257caa92-ccce-44ff-88ee-6a9c550ae5d9
# ╠═ef1bb6f5-00b8-405b-8702-244eca618644
# ╠═7ed5bcf5-dfc3-4e79-a608-d503124a1e96
# ╠═58caad92-a887-4527-ac84-be02fba5b23c
# ╠═4a0d2b79-be8e-460e-9cd2-f873d1af237a
# ╠═9a59505b-039b-41e1-8907-a8107ce68177
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═300aa043-3750-4d01-ae48-1bf765cd92b1
# ╠═7a50a176-96a5-4098-88d6-0fa2874d0f90
# ╠═77e7884c-0360-4b7f-b9ad-e81be2516552
# ╟─6964faa7-17d2-412c-b4a2-5f981f8c4b54
# ╠═d4d0a488-1c0e-4bc6-9a88-94b27d84e8ce
# ╠═9ed91d41-e650-4d3c-9e74-100ff3d57d82
# ╠═15f2a8e2-90df-48a9-a7bf-e86955f566ce
# ╠═8e0095e7-2198-439d-bd19-976bdb1d3766
# ╠═b7345279-4f80-47ad-a726-537571849eae
# ╟─bbf49122-11b1-4272-a660-0437c6aa2b3f
# ╠═fe6fa910-f41e-4657-836b-7eda2f0cddb2
# ╠═e472cbb6-258e-4303-85e2-56f26358c97b
# ╠═0d2dbc73-1ded-46e3-b142-9bc7b777728d
# ╠═ff644334-bf90-49e4-a3da-fc4ecc840929
# ╠═8f243944-3d0b-48ce-9f90-8d448c089239
# ╠═bd51fe42-7e39-4cd8-8065-58ab9814f966
# ╠═0c6c6d84-fb14-4eaf-858e-22fb90ab4d8d
# ╠═ab49efb3-2ab7-47d0-a3a5-a342c789aa9b
# ╠═de762a39-b430-4452-ba87-8b8cf1ad9852
# ╠═04ff1abd-c584-41de-8c83-7503482c3731
# ╠═0f50ae12-e207-4748-9c45-780c9215be73
# ╟─89552b84-d12e-4d88-a58e-8b89ad4b2569
# ╠═e6f2de3b-ce32-4d61-851f-4e42fcce95c0
# ╠═5b2fceff-9c3e-472d-9310-31e920137e41
# ╠═4b305b83-1a3b-48a6-b19f-6f3ebed0768f
# ╠═c218cfa8-2f55-4957-bcdd-8b3970fe639a
# ╠═f7213298-3dcb-49ac-a9f0-a129a03423aa
# ╟─f7ec8bba-9f45-435b-b67c-33182e992dfd
# ╠═93d185f2-1eaa-4d35-87dd-b84f385483de
# ╠═9213cc36-74d1-452f-bd9a-eb5c5cab1f87
# ╠═f9e1bcbd-2def-4e64-b69e-d1cfcf7ffedd
# ╠═609f55cd-aa49-4a5e-b871-413ec7ef0990
# ╠═f61debe4-8b23-4415-b77c-43c4465ccfcb
# ╠═0cf94e27-b200-4516-9495-dcac7f10d9a0
# ╠═102c73ef-2c95-4784-80df-ed0312511c00
# ╟─a488eeb2-bf24-408e-a2ba-7df1b305a2b8
# ╠═d3333b48-aa4e-42c1-9e0a-bbff98e3647d
# ╠═8bc140fa-6f4e-4ec5-bc95-989fc0ea48d1
# ╠═aee8f66c-ff93-496e-a21d-48faefd75a43
# ╠═d3ef6a98-8770-4359-9e37-94f08d2d9aa5
# ╠═78d39e37-44b5-4171-ab45-0672f607b517
# ╠═88ed48b4-baa9-4ac0-86e1-8348edcd59b4
# ╠═877706f1-08a0-496a-890d-622e3a2fd9ec
# ╠═73bd1553-e2ae-4bfb-aac1-0880346f5054
# ╠═222bb254-8b65-44d3-b3d2-b08fbcbe9950
# ╠═e3f05ee2-cc5f-437e-801d-3c7d842af709
# ╠═a13a33f1-65c4-4217-8636-639b1e14d109
# ╠═ee3c22db-6b6b-4c30-8d1f-86b103c018fc
# ╠═6bc02c4f-31a2-4e4e-8612-e66f8cc9c93e
# ╠═11fcf4f8-fcd5-4426-a4e9-b046138bde1b
# ╠═8c1457d6-32d4-4727-b45e-1a46244708a6
# ╠═63e28d21-c1ab-413e-ac7e-e6776bcf6ca5
# ╠═dbb51473-6506-4119-a11d-50df2ceb6dcf
# ╠═d913bf9b-c03b-4bce-a39d-7d23eca9b355
# ╠═e0905cda-92ab-4073-9632-ce6857ecb999
# ╠═24b2455e-1c9a-4b9a-9d88-19b98d4b2165
# ╠═70e75283-b6f2-4c37-95ea-1c0b2e3f2110
# ╠═01b3135f-7a72-4669-b586-4bc5894464ad
# ╠═d11edca7-b9ae-4269-9e1b-661d59bd965e
# ╠═22564a47-9b03-4778-b30c-d092581ec107
# ╟─4cda9764-f208-4532-b51f-5deb62992467
# ╠═9aafd6b9-6ec8-4b6c-908e-1faac4615e0b
# ╠═c099b81d-f987-47f7-85db-bdf883ac62f0
# ╠═d665895c-9a97-4c36-9ea9-c6e30d0a17ad
# ╠═15b1baf7-d9bc-4d3e-a9fa-9d8b8a4dbc6e
# ╠═73a41773-39bc-4bd5-9d99-39804000631a
# ╠═dcec54f8-27b4-405e-a386-6e07470dc073
# ╠═870571f3-f201-47fd-a215-995a416bc223
# ╠═5c86f548-4838-41e1-b9af-fdd93d900940
# ╠═c833ebb2-9f5e-4fb4-b812-dbc7fe57d8b5
# ╠═06b5c8d8-e531-4f00-a3cf-8d6202bb71f2
# ╟─db3177e7-7132-4f62-846a-f4416a804009
# ╠═7faa2813-e502-4187-855a-047a2f5dd48d
# ╠═a6bc7223-3267-4760-9b6a-d886ac6f4544
# ╠═02c55b13-5681-475b-bead-9e0e0b9e9656
# ╠═21daa18e-88c3-43ee-8dfa-ea8d116e5ff4
# ╠═0e098399-85ae-4cdf-9d82-654a9fd1cd35
# ╠═7567ad1b-3f31-49cc-8efa-ad5d9e0c3bf2
# ╠═2771f91a-a132-4fb1-ac93-36525690a3eb
# ╠═42fa6c2c-c9cd-4b0d-b1be-701a2db318ee
# ╠═6ea1eff9-14d4-459c-b086-199e0aaf84f0
# ╠═d7c091aa-a0d5-4c78-bd50-5e3c521fa1ed
# ╠═4d3cd232-81c0-4c32-bfc9-7ba8098d2e4d
# ╠═2b6a78d1-e69e-4661-85de-1f73c2ab1dcb
# ╠═6a406ddc-3fe9-40ff-adcf-9e2e429016a7
# ╠═ca3f437e-6c55-462f-b85d-68b19a241e8a
# ╠═2a07cc9f-e5a8-4baf-af1c-e69f1def7bc0
# ╠═18adfc17-ee1f-4cfb-b194-cd8494850200
# ╠═2ecbc309-2b17-4e75-a09b-6483b8b0a711
# ╠═520f0ae4-eae4-4ab1-86e6-5758a2d7ebad
# ╠═615f05fb-4907-40bc-9251-3065f565929b
# ╟─3655b6a0-ac9a-4a49-86fd-6a6484409819
# ╠═7f13339e-6a0c-4944-8a36-5ab136fd8415
# ╠═3920e27f-48fd-41a5-bf53-1f80edf3d56d
# ╠═a2370516-e7ec-4502-ae4b-b111bcf68d36
# ╠═aaaf5ba2-c9ed-41ec-a22a-d78ed96fd84e
# ╠═78e7aff0-3658-4d64-b117-58189a85307a
# ╠═92bfc52c-05f5-45e3-ae26-90c579ecfe2c
# ╠═f4f9dd06-1a1a-458b-be75-05d52623580c
# ╠═74152829-27ef-4d8d-8b32-ed30a18f30e4
# ╠═764b5306-20f9-4810-8188-1bdf9482260f
# ╠═5688e9c9-515c-4973-820f-1215031253f2
# ╟─ce3067d3-55e6-43a1-9b7b-4cf53c09ee88
# ╠═01e0cfcd-92af-4c04-a7d4-9c9d8fd5a9f1
# ╠═f15ad626-bb0d-4484-8004-d6f81ccf0825
# ╠═03462abb-32d6-433c-a9c4-65ef2dd58965
# ╠═82d4e049-7f47-4f8b-a3a0-7191ef82912b
# ╠═e343dcd0-078f-476d-abe0-fd84001b72ae
# ╠═0105a7ec-392e-426b-9bf9-a4b6144fe829
# ╠═7b22fa97-c311-42ae-ad5a-41d053c5334e
# ╠═e899faa9-580b-4aad-902e-382008048908
# ╠═c0a4e207-9778-4d43-8197-5fc8887b2b82
# ╠═69161135-1f19-4d9a-8ff6-ae63a79e3eb5
# ╠═ab12510f-9393-4acd-8ed0-dcffa06c65e6
# ╠═678a13d6-c892-43fe-a806-f3534661f785
# ╟─7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═30e3dc5b-3ce6-4dd7-9c2a-c82774909a8c
# ╠═33f54afc-cdb9-4eb8-887f-5a43281b837c
# ╠═dcb1bcf1-41e3-4389-b879-be6a53b5e670
# ╠═a70a2546-3aaa-4d78-a060-2a9a39940dc8
# ╠═382e7503-0e31-49d5-b593-a6a1a0f9c5e0
# ╠═0b7107f8-38d2-4b0c-aad1-e07d3927cf87
# ╠═e32f0d5a-9a6a-407c-a466-586c5d63fda8
# ╠═500815e1-f9e1-4401-ba2d-72f326cfa783
# ╠═62301344-869a-4ec0-8299-29f0ff5d1c15
# ╠═1deb1520-194a-40b5-8968-bf73b756ba3d
