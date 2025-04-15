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

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Radial Velocity CrossMatch

This notebook loads in radial velocity data preprocessed and matched to gaia from different sources and combines these into a signle catalogue.
"""

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, kurtosis, sem

# ╔═╡ 257caa92-ccce-44ff-88ee-6a9c550ae5d9
CairoMakie.activate!(type=:png)

# ╔═╡ ef1bb6f5-00b8-405b-8702-244eca618644
import DensityEstimators: histogram, calc_limits, make_bins

# ╔═╡ 58caad92-a887-4527-ac84-be02fba5b23c
import TOML

# ╔═╡ e3b19a39-fafb-4038-94da-f33400cf2a32
module RVUtils
	include("rv_utils.jl")
	not = !
end

# ╔═╡ 4a0d2b79-be8e-460e-9cd2-f873d1af237a
module GaiaFilters
	include(joinpath(ENV["DWARFS_ROOT"], "utils", "gaia_filters.jl"))
end

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "data/"

# ╔═╡ 300aa043-3750-4d01-ae48-1bf765cd92b1
md"""
we can load  stars without fbest, but this causes havoc later on since the probabilities are not defined.
"""

# ╔═╡ 7a50a176-96a5-4098-88d6-0fa2874d0f90
j24 = read_fits("processed/best_sample.fits")

# ╔═╡ 77e7884c-0360-4b7f-b9ad-e81be2516552
obs_properties = TOML.parsefile("observed_properties.toml")

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
@assert sum(.!ismissing.(all_stars.RV_apogee)) == length(apogee_all.RV)

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
	@info "mean difference = $μ ± $δμ"
	@info "std = $sigma"
	@info "mean uncertainty = $(mean(1 ./ sqrt.(w)))"
	lines!([75, 150], [75, 150], color=:black)
	return fig
end

# ╔═╡ f9e1bcbd-2def-4e64-b69e-d1cfcf7ffedd
compare_rv("t23", "w09")

# ╔═╡ f61debe4-8b23-4415-b77c-43c4465ccfcb
compare_rv("t23", "apogee")

# ╔═╡ 102c73ef-2c95-4784-80df-ed0312511c00
compare_rv("apogee", "w09")

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

# ╔═╡ 04323c34-f8dc-460b-add9-a5fcc7f9dd1b
md"""
# Likelihoods
"""

# ╔═╡ 6863d4cf-b19c-4383-9b97-ce9174199c46
μ_v = obs_properties["radial_velocity"]

# ╔═╡ 0326f962-f222-417f-b40b-eb5d847216b4
σ_v = obs_properties["sigma_v"]

# ╔═╡ c3ceb65f-bc86-4451-802c-b2e398f5bf11
import Distributions: pdf, Normal

# ╔═╡ 3ec085ba-a139-4bac-b3e1-fb699dc5bfd8
function L_rv(star)
	return pdf(Normal(μ_v, sqrt(σ_v^2 + star.RV_err^2)), star.RV)
end

# ╔═╡ 2934a2c6-2ff0-4b35-866c-4647389a24e6
function L_rv(RV, RV_err)
	return pdf(Normal(μ_v, sqrt(σ_v^2 + RV_err^2)), RV)
end

# ╔═╡ 685952a1-8a81-4618-ad2a-736e02da3d0a
sum(1 .- rv_meas.PSAT), sum(rv_meas.PSAT)

# ╔═╡ eebb35c0-0d29-488d-8f5b-c55dca624305
maximum(rv_meas.RV) - minimum(rv_meas.RV)

# ╔═╡ 181295fc-dd18-479b-a2b0-311cfc9fe53b
L_rv.(eachrow(rv_meas))

# ╔═╡ a79fe1ac-ba76-4ada-bce4-8d04750b51ec
rv_meas[!, :L_RV_SAT] = L_rv.(eachrow(rv_meas))

# ╔═╡ c401d944-10f0-483e-8837-c682a95578ff
LLR_rv = rv_meas.L_RV_SAT ./ rv_meas.L_RV_BKD

# ╔═╡ 7cb9c085-4cd1-438a-888c-c5ddb8e7ae22
sum(LLR_rv .< 1)

# ╔═╡ 44eec58d-1314-4664-aca7-c6c557036069
md"""
Find solar motion corrected RV
"""

# ╔═╡ d1d81751-7148-4186-b3c0-bf5c5a888610
gsr_static = LilGuys.GSR(ra=obs_properties["ra"], dec=obs_properties["dec"], distance = 100, pmra=0, pmdec=0, radial_velocity=0)

# ╔═╡ 07a758ab-7b0d-4f85-8778-cdc420367cdb
rv_offset = LilGuys.to_icrs(gsr_static).radial_velocity

# ╔═╡ c85af471-0923-4310-b9a6-6f42c9739eca
L_bg(rv) = pdf(Normal(rv_offset, 100), rv)

# ╔═╡ 452e261b-47d5-4799-9f3c-1ddbc99b6469
rv_meas[!, :L_RV_BKD] .= L_bg.(rv_meas.RV)

# ╔═╡ 9693cd53-26df-4b8b-9b17-94611da3e3c8
let 
	fig = Figure()
		
	ax = Axis(fig[1,1], xlabel = "radial velocity")
		
	
	stephist!(collect(skipmissing(rv_meas[(.!ismissing.(rv_meas.PSAT)) .& (rv_meas.PSAT .< 0.05), :RV])), normalization=:pdf)

	lines!(Normal(rv_offset, 100))

	fig
end

# ╔═╡ e948ea6f-fdba-4f39-a5f5-67eafd50b0a8
hist(asinh.(rv_meas.L_RV_SAT ./ rv_meas.L_RV_BKD ./ 10))

# ╔═╡ 27c22ea0-6359-43f6-9cf5-9ddf07ab6b77
f_sat = mean(skipmissing(all_stars.PSAT))

# ╔═╡ 4d8ecd12-1b4c-402b-9011-80d77fd9a527
L_SAT = @. rv_meas.L_CMD_SAT * rv_meas.L_PM_SAT * rv_meas.L_S_SAT * rv_meas.L_RV_SAT

# ╔═╡ b9e05adf-fa65-4a2d-b56b-24150f447abf
L_BKD = @. rv_meas.L_CMD_BKD * rv_meas.L_PM_BKD * rv_meas.L_S_BKD * rv_meas.L_RV_BKD

# ╔═╡ 878b633a-3f95-4b7c-bd60-439d1c28cfaa
rv_meas[!, :PSAT_RV] = @. f_sat * L_SAT / (f_sat * L_SAT + (1-f_sat) * L_BKD)

# ╔═╡ ae30a25e-6742-4650-ab0d-04ad07b4c5cd
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "PSAT + RV",
		ylabel = "PSAT",
		aspect = DataAspect()
			 )

	lines!([0,1], [0,1], color=:black)

	x = filter(isfinite, rv_meas.PSAT_RV)
	y = collect(skipmissing(rv_meas.PSAT))
	scatter!(x, y, alpha=0.3, markersize=3, color=rv_meas.RV)

	fig

end

# ╔═╡ e3584246-a99b-43f7-9302-97a28fb59da0
rv_meas.PSAT_RV

# ╔═╡ 2ec8db9f-f4d7-45b5-9a15-99cfced78765
hist(filter(x->isfinite(x) & !ismissing(x), rv_meas.PSAT_RV .- rv_meas.PSAT))

# ╔═╡ 5f17795b-ba5c-44d1-93ce-d57ba90550bd
hist(filter(isfinite, rv_meas.PSAT_RV), axis=(yscale=log10, yticks=Makie.automatic))

# ╔═╡ 35164ff4-5a83-4144-b74d-e024af17dc47
rv_meas[.!ismissing.(rv_meas.RV_gmos), [:source_id, :ra, :dec, :pmra, :pmdec, :LLR_S, :LLR_PM, :LLR_CMD, :PSAT, :PSAT_RV]]

# ╔═╡ ab27964f-b070-4328-81c2-6ecf4b8cec1e
md"""
# Corrected coordinate frame
"""

# ╔═╡ 57a77527-a14c-47a2-8949-1cbdafc49994
import Measurements: ±

# ╔═╡ 5930b133-3e3c-43be-9ed6-f125ea4418d4
import Measurements

# ╔═╡ 23d4b7c4-61f9-4748-b3c4-80eaad551914
distance = obs_properties["distance"] ± float.(obs_properties["distance_err"])

# ╔═╡ 6d6259fc-f8e6-4d65-b7a5-ae14ee1d5139
LilGuys.arcmin2kpc(30, obs_properties["distance"])

# ╔═╡ 7dab615b-e0a1-431e-88c0-e1e9d9c29568
ra0, dec0 = obs_properties["ra"], obs_properties["dec"]

# ╔═╡ 0e080c9f-d334-4dff-9041-4c83ba65506e
σ_pm = LilGuys.kms2pm(obs_properties["sigma_v"], obs_properties["distance"])

# ╔═╡ 46f6c9c6-20b3-49d7-8b8a-98eb2f1ffd4a
obs_properties["pmra_err"]

# ╔═╡ c4859522-5f50-4042-aa24-6bf9b214cb22
pmra = map(eachrow(rv_meas)) do row
	w1 = 1/row.pmra^2
	w2 = 1/(obs_properties["pmra_err"]^2 + σ_pm^2) * row.PSAT_RV
	m = (w1*row.pmra + w2 * obs_properties["pmra"] ) / (w1+w2)
	err = 1/sqrt(w1 + w2)
	return m ± err
end

# ╔═╡ ad7ae1a8-c5a4-4d79-bd86-35431fa71856
pmdec = map(eachrow(rv_meas)) do row
	w1 = 1/row.pmdec^2
	w2 = 1/(obs_properties["pmdec_err"]^2 + σ_pm^2) * row.PSAT_RV
	m = (w1*row.pmdec + w2 * obs_properties["pmdec"] ) / (w1+w2)
	err = 1/sqrt(w1 + w2)
	return m ± err

end

# ╔═╡ 61265c8d-5d16-497e-871a-72ee49b8b362
errorscatter(rv_meas.pmra, Measurements.value.(pmra), 
			 xerror=rv_meas.pmra, yerror=Measurements.uncertainty.(pmra) ,
			 color = (:black, 0.2),
			 axis = (;
				 limits = (-1, 1, -1, 1)
			 )
			)

# ╔═╡ f45540ff-658c-4ab2-bfa9-62ce3913cec7
hist(filter(x->isfinite(x) & !ismissing(x),
	(rv_meas.pmra .- Measurements.value.(pmra) ./ rv_meas.pmra_error)
		   ), bins=100)

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
vra = [lguys.pm2kms(o.pmra, distance) for o in gsr] 

# ╔═╡ 379bff0e-e65c-4abd-bce1-d6b776656bc8
vdec = [lguys.pm2kms(o.pmdec, distance) for o in gsr] 

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
		limits = (40, 100, -2, 2)
	)

	errorscatter!(Measurements.value.(rv), Measurements.value.(vz) .- Measurements.value.(rv), 
		#xerr=Measurements.uncertainty.(rv), yerr=Measurements.uncertainty.(vz),
		alpha=0.03, 
		color=:black
	)
	fig
end

# ╔═╡ 53d128cf-4af1-4cac-9bdc-6fed1c903524
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = L"\xi",
		ylabel = L"\eta",
		xreversed=true,
		aspect=DataAspect(),
	)
			  
			  
	p = scatter!(rv_meas.xi, rv_meas.eta, color = Measurements.value.(vz) .- Measurements.value.(rv), colorrange=(-1, 1), colormap=:redsblues, markersize=2,
	   )

	Colorbar(fig[1,2], p, label = "vz - rv")

	fig
end

# ╔═╡ 842d0ab5-7de7-424f-9da3-1acfa9282a4c
v_shift_outliers = rv_meas[abs.( Measurements.value.(vz) .- Measurements.value.(rv)) .> 1, :]

# ╔═╡ db5ccd85-e17a-4630-b5a1-b21df1ced791
sort(v_shift_outliers[:,  [:source_id, :xi, :eta, :pmra, :pmra_error, :pmdec, :pmdec_error, :LLR_PM, :LLR_CMD, :LLR_S, :RV, :RV_err, :PSAT, :PSAT_RV] ], :PSAT_RV)

# ╔═╡ 1acc759c-f28e-41ce-a26b-7f4301a45a65
md"""
Only star with large velocity shift and nonzero probability is 
- `5026130884816022016`
but this shift is only about 1 km/s and is because this is federico's star
"""

# ╔═╡ f5ac233f-eba0-49a9-8005-272e43446f95
filt_special_star = rv_meas.source_id .== 5026130884816022016

# ╔═╡ 675bb06f-723a-405d-bc77-5aa20dc4ea4d
rv[filt_special_star]

# ╔═╡ 7a20dd82-b545-4981-95b7-7d360bc9695b
vz[filt_special_star]

# ╔═╡ 0eecd302-185b-4467-8614-d25fa26a9b5d
let
	fig = Figure()
	ax = Axis(fig[1, 1], 
		xlabel="angular distance (degrees)",
		ylabel = L"$v_x - v_\textrm{los}$ / km\,s$^{-1}$",
		xgridvisible=false,
		ygridvisible=false,
			  limits=(nothing, nothing, -1, 1)
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
		limits = (40, 110, 0, 0.01)
		
	)

	errorscatter!(Measurements.value.(rv), log10.(Measurements.uncertainty.(vz) 
 ./ Measurements.uncertainty.(rv)), 
		alpha=0.03, 
		color=:black
	)
	fig
end

# ╔═╡ fad90f18-6797-4daa-9aed-c9d2c882a0a6
md"""
the increase in uncertainty by adding the projection is negligable.
"""

# ╔═╡ 750cf383-7004-4160-bbda-7c304690a0c7
rv_higher_uncertainty = rv_meas[Measurements.uncertainty.(vz)  ./ Measurements.uncertainty.(rv) .> 1.1, :]

# ╔═╡ dcd445bb-1c90-4e55-ab89-9a9697c3a360
rv_higher_uncertainty[:, [:xi, :eta, :LLR_S, :LLR_CMD, :LLR_PM, :RV, :RV_err, :PSAT_RV]]

# ╔═╡ 79ffc43f-be70-40c7-9419-8666d33b5947
3*9 * sind(2)

# ╔═╡ aa1edde8-2b95-4e9b-b844-ba3150e81a0a
hist(rv_meas.RV_count, axis=(;yscale=log10, yticks=Makie.automatic))

# ╔═╡ 9cd82f65-ad8e-4e32-a9f9-c8096be2becd
mean(rv_meas.RV_count)

# ╔═╡ 3b80b8ee-22b8-4fc7-b2f9-67274b494c03
  !ismissing(missing) && missing

# ╔═╡ 4adf3ac6-b11c-4391-8202-e5a67765e192
sqrt(1 / (1/0.65^2 + 1/1.0^2))

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
		ylabel = "count"
			 )

	
	hist!(filter(isfinite, log10.(sigma_sigma)), bins=120)
	vlines!(log10(5), color=COLORS[2])

	fig

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

# ╔═╡ 3920e27f-48fd-41a5-bf53-1f80edf3d56d
memb_stars = rv_meas[rv_meas.PSAT_RV .> 0.2 .&& F_qual, :]

# ╔═╡ a2370516-e7ec-4502-ae4b-b111bcf68d36
compare_rv_mean("apogee", memb_stars)

# ╔═╡ aaaf5ba2-c9ed-41ec-a22a-d78ed96fd84e
compare_rv_mean("w09", memb_stars)

# ╔═╡ 78e7aff0-3658-4d64-b117-58189a85307a
compare_rv_mean("t23", memb_stars)

# ╔═╡ 92bfc52c-05f5-45e3-ae26-90c579ecfe2c
hist(filter((!) ∘ ismissing, memb_stars.RV_err_t23 ./ memb_stars.RV_err))

# ╔═╡ f4f9dd06-1a1a-458b-be75-05d52623580c
compare_rv_mean("gmos", )

# ╔═╡ 74152829-27ef-4d8d-8b32-ed30a18f30e4
rv_meas.RV_gmos[.!ismissing.(rv_meas.RV_gmos)]

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

# ╔═╡ 5688e9c9-515c-4973-820f-1215031253f2
import StatsBase: var

# ╔═╡ ce3067d3-55e6-43a1-9b7b-4cf53c09ee88
md"""
# J24 purity checks
"""

# ╔═╡ 01e0cfcd-92af-4c04-a7d4-9c9d8fd5a9f1
rv_mean = obs_properties["radial_velocity"]; sigma_los = obs_properties["sigma_v"] # from mcmc modeling in next section

# ╔═╡ f246ba9f-6cf2-4793-9612-54f2a8475e30
⊕(a, b) = sqrt(a^2 + b^2)

# ╔═╡ 03462abb-32d6-433c-a9c4-65ef2dd58965
function is_rv_member(rv, err)
	return abs(rv - rv_mean) / (err ⊕ sigma_los) < 3
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
		xlabel = "PSAT",
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

	out_df[!, :vx] = Measurements.value.(vx)
	out_df[!, :vy] = Measurements.value.(vy)
	out_df[!, :vz] = Measurements.value.(vz)
	out_df[!, :vx_err] = Measurements.uncertainty.(vx)
	out_df[!, :vy_err] = Measurements.uncertainty.(vy)
	out_df[!, :vz_err] = Measurements.uncertainty.(vz)
	out_df[!, :F_RV] = F_qual

	
	
	select!(out_df, Not([:o_Target_w09, ]))

	out_df
end

# ╔═╡ 62301344-869a-4ec0-8299-29f0ff5d1c15
sum(sum(eachcol(ismissing.(out_df))))

# ╔═╡ 1deb1520-194a-40b5-8968-bf73b756ba3d
rv_meas[ .! ismissing.(rv_meas.RV_gmos), :PSAT_RV]

# ╔═╡ 615f05fb-4907-40bc-9251-3065f565929b
let
	filename = "processed/sculptor_all_rv.fits"
	if isfile(filename)
		rm(filename); 
	end
	write_fits(filename, out_df)
	@info "wrote data"
end

# ╔═╡ Cell order:
# ╟─811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═257caa92-ccce-44ff-88ee-6a9c550ae5d9
# ╠═ef1bb6f5-00b8-405b-8702-244eca618644
# ╠═7ed5bcf5-dfc3-4e79-a608-d503124a1e96
# ╠═58caad92-a887-4527-ac84-be02fba5b23c
# ╠═e3b19a39-fafb-4038-94da-f33400cf2a32
# ╠═4a0d2b79-be8e-460e-9cd2-f873d1af237a
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
# ╠═f61debe4-8b23-4415-b77c-43c4465ccfcb
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
# ╟─04323c34-f8dc-460b-add9-a5fcc7f9dd1b
# ╠═6863d4cf-b19c-4383-9b97-ce9174199c46
# ╠═0326f962-f222-417f-b40b-eb5d847216b4
# ╠═c3ceb65f-bc86-4451-802c-b2e398f5bf11
# ╠═3ec085ba-a139-4bac-b3e1-fb699dc5bfd8
# ╠═2934a2c6-2ff0-4b35-866c-4647389a24e6
# ╠═c85af471-0923-4310-b9a6-6f42c9739eca
# ╠═685952a1-8a81-4618-ad2a-736e02da3d0a
# ╠═eebb35c0-0d29-488d-8f5b-c55dca624305
# ╠═181295fc-dd18-479b-a2b0-311cfc9fe53b
# ╠═a79fe1ac-ba76-4ada-bce4-8d04750b51ec
# ╠═452e261b-47d5-4799-9f3c-1ddbc99b6469
# ╠═c401d944-10f0-483e-8837-c682a95578ff
# ╠═7cb9c085-4cd1-438a-888c-c5ddb8e7ae22
# ╠═9693cd53-26df-4b8b-9b17-94611da3e3c8
# ╠═44eec58d-1314-4664-aca7-c6c557036069
# ╠═d1d81751-7148-4186-b3c0-bf5c5a888610
# ╠═07a758ab-7b0d-4f85-8778-cdc420367cdb
# ╠═e948ea6f-fdba-4f39-a5f5-67eafd50b0a8
# ╠═27c22ea0-6359-43f6-9cf5-9ddf07ab6b77
# ╠═4d8ecd12-1b4c-402b-9011-80d77fd9a527
# ╠═b9e05adf-fa65-4a2d-b56b-24150f447abf
# ╠═878b633a-3f95-4b7c-bd60-439d1c28cfaa
# ╠═ae30a25e-6742-4650-ab0d-04ad07b4c5cd
# ╠═e3584246-a99b-43f7-9302-97a28fb59da0
# ╠═2ec8db9f-f4d7-45b5-9a15-99cfced78765
# ╠═5f17795b-ba5c-44d1-93ce-d57ba90550bd
# ╠═35164ff4-5a83-4144-b74d-e024af17dc47
# ╟─ab27964f-b070-4328-81c2-6ecf4b8cec1e
# ╠═57a77527-a14c-47a2-8949-1cbdafc49994
# ╠═5930b133-3e3c-43be-9ed6-f125ea4418d4
# ╠═23d4b7c4-61f9-4748-b3c4-80eaad551914
# ╠═6d6259fc-f8e6-4d65-b7a5-ae14ee1d5139
# ╠═7dab615b-e0a1-431e-88c0-e1e9d9c29568
# ╠═0e080c9f-d334-4dff-9041-4c83ba65506e
# ╠═46f6c9c6-20b3-49d7-8b8a-98eb2f1ffd4a
# ╠═c4859522-5f50-4042-aa24-6bf9b214cb22
# ╠═ad7ae1a8-c5a4-4d79-bd86-35431fa71856
# ╠═61265c8d-5d16-497e-871a-72ee49b8b362
# ╠═f45540ff-658c-4ab2-bfa9-62ce3913cec7
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
# ╠═53d128cf-4af1-4cac-9bdc-6fed1c903524
# ╠═842d0ab5-7de7-424f-9da3-1acfa9282a4c
# ╠═db5ccd85-e17a-4630-b5a1-b21df1ced791
# ╠═1acc759c-f28e-41ce-a26b-7f4301a45a65
# ╠═f5ac233f-eba0-49a9-8005-272e43446f95
# ╠═675bb06f-723a-405d-bc77-5aa20dc4ea4d
# ╠═7a20dd82-b545-4981-95b7-7d360bc9695b
# ╠═0eecd302-185b-4467-8614-d25fa26a9b5d
# ╠═e0111a22-8dd9-40d2-81e9-fcffe04adf4e
# ╠═fad90f18-6797-4daa-9aed-c9d2c882a0a6
# ╠═750cf383-7004-4160-bbda-7c304690a0c7
# ╠═dcd445bb-1c90-4e55-ab89-9a9697c3a360
# ╠═79ffc43f-be70-40c7-9419-8666d33b5947
# ╠═aa1edde8-2b95-4e9b-b844-ba3150e81a0a
# ╠═9cd82f65-ad8e-4e32-a9f9-c8096be2becd
# ╠═3b80b8ee-22b8-4fc7-b2f9-67274b494c03
# ╠═4adf3ac6-b11c-4391-8202-e5a67765e192
# ╠═22564a47-9b03-4778-b30c-d092581ec107
# ╠═4cda9764-f208-4532-b51f-5deb62992467
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
# ╠═ce3067d3-55e6-43a1-9b7b-4cf53c09ee88
# ╠═01e0cfcd-92af-4c04-a7d4-9c9d8fd5a9f1
# ╠═f246ba9f-6cf2-4793-9612-54f2a8475e30
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
# ╠═615f05fb-4907-40bc-9251-3065f565929b
