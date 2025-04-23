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
pace20_all = let
	df = read_fits("processed/rv_pace+20.fits")

	df = df[df.F_match, :]

	df[:, :RV_sigma] = replace(df[:, :RV_sigma], missing=>NaN)
	disallowmissing!(df, :RV_sigma)
	df
end

# ╔═╡ ed332de8-9524-43b6-8863-9742aa9be293
unique(pace20_all.source_id)

# ╔═╡ 5bcb34d8-e9a7-474c-9726-ca8b5f8f1310
spencer18_all_missing = read_fits("processed/rv_spencer+18.fits")

# ╔═╡ 15f2a8e2-90df-48a9-a7bf-e86955f566ce
spencer18_all = let
	df = filter(x->x.F_match,spencer18_all_missing)

	df[:, :RV_sigma] = replace(df[:, :RV_sigma], missing=>NaN)
	disallowmissing!(df, :RV_sigma)

	df
end

# ╔═╡ b7345279-4f80-47ad-a726-537571849eae
graces = let 
	graces = CSV.read("$data_dir/sestito+23_graces.csv", DataFrame)


	filt, idx = RVUtils.xmatch(graces, j24)
	println(filt)
	graces[!, :source_id] = j24.source_id[idx]

	graces[!, :RV_sigma] .= NaN
	graces[!, :RV_count] .= 1
	graces
end

# ╔═╡ bbf49122-11b1-4272-a660-0437c6aa2b3f
md"""
# Combined stars
"""

# ╔═╡ 38a57561-0b99-4cec-8e35-11811bef72a0
all_studies = ["graces", "apogee", "p20", "s18"]

# ╔═╡ fe6fa910-f41e-4657-836b-7eda2f0cddb2
function add_xmatch!(df, new, suffix)
	leftjoin!(df, rename(n->"$(n)_$suffix", new), on="source_id"=>"source_id_$suffix")
end

# ╔═╡ 00a539fb-e9e4-4f9a-8983-f35c3b73d761
md"""
### double checks
"""

# ╔═╡ 8f243944-3d0b-48ce-9f90-8d448c089239
md"""
double check we have all stars
"""

# ╔═╡ 04ff1abd-c584-41de-8c83-7503482c3731
md"""
- are any gaia RVs included?
Gaia stars are too bright to be members :(
"""

# ╔═╡ 22564a47-9b03-4778-b30c-d092581ec107
md"""
# Quality flag
"""

# ╔═╡ 4cda9764-f208-4532-b51f-5deb62992467
md"""
We want to remove stars with statistically large standard deviations.  If our interstudy std is > 5 times sqrt(n) times the standard error, we have a problem
"""

# ╔═╡ 40d26853-0d5e-4d55-9a5e-564da1722f33
function filter_qual_study(rv_meas, study)
	map(eachrow(rv_meas)) do row
		if ismissing(row["RV_sigma_$study"])
			return true
		elseif isnan(row["RV_sigma_$study"])
			return true
		else
			return row["F_scatter_$study"] 
		end

	end
end

# ╔═╡ db3177e7-7132-4f62-846a-f4416a804009
md"""
# write data
"""

# ╔═╡ 7faa2813-e502-4187-855a-047a2f5dd48d
σ_pm = lguys.kms2pm(obs_properties["sigma_v"], obs_properties["distance"])

# ╔═╡ 0e098399-85ae-4cdf-9d82-654a9fd1cd35
icrs0 = RVUtils.icrs(obs_properties)

# ╔═╡ 42fa6c2c-c9cd-4b0d-b1be-701a2db318ee
gsr0 = lguys.transform(lguys.GSR, icrs0)

# ╔═╡ 4d3cd232-81c0-4c32-bfc9-7ba8098d2e4d
icrs0.ra, icrs0.dec, icrs0.pmra, icrs0.pmdec, icrs0.distance, icrs0.radial_velocity

# ╔═╡ b48ef722-072d-4c9c-bcd2-6d7c07bd068d
gsr0.ra, gsr0.dec, gsr0.pmra, gsr0.pmdec, gsr0.distance, gsr0.radial_velocity

# ╔═╡ 24d67978-1013-40ed-b2ff-37eec77c00fc
md"""
- We have to remove a few problematic columns
"""

# ╔═╡ a23adf32-8854-422a-8b7b-2b625132b9f2
import Dates

# ╔═╡ da2b045f-e767-486c-a999-0838a2dbae87


# ╔═╡ c90a03c9-1da2-4bd3-bff8-c620b1bc98a2
md"""
# Numbers
"""

# ╔═╡ 89552b84-d12e-4d88-a58e-8b89ad4b2569
md"""
# Validation for xmatch
"""

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

# ╔═╡ 9e21f642-5315-4e33-a706-ac6cb84ad706
μ0 = obs_properties["radial_velocity"]

# ╔═╡ 47fcf9a4-742d-47b9-96fc-480a4a6b22a0
σ0 = obs_properties["sigma_v"]

# ╔═╡ 395118d9-f171-43e8-a1da-ed9aa72efbad
v_shift_p20 = -1.1

# ╔═╡ d06637b0-4167-40ff-8012-db7e66fe3394
v_shift_s18 = 1.1

# ╔═╡ 070de76a-e917-48ab-8efe-52f665ceca2f
v_shift_apogee = 0

# ╔═╡ e472cbb6-258e-4303-85e2-56f26358c97b
all_stars = let

	all_stars = copy(j24)
	add_xmatch!(all_stars, apogee_all, "apogee")
	add_xmatch!(all_stars, spencer18_all, "s18")
	add_xmatch!(all_stars, pace20_all, "p20")
	add_xmatch!(all_stars, graces, "graces")

	all_stars.RV_apogee .+= v_shift_apogee
	all_stars.RV_s18 .+= v_shift_s18
	all_stars.RV_p20 .+= v_shift_p20
	rename!(all_stars,
		:dr2_radial_velocity => "RV_gaia",
		:dr2_radial_velocity_error => "RV_err_gaia",
	)

	RVUtils.add_rv_means!(all_stars, all_studies)

	all_stars

end

# ╔═╡ 50bf826f-c50c-4465-8a85-2dad998d44cb
filt_is_meas = all_stars.RV_nstudy .> 0

# ╔═╡ 8ca865a2-6521-4e14-ae16-cf5b6a2622ec
rv_meas = all_stars[filt_is_meas, :]

# ╔═╡ 77ccf2bc-692c-42a9-a4fb-0efb7612294f
@assert length(pace20_all.source_id) == sum(pace20_all.F_match) == sum(.!ismissing.(rv_meas.RV_p20))

# ╔═╡ 0c6c6d84-fb14-4eaf-858e-22fb90ab4d8d
@assert length(spencer18_all.source_id) == sum(spencer18_all.F_match) ==  sum(.!ismissing.(rv_meas.RV_s18))

# ╔═╡ f6ebd104-84c7-4f86-ab58-acb65ebc789a
sum(skipmissing(rv_meas.F_scatter_s18))

# ╔═╡ ea333b79-abb0-4827-991f-e2c370bdb0be
sum(skipmissing(.!ismissing.(rv_meas.RV_err_gaia) .* (rv_meas.PSAT .> 001)))

# ╔═╡ 9aafd6b9-6ec8-4b6c-908e-1faac4615e0b
sigma_sigma = ifelse.(isfinite.(rv_meas.RV_sigma), 
	rv_meas.RV_sigma ./ (rv_meas.RV_err .* sqrt.(rv_meas.RV_nstudy)),
	0.				
	)

# ╔═╡ d665895c-9a97-4c36-9ea9-c6e30d0a17ad
mean(sigma_sigma .> 5)

# ╔═╡ 15b1baf7-d9bc-4d3e-a9fa-9d8b8a4dbc6e
F_qual_inter = sigma_sigma .< 5

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

# ╔═╡ dcec54f8-27b4-405e-a386-6e07470dc073
F_qual_study = filter_qual_study(rv_meas, "s18") .& filter_qual_study(rv_meas, "apogee") .& filter_qual_study(rv_meas, "p20")

# ╔═╡ 73a41773-39bc-4bd5-9d99-39804000631a
F_qual = F_qual_study .& F_qual_inter

# ╔═╡ 5c86f548-4838-41e1-b9af-fdd93d900940
sum(.!F_qual)

# ╔═╡ 870571f3-f201-47fd-a215-995a416bc223
sum(.!F_qual_study)

# ╔═╡ 21daa18e-88c3-43ee-8dfa-ea8d116e5ff4
sum(.!ismissing.(rv_meas.RV_graces))

# ╔═╡ 6ea1eff9-14d4-459c-b086-199e0aaf84f0
Δrv = RVUtils.rv_correction(rv_meas.ra, rv_meas.dec, gsr0)

# ╔═╡ 2ecbc309-2b17-4e75-a09b-6483b8b0a711
hist(Measurements.uncertainty.(Δrv) ./ rv_meas.RV_err, bins=100)

# ╔═╡ fd8a7210-ad46-400d-8ef1-d46c75b210c3
df_out = let
	df = copy(rv_meas)

	df[!, :F_scatter] = F_qual
	df[!, :delta_rv] = Measurements.value.(Δrv)
	df[!, :delta_rv_err] = Measurements.uncertainty.(Δrv)

	df
end

# ╔═╡ 57568618-1323-46bc-a8ef-53eb95fc2929
col_filt2 = eltype.(eachcol(df_out)) .∉ Ref(
	[Union{Missing, Dates.Time}, 
	 Missing,
	Union{Missing, String7},
	Union{Missing, String15}])

# ╔═╡ 615f05fb-4907-40bc-9251-3065f565929b
let
	filename = "processed/rv_combined.fits"
	
	write_fits(filename, df_out[:, col_filt2], overwrite=true)
	@info "wrote data"
end

# ╔═╡ 4e52deab-26c0-43e1-b2b0-5b7477737801
names(df_out)[.!col_filt2]

# ╔═╡ 6737a85e-55d5-4f0b-b14e-4866af3bc605
eltype.(eachcol(df_out[:, .!col_filt2]))

# ╔═╡ e08185ba-dba2-421c-801d-7acb387756ef
length(df_out.RV)

# ╔═╡ a7fd3ec2-964b-4cbb-a2dd-dee99cf33c03
median(df_out.RV_err)

# ╔═╡ 53fcc27a-700b-4b67-9ba6-3d4613edfb7b
sum(df_out.RV_count)

# ╔═╡ 5f3f4b6b-11c3-4c87-96e5-22bf6170f170
sum(df_out.F_scatter)

# ╔═╡ 0d2dbc73-1ded-46e3-b142-9bc7b777728d
all_stars[.!ismissing.(all_stars.RV_graces), [:ra, :dec, :source_id, :pmra, :pmdec, :PSAT]]

# ╔═╡ bd51fe42-7e39-4cd8-8065-58ab9814f966
 @assert sum(apogee_all.F_match) == sum(.!ismissing.(all_stars.RV_apogee)) 

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
plot_xmatch_radec("p20")

# ╔═╡ c218cfa8-2f55-4957-bcdd-8b3970fe639a
plot_xmatch_radec("s18")

# ╔═╡ f7213298-3dcb-49ac-a9f0-a129a03423aa
plot_xmatch_radec("graces")

# ╔═╡ 609f55cd-aa49-4a5e-b871-413ec7ef0990
import StatsBase as sb

# ╔═╡ 0cf94e27-b200-4516-9495-dcac7f10d9a0
sum(.!ismissing.(rv_meas.RV_apogee) .& ismissing.(rv_meas.RV_s18) .& ismissing.(rv_meas.RV_p20)) # most apogee stars are covered by other surveys.

# ╔═╡ 3655b6a0-ac9a-4a49-86fd-6a6484409819
md"""
## Check RV means make sense
"""

# ╔═╡ 7f13339e-6a0c-4944-8a36-5ab136fd8415
function compare_rv_mean(study1, rv_meas=rv_meas; low=μ0-3σ0, high=μ0+3σ0)
	rv1  = rv_meas[:, "RV_$study1"]
	rv1_err  = rv_meas[:, "RV_err_$study1"]
	
	rv2  = rv_meas[:, "RV"]
	rv2_err  = rv_meas[:, "RV_err"]

	filt = filt_missing(rv1, true; low=low, high=high) .& filt_missing(rv2, true;  low=low, high=high)

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
compare_rv_mean("graces", )

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

# ╔═╡ 03462abb-32d6-433c-a9c4-65ef2dd58965
function is_rv_member(rv, err)
	return abs(rv - rv_mean) / (err ⊕ sigma_los) < 3
end

# ╔═╡ 82d4e049-7f47-4f8b-a3a0-7191ef82912b
filt_rv = is_rv_member.(rv_meas.RV, rv_meas.RV_err) 

# ╔═╡ 3920e27f-48fd-41a5-bf53-1f80edf3d56d
memb_stars = rv_meas[rv_meas.PSAT .> 0.2 .&& F_qual .&& filt_rv , :]

# ╔═╡ 9213cc36-74d1-452f-bd9a-eb5c5cab1f87
function compare_rv(study1, study2)
	rv1  = memb_stars[:, "RV_$study1"]
	rv1_err  = memb_stars[:, "RV_err_$study1"]
	
	rv2  = memb_stars[:, "RV_$study2"]
	rv2_err  = memb_stars[:, "RV_err_$study2"]

	filt = filt_missing(rv1, true; low=μ0-3σ0, high=μ0+3σ0) .& filt_missing(rv2, true;  low=μ0-3σ0, high=μ0+3σ0)

	println("matched ", sum(filt), " stars")


	if sum(filt) == 0
		@error ("nothing to plot")
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
	μ = sb.median(rv1[filt] .- rv2[filt])
	δμ = sqrt(1/sum(w))
	sigma = sb.mad(rv1[filt] .- rv2[filt] .- μ)

	text!(0.1, 0.8, space=:relative, text="μ = $(μ ± δμ)\nσ = $(round(sigma, digits=2))", fontsize=10)

	
	lines!([μ0-3σ0, μ0+3σ0], [μ0-3σ0, μ0+3σ0], color=:black)
	return fig
end

# ╔═╡ f9e1bcbd-2def-4e64-b69e-d1cfcf7ffedd
@savefig "p20_s18_xmatch" compare_rv("p20", "s18")

# ╔═╡ f61debe4-8b23-4415-b77c-43c4465ccfcb
@savefig "s18_apogee" compare_rv("s18", "apogee")

# ╔═╡ 102c73ef-2c95-4784-80df-ed0312511c00
@savefig "apogee_w09_xmatch" compare_rv("apogee", "p20")

# ╔═╡ a2370516-e7ec-4502-ae4b-b111bcf68d36
compare_rv_mean("apogee", memb_stars)

# ╔═╡ aaaf5ba2-c9ed-41ec-a22a-d78ed96fd84e
compare_rv_mean("p20", memb_stars)

# ╔═╡ 78e7aff0-3658-4d64-b117-58189a85307a
compare_rv_mean("s18", memb_stars)

# ╔═╡ 92bfc52c-05f5-45e3-ae26-90c579ecfe2c
hist(filter((!) ∘ ismissing, memb_stars.RV_err_p20 ./ memb_stars.RV_err))

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
sum(.!ismissing.(memb_stars.RV_graces))

# ╔═╡ d0608300-496d-432a-9016-05e9f75df0a4
unique(memb_stars.RV_err_graces)

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
	scatter!(memb_stars.R_ell, memb_stars.RV_graces, label="Sestito+23", markersize=5)
	scatter!(memb_stars.R_ell, memb_stars.RV_p20, label="pace+20", markersize=3)
	scatter!(memb_stars.R_ell, memb_stars.RV_s18, label="spencer+18", markersize=3)

	Legend(fig[1,2], ax)
	#axislegend(ax, position=:rb)


	resize_to_layout!(fig)

	@savefig "rv_scatter_alstudy"
	fig
end

# ╔═╡ a70a2546-3aaa-4d78-a060-2a9a39940dc8
mean(.!ismissing.(memb_stars.RV_apogee))

# ╔═╡ 382e7503-0e31-49d5-b593-a6a1a0f9c5e0
mean(.!ismissing.(memb_stars.RV_p20))

# ╔═╡ 0b7107f8-38d2-4b0c-aad1-e07d3927cf87
mean(.!ismissing.(memb_stars.RV_s18))

# ╔═╡ e32f0d5a-9a6a-407c-a466-586c5d63fda8
mean(ismissing.(memb_stars.RV_s18) .& .!ismissing.(memb_stars.RV_p20))

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
	scatter!(rv_meas.R_ell, rv_meas.RV_graces, label="sestito + 23", markersize=3)
	scatter!(rv_meas.R_ell, rv_meas.RV_s18, label="spencer + 18", markersize=3)
	scatter!(rv_meas.R_ell, rv_meas.RV_p20, label="pace + 20", markersize=3)

	Legend(fig[1,2], ax)


	resize_to_layout!(fig)
	fig
end

# ╔═╡ 62301344-869a-4ec0-8299-29f0ff5d1c15
sum(sum(eachcol(ismissing.(rv_meas))))

# ╔═╡ 1deb1520-194a-40b5-8968-bf73b756ba3d
rv_meas[ .! ismissing.(rv_meas.RV_graces), :PSAT]

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
# ╟─300aa043-3750-4d01-ae48-1bf765cd92b1
# ╠═7a50a176-96a5-4098-88d6-0fa2874d0f90
# ╠═77e7884c-0360-4b7f-b9ad-e81be2516552
# ╟─6964faa7-17d2-412c-b4a2-5f981f8c4b54
# ╠═d4d0a488-1c0e-4bc6-9a88-94b27d84e8ce
# ╠═9ed91d41-e650-4d3c-9e74-100ff3d57d82
# ╠═ed332de8-9524-43b6-8863-9742aa9be293
# ╠═5bcb34d8-e9a7-474c-9726-ca8b5f8f1310
# ╠═15f2a8e2-90df-48a9-a7bf-e86955f566ce
# ╠═b7345279-4f80-47ad-a726-537571849eae
# ╟─bbf49122-11b1-4272-a660-0437c6aa2b3f
# ╠═38a57561-0b99-4cec-8e35-11811bef72a0
# ╠═e472cbb6-258e-4303-85e2-56f26358c97b
# ╠═fe6fa910-f41e-4657-836b-7eda2f0cddb2
# ╠═50bf826f-c50c-4465-8a85-2dad998d44cb
# ╠═8ca865a2-6521-4e14-ae16-cf5b6a2622ec
# ╟─00a539fb-e9e4-4f9a-8983-f35c3b73d761
# ╟─8f243944-3d0b-48ce-9f90-8d448c089239
# ╠═0d2dbc73-1ded-46e3-b142-9bc7b777728d
# ╠═bd51fe42-7e39-4cd8-8065-58ab9814f966
# ╠═77ccf2bc-692c-42a9-a4fb-0efb7612294f
# ╠═0c6c6d84-fb14-4eaf-858e-22fb90ab4d8d
# ╠═f6ebd104-84c7-4f86-ab58-acb65ebc789a
# ╟─04ff1abd-c584-41de-8c83-7503482c3731
# ╠═ea333b79-abb0-4827-991f-e2c370bdb0be
# ╟─22564a47-9b03-4778-b30c-d092581ec107
# ╟─4cda9764-f208-4532-b51f-5deb62992467
# ╠═9aafd6b9-6ec8-4b6c-908e-1faac4615e0b
# ╠═d665895c-9a97-4c36-9ea9-c6e30d0a17ad
# ╠═15b1baf7-d9bc-4d3e-a9fa-9d8b8a4dbc6e
# ╠═73a41773-39bc-4bd5-9d99-39804000631a
# ╠═dcec54f8-27b4-405e-a386-6e07470dc073
# ╠═870571f3-f201-47fd-a215-995a416bc223
# ╠═5c86f548-4838-41e1-b9af-fdd93d900940
# ╠═40d26853-0d5e-4d55-9a5e-564da1722f33
# ╠═06b5c8d8-e531-4f00-a3cf-8d6202bb71f2
# ╟─db3177e7-7132-4f62-846a-f4416a804009
# ╠═7faa2813-e502-4187-855a-047a2f5dd48d
# ╠═a6bc7223-3267-4760-9b6a-d886ac6f4544
# ╠═02c55b13-5681-475b-bead-9e0e0b9e9656
# ╠═21daa18e-88c3-43ee-8dfa-ea8d116e5ff4
# ╠═0e098399-85ae-4cdf-9d82-654a9fd1cd35
# ╠═42fa6c2c-c9cd-4b0d-b1be-701a2db318ee
# ╠═6ea1eff9-14d4-459c-b086-199e0aaf84f0
# ╠═4d3cd232-81c0-4c32-bfc9-7ba8098d2e4d
# ╠═b48ef722-072d-4c9c-bcd2-6d7c07bd068d
# ╠═2ecbc309-2b17-4e75-a09b-6483b8b0a711
# ╠═fd8a7210-ad46-400d-8ef1-d46c75b210c3
# ╠═615f05fb-4907-40bc-9251-3065f565929b
# ╠═24d67978-1013-40ed-b2ff-37eec77c00fc
# ╠═a23adf32-8854-422a-8b7b-2b625132b9f2
# ╠═da2b045f-e767-486c-a999-0838a2dbae87
# ╠═57568618-1323-46bc-a8ef-53eb95fc2929
# ╠═4e52deab-26c0-43e1-b2b0-5b7477737801
# ╠═6737a85e-55d5-4f0b-b14e-4866af3bc605
# ╠═c90a03c9-1da2-4bd3-bff8-c620b1bc98a2
# ╠═e08185ba-dba2-421c-801d-7acb387756ef
# ╠═a7fd3ec2-964b-4cbb-a2dd-dee99cf33c03
# ╠═53fcc27a-700b-4b67-9ba6-3d4613edfb7b
# ╠═5f3f4b6b-11c3-4c87-96e5-22bf6170f170
# ╟─89552b84-d12e-4d88-a58e-8b89ad4b2569
# ╠═e6f2de3b-ce32-4d61-851f-4e42fcce95c0
# ╠═5b2fceff-9c3e-472d-9310-31e920137e41
# ╠═4b305b83-1a3b-48a6-b19f-6f3ebed0768f
# ╠═c218cfa8-2f55-4957-bcdd-8b3970fe639a
# ╠═f7213298-3dcb-49ac-a9f0-a129a03423aa
# ╟─f7ec8bba-9f45-435b-b67c-33182e992dfd
# ╠═93d185f2-1eaa-4d35-87dd-b84f385483de
# ╠═9e21f642-5315-4e33-a706-ac6cb84ad706
# ╠═47fcf9a4-742d-47b9-96fc-480a4a6b22a0
# ╠═9213cc36-74d1-452f-bd9a-eb5c5cab1f87
# ╠═395118d9-f171-43e8-a1da-ed9aa72efbad
# ╠═d06637b0-4167-40ff-8012-db7e66fe3394
# ╠═070de76a-e917-48ab-8efe-52f665ceca2f
# ╠═f9e1bcbd-2def-4e64-b69e-d1cfcf7ffedd
# ╠═609f55cd-aa49-4a5e-b871-413ec7ef0990
# ╠═f61debe4-8b23-4415-b77c-43c4465ccfcb
# ╠═0cf94e27-b200-4516-9495-dcac7f10d9a0
# ╠═102c73ef-2c95-4784-80df-ed0312511c00
# ╟─3655b6a0-ac9a-4a49-86fd-6a6484409819
# ╠═7f13339e-6a0c-4944-8a36-5ab136fd8415
# ╠═3920e27f-48fd-41a5-bf53-1f80edf3d56d
# ╠═a2370516-e7ec-4502-ae4b-b111bcf68d36
# ╠═aaaf5ba2-c9ed-41ec-a22a-d78ed96fd84e
# ╠═78e7aff0-3658-4d64-b117-58189a85307a
# ╠═92bfc52c-05f5-45e3-ae26-90c579ecfe2c
# ╠═f4f9dd06-1a1a-458b-be75-05d52623580c
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
# ╠═d0608300-496d-432a-9016-05e9f75df0a4
# ╠═33f54afc-cdb9-4eb8-887f-5a43281b837c
# ╠═a70a2546-3aaa-4d78-a060-2a9a39940dc8
# ╠═382e7503-0e31-49d5-b593-a6a1a0f9c5e0
# ╠═0b7107f8-38d2-4b0c-aad1-e07d3927cf87
# ╠═e32f0d5a-9a6a-407c-a466-586c5d63fda8
# ╠═500815e1-f9e1-4401-ba2d-72f326cfa783
# ╠═62301344-869a-4ec0-8299-29f0ff5d1c15
# ╠═1deb1520-194a-40b5-8968-bf73b756ba3d
