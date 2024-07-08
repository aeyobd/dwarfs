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
	using FITSIO
	using Measurements

	import LilGuys as lguys
end

# ╔═╡ 9070c811-550c-4c49-9c58-0943b0f808b2
using Turing

# ╔═╡ e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
using PairPlots

# ╔═╡ d4123ffd-cb32-486e-a93d-f48d7112a831
include("density_fits/filter_utils.jl")

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

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median

# ╔═╡ dd29be70-4918-47e4-98cf-3608df14e88a
begin
	ell = 0.37
	PA = -94
	ra0 = 15.039170
	dec0 = -33.709180
end

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 49e728ee-e545-4ee5-bd9f-b0a38c4eaf15
function ang_dist(ra1, dec1, ra2, dec2)
	dra = (ra1 .- ra2) .* cosd.(dec1)
	ddec = dec1 .- dec2
	return @. sqrt(dra^2 + ddec^2)
end

# ╔═╡ 3d8b3c55-153b-4a4a-9289-78f34df07abc
"""
Cross matches 2 dataframes given angular seperation in arcmin
"""
function xmatch(df1::DataFrame, df2::DataFrame, max_sep=2)
	max_sep = max_sep / 3600
	dists = ang_dist.(df1.ra, df1.dec, Ref(df2.ra), Ref(df2.dec))
	filt = minimum.(dists) .< max_sep
	idxs = argmin.(dists)
	return filt, idxs
end

# ╔═╡ 248c2d5f-85cc-44be-bc63-43d4af470182
begin 
	params = read_file("density_fits/sculptor/fiducial.toml")
	params["filename"] = "../data/Sculptor.GAIASOURCE.RUWE.VELS.PROB.fits"
	params["PA"] = 92
	params["ellipticity"] = 0.36
	params["PSAT_min"] = nothing
	params = DensityParams(params)
end

# ╔═╡ 7a50a176-96a5-4098-88d6-0fa2874d0f90
begin #TODO: implement J24 crossmatching
	j24 = load_stars(params.filename, params) # does not filter
end

# ╔═╡ 7db47590-b82b-4822-8a01-eaff96c9389f
md"""
## Dart
"""

# ╔═╡ 23402e30-2064-494c-bd4a-8243cb474b61
begin 
	dart = CSV.read("../data/Sculptor_DARTS_BEST_psat40.csv", DataFrame)
	rename!(dart, 
		"vel"=>"RV",
		"evel"=>"RV_err",
		"feh" => "Fe_H",		
		# "feh_lo" => "Fe_H_lo",		
		# "feh_hi" => "Fe_H_hi",		
	)

	dart = dart[!, [:ra, :dec, :source_id, :RV, :RV_err, ]]
end

# ╔═╡ baacb491-9dff-4ff4-b62c-5842d79794da
hist(dart.RV)

# ╔═╡ c470c2b9-093d-42ab-b96d-9dac231ccabc
md"""
## Apogee
"""

# ╔═╡ 27063de7-01d4-48a8-a06f-cc24aec662d2
md"""
### Check apogee-dart xmatch against fed
"""

# ╔═╡ b6de4afb-9362-4618-a114-df460031e4f9
md"""
# GMOS alla Federico
"""

# ╔═╡ b7345279-4f80-47ad-a726-537571849eae
let
	global gmos = CSV.read("../data/targets_stellar_met.csv", DataFrame)
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
	FITS("../data/tolstoy_2023.fits", "r") do f
		global tolstoy23_all = DataFrame(f[2])
	end

	rename!(tolstoy23_all, 
		:Vlos => :RV,
		:e_Vlos => :RV_err,
		:RAJ2000 => :ra,
		:DEJ2000 => :dec,
		:GaiaDR3 => :source_id
	)


	tolstoy23_all = tolstoy23_all
end

# ╔═╡ 180fac98-678c-4a14-966c-385387c60ac3
md"""
## Walker et al. (2009)

Data downloaded from CDS associated with observation paper.
Should have 1365 scl members.
Weighted mean of measurements for repeated stars
"""

# ╔═╡ 77830e77-50a4-48b7-b070-8fbd7508c173
let 
	global walker09_single
	FITS("../data/walker_2009.fits", "r") do f
		global walker09_single = DataFrame(f[2])
	end
	walker09_single[!, "Galaxy"] = [s[1:3] for s in walker09_single.Target]

end

# ╔═╡ dee71790-ffeb-477d-adbe-112731dfe134
let 
	global walker09_averaged
	FITS("../data/walker+09.fits", "r") do f
		global walker09_averaged = DataFrame(f[2])
	end

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

# ╔═╡ 91ca7d42-4537-4777-b6d0-ebfd964e0171
unique(walker09_single.Target[walker09_single.Galaxy .== "Scl"])

# ╔═╡ fe6fa910-f41e-4657-836b-7eda2f0cddb2
function add_xmatch!(df, new, suffix)
	leftjoin!(df, rename(n->"$(n)_$suffix", new), on="source_id"=>"source_id_$suffix")
end

# ╔═╡ 89552b84-d12e-4d88-a58e-8b89ad4b2569
md"""
# Validation for xmatch
"""

# ╔═╡ f7ec8bba-9f45-435b-b67c-33182e992dfd
md"""
# Cross study RV
"""

# ╔═╡ d3333b48-aa4e-42c1-9e0a-bbff98e3647d
all_studies = ["dart", "gmos", "apogee", "w09", "t23"]

# ╔═╡ 222bb254-8b65-44d3-b3d2-b08fbcbe9950
all_studies

# ╔═╡ c724e720-18ca-4905-a4b6-39fc47abe39d
j24[j24.source_id .== 5006419626331394048, :]

# ╔═╡ 733fe42e-b7a5-4285-8c73-9a41e4488d40


# ╔═╡ c94b959d-2490-4a16-9ab7-c5c580316100


# ╔═╡ 3655b6a0-ac9a-4a49-86fd-6a6484409819
md"""
## Check RV means make sense
"""

# ╔═╡ 6734991c-16c0-4424-a2bb-84bfa811121f
md"""
## MCMC Priors
"""

# ╔═╡ 1b97d0e5-7a77-44a5-b609-ed8945cd959c
"""
Fits a normal (gaussian) distribution to 1d data with errors.
"""
@model function normal_dist(x, xerr; μ_min=90, μ_max=120)
	μ ~ Uniform(μ_min, μ_max)
	σ ~ LogNormal(2.5, 1) # approx 1 - 100 km / s, very broad but should cover all
	s = @. sqrt(σ^2 + xerr^2)

	x ~ MvNormal(fill(μ, length(x)), s)
end

# ╔═╡ abbd2a53-e077-4af7-a168-b571e1a906b8
xlabel = L"radial velocity / km s$^{-1}$"

# ╔═╡ 55504311-fbbd-4351-b9e2-ecab053ba207
function plot_samples!(samples, x;
		thin=10, color=:black, alpha=nothing, kwargs...)

	alpha = 1 / (size(samples, 1))^(1/3)
	for sample in eachrow(samples)[1:thin:end]
		y = lguys.gaussian.(x, sample.μ, sample.σ)
		lines!(x, y, color=color, alpha=alpha)
	end
end

# ╔═╡ 95ef12d2-174a-47c3-a106-9b4005b2a80d
md"""
# Full fits to data
"""

# ╔═╡ 9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
function fit_rv_sigma(rv, rv_err; N=3_000, p=0.16, burn=0.2)
	samples = DataFrame(sample(normal_dist(rv, rv_err), NUTS(0.65), N))
	burn = round(Int, burn * N)
	samples=samples[burn:end, :]

	μ = median(samples.μ)
	μ_p = quantile(samples.μ, [p, 1-p]) 
	μ_err = (μ - μ_p[1], μ_p[2] - μ)

	σ = median(samples.σ)
	σ_p = quantile(samples.σ, [p, 1-p])
	σ_err = (σ - σ_p[1], σ_p[2] - σ)

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

# ╔═╡ cfaafe3a-54bb-4e16-8b02-ef05fe7e4431
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

# ╔═╡ b5faa873-d5c9-447b-b90e-f693db26e6c2
md"""
### Comparing sample properties
"""

# ╔═╡ 7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
md"""
## Binned properties
"""

# ╔═╡ c2735c49-2892-46ac-bcf8-7cdcef409f44
function binned_mu_sigma(x, y, yerr, bins)

	if !issorted(bins)
		error("bins must be sorted")
	end

	N = length(bins) - 1
	μs = Vector{Float64}(undef, N)
	σs = Vector{Float64}(undef, N)
	μ_errs = Vector{Tuple{Float64, Float64}}(undef, N)
	σ_errs = Vector{Tuple{Float64, Float64}}(undef, N)

	
	for i in 1:N
		filt = x .>= bins[i]
		filt .&= x .< bins[i+1]
		@info "calculating bin $i"

		μs[i], σs[i], μ_errs[i], σ_errs[i] = fit_rv_sigma(y[filt], yerr[filt])
	end

	return μs, σs, μ_errs, σ_errs
end	

# ╔═╡ 82a0e58a-30a4-4e42-b9c1-cb184eb551aa
md"""
# Misc plots
"""

# ╔═╡ db78b22f-bee8-4d36-994e-6c7f9e96c9f2
value(x::Measurement) = x.val

# ╔═╡ 12ad8e67-0c7b-4029-91f6-5a2453ec799a
err(x::Measurement) = x.err

# ╔═╡ e976d103-9956-47f3-adb0-7c449214b9d9
import StatsBase as sb

# ╔═╡ 9b4a0a1f-4b0c-4c90-b871-2bd244f0a908
not = !

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

# ╔═╡ 1864bc12-cbc2-4fd4-84a8-6bb300f17c1b
filt_tolstoy = sigma_clip(tolstoy23_all.RV, 3)

# ╔═╡ b36c2a37-2359-4a1a-98fc-3a1cd17fd790
tolstoy23 = tolstoy23_all[filt_tolstoy, :]

# ╔═╡ 1d01f9a2-2b77-4ad0-9893-b31562de7924
filt_tolstoy .&= tolstoy23_all.Mem .== "m"

# ╔═╡ bb7f6769-ec92-460d-8423-449029175f79
begin 
	apogee_all = CSV.read("../data/sculptor_apogeeDR17_xmatch.csv", DataFrame)

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

# ╔═╡ 1bc2541c-69ac-40aa-94f5-743707ca5a66
hist(apogee.RV)

# ╔═╡ 48c62811-136f-4962-a42c-b1dd1fc74f8c
begin 
	apogee_notdart = CSV.read("../data/sculptor_apogeeDR17_xmatch_notDART.csv", DataFrame)

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

# ╔═╡ fcf50b2c-5703-4101-b3a7-45e5fb8e072f
size(a2, 1) + size(dart, 1) # should be 617

# ╔═╡ f8775eb1-2cb9-4d81-8c02-39c43ceb9b45
sort(a2.source_id) == sort(apogee_notdart.source_id)

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

# ╔═╡ 8c59d607-d484-497c-9ee7-751a4edd4992
hist(walker09.RV)

# ╔═╡ ebdf8cfc-ebff-4331-bab2-998a734f2cd1
datasets = Dict(
	:apogee => apogee,
	:dart => dart,
	:walker => walker09,
	:tolstoy => tolstoy23
)

# ╔═╡ ec4bc80f-22e6-49f9-a589-5c5bc9e50a8b
props = Dict(name => fit_rv_sigma(float.(data.RV), float.(data.RV_err)) for (name, data) in datasets)

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

# ╔═╡ eaf47845-95dd-4f7d-bf76-3d9b1154711a
for key in names(props)
	println(key)
	println(props[key])
	@info plot_sample_normal_fit(datasets[key], props[key], title=string(key))
end

# ╔═╡ e472cbb6-258e-4303-85e2-56f26358c97b
let
	global all_stars 

	all_stars = copy(j24)
	add_xmatch!(all_stars, dart, "dart")
	add_xmatch!(all_stars, apogee_all, "apogee")
	add_xmatch!(all_stars, walker09_all, "w09")
	add_xmatch!(all_stars, tolstoy23_all, "t23")
	add_xmatch!(all_stars, gmos, "gmos")


	rename!(all_stars,
		:dr2_radial_velocity => "RV_gaia",
		:dr2_radial_velocity_error => "RV_err_gaia",
	)

end

# ╔═╡ 88ed48b4-baa9-4ac0-86e1-8348edcd59b4
begin 
	rvs = [all_stars[:, "RV_$study"] for study in all_studies]
	rv_errs = [all_stars[:, "RV_err_$study"] for study in all_studies]

	rvs = hcat(rvs...)
	rv_errs = hcat(rv_errs...)
end

# ╔═╡ 6bc02c4f-31a2-4e4e-8612-e66f8cc9c93e
rvs

# ╔═╡ 11fcf4f8-fcd5-4426-a4e9-b046138bde1b
rv_errs

# ╔═╡ e3f05ee2-cc5f-437e-801d-3c7d842af709
all_stars[:, "RV_w09"]

# ╔═╡ 0d2dbc73-1ded-46e3-b142-9bc7b777728d
all_stars.RV_gmos[not.(ismissing.(all_stars.RV_gmos))]

# ╔═╡ 1510b6de-09ae-474e-91c4-eea4e2cacdce
sum(not.(ismissing.(all_stars.RV_dart))) == length(dart.RV)

# ╔═╡ bd51fe42-7e39-4cd8-8065-58ab9814f966
sum(not.(ismissing.(all_stars.RV_apogee))) == length(apogee_all.RV)

# ╔═╡ ab49efb3-2ab7-47d0-a3a5-a342c789aa9b
sum(not.(ismissing.(all_stars.RV_w09))) , length(walker09_all.RV)

# ╔═╡ de762a39-b430-4452-ba87-8b8cf1ad9852
sum(not.(ismissing.(all_stars.RV_t23))) == length(tolstoy23_all.RV)

# ╔═╡ e6f2de3b-ce32-4d61-851f-4e42fcce95c0
function plot_xmatch_radec(suffix)

	ra2 = all_stars[:, "ra_$suffix"]
	filt = not.(ismissing.(ra2))
	ra2 = ra2[filt]

	ra1 = all_stars.ra[filt]
	dec1 = all_stars.dec[filt]
	dec2 = all_stars[filt, "dec_$suffix"]

	fig, ax = FigAxis()

	scatter!(ra1, dec1)
	scatter!(float.(ra2), float.(dec2), markersize=5)

	fig
end

# ╔═╡ 5b2fceff-9c3e-472d-9310-31e920137e41
plot_xmatch_radec("apogee")

# ╔═╡ 7091dc6b-dd77-4f92-bd35-def8c7384f00
plot_xmatch_radec("dart")

# ╔═╡ 4b305b83-1a3b-48a6-b19f-6f3ebed0768f
plot_xmatch_radec("t23")

# ╔═╡ c218cfa8-2f55-4957-bcdd-8b3970fe639a
plot_xmatch_radec("w09")

# ╔═╡ f7213298-3dcb-49ac-a9f0-a129a03423aa
plot_xmatch_radec("gmos")

# ╔═╡ 93d185f2-1eaa-4d35-87dd-b84f385483de
function filt_missing(col, verbose=false; low=-Inf, high=Inf)
	filt =  not.(ismissing.(col)) .& not.(isnan.(col))
	filt1 = high .> col .> low
	if verbose
		println("excluding ", sum(not.(filt1)[filt]), " outliers")
	end
	return filt .& filt1
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

# ╔═╡ fae534ec-78a9-4f40-91d0-8511b0bab399
compare_rv("dart", "apogee")

# ╔═╡ 0df388a8-838b-46e6-a281-c48197459e37
compare_rv("dart", "t23")

# ╔═╡ fcb7c823-4584-4766-9709-1661aa8f87a7
compare_rv("dart", "w09")

# ╔═╡ f9e1bcbd-2def-4e64-b69e-d1cfcf7ffedd
compare_rv("t23", "w09")

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

# ╔═╡ 877706f1-08a0-496a-890d-622e3a2fd9ec
RV_a_std = [safe_weighted_mean(rvs[i, :], rv_errs[i, :]) for i in 1:size(rvs, 1)]

# ╔═╡ a13a33f1-65c4-4217-8636-639b1e14d109
sum(filt_missing(all_stars[:, "RV_err_w09"]))

# ╔═╡ ee3c22db-6b6b-4c30-8d1f-86b103c018fc
sum(filt_missing(walker09_all.RV_err))

# ╔═╡ 73bd1553-e2ae-4bfb-aac1-0880346f5054
begin
	filt_is_meas = not.(ismissing.(RV_a_std))
	rv_meas = all_stars[filt_is_meas, :]
	RV = value.(first.(RV_a_std[filt_is_meas]))
	RV_err = err.(first.(RV_a_std[filt_is_meas]))
	RV_std = (last.(RV_a_std[filt_is_meas]))

	rv_meas[!, :RV] = RV
	rv_meas[!, :RV_err] = RV_err
	rv_meas[!, :RV_std] = RV_std

end

# ╔═╡ 01b3135f-7a72-4669-b586-4bc5894464ad
sum(filt_is_meas)

# ╔═╡ 9e63c171-394e-4c05-b237-189fba0274e2
begin 
	memb_stars = rv_meas[(rv_meas.PSAT .> 0.2) .& (150 .> rv_meas.RV .> 60), :]

	# only Fed. data
	# memb_stars = memb_stars[not.(ismissing.(memb_stars.RV_apogee)) .| not.(ismissing.(memb_stars.RV_dart)), :]
end

# ╔═╡ 7f4a5254-ed6f-4faa-a71e-4b4986a99d45
hist(memb_stars.RV)

# ╔═╡ a1938588-ca40-4844-ab82-88c4254c435b
length(memb_stars.RV)

# ╔═╡ 5cf336f6-e3eb-4668-b074-18b396f027be
prior_samples = DataFrame(
	sample(normal_dist(memb_stars.RV, memb_stars.RV_err), Prior(), 10000)
)

# ╔═╡ ccbb694d-f973-40fd-bab7-a2aefbd9fb0b
pairplot(prior_samples[:, [:μ, :σ]])

# ╔═╡ 74ad07df-f15c-436f-b390-ce95b27f7fab
let
	fig, ax = FigAxis(
		xgridvisible=false,
		ygridvisible=false,
		xlabel=xlabel,
		ylabel="density",
		title="priors"
	)
	
	plot_samples!(prior_samples, LinRange(80, 130, 100))

	fig
end

# ╔═╡ 318d29b9-4c84-4d38-af6d-048518952970
samples = DataFrame(sample(normal_dist(memb_stars.RV, memb_stars.RV_err), NUTS(0.65), 10000))

# ╔═╡ b18e4622-41e0-4700-9e4b-3dbebeefea53
describe(samples)

# ╔═╡ bc7bd936-3f62-4430-8acd-8331ca3ee5ad
pairplot(samples[:, [:μ, :σ]])

# ╔═╡ 0f5a9d9e-c5ca-4eb6-a0d2-5bb39b81daf6
σ_m = median(samples.σ)

# ╔═╡ 319bd778-7e17-4bd7-856f-d6785b287219
quantile(samples.σ, [0.16, 0.84]) .- σ_m

# ╔═╡ 764b5306-20f9-4810-8188-1bdf9482260f
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = Arya.histogram(Float64.(memb_stars.RV), 30, normalization=:pdf)
	
	plot_samples!(samples, LinRange(70, 150, 100), thin=15)
	errscatter!(midpoints(h.bins), h.values, yerr=h.err, color=COLORS[6])

	fig
end

# ╔═╡ 61e15c47-c454-48db-95f9-02abe052676e
mean(memb_stars.RV)

# ╔═╡ d5938fc3-9c8a-4e33-8401-500b4201df14
sb.sem(memb_stars.RV)

# ╔═╡ 6514815e-e1f3-44a1-8c0b-a90f13e0077e
std(memb_stars.RV)

# ╔═╡ b654f5d5-ff48-4bae-bd69-a9db2673282b
sesd(memb_stars.RV)

# ╔═╡ 8b21cc49-ca17-4844-8238-e27e9752bee7
bins = Arya.bins_equal_number(memb_stars.r_ell, n=8)

# ╔═╡ 1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
bin_errs = diff(bins) / 2

# ╔═╡ c50f68d7-74c3-4c36-90c5-a5262982ed9f
μs, σs, μ_errs, σ_errs = binned_mu_sigma(memb_stars.r_ell, memb_stars.RV, memb_stars.RV_err, bins)

# ╔═╡ 614f3f09-1880-491d-b41e-4e229330d66f
let
	fig, ax = FigAxis(
		xlabel = "R / arcmin",
		ylabel = L"$\sigma_{v, \textrm{los}}$ / km s$^{-1}$"
	)

	errscatter!(midpoints(bins), σs, yerr=σ_errs, xerr=bin_errs, color=:black)
	hlines!(σ_m)

	fig
end

# ╔═╡ 33f54afc-cdb9-4eb8-887f-5a43281b837c
let
	fig, ax = FigAxis(
		xlabel = "R / arcmin",
		ylabel = L"RV / km s$^{-1}$",
		limits=(nothing, (60, 150))
	)

	scatter!(memb_stars.r_ell, memb_stars.RV_dart, label="DART")
	scatter!(memb_stars.r_ell, memb_stars.RV_apogee, label="APOGEE")
	scatter!(memb_stars.r_ell, memb_stars.RV_gmos, label="GMOS")
	scatter!(memb_stars.r_ell, memb_stars.RV_w09, label="Walker+09")
	scatter!(memb_stars.r_ell, memb_stars.RV_t23, label="tolstoy + 23")

	errscatter!(midpoints(bins), μs, yerr=μ_errs, color=:black)
	
	errorbars!(midpoints(bins), μs .+ σs, bin_errs, direction = :x, color=:black)
	errorbars!(midpoints(bins), μs .- σs, bin_errs, direction = :x, color=:black)

	Legend(fig[1,2], ax)


	fig
end

# ╔═╡ 86776e68-d47f-43ed-b37f-432c864050bb
let
	fig, ax = FigAxis(
		xlabel = "R / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)

	scatter!(memb_stars.r_ell, memb_stars.RV, color=COLORS[3], alpha=0.1)


	errscatter!(midpoints(bins), μs, yerr=μ_errs, color=:black)
	
	errorbars!(midpoints(bins), μs .+ σs, bin_errs, direction = :x, color=:black)
	errorbars!(midpoints(bins), μs .- σs, bin_errs, direction = :x, color=:black)


	fig
end

# ╔═╡ b7ebd916-bffc-4ffc-a7f7-44ee315e2528
stephist(memb_stars.r_ell, bins=bins)

# ╔═╡ 3a69f395-3c2d-4357-89af-5963d5fa79b8
let
	fig, ax = FigAxis(
		xlabel=L"$\log r_\textrm{ell}$ / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)
	
	scatter!(log10.(memb_stars.r_ell), memb_stars.RV)


	fig
end

# ╔═╡ 7f13339e-6a0c-4944-8a36-5ab136fd8415
function compare_rv_mean(study1)
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
compare_rv_mean("w09")

# ╔═╡ a2370516-e7ec-4502-ae4b-b111bcf68d36
compare_rv_mean("apogee")

# ╔═╡ 78e7aff0-3658-4d64-b117-58189a85307a
compare_rv_mean("t23")

# ╔═╡ f4f9dd06-1a1a-458b-be75-05d52623580c
compare_rv_mean("gmos")

# ╔═╡ 6b84c679-0532-465a-97dd-d62671077b61
compare_rv_mean("dart")

# ╔═╡ d11edca7-b9ae-4269-9e1b-661d59bd965e
all_stars[not.(ismissing.(all_stars.RV_gmos)), :].source_id

# ╔═╡ 74152829-27ef-4d8d-8b32-ed30a18f30e4
rv_meas.RV_gmos[not.(ismissing.(rv_meas.RV_gmos))]

# ╔═╡ 30e3dc5b-3ce6-4dd7-9c2a-c82774909a8c
sum(not.(ismissing.(memb_stars.RV_gmos)))

# ╔═╡ d688d2e5-faca-4b14-801a-d58b08fd6654
let
	fig, ax = FigAxis()

	filt = not.(ismissing.(memb_stars.RV))
	p = scatter!(memb_stars.ra[filt], memb_stars.dec[filt], color=memb_stars.RV[filt],
		colorrange=(90, 130),
		colormap=:bluesreds
	)

	#Colorbar(fig[1, 2], p)
	fig
end

# ╔═╡ 7178e5b9-cc42-4933-970a-4707ba69dbe9
let
	fig, ax = FigAxis()

	p = arrows!(memb_stars.ra, memb_stars.dec, memb_stars.pmra, memb_stars.pmdec, 				
		color=memb_stars.RV,
		colorrange=(90, 130),
		colormap=:bluesreds,
		lengthscale=0.1
	)

	#Colorbar(fig[1, 2], p)
	fig
end


# ╔═╡ Cell order:
# ╟─6f2359d8-332b-11ef-0db9-f1f06474c561
# ╠═fcf50b2c-5703-4101-b3a7-45e5fb8e072f
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═9070c811-550c-4c49-9c58-0943b0f808b2
# ╠═e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
# ╠═dd29be70-4918-47e4-98cf-3608df14e88a
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═49e728ee-e545-4ee5-bd9f-b0a38c4eaf15
# ╠═da5a3b57-72bc-46e1-b1a0-6c02eb101626
# ╠═3d8b3c55-153b-4a4a-9289-78f34df07abc
# ╠═d4123ffd-cb32-486e-a93d-f48d7112a831
# ╠═248c2d5f-85cc-44be-bc63-43d4af470182
# ╠═7a50a176-96a5-4098-88d6-0fa2874d0f90
# ╟─7db47590-b82b-4822-8a01-eaff96c9389f
# ╠═23402e30-2064-494c-bd4a-8243cb474b61
# ╠═baacb491-9dff-4ff4-b62c-5842d79794da
# ╟─c470c2b9-093d-42ab-b96d-9dac231ccabc
# ╠═bb7f6769-ec92-460d-8423-449029175f79
# ╠═5f71b8e9-5540-4431-8028-4ce14c8d7856
# ╠═fcc95615-f24a-40fa-a7d1-9e474df9f798
# ╠═4d63430e-f59c-4c68-97e1-7eaa5679e55f
# ╠═1bc2541c-69ac-40aa-94f5-743707ca5a66
# ╟─27063de7-01d4-48a8-a06f-cc24aec662d2
# ╠═48c62811-136f-4962-a42c-b1dd1fc74f8c
# ╠═2e7ce524-573e-45c9-b0f4-ce9fea68e026
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
# ╠═1510b6de-09ae-474e-91c4-eea4e2cacdce
# ╠═bd51fe42-7e39-4cd8-8065-58ab9814f966
# ╠═ab49efb3-2ab7-47d0-a3a5-a342c789aa9b
# ╠═de762a39-b430-4452-ba87-8b8cf1ad9852
# ╟─89552b84-d12e-4d88-a58e-8b89ad4b2569
# ╠═e6f2de3b-ce32-4d61-851f-4e42fcce95c0
# ╠═5b2fceff-9c3e-472d-9310-31e920137e41
# ╠═7091dc6b-dd77-4f92-bd35-def8c7384f00
# ╠═4b305b83-1a3b-48a6-b19f-6f3ebed0768f
# ╠═c218cfa8-2f55-4957-bcdd-8b3970fe639a
# ╠═f7213298-3dcb-49ac-a9f0-a129a03423aa
# ╠═f7ec8bba-9f45-435b-b67c-33182e992dfd
# ╠═93d185f2-1eaa-4d35-87dd-b84f385483de
# ╠═36d2d86f-6a75-46f4-b48f-36137e95e90d
# ╠═9213cc36-74d1-452f-bd9a-eb5c5cab1f87
# ╠═fae534ec-78a9-4f40-91d0-8511b0bab399
# ╠═0df388a8-838b-46e6-a281-c48197459e37
# ╠═fcb7c823-4584-4766-9709-1661aa8f87a7
# ╠═f9e1bcbd-2def-4e64-b69e-d1cfcf7ffedd
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
# ╠═9e63c171-394e-4c05-b237-189fba0274e2
# ╠═c94b959d-2490-4a16-9ab7-c5c580316100
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╠═a1938588-ca40-4844-ab82-88c4254c435b
# ╠═3655b6a0-ac9a-4a49-86fd-6a6484409819
# ╠═7f13339e-6a0c-4944-8a36-5ab136fd8415
# ╠═aaaf5ba2-c9ed-41ec-a22a-d78ed96fd84e
# ╠═a2370516-e7ec-4502-ae4b-b111bcf68d36
# ╠═78e7aff0-3658-4d64-b117-58189a85307a
# ╠═f4f9dd06-1a1a-458b-be75-05d52623580c
# ╠═74152829-27ef-4d8d-8b32-ed30a18f30e4
# ╠═6b84c679-0532-465a-97dd-d62671077b61
# ╠═6734991c-16c0-4424-a2bb-84bfa811121f
# ╠═1b97d0e5-7a77-44a5-b609-ed8945cd959c
# ╠═abbd2a53-e077-4af7-a168-b571e1a906b8
# ╠═5cf336f6-e3eb-4668-b074-18b396f027be
# ╠═55504311-fbbd-4351-b9e2-ecab053ba207
# ╠═ccbb694d-f973-40fd-bab7-a2aefbd9fb0b
# ╠═74ad07df-f15c-436f-b390-ce95b27f7fab
# ╟─95ef12d2-174a-47c3-a106-9b4005b2a80d
# ╠═9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
# ╠═9380d9d1-58b2-432d-8528-d247cf5724e9
# ╠═318d29b9-4c84-4d38-af6d-048518952970
# ╠═b18e4622-41e0-4700-9e4b-3dbebeefea53
# ╠═bc7bd936-3f62-4430-8acd-8331ca3ee5ad
# ╠═764b5306-20f9-4810-8188-1bdf9482260f
# ╠═61e15c47-c454-48db-95f9-02abe052676e
# ╠═d5938fc3-9c8a-4e33-8401-500b4201df14
# ╠═6514815e-e1f3-44a1-8c0b-a90f13e0077e
# ╠═b654f5d5-ff48-4bae-bd69-a9db2673282b
# ╠═cfaafe3a-54bb-4e16-8b02-ef05fe7e4431
# ╟─5bee3041-fe7f-4f67-a73d-fb60ed006959
# ╠═ebdf8cfc-ebff-4331-bab2-998a734f2cd1
# ╠═ec4bc80f-22e6-49f9-a589-5c5bc9e50a8b
# ╠═eaf47845-95dd-4f7d-bf76-3d9b1154711a
# ╟─b5faa873-d5c9-447b-b90e-f693db26e6c2
# ╠═46434fa6-09df-4b02-9270-cbdcc9648e38
# ╠═dc4c0453-fed0-4733-9f59-0b2df506b45c
# ╠═7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═c2735c49-2892-46ac-bcf8-7cdcef409f44
# ╠═8b21cc49-ca17-4844-8238-e27e9752bee7
# ╠═c50f68d7-74c3-4c36-90c5-a5262982ed9f
# ╠═30e3dc5b-3ce6-4dd7-9c2a-c82774909a8c
# ╠═33f54afc-cdb9-4eb8-887f-5a43281b837c
# ╠═86776e68-d47f-43ed-b37f-432c864050bb
# ╠═1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
# ╠═b7ebd916-bffc-4ffc-a7f7-44ee315e2528
# ╠═0f5a9d9e-c5ca-4eb6-a0d2-5bb39b81daf6
# ╠═614f3f09-1880-491d-b41e-4e229330d66f
# ╠═319bd778-7e17-4bd7-856f-d6785b287219
# ╠═82a0e58a-30a4-4e42-b9c1-cb184eb551aa
# ╠═3a69f395-3c2d-4357-89af-5963d5fa79b8
# ╠═db78b22f-bee8-4d36-994e-6c7f9e96c9f2
# ╠═12ad8e67-0c7b-4029-91f6-5a2453ec799a
# ╠═e976d103-9956-47f3-adb0-7c449214b9d9
# ╠═9b4a0a1f-4b0c-4c90-b871-2bd244f0a908
# ╠═d688d2e5-faca-4b14-801a-d58b08fd6654
# ╠═7178e5b9-cc42-4933-970a-4707ba69dbe9
