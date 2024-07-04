### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	
	using CSV, DataFrames
	using Arya
	using GLMakie
	using FITSIO
	using Measurements

	import LilGuys as lguys
end

# ╔═╡ 9070c811-550c-4c49-9c58-0943b0f808b2
using Turing

# ╔═╡ e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
using PairPlots

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

# ╔═╡ 09a0c8e9-7d46-4167-894f-d616e1bf3f5c
function add_xi_eta_rell!(df)
end

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
function xmatch(df1::DataFrame, df2::DataFrame, angle=2)
	angle = angle / 3600
	dists = ang_dist.(df1.ra, df1.dec, Ref(df2.ra), Ref(df2.dec))
	filt = minimum.(dists) .< angle
	idxs = argmin.(dists)
	return filt, idxs
end

# ╔═╡ 7db47590-b82b-4822-8a01-eaff96c9389f
md"""
## Dart
"""

# ╔═╡ 23402e30-2064-494c-bd4a-8243cb474b61
begin 
	dart = CSV.read("../data/Sculptor_DARTS_BEST_psat40.csv", DataFrame)
	rename!(dart, 
		"Elliptical radius"=>"r_h",
		"vel"=>"RV",
		"evel"=>"RV_err"
	)
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
# GMOS
"""

# ╔═╡ b7345279-4f80-47ad-a726-537571849eae
begin 
	gmos = CSV.read("data/targets_stellar_met.csv", DataFrame)
	rename!(gmos,
		"RA"=>"ra",
		"Dec"=>"dec"
	)

	gmos_rv_sys_err = 13.3
	gmos.RV_err .= @. sqrt(gmos.RV_err^2 + gmos_rv_sys_err^2)
end

# ╔═╡ c31f9e07-c520-443b-94dc-787519021d01
md"""
## Tolstoy
"""

# ╔═╡ 15f2a8e2-90df-48a9-a7bf-e86955f566ce

begin 
	FITS("../data/tolstoy_2023.fits", "r") do f
		global tolstoy23 = DataFrame(f[2])
	end

	rename!(tolstoy23, 
		:Vlos => :RV,
		:e_Vlos => :RV_err,
	)
end

# ╔═╡ bf919560-1349-492f-8335-04348f2ace1b
hist(tolstoy23.RV)

# ╔═╡ 180fac98-678c-4a14-966c-385387c60ac3
md"""
## Walker
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

# ╔═╡ 5cf336f6-e3eb-4668-b074-18b396f027be
prior_samples = DataFrame(
	sample(normal_dist(rv_meas.RV, rv_meas.RV_err), Prior(), 10000)
)

# ╔═╡ 55504311-fbbd-4351-b9e2-ecab053ba207
function plot_samples!(samples, x;
		thin=10, color=:black, alpha=nothing, kwargs...)

	alpha = 1 / (size(samples, 1))^(1/3)
	for sample in eachrow(samples)[1:thin:end]
		y = lguys.gaussian.(x, sample.μ, sample.σ)
		lines!(x, y, color=color, alpha=alpha)
	end
end

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
function plot_sample_normal_fit(sample, props)
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = Arya.histogram(Float64.(sample.RV), 10, normalization=:pdf)
	
	lines!(h)

	μ, σ, _, _ = props
	x_model = LinRange(80, 140, 1000)
	y_model = lguys.gaussian.(x_model, μ, σ)
	lines!(x_model, y_model)

	fig
end

# ╔═╡ 318d29b9-4c84-4d38-af6d-048518952970
samples = DataFrame(sample(normal_dist(rv_meas.RV, rv_meas.RV_err), NUTS(0.65), 10000))

# ╔═╡ b18e4622-41e0-4700-9e4b-3dbebeefea53
describe(samples)

# ╔═╡ bc7bd936-3f62-4430-8acd-8331ca3ee5ad
pairplot(samples[:, [:μ, :σ]])

# ╔═╡ 764b5306-20f9-4810-8188-1bdf9482260f
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = Arya.histogram(rv_meas.RV, 25, normalization=:pdf)
	
	plot_samples!(samples, LinRange(80, 140, 100))
	lines!(h)

	fig
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

# ╔═╡ 0f5a9d9e-c5ca-4eb6-a0d2-5bb39b81daf6
σ_m = median(samples.σ)

# ╔═╡ 319bd778-7e17-4bd7-856f-d6785b287219
quantile(samples.σ, [0.16, 0.84]) .- σ_m

# ╔═╡ 82a0e58a-30a4-4e42-b9c1-cb184eb551aa
md"""
# Misc plots
"""

# ╔═╡ 3a69f395-3c2d-4357-89af-5963d5fa79b8
let
	fig, ax = FigAxis(
		xlabel=L"$\log r_\textrm{ell}$ / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)
	
	scatter!(log10.(rv_meas.r_ell), rv_meas.RV)


	fig
end

# ╔═╡ c50f68d7-74c3-4c36-90c5-a5262982ed9f
μs, σs, μ_errs, σ_errs = binned_mu_sigma(rv_meas.r_ell, rv_meas.RV, rv_meas.RV_err, bins)

# ╔═╡ 1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
bin_errs = diff(bins) / 2

# ╔═╡ 33f54afc-cdb9-4eb8-887f-5a43281b837c
let
	fig, ax = FigAxis(
		xlabel = "R / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)

	scatter!(rv_meas.r_ell, rv_meas.RV)
	errscatter!(midpoints(bins), μs, yerr=μ_errs, color=:black)
	
	errorbars!(midpoints(bins), μs .+ σs, bin_errs, direction = :x, color=:black)
	errorbars!(midpoints(bins), μs .- σs, bin_errs, direction = :x, color=:black)


	fig
end

# ╔═╡ b7ebd916-bffc-4ffc-a7f7-44ee315e2528
stephist(rv_meas.r_ell, bins=bins)

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

# ╔═╡ db78b22f-bee8-4d36-994e-6c7f9e96c9f2
value(x::Measurement) = x.val

# ╔═╡ 12ad8e67-0c7b-4029-91f6-5a2453ec799a
err(x::Measurement) = x.err

# ╔═╡ e976d103-9956-47f3-adb0-7c449214b9d9
import StatsBase as sb

# ╔═╡ 2f5e5451-e653-43e8-8d2b-aeff4c575889
"""Standard error of variance"""
function sev(x)
	k = sb.kurtosis(x) + 3
	s = sb.std(x)
	n = length(x)

	var_s2 = (k - (n-3)/(n-1)) * s^4 / n
	return 1 / 2s * sqrt(var_s2)
end

# ╔═╡ f44d3c45-064c-4bf3-bf2c-e1f0d3504dd1
"""Standard error of variance"""
function sev_iid(x)
	s = sb.std(x)
	n = length(x)
	return s / sqrt(2n - 2)
end

# ╔═╡ 0fca9a64-347e-48f0-84fb-e7413e98d508
sev(rv_meas.RV)

# ╔═╡ e4294754-0931-4a0b-b9fe-0fa11878ea21
sev_iid(rv_meas.RV)

# ╔═╡ 0663e4c2-d68d-423c-9eed-1f47ab555621
rv_meas.RV_err

# ╔═╡ a5e55f3e-de53-443c-9d28-e539aa6206a8
errscatter(dart.r_h, dart.feh, yerr=collect(zip(dart.feh .- dart[:, "feh_lo"], dart[:, "fe_hi"] .- dart.feh)))

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

# ╔═╡ bb7f6769-ec92-460d-8423-449029175f79
begin 
	apogee_all = CSV.read("../data/sculptor_apogeeDR17_xmatch.csv", DataFrame)

	rename!(apogee_all, 
		"Elliptical half-light radii"=>"r_h",
		"VHELIO_AVG"=>"RV",
		"VERR"=>"RV_err",
		"GAIAEDR3_PMRA"=>"pmra",
		"GAIAEDR3_PMDEC"=>"pmdec",
		"RA (deg)"=>"ra",
		"Dec (deg)"=>"dec"
	)

	_apogee_all_filt = not.(ismissing.(apogee_all.RV))

	apogee_all = DataFrame(apogee_all[_apogee_all_filt, :])

	apogee_all[!, :RV] = float.(apogee_all.RV)
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
		"Dec (deg)"=>"dec"
	)

	_apogee_filt = not.(ismissing.(apogee_notdart.RV))

	apogee_notdart = DataFrame(apogee_notdart[_apogee_filt, :])
end

# ╔═╡ 2e7ce524-573e-45c9-b0f4-ce9fea68e026
a2 = apogee_all[not.(xmatch(apogee_all, dart)[1]), :];

# ╔═╡ fcf50b2c-5703-4101-b3a7-45e5fb8e072f
size(a2, 1) + size(dart, 1) # should be 617

# ╔═╡ f8775eb1-2cb9-4d81-8c02-39c43ceb9b45
sort(a2.GAIAEDR3_SOURCE_ID) == sort(apogee_notdart.GAIAEDR3_SOURCE_ID)

# ╔═╡ 77830e77-50a4-48b7-b070-8fbd7508c173
begin 
	FITS("../data/walker_2009.fits", "r") do f
		global walker09_all = DataFrame(f[2])
	end

	rename!(walker09_all,
		:HV => :RV,
		:e_HV => :RV_err,
	)

	_walker_filt = walker09_all.Mmb .> 0.9
	_walker_filt .&= not.(isnan.(walker09_all.RV))
	walker09_all = walker09_all[_walker_filt, :]
end

# ╔═╡ 20d11209-3d1e-4cef-9364-a9f3f6ee938d
filt_walker = 70 .< walker09_all.RV .< 150

# ╔═╡ e1c2e02e-05de-4faf-af2f-f93759ddadfe
filt_walker_2 = sigma_clip(walker09_all[filt_walker, :].RV, 2.5)

# ╔═╡ a5a95eba-d282-4881-a84e-a25d4c83f114
walker09 = walker09_all[filt_walker, :][filt_walker_2, :]

# ╔═╡ 9d1495e8-5f8a-4892-b5b5-b22f3eb6ab7c
let 
	fig, ax = FigAxis(xlabel=xlabel)
	bins = Arya.make_bins(walker09_all.RV, Arya.calc_limits(walker09_all.RV), bandwidth=1)


	
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
props = Dict(name => fit_rv_sigma(data.RV, data.RV_err) for (name, data) in datasets)

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
		ylabel=L"$\sigma_{v}$ / km s${^-1}$"
	)


	errscatter!(x, y, yerr=yerr)
	fig
end

# ╔═╡ eaf47845-95dd-4f7d-bf76-3d9b1154711a
for key in names(props)
	println(key)
	println(props[key])
	@info plot_sample_normal_fit(datasets[key], props[key])
end

# ╔═╡ 1e06b622-ef1c-43ed-9453-bceab92c1889
let 
	#global rv_meas = vcat(dart, apogee, gmos, cols=:union)

	filt = not.(ismissing.(rv_meas.RV))
	filt .&= rv_meas.RV .< 200
	filt .&= rv_meas.RV .> 50
	rv_meas = rv_meas[filt, :]

	rv_meas[!, :RV] = float.(rv_meas.RV)

	xi, eta = lguys.to_tangent(rv_meas.ra, rv_meas.dec, ra0, dec0)

	rv_meas[!, :xi] = xi
	rv_meas[!, :eta] = eta

	r_ell = 60lguys.calc_r_ell(xi, eta, ell, PA)

	rv_meas[!, :r_ell] = r_ell

	rv_meas

end

# ╔═╡ d688d2e5-faca-4b14-801a-d58b08fd6654
let
	fig, ax = FigAxis()

	filt = not.(ismissing.(rv_meas.RV))
	p = scatter!(rv_meas.ra[filt], rv_meas.dec[filt], color=rv_meas.RV[filt],
		colorrange=(90, 130),
		colormap=:bluesreds
	)

	#Colorbar(fig[1, 2], p)
	fig
end

# ╔═╡ 7178e5b9-cc42-4933-970a-4707ba69dbe9
let
	fig, ax = FigAxis()

	p = arrows!(rv_meas.ra, rv_meas.dec, rv_meas.pmra, rv_meas.pmdec, 				
		color=rv_meas.RV,
		colorrange=(90, 130),
		colormap=:bluesreds,
		lengthscale=0.1
	)

	#Colorbar(fig[1, 2], p)
	fig
end

# ╔═╡ 75442700-8532-4f86-8469-555c162edb97
Arya.std(rv_meas.RV)

# ╔═╡ 8b21cc49-ca17-4844-8238-e27e9752bee7
bins = Arya.bins_equal_number(rv_meas.r_ell, n=6)

# ╔═╡ 3741255c-d0f4-47c1-97dd-cc6823b547d9
# ╠═╡ disabled = true
#=╠═╡
bins = Arya.bins_equal_number(log10.(rv_meas.r_ell), n=10)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─6f2359d8-332b-11ef-0db9-f1f06474c561
# ╠═fcf50b2c-5703-4101-b3a7-45e5fb8e072f
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═9070c811-550c-4c49-9c58-0943b0f808b2
# ╠═e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
# ╠═dd29be70-4918-47e4-98cf-3608df14e88a
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═09a0c8e9-7d46-4167-894f-d616e1bf3f5c
# ╠═49e728ee-e545-4ee5-bd9f-b0a38c4eaf15
# ╠═da5a3b57-72bc-46e1-b1a0-6c02eb101626
# ╠═3d8b3c55-153b-4a4a-9289-78f34df07abc
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
# ╠═b6de4afb-9362-4618-a114-df460031e4f9
# ╠═b7345279-4f80-47ad-a726-537571849eae
# ╠═c31f9e07-c520-443b-94dc-787519021d01
# ╠═15f2a8e2-90df-48a9-a7bf-e86955f566ce
# ╠═bf919560-1349-492f-8335-04348f2ace1b
# ╠═180fac98-678c-4a14-966c-385387c60ac3
# ╠═77830e77-50a4-48b7-b070-8fbd7508c173
# ╠═9d1495e8-5f8a-4892-b5b5-b22f3eb6ab7c
# ╠═20d11209-3d1e-4cef-9364-a9f3f6ee938d
# ╠═e1c2e02e-05de-4faf-af2f-f93759ddadfe
# ╠═a5a95eba-d282-4881-a84e-a25d4c83f114
# ╠═8c59d607-d484-497c-9ee7-751a4edd4992
# ╠═1e06b622-ef1c-43ed-9453-bceab92c1889
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
# ╠═33f54afc-cdb9-4eb8-887f-5a43281b837c
# ╠═1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
# ╠═b7ebd916-bffc-4ffc-a7f7-44ee315e2528
# ╠═0f5a9d9e-c5ca-4eb6-a0d2-5bb39b81daf6
# ╠═614f3f09-1880-491d-b41e-4e229330d66f
# ╠═319bd778-7e17-4bd7-856f-d6785b287219
# ╠═82a0e58a-30a4-4e42-b9c1-cb184eb551aa
# ╠═3a69f395-3c2d-4357-89af-5963d5fa79b8
# ╠═3741255c-d0f4-47c1-97dd-cc6823b547d9
# ╠═db78b22f-bee8-4d36-994e-6c7f9e96c9f2
# ╠═12ad8e67-0c7b-4029-91f6-5a2453ec799a
# ╠═e976d103-9956-47f3-adb0-7c449214b9d9
# ╠═2f5e5451-e653-43e8-8d2b-aeff4c575889
# ╠═f44d3c45-064c-4bf3-bf2c-e1f0d3504dd1
# ╠═0fca9a64-347e-48f0-84fb-e7413e98d508
# ╠═e4294754-0931-4a0b-b9fe-0fa11878ea21
# ╠═0663e4c2-d68d-423c-9eed-1f47ab555621
# ╠═a5e55f3e-de53-443c-9d28-e539aa6206a8
# ╠═9b4a0a1f-4b0c-4c90-b871-2bd244f0a908
# ╠═d688d2e5-faca-4b14-801a-d58b08fd6654
# ╠═7178e5b9-cc42-4933-970a-4707ba69dbe9
# ╠═75442700-8532-4f86-8469-555c162edb97
