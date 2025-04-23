### A Pluto.jl notebook ###
# v0.20.6

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

# ╔═╡ 10120b1d-e5d1-4a07-b25f-35c79b85b3e9
using PyFITS

# ╔═╡ 2f5ca80c-d98d-45cc-a661-834002b620f6
using PairPlots

# ╔═╡ 7ed5bcf5-dfc3-4e79-a608-d503124a1e96
using LilGuys

# ╔═╡ e67ec69b-cfd1-4da7-8f0a-3c1a30f6f48c
using Turing

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Walker analy

This notebook analyzes the RV stars in walker+09 and creates a sample sutable to be combined with others for velocity analysis.

Data downloaded from CDS associated with observation paper.
Should have 1365 scl members.
Weighted mean of measurements for repeated stars.

Note that we combine two dataframes, one with individual measurements and a second using weighted sums. However, only the second one contains the stellar coordinates, so we crossmatch them.

"""

# ╔═╡ f2d0af5f-677a-4135-b7af-a6c3db43fdd4
md"""
about 11 stars are not xmatched well to gaia. We exclude these for simplicity.
"""

# ╔═╡ 0f370153-2fdc-4725-8d3f-63c1fa925b07
md"""
Creates:
- `processed/rv_walker+09.fits`
- ---`.csv`
Depends on:
- `../data/walker+09_observations.fits`
- `../data/walker+09_summary.fits`
- `../data/jensen+24_wide.fits`
- `../observed_properties.toml`

Assumes
- xmatch radius 3arcseconds
- Adds columns for RV sigma and count, and adds RV error for single star measurments.
"""

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, kurtosis, sem

# ╔═╡ 257caa92-ccce-44ff-88ee-6a9c550ae5d9
CairoMakie.activate!(type=:png)

# ╔═╡ dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
not = !

# ╔═╡ ef1bb6f5-00b8-405b-8702-244eca618644
import DensityEstimators: histogram, calc_limits, make_bins

# ╔═╡ 36634dea-21bc-4823-8a15-7bce20b6fc17
import TOML

# ╔═╡ 9a20ce08-79ce-4e23-ac4d-1c0af8de6ea7
module RVUtils
	include("../../rv_utils.jl")
	not = !
end

# ╔═╡ 79d7ebd9-2bbf-4d7a-a9db-c572372d53e8
⊕ = RVUtils.:⊕

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "../data/"

# ╔═╡ 77e7884c-0360-4b7f-b9ad-e81be2516552
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 7a50a176-96a5-4098-88d6-0fa2874d0f90
j24 = read_fits("../data/jensen+24_wide.fits")

# ╔═╡ ea398623-ebc2-420a-8377-65bef43c0154
begin 
	walker09_single = read_fits("$data_dir/walker+09_observations.fits") 
	
	walker09_single[:, "Galaxy"] = [s[1:3] for s in walker09_single.Target];

	walker09_single
end

# ╔═╡ 0e14808d-df92-419b-b272-f3a08f4b86b1
begin 
	walker09_averaged = read_fits("$data_dir/walker+09_summary.fits")

	rename!(walker09_averaged, 
		"RAJ2000"=>"ra",
		"DEJ2000" => "dec",
		)
	walker09_averaged[:, "Galaxy"] = [s[1:3] for s in walker09_averaged.Target]

	walker09_averaged
end

# ╔═╡ f8dd5663-910c-40c1-bd77-57621e2ad6c3
sum(walker09_averaged.Galaxy .== "Scl")

# ╔═╡ d4d0a488-1c0e-4bc6-9a88-94b27d84e8ce
begin
	walker09_all = copy(walker09_averaged)


	_walker_filt = walker09_all.Galaxy .== "Scl"
	walker09_all = walker09_all[_walker_filt, :]

	rename!(walker09_all, 
		:__HV_ => :RV,
		:e__HV_ => :RV_err
		   )
	
	walker09_all[:, :RV_count] .= 0
	walker09_all[:, :RV_sigma] .= 0.

	for (i, row) in enumerate(eachrow(walker09_all))
		xs = walker09_single[walker09_single.Target .== row.Target, :HV]
		x_errs = walker09_single[walker09_single.Target .== row.Target, :e_HV]

		ws = @. 1/x_errs^2
		walker09_all[i, :RV_count] = length(xs)
		walker09_all[i, :RV_sigma] = lguys.std(xs, ws)
		if ismissing(walker09_all.RV_err[i])
			walker09_all[i, :RV_err] = x_errs |> only
		end

		if length(xs) > 1 # will fail for single stars for some reason...
			@info "similar?", walker09_all.RV[i], lguys.mean(xs, ws), walker09_all.RV_err[i], RVUtils.sem_inv_var_weights(ws)
			@info "data:", xs, x_errs

			@assert isapprox(walker09_all.RV[i], lguys.mean(xs, ws), atol=max(walker09_all.RV_err[i], 0.1))
		end
	end

	j24_filt, j24_idx = RVUtils.xmatch(walker09_all, j24, 3) # catches 1 missed star

	walker09_all[!, :source_id] = j24.source_id[j24_idx]
	walker09_all[.!j24_filt, :source_id] .= missing

	walker09_all
end

# ╔═╡ 0ec29d83-a15e-4609-bfab-ba7366f034d2
sum(1 .- j24.F_BEST[j24_idx[j24_filt]])

# ╔═╡ a542dd17-0d01-4de9-83ee-1b3de933e781
F_match = j24_filt .& (j24.F_BEST[j24_idx] .== 1)

# ╔═╡ eb50e25e-c0a1-4c3c-af90-a89ffb2cff8b
p_chi2 =RVUtils.prob_chi2(walker09_all)

# ╔═╡ 034c138c-61f1-4e46-b312-6fcddc289721
F_scatter = @. isnan(p_chi2) || p_chi2 > 0.001

# ╔═╡ 39bd84c6-330b-4ddd-86f9-05136ba77316
(walker09_all.RV_sigma ./ walker09_all.RV_err)[.!F_scatter] # number of sigma for filtered stars

# ╔═╡ 0b9d4a5a-e448-4b73-86f1-4c73a02f8135
df_out = let
	df = copy(walker09_all)
	df[!, :F_scatter] = F_scatter
	df[!, :F_match] = F_match
	df
end

# ╔═╡ f71152ad-d576-4205-bece-92c85783c089
write_fits("processed/rv_walker+09.fits", df_out, overwrite=true)

# ╔═╡ 2c45110e-d323-4ef3-94c1-3e9922caa3c2
md"""
# Numbers
"""

# ╔═╡ b1120dfd-83f5-4b5a-a6b0-0fbd501fbbdb
sum(F_match)

# ╔═╡ fba308a6-b709-4037-963d-be75edd22c8b
sum(F_scatter)

# ╔═╡ 04cc2fc7-1550-4a36-9632-c7a14fd412d9
sum(F_scatter .& F_match)

# ╔═╡ d79ae7a9-7811-402d-a0c1-bce66ddf860a
sum(walker09_single.Galaxy .== "Scl")

# ╔═╡ a1601fa4-4851-41d0-87bd-6809c89a8faa
sum(walker09_all.RV_count)

# ╔═╡ a4d01dff-e5fe-49e2-90bb-f4929ede1fbd
sum(walker09_averaged.Galaxy .== "Scl"), length(walker09_all.Galaxy)

# ╔═╡ 66c42d45-0ac3-4613-a60c-99d3038d95d2
md"""
# Uncertainties
"""

# ╔═╡ 87a55d42-3000-4acf-a5de-087be9034da4
F_scatter_old = walker09_all.RV_sigma .< 3walker09_all.RV_err .* sqrt.(walker09_all.RV_count)

# ╔═╡ f9a69065-b4ff-4cbc-898a-8501dd3ff3ad
hist(p_chi2[.!F_scatter])

# ╔═╡ 8a0f11fd-86f9-4bab-9f5b-e5a8309d26d2
v_mean = leftjoin(walker09_single, walker09_averaged, on="Target", makeunique=true).__HV_

# ╔═╡ 5d447127-483f-4d2a-a0b2-8211aff7867e
v_mean_err = leftjoin(walker09_single, walker09_averaged, on="Target", makeunique=true).e__HV_

# ╔═╡ 393574f2-3c81-4815-80cd-b02de200d48f
filt_chi2 = walker09_all.RV_count .> 1

# ╔═╡ 78731218-3130-4445-b238-73d12f949fb2
Δv = @. walker09_all.RV_sigma ./ walker09_all.RV_err

# ╔═╡ e21b95df-92d1-4b25-b249-b06098f9375d
chi2 = Δv .^ 2

# ╔═╡ 24066b44-227f-49de-8375-c65929506889
md"""
If uncertanties are representative, then the plot below should be uniform with excess at zero (binarity).
"""

# ╔═╡ 0ec06d3c-ee56-4d91-978c-d29d2656da41
[1 - cdf(Chisq(n-1), 2n) for n in 2:30]

# ╔═╡ c9578ca4-dfd0-4ba8-b295-c3f25c41f22b
hist(filter(isfinite, p_chi2))

# ╔═╡ 7a959d88-95a0-4330-bdea-68df2574a902
md"""
## MCMC
"""

# ╔═╡ 3f31c73e-6953-498a-9de6-3f927f0dae87
md"""
## Exploratory
"""

# ╔═╡ 68ec5d7a-9665-405a-83eb-f2d49dcc254a
sum(walker09_single.Galaxy .== "Scl")

# ╔═╡ ee791fca-54f0-46b7-9589-c4984ad5e490
sum(walker09_averaged.Galaxy .== "Scl")

# ╔═╡ e8149a93-2069-4aac-8dc8-68a81c36c4ee
sum(walker09_all.RV_count)

# ╔═╡ 5a133aa4-cad9-4767-8dbd-ee37efc9319a
length(unique(walker09_all.Target)) == size(walker09_all, 1)

# ╔═╡ 656c34f6-ba9a-489b-9765-a27a8590a3e4
hist(walker09_all.Mmb)

# ╔═╡ 8faa8822-26a0-4425-95a3-4c99d535adb4
hist(walker09_all.RV_count)

# ╔═╡ d7f13c9e-67d7-412d-ae4b-d9064507a295
hist(filter(isfinite, walker09_all.RV_sigma))

# ╔═╡ d0813914-d8ee-4173-bb40-2e1c0e6fe6c4
hist(collect(skipmissing(walker09_all.RV_sigma ./ walker09_all.RV_err)))

# ╔═╡ 4098da43-4746-4f11-84b6-33ca5eb7a2c4
hist(collect(skipmissing(walker09_all.RV_err)))

# ╔═╡ df588ebc-80a1-4922-a907-a83750a8d24f
filt_walker = (walker09_all.Mmb .>= 0.5) .& (.!ismissing.(walker09_all.RV_err))

# ╔═╡ cf7c7083-5315-4d20-91e6-fe7718cfa07c
walker09 = walker09_all[filt_walker, :]

# ╔═╡ af916274-d081-45dd-ba62-626f7d409ebe
model = RVUtils.model_vel_1c(disallowmissing(walker09.RV), disallowmissing(walker09.RV_err))

# ╔═╡ cd34e5eb-c330-4476-b20a-fe74e67fe69c
chain = sample(model, NUTS(), MCMCThreads(), 1000, 16)

# ╔═╡ 68fd09a2-d2ca-44a1-b141-94fe9b968ce9
pairplot(chain)

# ╔═╡ 54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
df_summary = RVUtils.summarize(chain)

# ╔═╡ 60f45b81-8beb-4eb0-8f55-db04c5368eee
CSV.write("processed/rv_walker+09.csv", df_summary)

# ╔═╡ 0ae09680-cc34-492d-85bb-dec3ced936cf
samples = DataFrame(chain)

# ╔═╡ beb53bbc-0d2a-468d-a40c-2958df7e95ca
RVUtils.plot_samples(walker09, samples, bins=30)

# ╔═╡ 5ee4f9d3-21da-4323-b51f-c48c10193188
let 
	fig, ax = FigAxis(xlabel="radial velocity / km/s")
	bins = make_bins(walker09_all.RV, calc_limits(walker09_all.RV), bandwidth=3)


	
	h = histogram(walker09_all.RV, bins, normalization=:none)

	lines!(h, label="all")
	
	h = histogram(walker09.RV, bins, normalization=:none)
	lines!(h, label="selected")

	axislegend()

	fig
end

# ╔═╡ 0d45969a-1c24-4742-99df-5a3e13612e5d
hist(walker09.RV)

# ╔═╡ 0b02c16f-6478-457e-8947-a43058308a49
md"""
## Missed xmatches
May not be in gaia, closeish to magnitude limit and only 2 have a change of membership...
"""

# ╔═╡ 0df5f178-d644-4f5d-8e9d-05db24e9d535
walker09_missing = walker09_all[.!j24_filt, :]

# ╔═╡ 4dd94c3f-1eba-42ee-afdc-1d4f0a594584
_, matches_missing = RVUtils.xmatch(walker09_missing, j24, 3)

# ╔═╡ cd768c3d-7927-4b1c-9f67-90a2f1fc9ae5
let 
	fig = Figure()
	ax = Axis(fig[1,1])
		
	
	scatter!(walker09_missing.ra, walker09_missing.dec, color=Float64.(walker09_missing.Vmag))

	scatter!(j24.ra[matches_missing], j24.dec[matches_missing], markersize=2, color=COLORS[3])

	fig
end

# ╔═╡ 1cabc5ab-5088-455a-9f9e-93ae0ab94df3
let 
	fig = Figure()
	ax = Axis(fig[1,1])
		
	
	scatter!(walker09_missing.Vmag, j24.phot_g_mean_mag[matches_missing], color=walker09_missing.ra)

	fig
end

# ╔═╡ Cell order:
# ╟─811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╟─f2d0af5f-677a-4135-b7af-a6c3db43fdd4
# ╟─0f370153-2fdc-4725-8d3f-63c1fa925b07
# ╠═10120b1d-e5d1-4a07-b25f-35c79b85b3e9
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═2f5ca80c-d98d-45cc-a661-834002b620f6
# ╠═257caa92-ccce-44ff-88ee-6a9c550ae5d9
# ╠═dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
# ╠═ef1bb6f5-00b8-405b-8702-244eca618644
# ╠═7ed5bcf5-dfc3-4e79-a608-d503124a1e96
# ╠═36634dea-21bc-4823-8a15-7bce20b6fc17
# ╠═e67ec69b-cfd1-4da7-8f0a-3c1a30f6f48c
# ╠═9a20ce08-79ce-4e23-ac4d-1c0af8de6ea7
# ╠═79d7ebd9-2bbf-4d7a-a9db-c572372d53e8
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═f8dd5663-910c-40c1-bd77-57621e2ad6c3
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═77e7884c-0360-4b7f-b9ad-e81be2516552
# ╠═7a50a176-96a5-4098-88d6-0fa2874d0f90
# ╠═ea398623-ebc2-420a-8377-65bef43c0154
# ╠═0e14808d-df92-419b-b272-f3a08f4b86b1
# ╠═d4d0a488-1c0e-4bc6-9a88-94b27d84e8ce
# ╠═0ec29d83-a15e-4609-bfab-ba7366f034d2
# ╠═a542dd17-0d01-4de9-83ee-1b3de933e781
# ╠═eb50e25e-c0a1-4c3c-af90-a89ffb2cff8b
# ╠═034c138c-61f1-4e46-b312-6fcddc289721
# ╠═39bd84c6-330b-4ddd-86f9-05136ba77316
# ╠═0b9d4a5a-e448-4b73-86f1-4c73a02f8135
# ╠═f71152ad-d576-4205-bece-92c85783c089
# ╟─2c45110e-d323-4ef3-94c1-3e9922caa3c2
# ╠═b1120dfd-83f5-4b5a-a6b0-0fbd501fbbdb
# ╠═fba308a6-b709-4037-963d-be75edd22c8b
# ╠═04cc2fc7-1550-4a36-9632-c7a14fd412d9
# ╠═d79ae7a9-7811-402d-a0c1-bce66ddf860a
# ╠═a1601fa4-4851-41d0-87bd-6809c89a8faa
# ╠═a4d01dff-e5fe-49e2-90bb-f4929ede1fbd
# ╟─66c42d45-0ac3-4613-a60c-99d3038d95d2
# ╠═87a55d42-3000-4acf-a5de-087be9034da4
# ╠═f9a69065-b4ff-4cbc-898a-8501dd3ff3ad
# ╠═8a0f11fd-86f9-4bab-9f5b-e5a8309d26d2
# ╠═5d447127-483f-4d2a-a0b2-8211aff7867e
# ╠═393574f2-3c81-4815-80cd-b02de200d48f
# ╠═78731218-3130-4445-b238-73d12f949fb2
# ╠═e21b95df-92d1-4b25-b249-b06098f9375d
# ╟─24066b44-227f-49de-8375-c65929506889
# ╠═0ec06d3c-ee56-4d91-978c-d29d2656da41
# ╠═c9578ca4-dfd0-4ba8-b295-c3f25c41f22b
# ╟─7a959d88-95a0-4330-bdea-68df2574a902
# ╠═af916274-d081-45dd-ba62-626f7d409ebe
# ╠═cd34e5eb-c330-4476-b20a-fe74e67fe69c
# ╠═68fd09a2-d2ca-44a1-b141-94fe9b968ce9
# ╠═54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
# ╠═0ae09680-cc34-492d-85bb-dec3ced936cf
# ╠═beb53bbc-0d2a-468d-a40c-2958df7e95ca
# ╠═60f45b81-8beb-4eb0-8f55-db04c5368eee
# ╟─3f31c73e-6953-498a-9de6-3f927f0dae87
# ╠═68ec5d7a-9665-405a-83eb-f2d49dcc254a
# ╠═ee791fca-54f0-46b7-9589-c4984ad5e490
# ╠═e8149a93-2069-4aac-8dc8-68a81c36c4ee
# ╠═5a133aa4-cad9-4767-8dbd-ee37efc9319a
# ╠═656c34f6-ba9a-489b-9765-a27a8590a3e4
# ╠═8faa8822-26a0-4425-95a3-4c99d535adb4
# ╠═d7f13c9e-67d7-412d-ae4b-d9064507a295
# ╠═d0813914-d8ee-4173-bb40-2e1c0e6fe6c4
# ╠═4098da43-4746-4f11-84b6-33ca5eb7a2c4
# ╠═5ee4f9d3-21da-4323-b51f-c48c10193188
# ╠═df588ebc-80a1-4922-a907-a83750a8d24f
# ╠═cf7c7083-5315-4d20-91e6-fe7718cfa07c
# ╠═0d45969a-1c24-4742-99df-5a3e13612e5d
# ╟─0b02c16f-6478-457e-8947-a43058308a49
# ╠═0df5f178-d644-4f5d-8e9d-05db24e9d535
# ╠═4dd94c3f-1eba-42ee-afdc-1d4f0a594584
# ╠═cd768c3d-7927-4b1c-9f67-90a2f1fc9ae5
# ╠═1cabc5ab-5088-455a-9f9e-93ae0ab94df3
