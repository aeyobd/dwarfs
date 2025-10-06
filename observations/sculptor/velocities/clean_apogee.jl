### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	
	using CSV, DataFrames
	using Arya
	using CairoMakie

	import LilGuys as lguys
end

# ╔═╡ 05c93bf4-0f94-44cc-9aab-7c6193831806
using PyFITS

# ╔═╡ 2f5ca80c-d98d-45cc-a661-834002b620f6
using PairPlots

# ╔═╡ 7ed5bcf5-dfc3-4e79-a608-d503124a1e96
using LilGuys

# ╔═╡ e67ec69b-cfd1-4da7-8f0a-3c1a30f6f48c
using Turing

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Apogee RV analsysi.

This notebook analyzes the RV stars in APOGEE and creates a sample sutable to be combined with others for velocity analysis.

"""

# ╔═╡ d3c0f2b5-7812-4187-b2e1-b0f37114bd40
md"""
Creates:
- `processed/rv_apogee.fits`
- `mcmc_summary.---.csv`
Depends on:
- `processed/rv_apogee_xmatch.fits` (`apogee_xmatch.jl`)
"""

# ╔═╡ 7f3677f2-36e7-470e-92eb-5a33a02e53f6
md"""
xmatch by source_id
note that we use 
- `apogee_all.RV_sigma .< 3*apogee_all.RV_err`
(more restrictive...).
"""

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, kurtosis, sem

# ╔═╡ ef1bb6f5-00b8-405b-8702-244eca618644
import DensityEstimators: histogram, calc_limits, make_bins

# ╔═╡ 36634dea-21bc-4823-8a15-7bce20b6fc17
import TOML

# ╔═╡ 9a20ce08-79ce-4e23-ac4d-1c0af8de6ea7
module RVUtils
	include("../../rv_utils.jl")
	not = !
end

# ╔═╡ 257caa92-ccce-44ff-88ee-6a9c550ae5d9
CairoMakie.activate!(type=:png)

# ╔═╡ dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
not = !

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "../data/"

# ╔═╡ 77e7884c-0360-4b7f-b9ad-e81be2516552
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 7a50a176-96a5-4098-88d6-0fa2874d0f90
j24 = read_fits("..//data/jensen+24_wide_2c.fits")

# ╔═╡ d4d0a488-1c0e-4bc6-9a88-94b27d84e8ce
apogee_raw = read_fits("processed/apogee_xmatch.fits")

# ╔═╡ c470c2b9-093d-42ab-b96d-9dac231ccabc
md"""
## Apogee
"""

# ╔═╡ 45d6b4aa-ca44-4a71-afb8-ba6b2e674c7a
md"""
APOGEE DR 17 sample from federico's paper is apogee_raw_2. I have also xmatched the field against apogee which I use for self-consistency.
"""

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

	apogee_all = apogee_all[:, [:source_id, :ra, :dec, :RV, :RV_err, :RV_sigma, :RV_count, :fe_h, :fe_h_err, :STARFLAG, :RV_FLAG, :SNR]]
	apogee_all[ismissing.(apogee_all.fe_h), :fe_h] .= NaN
	apogee_all[ismissing.(apogee_all.fe_h_err), :fe_h_err] .= NaN

	apogee_all
end

# ╔═╡ 961e230d-e549-49b5-a98c-aaa330fa42da
@assert 0 == sum(apogee_all.RV_FLAG .> 0) # no need to include more flags

# ╔═╡ aa66b278-c408-49f4-b7e5-232341d98c3e
F_best = RVUtils.get_f_best.([j24], apogee_all.source_id)

# ╔═╡ 9c7aff9c-80d3-463d-9d38-585daf1feba6
p_chi2 = RVUtils.prob_chi2(apogee_all)

# ╔═╡ d1073372-b1bc-4c12-ab7f-e2bcd6ed47c3
F_scatter = RVUtils.filter_chi2(p_chi2)

# ╔═╡ 16559e75-f1b7-4078-a176-cb44e30942dd
F_match = .!ismissing.(apogee_all.source_id) .& (F_best .== 1.0)

# ╔═╡ 14094089-b9cd-464c-82aa-a41be8fc1bbf
df_out = let
	df = copy(apogee_all)
	df[!, :F_scatter] = F_scatter
	df[!, :F_match] = F_match
	df
end

# ╔═╡ f71152ad-d576-4205-bece-92c85783c089
write_fits("processed/rv_apogee.fits", df_out, overwrite=true)

# ╔═╡ bb5aefb9-7db5-4ae2-8f32-ef2a442e5c77
md"""
# Numbers
"""

# ╔═╡ f1742d6c-effc-465e-bbea-5c89ac575241
sum(F_match)

# ╔═╡ 013dcaf5-36f0-4bf1-bbd0-5facc465e30e
sum(F_scatter .& F_match)

# ╔═╡ 9c2f80c3-83ff-465b-9f1d-577661f1aa2f
sum(apogee_all.RV_count)

# ╔═╡ 59a20996-7b90-4afd-818a-fffd8029350f
md"""
# Plots and analysis
"""

# ╔═╡ bb57569e-0a4c-455e-b3c6-7bd2897bb843
sum(apogee_raw.RV_FLAG .> 0)

# ╔═╡ 0e14808d-df92-419b-b272-f3a08f4b86b1
apogee_raw_2 = CSV.read("$data_dir/sculptor_apogeeDR17_xmatch.csv", DataFrame)

# ╔═╡ bf088da8-a7c5-46e7-8e5f-ef5b91c10fb8
filt_apogee = RVUtils.sigma_clip(apogee_all.RV, 3) .& F_scatter .& F_match

# ╔═╡ 4d63430e-f59c-4c68-97e1-7eaa5679e55f
apogee = apogee_all[filt_apogee, :]

# ╔═╡ 5f71b8e9-5540-4431-8028-4ce14c8d7856
let 
	fig, ax = FigAxis(
		xlabel = "RV",
		ylabel = "count"
	)
	bins = make_bins(apogee_all.RV, calc_limits(apogee_all.RV), bandwidth=1)

	h = histogram(apogee.RV, bins, normalization=:none)
	lines!(h)
	
	h = histogram(apogee_all.RV, bins, normalization=:none)

	lines!(h)


	fig
end

# ╔═╡ 1bc2541c-69ac-40aa-94f5-743707ca5a66
hist(apogee.RV)

# ╔═╡ 63792b12-dde2-4361-88b2-1f5a789631ff
hist(filter(isfinite, apogee.fe_h))

# ╔═╡ 832e9183-13b5-4d67-8bec-b16af0ccd651
hist(filter(isfinite, (apogee.RV_sigma ./ apogee.RV_err) .^ 2 ./ sqrt.(apogee.RV_count)))

# ╔═╡ 5dafeaad-df9e-4558-99a7-721bdb91d28d
hist(apogee.RV_err)

# ╔═╡ 68960b6f-081b-46c1-8e27-8fda8a0521ee
hist(apogee.RV_err .* sqrt.(apogee.RV_count))

# ╔═╡ 949bdb6a-4554-4659-b0d3-7e669f470c10
hist(apogee.RV_sigma)

# ╔═╡ 4f230076-cad7-425a-8156-ca47e7f6c911
hist((apogee.RV_sigma) .^ 2 ./ (apogee.RV_err ) .^ 2)

# ╔═╡ 8bd988f3-884f-430b-b30c-3efc2ecbdd9c
hist(apogee.SNR)

# ╔═╡ 27603f1a-cce9-44c9-99e9-903cf375ee98
hist(filter(isfinite, p_chi2))

# ╔═╡ 2a514055-3f00-4dbc-86d4-882a0d0df77b
sum(p_chi2 .< 0.001)

# ╔═╡ 68cc4da9-b30c-4e51-9f96-b641dd89a7cf
sum(.!F_scatter)

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

# ╔═╡ 6009cc0e-6936-4589-85aa-dd2864948b7c
let
	fig = Figure()
    ax = Axis(fig[1,1])

	scatter!(apogee.ra, apogee.dec, markersize=5)
	scatter!(apogee_raw_2[:, "RA (deg)"], apogee_raw_2[:, "Dec (deg)"], markersize=3)

	fig
end

# ╔═╡ af916274-d081-45dd-ba62-626f7d409ebe
model = RVUtils.model_vel_1c(apogee.RV, apogee.RV_err)

# ╔═╡ cd34e5eb-c330-4476-b20a-fe74e67fe69c
chain = sample(model, NUTS(), MCMCThreads(), 1000, 16)

# ╔═╡ 68fd09a2-d2ca-44a1-b141-94fe9b968ce9
pairplot(chain)

# ╔═╡ 54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
df_summary = RVUtils.summarize(chain)

# ╔═╡ 5db61bab-a71e-4bcb-a706-f150f1126f4c
CSV.write("processed/mcmc_summary.rv_apogee.csv", df_summary)

# ╔═╡ ba337b2c-9a0c-40b3-867f-c9ebf683e1e3
apogee.RV

# ╔═╡ 0ae09680-cc34-492d-85bb-dec3ced936cf
samples = DataFrame(chain)

# ╔═╡ beb53bbc-0d2a-468d-a40c-2958df7e95ca
RVUtils.plot_samples(apogee, samples, bins=30)

# ╔═╡ 0ba4323f-8ffe-4ac7-999e-981f8c282708


# ╔═╡ Cell order:
# ╟─811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╟─d3c0f2b5-7812-4187-b2e1-b0f37114bd40
# ╟─7f3677f2-36e7-470e-92eb-5a33a02e53f6
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═05c93bf4-0f94-44cc-9aab-7c6193831806
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═2f5ca80c-d98d-45cc-a661-834002b620f6
# ╠═ef1bb6f5-00b8-405b-8702-244eca618644
# ╠═7ed5bcf5-dfc3-4e79-a608-d503124a1e96
# ╠═36634dea-21bc-4823-8a15-7bce20b6fc17
# ╠═e67ec69b-cfd1-4da7-8f0a-3c1a30f6f48c
# ╠═9a20ce08-79ce-4e23-ac4d-1c0af8de6ea7
# ╠═257caa92-ccce-44ff-88ee-6a9c550ae5d9
# ╠═dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═77e7884c-0360-4b7f-b9ad-e81be2516552
# ╠═7a50a176-96a5-4098-88d6-0fa2874d0f90
# ╠═d4d0a488-1c0e-4bc6-9a88-94b27d84e8ce
# ╟─c470c2b9-093d-42ab-b96d-9dac231ccabc
# ╟─45d6b4aa-ca44-4a71-afb8-ba6b2e674c7a
# ╠═bb7f6769-ec92-460d-8423-449029175f79
# ╠═961e230d-e549-49b5-a98c-aaa330fa42da
# ╠═d1073372-b1bc-4c12-ab7f-e2bcd6ed47c3
# ╠═aa66b278-c408-49f4-b7e5-232341d98c3e
# ╠═9c7aff9c-80d3-463d-9d38-585daf1feba6
# ╠═16559e75-f1b7-4078-a176-cb44e30942dd
# ╠═14094089-b9cd-464c-82aa-a41be8fc1bbf
# ╠═f71152ad-d576-4205-bece-92c85783c089
# ╟─bb5aefb9-7db5-4ae2-8f32-ef2a442e5c77
# ╠═f1742d6c-effc-465e-bbea-5c89ac575241
# ╠═013dcaf5-36f0-4bf1-bbd0-5facc465e30e
# ╠═9c2f80c3-83ff-465b-9f1d-577661f1aa2f
# ╟─59a20996-7b90-4afd-818a-fffd8029350f
# ╠═bb57569e-0a4c-455e-b3c6-7bd2897bb843
# ╠═0e14808d-df92-419b-b272-f3a08f4b86b1
# ╠═5f71b8e9-5540-4431-8028-4ce14c8d7856
# ╠═4d63430e-f59c-4c68-97e1-7eaa5679e55f
# ╠═bf088da8-a7c5-46e7-8e5f-ef5b91c10fb8
# ╠═1bc2541c-69ac-40aa-94f5-743707ca5a66
# ╠═63792b12-dde2-4361-88b2-1f5a789631ff
# ╠═832e9183-13b5-4d67-8bec-b16af0ccd651
# ╠═5dafeaad-df9e-4558-99a7-721bdb91d28d
# ╠═68960b6f-081b-46c1-8e27-8fda8a0521ee
# ╠═949bdb6a-4554-4659-b0d3-7e669f470c10
# ╠═4f230076-cad7-425a-8156-ca47e7f6c911
# ╠═8bd988f3-884f-430b-b30c-3efc2ecbdd9c
# ╠═27603f1a-cce9-44c9-99e9-903cf375ee98
# ╠═2a514055-3f00-4dbc-86d4-882a0d0df77b
# ╠═68cc4da9-b30c-4e51-9f96-b641dd89a7cf
# ╟─27063de7-01d4-48a8-a06f-cc24aec662d2
# ╠═48c62811-136f-4962-a42c-b1dd1fc74f8c
# ╠═6009cc0e-6936-4589-85aa-dd2864948b7c
# ╠═af916274-d081-45dd-ba62-626f7d409ebe
# ╠═cd34e5eb-c330-4476-b20a-fe74e67fe69c
# ╠═68fd09a2-d2ca-44a1-b141-94fe9b968ce9
# ╠═54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
# ╠═5db61bab-a71e-4bcb-a706-f150f1126f4c
# ╠═ba337b2c-9a0c-40b3-867f-c9ebf683e1e3
# ╠═0ae09680-cc34-492d-85bb-dec3ced936cf
# ╠═beb53bbc-0d2a-468d-a40c-2958df7e95ca
# ╠═0ba4323f-8ffe-4ac7-999e-981f8c282708
