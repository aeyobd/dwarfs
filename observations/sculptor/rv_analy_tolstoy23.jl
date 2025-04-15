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

# ╔═╡ 2f5ca80c-d98d-45cc-a661-834002b620f6
using PairPlots

# ╔═╡ 23b3766e-15b7-4b1e-bb41-d6af36a59caf
using PyFITS

# ╔═╡ 7ed5bcf5-dfc3-4e79-a608-d503124a1e96
using LilGuys

# ╔═╡ e67ec69b-cfd1-4da7-8f0a-3c1a30f6f48c
using Turing

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# TOLSTOY analy

This notebook analyzes the RV stars in tolstoy+23 and creates a sample sutable to be combined with others for velocity analysis.

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

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 9a20ce08-79ce-4e23-ac4d-1c0af8de6ea7
module RVUtils
	include("rv_utils.jl")
	not = !
end

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "data/"

# ╔═╡ 77e7884c-0360-4b7f-b9ad-e81be2516552
obs_properties = TOML.parsefile("observed_properties.toml")

# ╔═╡ c470c2b9-093d-42ab-b96d-9dac231ccabc
md"""
## Data loading
No need to filter further...
"""

# ╔═╡ bb7f6769-ec92-460d-8423-449029175f79

begin 
	tolstoy23_all = read_fits("$data_dir/tolstoy+23.fits")

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

	tolstoy23_all.RV_err = @. max((1.448 - 0.021*tolstoy23_all.S_N)*tolstoy23_all.RV_err, 0.65) # correction from arroyo-polonio 2024
	
	tolstoy23_all = tolstoy23_all
end

# ╔═╡ f71152ad-d576-4205-bece-92c85783c089
write_fits("processed/rv_tolstoy+23.fits", tolstoy23_all, overwrite=true)

# ╔═╡ 33d34e7e-0c0a-4c49-b0c1-0116c8ef7a38
filt_tolstoy = RVUtils.sigma_clip(tolstoy23_all.RV, 3) .& (tolstoy23_all.Mem .== "m")

# ╔═╡ 344f3c29-0873-4180-9025-51aeaeb2c681
tolstoy23 = tolstoy23_all[filt_tolstoy, :]

# ╔═╡ 4290a51e-961b-4b38-b954-1a65e946080c
md"""
# Plots
"""

# ╔═╡ 5f71b8e9-5540-4431-8028-4ce14c8d7856
let 
	fig, ax = FigAxis()
	bins = make_bins(tolstoy23_all.RV, calc_limits(tolstoy23_all.RV), bandwidth=1)

	h = histogram(tolstoy23_all.RV, bins, normalization=:none)
	lines!(h)
	
	h = histogram(tolstoy23.RV, bins, normalization=:none)

	lines!(h)


	fig
end

# ╔═╡ d86097ad-2feb-46d8-9ceb-a2743a74a1a9
hist((collect ∘ skipmissing)(tolstoy23_all.fe_h),
	 axis = (; xlabel = "[Fe/H]")
	)

# ╔═╡ af916274-d081-45dd-ba62-626f7d409ebe
model = RVUtils.model_vel_1c(tolstoy23.RV, tolstoy23.RV_err)

# ╔═╡ cd34e5eb-c330-4476-b20a-fe74e67fe69c
chain = sample(model, NUTS(), MCMCThreads(), 1000, 16)

# ╔═╡ 68fd09a2-d2ca-44a1-b141-94fe9b968ce9
pairplot(chain)

# ╔═╡ 54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
df_summary = Turing.summarize(chain)

# ╔═╡ 0ae09680-cc34-492d-85bb-dec3ced936cf
samples = DataFrame(chain)

# ╔═╡ beb53bbc-0d2a-468d-a40c-2958df7e95ca
RVUtils.plot_samples(tolstoy23, samples, bins=30)

# ╔═╡ 60f45b81-8beb-4eb0-8f55-db04c5368eee
CSV.write("processed/rv_tolstoy+23.csv", df_summary)

# ╔═╡ 8b909559-c7e9-4fac-9d29-bf3b654393fa
md"""
# Dart
"""

# ╔═╡ de449b47-2040-4514-aac3-9fe4729d6213
md"""
DART sample given from Federico. 

The sample is used in Battaglia + 2008 (Kinematic status of Scl). I cannot find the whole dataset publically available. 

Measurements taken using VLT/FLAMES. 

The paper reports $v_{\rm hel, sys} = 110.6\pm0.5$ and a velocity dispersion of $\sigma=10.1 \pm 0.3$
"""

# ╔═╡ ea398623-ebc2-420a-8377-65bef43c0154
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

# ╔═╡ e23e7d8c-e1ec-48db-a96d-83798da42e10
let
	filt, idx = RVUtils.xmatch(dart, tolstoy23)
	@info "matched $(sum(filt)) out of $(length(filt))"

	df_matched = copy(dart[filt, :])
	df_matched[:, :RV_tolstoy] = tolstoy23[idx[filt], :RV]
	df_matched[:, :RV_err_tolstoy] = tolstoy23[idx[filt], :RV_err]
	
	rv1  = df_matched.RV
	rv1_err  = df_matched.RV_err
	
	rv2  = df_matched.RV_tolstoy
	rv2_err  = df_matched.RV_err_tolstoy

	filt = RVUtils.filt_missing(rv1; low=75, high=150) .& RVUtils.filt_missing(rv2;  low=75, high=150)


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

# ╔═╡ ca774514-a1ea-4af2-8be4-cfef0f009574
model_dart = RVUtils.model_vel_1c(dart.RV, dart.RV_err)

# ╔═╡ 603b18e1-55d3-4a85-b281-0ae2a20222ea
chain_dart = sample(model_dart, NUTS(), MCMCThreads(), 1000, 16)

# ╔═╡ 1cca058a-57b8-4be2-9c81-e67b651a38fc
pairplot(chain_dart)

# ╔═╡ 3524c534-7d9d-41e2-8802-0526878ad498
df_summary_dart = Turing.summarize(chain_dart)

# ╔═╡ 155ac7dc-7ec6-43a1-bea4-3c079ab0f19c
samples_dart = DataFrame(chain_dart)

# ╔═╡ e9e30f25-76aa-4bfc-8f30-edad90efd47c
RVUtils.plot_samples(dart, samples_dart, bins=30)

# ╔═╡ db4256f6-23df-49e5-8a3b-d81c03257b5a
CSV.write("processed/rv_dart.csv", df_summary_dart)

# ╔═╡ e5e7a9ae-bec0-462a-97f9-262902fd41f5
md"""
# Gradient
"""

# ╔═╡ 05701206-01cc-4664-97ce-bc61000eeae2
"""
Fits a normal (gaussian) distribution to 3d data with errors (to include spatial gradient)
"""
@model function normal_gradient(x, xerr, ξ, η; μ_min=90, μ_max=120)
	μ ~ Uniform(μ_min, μ_max)
	A ~ Normal(0, 0.1)
	B ~ Normal(0, 0.1)
	m = @. μ + A*ξ + B*η
	σ ~ LogNormal(2.5, 1) # approx 1 - 100 km / s, very broad but should cover all
	s = @. sqrt(σ^2 + xerr^2)

	x ~ MvNormal(m, s)
end

# ╔═╡ cce6035b-6d83-4c6c-b81d-22fb83d16844


# ╔═╡ 1bdd359f-c111-43a8-9997-8c005e076588
samples_gradient = sample(normal_gradient(tolstoy23.RV, tolstoy23.RV_err, tolstoy23.xi, tolstoy23.eta, μ_min=30, μ_max=100), 
						  NUTS(0.65), MCMCThreads(), 1000, 16)

# ╔═╡ 589ac5e3-c5a1-4283-98e6-306647876051
pairplot(samples_gradient)

# ╔═╡ a08071ac-efaf-47b5-8688-bb3f30431f41
sqrt((0.073/0.017)^2 + (0.043 / 0.026)^2)

# ╔═╡ 47381bef-2d12-4813-be4a-cadc72f564b7
df_gradient = DataFrame(samples_gradient)

# ╔═╡ e32c52d1-9b45-4038-a3d7-94f4288aaacc
scatter(df_gradient.A, df_gradient.B, alpha=0.1, markersize=1, 
		axis=(;
			  limits=(-0.2, 0.2, -0.2, 0.2), 
			  aspect=DataAspect(),
			  xlabel = "A",
			  ylabel = "B",
			 )
	   )

# ╔═╡ aa868d3b-286c-45fc-b2d3-8525ba2fab58
import KernelDensity

# ╔═╡ 8b86d694-4acd-4a74-b235-c25a8c8679f9
kde = KernelDensity.kde((df_gradient.A, df_gradient.B))

# ╔═╡ 062528d2-12c3-4408-bd3a-77806b50730c
kde(0,0)

# ╔═╡ 37636050-1247-4ceb-b75c-2d3aeb87b364
sum(df_gradient.A .> 0)

# ╔═╡ b6933c35-8cb4-4171-9c28-315e3f787264
median(atand.(df_gradient.B ./ df_gradient.A))

# ╔═╡ a55a40a9-86b8-4af9-b259-309cd92600fa
median(sqrt.(df_gradient.B .^ 2 .+ df_gradient.A .^ 2)) * 60

# ╔═╡ 7fd48a81-159c-44ec-b3e4-5752eb709c0a
quantile(sqrt.(df_gradient.B .^ 2 .+ df_gradient.A .^ 2), [0.16, 0.5, 0.84]) * 60

# ╔═╡ 8964662b-6655-4281-b8e1-426293266382
quantile(atand.(df_gradient.B ./ df_gradient.A), [0.16, 0.5, 0.84])

# ╔═╡ d32f32a1-8b4f-4267-9176-82d30d623bf7
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.radial_velocity_gsr), 30, normalization=:pdf)
	
	plot_samples!(DataFrame(samples_gradient), LinRange(40, 110, 100), thin=15)
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

	fig
end

# ╔═╡ Cell order:
# ╠═811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═2f5ca80c-d98d-45cc-a661-834002b620f6
# ╠═257caa92-ccce-44ff-88ee-6a9c550ae5d9
# ╠═dfa3ccb0-66f7-49b9-bc6d-55e3f41070fe
# ╠═ef1bb6f5-00b8-405b-8702-244eca618644
# ╠═23b3766e-15b7-4b1e-bb41-d6af36a59caf
# ╠═7ed5bcf5-dfc3-4e79-a608-d503124a1e96
# ╠═36634dea-21bc-4823-8a15-7bce20b6fc17
# ╠═e67ec69b-cfd1-4da7-8f0a-3c1a30f6f48c
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═9a20ce08-79ce-4e23-ac4d-1c0af8de6ea7
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═77e7884c-0360-4b7f-b9ad-e81be2516552
# ╟─c470c2b9-093d-42ab-b96d-9dac231ccabc
# ╠═bb7f6769-ec92-460d-8423-449029175f79
# ╠═f71152ad-d576-4205-bece-92c85783c089
# ╠═344f3c29-0873-4180-9025-51aeaeb2c681
# ╠═33d34e7e-0c0a-4c49-b0c1-0116c8ef7a38
# ╟─4290a51e-961b-4b38-b954-1a65e946080c
# ╠═5f71b8e9-5540-4431-8028-4ce14c8d7856
# ╠═d86097ad-2feb-46d8-9ceb-a2743a74a1a9
# ╠═af916274-d081-45dd-ba62-626f7d409ebe
# ╠═cd34e5eb-c330-4476-b20a-fe74e67fe69c
# ╠═68fd09a2-d2ca-44a1-b141-94fe9b968ce9
# ╠═54c4c8b3-f9e3-4b85-afd6-5a5646f5eb79
# ╠═0ae09680-cc34-492d-85bb-dec3ced936cf
# ╠═beb53bbc-0d2a-468d-a40c-2958df7e95ca
# ╠═60f45b81-8beb-4eb0-8f55-db04c5368eee
# ╟─8b909559-c7e9-4fac-9d29-bf3b654393fa
# ╟─de449b47-2040-4514-aac3-9fe4729d6213
# ╠═ea398623-ebc2-420a-8377-65bef43c0154
# ╠═e23e7d8c-e1ec-48db-a96d-83798da42e10
# ╠═ca774514-a1ea-4af2-8be4-cfef0f009574
# ╠═603b18e1-55d3-4a85-b281-0ae2a20222ea
# ╠═1cca058a-57b8-4be2-9c81-e67b651a38fc
# ╠═3524c534-7d9d-41e2-8802-0526878ad498
# ╠═155ac7dc-7ec6-43a1-bea4-3c079ab0f19c
# ╠═e9e30f25-76aa-4bfc-8f30-edad90efd47c
# ╠═db4256f6-23df-49e5-8a3b-d81c03257b5a
# ╠═e5e7a9ae-bec0-462a-97f9-262902fd41f5
# ╠═05701206-01cc-4664-97ce-bc61000eeae2
# ╠═cce6035b-6d83-4c6c-b81d-22fb83d16844
# ╠═1bdd359f-c111-43a8-9997-8c005e076588
# ╠═589ac5e3-c5a1-4283-98e6-306647876051
# ╠═a08071ac-efaf-47b5-8688-bb3f30431f41
# ╠═47381bef-2d12-4813-be4a-cadc72f564b7
# ╠═e32c52d1-9b45-4038-a3d7-94f4288aaacc
# ╠═aa868d3b-286c-45fc-b2d3-8525ba2fab58
# ╠═8b86d694-4acd-4a74-b235-c25a8c8679f9
# ╠═062528d2-12c3-4408-bd3a-77806b50730c
# ╠═37636050-1247-4ceb-b75c-2d3aeb87b364
# ╠═b6933c35-8cb4-4171-9c28-315e3f787264
# ╠═a55a40a9-86b8-4af9-b259-309cd92600fa
# ╠═7fd48a81-159c-44ec-b3e4-5752eb709c0a
# ╠═8964662b-6655-4281-b8e1-426293266382
# ╠═d32f32a1-8b4f-4267-9176-82d30d623bf7
