### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 8e3f56fa-034b-11f0-1844-a38fa58e125c
begin
	using Pkg; Pkg.activate()
	using CairoMakie
	using Arya
	
	using Turing
end

# ╔═╡ 408d70ee-ab1c-4630-bd75-21358cd55489
using PlutoUI

# ╔═╡ 67ddb57a-9ee7-4f4c-848f-d0b77fba1855
using DataFrames, CSV

# ╔═╡ 6660b880-5cba-45aa-b572-e0d3588ccf17
using TuringBenchmarking

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	galaxy = ARGS[1]
	paramfile = ARGS[2]
else
	@bind galaxy confirm(TextField(default="leo2"))
	@bind paramfile confirm(TextField(default="hist_params.toml"))
end

# ╔═╡ 7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
begin 
	import PythonCall # enable fits
	using LilGuys

	FIGDIR = joinpath("..", galaxy, "figures"); FIGSUFFIX=".mcmc_spline"
end

# ╔═╡ 6014a851-525b-4565-b074-fbdd26a8ce2b


# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ 066c7b30-f818-4e4a-8db8-c8bac469f558
module MCMCUtils
	include("mcmc_utils.jl")
end

# ╔═╡ 0aa44389-f320-4274-abdd-0d7f88006a4d
log_r_label = L"$\log\,R_\textrm{ell}$ / arcmin"

# ╔═╡ 36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
log_Sigma_label = L"$\log\,\Sigma$"

# ╔═╡ 04053b71-bd55-40d7-885d-6df67035e3d6
md"""
# data loading
"""

# ╔═╡ 7e8124ea-7bbe-465b-a9dc-4b14d268c39e
obs_props = MCMCUtils.get_obs_props(galaxy)

# ╔═╡ 00380333-4c2f-4c32-9521-2764cffef265
stars = MCMCUtils.get_fits(galaxy, obs_props)

# ╔═╡ 1d505c9e-28c0-4ec4-b212-9292e25e735d
stars.L_CMD_BKD

# ╔═╡ 67b865b5-6d7d-4a62-84cd-82983a76f8ba
data = MCMCUtils.GaiaData(stars)

# ╔═╡ 46db429f-cdfb-4b6d-9484-22dd58627a8f
struct_params = MCMCUtils.StructuralParams(; LilGuys.dict_to_tuple(TOML.parsefile(joinpath("..", galaxy, "mcmc", "default_bins.toml")))...)

# ╔═╡ 0928e452-e137-4fe8-b471-d125858756a8
readdir("../leo2/mcmc/")

# ╔═╡ 43d159c4-1599-4138-96b9-c9447a8d6315
md"""
The idea with this model is to use the likelihood

``
{\cal L} \propto {\rm Prior}(\theta) \prod_i (f(\theta, x_i) P_{\rm memb,\ \it i} + (1-f(\theta, x_i)) P_{\rm bg,\it\ i})
``

where $f$ is a model for the fraction of stars in a small region belonging to the satellite.
"""

# ╔═╡ 1d9b3718-45d7-4765-ac7d-017dbf939ec8
md"""
# Histogram model
"""

# ╔═╡ 5eb1a414-86a1-44d3-856d-a0a7cc55720c
MCMCUtils.default_bins(stars)

# ╔═╡ 8e665feb-1a41-440c-a813-163cbf3be4f8
Nmemb = sum(stars.PSAT)

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Robust model
"""

# ╔═╡ e87b5ee8-e281-43e3-8172-4848018e2cd2
function read_error(key)
	if key * "_em" ∈ keys(obs_props)
		return max(obs_props[key * "_em"], obs_props[key * "_ep"])
	else
		return obs_props[key * "_err"]
	end
end

# ╔═╡ 82b3b170-729f-463e-9063-756a34840555
@time MCMCUtils.perturbed_radii(data, struct_params)

# ╔═╡ b13c0fc4-2ec0-4bed-a6b0-e87e3977c445
@time MCMCUtils.perturbed_radii(data, struct_params)

# ╔═╡ 22114c39-0e56-4665-918b-6fac012f19ae
log_Sigma_prior = Uniform(-12, 6)

# ╔═╡ 00a3b130-5ac4-464f-8e97-fceadefde4f9
import BSplineKit as BSK

# ╔═╡ 1d5de2ad-db30-4d18-845c-687bdf7ff2d1
md"""
Here, we increase the bin range on either end to allow shifts in parameter to matter less.
"""

# ╔═╡ ea7601be-e74c-4a52-8ea3-5afdfe7bc76d
sum(isnan.(MCMCUtils.Σ_hist(stars.R_ell, struct_params.bins, randn(length(struct_params.bins)-1))))

# ╔═╡ 3d25c533-6364-4c51-a1a1-a96a3f9b0dd5
Threads.nthreads()

# ╔═╡ 3ea22061-908c-4928-8fd1-d34142d20b1d
#using Profile

# ╔═╡ 1ef0a600-4c71-402a-b02a-90cc8d96b86a
#using LoggingExtras

# ╔═╡ d481d4a8-fe17-4e2b-9d56-fee13b683f3a
#import Mooncake

# ╔═╡ a5d4541b-c2a1-4d2a-9910-a48b91f667d3
# benchmark_model(mcmc_model, adbackends=[AutoForwardDiff(), AutoReverseDiff(), 
# 										#AutoMooncake(config=nothing)
# 									   ])

# ╔═╡ df3ab039-55ee-4431-bd68-1c77f843dd18
sampler = HMC(0.02, 10,
)
# 	NUTS(
# 	 100, 0.25,
# 	max_depth=10,
# 	Δ_max = 200.0,
# 	init_ϵ = 0.1,
# 	adtype = AutoReverseDiff(),
# )

# ╔═╡ ab2c3534-1325-4b1f-9405-27a43abbd0e7


# ╔═╡ dfc957ca-e66c-4b60-9f9c-4f71418ff91f
import Random

# ╔═╡ fa44af77-7031-4c16-8e5c-0d83dea63518
Nbins = length(struct_params.bins) - 1

# ╔═╡ 7fa5c2eb-5dd0-4e62-a728-7d8ab357ea20
struct_params

# ╔═╡ 36bb1876-5937-44fe-bb1c-add10470371d
pvalue = 0.16

# ╔═╡ c773432e-87fb-411e-bbed-32330b214367
bins = struct_params.bins

# ╔═╡ ea08987f-ff3e-4c9b-9590-3e057de70dea
spline_basis = BSK.BSplineBasis(BSK.BSplineOrder(3), copy(bins))

# ╔═╡ ef2785b8-91c7-4fb9-976f-7814f9c0bed6
BSK.Spline(spline_basis, randn(length(bins) +1))

# ╔═╡ 4874a13a-4095-421a-8acc-47f02132559c
function spline_log_Σ(radii, coeffs, x=bins)
	spline =  BSK.interpolate(x, coeffs, BSK.BSplineOrder(4), BSK.Natural())
	f = spline.(radii)
	return f
end

# ╔═╡ c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
@model function spline_model(data::MCMCUtils.GaiaData, params::MCMCUtils.StructuralParams=struct_params)
	Nb = length(params.bins) 

	dξ ~ Normal(0, params.position_err)
	dη ~ Normal(0, params.position_err)
	ellipticity ~ truncated(Normal(params.ellipticity, params.ellipticity_err), 
		lower=0, upper=0.99)
	position_angle ~ Normal(params.position_angle, params.position_angle_err)
	
	radii = LilGuys.calc_R_ell(data.xi .+ dξ, data.eta .+ dη, ellipticity, position_angle)
	
	coeffs ~ filldist(log_Sigma_prior, Nb)

	L_sat_space = 10 .^ spline_log_Σ(radii, coeffs)

	f_sat = @. L_sat_space  / (1 + L_sat_space)
	LL = sum(@. log10.(
		(1-f_sat) * data.P_bg 
		+ f_sat * data.P_sat
	))
	
	Turing.@addlogprob!(LL)
end

# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
begin 
	mcmc_model = spline_model(data)
	chain = sample(mcmc_model, sampler,  100, ) 
	chain
end

# ╔═╡ 115bb78a-ab2e-4ee5-bec2-b3054b42b482
begin 
	chain_summary = DataFrame(summarize(chain))
	Nvar = size(chain_summary, 1)
	log_Sigma_median = zeros(Nvar)
	log_Sigma_low = zeros(Nvar)
	log_Sigma_high = zeros(Nvar)
	for i in 1:Nvar
		log_Sigma_low[i], log_Sigma_high[i] = quantile(chain[:, i, :], [pvalue, 0.5, 1-pvalue])
	end

	chain_summary[!, :median] = log_Sigma_median
	chain_summary[!, :lower_error] = log_Sigma_median .- log_Sigma_low
	chain_summary[!, :upper_error] =  log_Sigma_high .- log_Sigma_median

	chain_summary[!, :parameters] = string.(chain_summary.parameters)
	chain_summary
end	

# ╔═╡ f3307345-6201-40af-af44-098629f84ffa
chain_summary_sigma = [startswith.(chain_summary.parameters, ["log_Sigma["]), :]

# ╔═╡ 9e8087d9-b61f-4669-9c9c-2fde55e63622
eachrow(chain_summary)[5]

# ╔═╡ ec5e0fff-0105-4c08-925a-2b771d7d78a2
let
	fontsize=4
	Nvar = size(chain_summary, 1)

	fig = Figure(size=(2*72, Nvar/4*72),
		yminorticksvisible=false
	)

	for i in 1:Nvar
		ax = Axis(fig[i, 1], 
			ylabelsize=fontsize, 
			xlabelsize=fontsize,
		  xticklabelsize=fontsize,
		  yticklabelsize=fontsize,
			ylabel=chain_summary.parameters[i],
			xlabel="step",
			yminorticksvisible=false,
		)
		
		for c in 1:size(chain.value, 3)
			y = chain[:, i, c]
			lines!((y), linewidth=0.1)
		end



		ax2 = Axis(fig[i, 2])
		hist!(vec(chain[:, i, :]), direction=:x)
		
		if i < Nvar
			hidexdecorations!(ax)
			
		end
		hidexdecorations!(ax2)
		hideydecorations!(ax2)
		linkyaxes!(ax, ax2)
	end
	rowgap!(fig.layout, 0)
	colgap!(fig.layout, 0)
	colsize!(fig.layout, 2, Relative(1/4))

	@savefig "chains"
	fig

end

# ╔═╡ 75930ff5-1add-4b1c-a495-64c50835ec4b
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		
		yticks = LogTicks(Arya.DefaultLinearTicks),
		xlabel = L"$\log\, r_\textrm{ell}$ / arcmin",
		ylabel = L"$\log \Sigma / \Sigma_\textrm{bg}$",
		#limits=(log10(bins_robust[2]), log10(bins_robust[end-1]), nothing, nothing,)
	)

	
	errorscatter!(log10.(midpoints(bins)), (chain_summary_sigma.median))
	

	@savefig "f_sat_vs_r_ell"
	fig
end

# ╔═╡ b8c73408-dffe-4dde-ae87-a427e1065251
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		
		yticks = LogTicks(Arya.DefaultLinearTicks),
		xlabel = L"$\log\, r_\textrm{ell}$ / arcmin",
		ylabel = L"$f_\textrm{sat}$",
		#limits=(log10(bins_robust[2]), log10(bins_robust[end-1]), nothing, nothing,)
	)

	xs = []
	ys = []
	yerrs = []
	
	for i in 1:length(bins)-1
		x = midpoints(bins)[i]
		y = 10 .^ median(chain[:, i, :])
		y = y ./ (1 .+ y)
		push!(xs, x)
		push!(ys, y)
		ye = 10 .^ quantile(chain[:, i, :], [pvalue, 1-pvalue])
		ye = (y-(ye[1] ./ (1 .+ ye[1])), (ye[2] ./ (1 .+ ye[2]))-y)
		push!(yerrs, ye)
	end

	@info sum(isnan.(first.(yerrs)))
	
	errorscatter!(log10.(xs), (ys), yerror=yerrs)

	@savefig "f_sat_vs_r_ell"
	fig
end

# ╔═╡ 9bf31971-1c10-4020-946e-d8cfebf6596a
md"""
# Outputs
"""

# ╔═╡ d1009491-62de-42d7-89ad-b41b8385aaaf
samplesout = joinpath("..", galaxy, "processed", "samples$FIGSUFFIX.csv")

# ╔═╡ a93579f2-d37d-492c-95e6-02bd7491dc8c
summaryout = joinpath("..", galaxy, "processed", "profile$FIGSUFFIX.toml")

# ╔═╡ a562c141-8626-4242-9f6a-c381a4da619b
function median_of(profiles, key)
	return median.(eachrow(
		hcat([getproperty(prof, Symbol(key)) for prof in profiles]...)
	))[filt_max]
end

# ╔═╡ ed16c378-2786-4f83-8152-144838dedd19
function err_of(profiles, key)
	m = median_of(profiles, key)
	h = quantile.(eachrow(
		hcat([getproperty(prof, Symbol(key)) for prof in profiles]...)
	), 1-pvalue)

	l = quantile.(eachrow(
		hcat([getproperty(prof, Symbol(key)) for prof in profiles]...)
	), pvalue)

	e = median.(eachrow(
		hcat([getproperty(prof, Symbol(key * "_err")) for prof in profiles]...)
	))

	return (@. max(h[filt_max]-m, m-l[filt_max]) + e[filt_max])
end

# ╔═╡ 36db7e1d-4b48-4510-99d5-7d567ac70d5d
df_out = DataFrame(chain)

# ╔═╡ aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
CSV.write(samplesout, df_out)

# ╔═╡ 4c0d5b06-d99b-41de-8a42-e29dbd6e0e53
CSV.write(summaryout, df_out)

# ╔═╡ Cell order:
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═6014a851-525b-4565-b074-fbdd26a8ce2b
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
# ╠═408d70ee-ab1c-4630-bd75-21358cd55489
# ╠═d1de613c-c3bb-4843-859e-7b8df54bafe0
# ╠═67ddb57a-9ee7-4f4c-848f-d0b77fba1855
# ╠═066c7b30-f818-4e4a-8db8-c8bac469f558
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╟─04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═7e8124ea-7bbe-465b-a9dc-4b14d268c39e
# ╠═00380333-4c2f-4c32-9521-2764cffef265
# ╠═1d505c9e-28c0-4ec4-b212-9292e25e735d
# ╠═67b865b5-6d7d-4a62-84cd-82983a76f8ba
# ╠═46db429f-cdfb-4b6d-9484-22dd58627a8f
# ╠═0928e452-e137-4fe8-b471-d125858756a8
# ╟─43d159c4-1599-4138-96b9-c9447a8d6315
# ╟─1d9b3718-45d7-4765-ac7d-017dbf939ec8
# ╠═5eb1a414-86a1-44d3-856d-a0a7cc55720c
# ╠═8e665feb-1a41-440c-a813-163cbf3be4f8
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═e87b5ee8-e281-43e3-8172-4848018e2cd2
# ╠═82b3b170-729f-463e-9063-756a34840555
# ╠═b13c0fc4-2ec0-4bed-a6b0-e87e3977c445
# ╠═22114c39-0e56-4665-918b-6fac012f19ae
# ╠═00a3b130-5ac4-464f-8e97-fceadefde4f9
# ╠═ea08987f-ff3e-4c9b-9590-3e057de70dea
# ╠═ef2785b8-91c7-4fb9-976f-7814f9c0bed6
# ╠═4874a13a-4095-421a-8acc-47f02132559c
# ╠═c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
# ╠═1d5de2ad-db30-4d18-845c-687bdf7ff2d1
# ╠═ea7601be-e74c-4a52-8ea3-5afdfe7bc76d
# ╠═3d25c533-6364-4c51-a1a1-a96a3f9b0dd5
# ╠═3ea22061-908c-4928-8fd1-d34142d20b1d
# ╠═1ef0a600-4c71-402a-b02a-90cc8d96b86a
# ╠═6660b880-5cba-45aa-b572-e0d3588ccf17
# ╠═d481d4a8-fe17-4e2b-9d56-fee13b683f3a
# ╠═a5d4541b-c2a1-4d2a-9910-a48b91f667d3
# ╠═df3ab039-55ee-4431-bd68-1c77f843dd18
# ╠═ab2c3534-1325-4b1f-9405-27a43abbd0e7
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═dfc957ca-e66c-4b60-9f9c-4f71418ff91f
# ╠═fa44af77-7031-4c16-8e5c-0d83dea63518
# ╠═7fa5c2eb-5dd0-4e62-a728-7d8ab357ea20
# ╠═115bb78a-ab2e-4ee5-bec2-b3054b42b482
# ╠═f3307345-6201-40af-af44-098629f84ffa
# ╠═ec5e0fff-0105-4c08-925a-2b771d7d78a2
# ╠═36bb1876-5937-44fe-bb1c-add10470371d
# ╠═75930ff5-1add-4b1c-a495-64c50835ec4b
# ╠═9e8087d9-b61f-4669-9c9c-2fde55e63622
# ╠═b8c73408-dffe-4dde-ae87-a427e1065251
# ╠═c773432e-87fb-411e-bbed-32330b214367
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═d1009491-62de-42d7-89ad-b41b8385aaaf
# ╠═a93579f2-d37d-492c-95e6-02bd7491dc8c
# ╠═a562c141-8626-4242-9f6a-c381a4da619b
# ╠═ed16c378-2786-4f83-8152-144838dedd19
# ╠═36db7e1d-4b48-4510-99d5-7d567ac70d5d
# ╠═aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
# ╠═4c0d5b06-d99b-41de-8a42-e29dbd6e0e53
