### A Pluto.jl notebook ###
# v0.20.18

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

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	galaxy = ARGS[1]
else
	@bind galaxy confirm(TextField(default="leo2"))
end

# ╔═╡ 05517bcc-7967-4bc7-9396-c271e420665d
import PairPlots

# ╔═╡ 6014a851-525b-4565-b074-fbdd26a8ce2b
outdir = joinpath("..", galaxy, "mcmc")

# ╔═╡ 7ff855b9-7a1e-422e-b8a5-1cb5ecfe368f
begin 
	# import PythonCall # enable fits
	using LilGuys

	FIGDIR = joinpath(outdir, "figures"); FIGSUFFIX=".mcmc_exp"
end

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

# ╔═╡ 84708283-6f48-4ab8-88f8-10b2f9376466
R_max = 100 * obs_props["R_h"]

# ╔═╡ 1d505c9e-28c0-4ec4-b212-9292e25e735d
stars.L_CMD_BKD

# ╔═╡ 67b865b5-6d7d-4a62-84cd-82983a76f8ba
data = MCMCUtils.GaiaData(stars)

# ╔═╡ 46db429f-cdfb-4b6d-9484-22dd58627a8f
struct_params = MCMCUtils.StructuralParams(; LilGuys.dict_to_tuple(TOML.parsefile(joinpath("..", galaxy, "mcmc", "default_bins.toml")))...)

# ╔═╡ 43d159c4-1599-4138-96b9-c9447a8d6315
md"""
The idea with this model is to use the likelihood

``
{\cal L} \propto {\rm Prior}(\theta) \prod_i (f(\theta, x_i) P_{\rm memb,\ \it i} + (1-f(\theta, x_i)) P_{\rm bg,\it\ i})
``

where $f$ is a model for the fraction of stars in a small region belonging to the satellite.
"""

# ╔═╡ 8e665feb-1a41-440c-a813-163cbf3be4f8
Nmemb = sum(stars.PSAT)

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Robust model
"""

# ╔═╡ 22114c39-0e56-4665-918b-6fac012f19ae
log_Sigma_prior = Uniform(-12, 6)

# ╔═╡ c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
@model function sersic_model(xi::Vector{Float64}, eta::Vector{Float64}, P_sat::Vector{Float64}, P_bg::Vector{Float64},)
	d_xi ~ Normal(0, 5)
	d_eta ~ Normal(0, 5)
	ellipticity ~ Uniform(0, 0.99)
	position_angle ~ Uniform(0, 180)
	
	R_s ~ LogNormal(0, 5)
	prof = LilGuys.Exp2D(R_s=R_s, M=1)

	R = LilGuys.calc_R_ell(xi .+ d_xi, eta .+ d_eta, ellipticity, position_angle)
	
	L_sat_space = @. LilGuys.surface_density.(prof, R)
	L_bg_space = 1 / (π*R_max^2)

	f_sat ~ Uniform(0, 1)
	LL = sum(@. log10.(
		(1-f_sat) * P_bg * L_bg_space
		+ f_sat * P_sat * L_sat_space
	))
	
	Turing.@addlogprob!(LL)
end

# ╔═╡ f96b8939-bf1b-4ae1-bf80-89c984198b2d
@model function sersic_model_simple(R::Vector{Float64}, P_sat::Vector{Float64}, P_bg::Vector{Float64},)
	R_s ~ LogNormal(2, 2)
	prof = LilGuys.Exp2D(R_s=R_s, M=1)

	L_sat_space = @. LilGuys.surface_density.(prof, R)
	Σ_bg = 1 / (π*R_max^2)
	f_sat ~ Uniform(0, 1)

	LL = sum(@. log10.(
		(1-f_sat) * Σ_bg * P_bg 
		+ f_sat * L_sat_space * P_sat
	))
	
	Turing.@addlogprob!(LL)
end

# ╔═╡ 3d25c533-6364-4c51-a1a1-a96a3f9b0dd5
Threads.nthreads()

# ╔═╡ df3ab039-55ee-4431-bd68-1c77f843dd18

sampler = 	NUTS(
 # 	 1000, 0.65,
	# max_depth = 5,
	# Δ_max = 200.0,
	# init_ϵ = 0.1,
	# adtype = AutoReverseDiff(),
)
#sampler = MH()

# ╔═╡ c9ba2348-7d79-4b86-bb0a-f6f7c4bef22e
R_ell = LilGuys.calc_R_ell(data.xi, data.eta, struct_params.ellipticity, struct_params.position_angle)

# ╔═╡ c6c8fc47-f607-41e9-ad97-2d5ecd4bfde5
sum(R_ell .> R_max)

# ╔═╡ 3db5c5fc-c1cd-4975-9b44-8c0966268c1f
filt = R_ell .< R_max

# ╔═╡ 7b1b4c0f-aa49-4ee0-b860-0bd927db8768
mcmc_model = sersic_model(data.xi[filt], data.eta[filt], data.P_sat[filt], 1data.P_bg[filt])

# ╔═╡ d357c1f6-559d-4c19-bb2a-9b86bc9a285e
mcmc_model_simple = sersic_model_simple(R_ell[filt], data.P_sat[filt], data.P_bg[filt])

# ╔═╡ 8550f4a6-a992-4504-b2f7-5049631faeb3
chain_simple = sample(mcmc_model_simple, sampler, MCMCThreads(), 1000, 16) 

# ╔═╡ 8f051c8a-def2-4a84-ab43-2ecc8b646b65
begin 
	chain = sample(mcmc_model, sampler, MCMCThreads(), 1000, 16) 
	chain
end

# ╔═╡ 00122357-4751-478c-af5e-743fed1d8505
df_prior = sample(mcmc_model, Prior(), 1000) |> DataFrame

# ╔═╡ d1600a1c-1a73-4a0e-a36f-3008a5e9ed23
prof_obs = LilGuys.SurfaceDensityProfile(R_ell[stars.PSAT .> 0.5]) |> LilGuys.filter_empty_bins

# ╔═╡ 51e22997-a5ec-48b3-9719-9ae0e26cc20c
prof_all = LilGuys.SurfaceDensityProfile(R_ell)

# ╔═╡ 6f6c7331-0875-419b-abdb-834057523ee8
Σ_bg = sum(R_ell[stars.PSAT .< 0.01] .< obs_props["R_h"]) / (π*obs_props["R_h"]^2)

# ╔═╡ ea2665ee-e24c-45df-8986-96ba52335882
log10(Σ_bg)

# ╔═╡ 428e8ecd-b28d-4c4e-b662-fb5d75e1ba4e
PairPlots.pairplot(chain)

# ╔═╡ 0b5bd493-8d3a-45d2-a00b-12b721ef6e78
PairPlots.pairplot(chain_simple)

# ╔═╡ 36bb1876-5937-44fe-bb1c-add10470371d
pvalue = 0.16

# ╔═╡ 115bb78a-ab2e-4ee5-bec2-b3054b42b482
begin 
	chain_summary = DataFrame(summarize(chain))
	Nvar = size(chain_summary, 1)
	log_Sigma_median = zeros(Nvar)
	log_Sigma_low = zeros(Nvar)
	log_Sigma_high = zeros(Nvar)
	for i in 1:Nvar
		log_Sigma_low[i],log_Sigma_median[i], log_Sigma_high[i] = quantile(chain[:, i, :], [pvalue, 0.5, 1-pvalue])
	end

	chain_summary[!, :median] = log_Sigma_median
	chain_summary[!, :lower_error] = log_Sigma_median .- log_Sigma_low
	chain_summary[!, :upper_error] =  log_Sigma_high .- log_Sigma_median

	chain_summary[!, :parameters] = string.(chain_summary.parameters)
	chain_summary
end	

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
			ylabelrotation=0,
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

# ╔═╡ c773432e-87fb-411e-bbed-32330b214367
bins = struct_params.bins

# ╔═╡ 9bf31971-1c10-4020-946e-d8cfebf6596a
md"""
# Outputs
"""

# ╔═╡ d1009491-62de-42d7-89ad-b41b8385aaaf
samplesout = joinpath(outdir, "samples$FIGSUFFIX.csv")

# ╔═╡ a93579f2-d37d-492c-95e6-02bd7491dc8c
summaryout = joinpath(outdir, "summary$FIGSUFFIX.csv")

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

# ╔═╡ d57d7605-3904-490b-b785-42320275b0c5
mean(df_out.acceptance_rate), mean(df_out.is_accept)

# ╔═╡ bcbcd4fb-e5d0-4e46-bfba-463e5e10d758
prof = LilGuys.Exp2D(R_s=median(df_out.R_s), M=10^median(df_out.log_M))

# ╔═╡ 1dae5f85-2720-4bcb-b413-8c5a6857b1c3
L_s_new = LilGuys.surface_density.(prof, R_ell)

# ╔═╡ a8dd3c6e-be4e-426f-ad46-0953e0c2149a
f_s_new = L_s_new ./ (1 .+ L_s_new)

# ╔═╡ b08ddf65-0b5e-40e1-ab57-459cce6745ba
PSAT_new = @. f_s_new * data.P_sat ./ (f_s_new * data.P_sat + (1-f_s_new) * data.P_bg)

# ╔═╡ 26bb9598-7150-4014-b789-a2b8ac82c48b
prof_new = LilGuys.SurfaceDensityProfile(R_ell[PSAT_new .> 0.5])

# ╔═╡ 2b780e99-54e3-4b90-beb9-85501449ff74
let
	fig = Figure()
	ax = Axis(fig[1,1])


	errorscatter!(prof_obs.log_R, (prof_obs.log_Sigma), yerror=error_interval.(prof_obs.log_Sigma))

	errorscatter!(prof_all.log_R, (prof_all.log_Sigma), yerror=error_interval.(prof_all.log_Sigma))
	# errorscatter!(prof_new.log_R, (prof_new.log_Sigma), yerror=error_interval.(prof_new.log_Sigma))

	
	x = LinRange(extrema(prof_obs.log_R)..., 1000)

	R = 10 .^ x

	df = df_out
	for i in 1:100
		M = size(stars, 1) * df.f_sat[i]
		prof = LilGuys.Exp2D(R_s=df.R_s[i], M=M)
		y = @. log10(LilGuys.surface_density(prof, R))

		lines!(x, y, color=COLORS[3], alpha=0.1)
	end

	ylims!(-6, 3)

	fig

end

# ╔═╡ aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
CSV.write(samplesout, df_out)

# ╔═╡ 4c0d5b06-d99b-41de-8a42-e29dbd6e0e53
CSV.write(summaryout, chain_summary)

# ╔═╡ Cell order:
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═05517bcc-7967-4bc7-9396-c271e420665d
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═6014a851-525b-4565-b074-fbdd26a8ce2b
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
# ╠═84708283-6f48-4ab8-88f8-10b2f9376466
# ╠═c6c8fc47-f607-41e9-ad97-2d5ecd4bfde5
# ╠═1d505c9e-28c0-4ec4-b212-9292e25e735d
# ╠═67b865b5-6d7d-4a62-84cd-82983a76f8ba
# ╠═46db429f-cdfb-4b6d-9484-22dd58627a8f
# ╟─43d159c4-1599-4138-96b9-c9447a8d6315
# ╠═8e665feb-1a41-440c-a813-163cbf3be4f8
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═22114c39-0e56-4665-918b-6fac012f19ae
# ╠═c38c2f54-1a6d-4bfd-966a-0ddf66ab94da
# ╠═f96b8939-bf1b-4ae1-bf80-89c984198b2d
# ╠═3d25c533-6364-4c51-a1a1-a96a3f9b0dd5
# ╠═df3ab039-55ee-4431-bd68-1c77f843dd18
# ╠═3db5c5fc-c1cd-4975-9b44-8c0966268c1f
# ╠═7b1b4c0f-aa49-4ee0-b860-0bd927db8768
# ╠═c9ba2348-7d79-4b86-bb0a-f6f7c4bef22e
# ╠═d357c1f6-559d-4c19-bb2a-9b86bc9a285e
# ╠═8550f4a6-a992-4504-b2f7-5049631faeb3
# ╠═8f051c8a-def2-4a84-ab43-2ecc8b646b65
# ╠═00122357-4751-478c-af5e-743fed1d8505
# ╠═d57d7605-3904-490b-b785-42320275b0c5
# ╠═115bb78a-ab2e-4ee5-bec2-b3054b42b482
# ╠═d1600a1c-1a73-4a0e-a36f-3008a5e9ed23
# ╠═51e22997-a5ec-48b3-9719-9ae0e26cc20c
# ╠═6f6c7331-0875-419b-abdb-834057523ee8
# ╠═ea2665ee-e24c-45df-8986-96ba52335882
# ╠═bcbcd4fb-e5d0-4e46-bfba-463e5e10d758
# ╠═1dae5f85-2720-4bcb-b413-8c5a6857b1c3
# ╠═a8dd3c6e-be4e-426f-ad46-0953e0c2149a
# ╠═b08ddf65-0b5e-40e1-ab57-459cce6745ba
# ╠═26bb9598-7150-4014-b789-a2b8ac82c48b
# ╠═2b780e99-54e3-4b90-beb9-85501449ff74
# ╠═428e8ecd-b28d-4c4e-b662-fb5d75e1ba4e
# ╠═0b5bd493-8d3a-45d2-a00b-12b721ef6e78
# ╠═ec5e0fff-0105-4c08-925a-2b771d7d78a2
# ╠═36bb1876-5937-44fe-bb1c-add10470371d
# ╠═c773432e-87fb-411e-bbed-32330b214367
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═d1009491-62de-42d7-89ad-b41b8385aaaf
# ╠═a93579f2-d37d-492c-95e6-02bd7491dc8c
# ╠═a562c141-8626-4242-9f6a-c381a4da619b
# ╠═ed16c378-2786-4f83-8152-144838dedd19
# ╠═36db7e1d-4b48-4510-99d5-7d567ac70d5d
# ╠═aa95112f-01d5-45cb-9c5c-b1c7e0ee7e45
# ╠═4c0d5b06-d99b-41de-8a42-e29dbd6e0e53
