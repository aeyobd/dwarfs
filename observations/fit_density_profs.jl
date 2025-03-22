### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ fa54e7aa-f477-45db-8dd9-56bd6b367604
begin
	import Pkg; Pkg.activate()

	using CairoMakie
	using Arya
end

# ╔═╡ cd1e43e6-65e1-4efe-9290-62633ee340a9
using PlutoUI

# ╔═╡ 51d52048-cd18-41ff-ac2b-41dfc097de16
using DataFrames

# ╔═╡ aae83e18-e739-4421-a763-37590cb7869a
using PairPlots

# ╔═╡ 60763b01-9c7b-4153-9b94-617820517927
using Optimization, OptimizationOptimJL

# ╔═╡ 37ee41d1-7b98-44f7-82ec-5b496622931b
using Turing

# ╔═╡ d0112dc0-8cc1-4fdf-b30d-a3b89660d9b4
begin
	using Measurements
	import LinearAlgebra: diag
	import NaNMath as nm
end

# ╔═╡ 27bd8a46-ee56-11ef-22e4-afe88d73af2d
md"""
The goal of this notebook is to fit the density profiles to e.g. King, Sersic, Exponential and save the best fit parameters for any dwarf + density profile combination.
"""

# ╔═╡ 9d3245ef-f774-4add-9b1d-60bb69b04c8d
md"""
The only inputs are the name of the galaxy and the name of the profile to load (in galaxyname/density_profiles/) 
"""

# ╔═╡ bd933a63-45e8-4e6c-9a0a-5d20c7d55eaf
md"""
# Setup
"""

# ╔═╡ 6cf3cb83-317e-40fa-be47-7a4e31a539ec
CairoMakie.activate!(pt_per_unit=2, type="svg")

# ╔═╡ 99315b4c-37c1-4287-bd71-c3339c8e76ad
import FillArrays: I, Eye

# ╔═╡ c3d9c4f3-6891-46e1-822e-5c217450f59b
import LinearAlgebra: diagm, identity

# ╔═╡ 56af9242-91cc-46db-b532-c127450adfeb
import Random: randperm

# ╔═╡ 62cce84c-e985-4bd2-8c2d-ab15ed239b99
function notebook_inputs(; kwargs...)
	return PlutoUI.combine() do Child
		
		user_inputs = [
			md""" $(string(name)): $(
				Child(name, obj)
			)"""
			
			for (name, obj) in kwargs
		]
		
		md"""
		#### Inputs
		$(user_inputs)
		"""
	end
end

# ╔═╡ 72524f04-bd5b-4f4b-a41d-8c038c180ce0
@bind inputs confirm(notebook_inputs(;
	galaxyname = TextField(default="ursa_minor"),
	profilename = TextField(default="../processed/profile.mcmc_hist.toml")
))

# ╔═╡ 1238d232-a2c0-44e5-936b-62fd9137d552
profilename = inputs.profilename

# ╔═╡ 25c3b74f-d8fa-4f4c-ae7c-71b137fb2ce7
galaxyname = inputs.galaxyname

# ╔═╡ 20b8c05c-91a6-4dc2-9558-1b02be4983f2
begin
	using LilGuys
	FIGDIR = joinpath(galaxyname, "figures")
	FIGSUFFIX=".fit_density_profs"
end

# ╔═╡ 11b84ff0-7f35-4716-ae30-e05a1c5f88ba
log_R_label = L"\log\, r\ /\ \mathrm{arcmin}"

# ╔═╡ 3eb9992d-425d-4bec-94d2-da4f8dc59a8d
log_sigma_label = L"$\log\, \Sigma $ /\ stars\ arcmin$^{-2}$"

# ╔═╡ db8a3705-f4ea-4588-9885-0e03cd42c47b
if !isdir(FIGDIR)
	mkdir(FIGDIR)
end

# ╔═╡ 8fdc3d18-6a3a-42eb-b0cb-d5cae1cd65d9
md"""
## Loading Data
"""

# ╔═╡ 20e017f2-c90b-46a5-a6f1-8e4aeb5ffde8
begin
	prof_in = LilGuys.StellarProfile(joinpath(galaxyname, "density_profiles/$(profilename)"))
	if prof_in.normalization == "mass"
		prof_in.log_Sigma .-= prof_in.log_m_scale
		prof_in.normalization = "none"
	end

	prof_in
end

# ╔═╡ 249fed2c-e86a-47ab-a7a6-06f7c1edc01a


# ╔═╡ 44a893d6-9a3d-4959-b76f-274635198aa2
md"""
not sure the bins are ddone right but should be okay???
"""

# ╔═╡ a0cb6ecd-59ab-4a7a-b95a-d55bfdd2d93b
begin
	_filt_prof = prof_in.log_Sigma_em .< 1
	_filt_prof .&= prof_in.log_Sigma_ep .< 1
	_filt_prof .&= isfinite.(prof_in.log_Sigma_em)
	_filt_prof .&= isfinite.(prof_in.log_Sigma)
	_filt_prof .&= isfinite.(prof_in.log_R)
	_filt_bins = [false; _filt_prof] .|| [_filt_prof; false]

	prof = (;
		log_R = prof_in.log_R[_filt_prof],
		log_Sigma = prof_in.log_Sigma[_filt_prof],
		log_Sigma_err = max.(prof_in.log_Sigma_em[_filt_prof] .+ prof_in.log_Sigma_ep[_filt_prof]),
		log_R_bins = prof_in.log_R_bins[_filt_bins],
	)
end

# ╔═╡ 995a9a15-38d0-47c3-ba12-58ba8209497c


# ╔═╡ 5026118a-85f6-4f71-bd0f-d25c4916c834
"a list of log radii to sample analytic models over"
log_r_model = LinRange(prof.log_R_bins[1], prof.log_R_bins[end], 1000)

# ╔═╡ 30154472-316b-4489-b1e5-f503447928cb
md"""
## Utility functions for ML optimization
"""

# ╔═╡ 4970ae98-a08d-4308-b67e-e9306d4193a3
md"""
These are not as clean, more as a double check
"""

# ╔═╡ 8d90dd5c-350a-45e6-b04c-9b2f85cf700b
function objective_function(x, u, analytic_profile, names; kwargs...)
	dlogx = u[1]
	dlogy = u[2]
	ana_kwargs = [Symbol.(sym) => ui for (sym, ui) in zip(names, u[3:end])]
	
	h = analytic_profile(; ana_kwargs..., kwargs...)
	
	y_m = nm.log10.(LilGuys.calc_Σ.(h, 10 .^ (x .- dlogx))) .+ dlogy
	return y_m
end

# ╔═╡ aa22fa06-8167-4d6e-aa05-399800370c1f
function filter_prof(prof, log_R_max)
	i_cut = findfirst(prof.log_R .> log_R_max)
	if i_cut === nothing
		return prof
	end
	i_cut = length(prof.log_R) - i_cut  + 1
	@info "cutting $i_cut"
	return (;
		log_R = prof.log_R[begin:end-i_cut],
		log_R_bins = prof.log_R_bins[begin:end-i_cut],
		log_Sigma = prof.log_Sigma[begin:end-i_cut],
		log_Sigma_err = prof.log_Sigma_err[begin:end-i_cut],
	)
end

# ╔═╡ 7ff0bee2-3622-4573-81d5-8b8fa96cf461
function fit_profile(prof, analytic_profile, names, p0; kwargs...)
	
	objective(x, u) = objective_function(x, u, analytic_profile, names; kwargs...)

	w = 1 ./prof.log_Sigma_err .^ 2

	local popt, covt, errs
	try
		popt, covt = LilGuys.curve_fit(objective, prof.log_R, prof.log_Sigma, w,  p0, )

	catch e
		@warn e
		popt = fill(NaN, length(names) + 2)
		errs = fill(NaN, length(names) + 2)
	else
		errs = sqrt.(diag(covt))
	end

	y_m = objective(prof.log_R, popt)
	res = prof.log_Sigma .- y_m
	chi2 = sum((res) .^ 2 .* w) ./ (length(w) - length(popt))


	result = OrderedDict()

	names_all = ["log_R_s"; "log_Sigma"; names]

	for (i, name) in enumerate(names_all)
		result[name] = popt[i]
		result[name * "_err"] = errs[i]
	end

	result["chi2"] = chi2
	result["popt"] = popt
	result["residuals"] = res
	result["log_Sigma_pred"] = y_m

	x_m = LinRange(-0.1 + minimum(prof.log_R), 0.1 + maximum(prof.log_R), 1000)
	y_m = objective(x_m, popt)
	result["log_Sigma_pred"] = y_m
	result["log_R_pred"] = x_m

	result["R_s"] = 10 .^ result["log_R_s"]
	result["R_s_err"] =  result["log_R_s_err"] * log(10) * result["R_s"]

	if analytic_profile <: LilGuys.Sersic
		result["Mtot"] = LilGuys.calc_M_2D(LilGuys.Sersic(), 1000) * 10 .^ (result["log_Sigma"] + 2result["log_R_s"])
	else
		
		kwargs = [Symbol.(sym) => ui for (sym, ui) in zip(names, popt[3:end])]
		h = analytic_profile(; M=1, kwargs...)	
	
		result["Mtot"] = LilGuys.get_M_tot(h) * 10 .^ (result["log_Sigma"]  + 2result["log_R_s"])
	end

	result["log_Mtot"] = log10.(result["Mtot"])
	result["log_Mtot_err"] = result["log_Sigma_err"] .+ 2result["log_R_s_err"]
	

	return result
end

# ╔═╡ 7a274a8d-f6eb-4d9b-97c5-0372dc9520fe
function plot_Σ_fit_res(obs, fit; res_max=1, nf=2, title="")
	fig = Figure(figsize=(3.25*72, 3.25*72))
    ax = Axis(fig[1, 1], 
        ylabel=log_sigma_label, 
		title=title,
	)
	
    errorscatter!(ax, obs.log_R, obs.log_Sigma, yerror=obs.log_Sigma_err)

	log_R = fit["log_R_pred"]
	pred = fit["log_Sigma_pred"]

    lines!(ax, log_R, pred, color=COLORS[2])
	
    lines!(ax, log_R, pred, color=COLORS[2], linestyle=:dash)
    
    ax2 = Axis(fig[2, 1],
        ylabel=L"\delta\log\Sigma", 
    	xlabel=log_R_label,
		limits = (nothing, (-res_max, res_max))
	)


	res = fit["residuals"]
    errorscatter!(ax2, obs.log_R, res, yerror=obs.log_Sigma_err, label="")

    hlines!(0, color=:black)
    
    rowsize!(fig.layout, 2, Relative(1/4))

    linkxaxes!(ax, ax2)
    hidexdecorations!(ax, grid=false)
    return fig
end

# ╔═╡ f91ae13d-4533-45aa-b835-d53fab594360
md"""
## Utility functions for MCMC optimization
"""

# ╔═╡ 3abb1fb9-1e3a-4251-abfd-959d1faa0f87
"""
	DensityModel(log_Σ, turing_model, argnames)

specifies a density model to run an MCMC fit on.
The attributes are 
- `log_Σ`. A function of `log_r` returning `log_Σ` using the kwargs specified in `argnames`
- `turing_model` a turing model specification (from `@model function...`) which takes four arguments: `log_Σ`, `log_r_obs`, `log_Σ_obs`, and `log_Σ_err_obs`
- `argnames`: kwargs in density profile model.
"""
struct DensityModel
	log_Σ::Function
	turing_model::Function
	argnames::Array{Symbol}
end

# ╔═╡ 95bd86c2-0222-4ef9-a8ef-d26e023fc171
"""
	MCMCFit(model, samples, summary, χ2, lp)

The MCMC result based on the specified model.
"""
struct MCMCFit
	model::DensityModel
	samples::Chains
	summary::DataFrame
	χ2::Float64
	lp::Float64
end

# ╔═╡ 3981b6f9-7ace-4de9-b72d-10168512689c
"""
	sample_model(model, log_Σ, prof; threads, chains, kwargs...)

Samples a model using the given number of threads and chains on each thread. Args passed to Turing.
"""
function sample_model(turing_model, log_Σ, prof; sampler=NUTS(), threads=4, chains=16, kwargs...)
	model = turing_model(log_Σ, prof.log_R, prof.log_Sigma, prof.log_Sigma_err)

	chains_per_thread = round(Int, chains / threads)
	samples = mapreduce(c -> sample(model, sampler, MCMCThreads(), 1000, threads; kwargs...), chainscat, 1:chains_per_thread)

	return samples
end

# ╔═╡ 3d7e4658-97d5-4aff-baff-2ee87c688871
"""
	predict(mcmc_fit, log_R)

Predicts log_Σ given log_R and a mcmc fit model (assuming median best parameters).
"""
function predict(mcmc_fit::MCMCFit, log_R::AbstractVector{<:Real})
	kwargs = Dict(arg => mcmc_fit.summary[only(findfirst(mcmc_fit.summary.parameters .== arg)), :median] for arg in mcmc_fit.model.argnames)

	return mcmc_fit.model.log_Σ(log_R; kwargs...)
end

# ╔═╡ 5e31ad85-fbc5-459c-b6d4-b3493d3fd29f
@doc raw"""
	predict(mcmc_fit, log_R)

Predicts $\log \Sigma$ given  $\log r$, an mcmc_model, and the specified iteration & chain.
"""
function predict(mcmc_fit::MCMCFit, log_R::AbstractVector{<:Real}, i::Integer, c::Integer)
	kwargs = Dict(arg => mcmc_fit.samples.value[i, arg, c] for arg in mcmc_fit.model.argnames)

	return mcmc_fit.model.log_Σ(log_R; kwargs...)
end

# ╔═╡ 9a02adf8-9b90-4b4b-afcb-3ec5170f538b
function get_χ2(fit, prof)
	res = prof.log_Sigma .- predict(fit, prof.log_R) 
	χ2 = sum(res .^ 2 ./ prof.log_Sigma_err .^ 2)

	χ2_red = χ2 / (length(prof.log_R) - length(fit.model.argnames))

	return χ2_red

end

# ╔═╡ 6c160b4e-177e-4b9e-9db4-3e2aa5465ce7
"""
	get_lp(samples::Chains)

Gets the best value of the log-posterior from MCMC samples
"""
function get_lp(samples::Chains)
	df = DataFrame(samples)
	lp = maximum(samples.value[:, :lp, :])
	return lp
end

# ╔═╡ 776d6ffd-9598-45cc-b8d9-6c6b4025583a
function summarize_chain(chain; p=0.16)
	df = DataFrame(chain)
	summary = DataFrame(Turing.summarize(chain))

	Nr = size(summary, 1)
	meds = zeros(Nr)
	err_low = zeros(Nr)
	err_high = zeros(Nr)
	
	for i in 1:Nr
		sym = summary[i, :parameters]
		meds[i] = median(df[!, sym])
		err_low[i] = meds[i] - quantile(df[!, sym], p)
		err_high[i] = quantile(df[!, sym], 1-p) - meds[i]
	end

	summary[!, :median] = meds
	summary[!, :err_low] = err_low
	summary[!, :err_high] = err_high

	select!(summary, :parameters, :median, :err_low, :err_high, Not([:parameters, :median, :err_low, :err_high]))
	return summary
end

# ╔═╡ 47ab55c6-e202-47a0-b1d5-684c9a1fef5d
"""
	sample_model(density_model; prof, kwargs...)

Samples the provided density model assuming observations in `prof`. 
Kwargs passed to turing.
"""
function sample_model(density_model::DensityModel; 
		prof=prof, 
		kwargs...
	)
	
	samples = sample_model(density_model.turing_model, density_model.log_Σ, prof; kwargs...)
	
	summary = summarize_chain(samples)

	# easier to reconstruct model here for chi2
	fit = MCMCFit(density_model, samples, summary, NaN, get_lp(samples))
	χ2 = get_χ2(fit, prof)
	return 	MCMCFit(density_model, samples, summary, χ2, get_lp(samples))
end

# ╔═╡ 6c59342a-36b8-458b-8b17-92b0f18afb7d
md"""
###  MCMC Plots
"""

# ╔═╡ 69b919f7-5985-4835-950a-ab7ad998d3fd
function subsample(samples, N=300)
	df = DataFrame(samples)
	Ntot = size(df, 1)
	if N > Ntot
		@error "more draws requested than samples"
	end

	idx = randperm(Ntot)[1:N]
	return df[idx, :]
end

# ╔═╡ 071356e5-53c4-462e-8448-7e5480dbfc0e
function plot_samples!(df, log_Σ, argnames; 
		x = log_R_model, 
		alpha = 0.01, 
		dy = zeros(length(x)),
		kwargs...
	)
	
	for row in eachrow(df)
		kwargs_prof = Dict(argname => row[argname] for argname in argnames)
		y = log_Σ(x; kwargs_prof...)

		lines!(x, y .+ dy, color=(:black, alpha); kwargs...)
	end
end

# ╔═╡ 566f2f87-a840-416e-b71e-bcfbb747e39b
function plot_chains(samples)
	fig = Figure(size=(3*72, 6*72))

	Nc = size(samples, 3)

	labels = samples.info[:varname_to_symbol]

	for (i, label) in enumerate(labels)
		ax = Axis(fig[i, 1],
			xlabel = "step",
			ylabel = last(label) |> string
		)
		
		for chain in 1:Nc
			x = samples.value[:, i, chain]
			lines!(x)
		end

		if i < length(labels)
			hidexdecorations!(ax, ticks=false, grid=false, minorticks=false)
		end
	end

	rowgap!(fig.layout, 0)
	linkxaxes!(fig.content...)
	fig
end

# ╔═╡ 2bc41ec4-6ced-474e-92e8-334551b38431


# ╔═╡ e31c35b1-9b37-4c95-8344-6177b685d5ab
function plot_samples_obs(mcmc_fit)
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_R_label, 
		ylabel = log_sigma_label,
	)

	
	plot_samples!(mcmc_fit.samples, mcmc_fit.log_Σ, mcmc_fit.model.argnames)
	errorscatter!(prof.log_R, prof.log_Sigma, yerror=prof.log_Sigma_err, color=COLORS[2])

	LilGuys.hide_grid!(ax)

	median =
	fig
end

# ╔═╡ 7779ec2e-1053-47e9-8d51-897d5063bb85
function plot_samples_residuals(mcmc_fit; res_ylims=nothing)
	samples = mcmc_fit.samples
	log_Σ = mcmc_fit.model.log_Σ
	argnames = mcmc_fit.model.argnames
	summary = mcmc_fit.summary
	
	x = LinRange(prof.log_R_bins[1], prof.log_R_bins[end], 1000)
	df = subsample(samples)

	fig = Figure(figsize=(3.25*72, 3.25*72))
	ax = Axis(fig[1,1],
		xlabel = log_R_label, 
		ylabel = log_sigma_label,
	)

	plot_samples!(df, log_Σ, argnames, x=x)
	errorscatter!(prof.log_R, prof.log_Sigma, yerror=prof.log_Sigma_err, color=COLORS[2])

	LilGuys.hide_grid!(ax)
	hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)

	# residuals
	ax_res = Axis(fig[2,1],
		xlabel = log_R_label, 
		ylabel = "residual",
	)


	kwargs_med = Dict(arg => summary[only(findfirst(summary.parameters .== arg)), :median] for arg in argnames)
	
	y_m = log_Σ.(x; kwargs_med...)

	plot_samples!(df, log_Σ, argnames, dy=-y_m)
	LilGuys.hide_grid!(ax_res)

	hlines!(0, color=:black)
	

	y_m = log_Σ.(prof.log_R; kwargs_med...)

	errorscatter!(prof.log_R, prof.log_Sigma .- y_m, yerror=prof.log_Sigma_err, color=COLORS[2])

	if isnothing(res_ylims)
		y_max = maximum(prof.log_Sigma .- y_m .+ prof.log_Sigma_err)
		y_min = minimum(prof.log_Sigma .- y_m .- prof.log_Sigma_err)

		y_abs_max = max(abs(y_max), abs(y_min))

		res_ylims = (-y_abs_max, y_abs_max)
	end

	ylims!(ax_res, res_ylims)

	rowsize!(fig.layout, 2, Relative(0.3))
	fig
end

# ╔═╡ 7115482e-0cc9-40ef-b63e-0ebad62bb949
md"""
# Fits
"""

# ╔═╡ a9e5aa14-6a11-436c-97ca-7afd5fd35478
md"""
## Sérsic
"""

# ╔═╡ a4d8323a-c668-452b-a958-17982177cac6
function log_Σ_sersic(log_R; log_M_s, log_R_h, n, log_Sigma_bg)
	r = 10 .^ log_R
	prof = LilGuys.Sersic(n=n, R_h=10 .^ log_R_h, _b_n=LilGuys.guess_b_n(n))
	Σ = LilGuys.calc_Σ.(prof, r)
	y = @. log10(Σ + 10^log_Sigma_bg)  .+ log_M_s

	return y
end

# ╔═╡ 2a3adaf5-2a0d-4674-a160-a7844b7f1463
@model function sersic_turing_model(model_log_Σ, log_R, log_Σ, log_Σ_err)
	n ~ Gamma(1)
	log_R_h ~ Normal(0.0, 1.0)
	log_M_s ~ Normal(4.0, 2.0)
	log_Sigma_bg ~ Normal(-4.0, 4.0)

	log_σ ~ Normal(-1, 1)
	
	σ = 10 .^ log_σ
	y_pred = model_log_Σ(log_R; log_M_s=log_M_s, log_R_h=log_R_h, n=n, log_Sigma_bg=log_Sigma_bg)
	
	log_Σ ~ MvNormal(y_pred, diagm(log_Σ_err.^2) + Eye(length(log_Σ)) * σ .^2)
end

# ╔═╡ 29ee50a7-9287-46e7-a804-1aada9cd88b8
sersic_model = DensityModel(log_Σ_sersic, sersic_turing_model, [:n, :log_M_s, :log_R_h, :log_Sigma_bg])

# ╔═╡ dd327a8e-a83d-4423-b507-54033c7934c4
mcmc_fit_sersic = sample_model(sersic_model)

# ╔═╡ efad6b83-73a5-4aea-944b-ddc1942ac71a
plot_chains(mcmc_fit_sersic.samples)

# ╔═╡ ebc0cbc2-d9f3-4cab-a3ee-6a21763f4d20
pairplot(mcmc_fit_sersic.samples)

# ╔═╡ 61b28598-11fd-4ef3-bb7b-527482815eaa
let
	fig = plot_samples_residuals(mcmc_fit_sersic)

	@savefig "sersic_density_mcmc_fits"
	
	fig
end

# ╔═╡ 85731957-d441-47ad-b67b-ef1ea5b56c77
ml_fit_sersic = fit_profile(prof, LilGuys.Sersic, ["n"], [0.5, -0.0, 1.0])

# ╔═╡ 6d7a1522-e61e-4440-a2ff-64e625b15906
plot_Σ_fit_res(prof, ml_fit_sersic, title="Sérsic", nf=3)

# ╔═╡ 5a54ce11-92b3-455d-a3b3-b147618585ae
md"""
## Exp2D
"""

# ╔═╡ cdc4390f-1efb-4550-b325-5c743890a5c6
function log_Σ_exp2d(log_r; log_M, log_R_s)
	r = 10 .^ log_r
	prof = LilGuys.Exp2D(
		R_s=10 .^ log_R_s, 
		M = 10 .^ log_M
	)
	
	Σ = LilGuys.calc_Σ.(prof, r)
	y = @. log10(Σ)

	return y
end

# ╔═╡ 7f367da6-4b49-4e74-9a93-d9fdf517b2bd
@model function turing_model_exp2d(exp2d_log_Σ, log_r, log_Σ, log_Σ_err)
	log_R_s ~ Normal(0.0, 1.0)
	log_M ~ Normal(4.0, 2.0)

	log_σ ~ Normal(-1, 1)
	
	σ = 10 .^ log_σ
	y_pred = exp2d_log_Σ(log_r; log_M=log_M, log_R_s=log_R_s)
	
	log_Σ ~ MvNormal(y_pred, diagm(log_Σ_err.^2) + Eye(length(log_Σ)) * σ .^2)
end

# ╔═╡ d7641f55-8f5f-44b7-81ba-16a7994242fc
model_exp2d = DensityModel(log_Σ_exp2d, turing_model_exp2d, [:log_M, :log_R_s])

# ╔═╡ 1ad2610a-7f62-4dbc-a22e-7209f6c3e881
mcmc_fit_exp2d = sample_model(model_exp2d)

# ╔═╡ 7239698d-f802-4be3-9c2e-c7102858998a
plot_chains(mcmc_fit_exp2d.samples)

# ╔═╡ d634a7c8-b35c-4237-95f2-6fb52c666b2f
pairplot(mcmc_fit_exp2d.samples)

# ╔═╡ 31725ed6-3e11-4526-9e99-d3ffda80f859
let
	fig = plot_samples_residuals(mcmc_fit_exp2d)

	@savefig "exp2d"
	fig
end

# ╔═╡ 415397c4-f4b0-4e21-9df4-cf8faf761648
ml_fit_exp2d  = fit_profile(prof, LilGuys.Exp2D, [], [0.0, 0.0])

# ╔═╡ 76326a80-b709-4f96-a379-8aaea4d3af44
plot_Σ_fit_res(prof, ml_fit_exp2d, title="Exp2D")

# ╔═╡ c2301858-2be3-4945-bef2-4525cd123030
md"""
### Inner exp2d
"""

# ╔═╡ 19f1f280-b36c-43e2-971b-2b2d3bdc7a4f
log_R_max_exp = let
	ml_fit_exp2d_iter = []

	ml_fit_exp2d_old = ml_fit_exp2d
	log_R_max = log10(3) + ml_fit_exp2d_old["log_R_s"]

	for i in 1:10
		log_R_max_old = log_R_max
		ml_fit_exp2d_new  = fit_profile(filter_prof(prof, log_R_max), LilGuys.Exp2D, [], [0.0, 0.0])
		ml_fit_exp2d_old = ml_fit_exp2d_new
		push!(ml_fit_exp2d_iter, ml_fit_exp2d_new)

		log_R_max = log10(3) + ml_fit_exp2d_old["log_R_s"]

		if log_R_max ≈ log_R_max_old rtol=1e-1
			@info "converged by iteration $i"
			break
		else
			@info "R max $log_R_max"
		end
	end

	ml_fit_exp2d
	log_R_max
end

# ╔═╡ 48f59ed4-209b-45c5-857f-3c309141e8a7
let
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_R_label,
		ylabel = log_sigma_label,
	)

	LilGuys.plot_density_prof!(ax, prof_in, label="input")
	scatter!(prof.log_R, prof.log_Sigma, color=COLORS[2], label="mcmc input")
	prof2 = filter_prof(prof, log_R_max_exp)
	scatter!(prof2.log_R, prof2.log_Sigma, color=COLORS[3], label="exp input")

	axislegend(position=:lb)
	fig
end

# ╔═╡ 6922839a-1b5d-4a5c-8f68-50cd57d986a8
model_exp2d_inner = DensityModel(log_Σ_exp2d, turing_model_exp2d, [:log_M, :log_R_s])

# ╔═╡ bdf90521-e8f3-40e5-a61f-ab8503cce9fe
mcmc_fit_exp2d_inner = sample_model(model_exp2d_inner, prof=filter_prof(prof, log_R_max_exp))

# ╔═╡ 8e81babf-1963-4187-a8c5-0b1e1a325ab3
plot_chains(mcmc_fit_exp2d_inner.samples)

# ╔═╡ 6518f711-4912-421f-b169-7b7d4fb365e5
pairplot(mcmc_fit_exp2d_inner.samples)

# ╔═╡ 44500a4e-9be6-486b-94e7-364b1c64749a
let
	fig = plot_samples_residuals(mcmc_fit_exp2d_inner)

	@savefig "exp2d_inner"
	fig
end

# ╔═╡ bde4f544-7b5b-4224-9ce3-090fa7fe8311
md"""
## Plummer
"""

# ╔═╡ fbd8d5b6-3d48-4525-b99d-483c946cc6a3
function log_Σ_plummer(log_r; log_M, log_R_s)
	r = 10 .^ log_r
	prof = LilGuys.Plummer(
		r_s=10 .^ log_R_s, 
		M = 10 .^ log_M
	)
	
	Σ = LilGuys.calc_Σ.(prof, r)
	y = @. log10(Σ)

	return y
end

# ╔═╡ 81e50d1f-609b-46f5-b9e2-873dfeba1b45
@model function turing_model_plummer(log_Σ_plummer, log_r, log_Σ, log_Σ_err)
	log_R_s ~ Normal(0.0, 1.0)
	log_M ~ Normal(4.0, 2.0)

	log_σ ~ Normal(-1, 1)
	
	σ = 10 .^ log_σ
	y_pred = log_Σ_plummer(log_r; log_M=log_M, log_R_s=log_R_s)
	
	log_Σ ~ MvNormal(y_pred, diagm(log_Σ_err.^2) + Eye(length(log_Σ)) * σ .^2)
end

# ╔═╡ 2b71aaa2-30a5-4c25-aa3a-3da5fdc75f10
model_plummer = DensityModel(log_Σ_plummer, turing_model_plummer, [:log_M, :log_R_s])

# ╔═╡ ad16e44b-4f56-4ff5-a85e-f619cbf0c24f
mcmc_fit_plummer = sample_model(model_plummer)

# ╔═╡ e8168f09-b8d1-42fc-93eb-16063417795c
plot_chains(mcmc_fit_plummer.samples)

# ╔═╡ 8ae13c5e-0bd0-4a37-977a-92b4730be75a
pairplot(mcmc_fit_plummer.samples)

# ╔═╡ 7d47a188-b26a-4beb-8bca-d6761bf21ee6
ml_fit_plummer = fit_profile(prof, LilGuys.Plummer, [], [1.0, 0.0])

# ╔═╡ 3407d0ad-856e-415b-948a-3592b3f6b28c
let 
	fig = plot_samples_residuals(mcmc_fit_plummer)

	@savefig "plummer_density_mcmc_fits"

	fig
end

# ╔═╡ 22b49e19-30b4-40dd-b641-6a8f477aa541
plot_Σ_fit_res(prof, ml_fit_plummer, title="Plummer")

# ╔═╡ 9be6bf5f-aa72-4f44-8a7d-42ca74e017f4
md"""
## King
"""

# ╔═╡ 1c1ff4cb-69b6-48a4-ab8d-df6b403ecab1
function log_Σ_king(log_r; log_M, log_R_s, c)
	r = 10 .^ log_r
	prof = LilGuys.KingProfile(
		R_s=10 .^ log_R_s, 
		M = 10 .^ log_M, 
		c=c
	)
	
	Σ = LilGuys.calc_Σ.(prof, r)
	y = @. log10(Σ)

	return y
end

# ╔═╡ 9a8b3d0b-783c-451d-811c-e25846d840a7
@model function turing_model_king(log_Σ_king, log_r, log_Σ, log_Σ_err)
	log_R_s ~ Normal(0.0, 1.0)
	log_M ~ Normal(4.0, 2.0)
	c ~ truncated(Gamma(10), 1, Inf)

	log_σ ~ Normal(-1, 1)
	
	σ = 10 .^ log_σ
	y_pred = log_Σ_king(log_r; log_M=log_M, log_R_s=log_R_s, c=c)
	
	log_Σ ~ MvNormal(y_pred, diagm(log_Σ_err.^2) + Eye(length(log_Σ)) * σ .^2)
end

# ╔═╡ 4a0a40cd-f408-4969-87be-6d0eae29c09a
model_king = DensityModel(log_Σ_king, turing_model_king, [:log_M, :log_R_s, :c])

# ╔═╡ 8e95dc0d-d020-4062-8a88-908055783bef
mcmc_fit_king = sample_model(model_king)

# ╔═╡ 288e87ab-b2da-4e1a-97fc-84c9841669fa
plot_chains(mcmc_fit_king.samples)

# ╔═╡ 4b316372-11e4-42dc-b12f-e74637a35824


# ╔═╡ 4d70a044-7954-48df-a8b2-07ac930ee0bb
ml_fit_king = fit_profile(prof, LilGuys.KingProfile, ["c"], [2.0, 2.0, 3.0], M=1)

# ╔═╡ 942a9c71-e568-4828-bfd5-6072d587bf47
pairplot(mcmc_fit_king.samples)

# ╔═╡ 6473e3e9-393b-4952-8600-a16506cd0537
let 
	fig = plot_samples_residuals(mcmc_fit_king)

	@savefig "king_density_mcmc_fits"

	fig
end

# ╔═╡ 9840c373-dbf2-4d26-85b5-15ecbbb2c657
plot_Σ_fit_res(prof, ml_fit_king, title="King")

# ╔═╡ de7eb968-6616-491e-b258-d3f8be1134e8
md"""
# Old Fits
"""

# ╔═╡ 6f101561-2262-4012-9e17-aae2bd873d95
all_fits = OrderedDict(
	"sersic" => mcmc_fit_sersic,
	"exp2d" => mcmc_fit_exp2d,
	"exp2d_inner" => mcmc_fit_exp2d_inner,
	"plummer" => mcmc_fit_plummer,
	"king" => mcmc_fit_king,
)

# ╔═╡ a406d50a-468e-4b8a-ac0d-2bc4ea8e4251
begin
	derived_props = OrderedDict()

	for (name, fit) in all_fits
		for argname in [fit.model.argnames; :log_σ]
			i = findfirst(fit.summary.parameters .== argname)

			derived_props["$(argname)_$(name)"] = fit.summary.median[i]
			derived_props["$(argname)_$(name)_em"] = fit.summary.err_low[i]
			derived_props["$(argname)_$(name)_ep"] = fit.summary.err_high[i]
		end
		derived_props["chi2_$name"] = fit.χ2
		derived_props["lp_$name"] = fit.lp
		derived_props["num_params_$name"] = length(fit.model.argnames)
	end

	derived_props
end

# ╔═╡ 84fe5bae-c470-4df8-86fa-51cdf4300f01
import TOML

# ╔═╡ 456d26b3-6ec0-4fdc-95cc-78d35cc2e7ca
open(joinpath(galaxyname, "density_profiles/$(profilename)_density_fits.toml"), "w") do f
	TOML.print(f, derived_props)
end

# ╔═╡ Cell order:
# ╟─27bd8a46-ee56-11ef-22e4-afe88d73af2d
# ╟─9d3245ef-f774-4add-9b1d-60bb69b04c8d
# ╠═72524f04-bd5b-4f4b-a41d-8c038c180ce0
# ╠═48f59ed4-209b-45c5-857f-3c309141e8a7
# ╟─bd933a63-45e8-4e6c-9a0a-5d20c7d55eaf
# ╠═fa54e7aa-f477-45db-8dd9-56bd6b367604
# ╠═6cf3cb83-317e-40fa-be47-7a4e31a539ec
# ╠═cd1e43e6-65e1-4efe-9290-62633ee340a9
# ╠═51d52048-cd18-41ff-ac2b-41dfc097de16
# ╠═99315b4c-37c1-4287-bd71-c3339c8e76ad
# ╠═c3d9c4f3-6891-46e1-822e-5c217450f59b
# ╠═aae83e18-e739-4421-a763-37590cb7869a
# ╠═60763b01-9c7b-4153-9b94-617820517927
# ╠═37ee41d1-7b98-44f7-82ec-5b496622931b
# ╠═d0112dc0-8cc1-4fdf-b30d-a3b89660d9b4
# ╠═56af9242-91cc-46db-b532-c127450adfeb
# ╠═20b8c05c-91a6-4dc2-9558-1b02be4983f2
# ╠═1238d232-a2c0-44e5-936b-62fd9137d552
# ╠═25c3b74f-d8fa-4f4c-ae7c-71b137fb2ce7
# ╠═62cce84c-e985-4bd2-8c2d-ab15ed239b99
# ╠═11b84ff0-7f35-4716-ae30-e05a1c5f88ba
# ╠═3eb9992d-425d-4bec-94d2-da4f8dc59a8d
# ╠═db8a3705-f4ea-4588-9885-0e03cd42c47b
# ╟─8fdc3d18-6a3a-42eb-b0cb-d5cae1cd65d9
# ╠═20e017f2-c90b-46a5-a6f1-8e4aeb5ffde8
# ╠═249fed2c-e86a-47ab-a7a6-06f7c1edc01a
# ╠═44a893d6-9a3d-4959-b76f-274635198aa2
# ╠═a0cb6ecd-59ab-4a7a-b95a-d55bfdd2d93b
# ╠═995a9a15-38d0-47c3-ba12-58ba8209497c
# ╠═5026118a-85f6-4f71-bd0f-d25c4916c834
# ╟─30154472-316b-4489-b1e5-f503447928cb
# ╟─4970ae98-a08d-4308-b67e-e9306d4193a3
# ╠═8d90dd5c-350a-45e6-b04c-9b2f85cf700b
# ╠═aa22fa06-8167-4d6e-aa05-399800370c1f
# ╠═7ff0bee2-3622-4573-81d5-8b8fa96cf461
# ╠═7a274a8d-f6eb-4d9b-97c5-0372dc9520fe
# ╟─f91ae13d-4533-45aa-b835-d53fab594360
# ╠═3abb1fb9-1e3a-4251-abfd-959d1faa0f87
# ╠═95bd86c2-0222-4ef9-a8ef-d26e023fc171
# ╠═47ab55c6-e202-47a0-b1d5-684c9a1fef5d
# ╠═3981b6f9-7ace-4de9-b72d-10168512689c
# ╠═3d7e4658-97d5-4aff-baff-2ee87c688871
# ╠═5e31ad85-fbc5-459c-b6d4-b3493d3fd29f
# ╠═9a02adf8-9b90-4b4b-afcb-3ec5170f538b
# ╠═6c160b4e-177e-4b9e-9db4-3e2aa5465ce7
# ╠═776d6ffd-9598-45cc-b8d9-6c6b4025583a
# ╟─6c59342a-36b8-458b-8b17-92b0f18afb7d
# ╠═69b919f7-5985-4835-950a-ab7ad998d3fd
# ╠═071356e5-53c4-462e-8448-7e5480dbfc0e
# ╠═566f2f87-a840-416e-b71e-bcfbb747e39b
# ╠═2bc41ec4-6ced-474e-92e8-334551b38431
# ╠═e31c35b1-9b37-4c95-8344-6177b685d5ab
# ╠═7779ec2e-1053-47e9-8d51-897d5063bb85
# ╟─7115482e-0cc9-40ef-b63e-0ebad62bb949
# ╟─a9e5aa14-6a11-436c-97ca-7afd5fd35478
# ╠═a4d8323a-c668-452b-a958-17982177cac6
# ╠═2a3adaf5-2a0d-4674-a160-a7844b7f1463
# ╠═29ee50a7-9287-46e7-a804-1aada9cd88b8
# ╠═dd327a8e-a83d-4423-b507-54033c7934c4
# ╠═efad6b83-73a5-4aea-944b-ddc1942ac71a
# ╠═ebc0cbc2-d9f3-4cab-a3ee-6a21763f4d20
# ╠═61b28598-11fd-4ef3-bb7b-527482815eaa
# ╠═85731957-d441-47ad-b67b-ef1ea5b56c77
# ╠═6d7a1522-e61e-4440-a2ff-64e625b15906
# ╟─5a54ce11-92b3-455d-a3b3-b147618585ae
# ╠═cdc4390f-1efb-4550-b325-5c743890a5c6
# ╠═7f367da6-4b49-4e74-9a93-d9fdf517b2bd
# ╠═d7641f55-8f5f-44b7-81ba-16a7994242fc
# ╠═1ad2610a-7f62-4dbc-a22e-7209f6c3e881
# ╠═7239698d-f802-4be3-9c2e-c7102858998a
# ╠═d634a7c8-b35c-4237-95f2-6fb52c666b2f
# ╠═31725ed6-3e11-4526-9e99-d3ffda80f859
# ╠═415397c4-f4b0-4e21-9df4-cf8faf761648
# ╠═76326a80-b709-4f96-a379-8aaea4d3af44
# ╠═c2301858-2be3-4945-bef2-4525cd123030
# ╠═19f1f280-b36c-43e2-971b-2b2d3bdc7a4f
# ╠═6922839a-1b5d-4a5c-8f68-50cd57d986a8
# ╠═bdf90521-e8f3-40e5-a61f-ab8503cce9fe
# ╠═8e81babf-1963-4187-a8c5-0b1e1a325ab3
# ╠═6518f711-4912-421f-b169-7b7d4fb365e5
# ╠═44500a4e-9be6-486b-94e7-364b1c64749a
# ╟─bde4f544-7b5b-4224-9ce3-090fa7fe8311
# ╠═fbd8d5b6-3d48-4525-b99d-483c946cc6a3
# ╠═81e50d1f-609b-46f5-b9e2-873dfeba1b45
# ╠═2b71aaa2-30a5-4c25-aa3a-3da5fdc75f10
# ╠═ad16e44b-4f56-4ff5-a85e-f619cbf0c24f
# ╠═e8168f09-b8d1-42fc-93eb-16063417795c
# ╠═8ae13c5e-0bd0-4a37-977a-92b4730be75a
# ╠═7d47a188-b26a-4beb-8bca-d6761bf21ee6
# ╟─3407d0ad-856e-415b-948a-3592b3f6b28c
# ╠═22b49e19-30b4-40dd-b641-6a8f477aa541
# ╟─9be6bf5f-aa72-4f44-8a7d-42ca74e017f4
# ╠═1c1ff4cb-69b6-48a4-ab8d-df6b403ecab1
# ╠═9a8b3d0b-783c-451d-811c-e25846d840a7
# ╠═4a0a40cd-f408-4969-87be-6d0eae29c09a
# ╠═8e95dc0d-d020-4062-8a88-908055783bef
# ╠═288e87ab-b2da-4e1a-97fc-84c9841669fa
# ╠═4b316372-11e4-42dc-b12f-e74637a35824
# ╠═4d70a044-7954-48df-a8b2-07ac930ee0bb
# ╠═942a9c71-e568-4828-bfd5-6072d587bf47
# ╠═6473e3e9-393b-4952-8600-a16506cd0537
# ╠═9840c373-dbf2-4d26-85b5-15ecbbb2c657
# ╟─de7eb968-6616-491e-b258-d3f8be1134e8
# ╠═6f101561-2262-4012-9e17-aae2bd873d95
# ╠═a406d50a-468e-4b8a-ac0d-2bc4ea8e4251
# ╠═84fe5bae-c470-4df8-86fa-51cdf4300f01
# ╠═456d26b3-6ec0-4fdc-95cc-78d35cc2e7ca
