### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ bfeedd02-4339-40b8-bf68-1ccdbeaa5245
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames, CSV
	
	using GLMakie
	using Measurements
	#using KernelDensity
	
	import SciPy
	using QuadGK
	
	import LinearAlgebra: diag
	
	import LilGuys as lguys
	using Arya
	
	using JSON

	import NaNMath as nm
end

# ╔═╡ 0d71de16-50d7-43da-b1fb-a477b949e15e
md"""
Given a set of density (histogram) observations, what profile is consistant?
"""

# ╔═╡ f93365e6-971d-4321-9d91-44e9e86610cb
md"""
# loading data...

We expect data as a
"""

# ╔═╡ d0992dc9-08f1-487a-a96a-90996f29cefd
name = "sculptor/fiducial"

# ╔═╡ aa23a8ab-3cff-400b-bbc9-592183f2e695
profile_name = name *  "_profile.json"

# ╔═╡ 88f02d81-927c-405a-ada4-064157c7dbf0
open(profile_name) do f
	global profile
	profile = JSON.parse(f)
end

# ╔═╡ f890490e-17ef-4e78-a18a-437c90724f86
md"""
# Fitting
"""

# ╔═╡ 333c27b4-1da5-4745-942a-961202399f6d
log_r_label = L"\log r / \mathrm{arcmin}"

# ╔═╡ 46d3af07-cf51-40eb-b129-8f3b8e0fdeed
function plot_Σ_fit_res(obs, pred, res)
    fig = Figure()
    ax = Axis(fig[1, 1], 
        ylabel=L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}")
	
    errscatter!(ax, obs["log_r"], obs["log_Sigma"], yerr=obs["log_Sigma_err"])


    lines!(ax, pred.log_r, log10.(pred.Σ), color=COLORS[2])
    
    ax2 = Axis(fig[2, 1],
        ylabel=L"\delta\log\Sigma", 
    	xlabel=log_r_label,
		limits = (nothing, (-1, 1))
	)

#     p2 = plot(ylabel=L"\Delta \log\Sigma", xlabel=log_r_label, ylim=(-2, 2))

	y = res
    errscatter!(ax2, obs["log_r"], y, yerr=obs["log_Sigma_err"], label="")


    hlines!(0, color=:black)
    
    rowsize!(fig.layout, 2, Relative(1/4))

    linkxaxes!(ax, ax2)
    hidexdecorations!(ax, grid=false)
#     return plot(p1, p2, layout=grid(2, 1, heights=(0.8, 0.2)), link=:x, bottom_margin=[-5Plots.mm 0Plots.mm])
    return fig
end

# ╔═╡ 6b909975-cdba-4961-bc0f-842c68f33ef9
function calc_Γ(log_rs, Σs, step=1)
	dx = lguys.gradient(log_rs)
	dρ = lguys.gradient(log10.(Σs))

	return dρ ./ dx #lguys.gradient(log10.(ρs), log_rs)
end

# ╔═╡ f6715279-c048-4260-b30a-8e0a4f7c5af5
function predict_properties(Σ_model; N=10_000, log_r_min=-2, log_r_max=2)
    log_r_bins = LinRange(log_r_min, log_r_max, 1000)
    log_r = lguys.midpoint(log_r_bins)
    r = 10 .^  log_r
    r_bins = 10 .^ log_r_bins
    
    Σ = Σ_model.(r)

    Γ = calc_Γ(log_r, Σ)
    M_in = [quadgk(rrr->2π*rrr*Σ_model(rrr), 0, rr)[1] for rr in r]
    Σ_m = M_in ./ (π * r .^ 2)
    Γ_max = 2*(1 .- Σ ./ Σ_m)
    counts =  Σ .* (2π *  r .* diff(r_bins) )
    
    return (log_r=log_r, log_r_bins=log_r_bins, counts=counts, Σ=Σ, Σ_m=Σ_m, Γ=Γ, Γ_max=Γ_max, M_in=M_in)

end

# ╔═╡ 347aee22-17bf-11ef-196f-0146bd88f688
function fit_profile(obs; r_max=Inf, N=10_000, profile=lguys.Exp2D, p0=[2, 0.3])
    r_val = 10 .^ obs["log_r"]
    log_Σ = obs["log_Sigma"] .± obs["log_Sigma_err"]
    filt = r_val .< r_max
    filt .&= map(x->isfinite(x), log_Σ)
    filt .&= @. !isnan(log_Σ)
    
    r_val = r_val[filt]
    log_Σ = log_Σ[filt]


    log_Σ_val = [s.val for s in log_Σ]
    log_Σ_e = [s.err for s in log_Σ]

	log_Σ_exp(r, popt...) = nm.log10.(lguys.calc_Σ.(profile(popt...), r))
	popt, covt = SciPy.optimize.curve_fit(log_Σ_exp, r_val, log_Σ_val, 
        sigma=log_Σ_e, p0=p0)
    
    popt_p = popt .± sqrt.(diag(covt))

	Σ_pred(r) = 10 .^ log_Σ_exp(r, popt...)
    props = predict_properties(Σ_pred, 
        N=N, log_r_min=obs["log_r_bins"][1], log_r_max=obs["log_r_bins"][end])
    
    log_Σ_pred = log_Σ_exp.(r_val, popt...)
    log_Σ_res = value.(log_Σ).- log_Σ_pred
    return popt_p, props, log_Σ_res
end

# ╔═╡ 3a0f2b8b-bb70-41bf-bcd7-aa241d4b7bdf
popt, pred, res = fit_profile(profile, p0=[2, 7])

# ╔═╡ 5cf332e8-47e4-4e9e-af94-199331a7aa3d
md"""
# Visual inspection
"""

# ╔═╡ d6b7764b-65b1-4a71-9759-0c8b0b8672b6
let
	fig, ax, p = errscatter(profile["log_r"], profile["counts"], yerr=sqrt.(profile["counts"]))

	ax.xlabel = log_r_label
	ax.ylabel = "count / bin"

	ax.yscale=Makie.pseudolog10

	fig
end

# ╔═╡ 98b776e5-824b-4db5-8455-b2433fba22b1
plot_Σ_fit_res(profile, pred, res)

# ╔═╡ 3c133454-d1c2-4aff-a27f-c3368bf06480
let
	fig = Figure(size=(700, 300))
	ax = Axis(fig[1, 1], limits=((-0.8, 3), nothing),
		xlabel=log_r_label,
		ylabel=L"\Gamma"
	)

	
	errscatter!(value.(profile["log_r"]), profile["Gamma"], yerr=profile["Gamma_err"])
	
	lines!(pred.log_r, pred.Γ, label="2D Exp", color=COLORS[2])


	ax_lin = Axis(fig[1, 2],
		xlabel="r / arcmin",
		yticklabelsvisible=false,
	)

	
	errscatter!(value.(10 .^ profile["log_r"]), profile["Gamma"], yerr=profile["Gamma_err"])
	
	lines!(10 .^ pred.log_r, pred.Γ, label="2D Exp", color=COLORS[2])

	linkyaxes!(ax, ax_lin)
	
	axislegend(ax)
	fig
end

# ╔═╡ ab2ffe07-0f81-4254-b1a0-f4e46532a077
let
	fig = Figure()
	ax = Axis(fig[1, 1], limits=((-0.8, 3), nothing),
		xlabel=log_r_label,
		ylabel=L"\Gamma_\textrm{max}"
	)

	
	errscatter!(value.(profile["log_r"]), profile["Gamma_max"], yerr=profile["Gamma_max_err"])
	
	lines!(pred.log_r, pred.Γ_max, label="2D Exp", color=COLORS[2])


	axislegend(ax)
	fig
end

# ╔═╡ e5cb1a72-90d9-45be-9450-a86c640b9420
md"""
# Comparting functional forms
"""

# ╔═╡ 0b79e7ec-b595-4dd0-9cfc-ab2eb0db9a12
let 
	popt, pred, res = fit_profile(profile, p0=[2, 7], profile=lguys.Exp3D)

	global popt_3d_exp = popt
	
	plot_Σ_fit_res(profile, pred, res)

end

# ╔═╡ 93cd5808-0e20-475a-8001-4e775d4ab4e7
let 
	popt, pred, res = fit_profile(profile, p0=[2, 7], profile=lguys.Exp2D)

	global popt_2d_exp = popt
	
	plot_Σ_fit_res(profile, pred, res)

end

# ╔═╡ 824fc065-8ad8-408c-a4f7-0c464ed7fa13
let 
	popt, pred, res = fit_profile(profile, profile=lguys.KingProfile, p0=[1, 35, 600])

	global popt_king = popt
	
	plot_Σ_fit_res(profile, pred, res)

end

# ╔═╡ fa96cf7f-248c-4765-aa3c-0ca9f656644d
popt_king

# ╔═╡ 1be3fa21-5945-4904-a23a-12c15cc4a485
let 
	popt, pred, res = fit_profile(profile, profile=lguys.LogCusp2D, p0=[1, 35])

	global popt_cusp = popt
	
	plot_Σ_fit_res(profile, pred, res)

end

# ╔═╡ 2f7bf886-5b1e-44d1-a16b-1d6214405a5f
md"""
# Saving fit parameters to json
"""

# ╔═╡ Cell order:
# ╠═0d71de16-50d7-43da-b1fb-a477b949e15e
# ╠═bfeedd02-4339-40b8-bf68-1ccdbeaa5245
# ╟─f93365e6-971d-4321-9d91-44e9e86610cb
# ╠═d0992dc9-08f1-487a-a96a-90996f29cefd
# ╠═aa23a8ab-3cff-400b-bbc9-592183f2e695
# ╠═88f02d81-927c-405a-ada4-064157c7dbf0
# ╠═f890490e-17ef-4e78-a18a-437c90724f86
# ╠═347aee22-17bf-11ef-196f-0146bd88f688
# ╠═f6715279-c048-4260-b30a-8e0a4f7c5af5
# ╠═46d3af07-cf51-40eb-b129-8f3b8e0fdeed
# ╠═333c27b4-1da5-4745-942a-961202399f6d
# ╠═6b909975-cdba-4961-bc0f-842c68f33ef9
# ╠═3a0f2b8b-bb70-41bf-bcd7-aa241d4b7bdf
# ╠═5cf332e8-47e4-4e9e-af94-199331a7aa3d
# ╠═d6b7764b-65b1-4a71-9759-0c8b0b8672b6
# ╠═98b776e5-824b-4db5-8455-b2433fba22b1
# ╠═3c133454-d1c2-4aff-a27f-c3368bf06480
# ╠═ab2ffe07-0f81-4254-b1a0-f4e46532a077
# ╟─e5cb1a72-90d9-45be-9450-a86c640b9420
# ╠═0b79e7ec-b595-4dd0-9cfc-ab2eb0db9a12
# ╠═93cd5808-0e20-475a-8001-4e775d4ab4e7
# ╠═824fc065-8ad8-408c-a4f7-0c464ed7fa13
# ╠═fa96cf7f-248c-4765-aa3c-0ca9f656644d
# ╠═1be3fa21-5945-4904-a23a-12c15cc4a485
# ╟─2f7bf886-5b1e-44d1-a16b-1d6214405a5f
