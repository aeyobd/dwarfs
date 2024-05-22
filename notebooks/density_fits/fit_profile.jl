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

# ╔═╡ fc892792-3f09-4bcc-84d1-eb9ad4f17b90
using Tables

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
name = "Scl"

# ╔═╡ aa23a8ab-3cff-400b-bbc9-592183f2e695
profile_name = name *  "_profile.fits"

# ╔═╡ b40e74a6-2546-4e23-9318-138e1599024e
f = FITS(profile_name)

# ╔═╡ 1a4538d3-31f4-4267-9344-68662d3445b7
df = DataFrame(f[2])

# ╔═╡ dd43be76-f3be-49a4-abb8-72f1ac9f4491
obs = Dict(name => read(f[2], name) for name in names(df))

# ╔═╡ ff486f44-b4e3-4e9e-b276-472c89b95761
read(obs[end, :])

# ╔═╡ 74520378-dabb-4539-ba1a-cc1aff23cae9
obs["log_r"] = obs["log_r"][1:end-1]

# ╔═╡ 39faf220-f090-4f31-a944-5566602b1fd8
obs["log_Sigma"] = obs["log_Sigma"][1:end-1]

# ╔═╡ f6af24d3-5d9d-488a-869f-35301beeb95b
obs["log_Sigma_err"] = obs["log_Sigma_err"][1:end-1]

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
    y = obs.log_Sigma  
	
    errscatter!(ax, obs["log_r"], obs["log_Sigma"], yerr=obs["log_Sigma_err"])


    lines!(ax, pred.log_r, log10.(pred.Σ), color=COLORS[2])
    
    ax2 = Axis(fig[2, 1],
        ylabel=L"\delta\log\Sigma", 
    	xlabel=log_r_label,
		limits = (nothing, (-1, 1))
	)

#     p2 = plot(ylabel=L"\Delta \log\Sigma", xlabel=log_r_label, ylim=(-2, 2))

	y = res
    errscatter!(ax2, obs.log_r, y, yerr=obs["log_Sigma_err"], label="")


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
popt, pred, res = fit_profile(obs, p0=[2, 7])

# ╔═╡ 1540f3ca-3279-47d8-876c-5544ed7fb759
obs["log_Sigma"]

# ╔═╡ 0efcb02f-29c6-4276-a9a4-3fb735eca3a5
obs

# ╔═╡ 5cf332e8-47e4-4e9e-af94-199331a7aa3d
md"""
# Visual inspection
"""

# ╔═╡ d6b7764b-65b1-4a71-9759-0c8b0b8672b6
let
	fig, ax, p = errscatter(obs["log_r"], obs["counts"], yerr=sqrt.(obs["counts"]))

	ax.xlabel = log_r_label
	ax.ylabel = "count / bin"

	ax.yscale=Makie.pseudolog10

	fig
end

# ╔═╡ 3c133454-d1c2-4aff-a27f-c3368bf06480
let
	fig = Figure()
	ax = Axis(fig[1, 1], limits=((-0.8, 2), (-6.5, -1.5)))
	errscatter!(value.(obs["log_r"]), obs["log_Sigma"], yerr=obs["log_Sigma_err"])
	
	lines!(pred.log_r, log10.(pred.Σ), label="2D Exp", color=COLORS[2])


	ax.xlabel = log_r_label
	ax.ylabel = L"\log\Sigma\quad [\textrm{fraction arcmin}^{-2}]"

	axislegend(ax)
	fig
end

# ╔═╡ 98b776e5-824b-4db5-8455-b2433fba22b1
plot_Σ_fit_res(obs, pred, res)

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
# ╠═ff486f44-b4e3-4e9e-b276-472c89b95761
# ╠═b40e74a6-2546-4e23-9318-138e1599024e
# ╠═1a4538d3-31f4-4267-9344-68662d3445b7
# ╠═fc892792-3f09-4bcc-84d1-eb9ad4f17b90
# ╠═dd43be76-f3be-49a4-abb8-72f1ac9f4491
# ╠═74520378-dabb-4539-ba1a-cc1aff23cae9
# ╠═39faf220-f090-4f31-a944-5566602b1fd8
# ╠═f6af24d3-5d9d-488a-869f-35301beeb95b
# ╠═f890490e-17ef-4e78-a18a-437c90724f86
# ╠═347aee22-17bf-11ef-196f-0146bd88f688
# ╠═f6715279-c048-4260-b30a-8e0a4f7c5af5
# ╠═46d3af07-cf51-40eb-b129-8f3b8e0fdeed
# ╠═333c27b4-1da5-4745-942a-961202399f6d
# ╠═6b909975-cdba-4961-bc0f-842c68f33ef9
# ╠═3a0f2b8b-bb70-41bf-bcd7-aa241d4b7bdf
# ╠═1540f3ca-3279-47d8-876c-5544ed7fb759
# ╠═0efcb02f-29c6-4276-a9a4-3fb735eca3a5
# ╠═5cf332e8-47e4-4e9e-af94-199331a7aa3d
# ╠═d6b7764b-65b1-4a71-9759-0c8b0b8672b6
# ╠═3c133454-d1c2-4aff-a27f-c3368bf06480
# ╠═98b776e5-824b-4db5-8455-b2433fba22b1
# ╟─2f7bf886-5b1e-44d1-a16b-1d6214405a5f
