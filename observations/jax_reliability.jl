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

# ╔═╡ bff50014-bfa9-11ee-33f0-0f67e543c2d4
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames 
	using CSV
	using CairoMakie

	import NaNMath as nm
	using KernelDensity
	using Measurements

end

# ╔═╡ 8acddd07-0eff-47a6-a364-f7f58680bd9c
begin 
	using Arya
	Arya.update_figsize!(7), Arya.update_fontsize!(24)
	update_theme!(
	CairoMakie = (; px_per_unit = 2.5,)
)
end

# ╔═╡ 8b3955e5-ab08-4d5a-abe5-8f40120ff8db
using PlutoUI

# ╔═╡ ae29bed0-6700-47f1-8952-35e867ce126b
using OrderedCollections

# ╔═╡ 47b8b3b0-0228-4f50-9da4-37d388ef9e9f
md"""
# Jensen et al. 2024 reliability

Some plots to understand the reliability of density profiles derived from the J+24 data sample.

We want to know if the choice of a spatial component affects the derived density profile and at what point these profile may become unreliable.
"""

# ╔═╡ 303ff946-ecdc-42c3-8289-c13e7a8cfe7c
import FileIO: FileIO, @format_str

# ╔═╡ eca9c1c5-e984-42d4-8854-b227fdec0a8a
md"""
galaxyname: $(@bind galaxyname confirm(TextField(default="galaxy")))
"""

# ╔═╡ 1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
begin 
	import LilGuys as lguys
	using LilGuys
	FIGDIR = joinpath(galaxyname, "figures")
	FIGSUFFIX = ".jax_observations"
end

# ╔═╡ da9ca7a7-18b8-49cb-a041-ab1c667920ff
import DensityEstimators: histogram2d

# ╔═╡ 489f6d21-9e9a-4b1e-b05c-c63a44ba1951
import StatsBase: percentile, mean

# ╔═╡ 0004f638-c57a-4dab-8b97-c77840cafbbf
import TOML

# ╔═╡ db1264b7-02c3-4b55-ae2b-9ce78aa1304a
import DensityEstimators: histogram

# ╔═╡ 81c677bf-8372-4fd3-9d60-af8377485f1c
import Statistics: cor

# ╔═╡ cf7aeb0a-d453-462a-b43d-a832567440fd
diverging_cmap = Reverse(:bluesreds)

# ╔═╡ fe7050d0-fede-44c0-8a8e-5d3dcaf0d1fb
CairoMakie.activate!(type="png")

# ╔═╡ ff41b9a8-e25c-4fc8-b00f-0d74ea6a1d31
FigAxis(px_)[1]

# ╔═╡ a955fb87-3674-48dc-8dbd-2c61926d195d
fig 

# ╔═╡ 2eb4aa78-0fea-460b-a18e-06a129c41504
md"""
# Inputs & Data loading
"""

# ╔═╡ 9ecf79a8-2ed3-40c6-b555-a102250ecbd4
observed_properties = TOML.parsefile(galaxyname * "/observed_properties.toml")

# ╔═╡ 4c8b9c30-f873-4756-a729-aa023b3a804e
readdir("$galaxyname/data")

# ╔═╡ 26cf1867-02be-4d36-8c35-6c58a1feca27
if isfile("$galaxyname/data/jensen+24_2c.fits")
	datafile = "$galaxyname/data/jensen+24_2c.fits"
else
	datafile = "$galaxyname/data/jensen+24_1c.fits"
end

# ╔═╡ ec227641-86e6-46b7-8019-9b02072ed9f7
begin 
	all_stars = lguys.read_fits(datafile)

	all_stars[:, :r_ell_old] = all_stars[:, :r_ell]
	all_stars.xi, all_stars.eta = lguys.to_tangent(all_stars.ra, all_stars.dec, observed_properties["ra"], observed_properties["dec"])

	all_stars.r_ell = 60lguys.calc_r_ell(all_stars.xi, all_stars.eta, observed_properties["ellipticity"], observed_properties["position_angle"])
	all_stars[!, :LL_S] = @. log10(all_stars.L_S_SAT / all_stars.L_S_BKD)
	all_stars[!, :LL_PM] = @. log10(all_stars.L_PM_SAT / all_stars.L_PM_BKD)
	all_stars[!, :LL_CMD] = @. log10(all_stars.L_CMD_SAT / all_stars.L_CMD_BKD)
	all_stars[!, :LL] = @. log10(all_stars.L_SAT / all_stars.L_BKD)

	all_stars[!, :LL_norm] = @. abs(all_stars.LL_S) + abs(all_stars.LL_PM) + abs(all_stars.LL_CMD)

	all_stars[!, :f_LL_S] = @. all_stars.LL_S / all_stars.LL
	all_stars[!, :f_LL_PM] = @. all_stars.LL_PM / all_stars.LL
	all_stars[!, :f_LL_CMD] = @. all_stars.LL_CMD / all_stars.LL

	all_stars[!, :f_LL_S_norm] = @. all_stars.LL_S / all_stars.LL_norm
	all_stars[!, :f_LL_PM_norm] = @. all_stars.LL_PM / all_stars.LL_norm
	all_stars[!, :f_LL_CMD_norm] = @. all_stars.LL_CMD / all_stars.LL_norm

	all_stars
end

# ╔═╡ 7ed70894-48a8-46af-8c78-8a75cb39d957
md"""
estimate of centring uncertanty in this dataset (standard error on mean and present difference between mean centre & assumed centre)

"""

# ╔═╡ 3abea230-2b95-40d3-9851-a91236f75c4a
md"""
each component of the cen_errs tuple below should be similar in magnitude
"""

# ╔═╡ 88fbdd09-30be-4fc3-95ae-acce6e0018e1
members = all_stars[all_stars.PSAT .> 0.2, :]

# ╔═╡ 223abc41-5158-49c2-96bf-df55b7be1114
cen = lguys.to_tangent(lguys.calc_centre2D(members.ra, members.dec, "mean")..., observed_properties["ra"], observed_properties["dec"])

# ╔═╡ 32e8e100-1fad-425f-b05c-2310e7ba559c
cen_errs =  (cen..., 
lguys.std(members.eta) / sqrt(length(members.eta)), 
lguys.std(members.xi) / sqrt(length(members.xi))
)

# ╔═╡ b1be5f97-0e61-4dda-87fd-7bb0f43147b6
cen_err = maximum(abs.(cen_errs))

# ╔═╡ 731ea468-5003-44e9-95b8-7fa7ef4b475b
Nmemb = size(members, 1)

# ╔═╡ 60d0e593-88fd-4b4c-9009-cc24a597c6d5
members_nospace = all_stars[all_stars.PSAT_NOSPACE .> 0.2, :]

# ╔═╡ 082a06dd-eeb5-4761-a233-1ee89e8cb819
best_stars = all_stars[all_stars.F_BEST .== 1.0, :]

# ╔═╡ a9d94121-ea6e-416a-bae8-aa93c16bde72
md"""
# Utilities
"""

# ╔═╡ f3d3f7b2-30c2-4f89-bef5-ac1a8720963f
logit(p) = log(p / (1-p))

# ╔═╡ 7cc17173-160c-45f1-a685-c79271679414
xi_label = L"\xi\,/\,degrees"

# ╔═╡ 965217b8-b2a5-485b-83de-cac887065b19
plot_labels = OrderedDict(
	:xi => L"$\xi$\,/\,degrees",
	:eta => L"$\eta$\,/\,degrees",
	:xi_am => L"$\xi$\,/\,arcmin",
	:eta_am => L"$\eta$\,/\,arcmin",
	:G => "G (mag)",
	:bp_rp => "BP – RP (mag)",
	:pmra => L"$\mu_{\alpha*}$ / mas yr$^{-1}$",
	:pmdec => L"$\mu_{\delta}$ / mas yr$^{-1}$",
	:r_ell => L"$r_\textrm{ell}$ / arcmin",
	:LL => L"\log\,\mathcal{L}_\textrm{sat} \, /\,\mathcal{L}_\textrm{MW}"
)

# ╔═╡ 2d4da56b-5d0a-49d3-83ae-a90f85192101
θmax = maximum(sqrt.(all_stars.xi.^2 .+ all_stars.eta .^ 2))

# ╔═╡ b0aec8b4-410f-4678-ad8c-cb055b743bd3
md"""
# Likelihood components
"""

# ╔═╡ 3c2ee4bb-d594-47e7-a7f3-5b7c9e38b7d4
md"""
## Likelihoods maps
"""

# ╔═╡ b191b7b8-28ed-4aee-b376-b2bdf22abc29
md"""
These maps are just for reference, to see how strongly each component of the likelihood varies in that parameter space.
"""

# ╔═╡ b06445ed-6290-42e8-bdb8-20f0372518c2
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi / degrees",
		ylabel = "eta / dgrees",
		aspect=DataAspect(),
		)
	
	p = scatter!((best_stars.xi), (best_stars.eta), color=(best_stars.L_S_BKD),
		colormap=:bluesreds,
		markersize=3
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_S_BKD", )
	
	fig
end

# ╔═╡ 2a8e13bd-b8ca-4cf3-b61b-096af218265b
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		)
	
	p = scatter!((best_stars.xi), (best_stars.eta), color=log10.(best_stars.L_S_SAT),
		markersize=3,
		#colorrange=(-1, 0)
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_S_SAT", )
	
	fig
end

# ╔═╡ 25dd2e38-1c78-4ef6-95ef-d38a0e80f0e4
let 
	f = Figure()
	ax =  Axis(f[1,1], yreversed=true,
	xlabel="bp - rp",
	ylabel = "G",
)
	LilGuys.hide_grid!(ax)

	
	h = scatter!(best_stars.bp_rp, best_stars.phot_g_mean_mag, color=log10.(best_stars.L_CMD_SAT), 
		markersize = 3
	)

	Colorbar(f[1, 2], h, label="log L_CMD_SAT")
	f
end

# ╔═╡ 1b63c1ee-0a46-49c6-880e-c1742055605e
let 
	f = Figure()
	ax =  Axis(f[1,1], yreversed=true,
	xlabel="bp - rp",
	ylabel = "G",
)

	
	h = scatter!(best_stars.bp_rp, best_stars.phot_g_mean_mag, color=log10.(best_stars.L_CMD_BKD), 
		markersize = 3
	)

	LilGuys.hide_grid!(ax)
	Colorbar(f[1, 2], h, label="log L_CMD_BKD")
	f
end

# ╔═╡ 643773f0-2b03-46f7-8df4-579b3c708909
let 
	f = Figure()
	ax =  Axis(f[1,1], yreversed=true,
		xlabel="bp - rp",
		ylabel = "G",
	)
	LilGuys.hide_grid!(ax)
	
	h = scatter!(best_stars.bp_rp, best_stars.phot_g_mean_mag, color=log10.(best_stars.L_CMD_SAT ./ best_stars.L_CMD_BKD), 
		colorrange=(-10, 1),
		markersize=3
	)

	Colorbar(f[1, 2], h, label="log L_CMD_SAT / L_CMD_BKD")
	f
end

# ╔═╡ cad88ab4-1ada-4cc3-8ef0-d4674500d57e
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "pmra",
		ylabel = "pmdec",
		aspect=DataAspect(),
		)
	
	p = scatter!((best_stars.pmra), (best_stars.pmdec), color=log10.(best_stars.L_PM_SAT),
		colorrange=(-30, 1),
		markersize=3
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_PM_SAT", )
	
	fig
end

# ╔═╡ 18252c3c-d4d6-46f2-8277-33a700ab4074
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "pmra",
		ylabel = "pmdec",
		aspect=DataAspect(),
		)
	
	p = scatter!((best_stars.pmra), (best_stars.pmdec), color=log10.(best_stars.L_PM_BKD),
		markersize=3
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_PM_BKD", )
	
	fig
end

# ╔═╡ d47a8d02-6c63-41fd-95b5-85550d1f5a85
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "pmra",
		ylabel = "pmdec",
		aspect=DataAspect(),
		)
	
	p = scatter!((best_stars.pmra), (best_stars.pmdec), color=log10.(best_stars.L_PM_SAT ./ best_stars.L_PM_BKD),
		markersize=3
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_PM_SAT / L_PM_BKD", )
	
	fig
end

# ╔═╡ b7e7e245-3732-453c-a2ec-3fe84c65c012
md"""
## Likelihood ratios
"""

# ╔═╡ 87f29963-2c7f-4d8f-8631-f70197d315c0
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		title = "LLR > -5"
		)

	df = best_stars[best_stars.LL .> -10, :]
	p = scatter!((df.xi), (df.eta), color=(df.LL_S) .- (df.LL_norm),
		markersize=3,
		colorrange=(-10, 0)
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_S / L_tot_norm", )
	
	fig
end

# ╔═╡ 72aa995b-8dd1-4d08-88c2-5fb527ff7f3e
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		title = "LLR > -10",
		)

	df = best_stars[best_stars.LL .> -10, :]
	p = scatter!((df.xi), (df.eta), color=(df.LL_S) .- (df.LL),
		markersize=3,
		colorrange=(-10, 10),
		colormap=diverging_cmap,
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L_S / L_tot", )
	
	fig
end

# ╔═╡ 6cdcd105-fa1d-48f5-9861-b722aafa28a3
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		)

	df = best_stars[best_stars.LL .> -10, :]
	p = scatter!((df.xi), (df.eta), color=(df.LL_PM) .- (df.LL),
		markersize=3,
		colorrange=(-10, 10),
		colormap=diverging_cmap,
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L pm / L tot", )
	
	fig
end

# ╔═╡ cf8e51a9-194c-4aa2-9766-39256a35ef02
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		)

	df = best_stars[best_stars.LL .> -5, :]
	p = scatter!((df.xi), (df.eta), color=(df.LL_PM) .- (df.LL_norm),
		markersize=3,
		colorrange=(-10, 0)
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L pm / L tot norm", )
	
	fig
end

# ╔═╡ e265ec56-d62d-4295-975f-98ba93cf5e16
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		)

	df = best_stars[best_stars.LL .> -10, :]
	p = scatter!((df.xi), (df.eta), color=(df.LL_CMD) .- (df.LL),
		markersize=3,
		colorrange=(-10, 10),
		colormap=diverging_cmap,
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L cmd / L tot", )
	
	fig
end

# ╔═╡ e2686165-6499-4558-9da9-91a78624b992
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "xi",
		ylabel = "eta",
		aspect=DataAspect(),
		)

	df = best_stars[best_stars.LL .> -5, :]
	p = scatter!((df.xi), (df.eta), color=(df.LL_CMD) .- (df.LL_norm),
		markersize=3,
		colorrange=(-10, 0)
	)
	LilGuys.hide_grid!(ax)

	Colorbar(fig[1,2], p, label="log L cmd / L tot norm", )
	
	fig
end

# ╔═╡ 1f2aa1ee-b395-4040-84cb-973d57b534d5
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "Log Likelihood",
		ylabel = "counts",
	)

	filt = isfinite.(best_stars.LL)

	stephist!(log10.(best_stars.L_S_SAT[filt]), 
			label="spatial")
	
	stephist!(log10.(best_stars.L_CMD_SAT[filt]),  label="cmd")

	stephist!(log10.(best_stars.L_PM_SAT[filt]),  label="pm")

	axislegend(position=:lt)
	fig
end

# ╔═╡ ca06c336-3950-4831-9f3f-7fe41406bc08
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log likilihood",
		ylabel = "counts",
		title = "members"
	)
	bins = LinRange(-10, 2, 30)
	
	stephist!(log10.(members.L_S_SAT), 
		bins=bins, label="spatial")
	
	stephist!(log10.(members.L_CMD_SAT), bins=bins, label="cmd")
	
	stephist!(log10.(members.L_PM_SAT), bins=bins,  label="pm")

	axislegend(position=:lt)
	fig
end

# ╔═╡ 63a3e02b-0ea0-472d-8c2e-962bd8e770c0
let 
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log likelihood sat / background",
		ylabel = "counts",
		title = "members"
	)
	
	stephist!(log10.(members.L_S_SAT ./ members.L_S_BKD), 
		label="spatial")
	
	stephist!(log10.(members.L_CMD_SAT ./ members.L_CMD_BKD), label="cmd")
	
	stephist!(log10.(members.L_PM_SAT ./ members.L_PM_BKD), label="pm")

	axislegend(position=:lt)
	fig
end

# ╔═╡ 3a52eba9-e07d-4445-a322-1c4b8bb33024
md"""
### Correlations
"""

# ╔═╡ f54c50ae-0e03-4917-899f-81a30f116b3f
function plot_corr(best_stars; markersize=3, title="")
	fig = Figure(
		size = (800, 400)
	)

	kwargs = (;
		markersize=markersize,
	)
	
	local p

	axs = []
	
	for i in 1:3
		col = [:LL_S, :LL_CMD, :LL_PM][i]
		label = ["space", "cmd", "pm"][i]

		ax = Axis(fig[1,i],
			xlabel = L"\log\mathcal{L}_\textrm{%$label}",
			ylabel = L"\log\mathcal{L}",
			aspect=DataAspect(),
			title = i == 2 ? title : ""
		)
		
		
		p = scatter!(best_stars[!, col], best_stars.LL;
			color= best_stars.PSAT ,
			colorrange=(0, 1),
			kwargs...
		)


		# correlation coefficient
		R = cor(best_stars[!, col], best_stars.LL)
		t = "ρ = $(round(R, digits=2))"
		text!(0.1, 0.9, text=t, space=:relative)

		# formatting
		if i > 1
			hideydecorations!(ax, grid=false, ticks=false, minorticks=false)
		end
		LilGuys.hide_grid!(ax)
		push!(axs, ax)

	end

	Colorbar(fig[1,4], p, label = "PSAT")

	linkaxes!(axs...)
	fig
end

# ╔═╡ b4f816f0-0abe-43c2-bf5b-48a2eb2b99c8
plot_corr(members, title="members")

# ╔═╡ dc280eae-4f10-4167-9545-c70d67af2d8f
let
	fig = plot_corr(best_stars[best_stars.LL .> 0, :], markersize=5, title = "LLR > 0")

	@savefig "LL_corr_trunc"

	fig
end

# ╔═╡ fb3a27be-62c7-41ab-ba5d-761efb74c044
let
	fig = plot_corr(best_stars[best_stars.LL .> -10, :], title="LLR > -10")

	fig
end

# ╔═╡ b4363369-4590-4446-9f4a-9f71ff31a870
plot_corr(best_stars[best_stars.LL .> -20, :], title = "LLR > -20")

# ╔═╡ 00d4679d-6c55-4929-bbd3-c796479ba70d
plot_corr(all_stars[all_stars.LL .> -80, :], title = "LLR > -80")

# ╔═╡ 6be3bd9d-c257-444f-87ab-2b3956fc83aa
md"""
## Fractional contributions
"""

# ╔═╡ 4ad0330e-f1d7-4d62-b6f9-f49be1aa0d80
md"""
Still not sure how useful these plots are, but they should represent an approximate metric for how substantially the log-likelihood is influenced by any component.
"""

# ╔═╡ 9a2ca944-dd5c-43f5-b70d-ab3d67b42743
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "LL",
		ylabel = "f spatial",
		limits=(-30, 10, -1.5, 3)
	)

	filt = best_stars.LL .|> isfinite
	scatter!(best_stars.LL[filt], (best_stars.f_LL_S[filt]), markersize=3, label = "best")
	scatter!(members.LL, (members.f_LL_S), label = "members")

	axislegend()
	fig
end

# ╔═╡ 79ac654d-7de9-4a2c-bc5d-2f8c8a7e057f
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "LL",
		ylabel = "f spatial normed",
		limits=(-30, 10, -1, 1)
	)

	filt = best_stars.LL .|> isfinite

	scatter!(best_stars.LL[filt], (best_stars.f_LL_S_norm[filt]), markersize=3)
	scatter!(members.LL, (members.f_LL_S_norm))

	fig
end

# ╔═╡ 479c9cad-8cae-4d9b-87d3-79a816b3af82
md"""
The three plots below verify that the normalized version of likelihood fractions makes sense as a diagnostic.
"""

# ╔═╡ a7dbb8bf-8c82-4ace-9f17-e9aa159a84ff
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "f spatial",
		ylabel = "f spatial normed",
		limits=(0, 1.5, 0, 1.5),
		aspect = DataAspect(),
	)
	
	scatter!(abs.(best_stars.f_LL_S), abs.(best_stars.f_LL_S_norm), markersize=3)
	scatter!(abs.(members.f_LL_S), abs.(members.f_LL_S_norm))

	fig
end

# ╔═╡ ec4a8ca8-a4b9-41d7-a369-8049b2a98795
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "f CMD",
		ylabel = "f CMD normed",
		limits=(0, 3, 0, 1.5),
		aspect = DataAspect(),
	)
	
	scatter!(abs.(best_stars.f_LL_CMD), abs.(best_stars.f_LL_CMD_norm), markersize=3)
	scatter!(abs.(members.f_LL_CMD), abs.(members.f_LL_CMD_norm))

	fig
end

# ╔═╡ 32fd5ed6-fb41-456c-bd05-d28ad013c2c7
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "f PM",
		ylabel = "f PM normed",
		limits=(0, 3, 0, 1.5),
		aspect = DataAspect(),
	)
	
	scatter!(abs.(best_stars.f_LL_PM), abs.(best_stars.f_LL_PM_norm), markersize=3)
	scatter!(abs.(members.f_LL_PM), abs.(members.f_LL_PM_norm))

	fig
end

# ╔═╡ 301f42ce-437f-4d45-a186-42df47303f25
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "f PM",
		ylabel = "f PM unnormed",
		limits=(-10, 10, -10, 10),
		aspect = DataAspect(),

	)

	filt = best_stars.LL .|> isfinite

	scatter!((best_stars.LL_PM .- best_stars.LL_norm)[filt], (best_stars.LL_PM .- best_stars.LL)[filt])

	fig
end

# ╔═╡ 74dddfdd-0195-446c-b499-9963d3b73059
extrema(best_stars.LL_PM .- best_stars.LL_norm)

# ╔═╡ 625795e7-aca7-4bb9-8db6-cfcfb87bddaf
import DensityEstimators as DE

# ╔═╡ 24e23b26-e5ca-4e5c-bea0-ecccb9a385ae
import StatsBase: quantile

# ╔═╡ 62a6fb8a-c49e-477f-a2a9-c21e9e006013
function medianscatter!(xs, ys; kwargs...)

	bins = DE.bins_equal_number(xs)
	bin_m = midpoints(bins)

	ym = DE.binned_statistic(xs, ys, bins, statistic=DE.median)
	yl = DE.binned_statistic(xs, ys, bins, statistic=x->quantile(x, 0.16))
	yh = DE.binned_statistic(xs, ys, bins, statistic=x->quantile(x, 1-0.16))

	ye = collect(zip(ym .- yl, yh .- ym))
	errorscatter!(bin_m, ym, yerror=ye; kwargs...)
	lines!(bin_m, ym; kwargs...)

end

# ╔═╡ 8bebb46c-f62f-4e4c-9334-a83a7d930c97
Cycled(2)

# ╔═╡ 84b10da6-92f2-4ccc-b548-5106f379e09f
function binhists!(xs, ys; limits=nothing, kwargs...)

	bins = DE.bins_equal_number(xs)

	idxs = DE.bin_indices(xs, bins)

	yh_max = 0.9 * minimum(diff(bins))

	h_ys= []
	h_xs = []
	dys = []
	for i in unique(idxs)
		filt = idxs .== i
		y = ys[filt]
		xm = midpoints(bins)[i]
		h = DE.histogram(y, limits=limits)
		push!(h_ys, h.values)
		push!(h_xs, midpoints(h.bins))
		push!(dys, xm)
	end

	h_scale = yh_max ./ maximum.(h_ys) 
	h_ys .*= h_scale

	for i in eachindex(h_ys)
		lines!(h_xs[i], h_ys[i] .+ dys[i]; kwargs...)
	end

end

# ╔═╡ 7ff47e82-efd7-45a3-abcc-9d8b0e764e8e
let
	fig = Figure()
	ax = Axis(fig[1,1],
		ylabel = "log r ell / arcmin",
		xlabel = "fsat space"
	)

	xs = members_nospace.r_ell .|> log10
	ys = members_nospace.f_LL_S_norm

	binhists!(xs, ys, label = "spatial", color=COLORS[1])

	ys = members_nospace.f_LL_PM_norm
	binhists!(xs, ys, label = "PM", color=COLORS[2])

	ys = members_nospace.f_LL_CMD_norm
	binhists!(xs, ys, label = "CMD", color=COLORS[3])
	fig
end

# ╔═╡ 184f312a-6e6c-43cc-9864-75b4846224c2
let
	fig = Figure()
	ax = Axis(fig[1,1],
		ylabel = "log r ell / arcmin",
		xlabel = "LL"
	)

	df = members_nospace
	
	xs = df.r_ell .|> log10
	ys = df.LL_S

	binhists!(xs, ys, label = "spatial", color=COLORS[1])

	ys = df.LL_PM
	binhists!(xs, ys, label = "PM", color=COLORS[2])

	ys = df.LL_CMD
	binhists!(xs, ys, label = "CMD", color=COLORS[3])

	ys = df.LL
	binhists!(xs, ys, label = "total", color=:black)

	vlines!(minimum(members.LL))
	fig
end

# ╔═╡ f229556e-4b94-41ad-a523-ad50509a8ed8
best_stars.L_S_SAT

# ╔═╡ b1b1b73a-b0ca-4776-84db-b9be3c94a949
let
	fig = Figure(size=(600, 600))
	ax = Axis(fig[1,1],
		ylabel = plot_labels[:r_ell],
		xlabel = plot_labels[:LL],
		title = "All stars"
	)

	limits = (-10, 5)

	df = best_stars[best_stars.LL .> -Inf, :]
	
	xs = df.r_ell .|> log10
	ys = df.LL_S

	binhists!(xs, ys, limits=limits, label = "spatial", color=COLORS[1])

	ys = df.LL_PM
	binhists!(xs, ys, limits=limits, label = "PM", color=COLORS[2])

	ys = df.LL_CMD
	binhists!(xs, ys,limits=limits,  label = "CMD", color=COLORS[3])

	ys = df.LL
	binhists!(xs, ys,limits=limits,  label = "total", color=:black)

	vlines!(minimum(best_stars.LL[best_stars.PSAT .> 0.2]))

	Legend(fig[1,2], ax, unique=true)
	fig
end

# ╔═╡ 69d072df-f522-454a-b4d5-bc22122af524
log10.(maximum(best_stars.L_S_SAT_Outer ./ best_stars.L_S_SAT_Inner))

# ╔═╡ 4768871c-1142-491c-9cd4-bb4c908a61b2
let
	fig = Figure(size=(600, 600))
	ax = Axis(fig[1,1],
		ylabel = plot_labels[:r_ell],
		xlabel = plot_labels[:LL],
		title = "All stars"
	)

	limits = (-10, 5)

	df = best_stars[best_stars.LL .> -Inf, :]
	
	xs = df.r_ell .|> log10
	ys = log10.(df.L_S_SAT_Inner ./ df.L_S_BKD)
	binhists!(xs, ys, limits=limits, label = "spatial inner", color=COLORS[1])

	ys = log10.(df.L_S_SAT_Outer ./ df.L_S_BKD)
	binhists!(xs, ys, limits=limits, label = "spatial outer", color=COLORS[2])
	ys = df.LL_S
	binhists!(xs, ys, limits=limits, label = "spatial", color=COLORS[3])

	vlines!(minimum(best_stars.LL[best_stars.PSAT .> 0.2]))

	Legend(fig[1,2], ax, unique=true)
	fig
end

# ╔═╡ 4793a9e8-7a58-42ab-997a-896c20edbbcc
let
	fig = Figure(size=(600, 600))
	ax = Axis(fig[1,1],
		ylabel = plot_labels[:r_ell],
		xlabel = plot_labels[:LL],
		title = "members"
	)

	limits = (-5, 5)

	df = members
	
	xs = df.r_ell .|> log10
	ys = df.LL_S

	binhists!(xs, ys, limits=limits, label = "spatial", color=COLORS[1])

	ys = df.LL_PM
	binhists!(xs, ys, limits=limits, label = "PM", color=COLORS[2])

	ys = df.LL_CMD
	binhists!(xs, ys,limits=limits,  label = "CMD", color=COLORS[3])

	ys = df.LL
	binhists!(xs, ys,limits=limits,  label = "total", color=:black)

	vlines!(minimum(best_stars.LL[best_stars.PSAT .> 0.2]))

	Legend(fig[1,2], ax, unique=true)
	fig
end

# ╔═╡ ed252f72-dff0-40c5-ac64-74ad6bcffd56
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log r ell / rh",
		ylabel = "fsat space"
	)

	xs = members_nospace.r_ell .|> log10
	xs .-= log10(observed_properties["r_h"])
	ys = members_nospace.f_LL_S_norm

	medianscatter!(xs, ys, label = "spatial")

	ys = members_nospace.f_LL_PM_norm
	medianscatter!(xs, ys, label = "PM")

	ys = members_nospace.f_LL_CMD_norm
	medianscatter!(xs, ys, label = "CMD")

	Legend(fig[1,2], ax, unique=true, merge=true)
	fig
end

# ╔═╡ 3b8ab87d-2226-4d30-9c45-ec0d8181bf45
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log r ell / arcmin",
		ylabel = "LLR"
	)

	df = best_stars[isfinite.(best_stars.LL), :]
	xs = df.r_ell .|> log10
	
	ys = df.LL_S .- df.LL
	medianscatter!(xs, ys, label = "spatial")

	ys = df.LL_PM .- df.LL
	medianscatter!(xs, ys, label = "spatial")

	
	ys = df.LL_CMD .- df.LL
	medianscatter!(xs, ys, label = "spatial")

	axislegend()
	fig
end

# ╔═╡ dbc3f4af-9e29-4437-a203-63d2f3d89ccc
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log r ell / arcmin",
		ylabel = "LLR"
	)

	df = members
	xs = df.r_ell .|> log10
	
	ys = df.LL_S 
	medianscatter!(xs, ys, label = "spatial")

	ys = df.LL_PM
	medianscatter!(xs, ys, label = "PM")

	
	ys = df.LL_CMD
	medianscatter!(xs, ys, label = "CMD")

	ys = df.LL
	medianscatter!(xs, ys, label = "total")

	hlines!(minimum(members.LL), label = "cutoff")
	Legend(fig[1,2], ax)
	fig
end

# ╔═╡ 59007a2d-8a11-4b1f-aedf-606f8de9fd12
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log r ell / arcmin",
		ylabel = "LLR",
		limits=(nothing, nothing, -10, 10),
	)
	jitter = 0.05
	df = best_stars[best_stars.LL .> -Inf, :]
	xs = df.r_ell .|> log10

	ys = df.LL
	medianscatter!(xs, ys, label = "total")

	
	ys = df.LL_S 
	medianscatter!(xs, ys, label = "spatial")

	ys = df.LL_PM
	medianscatter!(xs .+ jitter, ys, label = "PM")

	
	ys = df.LL_CMD
	medianscatter!(xs .- jitter, ys, label = "CMD")


	hlines!(minimum(members.LL), label = "cutoff")
	Legend(fig[1,2], ax)
	fig
end

# ╔═╡ e5cb3522-a756-4ca0-b843-b38d66ed6a7a


# ╔═╡ 6417beb2-eeaf-4c50-9fd5-ffeaf46bf64a
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log r ell / arcmin",
		ylabel = "LLR",
		limits=(nothing, nothing, -10, 10),
	)

	df = best_stars[best_stars.LL .> -Inf, :]
	xs = df.r_ell .|> log10
	
	ys = df.LL_S 
	medianscatter!(xs, ys, label = "spatial")

	ys = df.LL_PM
	medianscatter!(xs, ys, label = "PM")

	
	ys = df.LL_CMD
	medianscatter!(xs, ys, label = "CMD")

	ys = df.LL
	medianscatter!(xs, ys, label = "total")

	hlines!(minimum(members_nospace.LL), label = "cutoff")
	Legend(fig[1,2], ax)
	fig
end

# ╔═╡ 3c05db74-0cc1-4893-a468-b495c0121118
md"""
## Fractional contributions on the sky
"""

# ╔═╡ e3f951ed-3725-4ece-bf59-cb6998fdf170
function plot_f_l_tangent(best_stars; title="", markersize=5)
	fig = Figure(
		size=(800, 500)
	)

	filt = best_stars.PSAT .> 0.2
	df = best_stars[filt, :]


	for i in 1:3
		col = [:f_LL_S_norm, :f_LL_PM_norm, :f_LL_CMD_norm][i]
		label = ["space", "proper motion", "CMD"][i]
		
		ax = Axis(fig[1,i],
			title = label,
			xlabel = "f likelihood",
			ylabel = "counts",
			limits=(-1, 1, 0, nothing),
		)
		LilGuys.hide_grid!(ax)

		h = histogram(df[:, col])
		
		x = [h.bins[1]; midpoints(h.bins); h.bins[end]]
		y = [0; h.values; 0]
		lines!(x, y, color=x,
			colorrange=(-1, 1),
			colormap=diverging_cmap,
		)

		
		ax = Axis(fig[2,i],
			xlabel = "ξ / degree",
			ylabel = "η / degree",
			aspect=DataAspect(),
			#limits = (-2, 2, -2, 2)
			title = i == 2 ? title : ""
			)
		
		p = scatter!(df.xi, df.eta,
			markersize=markersize,
			color = df[:, col],
			colorrange=(-1, 1),
			colormap=diverging_cmap,
		)
		
		LilGuys.hide_grid!(ax)
	end



	rowsize!(fig.layout, 1, Relative(0.2))
	rowgap!(fig.layout, 1, 0)

	return fig
end

# ╔═╡ e4e6239f-9fdd-433b-8a39-ffaeaf96b59d
let
	fig = plot_f_l_tangent(members, title="members")

	@savefig "xi_eta_vs_f"
	fig
end

# ╔═╡ 363dfc70-ec52-44ed-ade4-f88322247871
let
	fig = plot_f_l_tangent(best_stars[best_stars.LL .> 0, :], title="LL > 0")

	@savefig "xi_eta_vs_f"
	fig
end

# ╔═╡ c68ff47f-97de-4bbb-8b0a-42d837da21f7
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = L"\log\,r_\textrm{ell}",
		ylabel = "f spatial",
	)
	
	p = scatter!(log10.(best_stars.r_ell), (best_stars.f_LL_S_norm), 
		color=(best_stars.LL),
		colorrange=(-10, 1),
		markersize=3,
	)
	LilGuys.hide_grid!(ax)
	
	Colorbar(fig[1,2], p, label="log Likelihood sat / background", )

	@savefig "f_space_vs_r_ell"
	fig
end

# ╔═╡ 8d497b75-a47d-42d1-9e54-e15bfb1db805
let
	fig = Figure()

	ax = Axis(fig[1,1],
		xlabel = "r_ell",
		ylabel = "f spatial",
	)
	
	scatter!((best_stars.r_ell), (best_stars.f_LL_S_norm), markersize=3)

	
	scatter!((members.r_ell), (members.f_LL_S_norm))

	fig
end

# ╔═╡ d2668f44-9cb7-4ab8-af49-1aef2baecdda
md"""
# Sanity checks
"""

# ╔═╡ 12bf64e9-ae49-4f93-8572-660dbc86ffc1
scatter(all_stars.r_ell, all_stars.r_ell_old * observed_properties["r_h"] * sqrt(1 - observed_properties["ellipticity"]), 
	axis = (;
	xlabel = "r ell",
	ylabel = "r ell jax",
	aspect = DataAspect(),
	)
)

# ╔═╡ Cell order:
# ╠═47b8b3b0-0228-4f50-9da4-37d388ef9e9f
# ╠═303ff946-ecdc-42c3-8289-c13e7a8cfe7c
# ╠═eca9c1c5-e984-42d4-8854-b227fdec0a8a
# ╠═bff50014-bfa9-11ee-33f0-0f67e543c2d4
# ╠═8acddd07-0eff-47a6-a364-f7f58680bd9c
# ╠═8b3955e5-ab08-4d5a-abe5-8f40120ff8db
# ╠═da9ca7a7-18b8-49cb-a041-ab1c667920ff
# ╠═489f6d21-9e9a-4b1e-b05c-c63a44ba1951
# ╠═0004f638-c57a-4dab-8b97-c77840cafbbf
# ╠═ae29bed0-6700-47f1-8952-35e867ce126b
# ╠═db1264b7-02c3-4b55-ae2b-9ce78aa1304a
# ╠═81c677bf-8372-4fd3-9d60-af8377485f1c
# ╠═1fbbd6cd-20d4-4025-829f-a2cc969b1cd7
# ╠═cf7aeb0a-d453-462a-b43d-a832567440fd
# ╠═fe7050d0-fede-44c0-8a8e-5d3dcaf0d1fb
# ╠═ff41b9a8-e25c-4fc8-b00f-0d74ea6a1d31
# ╠═a955fb87-3674-48dc-8dbd-2c61926d195d
# ╟─2eb4aa78-0fea-460b-a18e-06a129c41504
# ╠═9ecf79a8-2ed3-40c6-b555-a102250ecbd4
# ╠═4c8b9c30-f873-4756-a729-aa023b3a804e
# ╠═26cf1867-02be-4d36-8c35-6c58a1feca27
# ╠═ec227641-86e6-46b7-8019-9b02072ed9f7
# ╟─7ed70894-48a8-46af-8c78-8a75cb39d957
# ╠═223abc41-5158-49c2-96bf-df55b7be1114
# ╟─3abea230-2b95-40d3-9851-a91236f75c4a
# ╠═32e8e100-1fad-425f-b05c-2310e7ba559c
# ╠═b1be5f97-0e61-4dda-87fd-7bb0f43147b6
# ╠═731ea468-5003-44e9-95b8-7fa7ef4b475b
# ╠═88fbdd09-30be-4fc3-95ae-acce6e0018e1
# ╠═60d0e593-88fd-4b4c-9009-cc24a597c6d5
# ╠═082a06dd-eeb5-4761-a233-1ee89e8cb819
# ╟─a9d94121-ea6e-416a-bae8-aa93c16bde72
# ╠═f3d3f7b2-30c2-4f89-bef5-ac1a8720963f
# ╠═7cc17173-160c-45f1-a685-c79271679414
# ╠═965217b8-b2a5-485b-83de-cac887065b19
# ╠═2d4da56b-5d0a-49d3-83ae-a90f85192101
# ╟─b0aec8b4-410f-4678-ad8c-cb055b743bd3
# ╟─3c2ee4bb-d594-47e7-a7f3-5b7c9e38b7d4
# ╟─b191b7b8-28ed-4aee-b376-b2bdf22abc29
# ╟─b06445ed-6290-42e8-bdb8-20f0372518c2
# ╟─2a8e13bd-b8ca-4cf3-b61b-096af218265b
# ╟─25dd2e38-1c78-4ef6-95ef-d38a0e80f0e4
# ╟─1b63c1ee-0a46-49c6-880e-c1742055605e
# ╟─643773f0-2b03-46f7-8df4-579b3c708909
# ╟─cad88ab4-1ada-4cc3-8ef0-d4674500d57e
# ╟─18252c3c-d4d6-46f2-8277-33a700ab4074
# ╟─d47a8d02-6c63-41fd-95b5-85550d1f5a85
# ╟─b7e7e245-3732-453c-a2ec-3fe84c65c012
# ╟─87f29963-2c7f-4d8f-8631-f70197d315c0
# ╟─72aa995b-8dd1-4d08-88c2-5fb527ff7f3e
# ╟─6cdcd105-fa1d-48f5-9861-b722aafa28a3
# ╟─cf8e51a9-194c-4aa2-9766-39256a35ef02
# ╟─e265ec56-d62d-4295-975f-98ba93cf5e16
# ╟─e2686165-6499-4558-9da9-91a78624b992
# ╟─1f2aa1ee-b395-4040-84cb-973d57b534d5
# ╟─ca06c336-3950-4831-9f3f-7fe41406bc08
# ╟─63a3e02b-0ea0-472d-8c2e-962bd8e770c0
# ╟─3a52eba9-e07d-4445-a322-1c4b8bb33024
# ╠═f54c50ae-0e03-4917-899f-81a30f116b3f
# ╠═b4f816f0-0abe-43c2-bf5b-48a2eb2b99c8
# ╟─dc280eae-4f10-4167-9545-c70d67af2d8f
# ╟─fb3a27be-62c7-41ab-ba5d-761efb74c044
# ╟─b4363369-4590-4446-9f4a-9f71ff31a870
# ╟─00d4679d-6c55-4929-bbd3-c796479ba70d
# ╟─6be3bd9d-c257-444f-87ab-2b3956fc83aa
# ╟─4ad0330e-f1d7-4d62-b6f9-f49be1aa0d80
# ╠═9a2ca944-dd5c-43f5-b70d-ab3d67b42743
# ╟─79ac654d-7de9-4a2c-bc5d-2f8c8a7e057f
# ╟─479c9cad-8cae-4d9b-87d3-79a816b3af82
# ╟─a7dbb8bf-8c82-4ace-9f17-e9aa159a84ff
# ╟─ec4a8ca8-a4b9-41d7-a369-8049b2a98795
# ╠═32fd5ed6-fb41-456c-bd05-d28ad013c2c7
# ╟─301f42ce-437f-4d45-a186-42df47303f25
# ╠═74dddfdd-0195-446c-b499-9963d3b73059
# ╠═625795e7-aca7-4bb9-8db6-cfcfb87bddaf
# ╠═24e23b26-e5ca-4e5c-bea0-ecccb9a385ae
# ╠═62a6fb8a-c49e-477f-a2a9-c21e9e006013
# ╠═8bebb46c-f62f-4e4c-9334-a83a7d930c97
# ╠═84b10da6-92f2-4ccc-b548-5106f379e09f
# ╠═7ff47e82-efd7-45a3-abcc-9d8b0e764e8e
# ╠═184f312a-6e6c-43cc-9864-75b4846224c2
# ╠═f229556e-4b94-41ad-a523-ad50509a8ed8
# ╠═b1b1b73a-b0ca-4776-84db-b9be3c94a949
# ╠═69d072df-f522-454a-b4d5-bc22122af524
# ╠═4768871c-1142-491c-9cd4-bb4c908a61b2
# ╠═4793a9e8-7a58-42ab-997a-896c20edbbcc
# ╠═ed252f72-dff0-40c5-ac64-74ad6bcffd56
# ╠═3b8ab87d-2226-4d30-9c45-ec0d8181bf45
# ╠═dbc3f4af-9e29-4437-a203-63d2f3d89ccc
# ╠═59007a2d-8a11-4b1f-aedf-606f8de9fd12
# ╠═e5cb3522-a756-4ca0-b843-b38d66ed6a7a
# ╠═6417beb2-eeaf-4c50-9fd5-ffeaf46bf64a
# ╟─3c05db74-0cc1-4893-a468-b495c0121118
# ╠═e3f951ed-3725-4ece-bf59-cb6998fdf170
# ╠═e4e6239f-9fdd-433b-8a39-ffaeaf96b59d
# ╠═363dfc70-ec52-44ed-ade4-f88322247871
# ╟─c68ff47f-97de-4bbb-8b0a-42d837da21f7
# ╟─8d497b75-a47d-42d1-9e54-e15bfb1db805
# ╟─d2668f44-9cb7-4ab8-af49-1aef2baecdda
# ╠═12bf64e9-ae49-4f93-8572-660dbc86ffc1
