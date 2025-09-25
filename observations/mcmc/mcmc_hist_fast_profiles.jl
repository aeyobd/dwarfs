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

	using LilGuys
	using PyFITS

	import TOML
	using OrderedCollections
	using DataFrames, CSV
end

# ╔═╡ e997f77b-02d5-406a-8f78-b9e115baaf88
using PlutoUI

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	galaxy = ARGS[1]
else
	@bind galaxy confirm(TextField(default="leo2"))
end

# ╔═╡ 1f497448-9ea2-460f-b403-618b78b47565
module MCMCUtils
	include("mcmc_utils.jl")
end

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

# ╔═╡ 11d9d9f5-0fb7-4b06-bd4a-36bdb1f00188
outdir = joinpath("..", galaxy, "mcmc")

# ╔═╡ b94c3346-bd31-409e-ad6f-5e7afb891ad1
FIGDIR = joinpath(outdir, "figures"); FIGSUFFIX=".mcmc_hist_fast"

# ╔═╡ 0aa44389-f320-4274-abdd-0d7f88006a4d
log_r_label = L"$\log\,R_\textrm{ell}$ / arcmin"

# ╔═╡ 36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
log_Sigma_label = L"$\log\,\Sigma$"

# ╔═╡ 04053b71-bd55-40d7-885d-6df67035e3d6
md"""
# data loading
"""

# ╔═╡ b7b4361a-1ee8-4b1d-a62e-9f7bd0363160
obs_props = MCMCUtils.get_obs_props(galaxy)

# ╔═╡ aae56f16-2433-4847-83db-9184a74f28f8
Threads.nthreads()

# ╔═╡ 4b77612f-e124-4d77-974c-d40c2f5a37ff
struct_params = MCMCUtils.StructuralParams(; LilGuys.dict_to_tuple(TOML.parsefile(joinpath("..", galaxy, "mcmc", "default_bins.toml")))...)

# ╔═╡ 9726210e-278f-4478-9c4a-8be8975a57d6
bins = struct_params.bins

# ╔═╡ c58ed939-07c1-4d34-90a4-3c5d9cc9910b
stars = MCMCUtils.get_fits(galaxy, obs_props)

# ╔═╡ 6f016a8e-38ae-4f05-a7ee-c292ac0e5741
df_chains = CSV.read(joinpath(outdir, "samples.mcmc_hist_fast.csv"), DataFrame)

# ╔═╡ 3e99a9b9-ba06-47ee-bce7-5720a4eb1944
df_summary = CSV.read(joinpath(outdir, "summary.mcmc_hist_fast.csv"), DataFrame)

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Analysis
"""

# ╔═╡ 3067ae10-854f-48ce-8285-0d3f3fa7f25c
md"""
We need the counts to turn the log densities in the MCMC hist into true densities
"""

# ╔═╡ d254b67b-8017-4bbf-ba04-50b5b76f48f0
_prof_counts = LilGuys.SurfaceDensityProfile(stars.R_ell, bins=log10.(bins),  errors=:weighted)

# ╔═╡ 4b644aaf-76ab-4e30-8a78-7e127c2ba7b5
md"""
Consistency checks:
"""

# ╔═╡ 85f59eb3-4f90-498a-923e-36cd679af1c0
@assert df_summary.log_R ≈ midpoints(log10.(bins))

# ╔═╡ 575a31e1-df22-42e9-9a9e-1dc4b678a5b6
@assert isapprox(df_summary.N_stars, _prof_counts.counts, atol=2)

# ╔═╡ f0afb784-c1f4-40f6-84e6-a5c8f12ac67b
log_Sigmas_fast = let
	Nbins = length(bins) - 1

	Nc = size(df_chains, 1)
	Ns = size(stars, 1)
	
	log_Sigmas_fast = Matrix{Float64}(undef, Nbins, Nc)
	areas = diff(π * bins .^2)
	counts = df_summary.N_stars

	for i in 1:Nc
		if i % 200 == 0
			@info "calculated $i / $Nc"
		end
		
		params = [df_chains[i, "params[$j]"] for j in 1:Nbins]
		f_sat = exp10.(params) ./(1 .+ exp10.(params))
		log_Sigmas_fast[:, i] .= log10.( f_sat .* counts ./ areas)
	end

	log_Sigmas_fast

end

# ╔═╡ fe3d955f-c96b-4202-9218-f749a5e9d41f
fsats, psats = let
	Nbins = length(bins) - 1

	Nc = size(df_chains, 1)
	Ns = size(stars, 1)

	idxs = 1:Nc
	fsats = Matrix{Float64}(undef, Ns, length(idxs))
	psats = Matrix{Float64}(undef, Ns, length(idxs))

	Threads.@threads for i in eachindex(idxs)
		idx = idxs[i]
		if i % 10 == 0
			@info "sample $i"
		end
		
		params = [df_chains[idx, "params[$j]"] for j in 1:Nbins]
		radii = stars.R_ell

		Σ = MCMCUtils.Σ_hist(radii, bins, params)
		f = 10 .^ Σ ./ (1 .+ 10 .^ Σ)
		fsats[:, i] .= f

		Ls = stars.L_CMD_SAT .* stars.L_PM_SAT
		Lb = stars.L_CMD_BKD .* stars.L_PM_BKD
		psat =  @. f*Ls / (f*Ls + (1-f) * Lb)
		psats[:, i] .= psat

	end

	fsats, psats
end


# ╔═╡ 5feaa9b0-ca44-48d6-80d1-16303a34005f
prof_jax_cut = LilGuys.SurfaceDensityProfile(stars.R_ell[stars.PSAT .> 0.2], bins=log10.(bins), errors=:weighted, normalization=:none)

# ╔═╡ ab124afe-c94a-44d9-8869-bd4d4cee3fbd
prof_jax_weighted = LilGuys.SurfaceDensityProfile(stars.R_ell[stars.PSAT .> 0], weights=stars.PSAT[stars.PSAT .> 0] |> disallowmissing, bins=log10.(bins), errors=:weighted,  normalization=:none)

# ╔═╡ 6685cdec-7c99-45ca-98c8-f6439dfd7290
md"""
## MCMC Profiles
"""

# ╔═╡ 6470edfd-bd06-43b7-aa5e-a80645af78fb
psats_pm = MCMCUtils.to_measurements(psats)

# ╔═╡ c6a8bf2b-8a65-48ff-a1ca-7790f7806a94
log_Sigma_mean = MCMCUtils.to_measurements(log_Sigmas_fast)

# ╔═╡ f80631b5-aa7e-454a-937f-c2e66b042105
profiles = [
	LilGuys.SurfaceDensityProfile(
		R_units = "arcmin",
		log_Sigma = log_Sigma,
		log_R_bins = log10.(bins),
		log_R = midpoints(log10.(bins)),
		Gamma = LilGuys.gradient(log_Sigma, midpoints(log10.(bins)))
	) 
	for log_Sigma in eachcol(log_Sigmas_fast)
]

# ╔═╡ 93fe93ec-b98c-4037-9504-bf584e31b420
Gamma_mean = MCMCUtils.to_measurements(
	hcat((middle.(prof.Gamma) for prof in profiles)...))

# ╔═╡ 4106ea1f-948f-42fb-9a19-64580a94ce9a
prof_mean = LilGuys.SurfaceDensityProfile(
		R_units = "arcmin",
		log_Sigma=log_Sigma_mean,
		log_R_bins=log10.(bins), # use normal bins for these
		log_R=midpoints(log10.(bins)),
		Gamma = Gamma_mean
	)

# ╔═╡ b8ba7f13-6a1a-43d5-90db-930196dfbbe4
md"""
# Plots
"""

# ╔═╡ da737b3e-f28a-4917-91d6-eea97458ddb0
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "PSAT J+24", ylabel = "PSAT MCMC")

	skip = 1000

	for i in 1:skip:size(psats, 2)
		y = psats[:, i]
		filt = 1e-6 .< y .< 1-1e-6
		scatter!(stars.PSAT[filt], y[filt], color=:black, alpha=0.01, markersize=2, rasterize=true)
	end
	@savefig "psat_j+24_vs_mcmc"
	
	fig
end

# ╔═╡ ad85b789-8f43-489b-a41c-a966e08d78ad
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log r_ell", ylabel = "psat me - jax")

	skip = 100

	scatter!(log10.(stars.R_ell), middle.(psats_pm) .- stars.PSAT, color=:black, alpha=0.03, markersize=1, rasterize=true)

	@savefig "j+24_vs_mcmc_diff"

	fig
end

# ╔═╡ f2c43c30-b069-4d17-8d61-be801d85b245
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel =log_Sigma_label,
		limits = (nothing, nothing, -8, 3)
	)
	skip = 100
	jitter = 0.005
	
	for h in profiles[1:skip:end]
		filt = isfinite.(h.log_Sigma)
		scatter!(h.log_R[filt] .+ jitter * randn(sum(filt)), h.log_Sigma[filt], color=:black, alpha=0.03, markersize=1)
	end

	@savefig "scatter_profiles"

	fig
end

# ╔═╡ c507975a-8041-48f0-b46e-ac1ebccf3a94
let
	jitter = 0.005
	fig = Figure(size=(5*72, 3*72))
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = L"\Gamma",
		limits = (nothing, nothing, -10, 2)
	)
	skip = 100
	
	for h in profiles[1:skip:end]

		Gamma = LilGuys.gradient(middle.(h.log_Sigma), h.log_R)
		filt = isfinite.(Gamma)

		scatter!(h.log_R[filt].+ jitter*randn(sum(filt)), Gamma[filt], color=:black, alpha=0.1, markersize=1)
	end

	ax2 = Axis(fig[1,2],
		xlabel = L"$R_\textrm{ell}$ / arcmin",
		ylabel = L"\Gamma",
		limits = (nothing, nothing, -10, 2)
	)
	skip = 10
	
	for h in profiles[1:skip:end]

		Gamma = LilGuys.gradient(middle.(h.log_Sigma), h.log_R)
		filt = isfinite.(Gamma)

		scatter!(10 .^ (h.log_R[filt].+ jitter*randn(sum(filt))), Gamma[filt], color=:black, alpha=0.03, markersize=1)
	end

	linkyaxes!(ax, ax2)
	hideydecorations!(ax2, ticks=false, minorticks=false)
	colgap!(fig.layout, 0)
	
	@savefig "scatter_slope"


	fig
end

# ╔═╡ bab50892-fab8-44fe-931b-66ff466b6974
let
	jitter = 0.005
	fig = Figure(size=(5*72, 3*72))
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = L"\Gamma",
		limits = (nothing, nothing, -10, 2)
	)


	profs = OrderedDict(
		"J+24 probability cut" => prof_jax_cut,
		"J+24 weighted" => prof_jax_weighted,
		"MCMC hist" => prof_mean,

	)

	
	for (i, (label, prof)) in enumerate(profs)
		errorscatter!(prof.log_R .+ 0.01*i, prof.Gamma, yerror=error_interval.(prof.Gamma),
					  label=label, markersize=3, errorcolor=(COLORS[i], 0.5))
	end

	axislegend(position=:lb)
	
	@savefig "Gamma_j+24_vs_mcmc"

	fig

end

# ╔═╡ 7c02d24f-7d59-410b-99bc-93513afa4c5d
let
	fig = Figure(size=(4*72, 3*72))
	y_max = maximum(middle.(prof_jax_cut.log_Sigma))
	
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
		limits = (nothing, nothing, y_max + 0.5 - 8, y_max + 0.5)
	)
	jitter=0.04

	profs = OrderedDict(
		"J+24 probability cut" => prof_jax_cut,
		"J+24 weighted" => prof_jax_weighted,
		"MCMC hist" => prof_mean,
	)

	
	for (label, prof) in profs

		filt = isfinite.(prof.log_Sigma)
		x = prof.log_R[filt]
		y = prof.log_Sigma[filt] 

		yerr = error_interval.(y)
		errorscatter!(x .+ jitter * LilGuys.randu(-1, 1), y, yerror=yerr, label=label, )	
	end

	axislegend(ax, position=:lb)

	
	ax_res = Axis(fig[2,1],
		xlabel = log_r_label,
		ylabel = "residual",
		limits=(nothing, nothing, -0.5, 0.5)
	)

	prof_ref = prof_mean
	

	for (label, prof) in profs
		filt = isfinite.(prof_ref.log_Sigma)

		filt .&= isfinite.(prof.log_Sigma)
		x = prof.log_R[filt]
		y = prof.log_Sigma[filt] - prof_ref.log_Sigma[filt]
		yerr = LilGuys.error_interval.(prof.log_Sigma[filt])

		errorscatter!(x .+ jitter * LilGuys.randu(-1, 1), y, yerror=yerr, label=label)
	end

	rowsize!(fig.layout, 2, Relative(0.3))
	rowgap!(fig.layout, 0)
	hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
	linkxaxes!(ax_res, ax)

	@savefig "density_j+24_vs_mcmc"
end

# ╔═╡ 9bf31971-1c10-4020-946e-d8cfebf6596a
md"""
# Outputs
"""

# ╔═╡ 8c565a28-84bc-4bc7-8a0e-b2e0dff76665
if !isdir(outdir)
	mkdir(outdir)
end

# ╔═╡ a93579f2-d37d-492c-95e6-02bd7491dc8c
profout = joinpath(outdir, "hist_fast_profile.toml")

# ╔═╡ 28b67692-45a9-435d-8862-882b0f06f947
starsout = joinpath(outdir, "stars$FIGSUFFIX.fits")

# ╔═╡ 6f2a4b25-b9ea-4560-ba6e-fc22af9f23e6
open(profout, "w") do f
	print(f, prof_mean)
end

# ╔═╡ 0e6fe0e9-51d9-4948-9eca-278b82611dea
import Statistics: median

# ╔═╡ 8b05449b-6d03-40cc-b4d8-9230303c15c1
let
	stars_new = copy(stars)
	stars_new[!, :PSAT] = middle.(psats_pm)
	stars_new[!, :PSAT_em] = lower_error.(psats_pm)
	stars_new[!, :PSAT_ep] =  upper_error.(psats_pm)
	stars_new[!, :f_sat] = median.(eachrow(fsats))

	write_fits(starsout, stars_new, overwrite=true)

	stars_new
end	

# ╔═╡ fcea9627-7150-4c6b-86c7-c73671554926
md"""
# Additional checks
"""

# ╔═╡ 7a490042-abb2-424e-b61b-e78ed885034a
scatter(df_summary.N_stars, _prof_counts.counts .- df_summary.N_stars)

# ╔═╡ Cell order:
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═e997f77b-02d5-406a-8f78-b9e115baaf88
# ╠═1f497448-9ea2-460f-b403-618b78b47565
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╠═b94c3346-bd31-409e-ad6f-5e7afb891ad1
# ╠═11d9d9f5-0fb7-4b06-bd4a-36bdb1f00188
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╟─04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═b7b4361a-1ee8-4b1d-a62e-9f7bd0363160
# ╠═aae56f16-2433-4847-83db-9184a74f28f8
# ╠═4b77612f-e124-4d77-974c-d40c2f5a37ff
# ╠═9726210e-278f-4478-9c4a-8be8975a57d6
# ╠═c58ed939-07c1-4d34-90a4-3c5d9cc9910b
# ╠═6f016a8e-38ae-4f05-a7ee-c292ac0e5741
# ╠═3e99a9b9-ba06-47ee-bce7-5720a4eb1944
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╟─3067ae10-854f-48ce-8285-0d3f3fa7f25c
# ╠═d254b67b-8017-4bbf-ba04-50b5b76f48f0
# ╟─4b644aaf-76ab-4e30-8a78-7e127c2ba7b5
# ╠═85f59eb3-4f90-498a-923e-36cd679af1c0
# ╠═575a31e1-df22-42e9-9a9e-1dc4b678a5b6
# ╠═f0afb784-c1f4-40f6-84e6-a5c8f12ac67b
# ╠═fe3d955f-c96b-4202-9218-f749a5e9d41f
# ╠═5feaa9b0-ca44-48d6-80d1-16303a34005f
# ╠═ab124afe-c94a-44d9-8869-bd4d4cee3fbd
# ╟─6685cdec-7c99-45ca-98c8-f6439dfd7290
# ╠═6470edfd-bd06-43b7-aa5e-a80645af78fb
# ╠═c6a8bf2b-8a65-48ff-a1ca-7790f7806a94
# ╠═93fe93ec-b98c-4037-9504-bf584e31b420
# ╠═f80631b5-aa7e-454a-937f-c2e66b042105
# ╠═4106ea1f-948f-42fb-9a19-64580a94ce9a
# ╟─b8ba7f13-6a1a-43d5-90db-930196dfbbe4
# ╠═da737b3e-f28a-4917-91d6-eea97458ddb0
# ╠═ad85b789-8f43-489b-a41c-a966e08d78ad
# ╠═f2c43c30-b069-4d17-8d61-be801d85b245
# ╠═c507975a-8041-48f0-b46e-ac1ebccf3a94
# ╠═bab50892-fab8-44fe-931b-66ff466b6974
# ╠═7c02d24f-7d59-410b-99bc-93513afa4c5d
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═8c565a28-84bc-4bc7-8a0e-b2e0dff76665
# ╠═a93579f2-d37d-492c-95e6-02bd7491dc8c
# ╠═28b67692-45a9-435d-8862-882b0f06f947
# ╠═6f2a4b25-b9ea-4560-ba6e-fc22af9f23e6
# ╠═0e6fe0e9-51d9-4948-9eca-278b82611dea
# ╠═8b05449b-6d03-40cc-b4d8-9230303c15c1
# ╟─fcea9627-7150-4c6b-86c7-c73671554926
# ╠═7a490042-abb2-424e-b61b-e78ed885034a
