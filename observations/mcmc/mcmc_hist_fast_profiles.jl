### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 8e3f56fa-034b-11f0-1844-a38fa58e125c
begin
	using Pkg; Pkg.activate()
	using CairoMakie
	using Arya

	using LilGuys
end

# ╔═╡ 67ddb57a-9ee7-4f4c-848f-d0b77fba1855
using DataFrames, CSV

# ╔═╡ e74e6a96-a18e-40bf-ade8-07ce0e30e5c4
using Distributions

# ╔═╡ 6c4427b6-fe6d-4ad7-812a-21c9788591d5
using OrderedCollections

# ╔═╡ 08d97b62-2760-47e4-b891-8f446e858c88
if !@isdefined(PlutoRunner)
	galaxy = ARGS[1]
else
	galaxy = "antlia2"
	skip = 10
	all_profiles = true
end

# ╔═╡ d1de613c-c3bb-4843-859e-7b8df54bafe0
import TOML

# ╔═╡ 1f497448-9ea2-460f-b403-618b78b47565
module MCMCUtils
	include("mcmc_utils.jl")
end

# ╔═╡ 11d9d9f5-0fb7-4b06-bd4a-36bdb1f00188
outdir = joinpath("..", galaxy, "mcmc")

# ╔═╡ b94c3346-bd31-409e-ad6f-5e7afb891ad1
FIGDIR = joinpath(outdir, "figures"); FIGSUFFIX=".mcmc_hist_fast"

# ╔═╡ 133a025f-407f-49eb-9e02-0c620d5b77ba
CairoMakie.activate!(type=:png)

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

# ╔═╡ c58ed939-07c1-4d34-90a4-3c5d9cc9910b
stars = MCMCUtils.get_fits(galaxy, obs_props)

# ╔═╡ f078d929-62be-48f2-acdb-808357801a7b
md"""
note, we do change the first and last bins to be something more workable (just double the last bin)
"""

# ╔═╡ 4b77612f-e124-4d77-974c-d40c2f5a37ff
struct_params = MCMCUtils.StructuralParams(; LilGuys.dict_to_tuple(TOML.parsefile(joinpath("..", galaxy, "mcmc", "default_bins.toml")))...)

# ╔═╡ 9726210e-278f-4478-9c4a-8be8975a57d6
bins = struct_params.bins

# ╔═╡ 6f016a8e-38ae-4f05-a7ee-c292ac0e5741
df_chains = CSV.read(joinpath(outdir, "samples.mcmc_hist_fast.csv"), DataFrame)

# ╔═╡ 3e99a9b9-ba06-47ee-bce7-5720a4eb1944
df_summary = CSV.read(joinpath(outdir, "summary.mcmc_hist_fast.csv"), DataFrame)

# ╔═╡ 24a65d65-6e0b-4108-8041-79fee06cd28a
md"""
# Analysis
"""

# ╔═╡ 7520415d-a482-4a9a-bf8a-7f8443d3b613
md"""
### Reference profiles
"""

# ╔═╡ 5feaa9b0-ca44-48d6-80d1-16303a34005f
prof_simple = LilGuys.StellarDensityProfile(stars.R_ell[stars.PSAT .> 0.2], bins=log10.(bins), errors=:weighted, normalization=:none)

# ╔═╡ ab124afe-c94a-44d9-8869-bd4d4cee3fbd
prof_weighted = LilGuys.StellarDensityProfile(stars.R_ell, weights=stars.PSAT, bins=log10.(bins), errors=:weighted,  normalization=:none)

# ╔═╡ babcaaa8-8d95-4dcb-b22b-5e36c10c57dd
prof_bernoulli = LilGuys.StellarDensityProfile(stars.R_ell, weights=stars.PSAT, bins=log10.(bins), errors=:bernoulli,  normalization=:none)

# ╔═╡ e60ea082-818d-46f8-bf50-30a932f97ef2
prof_auto = LilGuys.StellarDensityProfile(stars.R_ell[stars.PSAT .> 0.2], errors=:weighted, normalization=:none)

# ╔═╡ 82f20e46-90c7-4925-83c9-3f49a909664d
sum(prof_simple.counts)

# ╔═╡ 6685cdec-7c99-45ca-98c8-f6439dfd7290
md"""
## MCMC Profiles
"""

# ╔═╡ d254b67b-8017-4bbf-ba04-50b5b76f48f0
prof_counts = LilGuys.StellarDensityProfile(stars.R_ell, bins=log10.(bins),  errors=:weighted)

# ╔═╡ f0afb784-c1f4-40f6-84e6-a5c8f12ac67b
let
	skip = 1
	Nbins = length(bins) - 1

	Nc = size(df_chains, 1)
	Ns = size(stars, 1)
	
	global log_Sigmas_fast = Matrix{Float64}(undef, Nbins, Nc)
	areas = diff(π * bins .^2)
	counts = prof_counts.counts

	for i in 1:skip:Nc
		if i % 200 == 0
			@info "calculated $i / $Nc"
		end
		
		params = [df_chains[i, "params[$j]"] for j in 1:Nbins]
		f_sat = exp10.(params) ./(1 .+ exp10.(params))
		log_Sigmas_fast[:, i] .= log10.( f_sat .* counts ./ areas)
	end

end

# ╔═╡ 2027b636-bb17-4004-941e-b715ba8d0216
psat_min = 1e-10

# ╔═╡ 8abe3d92-f9e8-4e89-96ec-14eaf9482e64
Nsteps = size(df_chains, 1)

# ╔═╡ 939aa449-26ea-4d2d-8773-42a73c3d0410
if all_profiles
	Nbins = length(bins) - 1

	Nc = size(df_chains, 1)
	Ns = size(stars, 1)

	idxs = 1:skip:Nc
	fsats = Matrix{Float64}(undef, Ns, length(idxs))
	psats = Matrix{Float64}(undef, Ns, length(idxs))
	profiles = Vector{LilGuys.StellarDensityProfile}(undef, length(idxs))

	Threads.@threads for i in eachindex(idxs)
		idx = idxs[i]
		if i % 10 == 0
			@info "sample $i"
		end
		
		params = [df_chains[idx, "params[$j]"] for j in 1:Nbins]
		radii = MCMCUtils.perturbed_radii(stars, struct_params)

		Σ = MCMCUtils.Σ_hist(radii, bins, params)
		f = 10 .^ Σ ./ (1 .+ 10 .^ Σ)
		fsats[:, i] .= f

		Ls = stars.L_CMD_SAT .* stars.L_PM_SAT
		Lb = stars.L_CMD_BKD .* stars.L_PM_BKD
		psat =  @. f*Ls / (f*Ls + (1-f) * Lb)
		psats[:, i] .= psat

		weights = zeros(length(psat))
		filt = psat .> 0
		weights[filt] = rand(Dirichlet(psat[filt]))
		weights .*= sum(psat)
		prof = LilGuys.StellarDensityProfile(radii, bins=log10.(bins), 
			weights=weights, errors=:weighted)

		profiles[i] = prof
	end

end

# ╔═╡ ba280361-c910-483a-9108-fd4be4906606
sum(middle.(10 .^ profiles[1].log_Sigma) .* diff(exp10.(profiles[1].log_R_bins) .^ 2) .* π)

# ╔═╡ 36f7e33d-0ee3-4ea3-bb0f-3cedd34da773
sum(psats[:, 1])

# ╔═╡ 36bb1876-5937-44fe-bb1c-add10470371d
pvalue = 0.16

# ╔═╡ c48f8ad0-a8c8-46ef-9b56-6f4ff5ac7a83
all_Sigmas = hcat([prof.log_Sigma for prof in profiles]...)

# ╔═╡ 084a491c-e53a-4371-9c5b-f9450ec5deec
sum(isnan.(all_Sigmas))

# ╔═╡ 428d53fd-39cd-40f3-ac6b-389f637e351e
all_Sigmas[isnan.(all_Sigmas)] .= -Inf

# ╔═╡ f4d9ddee-9ca9-4e5a-9d6c-a18e866ba71a
all_Sigmas_m = middle.(all_Sigmas)

# ╔═╡ 053d441b-9b1b-4c1a-bc96-fa5b1e9f059e
function to_measurements(A::Matrix; pvalue=pvalue)
	m = dropdims(median(A, dims=2), dims=2)

	h = quantile.(eachrow(A), 1-pvalue)

	l = quantile.(eachrow(A), pvalue)

	return Measurement.(m, m .- l, h .- m)
end

# ╔═╡ 6470edfd-bd06-43b7-aa5e-a80645af78fb
psats_pm = to_measurements(psats)

# ╔═╡ 5f5bcf71-f306-45db-ba76-af69e88675c2
log_Sigma_slow = to_measurements(middle.(all_Sigmas))

# ╔═╡ c6a8bf2b-8a65-48ff-a1ca-7790f7806a94
log_Sigma_fast = to_measurements(log_Sigmas_fast)

# ╔═╡ 81e5a128-c1ed-4a8a-8bb6-e4efd0593c80
sum(stars.PSAT)

# ╔═╡ da737b3e-f28a-4917-91d6-eea97458ddb0
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "PSAT J+24", ylabel = "PSAT MCMC")

	skip = 100

	for i in 1:skip:size(psats, 2)
		y = psats[:, i]
		filt = 1e-6 .< y .< 1-1e-6
		scatter!(stars.PSAT[filt], y[filt], color=:black, alpha=0.1, markersize=2, rasterize=true)
	end
	@savefig "j+24_vs_mcmc"
	
	fig
end

# ╔═╡ b2b89764-7990-4441-bd3c-5f369fbbaceb
size(psats)

# ╔═╡ ad85b789-8f43-489b-a41c-a966e08d78ad
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log r_ell", ylabel = "psat me - jax")

	skip = 10

	scatter!(log10.(stars.R_ell), middle.(psats_pm) .- stars.PSAT, color=:black, alpha=0.1, markersize=2, rasterize=true)

	
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
	skip = 10
	
	for h in profiles[1:skip:end]
		filt = isfinite.(h.log_Sigma)
		scatter!(h.log_R[filt] .+ 0.01*randn(sum(filt)), h.log_Sigma[filt], color=:black, alpha=0.1, markersize=1)
	end

	@savefig "scatter_profiles"

	fig
end

# ╔═╡ 08424d6a-6a21-45d2-9b1a-d28e129c8158
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = L"\Gamma",
		limits = (nothing, nothing, -10, 2)
	)
	skip = 10
	
	for h in profiles[1:skip:end]

		filt = isfinite.(h.Gamma)
		scatter!(h.log_R[filt].+ 0.01*randn(sum(filt)), h.Gamma[filt], color=:black, alpha=0.1, markersize=1)
	end


	fig
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

# ╔═╡ 4482bead-75a5-413f-b442-e60ab3d3ec88
prof_mc_slow = LilGuys.StellarDensityProfile(
	R_units = "arcmin",
	log_Sigma=log_Sigma_slow,
	log_R_bins=log10.(bins),
	log_R=midpoints(log10.(bins)),
)

# ╔═╡ 94d6a13f-1165-4d8f-a40a-a4cfa9224093
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel =log_Sigma_label,
		limits = (nothing, nothing, -8, 3)
	)
	skip = 10
	
	for i in 1:skip:Nsteps
		y = log_Sigmas_fast[2:end-1, i]
		x = prof_counts.log_R[2:end-1] .+ 0.01*randn(length(prof_counts.log_R) - 2)
		
		filt = isfinite.(y)
		scatter!(x[filt], y[filt], color=:black, alpha=0.1, markersize=1)
	end
	
	if all_profiles
		LilGuys.plot_log_Σ!(ax, prof_mc_slow, label="mc hist", size=2, color=(COLORS[2], 0.3))
	end
	
	@savefig "scatter_profiles_approx"

	fig
end

# ╔═╡ bf8d5936-1a4e-47c2-bb22-531ab344b8ad
prof_mc = LilGuys.StellarDensityProfile(
		R_units = "arcmin",
		log_Sigma=log_Sigma_fast,
		log_R_bins=log10.(bins), # use normal bins for these
		log_R=midpoints(log10.(bins)),
	)

# ╔═╡ 7c02d24f-7d59-410b-99bc-93513afa4c5d
let
	fig = Figure(size=(5*72, 3*72))
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
		limits = (nothing, nothing, -12, 4)
	)
	jitter=0.01

	profs = OrderedDict(
		"mc hist" => prof_mc_slow,
		"mc hist approx" => prof_mc,
		"weighted" => prof_weighted,
		"fiducial" => prof_simple,
	)

	
	for (label, prof) in profs
		LilGuys.plot_log_Σ!(ax, prof, label=label)
	end

	Legend(fig[1,2], ax)

	
	ax_res = Axis(fig[2,1],
		xlabel = log_r_label,
		ylabel = "residual",
		limits=(nothing, nothing, -0.5, 0.5)
	)

	prof_ref = prof_mc_slow
	

	for (label, prof) in profs
		filt = isfinite.(prof_ref.log_Sigma)

		filt .&= isfinite.(prof.log_Sigma)
		x = prof.log_R[filt]
		y = prof.log_Sigma[filt] - prof_ref.log_Sigma[filt]
		yerr = LilGuys.error_interval.(prof.log_Sigma[filt])

		errorscatter!(x, y, yerror=yerr, label=label)
	end

	rowsize!(fig.layout, 2, Relative(0.3))
	rowgap!(fig.layout, 0)
	hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
	linkxaxes!(ax_res, ax)

	@savefig "density_versus_simple"
end

# ╔═╡ 6f2a4b25-b9ea-4560-ba6e-fc22af9f23e6
open(profout, "w") do f
	print(f, prof_mc_slow)
end

# ╔═╡ 17b68535-445d-4203-861f-f278090de89a
open(joinpath(outdir, "hist_fast_approx_profile.toml"), "w") do f
	print(f, prof_mc)
end

# ╔═╡ 8b05449b-6d03-40cc-b4d8-9230303c15c1
let
	stars_new = copy(stars)
	stars_new[!, :PSAT] = middle.(psats_pm)
	stars_new[!, :PSAT_em] = lower_error.(psats_pm)
	stars_new[!, :PSAT_ep] =  upper_error.(psats_pm)
	stars_new[!, :f_sat] = median.(eachrow(fsats))

	LilGuys.write_fits(starsout, stars_new, overwrite=true)
end	

# ╔═╡ Cell order:
# ╠═08d97b62-2760-47e4-b891-8f446e858c88
# ╠═8e3f56fa-034b-11f0-1844-a38fa58e125c
# ╠═d1de613c-c3bb-4843-859e-7b8df54bafe0
# ╠═67ddb57a-9ee7-4f4c-848f-d0b77fba1855
# ╠═e74e6a96-a18e-40bf-ade8-07ce0e30e5c4
# ╠═6c4427b6-fe6d-4ad7-812a-21c9788591d5
# ╠═1f497448-9ea2-460f-b403-618b78b47565
# ╠═11d9d9f5-0fb7-4b06-bd4a-36bdb1f00188
# ╠═b94c3346-bd31-409e-ad6f-5e7afb891ad1
# ╠═133a025f-407f-49eb-9e02-0c620d5b77ba
# ╠═0aa44389-f320-4274-abdd-0d7f88006a4d
# ╠═36ef5f0d-3a14-4e5d-a2fb-f6626a941d5f
# ╟─04053b71-bd55-40d7-885d-6df67035e3d6
# ╠═b7b4361a-1ee8-4b1d-a62e-9f7bd0363160
# ╠═c58ed939-07c1-4d34-90a4-3c5d9cc9910b
# ╠═f078d929-62be-48f2-acdb-808357801a7b
# ╠═4b77612f-e124-4d77-974c-d40c2f5a37ff
# ╠═9726210e-278f-4478-9c4a-8be8975a57d6
# ╠═6f016a8e-38ae-4f05-a7ee-c292ac0e5741
# ╠═3e99a9b9-ba06-47ee-bce7-5720a4eb1944
# ╟─24a65d65-6e0b-4108-8041-79fee06cd28a
# ╠═7520415d-a482-4a9a-bf8a-7f8443d3b613
# ╠═5feaa9b0-ca44-48d6-80d1-16303a34005f
# ╠═ab124afe-c94a-44d9-8869-bd4d4cee3fbd
# ╠═babcaaa8-8d95-4dcb-b22b-5e36c10c57dd
# ╠═e60ea082-818d-46f8-bf50-30a932f97ef2
# ╠═82f20e46-90c7-4925-83c9-3f49a909664d
# ╟─6685cdec-7c99-45ca-98c8-f6439dfd7290
# ╠═d254b67b-8017-4bbf-ba04-50b5b76f48f0
# ╠═f0afb784-c1f4-40f6-84e6-a5c8f12ac67b
# ╠═2027b636-bb17-4004-941e-b715ba8d0216
# ╠═8abe3d92-f9e8-4e89-96ec-14eaf9482e64
# ╠═939aa449-26ea-4d2d-8773-42a73c3d0410
# ╠═ba280361-c910-483a-9108-fd4be4906606
# ╠═36f7e33d-0ee3-4ea3-bb0f-3cedd34da773
# ╠═36bb1876-5937-44fe-bb1c-add10470371d
# ╠═c48f8ad0-a8c8-46ef-9b56-6f4ff5ac7a83
# ╠═084a491c-e53a-4371-9c5b-f9450ec5deec
# ╠═428d53fd-39cd-40f3-ac6b-389f637e351e
# ╠═f4d9ddee-9ca9-4e5a-9d6c-a18e866ba71a
# ╠═053d441b-9b1b-4c1a-bc96-fa5b1e9f059e
# ╠═6470edfd-bd06-43b7-aa5e-a80645af78fb
# ╠═5f5bcf71-f306-45db-ba76-af69e88675c2
# ╠═c6a8bf2b-8a65-48ff-a1ca-7790f7806a94
# ╠═81e5a128-c1ed-4a8a-8bb6-e4efd0593c80
# ╠═da737b3e-f28a-4917-91d6-eea97458ddb0
# ╠═b2b89764-7990-4441-bd3c-5f369fbbaceb
# ╠═ad85b789-8f43-489b-a41c-a966e08d78ad
# ╠═f2c43c30-b069-4d17-8d61-be801d85b245
# ╠═94d6a13f-1165-4d8f-a40a-a4cfa9224093
# ╠═08424d6a-6a21-45d2-9b1a-d28e129c8158
# ╠═7c02d24f-7d59-410b-99bc-93513afa4c5d
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═8c565a28-84bc-4bc7-8a0e-b2e0dff76665
# ╠═a93579f2-d37d-492c-95e6-02bd7491dc8c
# ╠═28b67692-45a9-435d-8862-882b0f06f947
# ╠═4482bead-75a5-413f-b442-e60ab3d3ec88
# ╠═bf8d5936-1a4e-47c2-bb22-531ab344b8ad
# ╠═6f2a4b25-b9ea-4560-ba6e-fc22af9f23e6
# ╠═17b68535-445d-4203-861f-f278090de89a
# ╠═8b05449b-6d03-40cc-b4d8-9230303c15c1
