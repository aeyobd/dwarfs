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
	@bind galaxy confirm(TextField(default="galaxyname"))
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
# stars = MCMCUtils.get_fits(galaxy, obs_props)

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

# ╔═╡ 4b644aaf-76ab-4e30-8a78-7e127c2ba7b5
md"""
Consistency checks:
"""

# ╔═╡ 85f59eb3-4f90-498a-923e-36cd679af1c0
@assert df_summary.log_R ≈ midpoints(log10.(bins))

# ╔═╡ f0afb784-c1f4-40f6-84e6-a5c8f12ac67b
log_Sigmas_fast = let
	Nbins = length(bins) - 1

	Nc = size(df_chains, 1)
	
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

# ╔═╡ 5feaa9b0-ca44-48d6-80d1-16303a34005f
prof_jax_cut = let
	local profile = nothing
	
	for filename in ["jax_eqw_profile", "jax_profile"]
		path = joinpath("..", galaxy, "density_profiles/", filename * ".toml")
		if isfile(path)
			profile = LilGuys.SurfaceDensityProfile(path)
			break
		end
	end
	
	profile
end

# ╔═╡ 55624f98-a60f-462c-8c76-394edb27cf51
prof_jax_2c = LilGuys.SurfaceDensityProfile(joinpath("..", galaxy, "density_profiles/jax_2c_eqw_profile.toml"))

# ╔═╡ ab124afe-c94a-44d9-8869-bd4d4cee3fbd
prof_jax_weighted = LilGuys.SurfaceDensityProfile(joinpath("..", galaxy, "density_profiles/weighted_profile.toml"))

# ╔═╡ 6685cdec-7c99-45ca-98c8-f6439dfd7290
md"""
## MCMC Profiles
"""

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

# ╔═╡ 43787c14-2d85-4f5a-8377-5999eaa49486
begin 
	profs = OrderedDict(
		"J+24 probability cut" => prof_jax_cut,
		"MCMC hist" => prof_mean,
	)

	if @isdefined prof_jax_weighted
		profs["J+24 weighted"] = prof_jax_weighted
	end

	if @isdefined prof_jax_2c
		profs["J+24 2comp"] = prof_jax_2c
	end
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


	if @isdefined prof_jax_weighted
		profs["J+24 weighted"] => prof_jax_weighted
	end

	
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

	prof_ref = LilGuys.Exp2D(R_s=obs_props["R_h"] / LilGuys.R_h(LilGuys.Exp2D()), M=sum(prof_jax_cut.counts))
	
	f(log_R) = log10(LilGuys.surface_density(prof_ref, 10^log_R))

	
	for (label, prof) in profs

		filt = isfinite.(prof.log_Sigma)
		x = prof.log_R[filt]
		y = prof.log_Sigma[filt] 

		yerr = error_interval.(y)
		errorscatter!(x .+ jitter * LilGuys.randu(-1, 1), y, yerror=yerr, label=label, )	
	end

	lines!(LinRange(-0.5, 2, 100), f, color=:black)

	axislegend(ax, position=:lb)

	
	ax_res = Axis(fig[2,1],
		xlabel = log_r_label,
		ylabel = "residual",
		limits=(nothing, nothing, -2, 2)
	)


	lines!(LinRange(-0.5, 2, 100), x->0, color=:black)


	for (label, prof) in profs
		filt = isfinite.(prof.log_Sigma)
		x = prof.log_R[filt]
		y = prof.log_Sigma[filt] - f.(x)
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

# ╔═╡ 6f2a4b25-b9ea-4560-ba6e-fc22af9f23e6
open(profout, "w") do f
	print(f, prof_mean)
end

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
# ╟─4b644aaf-76ab-4e30-8a78-7e127c2ba7b5
# ╠═85f59eb3-4f90-498a-923e-36cd679af1c0
# ╠═f0afb784-c1f4-40f6-84e6-a5c8f12ac67b
# ╠═5feaa9b0-ca44-48d6-80d1-16303a34005f
# ╠═55624f98-a60f-462c-8c76-394edb27cf51
# ╠═ab124afe-c94a-44d9-8869-bd4d4cee3fbd
# ╟─6685cdec-7c99-45ca-98c8-f6439dfd7290
# ╠═c6a8bf2b-8a65-48ff-a1ca-7790f7806a94
# ╠═93fe93ec-b98c-4037-9504-bf584e31b420
# ╠═f80631b5-aa7e-454a-937f-c2e66b042105
# ╠═4106ea1f-948f-42fb-9a19-64580a94ce9a
# ╟─b8ba7f13-6a1a-43d5-90db-930196dfbbe4
# ╠═f2c43c30-b069-4d17-8d61-be801d85b245
# ╠═c507975a-8041-48f0-b46e-ac1ebccf3a94
# ╠═43787c14-2d85-4f5a-8377-5999eaa49486
# ╠═bab50892-fab8-44fe-931b-66ff466b6974
# ╠═7c02d24f-7d59-410b-99bc-93513afa4c5d
# ╟─9bf31971-1c10-4020-946e-d8cfebf6596a
# ╠═8c565a28-84bc-4bc7-8a0e-b2e0dff76665
# ╠═a93579f2-d37d-492c-95e6-02bd7491dc8c
# ╠═6f2a4b25-b9ea-4560-ba6e-fc22af9f23e6
