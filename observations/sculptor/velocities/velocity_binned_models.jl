### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	using Arya
	using CairoMakie

	import LilGuys as lguys
end

# ╔═╡ 93838644-cad6-4df3-b554-208b7afeb3b8
using PyFITS

# ╔═╡ 72f1febc-c6ea-449a-8cec-cd0e49c4e20c
using DataFrames

# ╔═╡ 9070c811-550c-4c49-9c58-0943b0f808b2
using Turing

# ╔═╡ bd6dfd17-02ee-4855-be37-fecfdab6776f
using LilGuys; FIGDIR = "figures"

# ╔═╡ 3ed8c28f-5908-42dc-a56b-24a9b2685a07
md"""
## Inputs
"""

# ╔═╡ 50488b8f-6886-4191-8778-af66929f1445
rv_file = "rv_combined_x_wide_2c_psat_0.2.fits"

# ╔═╡ 4827c99d-6326-42fb-a329-d56992f6d20d
j24_sample = "wide_2c"

# ╔═╡ 8b3ad5b9-0ab3-4349-90d0-013ac96ff6b1
n_samples = 10000

# ╔═╡ 0b8f0193-f1e5-471e-8c34-498ac827be63
n_threads = 16

# ╔═╡ ad79710a-0392-4b75-94ac-23881376a70b
sampler = NUTS(0.65)

# ╔═╡ 680e7f76-cb4d-40d6-9a9f-d4672427a633
md"""
## derived
"""

# ╔═╡ 86fe351f-ef12-474a-85cc-c10c22a65e77
FIGSUFFIX  = "." * splitext(basename(rv_file))[1]

# ╔═╡ 7330c75e-1bf9-476a-8274-ebc86d555e6f
md"""
# RV sample models
"""

# ╔═╡ 6bec7416-40c8-4e2b-9d3d-14aa19e5642d
md"""
This notebook takes the dataset created from velocity_xmatch.jl and analyzes it using MCMC to estimate the velocity dispersion and search for any possible velocity gradients.
"""

# ╔═╡ 34e43f4a-bcff-41cb-92c4-0c8d600fd053
import CSV

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, sem

# ╔═╡ c3298b45-7c5d-4937-8e4c-c87de36a1354
import DensityEstimators: histogram, bins_equal_number

# ╔═╡ d2888213-61e3-4a6f-872b-48a075640ef5
import TOML

# ╔═╡ b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
CairoMakie.activate!(type=:png)

# ╔═╡ 2ab0018f-1628-4f5e-b7da-370eb20c00d0
module RVUtils
	include("../../rv_utils.jl")
end

# ╔═╡ 5ec475a1-14bb-40f6-856a-69fa9efe087a
⊕ = RVUtils.:⊕

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "processed"

# ╔═╡ 1bc6b7f5-4884-479d-b4be-54f28c2e0a8a
f_sat = TOML.parsefile("../j+24_fsat.toml")[j24_sample]["f_sat"]

# ╔═╡ 3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
obs_properties = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml")

# ╔═╡ 66c35421-f6d3-4b2d-86e4-319f5476b222
σv = obs_properties["sigma_v"]

# ╔═╡ 2d110151-f7f3-4b09-8684-a2fc4327814b
Δv_gsr = RVUtils.rv_gsr_shift(obs_properties["ra"], obs_properties["dec"])

# ╔═╡ 84509e42-8484-410a-8a76-38473b9f4b71
rv0 = obs_properties["radial_velocity"] .- Δv_gsr

# ╔═╡ 9e2420ea-8d47-4eab-a4bd-0caeb09d9ebb
θ_orbit = obs_properties["theta_pm_gsr"]

# ╔═╡ 4dac920b-8252-48a7-86f5-b9f96de6aaa0
memb_stars = let
	rv_meas = read_fits(joinpath(data_dir, rv_file))

	RVUtils.add_gsr!(rv_meas, 
					 distance=obs_properties["distance"],
					pmra=obs_properties["pmra"], pmdec=obs_properties["pmdec"])


	rv_meas[:, :vz] .= rv_meas.radial_velocity_gsr .+ rv_meas.delta_rv
	rv_meas[:, :vz_err] .= rv_meas.RV_err .⊕ rv_meas.delta_rv_err

	RVUtils.add_PSAT_RV!(rv_meas; sigma_v=σv, radial_velocity_gsr=rv0, f_sat=f_sat)
	
	rv_meas
end

# ╔═╡ c28071fe-6077-43ad-b930-604483d5eb28
md"""
## Membership
"""

# ╔═╡ 7f4a5254-ed6f-4faa-a71e-4b4986a99d45
hist(memb_stars.RV)

# ╔═╡ abbd2a53-e077-4af7-a168-b571e1a906b8
xlabel = L"radial velocity / km s$^{-1}$"

# ╔═╡ 95ef12d2-174a-47c3-a106-9b4005b2a80d
md"""
# Full fits to data
"""

# ╔═╡ 9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
function fit_rv_sigma(rv, rv_err; μ_0_prior=0, N=3_000, p=0.16)
	samples = DataFrame(sample(RVUtils.model_vel_1c(rv, rv_err, μ_0_prior=μ_0_prior), sampler, MCMCThreads(), N, n_threads))

	μ = median(samples.μ)
	μ_p = quantile(samples.μ, [p, 1-p]) 
	μ_err = (μ - μ_p[1], μ_p[2] - μ)

	σ = median(samples.σ)
	σ_p = quantile(samples.σ, [p, 1-p])
	σ_err = (σ - σ_p[1], σ_p[2] - σ)

	return μ, σ, μ_err, σ_err
end

# ╔═╡ 965a9390-3088-42dd-8a2c-328f366ee2b0
df_Rell = CSV.read("processed/mcmc_samples_Rell_sigma$FIGSUFFIX.csv", DataFrame)

# ╔═╡ 8555e608-53c1-40d3-b21e-413af8953c30
md"""
code below validates induced PM gradient (should be approx 2).
"""

# ╔═╡ 43175a56-4c8e-4274-baef-db76f6d66bea
df_gradient = CSV.read("processed/mcmc_samples_gradient$FIGSUFFIX.csv", DataFrame)

# ╔═╡ 7a9d2e0f-4dcb-4b69-9f68-d43d6dde8bf2
θ_m = median(df_gradient.Θ_grad)

# ╔═╡ 3195286d-d85f-43a3-aa25-dae2134f570b
xi_rot, eta_rot = lguys.to_orbit_coords(memb_stars.ra, memb_stars.dec, obs_properties["ra"], obs_properties["dec"], θ_m) .* 60


# ╔═╡ 88f2918e-e126-420a-96a2-5746a8010f73
icrs0 = lguys.ICRS(obs_properties)

# ╔═╡ c48bb30d-1186-4940-b061-91f53e8335e1
vec_pm = LilGuys.pm2kms.([icrs0.pmra, icrs0.pmdec], icrs0.distance) /(180/π)

# ╔═╡ af0d2050-b42e-4a7f-aabb-5b08d23381e9
lguys.transform(ICRS, lguys.GSR(ra=icrs0.ra + 2/vec_pm[1] *cos(icrs0.dec), dec=icrs0.dec + 2/vec_pm[2], distance=icrs0.distance, radial_velocity=0)).radial_velocity  .- Δv_gsr

# ╔═╡ 4a473039-79f0-4d77-aa0c-681e2fba4f4c
gsr0 = lguys.transform(lguys.GSR, icrs0)

# ╔═╡ 7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
md"""
## Binned properties along orbit
"""

# ╔═╡ c2735c49-2892-46ac-bcf8-7cdcef409f44
function calc_binned_mu_sigma(x, y, yerr, bins; kwargs...)

	if !issorted(bins)
		error("bins must be sorted")
	end

	N = length(bins) - 1
	μs = Vector{Float64}(undef, N)
	σs = Vector{Float64}(undef, N)
	μ_errs = Vector{Tuple{Float64, Float64}}(undef, N)
	σ_errs = Vector{Tuple{Float64, Float64}}(undef, N)
	xs = Vector{Float64}(undef, N)

	
	for i in 1:N
		filt = x .>= bins[i]
		filt .&= x .< bins[i+1]
		@info "calculating bin $i"

		μs[i], σs[i], μ_errs[i], σ_errs[i] = fit_rv_sigma(y[filt], yerr[filt]; kwargs...)
		xs[i] = median(x[filt])
	end

	return DataFrame(
		x=xs,
		x_low = bins[1:end-1],
		x_high = bins[2:end],
		μ=μs, 
		σ = σs, 
		μ_em = first.(μ_errs), 
		μ_ep = last.(μ_errs), 
		σ_em = first.(σ_errs),
		σ_ep = last.(σ_errs),
		
	)
end	

# ╔═╡ 8b21cc49-ca17-4844-8238-e27e9752bee7
bins = bins_equal_number(memb_stars.R_ell, n=10)

# ╔═╡ f6d0dc0a-3ae2-4382-8aef-bfc816cdb721
df_r_ell_z = calc_binned_mu_sigma(memb_stars.R_ell, memb_stars.vz, memb_stars.vz_err, bins)

# ╔═╡ 03928d42-ee61-4480-82d3-a482842f8521
CSV.write("processed/vz_r_ell_binned$FIGSUFFIX.csv", df_r_ell_z)

# ╔═╡ 38da4da1-74f5-4661-89f2-4b25562a1faf
function scatter_range!(df_r_ell)
	errorscatter!(df_r_ell.x, df_r_ell.μ, yerror=collect(zip(df_r_ell.μ_em, df_r_ell.μ_ep)), color=:black)
	
	errorbars!(df_r_ell.x, df_r_ell.μ .+ df_r_ell.σ, df_r_ell.x .- df_r_ell.x_low, df_r_ell.x_high .- df_r_ell.x,  direction = :x, color=:black)
	errorbars!(df_r_ell.x, df_r_ell.μ .- df_r_ell.σ, df_r_ell.x .- df_r_ell.x_low, df_r_ell.x_high .- df_r_ell.x, direction = :x, color=:black)
end

# ╔═╡ e05ec20a-3165-4360-866e-3e8cae8665e5
let
	fig, ax = FigAxis(
		xlabel = "R / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)

	scatter!(memb_stars.R_ell, memb_stars.vz, color=COLORS[3], alpha=0.1)

	scatter_range!(df_r_ell_z)

	fig
end

# ╔═╡ e50a1968-9f91-4751-a7a9-0a4768c9ab7e
let
	fig, ax = FigAxis(
		xlabel = L"$\log\ R$ / arcmin",
		ylabel = L"$\log\ \sigma_{v, \textrm{los}}$ / km s$^{-1}$"
	)

	x_mid = midpoints(log10.(bins))
	log_bin_errs = [(x_mid[i] - log10(bins[i]), log10(bins[i+1]) - x_mid[i]) for i in eachindex(x_mid)]
	@info log_bin_errs
	
	errorscatter!(x_mid, log10.(df_r_ell_z.σ), yerror=max.(df_r_ell_z.σ_em, df_r_ell_z.σ_ep) ./ df_r_ell_z.σ ./ log(10), xerror=log_bin_errs, color=:black)
	# hlines!(log10.(σ_m), color=:black)

	for i in 1:400:size(df_Rell, 1)
		x=LinRange(-0.5, 2, 100)
		y = df_Rell.σ[i] .* 10 .^ ((x .- 1.0) .* df_Rell.dlσ_dlR[i])
		lines!(x, log10.(y), alpha=0.03, color=COLORS[1])
	end
	@savefig "log_sigma_log_R"
	fig
end

# ╔═╡ 31aa8fc5-1415-4c44-9b92-a7d097181639
md"""
# Binned properties with radius
"""

# ╔═╡ 2f5db664-fc99-41ac-a726-f21dd5d88ad4
df_xi_p = calc_binned_mu_sigma(xi_rot, memb_stars.vz, memb_stars.vz_err, bins_equal_number(xi_rot, n=10))

# ╔═╡ 301ad1d5-83a1-418d-aa06-583f7212541d
CSV.write("processed/vz_xi_p_binned$FIGSUFFIX.csv", df_xi_p)

# ╔═╡ 090acae4-1209-49e6-882c-20ac2c972dd5
df_eta_p = calc_binned_mu_sigma(eta_rot, memb_stars.vz, memb_stars.vz_err, bins_equal_number(eta_rot, n=10))

# ╔═╡ fd0a74a1-6513-4612-8181-745d5b7c3f4c
let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$ / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)

	scatter!(memb_stars.xi, memb_stars.radial_velocity_gsr, color=COLORS[3], alpha=0.1)

	scatter_range!(df_xi_p)

	fig
end

# ╔═╡ 3b411c40-2436-4d1f-bf41-8f2c3f0bf3a4
@savefig "vel_gradient_binned" let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$ / degree",
		ylabel = L"$\mu_{v,z}$ / km s$^{-1}$"
	)

	errorscatter!(df_xi_p.x ./ 60, df_xi_p.μ, yerror=collect(zip(df_xi_p.μ_em, df_xi_p.μ_ep)), color=:black)

	for i in 1:400:size(df_gradient, 1)
		M = 60df_gradient.B[i] .* sind(θ_m) .+ 60df_gradient.A[i] .* cosd(θ_m)
		x=LinRange(-0.4, 0.4, 100)
		y = M*x  .+ df_gradient.μ[i]
		lines!(x, y, alpha=0.03, color=COLORS[1])
	end
			
	fig
end

# ╔═╡ bf8a9940-1c8b-4fdd-a231-1861721ea8e9
sind(θ_m), cosd(θ_m)

# ╔═╡ 14e81f66-8ad9-48d5-aa0b-a09bc2a3bf52
let
	fig, ax = FigAxis(
		xlabel = L"$\eta'$ / arcmin",
		ylabel = L"$\mu_{v, \textrm{gsr}}$ / km s$^{-1}$"
	)

	errorscatter!(df_eta_p.x, df_eta_p.μ, yerror=collect(zip(df_xi_p.μ_em, df_xi_p.μ_ep)), color=:black)

	fig
end

# ╔═╡ d8d85c45-67ca-4c7a-92c1-a63bbf873c3d
let
	fig, ax = FigAxis(
		xlabel = L"$\xi'$ / arcmin",
		ylabel = L"$\sigma_{v, \textrm{gsr}}$ / km s$^{-1}$"
	)

	errorscatter!(df_xi_p.x, df_xi_p.σ, yerror=collect(zip(df_xi_p.σ_em, df_xi_p.σ_ep)), color=:black)


		
	@savefig "sigma_v_vs_eta_p"
	fig
end

# ╔═╡ 1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
bin_errs = diff(bins) / 2

# ╔═╡ f82d2ff7-7a7f-4520-811e-126f3f4f5349
let
	fig, ax = FigAxis(
		xlabel = L"$R_\textrm{ell}$ / arcmin",
		ylabel = L"$\sigma_{v, \textrm{los}}$ / km s$^{-1}$"
	)

	errorscatter!(midpoints(bins), df_r_ell_z.σ, yerror=collect(zip(df_r_ell_z.σ_em, df_r_ell_z.σ_ep)), xerror=bin_errs, color=:black)
	# hlines!(σ_m)

	@savefig "sigma_v_vs_rell"
	fig
end

# ╔═╡ Cell order:
# ╟─3ed8c28f-5908-42dc-a56b-24a9b2685a07
# ╠═50488b8f-6886-4191-8778-af66929f1445
# ╠═4827c99d-6326-42fb-a329-d56992f6d20d
# ╠═8b3ad5b9-0ab3-4349-90d0-013ac96ff6b1
# ╠═0b8f0193-f1e5-471e-8c34-498ac827be63
# ╠═ad79710a-0392-4b75-94ac-23881376a70b
# ╟─680e7f76-cb4d-40d6-9a9f-d4672427a633
# ╠═86fe351f-ef12-474a-85cc-c10c22a65e77
# ╟─7330c75e-1bf9-476a-8274-ebc86d555e6f
# ╟─6bec7416-40c8-4e2b-9d3d-14aa19e5642d
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═34e43f4a-bcff-41cb-92c4-0c8d600fd053
# ╠═93838644-cad6-4df3-b554-208b7afeb3b8
# ╠═72f1febc-c6ea-449a-8cec-cd0e49c4e20c
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═9070c811-550c-4c49-9c58-0943b0f808b2
# ╠═c3298b45-7c5d-4937-8e4c-c87de36a1354
# ╠═d2888213-61e3-4a6f-872b-48a075640ef5
# ╠═b00b7e2b-3a72-466f-ac09-86cdda1e4a9c
# ╠═bd6dfd17-02ee-4855-be37-fecfdab6776f
# ╠═2ab0018f-1628-4f5e-b7da-370eb20c00d0
# ╠═5ec475a1-14bb-40f6-856a-69fa9efe087a
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═1bc6b7f5-4884-479d-b4be-54f28c2e0a8a
# ╠═66c35421-f6d3-4b2d-86e4-319f5476b222
# ╠═2d110151-f7f3-4b09-8684-a2fc4327814b
# ╠═84509e42-8484-410a-8a76-38473b9f4b71
# ╠═9e2420ea-8d47-4eab-a4bd-0caeb09d9ebb
# ╠═3eb74a2e-ca74-4145-a2a4-7ffbe5fffe94
# ╠═4dac920b-8252-48a7-86f5-b9f96de6aaa0
# ╟─c28071fe-6077-43ad-b930-604483d5eb28
# ╠═7f4a5254-ed6f-4faa-a71e-4b4986a99d45
# ╠═abbd2a53-e077-4af7-a168-b571e1a906b8
# ╟─95ef12d2-174a-47c3-a106-9b4005b2a80d
# ╠═9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
# ╠═965a9390-3088-42dd-8a2c-328f366ee2b0
# ╠═c48bb30d-1186-4940-b061-91f53e8335e1
# ╟─8555e608-53c1-40d3-b21e-413af8953c30
# ╠═af0d2050-b42e-4a7f-aabb-5b08d23381e9
# ╠═43175a56-4c8e-4274-baef-db76f6d66bea
# ╠═7a9d2e0f-4dcb-4b69-9f68-d43d6dde8bf2
# ╠═3195286d-d85f-43a3-aa25-dae2134f570b
# ╠═4a473039-79f0-4d77-aa0c-681e2fba4f4c
# ╠═88f2918e-e126-420a-96a2-5746a8010f73
# ╟─7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═c2735c49-2892-46ac-bcf8-7cdcef409f44
# ╠═8b21cc49-ca17-4844-8238-e27e9752bee7
# ╠═f6d0dc0a-3ae2-4382-8aef-bfc816cdb721
# ╠═03928d42-ee61-4480-82d3-a482842f8521
# ╠═38da4da1-74f5-4661-89f2-4b25562a1faf
# ╠═e05ec20a-3165-4360-866e-3e8cae8665e5
# ╠═f82d2ff7-7a7f-4520-811e-126f3f4f5349
# ╠═e50a1968-9f91-4751-a7a9-0a4768c9ab7e
# ╟─31aa8fc5-1415-4c44-9b92-a7d097181639
# ╠═2f5db664-fc99-41ac-a726-f21dd5d88ad4
# ╠═301ad1d5-83a1-418d-aa06-583f7212541d
# ╠═090acae4-1209-49e6-882c-20ac2c972dd5
# ╠═fd0a74a1-6513-4612-8181-745d5b7c3f4c
# ╠═3b411c40-2436-4d1f-bf41-8f2c3f0bf3a4
# ╠═bf8a9940-1c8b-4fdd-a231-1861721ea8e9
# ╠═14e81f66-8ad9-48d5-aa0b-a09bc2a3bf52
# ╠═d8d85c45-67ca-4c7a-92c1-a63bbf873c3d
# ╠═1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
