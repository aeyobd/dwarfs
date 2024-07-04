### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	
	using CSV, DataFrames
	using Arya
	using GLMakie
	using FITSIO
	using Measurements

	import LilGuys as lguys
end

# ╔═╡ 9070c811-550c-4c49-9c58-0943b0f808b2
using Turing

# ╔═╡ e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
using PairPlots

# ╔═╡ 6f2359d8-332b-11ef-0db9-f1f06474c561
md"""
From Fed

Here attached some tables that I used for my paper on Sculptor. The literature compilations are from the DART survey and APOGEE DR17. There are some member stars in common between the two compilations. Therefore, I have divided into 4 files, I guess that depending on your needs, they might give you the same results.  I think all of them has been crossmatched with the Jaclyn’s member list and probabilities. I selected stars with a probability of being members of > 40 percent.

- The `Sculptor_DARTS_BEST_psat40.csv` contains all the DART members. 
- The table `Sculptor_DARTS_BEST_psat40_notAPOGEE` is the same but without the stars that are also in APOGEE. 

- `sculptor_apogeeDR17_xmatch.csv` has all the APOGEE stars,
- while DARt stars are removed in `sculptor_apogeeDR17_xmatch_notDART.csv`. 

The total members from DART and APOGEE (removing duplicates) should be 617 stars, either from `Sculptor_DARTS_BEST_psat40` + `sculptor_apogeeDR17_xmatch_notDART`  or Sculptor_DARTS_BEST_psat40_notAPOGEE + `sculptor_apogeeDR17_xmatch`.

Then, you should also include the stars observed with GMOS in the paper, they are in the table `targets_stellar_met.csv`. Since GMOS has relatively large systematics uncertainties on the RV, you should sum in quadrature to the RV_err the value 13.3 km/s (as in the paper).

Let me know if you have questions on the tables or things in Sculptor.
"""

# ╔═╡ dd29be70-4918-47e4-98cf-3608df14e88a
begin
	ell = 0.37
	PA = -94
	ra0 = 15.039170
	dec0 = -33.709180
end

# ╔═╡ 23402e30-2064-494c-bd4a-8243cb474b61
begin 
	dart = CSV.read("data/Sculptor_DARTS_BEST_psat40.csv", DataFrame)
	rename!(dart, 
		"Elliptical radius"=>"r_h",
		"vel"=>"RV",
		"evel"=>"RV_err"
	)
end

# ╔═╡ 48c62811-136f-4962-a42c-b1dd1fc74f8c
begin 
	apogee = CSV.read("data/sculptor_apogeeDR17_xmatch_notDART.csv", DataFrame)

	rename!(apogee, 
		"Elliptical half-light radii"=>"r_h",
		"VHELIO_AVG"=>"RV",
		"VERR"=>"RV_err",
		"GAIAEDR3_PMRA"=>"pmra",
		"GAIAEDR3_PMDEC"=>"pmdec",
		"RA (deg)"=>"ra",
		"Dec (deg)"=>"dec"
	)
end

# ╔═╡ b7345279-4f80-47ad-a726-537571849eae
begin 
	gmos = CSV.read("data/targets_stellar_met.csv", DataFrame)
	rename!(gmos,
		"RA"=>"ra",
		"Dec"=>"dec"
	)
end

# ╔═╡ a141ae2b-debd-453e-a7ea-d5a7287960a5
gmos_rv_sys_err = 13.3

# ╔═╡ 2513e518-cb84-4fca-90d6-2f169e14f4f1
gmos.RV_err .= @. sqrt(gmos.RV_err^2 + gmos_rv_sys_err^2)

# ╔═╡ 6734991c-16c0-4424-a2bb-84bfa811121f
md"""
## MCMC Priors
"""

# ╔═╡ 1b97d0e5-7a77-44a5-b609-ed8945cd959c
"""
Fits a normal (gaussian) distribution to 1d data with errors.
"""
@model function normal_dist(x, xerr; μ_min=90, μ_max=120)
	μ ~ Uniform(μ_min, μ_max)
	σ ~ LogNormal(2.5, 1) # approx 1 - 100 km / s, very broad but should cover all
	s = @. sqrt(σ^2 + xerr^2)

	x ~ MvNormal(fill(μ, length(x)), s)
end

# ╔═╡ 2f1d6fad-7c99-4dae-822a-702d9a74bb94
plot(Normal(2.5/log(10), 1/log(10)))

# ╔═╡ abbd2a53-e077-4af7-a168-b571e1a906b8
xlabel = L"radial velocity / km s$^{-1}$"

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile

# ╔═╡ 9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
function fit_rv_sigma(rv, rv_err; N=3_000, p=0.16)
	samples = DataFrame(sample(normal_dist(rv, rv_err), NUTS(0.65), N))

	μ = median(samples.μ)
	μ_p = quantile(samples.μ, [p, 1-p]) 
	μ_err = (μ - μ_p[1], μ_p[2] - μ)

	σ = median(samples.σ)
	σ_p = quantile(samples.σ, [p, 1-p])
	σ_err = (σ - σ_p[1], σ_p[2] - σ)

	return μ, σ, μ_err, σ_err
end

# ╔═╡ 2474644e-06cc-4ee5-9318-dc7562de6b6c
tuple([1,2]...)

# ╔═╡ d8f07585-ce89-4741-88d6-5bd3d5e4c71c
lguys.gaussian

# ╔═╡ 55504311-fbbd-4351-b9e2-ecab053ba207
function plot_samples!(samples, x;
		thin=10, color=:black, alpha=nothing, kwargs...)

	alpha = 1 / (size(samples, 1))^(1/3)
	for sample in eachrow(samples)[1:thin:end]
		y = lguys.gaussian.(x, sample.μ, sample.σ)
		lines!(x, y, color=color, alpha=alpha)
	end
end

# ╔═╡ 7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
md"""
## Binned properties
"""

# ╔═╡ c2735c49-2892-46ac-bcf8-7cdcef409f44
function binned_mu_sigma(x, y, yerr, bins)

	if !issorted(bins)
		error("bins must be sorted")
	end

	N = length(bins) - 1
	μs = Vector{Float64}(undef, N)
	σs = Vector{Float64}(undef, N)
	μ_errs = Vector{Tuple{Float64, Float64}}(undef, N)
	σ_errs = Vector{Tuple{Float64, Float64}}(undef, N)

	
	for i in 1:N
		filt = x .>= bins[i]
		filt .&= x .< bins[i+1]
		@info "calculating bin $i"

		μs[i], σs[i], μ_errs[i], σ_errs[i] = fit_rv_sigma(y[filt], yerr[filt])
	end

	return μs, σs, μ_errs, σ_errs
end	

# ╔═╡ 82a0e58a-30a4-4e42-b9c1-cb184eb551aa
md"""
# Misc plots
"""

# ╔═╡ db78b22f-bee8-4d36-994e-6c7f9e96c9f2
value(x::Measurement) = x.val

# ╔═╡ 12ad8e67-0c7b-4029-91f6-5a2453ec799a
err(x::Measurement) = x.err

# ╔═╡ e976d103-9956-47f3-adb0-7c449214b9d9
import StatsBase as sb

# ╔═╡ 2f5e5451-e653-43e8-8d2b-aeff4c575889
"""Standard error of variance"""
function sev(x)
	k = sb.kurtosis(x) + 3
	s = sb.std(x)
	n = length(x)

	var_s2 = (k - (n-3)/(n-1)) * s^4 / n
	return 1 / 2s * sqrt(var_s2)
end

# ╔═╡ f44d3c45-064c-4bf3-bf2c-e1f0d3504dd1
"""Standard error of variance"""
function sev_iid(x)
	s = sb.std(x)
	n = length(x)
	return s / sqrt(2n - 2)
end

# ╔═╡ a5e55f3e-de53-443c-9d28-e539aa6206a8
errscatter(dart.r_h, dart.feh, yerr=collect(zip(dart.feh .- dart[:, "feh_lo"], dart[:, "fe_hi"] .- dart.feh)))

# ╔═╡ 9b4a0a1f-4b0c-4c90-b871-2bd244f0a908
not = !

# ╔═╡ 1e06b622-ef1c-43ed-9453-bceab92c1889
let 
	global rv_meas = vcat(dart, apogee, gmos, cols=:union)

	filt = not.(ismissing.(rv_meas.RV))
	filt .&= rv_meas.RV .< 200
	filt .&= rv_meas.RV .> 50
	rv_meas = rv_meas[filt, :]

	rv_meas[!, :RV] = float.(rv_meas.RV)

	xi, eta = lguys.to_tangent(rv_meas.ra, rv_meas.dec, ra0, dec0)

	rv_meas[!, :xi] = xi
	rv_meas[!, :eta] = eta

	r_ell = 60lguys.calc_r_ell(xi, eta, ell, PA)

	rv_meas[!, :r_ell] = r_ell

	rv_meas

end

# ╔═╡ fcf50b2c-5703-4101-b3a7-45e5fb8e072f
size(apogee, 1) + size(rv_meas, 1) # should be 617

# ╔═╡ 5cf336f6-e3eb-4668-b074-18b396f027be
prior_samples = DataFrame(
	sample(normal_dist(rv_meas.RV, rv_meas.RV_err), Prior(), 10000)
)

# ╔═╡ 74ad07df-f15c-436f-b390-ce95b27f7fab
let
	fig, ax = FigAxis(
		xgridvisible=false,
		ygridvisible=false,
		xlabel=xlabel,
		ylabel="density",
		title="priors"
	)
	
	plot_samples!(prior_samples, LinRange(80, 130, 100))

	fig
end

# ╔═╡ f9733c50-79c2-4fa5-8d9f-a259e553b6e9
hist(prior_samples.σ)

# ╔═╡ ccbb694d-f973-40fd-bab7-a2aefbd9fb0b
pairplot(prior_samples[:, [:μ, :σ]])

# ╔═╡ 87601bb8-4c90-4433-ab4c-275f3afc8d6a
fit_rv_sigma(rv_meas.RV, rv_meas.RV_err)

# ╔═╡ 318d29b9-4c84-4d38-af6d-048518952970
samples = DataFrame(sample(normal_dist(rv_meas.RV, rv_meas.RV_err), NUTS(0.65), 10000))

# ╔═╡ b18e4622-41e0-4700-9e4b-3dbebeefea53
describe(samples)

# ╔═╡ bc7bd936-3f62-4430-8acd-8331ca3ee5ad
pairplot(samples[:, [:μ, :σ]])

# ╔═╡ 0f5a9d9e-c5ca-4eb6-a0d2-5bb39b81daf6
σ_m = median(samples.σ)

# ╔═╡ 764b5306-20f9-4810-8188-1bdf9482260f
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = Arya.histogram(rv_meas.RV, 25, normalization=:pdf)
	
	plot_samples!(samples, LinRange(80, 140, 100))
	lines!(h)

	fig
end

# ╔═╡ 3a69f395-3c2d-4357-89af-5963d5fa79b8
let
	fig, ax = FigAxis(
		xlabel=L"$\log r_\textrm{ell}$ / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)
	
	scatter!(log10.(rv_meas.r_ell), rv_meas.RV)


	fig
end

# ╔═╡ c50f68d7-74c3-4c36-90c5-a5262982ed9f
μs, σs, μ_errs, σ_errs = binned_mu_sigma(rv_meas.r_ell, rv_meas.RV, rv_meas.RV_err, bins)

# ╔═╡ 33f54afc-cdb9-4eb8-887f-5a43281b837c
let
	fig, ax = FigAxis(
		xlabel = "R / arcmin",
		ylabel = L"RV / km s$^{-1}$"
	)

	scatter!(rv_meas.r_ell, rv_meas.RV)
	errscatter!(midpoints(bins), μs, yerr=μ_errs, color=:black)
	lines!(midpoints(bins), μs .+ σs, color=:black)
	lines!(midpoints(bins), μs .- σs, color=:black)

	fig
end

# ╔═╡ 1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
bin_errs = diff(bins) / 2

# ╔═╡ b7ebd916-bffc-4ffc-a7f7-44ee315e2528
stephist(rv_meas.r_ell, bins=bins)

# ╔═╡ 614f3f09-1880-491d-b41e-4e229330d66f
let
	fig, ax = FigAxis(
		xlabel = "R / arcmin",
		ylabel = L"$\sigma_{v, \textrm{los}}$ / km s$^{-1}$"
	)

	errscatter!(midpoints(bins), σs, yerr=σ_errs, xerr=bin_errs, color=:black)
	hlines!(σ_m)

	fig
end

# ╔═╡ 0fca9a64-347e-48f0-84fb-e7413e98d508
sev(rv_meas.RV)

# ╔═╡ e4294754-0931-4a0b-b9fe-0fa11878ea21
sev_iid(rv_meas.RV)

# ╔═╡ 0663e4c2-d68d-423c-9eed-1f47ab555621
rv_meas.RV_err

# ╔═╡ d688d2e5-faca-4b14-801a-d58b08fd6654
let
	fig, ax = FigAxis()

	filt = not.(ismissing.(rv_meas.RV))
	p = scatter!(rv_meas.ra[filt], rv_meas.dec[filt], color=rv_meas.RV[filt],
		colorrange=(90, 130),
		colormap=:bluesreds
	)

	#Colorbar(fig[1, 2], p)
	fig
end

# ╔═╡ 7178e5b9-cc42-4933-970a-4707ba69dbe9
let
	fig, ax = FigAxis()

	p = arrows!(rv_meas.ra, rv_meas.dec, rv_meas.pmra, rv_meas.pmdec, 				
		color=rv_meas.RV,
		colorrange=(90, 130),
		colormap=:bluesreds,
		lengthscale=0.1
	)

	#Colorbar(fig[1, 2], p)
	fig
end

# ╔═╡ 75442700-8532-4f86-8469-555c162edb97
Arya.std(rv_meas.RV)

# ╔═╡ 8b21cc49-ca17-4844-8238-e27e9752bee7
bins = Arya.bins_equal_number(rv_meas.r_ell, n=6)

# ╔═╡ 3741255c-d0f4-47c1-97dd-cc6823b547d9
# ╠═╡ disabled = true
#=╠═╡
bins = Arya.bins_equal_number(log10.(rv_meas.r_ell), n=10)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─6f2359d8-332b-11ef-0db9-f1f06474c561
# ╠═fcf50b2c-5703-4101-b3a7-45e5fb8e072f
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═9070c811-550c-4c49-9c58-0943b0f808b2
# ╠═e1cdc7ac-b1a4-45db-a363-2ea5b5ad9990
# ╠═dd29be70-4918-47e4-98cf-3608df14e88a
# ╠═23402e30-2064-494c-bd4a-8243cb474b61
# ╠═48c62811-136f-4962-a42c-b1dd1fc74f8c
# ╠═b7345279-4f80-47ad-a726-537571849eae
# ╠═a141ae2b-debd-453e-a7ea-d5a7287960a5
# ╠═2513e518-cb84-4fca-90d6-2f169e14f4f1
# ╠═1e06b622-ef1c-43ed-9453-bceab92c1889
# ╠═6734991c-16c0-4424-a2bb-84bfa811121f
# ╠═1b97d0e5-7a77-44a5-b609-ed8945cd959c
# ╠═2f1d6fad-7c99-4dae-822a-702d9a74bb94
# ╠═abbd2a53-e077-4af7-a168-b571e1a906b8
# ╠═74ad07df-f15c-436f-b390-ce95b27f7fab
# ╠═5cf336f6-e3eb-4668-b074-18b396f027be
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═9a99b3cb-90c0-4e5b-82c8-ae567ef6f7fa
# ╠═2474644e-06cc-4ee5-9318-dc7562de6b6c
# ╠═87601bb8-4c90-4433-ab4c-275f3afc8d6a
# ╠═318d29b9-4c84-4d38-af6d-048518952970
# ╠═d8f07585-ce89-4741-88d6-5bd3d5e4c71c
# ╠═55504311-fbbd-4351-b9e2-ecab053ba207
# ╠═764b5306-20f9-4810-8188-1bdf9482260f
# ╠═f9733c50-79c2-4fa5-8d9f-a259e553b6e9
# ╠═b18e4622-41e0-4700-9e4b-3dbebeefea53
# ╠═ccbb694d-f973-40fd-bab7-a2aefbd9fb0b
# ╠═bc7bd936-3f62-4430-8acd-8331ca3ee5ad
# ╠═7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═c2735c49-2892-46ac-bcf8-7cdcef409f44
# ╠═8b21cc49-ca17-4844-8238-e27e9752bee7
# ╠═c50f68d7-74c3-4c36-90c5-a5262982ed9f
# ╠═33f54afc-cdb9-4eb8-887f-5a43281b837c
# ╠═1eeb1572-4b97-4ccf-ad7a-dfd1e353bda7
# ╠═b7ebd916-bffc-4ffc-a7f7-44ee315e2528
# ╠═0f5a9d9e-c5ca-4eb6-a0d2-5bb39b81daf6
# ╠═614f3f09-1880-491d-b41e-4e229330d66f
# ╠═82a0e58a-30a4-4e42-b9c1-cb184eb551aa
# ╠═3a69f395-3c2d-4357-89af-5963d5fa79b8
# ╠═3741255c-d0f4-47c1-97dd-cc6823b547d9
# ╠═db78b22f-bee8-4d36-994e-6c7f9e96c9f2
# ╠═12ad8e67-0c7b-4029-91f6-5a2453ec799a
# ╠═e976d103-9956-47f3-adb0-7c449214b9d9
# ╠═2f5e5451-e653-43e8-8d2b-aeff4c575889
# ╠═f44d3c45-064c-4bf3-bf2c-e1f0d3504dd1
# ╠═0fca9a64-347e-48f0-84fb-e7413e98d508
# ╠═e4294754-0931-4a0b-b9fe-0fa11878ea21
# ╠═0663e4c2-d68d-423c-9eed-1f47ab555621
# ╠═a5e55f3e-de53-443c-9d28-e539aa6206a8
# ╠═9b4a0a1f-4b0c-4c90-b871-2bd244f0a908
# ╠═d688d2e5-faca-4b14-801a-d58b08fd6654
# ╠═7178e5b9-cc42-4933-970a-4707ba69dbe9
# ╠═75442700-8532-4f86-8469-555c162edb97
