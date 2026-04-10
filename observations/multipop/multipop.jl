### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 7d781874-b684-11f0-a926-e3f86a82519c
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ 77cc624b-77e7-4f7b-b72e-7309a2be8dad
using PyFITS

# ╔═╡ 5fabb947-3f85-4b81-b92d-43f58ae55973
using PairPlots

# ╔═╡ e070940a-d1e5-45f6-8737-a61431f29d58
using CSV, DataFrames

# ╔═╡ fc38836b-9f74-4ba4-81b6-5511a8bb7534
using KernelDensity

# ╔═╡ 7daee39d-f973-4a0b-aaf7-5e308e79b45c
using Turing

# ╔═╡ a662b1b4-5314-40d2-990f-28db38d21bf3
galaxyname = "sculptor"

# ╔═╡ 902876ac-2f5c-4fe3-9233-9ddace247d41
import TOML

# ╔═╡ 25e3c178-8771-427b-97c7-d822d594ae6c
import Random

# ╔═╡ 436f47c7-c161-4e85-ba87-462377b1705c
import Interpolations

# ╔═╡ f5392cbf-93a9-404b-bd4d-a26b080de5c6
md"""
# KDE distributions
"""

# ╔═╡ 0c5b69c6-9656-4e2d-ac60-5562dd2fb727
begin 
	struct InterpKDEDistribution{T<:Real,K<:KernelDensity.InterpKDE} <: ContinuousUnivariateDistribution
	    kde::K
	end
	function InterpKDEDistribution(k::KernelDensity.InterpKDE)
	    T = eltype(k.kde.x)
	    return InterpKDEDistribution{T,typeof(k)}(k)
	end
	function InterpKDEDistribution(k::KernelDensity.UnivariateKDE)
	    return InterpKDEDistribution(KernelDensity.InterpKDE(k))
	end
	
	function Distributions.minimum(d::InterpKDEDistribution)
	    return first(only(Interpolations.bounds(d.kde.itp.itp)))
	end
	
	function Distributions.maximum(d::InterpKDEDistribution)
	    return last(only(Interpolations.bounds(d.kde.itp.itp)))
	end
	
	function Distributions.pdf(d::InterpKDEDistribution, x::Real)
	    return pdf(d.kde, x)
	end
	
	function Distributions.logpdf(d::InterpKDEDistribution, x::Real)
	    return log(pdf(d, x))
	end
	
	# need to at least have a very rough implementation of this
	# much better/more efficient implementation is possible,
	# but this will only be called once when sampling starts 
	function Random.rand(rng::Random.AbstractRNG, d::InterpKDEDistribution)
	    (; kde) = d
	    knots = Interpolations.knots(kde.itp.itp).knots
	    cdf = cumsum(pdf.(Ref(kde), knots) .* LilGuys.gradient(knots))
	    u = rand(rng)
		if u >= maximum(cdf)
			return knots[end]
		elseif u <= minimum(cdf)
			return knots[1]
		end
	    return knots[findlast(u .> cdf)]
	end

end

# ╔═╡ d32e12a3-1bd6-42d7-8614-68f5ad4f69c9
function get_kde(samples, param)
	x = samples[!, param]
	k = kde(x, boundary=extrema(x))
	@info InterpKDE(k).kde.x
	@info Interpolations.bounds(InterpKDE(k).itp.itp)
	return InterpKDEDistribution(k)
end

# ╔═╡ b92892ce-3bf8-4e60-95da-9776e59f75cb
md"""
# Data loading
"""

# ╔═╡ 9c369fe4-d8c6-4890-95ec-fbc43eb80769
obs_props_scl = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxyname/observed_properties.toml"))

# ╔═╡ e49b4ca6-1376-4cc5-a43d-859b4db6482d
starsfile = Dict(
	"sculptor" => "rv_combined_x_wide_2c_psat_0.2.fits",
	"ursa_minor" => "rv_combined_x_2c_psat_0.2.fits",
)[galaxyname]

# ╔═╡ 06ec245d-129c-4187-aec2-875da4faa2ca
stars = read_fits(ENV["DWARFS_ROOT"] * "/observations/$galaxyname/velocities/processed/" * starsfile)

# ╔═╡ c0b8d74a-fd7b-467e-b891-c5cb28f6ed37
exp_samples = CSV.read(ENV["DWARFS_ROOT"] * "/observations/$galaxyname/mcmc/samples.mcmc_2exp.csv", DataFrame)

# ╔═╡ ec5a8eb5-3c3f-4a81-b3b3-9f4aca18bd69
function combine_fe_h(stars)
	cols = filter(x->startswith(x, "fe_h") && !contains(x, "_err"),
				 collect(names(stars)))

	Nr = size(stars, 1)
	fe_h = Vector{Union{Missing, Float64}}(undef, Nr)
	fe_h_err = zeros(Nr)
	filt = falses(Nr)

	for i in 1:Nr
		for col in cols
			if !ismissing(stars[i, col])
				fe_h[i] = stars[i, col]
				filt[i] = true

				# study =
				# fe_h_err[i] = stars[i, "fe_h"]
				break
			end
		end

	end

	fe_h[.!filt] .= missing
	return fe_h[filt], fe_h_err[filt], stars.xi[filt], stars.eta[filt]

end

# ╔═╡ ce4a09ec-8a88-48a0-bde2-8ceec99c7431
fe_h, fe_h_err, xi_fe, eta_fe = combine_fe_h(stars)

# ╔═╡ fd5a1fd7-49a9-4cba-94d7-610a4333d1a4
md"""
# Model fit
"""

# ╔═╡ bb23fb9c-d8c5-4119-9223-f0aad311b9fc
@model function two_pop_model( fe_h, fe_h_err, xi_fe, eta_fe;
							 prior_R_s = get_kde(exp_samples, "R_s"),
							 prior_R_s_outer = get_kde(exp_samples, "R_s_outer"),
							 prior_f_outer = get_kde(exp_samples, "f_outer"),
							 prior_ell = get_kde(exp_samples, "ellipticity"),
							 prior_ell_outer = get_kde(exp_samples, "ellipticity_outer"),
							 prior_PA = get_kde(exp_samples, "position_angle"),
							 prior_PA_outer = get_kde(exp_samples, "position_angle_outer"),
							 )
	mu_fe_a ~ Normal(-2, 1)
	sigma_fe_a ~ Uniform(0, 3)
	mu_fe_b ~ Normal(-2, 1)
	sigma_fe_b ~ Uniform(0, 3)


	f_outer ~ prior_f_outer
	R_s ~ prior_R_s
	R_s_outer ~ prior_R_s_outer
	ell ~ prior_ell
	ell_outer ~ prior_ell_outer
	PA ~ prior_PA
	PA_outer ~ prior_PA_outer

	R_ell = LilGuys.calc_R_ell.(xi_fe, eta_fe, ell, PA)
	R_ell_outer = LilGuys.calc_R_ell.(xi_fe, eta_fe, ell_outer, PA_outer)
	@assert all(R_ell .> 0) && all(R_ell_outer .> 0)

	
	prof = LilGuys.Exp2D(R_s=R_s, M=1 - f_outer)
	prof_outer = LilGuys.Exp2D(R_s=R_s_outer, M=f_outer)

	Σ_out = surface_density.(prof_outer, R_ell_outer)
	Σ_in = surface_density.(prof, R_ell)
	
	f_in = @.  Σ_in/(Σ_in + Σ_out)

	dist_fe_a = Normal(mu_fe_a, sigma_fe_a)
	dist_fe_b = Normal(mu_fe_b, sigma_fe_b)

	for i in eachindex(fe_h)
		fe_h[i] ~ MixtureModel([dist_fe_a, dist_fe_b], [f_in[i], 1 - f_in[i]])
	end
	
end

# ╔═╡ 62182131-faae-4b8b-ac15-697895c1e3c2
model = two_pop_model(fe_h, fe_h_err, xi_fe, eta_fe)

# ╔═╡ 6ce1df6c-a64c-424c-83d5-d6e2443d8a98
samples = sample(model, NUTS(), MCMCThreads(), 1000, 2)

# ╔═╡ c33a2b82-596c-423a-b4b0-abf2f3df941d
pairplot(samples)

# ╔═╡ 6d4b2d04-7083-44ab-adc0-4cefaab8c6f4
PVALUE = 0.16

# ╔═╡ 041f2040-dfce-4b8c-93ef-d0709d0dbfd8
"""
    summarize(samples)

Sumerizes a Turing chain. This loosely wraps Turings version
except adds median and quantile uncertainties.
"""
function summarize(samples)
	chain_summary = DataFrame(Turing.summarize(samples))
	Nvar = size(chain_summary, 1)

	medians = zeros(Nvar)
	err_lows = zeros(Nvar)
	err_highs = zeros(Nvar)
	for i in 1:Nvar
        l, m, h = quantile(samples[:, i, :], [PVALUE, 0.5, 1-PVALUE])
		medians[i] = m
		err_lows[i] = m - l
		err_highs[i] = h - m
	end

	chain_summary[!, :parameters] = string.(chain_summary.parameters)
	chain_summary[!, :median] = medians
	chain_summary[!, :lower_error] = err_lows
	chain_summary[!, :upper_error] = err_highs

	chain_summary
end

# ╔═╡ e0544f33-76ef-4d43-83bc-914ad3c37721
df = DataFrame(samples)

# ╔═╡ a3f7bbaa-13c2-4caf-9636-59fd4e311c43
summary = summarize(samples)

# ╔═╡ 9f6f33bf-5a87-4c39-832e-eb712e8df504
filename_out = "processed/$galaxyname.mcmc_2pop_vel_fe"

# ╔═╡ 7246570d-a920-4b7d-8e80-453310c27881
CSV.write(filename_out * ".summary.csv", summary)

# ╔═╡ 7c87c114-3e83-4382-8975-1eab0205329f
write_fits(filename_out * ".samples.fits", df, overwrite=true)

# ╔═╡ bd7dc274-a1aa-49e0-a517-160c63c9de63
md"""
# Additional plots
"""

# ╔═╡ 965cbbd8-f79c-48a0-a168-1a383d9d31dd
let
	fig=Figure()
	ax = Axis(fig[1,1],
			 xlabel = "[Fe/H] (A)",
			 ylabel = "Δ [Fe/H] (B - A)",
			 )

	scatter!(df.mu_fe_a, df.mu_fe_b .- df.mu_fe_a)
	hlines!(0, color=:black)

	fig
end

# ╔═╡ 250674e8-f676-40a2-b66a-59d8e5392e1c
let
	var = :f_outer
	fig = Figure()
	dist = get_kde(exp_samples, var)

	ax = Axis(fig[1,1])

	lines!(dist)
	stephist!(exp_samples[!, var], bins=100, normalization=:pdf)

	x = rand(dist, 10_000)
	stephist!(x, bins=100, normalization=:pdf)
	fig
end

# ╔═╡ Cell order:
# ╠═a662b1b4-5314-40d2-990f-28db38d21bf3
# ╠═7d781874-b684-11f0-a926-e3f86a82519c
# ╠═77cc624b-77e7-4f7b-b72e-7309a2be8dad
# ╠═902876ac-2f5c-4fe3-9233-9ddace247d41
# ╠═5fabb947-3f85-4b81-b92d-43f58ae55973
# ╠═e070940a-d1e5-45f6-8737-a61431f29d58
# ╠═fc38836b-9f74-4ba4-81b6-5511a8bb7534
# ╠═25e3c178-8771-427b-97c7-d822d594ae6c
# ╠═436f47c7-c161-4e85-ba87-462377b1705c
# ╠═7daee39d-f973-4a0b-aaf7-5e308e79b45c
# ╟─f5392cbf-93a9-404b-bd4d-a26b080de5c6
# ╠═0c5b69c6-9656-4e2d-ac60-5562dd2fb727
# ╠═d32e12a3-1bd6-42d7-8614-68f5ad4f69c9
# ╟─b92892ce-3bf8-4e60-95da-9776e59f75cb
# ╠═9c369fe4-d8c6-4890-95ec-fbc43eb80769
# ╠═e49b4ca6-1376-4cc5-a43d-859b4db6482d
# ╠═06ec245d-129c-4187-aec2-875da4faa2ca
# ╠═c0b8d74a-fd7b-467e-b891-c5cb28f6ed37
# ╠═ec5a8eb5-3c3f-4a81-b3b3-9f4aca18bd69
# ╠═ce4a09ec-8a88-48a0-bde2-8ceec99c7431
# ╟─fd5a1fd7-49a9-4cba-94d7-610a4333d1a4
# ╠═bb23fb9c-d8c5-4119-9223-f0aad311b9fc
# ╠═62182131-faae-4b8b-ac15-697895c1e3c2
# ╠═6ce1df6c-a64c-424c-83d5-d6e2443d8a98
# ╠═c33a2b82-596c-423a-b4b0-abf2f3df941d
# ╠═6d4b2d04-7083-44ab-adc0-4cefaab8c6f4
# ╠═041f2040-dfce-4b8c-93ef-d0709d0dbfd8
# ╠═e0544f33-76ef-4d43-83bc-914ad3c37721
# ╠═a3f7bbaa-13c2-4caf-9636-59fd4e311c43
# ╠═9f6f33bf-5a87-4c39-832e-eb712e8df504
# ╠═7246570d-a920-4b7d-8e80-453310c27881
# ╠═7c87c114-3e83-4382-8975-1eab0205329f
# ╟─bd7dc274-a1aa-49e0-a517-160c63c9de63
# ╠═965cbbd8-f79c-48a0-a168-1a383d9d31dd
# ╠═250674e8-f676-40a2-b66a-59d8e5392e1c
