### A Pluto.jl notebook ###
# v0.20.20

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

# ╔═╡ 7daee39d-f973-4a0b-aaf7-5e308e79b45c
using Turing

# ╔═╡ a662b1b4-5314-40d2-990f-28db38d21bf3
galaxyname = "sculptor"

# ╔═╡ 902876ac-2f5c-4fe3-9233-9ddace247d41
import TOML

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

# ╔═╡ ad2bed21-f45b-48fa-9dc9-e9bcff981d56
exp_fit = CSV.read(ENV["DWARFS_ROOT"] * "/observations/$galaxyname/mcmc/summary.mcmc_2exp.csv", DataFrame)

# ╔═╡ 70632fc8-71d9-4dc8-acde-381cb5fd9dde
function get_param(exp_fit, param)
	idx = only(findall(exp_fit.parameters .== [param]))
	return exp_fit.median[idx]
end

# ╔═╡ 4cbf1e33-7f7a-4d26-ae39-83cd28c52ff5
R_s_a = get_param(exp_fit, "R_s")

# ╔═╡ a3a48158-2930-46d2-815f-3855e62a709a
R_s_b = get_param(exp_fit, "R_s_outer")

# ╔═╡ af6f3c41-4c47-441a-94fc-63dcec90cc35
f_outer = get_param(exp_fit, "f_outer")

# ╔═╡ 080c0ea5-4ab5-4499-9909-a4338f9fdacb
minimum(stars.PSAT_RV)

# ╔═╡ 8483b37a-4751-427d-9e66-4ff2fd6ccfdc
hist(stars.vz)

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
	return fe_h[filt], fe_h_err[filt], stars.R_ell[filt]

end

# ╔═╡ ce4a09ec-8a88-48a0-bde2-8ceec99c7431
fe_h, fe_h_err, R_ell_fe = combine_fe_h(stars)

# ╔═╡ 65a50351-ae87-4675-91ff-ed23c0261c9b
hist(collect(skipmissing(fe_h)))

# ╔═╡ cf87677d-6ba1-4a0b-b51c-ff7dc3cccef5
scatter(stars.xi, stars.eta)

# ╔═╡ 25154575-f5ee-42b5-a35a-314770c158e7
hist(stars.R_ell)

# ╔═╡ fd5a1fd7-49a9-4cba-94d7-610a4333d1a4
md"""
# Model fit
"""

# ╔═╡ bb23fb9c-d8c5-4119-9223-f0aad311b9fc
@model function two_pop_model(v, v_err, R_ell, fe_h, fe_h_err, R_ell_fe)
	mu_fe_a ~ Normal(-2, 1)
	sigma_fe_a ~ Uniform(0, 3)
	mu_fe_b ~ Normal(-2, 1)
	sigma_fe_b ~ Uniform(0, 3)


	mu_vel_a ~ Normal(0, 100)
	sigma_vel_a ~ Uniform(0, 30)
	mu_vel_b ~ Normal(0, 100)
	sigma_vel_b ~ Uniform(0, 30)


	prof = LilGuys.Exp2D(R_s=R_s_a, M=1 - f_outer)
	prof_outer = LilGuys.Exp2D(R_s=R_s_b, M=f_outer)

	f_b = @. surface_density(prof_outer, R_ell) / (surface_density(prof, R_ell) + surface_density(prof_outer, R_ell))

	dist_vel_a = Normal(mu_vel_a, sigma_vel_a)
	dist_vel_b = Normal(mu_vel_b, sigma_vel_b)
	dist_fe_a = Normal(mu_fe_a, sigma_fe_a)
	dist_fe_b = Normal(mu_fe_b, sigma_fe_b)

	for i in eachindex(v)
		v[i] ~ MixtureModel([dist_vel_a, dist_vel_b], [1-f_b[i], f_b[i]])
	end


	f_b_fe = @. surface_density(prof_outer, R_ell_fe) / (surface_density(prof, R_ell_fe) + surface_density(prof_outer, R_ell_fe))

	for i in eachindex(fe_h)
		fe_h[i] ~ MixtureModel([dist_fe_a, dist_fe_b], [1-f_b_fe[i], f_b_fe[i]])
	end
	
end

# ╔═╡ 62182131-faae-4b8b-ac15-697895c1e3c2
model = two_pop_model(stars.vz, stars.vz_err, stars.R_ell, 
					  fe_h, fe_h_err, R_ell_fe)

# ╔═╡ 6ce1df6c-a64c-424c-83d5-d6e2443d8a98
samples = sample(model, NUTS(), 1000)

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

# ╔═╡ db418a63-7aaf-496f-ab08-272ffde18079
let
	fig=Figure()
	ax = Axis(fig[1,1],
			 xlabel = "μv (A)",
			 ylabel = "Δ μv (B - A)",
			 )

	scatter!(df.mu_vel_a, df.mu_vel_b .- df.mu_vel_a)
	hlines!(0, color=:black)

	fig
end

# ╔═╡ 0e42f462-eec5-4d56-b68d-ec5dc09eb868
let
	fig=Figure()
	ax = Axis(fig[1,1],
			 xlabel = "σv (A)",
			 ylabel = "Δ σv (B - A)",

			 )

	scatter!(df.sigma_vel_a, df.sigma_vel_b .- df.sigma_vel_a)
	hlines!(0, color=:black)
	fig
end

# ╔═╡ Cell order:
# ╠═a662b1b4-5314-40d2-990f-28db38d21bf3
# ╠═7d781874-b684-11f0-a926-e3f86a82519c
# ╠═77cc624b-77e7-4f7b-b72e-7309a2be8dad
# ╠═902876ac-2f5c-4fe3-9233-9ddace247d41
# ╠═5fabb947-3f85-4b81-b92d-43f58ae55973
# ╠═e070940a-d1e5-45f6-8737-a61431f29d58
# ╠═b92892ce-3bf8-4e60-95da-9776e59f75cb
# ╠═9c369fe4-d8c6-4890-95ec-fbc43eb80769
# ╠═e49b4ca6-1376-4cc5-a43d-859b4db6482d
# ╠═06ec245d-129c-4187-aec2-875da4faa2ca
# ╠═ad2bed21-f45b-48fa-9dc9-e9bcff981d56
# ╠═70632fc8-71d9-4dc8-acde-381cb5fd9dde
# ╠═4cbf1e33-7f7a-4d26-ae39-83cd28c52ff5
# ╠═a3a48158-2930-46d2-815f-3855e62a709a
# ╠═af6f3c41-4c47-441a-94fc-63dcec90cc35
# ╠═080c0ea5-4ab5-4499-9909-a4338f9fdacb
# ╠═8483b37a-4751-427d-9e66-4ff2fd6ccfdc
# ╠═ec5a8eb5-3c3f-4a81-b3b3-9f4aca18bd69
# ╠═ce4a09ec-8a88-48a0-bde2-8ceec99c7431
# ╠═65a50351-ae87-4675-91ff-ed23c0261c9b
# ╠═cf87677d-6ba1-4a0b-b51c-ff7dc3cccef5
# ╠═25154575-f5ee-42b5-a35a-314770c158e7
# ╟─fd5a1fd7-49a9-4cba-94d7-610a4333d1a4
# ╠═7daee39d-f973-4a0b-aaf7-5e308e79b45c
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
# ╠═db418a63-7aaf-496f-ab08-272ffde18079
# ╠═0e42f462-eec5-4d56-b68d-ec5dc09eb868
