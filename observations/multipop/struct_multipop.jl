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

# ╔═╡ 05ae7aa5-737f-401f-a844-1582bc1657ee
module MCMCUtils
	include("../mcmc/mcmc_utils.jl")
end

# ╔═╡ b92892ce-3bf8-4e60-95da-9776e59f75cb
md"""
# Data loading
"""

# ╔═╡ ca90b810-452e-49dc-9ab4-18878069d250
obs_props = MCMCUtils.get_obs_props(galaxyname)

# ╔═╡ 6533ff08-2fba-49a1-896e-d65880b7b69c
stars = MCMCUtils.get_fits(galaxyname, obs_props)

# ╔═╡ 98b6576b-b0ea-4e31-90fa-7b46c70dbbd9
R_max = maximum(stars.R_ell)

# ╔═╡ 299a6506-4d89-4f7a-9de8-46eb13181959


# ╔═╡ 9c369fe4-d8c6-4890-95ec-fbc43eb80769
obs_props_scl = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations/$galaxyname/observed_properties.toml"))

# ╔═╡ e49b4ca6-1376-4cc5-a43d-859b4db6482d
starsfile = Dict(
	"sculptor" => "rv_combined_x_wide_2c_psat_0.2.fits",
	"ursa_minor" => "rv_combined_x_2c_psat_0.2.fits",
)[galaxyname]

# ╔═╡ 06ec245d-129c-4187-aec2-875da4faa2ca
rv_stars = read_fits(ENV["DWARFS_ROOT"] * "/observations/$galaxyname/velocities/processed/" * starsfile)

# ╔═╡ 70632fc8-71d9-4dc8-acde-381cb5fd9dde
function get_param(exp_fit, param)
	idx = only(findall(exp_fit.parameters .== [param]))
	return exp_fit.median[idx]
end

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
fe_h, fe_h_err, R_ell_fe = combine_fe_h(rv_stars)

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
@model function two_pop_model(data::MCMCUtils.GaiaData, 
  		v::Vector{Float64}, v_err::Vector{Float64}, R_ell::Vector{Float64}, 
  		fe_h::Vector{Float64}, fe_h_err::Vector{Float64}, R_ell_fe::Vector{Float64})
	
	# structural 
	d_xi = 0# ~ Normal(0, 5)
	d_eta = 0# ~ Normal(0, 5)
	ellipticity ~ Uniform(0, 0.99)
	position_angle ~ Uniform(0, 180)
	R_s ~ LogUniform(0.1, 1e3)

	f_outer ~ Uniform(0, 1)

	R_s_outer ~ LogUniform(R_s * 1.5, 1.5e3)
	ellipticity_outer ~ Uniform(0, 0.99)
	position_angle_outer ~ Uniform(0, 180)


	prof = LilGuys.Exp2D(R_s=R_s, M=1 - f_outer)
	prof_outer = LilGuys.Exp2D(R_s=R_s_outer, M=f_outer)


	
	R = LilGuys.calc_R_ell(data.xi .+ d_xi, data.eta .+ d_eta, ellipticity, position_angle)
	R_outer = LilGuys.calc_R_ell(data.xi .+ d_xi, data.eta .+ d_eta, ellipticity_outer, position_angle_outer)

	
	L_sat_space = @. LilGuys.surface_density.(prof, R) + LilGuys.surface_density(prof_outer, R_outer)
	L_bg_space = 1 / (π*R_max^2)

	f_sat ~ Uniform(0, 1)
	LL = sum(@. log10.(
		(1-f_sat) * data.L_bg * L_bg_space
		+ f_sat * data.L_sat * L_sat_space
	))
	
	Turing.@addlogprob!(LL)

	#

	mu_fe_a ~ Normal(-2, 1)
	sigma_fe_a ~ Uniform(0, 3)
	mu_fe_b ~ Normal(-2, 1)
	sigma_fe_b ~ Uniform(0, 3)


	# mu_vel_a ~ Normal(0, 100)
	# sigma_vel_a ~ Uniform(0, 30)
	# mu_vel_b ~ Normal(0, 100)
	# sigma_vel_b ~ Uniform(0, 30)

	# f_b = @. surface_density(prof_outer, R_ell) / (surface_density(prof, R_ell) + surface_density(prof_outer, R_ell))

	# dist_vel_a = Normal(mu_vel_a, sigma_vel_a)
	# dist_vel_b = Normal(mu_vel_b, sigma_vel_b)
	dist_fe_a = Normal(mu_fe_a, sigma_fe_a)
	dist_fe_b = Normal(mu_fe_b, sigma_fe_b)

	# for i in eachindex(v)
	# 	v[i] ~ MixtureModel([dist_vel_a, dist_vel_b], [1-f_b[i], f_b[i]])
	# end


	f_b_fe = @. surface_density(prof_outer, R_ell_fe) / (surface_density(prof, R_ell_fe) + surface_density(prof_outer, R_ell_fe))

	for i in eachindex(fe_h)
		fe_h[i] ~ MixtureModel([dist_fe_a, dist_fe_b], [1-f_b_fe[i], f_b_fe[i]])
	end
	
end

# ╔═╡ 658da40d-355b-4efe-a6ea-8b6553398165
x = rand(10)

# ╔═╡ f0a93bb7-87c5-488b-ac7f-37f214f165e6
MixtureModel([filldist(Normal(0, 1), 10), filldist(Normal(2, 1), 10)], [x, 1 .-x])

# ╔═╡ b1428f5d-cc6e-463b-a056-7f7082ba82c8
@model function struct_model(data::MCMCUtils.GaiaData)
	
	# structural 
	d_xi = 0# ~ Normal(0, 5)
	d_eta = 0# ~ Normal(0, 5)
	ellipticity ~ Uniform(0, 0.99)
	position_angle ~ Uniform(0, 180)
	R_s ~ LogUniform(0.1, 1e3)

	f_outer ~ Uniform(0, 1)

	R_s_outer ~ LogUniform(R_s * 1.5, 2e3)
	ellipticity_outer ~ Uniform(0, 0.99)
	position_angle_outer ~ Uniform(0, 180)


	prof = LilGuys.Exp2D(R_s=R_s, M=1 - f_outer)
	prof_outer = LilGuys.Exp2D(R_s=R_s_outer, M=f_outer)


	
	R = LilGuys.calc_R_ell(data.xi .+ d_xi, data.eta .+ d_eta, ellipticity, position_angle)
	R_outer = LilGuys.calc_R_ell(data.xi .+ d_xi, data.eta .+ d_eta, ellipticity_outer, position_angle_outer)

	
	L_sat_space = @. LilGuys.surface_density.(prof, R) + LilGuys.surface_density(prof_outer, R_outer)
	L_bg_space = 1 / (π*R_max^2)

	f_sat ~ Uniform(0, 1)
	LL = sum(@. log10.(
		(1-f_sat) * data.L_bg * L_bg_space
		+ f_sat * data.L_sat * L_sat_space
	))
	
	Turing.@addlogprob!(LL)
end

# ╔═╡ 13310eae-a7e6-4994-b698-91d976e9b53a
gaia_data = MCMCUtils.GaiaData(stars)

# ╔═╡ 62182131-faae-4b8b-ac15-697895c1e3c2
model = two_pop_model(gaia_data,
					  rv_stars.vz, rv_stars.vz_err, rv_stars.R_ell, 
					  disallowmissing(fe_h), fe_h_err, R_ell_fe)

# ╔═╡ e4d79a2e-bdca-4121-8dee-35d9ec14a596
samples = sample(model, NUTS(), 100)

# ╔═╡ 63159af8-360a-432b-96dd-c0c4f9ee9f4c
model_st = struct_model(gaia_data)

# ╔═╡ 6ce1df6c-a64c-424c-83d5-d6e2443d8a98
samples_st = sample(model_st, NUTS(), 100)

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
# ╠═05ae7aa5-737f-401f-a844-1582bc1657ee
# ╠═b92892ce-3bf8-4e60-95da-9776e59f75cb
# ╠═ca90b810-452e-49dc-9ab4-18878069d250
# ╠═6533ff08-2fba-49a1-896e-d65880b7b69c
# ╠═98b6576b-b0ea-4e31-90fa-7b46c70dbbd9
# ╠═299a6506-4d89-4f7a-9de8-46eb13181959
# ╠═9c369fe4-d8c6-4890-95ec-fbc43eb80769
# ╠═e49b4ca6-1376-4cc5-a43d-859b4db6482d
# ╠═06ec245d-129c-4187-aec2-875da4faa2ca
# ╠═70632fc8-71d9-4dc8-acde-381cb5fd9dde
# ╠═ec5a8eb5-3c3f-4a81-b3b3-9f4aca18bd69
# ╠═ce4a09ec-8a88-48a0-bde2-8ceec99c7431
# ╠═65a50351-ae87-4675-91ff-ed23c0261c9b
# ╠═cf87677d-6ba1-4a0b-b51c-ff7dc3cccef5
# ╠═25154575-f5ee-42b5-a35a-314770c158e7
# ╟─fd5a1fd7-49a9-4cba-94d7-610a4333d1a4
# ╠═7daee39d-f973-4a0b-aaf7-5e308e79b45c
# ╠═bb23fb9c-d8c5-4119-9223-f0aad311b9fc
# ╠═f0a93bb7-87c5-488b-ac7f-37f214f165e6
# ╠═658da40d-355b-4efe-a6ea-8b6553398165
# ╠═b1428f5d-cc6e-463b-a056-7f7082ba82c8
# ╠═13310eae-a7e6-4994-b698-91d976e9b53a
# ╠═62182131-faae-4b8b-ac15-697895c1e3c2
# ╠═e4d79a2e-bdca-4121-8dee-35d9ec14a596
# ╠═63159af8-360a-432b-96dd-c0c4f9ee9f4c
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
