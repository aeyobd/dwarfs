### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ e9d818c0-1897-11ef-1913-1b38c9d4b4d6
begin 
	import Pkg; Pkg.activate()
	
	using CairoMakie
	
	import LilGuys as lguys
	using Arya
end

# ╔═╡ 49dfa506-3ab2-44ad-949e-2ce709d393c4
begin 
	profiles = Dict(
		"0.03" => lguys.ObsProfile("exp2d_rs0.03_stars_today_profile.toml"),
		"0.05" => lguys.ObsProfile("exp2d_rs0.05_stars_today_profile.toml"),
		"0.1" => lguys.ObsProfile("exp2d_rs0.1_stars_today_profile.toml"),
		"0.2" => lguys.ObsProfile("exp2d_rs0.2_stars_today_profile.toml"),
	) 
end

# ╔═╡ 787131ca-5d40-43f9-b48e-e8d19195f32a
begin 
	profiles_i = Dict(
		"0.03" => lguys.ObsProfile("exp2d_rs0.03_stars_i_today_profile.toml"),

		"0.05" => lguys.ObsProfile("exp2d_rs0.05_stars_i_today_profile.toml"),
		"0.1" => lguys.ObsProfile("exp2d_rs0.1_stars_i_today_profile.toml"),
		"0.2" => lguys.ObsProfile("exp2d_rs0.2_stars_i_today_profile.toml"),
		"0.3" => lguys.ObsProfile("exp2d_rs0.3_stars_i_today_profile.toml"),
		"0.5" => lguys.ObsProfile("exp2d_rs0.5_stars_i_today_profile.toml"),
	)
end

# ╔═╡ 91ba945f-3b51-4e67-ba96-922b3a880d3c
ks = ["0.03", "0.05", "0.1", "0.2"]

# ╔═╡ 6d207751-695f-4247-94da-ced5a146092f
prof_expected = lguys.ObsProfile("/cosma/home/durham/dc-boye1/sculptor/fiducial_sample_profile.toml")

# ╔═╡ de69d265-cda1-4543-ad02-2ee3091964d6
log_r_label = "log r / arcmin"

# ╔═╡ 932c4fef-992b-4518-80d0-59c8e126ccb5
let 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=log_r_label,
		ylabel = L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
		limits=((-1.5, 3), (-10, 0))
	)

	for k in ks
		profile = profiles[k]

		label = "$k"
		lines!(profile.log_r, profile.log_Sigma, 
			label=label)
	end
	errscatter!(prof_expected.log_r, prof_expected.log_Sigma,
		yerr=prof_expected.log_Sigma_err,
		label="J+24"
	)
	
	
	vlines!(log10(295))
	axislegend("Rs / kpc")

	fig
end

# ╔═╡ c843b157-c95c-442d-8166-503bf0484fee
let 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=log_r_label,
		ylabel = L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
		limits=((-1.5, 3), (-10, 0))
	)

	for k in ks
		profile = profiles_i[k]

		label = "$k"
		lines!(profile.log_r, profile.log_Sigma, 
			label=label)
	end
	
	errscatter!(prof_expected.log_r, prof_expected.log_Sigma,
		yerr=prof_expected.log_Sigma_err,
		label="J+24"
	)
	
	
	vlines!(log10(295))
	axislegend(ax)

	fig
end

# ╔═╡ ddf6a3f5-a6ca-4274-89fa-76ff940c1abd
arcmin_to_rad = 60/206265

# ╔═╡ d3e934c2-cc73-4873-8157-8303161ceecf
distance=0.3*86

# ╔═╡ e0db88b5-b601-48d3-8c59-a8158736acfb
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel=log_r_label,
		ylabel = L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
		limits=((-1.5, 3), (-10, 0))
	)

	prof = profiles_i["0.1"]
	
	errscatter!(prof.log_r, prof.log_Sigma,
		yerr=prof.log_Sigma_err,
	)


	R_s_kpc = 0.1

	R_s = R_s_kpc / ( distance*arcmin_to_rad )
	prof_a = lguys.Exp2D(R_s=R_s)
	println(R_s)

	x_model = LinRange(-0.5, 3, 1000)
	r = 10 .^ x_model
	Σ(r) = lguys.calc_Σ(prof_a, r)
	y_model = log10.(Σ.(r))

	lines!(x_model, y_model)


	fig
end

# ╔═╡ eebcdc5e-a755-45fe-a740-d2c77d647fa7
lguys.calc_Σ_from_ρ(lguys.Exp2D(R_s=0.1), 0.3)

# ╔═╡ b23c5630-fb8d-4132-a29d-d2d408e247ab
lguys.calc_Σ(lguys.Exp2D(R_s=0.1), 0.3)

# ╔═╡ 5d03ccd2-2dd6-45d9-8ec1-1af9bea3475b
md"""
# Scratch
"""

# ╔═╡ Cell order:
# ╠═e9d818c0-1897-11ef-1913-1b38c9d4b4d6
# ╠═49dfa506-3ab2-44ad-949e-2ce709d393c4
# ╠═787131ca-5d40-43f9-b48e-e8d19195f32a
# ╠═91ba945f-3b51-4e67-ba96-922b3a880d3c
# ╠═6d207751-695f-4247-94da-ced5a146092f
# ╠═de69d265-cda1-4543-ad02-2ee3091964d6
# ╠═932c4fef-992b-4518-80d0-59c8e126ccb5
# ╠═c843b157-c95c-442d-8166-503bf0484fee
# ╠═ddf6a3f5-a6ca-4274-89fa-76ff940c1abd
# ╠═d3e934c2-cc73-4873-8157-8303161ceecf
# ╠═e0db88b5-b601-48d3-8c59-a8158736acfb
# ╠═eebcdc5e-a755-45fe-a740-d2c77d647fa7
# ╠═b23c5630-fb8d-4132-a29d-d2d408e247ab
# ╟─5d03ccd2-2dd6-45d9-8ec1-1af9bea3475b
