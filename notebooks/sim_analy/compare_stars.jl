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
		"0.05" => lguys.ObsProfile("exp2d_rs0.05_today_profile.toml"),
		"0.1" => lguys.ObsProfile("exp2d_rs0.1_today_profile.toml"),
		"0.2" => lguys.ObsProfile("exp2d_rs0.2_today_profile.toml"),
	) 
end

# ╔═╡ 787131ca-5d40-43f9-b48e-e8d19195f32a
begin 
	profiles_i = Dict(
		"0.05" => lguys.ObsProfile("exp2d_rs0.05_i_today_profile.toml"),
		"0.1" => lguys.ObsProfile("exp2d_rs0.1_i_today_profile.toml"),
		"0.2" => lguys.ObsProfile("exp2d_rs0.2_i_today_profile.toml"),
	)
end

# ╔═╡ 91ba945f-3b51-4e67-ba96-922b3a880d3c
ks = ["0.05", "0.1", "0.2"]

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

	errscatter!(prof_expected.log_r, prof_expected.log_Sigma,
		yerr=prof_expected.log_Sigma_err,
		label="J+24",
		color=:black
	)

	
	for k in ks
		profile = profiles[k]

		label = "$k"
		lines!(profile.log_r, profile.log_Sigma, 
			label=label)
	end
	
	for i in 1:length(ks)
		k = ks[i]
		profile = profiles_i[k]

		label = "$k"
		lines!(profile.log_r, profile.log_Sigma,
			linestyle=:dash, color=COLORS[i]
		)
	end

	
	vlines!(log10(295))
	axislegend("Rs / kpc")

	fig
end

# ╔═╡ 88f31bfc-a67d-4654-af8b-46dc91500558
r_b = 296

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
	
	
	vlines!(log10(r_b))
	axislegend(ax)

	fig
end

# ╔═╡ 09207e06-3ab5-41b2-a9af-e77d89b34f59
let
	profile = profiles["0.1"]
	fig = Figure(size=(700, 300))
	ax = Axis(fig[1, 1], limits=((-0.8, 2), (-12, 10)),
		xlabel=log_r_label,
		ylabel=L"\Gamma"
	)

	
	errscatter!(profile.log_r, profile.Gamma, yerr=profile.Gamma_err)
	
	vlines!(log10(r_b))

	ax_lin = Axis(fig[1, 2],
		xlabel="r / arcmin",
		yticklabelsvisible=false,
	)

	
	errscatter!(value.(10 .^ profile.log_r), profile.Gamma, yerr=profile.Gamma_err)
	vlines!((r_b))

	linkyaxes!(ax, ax_lin)

	fig
end

# ╔═╡ d3e934c2-cc73-4873-8157-8303161ceecf
distance=86

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

	R_s_arcmin = lguys.kpc_to_arcmin(R_s_kpc, distance)
	println(R_s_arcmin)

	profile2 = lguys.Exp2D(R_s = R_s_arcmin)

	log_Σ(r) = log10(lguys.calc_Σ(profile2, r))

	log_R = LinRange(-2, 2, 1000)
	y = log_Σ.(10 .^ log_R)
	
	lines!(log_R, y)


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
# ╠═88f31bfc-a67d-4654-af8b-46dc91500558
# ╠═09207e06-3ab5-41b2-a9af-e77d89b34f59
# ╠═d3e934c2-cc73-4873-8157-8303161ceecf
# ╠═e0db88b5-b601-48d3-8c59-a8158736acfb
# ╠═eebcdc5e-a755-45fe-a740-d2c77d647fa7
# ╠═b23c5630-fb8d-4132-a29d-d2d408e247ab
# ╟─5d03ccd2-2dd6-45d9-8ec1-1af9bea3475b