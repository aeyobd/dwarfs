### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ e9d818c0-1897-11ef-1913-1b38c9d4b4d6
begin 
	import Pkg; Pkg.activate()
	
	using CairoMakie
	
	import LilGuys as lguys
	using Arya
end

# ╔═╡ 0f7680b9-4529-4438-8a5e-a9d6e6530eaf
import TOML

# ╔═╡ 96ae7f95-886b-4fd6-a1a0-3d4f0790cc34
cd("/astro/dboyea/sculptor/orbits/1e6/orbit1//V70_r0.4/stars")

# ╔═╡ c64b9fe9-b3ff-498f-85dd-55e227443b63
name = "exp2d_rs0.07"

# ╔═╡ 6d207751-695f-4247-94da-ced5a146092f
prof_expected = lguys.ObsProfile("/arc7/home/dboyea/dwarfs/notebooks/density_fits/sculptor/fiducial_profile.toml")

# ╔═╡ 88f31bfc-a67d-4654-af8b-46dc91500558
r_b = 68

# ╔═╡ de26b439-d9dd-449e-9289-9d3daed87cc7
fig_dir = "figures"

# ╔═╡ b00d8ca9-f080-46a5-883d-1aa45a683e3f
mkpath(fig_dir)

# ╔═╡ 9e8f7fab-fdbf-4af3-afec-09b65da4d019
md"""
# loading data
"""

# ╔═╡ e032f8e2-28b5-46c4-896e-27da9dfff22f
params = TOML.parsefile("/astro/dboyea/sculptor/isolation/1e6/halos/V70_r0.4/stars/$name.toml")

# ╔═╡ 41d07a5b-32a7-43f2-9a34-eab54c8ab4e0
R_s = params["profile_kwargs"]["R_s"]

# ╔═╡ 3b3a4fc6-5bf3-4b81-8621-8f425c7dd517
R_s_arcmin = lguys.kpc_to_arcmin(R_s, 82.3)

# ╔═╡ f700b0bb-64b5-4be0-bf8a-88e72a7500e1
prof_f = lguys.ObsProfile("$(name)_today_profile.toml")

# ╔═╡ 2e81adf7-47f4-4f61-857e-8997c15dc943
prof_i = lguys.ObsProfile("$(name)_i_today_profile.toml")

# ╔═╡ 2e69a41e-1b32-43e4-a23c-4550539275e6
prof_a_sky = lguys.Exp2D(; R_s=R_s_arcmin)

# ╔═╡ 49dfa506-3ab2-44ad-949e-2ce709d393c4
begin 
	profiles = Dict(
		"initial" => prof_i,
		"final" => prof_f
	) 
end

# ╔═╡ 2be5fd09-3b4c-4749-b161-fdc6c4d10142
begin 
	ks = collect(keys(profiles))
end

# ╔═╡ 10792def-f5cc-4963-9822-88ff8eabed95
function plot_model_and_exp!(prof, R_s_kpc; 
		y_offset=0, distance=distance, kwargs...)
	R_s_arcmin = lguys.kpc_to_arcmin(R_s_kpc, distance)
	println(R_s_arcmin)

	profile2 = lguys.Exp2D(R_s = R_s_arcmin)

	log_Σ(r) = log10(lguys.calc_Σ(profile2, r))

	log_R = LinRange(-2, maximum(prof.log_r), 1000)
	y = log_Σ.(10 .^ log_R)
	y .-= y[1]
	
	lines!(log_R, y .+ y_offset; kwargs...)

	errscatter!(prof.log_r, prof.log_Sigma .+ y_offset,
		yerr=prof.log_Sigma_err,
	)
end

# ╔═╡ de69d265-cda1-4543-ad02-2ee3091964d6
log_r_label = "log R / arcmin"

# ╔═╡ bfab4ae8-a94b-4a81-aad3-9b706f2474bb
function sigma_axis(; kwargs...) 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=log_r_label,
		ylabel = L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}"
		;kwargs...
	)

	return fig, ax
end

# ╔═╡ 932c4fef-992b-4518-80d0-59c8e126ccb5
let 
	fig = Figure(
		backgroundcolor=:transparent
	)
	ax = Axis(fig[1,1], 
		xlabel=log_r_label,
		ylabel = L"\log\, \Sigma\,/\,\Sigma_0",
		limits=((-1, 2.1), (-4.5, 0.5)),
		backgroundcolor=(:black, 0.)
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
	

	
	vlines!(log10(r_b), color=:grey, label="break radius")
	axislegend(position=:lb,
		backgroundcolor=:transparent
	)

	save(joinpath(fig_dir, "density_i_f.pdf"), fig)

	fig
end

# ╔═╡ 693cb85b-8f19-414f-8a26-da2035232df0
abspath(fig_dir)

# ╔═╡ 09207e06-3ab5-41b2-a9af-e77d89b34f59
let
	fig = Figure(size=(700, 300))
	ax = Axis(fig[1, 1], limits=((-0.8, 2), (-12, 10)),
		xlabel=log_r_label,
		ylabel=L"\Gamma"
	)

	for (name, profile) in profiles
		errscatter!(profile.log_r, profile.Gamma, yerr=profile.Gamma_err)
	end

	
	vlines!(log10(r_b))

	ax_lin = Axis(fig[1, 2],
		xlabel="r / arcmin",
		yticklabelsvisible=false,
		limits=((0, 200), nothing)
	)
	
	for (name, profile) in profiles
		errscatter!(10 .^ profile.log_r, profile.Gamma, yerr=profile.Gamma_err)
	end
	vlines!((r_b))

	linkyaxes!(ax, ax_lin)

	fig
end

# ╔═╡ b23c5630-fb8d-4132-a29d-d2d408e247ab
md"""
## Initial conditions
"""

# ╔═╡ 67f4b666-184a-4276-9306-e85a72b399a3
let 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=log_r_label,
		ylabel = L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
		limits=((-1.5, 3), (-5, 1))
	)

	
	
	scatter!(prof_i.log_r, prof_i.log_Sigma, 
			label="initial conditions")


	log_r = LinRange(-1, 2, 1000)
	r = 10 .^ log_r
	sigma = lguys.calc_Σ.(prof_a_sky, r)
	sigma ./= sigma[1] / 1.5

	lines!(log_r, log10.(sigma))

	fig
end

# ╔═╡ 5d03ccd2-2dd6-45d9-8ec1-1af9bea3475b


# ╔═╡ Cell order:
# ╠═e9d818c0-1897-11ef-1913-1b38c9d4b4d6
# ╠═0f7680b9-4529-4438-8a5e-a9d6e6530eaf
# ╠═96ae7f95-886b-4fd6-a1a0-3d4f0790cc34
# ╠═c64b9fe9-b3ff-498f-85dd-55e227443b63
# ╠═6d207751-695f-4247-94da-ced5a146092f
# ╠═88f31bfc-a67d-4654-af8b-46dc91500558
# ╠═de26b439-d9dd-449e-9289-9d3daed87cc7
# ╠═b00d8ca9-f080-46a5-883d-1aa45a683e3f
# ╠═9e8f7fab-fdbf-4af3-afec-09b65da4d019
# ╠═e032f8e2-28b5-46c4-896e-27da9dfff22f
# ╠═41d07a5b-32a7-43f2-9a34-eab54c8ab4e0
# ╠═3b3a4fc6-5bf3-4b81-8621-8f425c7dd517
# ╠═f700b0bb-64b5-4be0-bf8a-88e72a7500e1
# ╠═2e81adf7-47f4-4f61-857e-8997c15dc943
# ╠═2e69a41e-1b32-43e4-a23c-4550539275e6
# ╠═49dfa506-3ab2-44ad-949e-2ce709d393c4
# ╠═2be5fd09-3b4c-4749-b161-fdc6c4d10142
# ╠═10792def-f5cc-4963-9822-88ff8eabed95
# ╠═bfab4ae8-a94b-4a81-aad3-9b706f2474bb
# ╠═de69d265-cda1-4543-ad02-2ee3091964d6
# ╠═932c4fef-992b-4518-80d0-59c8e126ccb5
# ╠═693cb85b-8f19-414f-8a26-da2035232df0
# ╠═09207e06-3ab5-41b2-a9af-e77d89b34f59
# ╟─b23c5630-fb8d-4132-a29d-d2d408e247ab
# ╠═67f4b666-184a-4276-9306-e85a72b399a3
# ╠═5d03ccd2-2dd6-45d9-8ec1-1af9bea3475b
