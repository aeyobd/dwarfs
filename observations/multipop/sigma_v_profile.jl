### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 8a84cd7a-b75e-11f0-ac6e-5386c69b3015
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using CairoMakie
	using Arya
end

# ╔═╡ 9850dd29-11ec-4872-90b3-7161e61a2cac
md"""
# Theoretical modeling
"""

# ╔═╡ 19866432-ee61-4891-88b8-6c2ab3f20c64
halo = NFW(r_circ_max = 3.2, v_circ_max=31/V2KMS)

# ╔═╡ bb15d2b6-313c-4f2c-877f-4a4aceb683f0
prof = LilGuys.Exp2D(R_s=0.10)

# ╔═╡ 0d6a1b04-7f20-42c3-870c-7f6e76ceea17
md"""
# Snapshot tyme
"""

# ╔═╡ 6775faf6-8bba-4285-9414-4c5c8263e6be
f_outer = 0.32

# ╔═╡ 7e12fb10-688f-4305-8922-a9f48db7680e
snap_double = let
	galaxyname, modelname, starsname = "sculptor", "1e6_new_v31_r3.2/orbit_smallperi", "exp2d_rs0.10"
	starsname_outer = "exp2d_rs0.22"

	haloname = joinpath(modelname, "..")
	stars_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, "stars")

	df_probs = LilGuys.read_hdf5_table(joinpath(stars_dir, starsname, "probabilities_stars.hdf5"))
	df_probs_outer = LilGuys.read_hdf5_table(joinpath(stars_dir, starsname_outer, "probabilities_stars.hdf5"))

	 
	snap_i = Snapshot(joinpath(stars_dir, "iso_paint.hdf5"))

	probs = @. df_probs.probability * (1-f_outer) + df_probs_outer.probability * (f_outer)
	LilGuys.add_stars!(snap_i, probs)


	snap_i
end

# ╔═╡ de504a7f-a151-490b-9535-4122a630a714
function load_snaps(galaxyname, modelname, starsname)
	haloname = joinpath(modelname, "..")
	stars_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, haloname, "stars")

	df_probs = LilGuys.read_hdf5_table(joinpath(stars_dir, starsname, "probabilities_stars.hdf5"))

	snap_i = Snapshot(joinpath(stars_dir, "iso_paint.hdf5"))
	LilGuys.add_stars!(snap_i, df_probs.probability)

	snap_f = Snapshot(joinpath(stars_dir, "iso_final.hdf5"))
	LilGuys.add_stars!(snap_f, df_probs.probability)

	return snap_i, snap_f
end

# ╔═╡ 90f7565d-6ee8-4ed8-bfcb-f2e97cfa9252
snap_i = load_snaps("sculptor", "1e6_new_v31_r3.2/orbit_smallperi", "exp2d_rs0.10")[1]

# ╔═╡ 39f091fe-c7c3-4427-8b8e-0dc7245db018
snap_plummer = load_snaps("sculptor", "1e6_new_v31_r3.2/orbit_smallperi", "plummer_rs0.20")[1]

# ╔═╡ aa68ef1a-0f77-4a2f-b6d2-eb6e8deb01c5
function calc_velocity_profile(snap)
	rs = LilGuys.radii(snap.positions[1:2, :])
	
	bins = LilGuys.bins_both(log10.(rs), snap.weights, bin_width=0.05, num_per_bin=2)

	N_bins = length(bins) - 1

	σv = zeros(N_bins)
	Ns = zeros(N_bins)
	
	for i in 1:N_bins
		filt = bins[i] .<= log10.(rs) .< bins[i+1]
		σ = LilGuys.std(snap.velocities[3, filt], snap.weights[filt])

		Ns[i] = sum(filt)
		σv[i] = σ
	end

	10 .^ midpoints(bins), σv, bins, Ns
end

# ╔═╡ c1a6ddbc-d45d-4a14-b890-8ec95f0a64d6
sigma_prof = calc_velocity_profile(snap_i)

# ╔═╡ 85650542-649f-41f6-ba19-3fdd64601c4a
md"""
# Comparisons
"""

# ╔═╡ 7f9c948b-4ee4-4ecb-834c-3c165a1d8c83
Base.@kwdef struct DoubleExp2D <: LilGuys.SphericalProfile
	R_s::Float64
	M::Float64

	R_s_outer::Float64
	M_outer::Float64
end

# ╔═╡ 33db3986-0db3-4717-b635-b58845f4f837
function LilGuys.density(prof::DoubleExp2D, r::Real)
	return LilGuys.density(LilGuys.Exp2D(prof.M, prof.R_s, ), r) + LilGuys.density(LilGuys.Exp2D(prof.M_outer, prof.R_s_outer), r)
end

# ╔═╡ 13ebcd2c-81c2-49ae-815f-955a62e2932e
function rho_sigma2_r(halo, prof, r)
	integrand(r) = LilGuys.density(prof, r) * LilGuys.mass(halo, r) / r^2

	return LilGuys.integrate(integrand, r, Inf)
end

# ╔═╡ 35a5458b-0ceb-49aa-8cae-d5212bb4a329
function lerp_rho_sigma2(halo, prof, rs)
	ys = rho_sigma2_r.(halo, prof, rs)
	return LilGuys.lerp(rs, ys)
end

# ╔═╡ eeb26b49-5377-4ba6-8362-5c2d98ccb63c
rho_sigma2_r(halo, prof, 0)

# ╔═╡ fdf65715-ae1d-465e-ab72-811568f42ab1
function LilGuys.surface_density(prof::DoubleExp2D, r::Real)
	return LilGuys.surface_density(LilGuys.Exp2D(prof.M, prof.R_s), r) + LilGuys.surface_density(LilGuys.Exp2D(prof.M_outer, prof.R_s_outer), r)
end

# ╔═╡ 290f1743-a8c6-4932-8333-7ecd88d564a0
function sigma_los(halo::LilGuys.SphericalProfile, prof, R)
	integrand(r) = rho_sigma2_r(halo, prof, r) * r / sqrt(r^2 - R^2)

	Sigma_sigma2 = 2*LilGuys.integrate(integrand, R* (1+1e-10), Inf)

	return sqrt(Sigma_sigma2 / LilGuys.surface_density(prof, R))
end

# ╔═╡ faf98802-44cb-4cda-a066-ecc9a92437a0
function sigma_los(rho_sigma2, prof, R)
	integrand(r) = rho_sigma2(r) * r / sqrt(r^2 - R^2)

	Sigma_sigma2 = 2*LilGuys.integrate(integrand, R * (1+1e-2), Inf)

	return sqrt(Sigma_sigma2 / LilGuys.surface_density(prof, R))
end

# ╔═╡ 94209581-6029-446c-8538-b58dfe67849d
sigma_los(halo, prof, 20) * V2KMS

# ╔═╡ 385892d5-d9e1-4ea3-9b26-02893c7aaca0
prof_double = DoubleExp2D(M=1-f_outer, R_s=0.10, M_outer=f_outer, R_s_outer=0.22)

# ╔═╡ a1622346-d34c-4553-b21d-e59ceda6c947
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xscale = log10,
			 yscale =  log10,
			  xlabel = "R / kpc",
			  ylabel = L"$\sigma_\textrm{v}$ / km s$^{-1}$",
			  xticks = Makie.automatic,
			  yticks = 5:15,
			 )


	lines!(sigma_prof[1], sigma_prof[2] * V2KMS)
	
	x = logrange(0.01, 100, 100)
	v = sigma_los.(halo, prof, x)
	lines!(x, v*V2KMS)


	xlims!(0.01, 3)
	ylims!(5, 15)

	fig
end

# ╔═╡ 4a83e9d5-b6b1-454b-a870-2ac6b5e022ec
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xscale = log10,
			 yscale =  log10,
			  xlabel = "R / kpc",
			  ylabel = L"$\sigma_\textrm{v}$ / km s$^{-1}$",
			  xticks = Makie.automatic,
			  yticks = 5:15,
			 )

	
	sigma_prof = calc_velocity_profile(snap_plummer)
	lines!(sigma_prof[1], sigma_prof[2] * V2KMS)
	
	x = logrange(0.01, 100, 100)
	v = sigma_los.(halo, LilGuys.Plummer(r_s=0.20), x)
	lines!(x, v*V2KMS)


	xlims!(0.01, 3)
	ylims!(5, 15)

	fig
end

# ╔═╡ c27ba154-10bc-421d-b76a-e611fb656406
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xscale = log10,
			 yscale =  log10,
			  xlabel = "R / kpc",
			  ylabel = L"$\sigma_\textrm{v}$ / km s$^{-1}$",
			  xticks = Makie.automatic,
			  yticks = 5:15,
			 )

	
	sigma_prof = calc_velocity_profile(snap_double)
	lines!(sigma_prof[1], sigma_prof[2] * V2KMS)
	
	x = logrange(0.01, 100, 100)
	v = sigma_los.(halo, prof_double, x)
	lines!(x, v*V2KMS)


	xlims!(0.01, 3)
	ylims!(5, 15)

	fig
end

# ╔═╡ 0e239ccd-507b-456e-bf9f-94235cf9d0fb
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xscale = log10,
			 yscale =  log10,
			  xlabel = "R / kpc",
			  ylabel = L"$\Sigma$",
			  xticks = Makie.automatic,
			  # yticks=5:5:20
			 )


	
	x = logrange(0.01, 10, 100)
	y = LilGuys.density.(prof_double, x)
	lines!(x, y)

	fig

end

# ╔═╡ fb250f85-2084-4f61-97cc-d0aaf051ed4a
let
	fig = Figure(size=(4, 3) .* 72)
	ax = Axis(fig[1,1], 
			 xscale = log10,
			 yscale =  log10,
			  xlabel = "R / kpc",
			  ylabel = L"$\sigma_\textrm{v}$ / km s$^{-1}$",
			  yticks=5:5:20,
			 )


	
	x = logrange(0.01, 100, 100)
	v = sigma_los.(halo, prof, x)
	lines!(x, v*V2KMS, label="exp inner")

	r_scale = 2.2

	v = sigma_los.(halo, LilGuys.scale(prof, r_scale, 1), x)
	lines!(x, v*V2KMS, label="exp outer")


	v = sigma_los.(halo, prof_double, x)
	lines!(x, v*V2KMS, label="double exp")

	v = sigma_los.(halo, LilGuys.Plummer(r_s=0.20), x)
	lines!(x, v*V2KMS, label="Plummer")


	# Legend(fig[1, 2],ax)
	axislegend()
	xlims!(0.01, 3)
	ax.xticks[] = [0.01, 0.1, 1]
	ylims!(8, 20)

	fig
end

# ╔═╡ 9fed0d38-97f1-4f8b-b86b-2fd5fb2fb579
x = LinRange(0.01, 1.9, 100)

# ╔═╡ 308bd799-12c4-45c5-87b8-1ff9e81511bd
l = lerp_rho_sigma2(halo, prof, logrange(0.01, 10, 1000))

# ╔═╡ 4fa05a52-a10a-4cb5-8256-ae5789a51a65
sigma_los(l, prof, 1.9)

# ╔═╡ c1f1e0ae-518d-4215-8a88-9a1fcc0d8055
v = sigma_los.(halo, prof, x)

# ╔═╡ ba778dea-e5a5-4b9c-ad02-445da29708f2
md"""
# Double Exponential
"""

# ╔═╡ b1436cbd-eeba-4483-93d3-64e8807997b1
import Distributions

# ╔═╡ d61eb615-84be-4ca0-98f7-7e3e8a4720eb
let
	fig = Figure()

	ax = Axis(fig[1,1])

	x = snap_double.velocities[1, :] / LilGuys.std(snap_double.velocities[1, :], snap_double.weights)
	
	bins, values, _ = LilGuys.histogram(x, weights=snap_double.weights, errors=:weighted, normalization=:pdf)

	lines!(midpoints(bins), values)


	x = snap_i.velocities[1, :] / LilGuys.std(snap_i.velocities[1, :], snap_i.weights)

	bins, values, _ = LilGuys.histogram(x,  weights=snap_i.weights, errors=:weighted, normalization=:pdf)

	lines!(midpoints(bins), values)

	lines!(Distributions.Normal(), color=COLORS[3])

	fig

end

# ╔═╡ Cell order:
# ╠═8a84cd7a-b75e-11f0-ac6e-5386c69b3015
# ╟─9850dd29-11ec-4872-90b3-7161e61a2cac
# ╠═19866432-ee61-4891-88b8-6c2ab3f20c64
# ╠═bb15d2b6-313c-4f2c-877f-4a4aceb683f0
# ╠═13ebcd2c-81c2-49ae-815f-955a62e2932e
# ╠═35a5458b-0ceb-49aa-8cae-d5212bb4a329
# ╠═eeb26b49-5377-4ba6-8362-5c2d98ccb63c
# ╠═290f1743-a8c6-4932-8333-7ecd88d564a0
# ╠═faf98802-44cb-4cda-a066-ecc9a92437a0
# ╠═94209581-6029-446c-8538-b58dfe67849d
# ╠═0d6a1b04-7f20-42c3-870c-7f6e76ceea17
# ╠═90f7565d-6ee8-4ed8-bfcb-f2e97cfa9252
# ╠═39f091fe-c7c3-4427-8b8e-0dc7245db018
# ╠═6775faf6-8bba-4285-9414-4c5c8263e6be
# ╠═7e12fb10-688f-4305-8922-a9f48db7680e
# ╠═de504a7f-a151-490b-9535-4122a630a714
# ╠═aa68ef1a-0f77-4a2f-b6d2-eb6e8deb01c5
# ╠═c1a6ddbc-d45d-4a14-b890-8ec95f0a64d6
# ╠═85650542-649f-41f6-ba19-3fdd64601c4a
# ╠═7f9c948b-4ee4-4ecb-834c-3c165a1d8c83
# ╠═33db3986-0db3-4717-b635-b58845f4f837
# ╠═fdf65715-ae1d-465e-ab72-811568f42ab1
# ╠═385892d5-d9e1-4ea3-9b26-02893c7aaca0
# ╠═a1622346-d34c-4553-b21d-e59ceda6c947
# ╠═4a83e9d5-b6b1-454b-a870-2ac6b5e022ec
# ╠═c27ba154-10bc-421d-b76a-e611fb656406
# ╠═0e239ccd-507b-456e-bf9f-94235cf9d0fb
# ╠═fb250f85-2084-4f61-97cc-d0aaf051ed4a
# ╠═9fed0d38-97f1-4f8b-b86b-2fd5fb2fb579
# ╠═308bd799-12c4-45c5-87b8-1ff9e81511bd
# ╠═4fa05a52-a10a-4cb5-8256-ae5789a51a65
# ╠═c1f1e0ae-518d-4215-8a88-9a1fcc0d8055
# ╟─ba778dea-e5a5-4b9c-ad02-445da29708f2
# ╠═b1436cbd-eeba-4483-93d3-64e8807997b1
# ╠═d61eb615-84be-4ca0-98f7-7e3e8a4720eb
