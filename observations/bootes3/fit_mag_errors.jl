### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 8c8fb8de-35e4-4efe-aa0e-aa8d41743752
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya, CairoMakie
	FIGDIR = "./figures"
end

# ╔═╡ c6774d3d-96cc-468f-a282-24a9420ccd68
using PyFITS

# ╔═╡ d93400fb-9e35-4337-ae37-8f5b67d2f39c
md"""
The goal of this notebook is to fit the photometric uncertainties with simple functional forms
"""

# ╔═╡ ea62c9fc-7787-4e0d-b35d-f2665a8ba460
CairoMakie.activate!(type=:png)

# ╔═╡ 66ccecca-1db1-4e8f-b1df-789505b0b028
module Utils
	include("delve_utils.jl")
end

# ╔═╡ 2e34c0d6-a3d2-4a2c-81d3-0ba38624623b
⊕(x, y) = sqrt(x^2 + y^2)

# ╔═╡ 9c25bb60-a514-41c6-974d-3e80696acc69
delve = let
	df = read_fits("data/delve_dr2_good_6deg.fits")
	df
end

# ╔═╡ 07878b62-cb99-4a64-9d08-fb5df7f0403f
hist(delve.magerr_psf_g .⊕ delve.magerr_psf_r)

# ╔═╡ 0707b3fd-191d-486a-b583-616067f83c52
delve_gr_err = Utils.fit_color_err(delve.gmag, delve.magerr_psf_g .⊕ delve.magerr_psf_r)

# ╔═╡ 8ba4c6e6-b42b-4557-b042-f01259a3fdaa
delve_g_err = Utils.fit_color_err(delve.mag_psf_g, delve.magerr_psf_g)

# ╔═╡ 1bfb0b4d-8d73-4da6-bed2-7418698aab35
delve_g_err(29)

# ╔═╡ f5e08f1a-6b59-41a3-a2c8-378c1a8a676e
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "DELVE gmag", 
			 ylabel = "delve g-r err")
	scatter!(delve.gmag, delve.magerr_psf_g .⊕ delve.magerr_psf_r, markersize=0.3, color=:black, alpha=0.03)

	x = LinRange(extrema(delve.gmag)..., 1000)
	lines!(x, Utils.delve_gr_err.(x))
	fig
end

# ╔═╡ dcea0cbe-97d7-49a7-8e17-db536d231b23
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "DELVE gmag", 
			 ylabel = "delve g err")
	scatter!(delve.gmag, delve.magerr_psf_g , markersize=0.3, color=:black, alpha=0.03)

	x = LinRange(extrema(delve.gmag)..., 1000)
	lines!(x, Utils.delve_g_err.(x))
	fig
end

# ╔═╡ 78c20ae2-384c-49aa-a7b5-78f1bd285d34
gaia = let
	df = read_fits("data/j24_1c.fits")
	df[df.F_BEST .== 1, :]
end

# ╔═╡ 6b754609-26b3-4ba3-b087-baac8b23b5e3
gaia_bp_rp_err = Utils.fit_color_err(gaia.phot_g_mean_mag,  gaia.dBP .⊕ gaia.dRP)

# ╔═╡ 1406c198-60e4-491a-830e-bb3a6e36852d
gaia_G_err = Utils.fit_color_err(gaia.phot_g_mean_mag, gaia.dG)

# ╔═╡ 2ca9f103-0b24-4dc0-ad8e-81f1d64e7970
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 xlabel = "gaia G",
			 ylabel = "gaia d G",
			 limits=(nothing, (-0.001, 0.015))
	)
	
	scatter!(gaia.phot_g_mean_mag, gaia.dG, markersize=0.5, color=:black, alpha=0.03)

	x = LinRange(extrema(gaia.phot_g_mean_mag)..., 1000)
	lines!(x, Utils.gaia_G_err.(x))
	fig
end

# ╔═╡ 1d1645d1-65b3-4574-b584-a16ac112d630
let
	fig = Figure()
	ax = Axis(fig[1,1], 
			 limits=(nothing, (-0.01, 0.4)),
			 xlabel = "gaia G",
			 ylabel = "gaia d bp-rp")
	
	scatter!(gaia.phot_g_mean_mag, gaia.dBP .⊕ gaia.dRP, markersize=0.5, color=:black, alpha=0.03)

	x = LinRange(extrema(gaia.phot_g_mean_mag)..., 1000)
	lines!(x, Utils.gaia_bp_rp_err.(x))
	fig
end

# ╔═╡ Cell order:
# ╠═d93400fb-9e35-4337-ae37-8f5b67d2f39c
# ╠═8c8fb8de-35e4-4efe-aa0e-aa8d41743752
# ╠═ea62c9fc-7787-4e0d-b35d-f2665a8ba460
# ╠═c6774d3d-96cc-468f-a282-24a9420ccd68
# ╠═66ccecca-1db1-4e8f-b1df-789505b0b028
# ╠═2e34c0d6-a3d2-4a2c-81d3-0ba38624623b
# ╠═07878b62-cb99-4a64-9d08-fb5df7f0403f
# ╠═0707b3fd-191d-486a-b583-616067f83c52
# ╠═8ba4c6e6-b42b-4557-b042-f01259a3fdaa
# ╠═6b754609-26b3-4ba3-b087-baac8b23b5e3
# ╠═1406c198-60e4-491a-830e-bb3a6e36852d
# ╠═f5e08f1a-6b59-41a3-a2c8-378c1a8a676e
# ╠═1bfb0b4d-8d73-4da6-bed2-7418698aab35
# ╠═dcea0cbe-97d7-49a7-8e17-db536d231b23
# ╠═2ca9f103-0b24-4dc0-ad8e-81f1d64e7970
# ╠═1d1645d1-65b3-4574-b584-a16ac112d630
# ╠═9c25bb60-a514-41c6-974d-3e80696acc69
# ╠═78c20ae2-384c-49aa-a7b5-78f1bd285d34
