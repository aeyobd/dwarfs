### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	
	using CSV, DataFrames
	using Arya
	using CairoMakie

	using PyFITS
	import LilGuys as lguys
end

# ╔═╡ 7ed5bcf5-dfc3-4e79-a608-d503124a1e96
using LilGuys

# ╔═╡ 68da1f4c-93d8-4125-8792-188bcf15b306
using Measurements

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Radial Velocity CrossMatch

This notebook loads in radial velocity data preprocessed and matched to gaia from different sources and combines these into a signle catalogue.
"""

# ╔═╡ 9485165e-0cb3-452d-a0ec-953b790c9d7f
FIGDIR = "figures"

# ╔═╡ 9e9ba645-b780-4afa-b305-a2b1d8a97220
import StatsBase: quantile, mean, std, median, kurtosis, sem

# ╔═╡ 257caa92-ccce-44ff-88ee-6a9c550ae5d9
CairoMakie.activate!(type=:png)

# ╔═╡ ef1bb6f5-00b8-405b-8702-244eca618644
import DensityEstimators: histogram, calc_limits, make_bins

# ╔═╡ 58caad92-a887-4527-ac84-be02fba5b23c
import TOML

# ╔═╡ 9a59505b-039b-41e1-8907-a8107ce68177
module RVUtils
	include("../../rv_utils.jl")
end

# ╔═╡ d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
md"""
# Loading data tables
"""

# ╔═╡ 3e0eb6d1-6be4-41ec-98a5-5e9167506e61
data_dir = "../data/"

# ╔═╡ 77e7884c-0360-4b7f-b9ad-e81be2516552
obs_properties = TOML.parsefile("../observed_properties.toml")

# ╔═╡ 7a50a176-96a5-4098-88d6-0fa2874d0f90
j24 = read_fits("../samples/best_sample.fits")

# ╔═╡ d4d0a488-1c0e-4bc6-9a88-94b27d84e8ce
li_all = let
	df = read_fits("processed/rv_li.fits")

	df[!, "RV_sigma"] .= 0
	df[!, "RV_count"] .= 1

	df
end

# ╔═╡ 9ed91d41-e650-4d3c-9e74-100ff3d57d82
geha_all = let
	df = read_fits("processed/rv_deimos.fits")
	filt = df.system .== Ref("Boo3")
	filt .&= .!ismissing.(df.source_id )
	df[!, "RV_sigma"] .= 0
	df[!, "RV_count"] .= 1
	df[filt, :]
end

# ╔═╡ bbf49122-11b1-4272-a660-0437c6aa2b3f
md"""
# Combined stars
"""

# ╔═╡ d3333b48-aa4e-42c1-9e0a-bbff98e3647d
all_studies = ["geha26", "li26"]

# ╔═╡ fe6fa910-f41e-4657-836b-7eda2f0cddb2
function add_xmatch!(df, new, suffix)
	leftjoin!(df, rename(n->"$(n)_$suffix", new), on="source_id"=>"source_id_$suffix")
end

# ╔═╡ e472cbb6-258e-4303-85e2-56f26358c97b
all_stars = let

	all_stars = copy(j24)[:, [:source_id, :ra, :dec, :dr2_radial_velocity, :dr2_radial_velocity_error]]
	add_xmatch!(all_stars, li_all, "li26")
	add_xmatch!(all_stars, geha_all, "geha26")

	rename!(all_stars,
		:dr2_radial_velocity => "RV_gaia",
		:dr2_radial_velocity_error => "RV_err_gaia",
	)

	RVUtils.add_rv_means!(all_stars, all_studies)
	all_stars

end

# ╔═╡ 3e2a9c87-e464-4c84-95cf-d70f0cece9e2
filt_is_meas = all_stars.RV_nstudy .> 0

# ╔═╡ 73bd1553-e2ae-4bfb-aac1-0880346f5054
rv_meas = all_stars[filt_is_meas, :]

# ╔═╡ 76d38672-0a84-449e-bc21-83f3bba01d26
md"""
## double checks
"""

# ╔═╡ 8f243944-3d0b-48ce-9f90-8d448c089239
md"""
- double check we have all stars
"""

# ╔═╡ ab49efb3-2ab7-47d0-a3a5-a342c789aa9b
@assert length(li_all.source_id) == sum(.!ismissing.(rv_meas.RV_li26))

# ╔═╡ f505b96f-d134-4b9f-a63c-b275f16be5d6
@assert length(geha_all.source_id) == sum(.!ismissing.(rv_meas.RV_geha26))

# ╔═╡ 8862360f-0888-4d3c-8ae0-6274c32d1818
md"""
- double check that the F_match corresponds to actual kept xmatches
"""

# ╔═╡ de762a39-b430-4452-ba87-8b8cf1ad9852
@assert sum(geha_all.F_match) ==  sum(.!ismissing.(rv_meas.RV_geha26))

# ╔═╡ 04ff1abd-c584-41de-8c83-7503482c3731
md"""
- Are any gaia RVs included?
Gaia stars are too bright to be members :(
"""

# ╔═╡ bebba70b-3d89-48a6-9fe8-53ca4428f241
PSAT = [only(j24.PSAT[j24.source_id .== row.source_id]) for row in eachrow(rv_meas)]

# ╔═╡ 0f50ae12-e207-4748-9c45-780c9215be73
sum(skipmissing(.!ismissing.(rv_meas.RV_err_gaia) .* (PSAT .> 001)))

# ╔═╡ 103f7c18-fefa-4d3c-a1b2-3324cc2239d3
R_ell = [only(j24.R_ell[j24.source_id .== row.source_id]) for row in eachrow(rv_meas)]

# ╔═╡ db3177e7-7132-4f62-846a-f4416a804009
md"""
# write data
"""

# ╔═╡ 615f05fb-4907-40bc-9251-3065f565929b
let
	filename = "processed/rv_combined.fits"
	if isfile(filename)
		rm(filename); 
	end
	rv_meas[!, :F_scatter] .= true
	rv_meas[!, :F_match] .= true
	write_fits(filename, rv_meas)
	@info "wrote data"
end

# ╔═╡ 23d894e1-d976-45d0-bacb-286cc0d6915c
md"""
# Counts
"""

# ╔═╡ ad832280-d456-411d-9199-47a1a2909c79
sum(rv_meas.RV_count)

# ╔═╡ 2ce6f17b-38f4-4a60-8e11-510650b83f04
length(rv_meas.RV)

# ╔═╡ 8162e31f-9b76-4348-8ae4-d1e66d6ac429
median(rv_meas.RV_err)

# ╔═╡ f7ec8bba-9f45-435b-b67c-33182e992dfd
md"""
# Cross study RV
"""

# ╔═╡ 93d185f2-1eaa-4d35-87dd-b84f385483de
function filt_missing(col, verbose=false; low=-Inf, high=Inf)
	filt =  @. !ismissing(col) & !isnan(col)
	filt1 = high .> col .> low
	if verbose
		println("excluding ", sum(.!(filt1)[filt]), " outliers")
	end
	return filt .& filt1
end

# ╔═╡ f69aeefd-5834-4b3b-b332-4d89801f2763
rv_low = 150

# ╔═╡ e5b82f39-a45b-426d-9a9e-c203670959c5
rv_high = 250

# ╔═╡ 609f55cd-aa49-4a5e-b871-413ec7ef0990
import StatsBase as sb

# ╔═╡ 9213cc36-74d1-452f-bd9a-eb5c5cab1f87
function compare_rv(study1, study2)
	rv1  = all_stars[:, "RV_$study1"]
	rv1_err  = all_stars[:, "RV_err_$study1"]
	
	rv2  = all_stars[:, "RV_$study2"]
	rv2_err  = all_stars[:, "RV_err_$study2"]

	filt = filt_missing(rv1, true) .& filt_missing(rv2, true)

	println("matched ", sum(filt), " stars")

	println(rv1[filt], " pm ", rv1_err[filt])
	println(rv2[filt], " pm ", rv2_err[filt])

	if sum(filt) == 0
		println("nothing to plot")
		return
	end
	
	println("plotting ", sum(filt), " stars")


	fig, ax = FigAxis(
		xlabel = "RV $study1",
		ylabel = "RV $study2",
		aspect=DataAspect()
	)

	errorscatter!(disallowmissing(rv1[filt]), disallowmissing(rv2[filt]), xerror=rv1_err[filt], yerror=rv2_err[filt])

	w = @. 1/(rv1_err[filt]^2 + rv2_err[filt]^2)
	μ = median(rv1[filt] .- rv2[filt])
	sigma = sb.mad(rv1[filt] .- rv2[filt] .- μ)
	δμ = sigma / sqrt(length(w))

	text!(0.1, 0.8, space=:relative, text="μ = $(μ ± δμ)\nσ = $(round(sigma, digits=2))", fontsize=10)


	scatter!(rv1[filt], rv1[filt], color=:black)
	return fig
end

# ╔═╡ f61debe4-8b23-4415-b77c-43c4465ccfcb
compare_rv("geha26", "li26")

# ╔═╡ 7c5cee19-8ff4-4a77-85c4-115273f08afe
print(intersect(Set(geha_all.source_id), Set(li_all.source_id)))

# ╔═╡ a784d463-891a-4ccc-971f-2792fcdfd722
geha_all[geha_all.source_id .== [1450782713560551808], :]

# ╔═╡ f73dd9b4-a52e-497d-93db-2a9f59fdf058
li_all[li_all.source_id .== [1450782713560551808], :]

# ╔═╡ 102c73ef-2c95-4784-80df-ed0312511c00


# ╔═╡ 3655b6a0-ac9a-4a49-86fd-6a6484409819
md"""
## Check RV means make sense
"""

# ╔═╡ 7f13339e-6a0c-4944-8a36-5ab136fd8415
function compare_rv_mean(study1, rv_meas=rv_meas)
	rv1  = rv_meas[:, "RV_$study1"]
	rv1_err  = rv_meas[:, "RV_err_$study1"]
	
	rv2  = rv_meas[:, "RV"]
	rv2_err  = rv_meas[:, "RV_err"]

	filt = filt_missing(rv1, true; low=rv_low, high=rv_high) .& filt_missing(rv2, true;  low=rv_low, high=rv_high)

	println("matched ", sum(filt), " stars")


	if sum(filt) == 0
		println("nothing to plot")
		return
	end
	
	println("plotting ", sum(filt), " stars")


	fig, ax = FigAxis(
		xlabel = "RV mean",
		ylabel = "RV $study1 - RV mean",
		#aspect=DataAspect()
	)

	x = float.(rv2[filt])
	y = float.(rv1[filt] .- rv2[filt])
	xerr = float.(rv2_err[filt])
	yerr = float.(rv1_err[filt])
	errorscatter!(x, y, xerror=xerr, yerror=yerr, color=(COLORS[1], 0.2), markersize=1)

	w = 1 ./ xerr .^2
	μ = LilGuys.mean(y, w)
	δμ = sqrt( 1 / sum(w))
	s = LilGuys.std(y, w)
	text!(0.1, 0.1, space=:relative, text="μ = $(μ ± δμ)\nσ = $(round(s, digits=2))", fontsize=10)

	hlines!(0, color=:black)
	return fig
end

# ╔═╡ 5688e9c9-515c-4973-820f-1215031253f2
import StatsBase: var

# ╔═╡ ce3067d3-55e6-43a1-9b7b-4cf53c09ee88
md"""
# J24 purity checks
"""

# ╔═╡ 01e0cfcd-92af-4c04-a7d4-9c9d8fd5a9f1
rv_mean = obs_properties["radial_velocity"]; sigma_los = obs_properties["sigma_v"] # from mcmc modeling in next section

# ╔═╡ 92a32c4c-59d8-4d3f-8709-d66fae10bfe9
rv_mean

# ╔═╡ f15ad626-bb0d-4484-8004-d6f81ccf0825
⊕ = RVUtils.:⊕

# ╔═╡ 03462abb-32d6-433c-a9c4-65ef2dd58965
function is_rv_member(rv, err)
	return abs(rv - rv_mean) / (err ⊕ sigma_los) < 3
end

# ╔═╡ 8ce732f4-3678-4cdd-bdf3-c6c2c34ee404
PSAT

# ╔═╡ 82d4e049-7f47-4f8b-a3a0-7191ef82912b
filt_rv = is_rv_member.(rv_meas.RV, rv_meas.RV_err) 

# ╔═╡ daa1295f-3e76-4794-ab6d-dba1b3b9f779
memb_filt = (abs.(rv_meas.RV .- rv_mean) .< 3*obs_properties["sigma_v"]) .&& filt_rv

# ╔═╡ 3920e27f-48fd-41a5-bf53-1f80edf3d56d
memb_stars = rv_meas[memb_filt, :]

# ╔═╡ aaaf5ba2-c9ed-41ec-a22a-d78ed96fd84e
compare_rv_mean("geha26", memb_stars)

# ╔═╡ 78e7aff0-3658-4d64-b117-58189a85307a
compare_rv_mean("li26", memb_stars)

# ╔═╡ 764b5306-20f9-4810-8188-1bdf9482260f
let
	fig, ax = FigAxis(
		xlabel=L"radial velocity / km s$^{-1}$",
		ylabel="density"
	)
	h = histogram(Float64.(memb_stars.RV), 30, normalization=:pdf)
	
	errorscatter!(midpoints(h.bins), h.values, yerror=h.err, color=COLORS[6])

	fig
end

# ╔═╡ e343dcd0-078f-476d-abe0-fd84001b72ae
mean(filt_rv[.!isnan.(rv_meas.RV_err)])

# ╔═╡ 0105a7ec-392e-426b-9bf9-a4b6144fe829
rv_meas.RV .- rv_mean

# ╔═╡ 7b22fa97-c311-42ae-ad5a-41d053c5334e
rv_meas.RV

# ╔═╡ e899faa9-580b-4aad-902e-382008048908
function purity_given_p(psatlow, psathigh)
	filt = psathigh .> PSAT .>= psatlow
	filt .&= .!isnan.(PSAT)
	filt .&= .!ismissing.(filt_rv)
	return sum(filt_rv[filt]) / sum(filt)
end

# ╔═╡ c0a4e207-9778-4d43-8197-5fc8887b2b82
function number_satisfying(psatlow, psathigh)
	filt = psathigh .> PSAT .>= psatlow
	filt .&= .!isnan.(PSAT)
	filt .&= .!ismissing.(filt_rv)

	return sum(filt)
end

# ╔═╡ 69161135-1f19-4d9a-8ff6-ae63a79e3eb5
purity_given_p(0.0, 1)

# ╔═╡ ab12510f-9393-4acd-8ed0-dcffa06c65e6
purity_given_p(0.0, 1.1)

# ╔═╡ 678a13d6-c892-43fe-a806-f3534661f785
let
	fig, ax = FigAxis(
		ylabel = "purity",
		limits=(0, 1, 0, 1.05)
	)

	ps = [0, 0.1, 0.3, 0.5, 0.7, 0.9, 1]

	y = purity_given_p.(ps[1:end-1], ps[2:end])
	lines!([0, 1], [0, 1], color=:black, linestyle=:dot)

	scatterlines!(midpoints(ps), y)

	hidexdecorations!(ax, ticks=false, minorticks=false)

	ax2 = Axis(fig[2, 1], ylabel = "# in bin", xlabel="PSAT",
			  yscale=log10, yticks=[1, 10],
			  )
	ylims!(1, nothing)

	y = number_satisfying.(ps[1:end-1], ps[2:end])

	lines!(midpoints(ps), y)

	rowgap!(fig.layout, 0)
	rowsize!(fig.layout, 2, Relative(1/3))
	
	@savefig "purity_vs_psat"
end

# ╔═╡ 7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
md"""
## Quick properties
"""

# ╔═╡ 33f54afc-cdb9-4eb8-887f-5a43281b837c
let
	fig = Figure(
		size=(5*72, 3*72)
	)
	ax = Axis(fig[1,1],
		xlabel = L"$R_\textrm{ell}$ / arcmin",
		ylabel = L"RV / km s$^{-1}$",
	)

	scatter!(R_ell[memb_filt], memb_stars.RV_li26, label="Li+26", markersize=3)
	scatter!(R_ell[memb_filt], memb_stars.RV_geha26, label="Geha+26", markersize=5)

	axislegend(ax, position=:rb)


	resize_to_layout!(fig)

	@savefig "rv_scatter_alstudy"
	fig
end

# ╔═╡ e000571c-c471-46af-92c4-2d16e39e66e2
let
	fig = Figure(size=(4, 4) .* 72)
	ax = Axis(fig[1,1], 
			 xlabel = "RA", 
			 ylabel = "Dec", 
			  title = "Boo III spectroscopic membs"
	)


	filt = .!ismissing.(memb_stars.RV_geha26)
	scatter!(memb_stars.ra[filt], memb_stars.dec[filt], alpha=0.5, label="Keck", marker=:utriangle)


	filt = .!ismissing.(memb_stars.RV_li26)
	scatter!(memb_stars.ra[filt], memb_stars.dec[filt], alpha=0.5, label="S5")

	filt = j24.PSAT .> 0.2
	scatter!(j24.ra[filt], j24.dec[filt], markersize=2, label="Gaia only", marker=:circle, color=:black, alpha=0.5)


	Legend(fig[2, 1], ax, tellwidth=false, tellheight=true)

	@savefig "ra_dec_spectroscopic"
	fig
end

# ╔═╡ Cell order:
# ╟─811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═9485165e-0cb3-452d-a0ec-953b790c9d7f
# ╠═9e9ba645-b780-4afa-b305-a2b1d8a97220
# ╠═257caa92-ccce-44ff-88ee-6a9c550ae5d9
# ╠═ef1bb6f5-00b8-405b-8702-244eca618644
# ╠═7ed5bcf5-dfc3-4e79-a608-d503124a1e96
# ╠═68da1f4c-93d8-4125-8792-188bcf15b306
# ╠═58caad92-a887-4527-ac84-be02fba5b23c
# ╠═9a59505b-039b-41e1-8907-a8107ce68177
# ╟─d4eb6d0f-4fe0-4e9d-b617-7a41f78da940
# ╠═3e0eb6d1-6be4-41ec-98a5-5e9167506e61
# ╠═77e7884c-0360-4b7f-b9ad-e81be2516552
# ╠═7a50a176-96a5-4098-88d6-0fa2874d0f90
# ╠═d4d0a488-1c0e-4bc6-9a88-94b27d84e8ce
# ╠═9ed91d41-e650-4d3c-9e74-100ff3d57d82
# ╟─bbf49122-11b1-4272-a660-0437c6aa2b3f
# ╠═d3333b48-aa4e-42c1-9e0a-bbff98e3647d
# ╠═e472cbb6-258e-4303-85e2-56f26358c97b
# ╠═fe6fa910-f41e-4657-836b-7eda2f0cddb2
# ╠═3e2a9c87-e464-4c84-95cf-d70f0cece9e2
# ╠═73bd1553-e2ae-4bfb-aac1-0880346f5054
# ╟─76d38672-0a84-449e-bc21-83f3bba01d26
# ╟─8f243944-3d0b-48ce-9f90-8d448c089239
# ╠═ab49efb3-2ab7-47d0-a3a5-a342c789aa9b
# ╠═f505b96f-d134-4b9f-a63c-b275f16be5d6
# ╟─8862360f-0888-4d3c-8ae0-6274c32d1818
# ╠═de762a39-b430-4452-ba87-8b8cf1ad9852
# ╟─04ff1abd-c584-41de-8c83-7503482c3731
# ╠═0f50ae12-e207-4748-9c45-780c9215be73
# ╠═bebba70b-3d89-48a6-9fe8-53ca4428f241
# ╠═103f7c18-fefa-4d3c-a1b2-3324cc2239d3
# ╟─db3177e7-7132-4f62-846a-f4416a804009
# ╠═615f05fb-4907-40bc-9251-3065f565929b
# ╟─23d894e1-d976-45d0-bacb-286cc0d6915c
# ╠═ad832280-d456-411d-9199-47a1a2909c79
# ╠═2ce6f17b-38f4-4a60-8e11-510650b83f04
# ╠═8162e31f-9b76-4348-8ae4-d1e66d6ac429
# ╟─f7ec8bba-9f45-435b-b67c-33182e992dfd
# ╠═93d185f2-1eaa-4d35-87dd-b84f385483de
# ╠═f69aeefd-5834-4b3b-b332-4d89801f2763
# ╠═e5b82f39-a45b-426d-9a9e-c203670959c5
# ╠═9213cc36-74d1-452f-bd9a-eb5c5cab1f87
# ╠═609f55cd-aa49-4a5e-b871-413ec7ef0990
# ╠═f61debe4-8b23-4415-b77c-43c4465ccfcb
# ╠═7c5cee19-8ff4-4a77-85c4-115273f08afe
# ╠═a784d463-891a-4ccc-971f-2792fcdfd722
# ╠═f73dd9b4-a52e-497d-93db-2a9f59fdf058
# ╠═102c73ef-2c95-4784-80df-ed0312511c00
# ╟─3655b6a0-ac9a-4a49-86fd-6a6484409819
# ╠═7f13339e-6a0c-4944-8a36-5ab136fd8415
# ╠═92a32c4c-59d8-4d3f-8709-d66fae10bfe9
# ╠═daa1295f-3e76-4794-ab6d-dba1b3b9f779
# ╠═3920e27f-48fd-41a5-bf53-1f80edf3d56d
# ╠═aaaf5ba2-c9ed-41ec-a22a-d78ed96fd84e
# ╠═78e7aff0-3658-4d64-b117-58189a85307a
# ╠═764b5306-20f9-4810-8188-1bdf9482260f
# ╠═5688e9c9-515c-4973-820f-1215031253f2
# ╟─ce3067d3-55e6-43a1-9b7b-4cf53c09ee88
# ╠═01e0cfcd-92af-4c04-a7d4-9c9d8fd5a9f1
# ╠═f15ad626-bb0d-4484-8004-d6f81ccf0825
# ╠═03462abb-32d6-433c-a9c4-65ef2dd58965
# ╠═8ce732f4-3678-4cdd-bdf3-c6c2c34ee404
# ╠═82d4e049-7f47-4f8b-a3a0-7191ef82912b
# ╠═e343dcd0-078f-476d-abe0-fd84001b72ae
# ╠═0105a7ec-392e-426b-9bf9-a4b6144fe829
# ╠═7b22fa97-c311-42ae-ad5a-41d053c5334e
# ╠═e899faa9-580b-4aad-902e-382008048908
# ╠═c0a4e207-9778-4d43-8197-5fc8887b2b82
# ╠═69161135-1f19-4d9a-8ff6-ae63a79e3eb5
# ╠═ab12510f-9393-4acd-8ed0-dcffa06c65e6
# ╠═678a13d6-c892-43fe-a806-f3534661f785
# ╟─7a1a920e-45e7-4d6f-925c-88dfb77f6dfb
# ╠═33f54afc-cdb9-4eb8-887f-5a43281b837c
# ╠═e000571c-c471-46af-92c4-2d16e39e66e2
