### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 9a85e4f6-2960-11f1-a0af-dd4aa796aa14
begin
	import Pkg; Pkg.activate()

	using Arya, CairoMakie; CairoMakie.activate!(type=:png)
	using LilGuys
	import TOML
end

# ╔═╡ f3c25159-2321-4265-8a43-3c93ee75a5e1
using PyFITS

# ╔═╡ 83ad8ecc-b8bd-4f79-9380-393d72a2b612
using Measurements: ±

# ╔═╡ b01546fe-c3f9-4b68-a86e-0195762e4c3f
using Turing

# ╔═╡ a9480d56-d6f8-45eb-96d7-2e57da8d72a3
using PairPlots

# ╔═╡ 9d66e7e9-004d-4bf1-bd34-a689b07b9426
using DataFrames

# ╔═╡ 6862f562-bb38-40b0-83b6-205037de6b40
using KernelDensity

# ╔═╡ 8649ca7f-9e6a-4c23-bed8-a9181c42be0e
md"""
# Finding a stream

This notebook aims to search for a possible stream around Boo III, as reported in Grillmair2009. 

While there is a hint of a similar stream, we find that the stream is likely explained by the pattern of extinction nearby to Boo III. Even though the stream roughly aligns with the proper moteion of Boo III (``\theta \approx 80^\circ``), we only find a 1-2``\sigma`` overdensity. As such, we conclude the strema is likely not visible in the current DELVE dataset.
"""

# ╔═╡ e1911072-a189-4577-a817-d9ffadd0afa8
import StatsBase: median, mad, std

# ╔═╡ 6cd21893-a5f3-4f8a-8036-82f455c0fde8
FIGDIR = "./figures"

# ╔═╡ 6c28b525-dd73-4a6b-8deb-933e1d76ab31
stars = read_fits("samples/delve_matched_filter.fits")

# ╔═╡ 9903762b-0cb9-408d-aa84-d0ee4c924f1b
icrs = ICRS(TOML.parsefile("observed_properties.toml"))

# ╔═╡ 477c460e-bd11-49a6-8105-708471845fd4
md"""
To better produce histograms of the across-stream stellar density, we mask out the regions of the stream nearby to the Globular Clusters NGC 5466 and NGC 5272
"""

# ╔═╡ ac40e18b-1d33-4d87-9837-72f263641018
props_ngc5466 = (TOML.parsefile("observed_properties_ngc5466.toml"))

# ╔═╡ 455e1f2a-bb73-4f89-9134-986d34cf0e84
props_ngc5272 = (TOML.parsefile("observed_properties_ngc5272.toml"))

# ╔═╡ 23e8fe64-1fdb-47bb-8ed2-5db7845201df
R_h_ngc5272 = props_ngc5272["R_h"] / 60

# ╔═╡ 5f343cf5-29db-48e0-863c-cb7a6a2d7209
R_h_ngc5466 = props_ngc5466["R_h"] / 60

# ╔═╡ cf204b42-19cc-4ac6-8e9e-a1dc41ddaf03
md"""
As a quick check, we compute the 
"""

# ╔═╡ 579d69fc-1938-4fb0-89f6-bbef01c8a1e6
pmra_err, pmdec_err = [TOML.parsefile("observed_properties.toml")[k] for k in ["pmra_em", "pmdec_em"]]

# ╔═╡ fcb870e5-d831-4398-8754-932f1fe456eb
gsr = LilGuys.transform(GSR, icrs);

# ╔═╡ 7ea273b6-52d6-41b3-bfc4-7dd8858e1ce5
θ_pm = atand(gsr.pmra ± pmra_err, gsr.pmdec ± pmdec_err)

# ╔═╡ 42d25b54-345a-4cf9-a5be-35ea379a68b0
w_all = (stars.L_sat ./ stars.L_bg)

# ╔═╡ ccaa7118-7cf5-458a-8ed3-70e02c598e35
md"""
# The best case

Below, we explore the best candidate for a stream direction. 
We define stream coordinates ``\xi'`` and ``\eta'`` which are merely rotated versions of the ordinay tangent plane coordinates, according to some angle.


"""

# ╔═╡ 9def7bd4-6f54-417c-8343-7282af3a4f63
function tangent_axis(gs)
	Axis(gs, 
		xlabel = "xi / arcmin",
		ylabel = "eta / arcmin",
		 xreversed = true,
		 aspect = DataAspect()
		)
end

# ╔═╡ f9f3e7f7-c5f0-4ecb-af1c-26f75cd98d21
function stream_axis(gs)
	Axis(gs, 
		xlabel = L"$\xi'$ / deg",
		ylabel = L"$\eta'$ / deg",
		 xreversed = true,
		 aspect = DataAspect()
		)
end

# ╔═╡ c5928b42-7a8d-4a73-894f-f00368c22298
R_h = 30 / 60

# ╔═╡ 7a945aef-0a3d-41a3-baa5-a8778d21e08a
function rotated_tangent_sample(θ_stream; η_range=(-4, 4), ξ_range = (-6, 6), filter=true, r_mag_max=24)
	transform(ra, dec) = LilGuys.to_orbit_coords(ra, dec, icrs.ra, icrs.dec, θ_stream)
	
	ξ_all, η_all = transform(stars.ra, stars.dec)
	ξ_ngc5466, η_ngc5466 = transform(props_ngc5466["ra"], props_ngc5466["dec"])
	
	ξ_ngc5272, η_ngc5272 = transform(props_ngc5272["ra"], props_ngc5272["dec"])

	if !filter
		return ξ_all, η_all, w_all
	end

	filt = (ξ_range[1] .< ξ_all .< ξ_range[2]) .& (η_range[1] .< η_all .< η_range[2])
	filt .&= @. !(-10R_h < ξ_all < 10R_h)
	filt .&= @. !(ξ_ngc5272 - 10R_h_ngc5272 < ξ_all < ξ_ngc5272 + 10R_h_ngc5272)
	filt .&= @. !(ξ_ngc5466 - 10R_h_ngc5466 < ξ_all < ξ_ngc5466 + 10R_h_ngc5466)
	
	filt .&= isfinite.(w_all)
	filt .&= stars.mag_psf_r .< r_mag_max

	return ξ_all[filt], η_all[filt], w_all[filt]
end
	

# ╔═╡ 072a5cb6-3139-4b1c-a77f-dc632ab51013
θ_example = 70

# ╔═╡ 0dfaea3d-8068-4a36-8a84-a2e5cbc79e51
ξ_all, η_all, _ = rotated_tangent_sample(θ_example, filter=false)

# ╔═╡ 5a80f035-7e11-4cc1-8a8d-c92843995c86
ξ, η, w = rotated_tangent_sample(θ_example, η_range=(-6, 5.5), ξ_range=(0, 10), r_mag_max=22)

# ╔═╡ 266493ce-8b6b-49b5-9a97-0fa233523aa0
md"""
The left plot below scatter stars, with markersize weighted by the satellite likelihood. The right plot shows the SFD extinction in the region (linear scaling of A(g))

The bottom two are identical, except showing the prefered direction of the example stream in orange
"""

# ╔═╡ 49d1d9c0-d350-465e-a129-23ac468a5438
let
	fig = Figure(size=(6, 6) .* 72, rasterize=true)
	ax = tangent_axis(fig[1,1])
	filt = abs.(η_all) .< 0.2
		
	scatter!(stars.xi, stars.eta, markersize=sqrt.(w_all) ./ 6, alpha=1, color=:black, rasterize = true)
	hidexdecorations!(ticks=false, minorticks=false)


	ax = tangent_axis(fig[1,2])
	scatter!(stars.xi, stars.eta, markersize=1, alpha=1, color=stars.A_g, rasterize = true)
	hidedecorations!(ticks=false, minorticks=false, )

	
	ax = tangent_axis(fig[2,1])
		
	scatter!(stars.xi, stars.eta, markersize=sqrt.(w_all) ./ 6, alpha=1, color=:black, rasterize = true)

	scatter!(stars.xi[filt], stars.eta[filt], markersize=0.5, alpha=0.03, rasterize = true)

	ax = tangent_axis(fig[2,2])
	scatter!(stars.xi, stars.eta, markersize=1, alpha=1, color=stars.A_g, rasterize = true)

	scatter!(stars.xi[filt], stars.eta[filt], markersize=0.5, alpha=0.03, rasterize = true)

	hideydecorations!(ticks=false, minorticks=false)

	# @savefig "stream_tangent_w_extinction"

	fig
end

# ╔═╡ afd1f0d3-fc0a-446f-ab5c-52d9f5776145
let
	fig = Figure(size=(6, 3) .* 72)
	ax = stream_axis(fig[1,1])
	ax.title = "All"
					 
	scatter!(ξ_all, η_all, markersize=sqrt.(w_all) ./ 5, alpha=0.3, color=:black)


	ax2 = stream_axis(fig[1,2])
	ax2.title = "selected for hist"

	scatter!(ξ, η, markersize=sqrt.(w) ./ 5, alpha=0.3, color=:black)
	hideydecorations!(ticks=false, minorticks=false)

	linkaxes!(ax, ax2)
	fig

end

# ╔═╡ d02597d7-95e0-40fc-bb2d-57cfc1a1d008
md"""
The below plot illustrates the density of stars in each bin on the off-stream cordinate ``\eta'``. Note that the unweighted sample gently slopes across the field (with the decrease starting at 4 deg due to the shape of the selected stars). While there is a mild overdensity from about -2 to 1 deg in the samples weighted by CMD likelihood or with CMD likelihood > 3, the natural variation nearby means this peak is likely only marginally significant. 

Given that the peak corresponds to the reduction in extinction, we conclude this overdensity likely arises from extinction-based incompleteness effects near the expected main sequence of Boo III. 
"""

# ╔═╡ a743fc66-782e-4cc2-929c-e13152d7ce72
let
	fig = Figure()
	ax = Axis(fig[1,1], ylabel="density", xlabel=L"$\eta'$ / deg")

	stephist!(η, bins=40, weights=w, normalization=:pdf, label="weighted")
	stephist!(η, bins=40, weights=w .> 3, normalization=:pdf, label="CMD L >3")
	stephist!(η, bins=40, normalization=:pdf, label="unweighted")


	axislegend(position=:lb)
	fig
end

# ╔═╡ 7ea39641-d7bb-4ebf-8d45-d93f991edef5
md"""
# Quick likelihood fit
"""

# ╔═╡ 085e3621-1bf5-43ef-875b-42ea1cdd66c5
η_max = 3.5

# ╔═╡ 95902f07-87a2-458f-a466-5814d1758f8c
@model function linear_gaussian_model(η, w)
	a ~ Normal(0, 0.1)
	L_linear = @. (1+a*η) / (2*η_max) 

	μ ~ Uniform(-4, 4)
	σ ~ LogUniform(0.1, 5)
	f ~ Uniform(-1, 1)
	# L_gaus = 0 #@. A*LilGuys.gaussian.(η, μ, σ)
	L_gaus = @. 1/sqrt(2π) /σ * exp(-(η -μ)^2/(2σ^2))

	L_tot = L_linear .* (1-f) .+ L_gaus .* w .* f
	Turing.@addlogprob!(sum(log10.(max.(L_tot, 1e-12))))
end

# ╔═╡ c029b4ec-d00b-4871-8573-c9fdff6af379
filt_model = (w .> 1e-4) .& (abs.(η) .< η_max)

# ╔═╡ 04168d83-b6e6-4912-ad10-19945ea29f0a
g_model = linear_gaussian_model(η[filt_model], w[filt_model])

# ╔═╡ 18ea911f-81c4-407b-b277-52bec9c45958
samples = sample(g_model, NUTS(), 1000)

# ╔═╡ 13d75431-de77-4588-b76d-23a36f96cd71
pairplot(samples)

# ╔═╡ e8c072db-2c48-48ab-ac64-3515845a9c9e
df_samples = DataFrame(samples)

# ╔═╡ 0c84fca8-8d3d-46c5-b6c8-428f7e072a26
mean(df_samples.f * sum(filt_model) .< 300) 

# ╔═╡ 5205f4ac-9416-4f47-9802-dc63846cda9e
sum(filt_model)

# ╔═╡ 5401956e-9503-4a2b-8c79-54ae24fca1fd
sum(filt_model) .* median(df_samples.f)

# ╔═╡ 5170c376-2a3f-457a-9cf0-0beab95c4c0c
median(df_samples.f)

# ╔═╡ a42ea27b-a9cb-4cf1-a222-97b1eb9c57f5


# ╔═╡ 971f77bf-1a95-4f33-8d3d-263e9fb22741
let
	fig = Figure()
	ax = Axis(fig[1,1])

	for i in 1:100
		x = LinRange(-3.5, 3.5, 1000)
		L_linear = @. 1 + df_samples.a[i]*x
		
		L_gaus = @. 1 /sqrt(2π) /df_samples.σ[i] * exp(-(x -df_samples.μ[i])^2/(2df_samples.σ[i]^2))

		L_tot = L_linear.* (1 - df_samples.f[i]) .+ L_gaus*df_samples.f[i]
		norm = sum(L_tot .* LilGuys.gradient(x))

		lines!(x, L_tot./ norm)
		lines!(x, L_gaus*df_samples.f[i] ./ norm)

	end
	stephist!(η[filt_model], bins=60, weights=w[filt_model], normalization=:pdf, label="weighted")

	fig
end


# ╔═╡ 77cfde22-a64a-401b-8510-0870ee86c930
md"""
# Changing the angle to search for streams

When changing the angle systematically, and creating density plots as above (either by likelihood weighting or cutting in likelihood), the strongest overdensity isthe one near 65-70 degrees we explored earlier. 
"""

# ╔═╡ ded6e443-2d09-45d0-b544-79542e95b3dc
let
	fig = Figure(size=(2, 8) .* 72)
	ax = Axis(fig[1,1])

	dy = 0.1
	
	for (i, θ) in enumerate(0:5:180)
		ξ, η, w = rotated_tangent_sample(θ)

		h =  KernelDensity.kde(η, weights=w, bandwidth=6/60, boundary=(-4, 4))

		lines!(h.x, h.density .+ dy*i)
		μ = median(h.density[abs.(h.x) .> 1])
		σ = std(h.density[abs.(h.x) .> 1])
		hlines!(μ + dy*i)
		hspan!(μ + dy*i - σ, μ + dy*i + σ, alpha=0.2)
		text!(-4, median(h.density) + dy*i, text="$θ", color=COLORS[(i-1) % 10 + 1 ])
	end

	fig
end

# ╔═╡ b3c1b2dc-57e7-4bf3-823d-468c37f30545
let
	fig = Figure(size=(2, 8) .* 72)
	ax = Axis(fig[1,1])

	dy = 0.1
	
	for (i, θ) in enumerate(0:5:180)
		ξ, η, w = rotated_tangent_sample(θ)

		h =  KernelDensity.kde(η, weights=w .> 3, bandwidth=6/60, boundary=(-4, 4))

		lines!(h.x, h.density .+ dy*i)
		μ = median(h.density[abs.(h.x) .> 1])
		σ = std(h.density[abs.(h.x) .> 1])
		hlines!(μ + dy*i)
		hspan!(μ + dy*i - σ, μ + dy*i + σ, alpha=0.2)
		text!(-4, median(h.density) + dy*i, text="$θ", color=COLORS[(i-1) % 10 + 1 ])
	end

	fig
end

# ╔═╡ fec6a72f-9591-4bf7-85f8-5ba6817371de
md"""
These plots simply verify that the masking technique doesn't leave the GCs in the selected region
"""

# ╔═╡ 835df023-1f5a-455d-a415-e0e93ff1966f
let
	fig = Figure(size=(3, 10) .* 72)

	dy = 100
	
	for (i, θ) in enumerate((0:45:360) .+ 25)
		ξ, η, w = rotated_tangent_sample(θ)

		ax = Axis(fig[i, 1], title="$θ")
		scatter!(ξ, η, markersize=sqrt.(w) ./ 3, alpha=0.3, color=:black)
	end

	fig
end

# ╔═╡ Cell order:
# ╟─8649ca7f-9e6a-4c23-bed8-a9181c42be0e
# ╠═9a85e4f6-2960-11f1-a0af-dd4aa796aa14
# ╠═e1911072-a189-4577-a817-d9ffadd0afa8
# ╠═6cd21893-a5f3-4f8a-8036-82f455c0fde8
# ╠═f3c25159-2321-4265-8a43-3c93ee75a5e1
# ╠═6c28b525-dd73-4a6b-8deb-933e1d76ab31
# ╠═9903762b-0cb9-408d-aa84-d0ee4c924f1b
# ╠═477c460e-bd11-49a6-8105-708471845fd4
# ╠═ac40e18b-1d33-4d87-9837-72f263641018
# ╠═455e1f2a-bb73-4f89-9134-986d34cf0e84
# ╠═23e8fe64-1fdb-47bb-8ed2-5db7845201df
# ╠═5f343cf5-29db-48e0-863c-cb7a6a2d7209
# ╠═cf204b42-19cc-4ac6-8e9e-a1dc41ddaf03
# ╠═579d69fc-1938-4fb0-89f6-bbef01c8a1e6
# ╠═fcb870e5-d831-4398-8754-932f1fe456eb
# ╠═83ad8ecc-b8bd-4f79-9380-393d72a2b612
# ╠═7ea273b6-52d6-41b3-bfc4-7dd8858e1ce5
# ╠═42d25b54-345a-4cf9-a5be-35ea379a68b0
# ╟─ccaa7118-7cf5-458a-8ed3-70e02c598e35
# ╠═9def7bd4-6f54-417c-8343-7282af3a4f63
# ╠═f9f3e7f7-c5f0-4ecb-af1c-26f75cd98d21
# ╠═7a945aef-0a3d-41a3-baa5-a8778d21e08a
# ╠═c5928b42-7a8d-4a73-894f-f00368c22298
# ╠═072a5cb6-3139-4b1c-a77f-dc632ab51013
# ╠═0dfaea3d-8068-4a36-8a84-a2e5cbc79e51
# ╠═5a80f035-7e11-4cc1-8a8d-c92843995c86
# ╟─266493ce-8b6b-49b5-9a97-0fa233523aa0
# ╠═49d1d9c0-d350-465e-a129-23ac468a5438
# ╠═afd1f0d3-fc0a-446f-ab5c-52d9f5776145
# ╟─d02597d7-95e0-40fc-bb2d-57cfc1a1d008
# ╠═a743fc66-782e-4cc2-929c-e13152d7ce72
# ╠═7ea39641-d7bb-4ebf-8d45-d93f991edef5
# ╠═b01546fe-c3f9-4b68-a86e-0195762e4c3f
# ╠═95902f07-87a2-458f-a466-5814d1758f8c
# ╠═085e3621-1bf5-43ef-875b-42ea1cdd66c5
# ╠═c029b4ec-d00b-4871-8573-c9fdff6af379
# ╠═04168d83-b6e6-4912-ad10-19945ea29f0a
# ╠═18ea911f-81c4-407b-b277-52bec9c45958
# ╠═a9480d56-d6f8-45eb-96d7-2e57da8d72a3
# ╠═13d75431-de77-4588-b76d-23a36f96cd71
# ╠═9d66e7e9-004d-4bf1-bd34-a689b07b9426
# ╠═e8c072db-2c48-48ab-ac64-3515845a9c9e
# ╠═0c84fca8-8d3d-46c5-b6c8-428f7e072a26
# ╠═5205f4ac-9416-4f47-9802-dc63846cda9e
# ╠═5401956e-9503-4a2b-8c79-54ae24fca1fd
# ╠═5170c376-2a3f-457a-9cf0-0beab95c4c0c
# ╠═a42ea27b-a9cb-4cf1-a222-97b1eb9c57f5
# ╠═971f77bf-1a95-4f33-8d3d-263e9fb22741
# ╟─77cfde22-a64a-401b-8510-0870ee86c930
# ╠═6862f562-bb38-40b0-83b6-205037de6b40
# ╠═ded6e443-2d09-45d0-b544-79542e95b3dc
# ╠═b3c1b2dc-57e7-4bf3-823d-468c37f30545
# ╠═fec6a72f-9591-4bf7-85f8-5ba6817371de
# ╠═835df023-1f5a-455d-a415-e0e93ff1966f
