### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ e189b810-3436-11f1-bacc-7b98308ba75d
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya, CairoMakie

	import CSV
	import DataFrames: DataFrame
	import TOML

	import StatsBase: median
end

# ╔═╡ c10ae5af-d203-4d72-aa2f-28cfeba0522e
module Utils # only used for obs_dir and kroupa_imf
	obs_dir = joinpath(ENV["DWARFS_ROOT"], "observations/bootes3")
	include(joinpath(obs_dir, "delve_utils.jl"))
end

# ╔═╡ 1ba72e18-5744-4d3f-9160-3c51ff0f608d
CairoMakie.activate!(type=:png)

# ╔═╡ 71395e6f-ddd5-43ac-8a80-fbef8be9bf46
obs_dir = Utils.obs_dir

# ╔═╡ 035dad86-95ef-4784-9392-4cc88dcebf77
all_isochrones_gaia = Dict(
	age => Utils.read_mist_file(joinpath(obs_dir, "../padova", "isochrone.gaiadr3.$(age)Gyrs.dat")) for age in [8, 9, 10, 12]
);

# ╔═╡ 3e809fef-e29a-4440-a3c9-8c2369ef0e8c
function get_isochrone(isochrones, age, M_H; stage_max=5)

	if age ∉ keys(isochrones)
		throw(KeyError("$age not in available ages"))
	end
	
	isos = isochrones[age]
	@assert isapprox(isos.logAge[1],  log10(age) + 9, atol=0.01)
	if (M_H < minimum(isos.MH)) || (M_H > maximum(isos.MH))
		throw(DomainError(M_H, "metallicity out of isochrone range: $(extrema(isos.MH))"))
	end

	M_Hs = unique(isos.MH)
	M_H_adopted = M_Hs[argmin(abs.(M_H .- M_Hs))]

	filt = isapprox.(isos.MH, M_H_adopted)
	filt .&= isos.label .<= stage_max # only keep through HB

	return isos[filt, :]
end

# ╔═╡ ed6cdd98-2df6-4e81-9f3e-d3b892127110
iso = get_isochrone(all_isochrones_gaia, 12, -2.19)

# ╔═╡ d74fe6d3-356c-4c1c-996d-333db8c1180e
obs_props = TOML.parsefile(joinpath(obs_dir, "observed_properties.toml"))

# ╔═╡ fa4fb6e3-ba1b-496e-9ad2-ad60cba86583
dm_0 = obs_props["distance_modulus"]

# ╔═╡ 6a5bb2e3-ea16-4723-af51-7473619cbd88
mass_sampler = Utils.create_kroupa_sampler(iso.Mini[end])

# ╔═╡ b3267dee-0509-4a63-9232-72cc57656c5a
function interpolate_magnitude(iso, mag_col)
	return LilGuys.lerp(iso.Mini, iso[!, mag_col])
end

# ╔═╡ 822b8561-2dfe-4c66-8735-7765501b7b0c
let
	fig = Figure(size=(4, 3) .* 72)
	ax = Axis(fig[1,1],
			 yreversed=true, 
			 limits=(nothing, nothing, 15, 22),
			 ylabel = "G", 
			 xlabel = "BP - RP")

	lines!(iso.G_BPmag .- iso.G_RPmag, iso.Gmag .+ dm_0)


	ax_hist = Axis(fig[1, 2], yreversed=true,  limits=(nothing, nothing, 15, 23),
				  xlabel = "density")


	masses = [mass_sampler() for _ in 1:1_000_000]
	mags = interpolate_magnitude(iso, "Gmag").(masses) .+ dm_0

	M_mean = LilGuys.mean(masses)
	@info M_mean

	@info LilGuys.mean((mags .< 21))
	hist!(mags, direction=:x, bins=100, normalization=:pdf)

	m = LinRange(0.08, maximum(iso.Mini), 10000)
	w = Utils.kroupa_imf.(m)

	f = LilGuys.lerp(iso.Mini, iso.Gmag .+ dm_0)

	mags = f.(m)
	@info sum(w[(mags .< 21)]) / sum(w)

	hist!(mags, bins=100, weights=w, normalization=:pdf, color=:transparent,  label = "isochrone mag(kroupa IMF)", direction=:x, strokewidth=1)


	
	xlims!(0, 0.003)
	linkyaxes!(ax, ax_hist)



	fig
end

# ╔═╡ 6062be48-e24a-4b71-b7cb-e7961fe7e7c9
110 * 120

# ╔═╡ Cell order:
# ╠═e189b810-3436-11f1-bacc-7b98308ba75d
# ╠═c10ae5af-d203-4d72-aa2f-28cfeba0522e
# ╠═1ba72e18-5744-4d3f-9160-3c51ff0f608d
# ╠═71395e6f-ddd5-43ac-8a80-fbef8be9bf46
# ╠═035dad86-95ef-4784-9392-4cc88dcebf77
# ╠═3e809fef-e29a-4440-a3c9-8c2369ef0e8c
# ╠═ed6cdd98-2df6-4e81-9f3e-d3b892127110
# ╠═d74fe6d3-356c-4c1c-996d-333db8c1180e
# ╠═fa4fb6e3-ba1b-496e-9ad2-ad60cba86583
# ╠═6a5bb2e3-ea16-4723-af51-7473619cbd88
# ╠═b3267dee-0509-4a63-9232-72cc57656c5a
# ╠═822b8561-2dfe-4c66-8735-7765501b7b0c
# ╠═6062be48-e24a-4b71-b7cb-e7961fe7e7c9
