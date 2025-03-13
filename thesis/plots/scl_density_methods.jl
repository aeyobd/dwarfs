### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ e7ae2044-a31e-4a98-b496-3b1b3be056f0
using OrderedCollections

# ╔═╡ 53a18bc2-3bbd-4fe3-8a77-d9c49d2ddea9
include("paper_style.jl")

# ╔═╡ 0ac92e62-f8fd-4e17-8b2f-feecdab75497
import TOML

# ╔═╡ b11f3588-9904-46ef-8f35-a7c1d623e020
log_r_ell_label = L"\log\ r_\textrm{ell}\,/\,\textrm{arcmin}"

# ╔═╡ d42a0cd3-cc8e-4a24-8887-7100f3927961
log_Σ_label = L"$\log \Sigma$\,/\,stars\ arcmin$^{-2}$"

# ╔═╡ dedba1c4-b769-4bd1-a005-0385a1a1ffa8
function get_num_per_bin(x)
	num_per_bin = round(Int64, LilGuys.Interface.default_n_per_bin(x))
	if length(x) < 2
		num_per_bin = 1
	end
	return num_per_bin
end


# ╔═╡ 972f8fe3-a819-402d-97ec-c0dcbdc07f5f
datafile = ENV["DWARFS_ROOT"] * "/observations/sculptor/data/jensen+24_wide.fits"

# ╔═╡ 389fefd6-60fd-4dd8-ba77-29e87b4ed846
observed_properties = TOML.parsefile(ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml")

# ╔═╡ 412b80f4-106e-44f6-bbf9-8f29b62e2fd5
begin 
	all_stars = LilGuys.read_fits(datafile)

	all_stars[:, :r_ell_old] = all_stars[:, :r_ell]
	all_stars.xi, all_stars.eta = LilGuys.to_tangent(all_stars.ra, all_stars.dec, observed_properties["ra"], observed_properties["dec"])

	all_stars.r_ell = 60LilGuys.calc_r_ell(all_stars.xi, all_stars.eta, observed_properties["ellipticity"], observed_properties["position_angle"])
	all_stars
end

# ╔═╡ 8b51c994-b4a2-41dd-85f6-c5f96bdd1ea8
best_stars = all_stars[all_stars.F_BEST .== 1., :]

# ╔═╡ db45fd11-ca8b-46e1-a4ac-af426410a1ba
members = best_stars[best_stars.PSAT .> 0.2, :]

# ╔═╡ fb319d6b-57e1-41af-af33-46931d16818f
r_ell_max = (1-observed_properties["ellipticity"]) * maximum(all_stars.r_ell)

# ╔═╡ d63dd2d8-b4f4-46ef-8745-d06d5cbfeb19
function get_bins(r_ell; bin_width = 0.05, num_per_bin=nothing)

	if num_per_bin === nothing
		num_per_bin = get_num_per_bin(log10.(r_ell))
	end

	num_per_bin = round(Int, num_per_bin)
	bins = LilGuys.Interface.bins_both(log10.(r_ell), nothing, bin_width=bin_width, num_per_bin=num_per_bin)
	bins[end] = min(bins[end], log10(r_ell_max))

	return bins
end

# ╔═╡ 16a1e52a-2ce4-45d6-8c25-bb7c82ff718f
function plot_density_for_filter(filters; sequential=false, normalization=:none, axislabel="")
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel = log_r_ell_label,
		ylabel = log_Σ_label,
	)

	names = collect(keys(filters))
	for i in 1:length(filters)
		name = names[i]
		x = filters[name]
		
		x = x[r_ell_max .> x .> 0]
		if length(x) < 2
			println("skipping with < 2 members: ", "\t", name)
			continue
		end
		
		println(name, "\t", length(x))
		bins = get_bins(x)

		if length(bins) < 2
			println("skipping with < 2 bins: ", "\t", name)
			continue
		end
		
		prof = LilGuys.StellarProfile(x, normalization=normalization, bins=(bins))

		if sequential
			 kwargs = (; color=i, colorrange=(1, length(names)))
			else
			 kwargs = (; )
		end

		if name == "bright"
			prof.log_Sigma .+= log10(2)
		end
		
		scatterlines!(prof.log_r, prof.log_Sigma, label=name; kwargs...)
	end


	return fig
end

# ╔═╡ c51c6f13-b5ed-4671-9f08-4adb4f7b38b7
import Statistics: median

# ╔═╡ 82aee929-129d-490c-99c3-e504d0ac5309
LL_NOSPACE = log10.(best_stars.L_CMD_SAT ./ best_stars.L_CMD_BKD .* best_stars.L_PM_SAT ./ best_stars.L_PM_BKD)

# ╔═╡ 2b83d8d5-ded6-4bfd-b77c-726484203988
LL = log10.(best_stars.L_SAT ./ best_stars.L_BKD)

# ╔═╡ 72895e1e-afab-468a-89bb-9aa1e5665269
L_min = minimum(LL[best_stars.PSAT .> 0.2])

# ╔═╡ a567a789-2116-4d6e-a04f-c1edfe74a4aa
minimum(LL_NOSPACE[best_stars.PSAT_NOSPACE .> 0.2])

# ╔═╡ 95695b01-d0a0-4bac-be51-e00798363c4a
sum(LL_NOSPACE .> L_min)

# ╔═╡ ca32f73a-f43a-4b48-a17d-66f0e5f1edc9
sum(best_stars.PSAT_NOSPACE .> 0.0001)

# ╔═╡ e7f92f8e-fdc3-4970-af5d-55b002f16038
hist(LL_NOSPACE[LL_NOSPACE .> -10])

# ╔═╡ ce4c3ccf-0f54-40dc-a6b7-e2a93f9f2247
members_nospace = best_stars[LL_NOSPACE .> L_min, :]

# ╔═╡ 30a60b7d-630a-4fa3-ac57-e873655ad754
let

	fig = plot_density_for_filter(OrderedDict(
		"fiducial" => members.r_ell,
		"circular" => @.(60 * sqrt(members.xi^2 + members.eta^2)),
		"CMD+PM" => members_nospace.r_ell,
		"bright" => members.r_ell[members.phot_g_mean_mag .< median(members.phot_g_mean_mag)],
		"PSAT > 0.8" => members.r_ell[members.PSAT .> 0.8],
		),
		#normalization=:mass,
	)
	Legend(fig[1,1], fig.content[1], halign=0.1, valign=0.1, tellwidth=false)

	@savefig "scl_density_methods"

	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═0ac92e62-f8fd-4e17-8b2f-feecdab75497
# ╠═53a18bc2-3bbd-4fe3-8a77-d9c49d2ddea9
# ╠═b11f3588-9904-46ef-8f35-a7c1d623e020
# ╠═d42a0cd3-cc8e-4a24-8887-7100f3927961
# ╠═d63dd2d8-b4f4-46ef-8745-d06d5cbfeb19
# ╠═dedba1c4-b769-4bd1-a005-0385a1a1ffa8
# ╠═972f8fe3-a819-402d-97ec-c0dcbdc07f5f
# ╠═389fefd6-60fd-4dd8-ba77-29e87b4ed846
# ╠═412b80f4-106e-44f6-bbf9-8f29b62e2fd5
# ╠═8b51c994-b4a2-41dd-85f6-c5f96bdd1ea8
# ╠═db45fd11-ca8b-46e1-a4ac-af426410a1ba
# ╠═fb319d6b-57e1-41af-af33-46931d16818f
# ╠═16a1e52a-2ce4-45d6-8c25-bb7c82ff718f
# ╠═e7ae2044-a31e-4a98-b496-3b1b3be056f0
# ╠═c51c6f13-b5ed-4671-9f08-4adb4f7b38b7
# ╠═82aee929-129d-490c-99c3-e504d0ac5309
# ╠═2b83d8d5-ded6-4bfd-b77c-726484203988
# ╠═72895e1e-afab-468a-89bb-9aa1e5665269
# ╠═a567a789-2116-4d6e-a04f-c1edfe74a4aa
# ╠═95695b01-d0a0-4bac-be51-e00798363c4a
# ╠═ca32f73a-f43a-4b48-a17d-66f0e5f1edc9
# ╠═e7f92f8e-fdc3-4970-af5d-55b002f16038
# ╠═ce4c3ccf-0f54-40dc-a6b7-e2a93f9f2247
# ╠═30a60b7d-630a-4fa3-ac57-e873655ad754
