### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ e9d818c0-1897-11ef-1913-1b38c9d4b4d6
begin 
	import Pkg; Pkg.activate()
	
	using CairoMakie
	
	import LilGuys as lguys
	using Arya
end

# ╔═╡ 2fbbea8f-1a5a-46eb-89b8-8e64265cd66d
cd("/astro/dboyea/sculptor/isolation/1e4/stars")

# ╔═╡ 49dfa506-3ab2-44ad-949e-2ce709d393c4
begin 
	profiles = Dict(
		#"0.05" => lguys.ObsProfile("exp2d_rs0.05_today_profile.toml"),
		"0.1" => lguys.ObsProfile("exp2d_rs0.1_mock_stars_profile.toml"),
		"0.15" => lguys.ObsProfile("exp2d_rs0.15_mock_stars_profile.toml"),
		"0.2" => lguys.ObsProfile("exp2d_rs0.2_mock_stars_profile.toml"),
		"0.5" => lguys.ObsProfile("exp2d_rs0.5_mock_stars_profile.toml"),
		"1" => lguys.ObsProfile("exp2d_rs1_mock_stars_profile.toml"),
	) 
end

# ╔═╡ 421f707a-94de-4791-8f16-b9ea38fbfbf1
begin 
	ks = collect(keys(profiles))
	Rss = parse.(Float64, ks)
	ks = ks[sortperm(Rss)]
	Rss = sort(Rss)
end

# ╔═╡ 6d207751-695f-4247-94da-ced5a146092f
prof_expected = lguys.ObsProfile("/astro/dboyea/sculptor/fiducial_sample_profile.toml")

# ╔═╡ de69d265-cda1-4543-ad02-2ee3091964d6
log_r_label = "log r / arcmin"

# ╔═╡ d3e934c2-cc73-4873-8157-8303161ceecf
distance=86

# ╔═╡ 76ca8391-6f20-4d50-b0fc-3fd9ff0edd1b
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
	

# ╔═╡ c3916053-02ed-4b84-afed-faee4190e4be
let
	fig, ax = FigAxis(
		xlabel=log_r_label,
		ylabel=L"\Sigma/\Sigma_0",
		limits=(-1, 3.5, nothing, nothing)
	)

	for i in 1:length(Rss)
		Rs = Rss[i]
		
		plot_model_and_exp!(profiles[ks[i]], Rs, y_offset=i, label=ks[i])
	end

	axislegend(L"$R_s$ / kpc", position=:lb)
	fig
end

# ╔═╡ Cell order:
# ╠═e9d818c0-1897-11ef-1913-1b38c9d4b4d6
# ╠═2fbbea8f-1a5a-46eb-89b8-8e64265cd66d
# ╠═49dfa506-3ab2-44ad-949e-2ce709d393c4
# ╠═76ca8391-6f20-4d50-b0fc-3fd9ff0edd1b
# ╠═c3916053-02ed-4b84-afed-faee4190e4be
# ╠═421f707a-94de-4791-8f16-b9ea38fbfbf1
# ╠═6d207751-695f-4247-94da-ced5a146092f
# ╠═de69d265-cda1-4543-ad02-2ee3091964d6
# ╠═d3e934c2-cc73-4873-8157-8303161ceecf
