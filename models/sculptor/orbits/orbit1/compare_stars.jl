### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ e9d818c0-1897-11ef-1913-1b38c9d4b4d6
begin 
	import Pkg; Pkg.activate()
	
	using CairoMakie
	using Measurements

	
	import LilGuys as lguys
	using Arya
	
	using JSON
end

# ╔═╡ 1541f287-a351-41e1-8ea5-d78316202804
function load_profile(filename)
	local x
	open(filename, "r") do f
		x = JSON.parse(f)
	end

	for key in keys(x)
		col = x[key]
		if typeof(col) <: AbstractArray
			x[key][x[key] .=== nothing] .= NaN
			x[key] = Float64.(x[key])
		end
	end
	
	return x
end

# ╔═╡ 49dfa506-3ab2-44ad-949e-2ce709d393c4
begin 
	profiles = Dict(
		"big" => load_profile("exp2d_big_stars_today_profile.json"),
		"simulation" => load_profile("exp2d_stars_today_profile.json"),
		"small" => load_profile("exp2d_small_stars_today_profile.json"),
	)
end

# ╔═╡ 6d207751-695f-4247-94da-ced5a146092f
prof_expected = load_profile("/cosma/home/durham/dc-boye1/sculptor/fiducial_sample_profile.json")

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

	for (label, profile) in profiles
	
		lines!(profile["log_r"], profile["log_Sigma"], 
			label=label)
	end
	errscatter!(prof_expected["log_r"], prof_expected["log_Sigma"],
		yerr=prof_expected["log_Sigma_err"],
		label="J+24"
	)
	
	
	vlines!(log10(295))
	axislegend(ax)

	fig
end

# ╔═╡ Cell order:
# ╠═e9d818c0-1897-11ef-1913-1b38c9d4b4d6
# ╠═1541f287-a351-41e1-8ea5-d78316202804
# ╠═49dfa506-3ab2-44ad-949e-2ce709d393c4
# ╠═6d207751-695f-4247-94da-ced5a146092f
# ╠═de69d265-cda1-4543-ad02-2ee3091964d6
# ╠═932c4fef-992b-4518-80d0-59c8e126ccb5
