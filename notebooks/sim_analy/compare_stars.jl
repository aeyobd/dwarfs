### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ e9d818c0-1897-11ef-1913-1b38c9d4b4d6
begin 
	import Pkg; Pkg.activate()
	
	using CairoMakie
	
	import LilGuys as lguys
	using Arya
end

# ╔═╡ 96ae7f95-886b-4fd6-a1a0-3d4f0790cc34
cd("/astro/dboyea/sculptor/orbits/1e6/orbit1/V70_r0.4/stars/")

# ╔═╡ 88f31bfc-a67d-4654-af8b-46dc91500558
r_b = 76

# ╔═╡ 658ea28e-b0c9-433f-80e6-71d868f882a8
Rs = [
	"0.02",
	"0.07", 
	"0.13"
]

# ╔═╡ 51af74d4-2643-4d20-a0da-01f3cd0ff264
names_f = ["exp2d_rs$(r)_stars_profile.toml" for r in Rs]

# ╔═╡ c441f04d-9793-4315-a29e-3c248dd66fd3
names_i = ["exp2d_rs$(r)_stars_1_profile.toml" for r in Rs]

# ╔═╡ 49dfa506-3ab2-44ad-949e-2ce709d393c4
profiles = lguys.ObsProfile.(names_f)

# ╔═╡ bd5c7170-ed86-4273-a775-4f574c28932c
profiles_i = lguys.ObsProfile.(names_i)

# ╔═╡ 6d207751-695f-4247-94da-ced5a146092f
prof_expected = lguys.ObsProfile("/astro/dboyea/sculptor/fiducial_sample_profile.toml")

# ╔═╡ de69d265-cda1-4543-ad02-2ee3091964d6
log_r_label = "log r / arcmin"

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
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=log_r_label,
		ylabel = L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
		limits=((-1.5, 3), (-5, 1))
	)

	errscatter!(prof_expected.log_r, prof_expected.log_Sigma,
		yerr=prof_expected.log_Sigma_err,
		label="J+24",
		color=:black
	)

	
	for i in eachindex(Rs)
		profile = profiles[i]

		label = "$(Rs[i])"
		lines!(profile.log_r, profile.log_Sigma, 
			label=label)
	end
	

	
	vlines!(log10(r_b), color=:grey, label="break radius")
	axislegend("Rs / kpc")

	fig
end

# ╔═╡ 91dcf440-0660-445c-874d-9ebdef4b36ab
let 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=log_r_label,
		ylabel = L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
		limits=((-1.5, 3), (-5, 1))
	)

	
	for i in eachindex(Rs)
		profile = profiles_i[i]

		label = "$(Rs[i])"
		lines!(profile.log_r, profile.log_Sigma .+ 0i, 
			label=label)
	end
	

	
	axislegend("Rs / kpc")

	fig
end

# ╔═╡ 64baeb89-4cbf-421c-97f6-d2a624b30882
let 
	for i in eachindex(Rs)
		fig = Figure()
		ax = Axis(fig[1,1], 
			xlabel=log_r_label,
			ylabel = L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
			limits=((-1.5, 3), (-5, 1)),
			title = "Rs = $(Rs[i]) kpc"
		)
	

		lines!(profiles_i[i].log_r, profiles_i[i].log_Sigma, 
				label="initial", linestyle=:dot)
		
		lines!(profiles[i].log_r, profiles[i].log_Sigma, 
				label="final")

		
	
		
		axislegend()
	
		fig
		@info fig
	end
end

# ╔═╡ 5d03ccd2-2dd6-45d9-8ec1-1af9bea3475b
md"""
# Scratch
"""

# ╔═╡ Cell order:
# ╠═e9d818c0-1897-11ef-1913-1b38c9d4b4d6
# ╠═96ae7f95-886b-4fd6-a1a0-3d4f0790cc34
# ╠═88f31bfc-a67d-4654-af8b-46dc91500558
# ╠═658ea28e-b0c9-433f-80e6-71d868f882a8
# ╠═51af74d4-2643-4d20-a0da-01f3cd0ff264
# ╠═c441f04d-9793-4315-a29e-3c248dd66fd3
# ╠═49dfa506-3ab2-44ad-949e-2ce709d393c4
# ╠═bd5c7170-ed86-4273-a775-4f574c28932c
# ╠═bfab4ae8-a94b-4a81-aad3-9b706f2474bb
# ╠═6d207751-695f-4247-94da-ced5a146092f
# ╠═de69d265-cda1-4543-ad02-2ee3091964d6
# ╠═932c4fef-992b-4518-80d0-59c8e126ccb5
# ╠═91dcf440-0660-445c-874d-9ebdef4b36ab
# ╠═64baeb89-4cbf-421c-97f6-d2a624b30882
# ╟─5d03ccd2-2dd6-45d9-8ec1-1af9bea3475b
