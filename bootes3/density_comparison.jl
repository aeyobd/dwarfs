### A Pluto.jl notebook ###
# v0.20.24

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

# ╔═╡ d3c09dc2-db77-471c-9384-3ffeb0df39dc
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
import TOML

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ e93918f6-0c72-4dec-9cdf-97f185c0bceb
module ModelUtils
	include("model_utils.jl")
end

# ╔═╡ 437ac44b-5440-4b3c-aa58-4bbd93be7174
function plot_prof!(prof; kwargs...)
	lines!(LilGuys.log_radii(prof), LilGuys.log_surface_density(prof); kwargs...)
end

# ╔═╡ 56365cfd-16e0-4567-a85c-f59372a611c7
function compare_density(gs, modelname_exp; modelname_plummer=nothing, limits=(-0.5, 3, -8, -1.0), lmc=false, y_low_R_h=nothing, y_low=nothing, color=COLORS[3],
						break_label_visible = true, jacobi_label_visible = true, norm_shift_exp=0)

	galaxy = modelname_exp[1]
	ax = Axis(gs,
		xlabel = "log Radii / arcmin",
		ylabel = "log surface density",
		limits = limits,
			 )
	ModelUtils.plot_expected_profile!(galaxy)

	prof_i, prof_f, _ = ModelUtils.load_stellar_profiles(modelname_exp..., norm_shift=norm_shift_exp)

	plot_prof!(prof_i, color=color, linewidth=1.5, linestyle=:dot, label="stars initial")
	plot_prof!(prof_f, color=color, linewidth=1.5, linestyle=:solid, label="stars final")


	if !isnothing(modelname_plummer)
		prof_i, prof_f, _ = ModelUtils.load_stellar_profiles(modelname_plummer...)
	
		plot_prof!(prof_i, color=COLORS[2], linestyle=:dot, label="2exp initial")
		plot_prof!(prof_f, color=COLORS[2], linestyle=:solid, label="2exp final")
	end


	limits!(limits...)
	
	if isnothing(y_low)
		y_low = limits[3]
	end
	r_b = ModelUtils.get_r_b(modelname_exp..., lmc=lmc)
	ModelUtils.plot_r_break_arrow!(r_b, y_low, break_label_visible)

	r_j = ModelUtils.get_r_j(modelname_exp[1:2]..., lmc=lmc)
	ModelUtils.plot_r_jacobi_arrow!(r_j, y_low, jacobi_label_visible)

	R_h = ModelUtils.get_R_h(galaxy)
	if isnothing(y_low_R_h)
		y_low_R_h = limits[3]
	end
	ModelUtils.plot_R_h_arrow!(R_h, y_low_R_h)
	
	ax
end

# ╔═╡ fc0a118f-f8b7-425e-aacf-c8c7630d61a5
function model_label!(text)
	text!(0, 0, space=:relative, text=text, align=(:left, :bottom), offset=(4, 4), fontsize=8)
end

# ╔═╡ f38390b2-b6dd-4c42-a548-fab9b6ac2a73
function compare_density(modelname; title="", label=nothing, kwargs...)
	fig = Figure()
	ax = compare_density(fig[1,1], modelname; kwargs...)

	ax.title = title
	if !isnothing(label)
		model_label!(label)
	end
	
	fig
end

# ╔═╡ 8c2d254f-7443-4ff1-ae5e-55b3bd0df766
modelnames = Dict(
	"1x12kpc" => ["bootes3", "1e6_v22_r3.9/1_peri_12kpc", "exp2d_rs0.20"],
	"1x1.5kpc" => ["bootes3", "1e6_v30_r3.0/1_peri_1.5kpc", "exp2d_rs0.20"],
	"3x26kpc" => ["bootes3", "1e6_v22_r3.9/3_peri_26kpc", "exp2d_rs0.20"],
	"5x18kpc" => ["bootes3", "1e6_v30_r3.0/5_peri_18kpc", "exp2d_rs0.20"],
	"5x18kpc-P" => ["bootes3", "1e6_v30_r3.0/5_peri_18kpc", "plummer_rs0.30"],
)

# ╔═╡ f1443328-05ae-41be-8242-d92ed430908b
compare_density(modelnames["1x12kpc"], norm_shift_exp=0.2, title="1x12kpc")

# ╔═╡ 080c9586-bb3d-4f80-b4ea-5bae53f174f4
compare_density(modelnames["5x18kpc-P"], norm_shift_exp=0.2, title="5x18kpc Plummer")

# ╔═╡ d26db585-3e68-4870-99da-bb5d9182c923
compare_density(modelnames["5x18kpc"], title="5x18kpc")

# ╔═╡ dbed5ec4-6030-4667-a87b-e8100deebe4a
compare_density(modelnames["1x1.5kpc"], title="1x1.5kpc")

# ╔═╡ 9e8271b9-173a-47d6-a144-7dbaeb3badc7
compare_density(modelnames["3x26kpc"], norm_shift_exp=0.3, title="3x26kpc")

# ╔═╡ 1f026ba8-0c58-4cf1-8b4e-c4344f7c8c72
md"""
# Boundmass
"""

# ╔═╡ 713590cd-0446-4096-b3b8-c0b8df850bea
function get_scalars(modelname)
	dir = ModelUtils.get_starsdir_out(modelname...)
	t_f_Gyr = TOML.parsefile(joinpath(dir, "../../orbital_properties.toml"))["t_f_gyr"]
	@info t_f_Gyr
	df = read_fits(joinpath(dir, "stellar_profiles_3d_scalars.fits"))
	df.time .-= t_f_Gyr ./ T2GYR
	return df
end

# ╔═╡ 604a69dd-1768-4156-9e65-9e4d3f8a120e
scalars = Dict(
	label => get_scalars(modelname) for (label, modelname) in modelnames
)

# ╔═╡ 898d0210-1da9-4dfd-84b8-114bd2511e62
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "time / Gyr",
			 ylabel = "bound mass fraction",
			 )

	for model in keys(modelnames)
		df = scalars[model]
		lines!(df.time * T2GYR, (df.bound_mass ./ df.bound_mass[1]), label=model)
	end

	ylims!(0, 1.05)
	vlines!(0, color=:black, linestyle=:dot)

	axislegend(position=:lb)
	fig
end

# ╔═╡ e66a81b1-36ba-48ca-be63-1658ec3804f6
let
	fig = Figure()
	ax = Axis(fig[1,1],
			 xlabel = "time / Gyr",
			 ylabel = "sigmav",
			 )

	for model in keys(modelnames)
		df = scalars[model]
		lines!(df.time * T2GYR, df.sigma_v * V2KMS, label=model)
	end
	vlines!(0, color=:black, linestyle=:dot)

	hlines!(7.5, color=:black,)
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═3a8954e9-0f5b-4f2c-8a9e-6e66f0f20ccc
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═e93918f6-0c72-4dec-9cdf-97f185c0bceb
# ╠═f38390b2-b6dd-4c42-a548-fab9b6ac2a73
# ╠═56365cfd-16e0-4567-a85c-f59372a611c7
# ╠═437ac44b-5440-4b3c-aa58-4bbd93be7174
# ╠═fc0a118f-f8b7-425e-aacf-c8c7630d61a5
# ╠═8c2d254f-7443-4ff1-ae5e-55b3bd0df766
# ╠═f1443328-05ae-41be-8242-d92ed430908b
# ╠═080c9586-bb3d-4f80-b4ea-5bae53f174f4
# ╠═d26db585-3e68-4870-99da-bb5d9182c923
# ╠═dbed5ec4-6030-4667-a87b-e8100deebe4a
# ╠═9e8271b9-173a-47d6-a144-7dbaeb3badc7
# ╠═1f026ba8-0c58-4cf1-8b4e-c4344f7c8c72
# ╠═d3c09dc2-db77-471c-9384-3ffeb0df39dc
# ╠═713590cd-0446-4096-b3b8-c0b8df850bea
# ╠═604a69dd-1768-4156-9e65-9e4d3f8a120e
# ╠═898d0210-1da9-4dfd-84b8-114bd2511e62
# ╠═e66a81b1-36ba-48ca-be63-1658ec3804f6
