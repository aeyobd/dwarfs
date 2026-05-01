### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 3a02cce0-3430-11f1-abe3-ef9766dea525
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya, CairoMakie
	
	using PyFITS
end

# ╔═╡ 84d23982-3393-4319-95b8-ded7e66be1c8
using OrderedCollections

# ╔═╡ f912dcd3-51fe-456f-8189-b960dda0cdc8
import TOML

# ╔═╡ 84e4e3a4-63c0-4567-869a-8aa889913a6f
CairoMakie.activate!(type=:png)

# ╔═╡ b16da089-2324-4dac-8121-353013fb9227
md"""
# Data loading
"""

# ╔═╡ 4f5c53c0-3552-4ab3-8bc9-86042b738622
model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", "bootes3")

# ╔═╡ a2daf333-f3ab-45bf-9160-a37d27bcf0b7
function get_time_f(name)
	props = TOML.parsefile(joinpath(model_dir, name, "orbital_properties.toml"))

	props["t_f_gyr"] / T2GYR
end


# ╔═╡ d8e2a511-8209-4a14-9692-0cd90d714d76
function read_scalars(name)
	df= read_fits(joinpath(model_dir, name, "profiles_scalars.fits"))
	df.time .-= get_time_f(name)
	df
end

# ╔═╡ 7af7cf5a-d149-4fe0-9566-d23b293af3ec
scalars = OrderedDict(
	"compact: 1x1.5kpc" => read_scalars("1e6_v30_r3.0/1_peri_1.5kpc"),
	"compact: 2x7kpc" => read_scalars("1e6_v30_r3.0/2_peri_7kpc"),
	"compact: 5x18kpc" => read_scalars("1e6_v30_r3.0/5_peri_18kpc"),
	"mean: 1x12kpc" => read_scalars("1e6_v22_r3.9/1_peri_12kpc"),
	"mean: 3x26kpc" => read_scalars("1e6_v22_r3.9/3_peri_26kpc"),

)

# ╔═╡ 43b5f031-08ac-45ae-a00b-5cdb375f979e
scalars_cored = OrderedDict(
	"NFW" => read_scalars("1e6_v30_r3.0/5_peri_18kpc"),
	"0.1 core" => read_scalars("1e5_v30_r3.0_c0.1/5_peri_18kpc"),
	"0.1 core-compact" => read_scalars("1e5_v30_r2.2_c0.1/5_peri_18kpc"),
	"1.0 core" => read_scalars("1e5_v30_r2.2_c1/5_peri_18kpc"),

)

# ╔═╡ 62b21b67-f511-473b-bae8-9cdf3c77e03a
md"""
# Plot
"""

# ╔═╡ 76789b1e-022a-425c-ab8f-0e42f0e422f2
function plot_rmax_vmax(scalars)
	fig = Figure(size=(6, 3) .* 72)
	
	ax = Axis(fig[1,1],
			 xlabel = L"r_\text{max}\ / \ r_\text{max, 0}", 
			 ylabel = L"\text{v}_\text{max}\ / \ \text{v}_\text{max, 0}",
			 xscale=log10,
			 yscale=log10)

	for (label,df) in scalars
		filt = df.time .< 0
		lines!(df.r_circ_max[filt] ./ df.r_circ_max[1], df.v_circ_max[filt] / df.v_circ_max[1], label=label)
	end

	x, y = LilGuys.EN21_tidal_track(1, 1, x_min=0.2)
	lines!(x, y, color=:black)
	

	ax = Axis(fig[1,2],
			 xlabel = "time / Gyr", 
			 ylabel = L"\text{v}_\text{max}\ /\ \text{km}\,\text{s}^{-1}",
			 )

	for (label,df) in scalars
		lines!(df.time * T2GYR, df.v_circ_max*V2KMS, label=label, alpha=1)
	end

	axislegend(position=:lb)
	fig

end

# ╔═╡ c98463d0-b895-473f-95c6-64e237334801
plot_rmax_vmax(scalars)

# ╔═╡ ed0ebba7-5ec8-43ad-a259-419342ac00e8
plot_rmax_vmax(scalars_cored)

# ╔═╡ Cell order:
# ╠═3a02cce0-3430-11f1-abe3-ef9766dea525
# ╠═f912dcd3-51fe-456f-8189-b960dda0cdc8
# ╠═84d23982-3393-4319-95b8-ded7e66be1c8
# ╠═84e4e3a4-63c0-4567-869a-8aa889913a6f
# ╠═b16da089-2324-4dac-8121-353013fb9227
# ╠═4f5c53c0-3552-4ab3-8bc9-86042b738622
# ╠═a2daf333-f3ab-45bf-9160-a37d27bcf0b7
# ╠═d8e2a511-8209-4a14-9692-0cd90d714d76
# ╠═7af7cf5a-d149-4fe0-9566-d23b293af3ec
# ╠═43b5f031-08ac-45ae-a00b-5cdb375f979e
# ╠═62b21b67-f511-473b-bae8-9cdf3c77e03a
# ╠═76789b1e-022a-425c-ab8f-0e42f0e422f2
# ╠═c98463d0-b895-473f-95c6-64e237334801
# ╠═ed0ebba7-5ec8-43ad-a259-419342ac00e8
