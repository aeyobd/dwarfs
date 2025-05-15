### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 4248b700-31d9-11f0-2ec1-9b80c5cbb7d2
begin 
	import Pkg; Pkg.activate()
	using Arya, CairoMakie, PyFITS
end

# ╔═╡ 5c3a6faa-fc01-44cf-a0fc-6969df95aef1
using LilGuys

# ╔═╡ bc5b5c00-ac45-4ea4-ade1-6b2130ff1cf5
df = read_fits(joinpath(ENV["DWARFS_ROOT"], "observations/ursa_minor/data/jensen+24_2c.fits"))

# ╔═╡ 8b512217-a471-4ca5-99b8-e33632adc71e
nextfloat(1.0) ==1.0f0

# ╔═╡ 434582b0-cddc-4c33-b44d-6117fb3b3cca
CairoMakie.activate!(type=:png)

# ╔═╡ c7e2f35a-20b6-4733-82b5-233e0ca79ff7
@savefig "figures/umi_gaia_memb" let
	fig = Figure(
		backgroundcolor=:transparent
	)
	ax = Axis(fig[1,1],
			 xreversed=true,
			 aspect=DataAspect(),
			 limits = (-1., 1., -1., 1.) ,
			  backgroundcolor = :transparent
			 )

	hidespines!()
	hidexdecorations!()
	hideydecorations!()

	filt = (df.PSAT .> 0.2) .& (df.F_BEST .== 1.)
	#filt .&= .!ismissing.(filt)
	scatter!(df.xi[filt], df.eta[filt], markersize=1, color=COLORS[2])

	resize_to_layout!(fig)
	fig
end

# ╔═╡ d07e26e0-db92-426e-9c8f-a0c2c3a6b247
@savefig "figures/umi_gaia_all" let
	fig = Figure(
		backgroundcolor=:transparent
	)
	ax = Axis(fig[1,1],
			 xreversed=true,
			 aspect=DataAspect(),
			 limits = (-1., 1., -1., 1.),
			  backgroundcolor = :transparent
			 )

	hidespines!()
	hidexdecorations!()
	hideydecorations!()

	filt = .!ismissing.(df.phot_g_mean_mag) .& (df.phot_g_mean_mag .< 22)
	scatter!(df.xi[filt], df.eta[filt], markersize=200 ./ df.phot_g_mean_mag[filt] .^ 2, color=COLORS[1])

	resize_to_layout!(fig)
	fig
end

# ╔═╡ Cell order:
# ╠═4248b700-31d9-11f0-2ec1-9b80c5cbb7d2
# ╠═bc5b5c00-ac45-4ea4-ade1-6b2130ff1cf5
# ╠═8b512217-a471-4ca5-99b8-e33632adc71e
# ╠═5c3a6faa-fc01-44cf-a0fc-6969df95aef1
# ╠═434582b0-cddc-4c33-b44d-6117fb3b3cca
# ╠═c7e2f35a-20b6-4733-82b5-233e0ca79ff7
# ╠═d07e26e0-db92-426e-9c8f-a0c2c3a6b247
