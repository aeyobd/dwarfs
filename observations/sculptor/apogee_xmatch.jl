### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 04bbc735-e0b4-4f0a-9a83-e50c8b923caf
begin 
	import Pkg; Pkg.activate()
	
	using Arya
	using CairoMakie
	
	import LilGuys as lguys
end

# ╔═╡ 399bc588-515c-4a5f-88c4-cae5272d20a2
using PyFITS

# ╔═╡ 08fc8a65-fa92-4739-879e-a32d0228d181
using DataFrames: innerjoin, rename!

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Radial Velocity CrossMatch

This notebook loads in the complete APOGEE radial velocity data and crossmatches with the J+24.

With this catalogue, detailed radial velocity analysis can then be calculated.
"""

# ╔═╡ 0efc36d3-7958-4873-a15a-244778e94c61
import DataFrames: disallowmissing

# ╔═╡ ab07e08f-0d80-409a-98cd-35cedf341ba3
CairoMakie.activate!(type=:png)

# ╔═╡ 4e6691d3-42f2-4929-9351-60be5d3235d5
md"""
## load datas
"""

# ╔═╡ f26e6eb1-640b-40fa-8cbc-e87092f132f7
j24 = read_fits("data/jensen+24_wide.fits")

# ╔═╡ 5a9c76ed-eccf-42d7-96f1-814471f66d89
 apogee_f = read_fits("../all/data/allStarLite-dr17-synspec_rev1.fits")

# ╔═╡ 0587bd4e-459e-43f6-b0b3-ef39c393d55b
md"""
## Proces datums
"""

# ╔═╡ a1dd09f8-b635-48e5-8d0c-afd58dff4641
apogee = rename!(apogee_f, 
	:GAIAEDR3_SOURCE_ID => :source_id,
	:VHELIO_AVG => :RV,
	:VERR => :RV_err,
	:VSCATTER => :RV_sigma,
	:NVISITS => :RV_count,
	:RA => :RA_apogee,
	:DEC => :DEC_apogee,
)

# ╔═╡ bc8791ad-dfc3-49c2-9849-3d7bcf448978
joined = innerjoin(j24, apogee, on=:source_id)

# ╔═╡ a620300d-1078-4c4c-9803-34f8ef10ec4d
good = joined[.!ismissing.(joined.RV), :]

# ╔═╡ ff98fe21-fc5a-4bb8-bb0e-3cd2f4725fa3
write_fits( "data/apogee_xmatch.fits", good, overwrite=true)

# ╔═╡ 62c90fb6-1310-47d7-a0ca-b8505f9ad129
md"""
## Plots
"""

# ╔═╡ fdb6e310-070c-41e7-b81d-b239e85c961c
hist(filter(x -> !ismissing(x) && x<30, apogee.GAIAEDR3_PHOT_G_MEAN_MAG))

# ╔═╡ 209031d1-9305-40a2-817a-8eb5a904af76
let
	fig = Figure()
    ax = Axis(fig[1,1])

	scatter!(joined.ra, joined.dec, markersize=5)
	scatter!(disallowmissing(joined.RA_apogee), disallowmissing(joined.DEC_apogee), markersize=3)

	fig
end

# ╔═╡ 5488c919-4a20-4918-af26-647e0f1ca716
hist(good.RV)

# ╔═╡ f7901e49-8ac6-40ea-8128-77785e72dbf7
hist(filter(x -> !ismissing(x) && x<30, apogee.K))

# ╔═╡ Cell order:
# ╠═811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═399bc588-515c-4a5f-88c4-cae5272d20a2
# ╠═08fc8a65-fa92-4739-879e-a32d0228d181
# ╠═0efc36d3-7958-4873-a15a-244778e94c61
# ╠═ab07e08f-0d80-409a-98cd-35cedf341ba3
# ╟─4e6691d3-42f2-4929-9351-60be5d3235d5
# ╠═f26e6eb1-640b-40fa-8cbc-e87092f132f7
# ╠═5a9c76ed-eccf-42d7-96f1-814471f66d89
# ╟─0587bd4e-459e-43f6-b0b3-ef39c393d55b
# ╠═a1dd09f8-b635-48e5-8d0c-afd58dff4641
# ╠═bc8791ad-dfc3-49c2-9849-3d7bcf448978
# ╠═a620300d-1078-4c4c-9803-34f8ef10ec4d
# ╠═ff98fe21-fc5a-4bb8-bb0e-3cd2f4725fa3
# ╠═62c90fb6-1310-47d7-a0ca-b8505f9ad129
# ╠═fdb6e310-070c-41e7-b81d-b239e85c961c
# ╠═209031d1-9305-40a2-817a-8eb5a904af76
# ╠═5488c919-4a20-4918-af26-647e0f1ca716
# ╠═f7901e49-8ac6-40ea-8128-77785e72dbf7
