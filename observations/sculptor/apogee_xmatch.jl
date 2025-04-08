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
	import PythonCall # fits
end

# ╔═╡ 08fc8a65-fa92-4739-879e-a32d0228d181
using DataFrames: innerjoin, rename!

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Radial Velocity CrossMatch

This notebook loads in the complete APOGEE radial velocity data and crossmatches with the J+24.

With this catalogue, detailed radial velocity analysis can then be calculated.
"""

# ╔═╡ f26e6eb1-640b-40fa-8cbc-e87092f132f7
j24 = lguys.read_fits("data/jensen+24_wide.fits")

# ╔═╡ 5a9c76ed-eccf-42d7-96f1-814471f66d89
 apogee_f = lguys.read_fits("../all/data/allStarLite-dr17-synspec_rev1.fits")

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
good = joined[.!isnan.(joined.RV), :]

# ╔═╡ ab07e08f-0d80-409a-98cd-35cedf341ba3
CairoMakie.activate!(type=:png)

# ╔═╡ 209031d1-9305-40a2-817a-8eb5a904af76
let
	fig = Figure()
    ax = Axis(fig[1,1])

	scatter!(joined.ra, joined.dec, markersize=5)
	scatter!(joined.RA_apogee, joined.DEC_apogee, markersize=3)

	fig
end

# ╔═╡ 5488c919-4a20-4918-af26-647e0f1ca716
hist(good.RV)

# ╔═╡ 8375c521-5834-48f5-acaf-51c28dd4c5d4


# ╔═╡ ff98fe21-fc5a-4bb8-bb0e-3cd2f4725fa3
lguys.write_fits( "data/apogee_xmatch.fits", good, overwrite=true)

# ╔═╡ Cell order:
# ╠═811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═f26e6eb1-640b-40fa-8cbc-e87092f132f7
# ╠═5a9c76ed-eccf-42d7-96f1-814471f66d89
# ╠═a1dd09f8-b635-48e5-8d0c-afd58dff4641
# ╠═08fc8a65-fa92-4739-879e-a32d0228d181
# ╠═bc8791ad-dfc3-49c2-9849-3d7bcf448978
# ╠═a620300d-1078-4c4c-9803-34f8ef10ec4d
# ╠═ab07e08f-0d80-409a-98cd-35cedf341ba3
# ╠═209031d1-9305-40a2-817a-8eb5a904af76
# ╠═5488c919-4a20-4918-af26-647e0f1ca716
# ╠═8375c521-5834-48f5-acaf-51c28dd4c5d4
# ╠═ff98fe21-fc5a-4bb8-bb0e-3cd2f4725fa3
