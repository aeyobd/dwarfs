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
using DataFrames: innerjoin

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Radial Velocity CrossMatch

This notebook loads in the complete APOGEE radial velocity data and crossmatches with the J+24.

With this catalogue, detailed radial velocity analysis can then be calculated.
"""

# ╔═╡ f26e6eb1-640b-40fa-8cbc-e87092f132f7
j24 = lguys.read_fits("data/jensen+24_2c.fits")

# ╔═╡ 5a9c76ed-eccf-42d7-96f1-814471f66d89
 apogee_f = lguys.read_fits("../all/data/allStarLite-dr17-synspec_rev1.fits", columns=["GAIAEDR3_SOURCE_ID", "VHELIO_AVG", "VERR", "VSCATTER", "NVISITS", "RV_FLAG", "RA", "DEC"])

# ╔═╡ a1dd09f8-b635-48e5-8d0c-afd58dff4641
apogee = lguys.DataFrame(
	:source_id => apogee_f.GAIAEDR3_SOURCE_ID,
	:RV => apogee_f.VHELIO_AVG,
	:RV_err => apogee_f.VERR,
	:RV_sigma => apogee_f.VSCATTER,
	:RV_count => apogee_f.NVISITS,
	:RV_flag => apogee_f.RV_FLAG,
	:RA_apogee => apogee_f.RA,
	:DEC_apogee => apogee_f.DEC,
)

# ╔═╡ bc8791ad-dfc3-49c2-9849-3d7bcf448978
joined = innerjoin(j24, apogee, on=:source_id)

# ╔═╡ a620300d-1078-4c4c-9803-34f8ef10ec4d
good = joined[.!isnan.(joined.RV), :]

# ╔═╡ adbbadd8-a208-4001-979d-1b1c279e89a7
scatter(joined.ra, joined.RA_apogee)

# ╔═╡ 9479f3e0-4b35-4c1e-ab13-0ab0af19728f
scatter(joined.dec, joined.DEC_apogee)

# ╔═╡ 5488c919-4a20-4918-af26-647e0f1ca716
hist(good.RV)

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
# ╠═adbbadd8-a208-4001-979d-1b1c279e89a7
# ╠═9479f3e0-4b35-4c1e-ab13-0ab0af19728f
# ╠═5488c919-4a20-4918-af26-647e0f1ca716
# ╠═ff98fe21-fc5a-4bb8-bb0e-3cd2f4725fa3
