### A Pluto.jl notebook ###
# v0.20.8

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

# ╔═╡ 1142d331-52c9-43f2-95f4-9983dae8aee0
using DataFrames

# ╔═╡ 811c5da0-7e70-4393-b59d-c0fdb89523ca
md"""
# Radial Velocity CrossMatch

This notebook loads in the complete APOGEE radial velocity data and crossmatches with the J+24.

With this catalogue, detailed radial velocity analysis can then be calculated.
"""

# ╔═╡ 72fe0dab-a127-4dd1-9bb5-d3868edccf59
md"""
Creates:
- `processed/apogee_xmatch.fits`
- `processed/apogee_allvisit_xmatch.fits`

Depends on:
- `../../all/data/allVisit-dr17-synspec_rev1.fits`
- `../../all/data/allStarLite-dr17-synspec_rev1.fits`
- `../data/jensen+24_wide.fits`
"""

# ╔═╡ 0efc36d3-7958-4873-a15a-244778e94c61
import DataFrames: disallowmissing, leftjoin, groupby

# ╔═╡ ab07e08f-0d80-409a-98cd-35cedf341ba3
CairoMakie.activate!(type=:png)

# ╔═╡ 4e6691d3-42f2-4929-9351-60be5d3235d5
md"""
## load datas
"""

# ╔═╡ f26e6eb1-640b-40fa-8cbc-e87092f132f7
j24 = read_fits("../data/jensen+24_wide_2c.fits")

# ╔═╡ 5a9c76ed-eccf-42d7-96f1-814471f66d89
 apogee_f = read_fits("../../all/data/allStarLite-dr17-synspec_rev1.fits")

# ╔═╡ 0587bd4e-459e-43f6-b0b3-ef39c393d55b
md"""
## Proces datums
"""

# ╔═╡ a1dd09f8-b635-48e5-8d0c-afd58dff4641
apogee = rename(apogee_f, 
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
write_fits( "processed/apogee_xmatch.fits", good, overwrite=true)

# ╔═╡ 80cf2f44-1e4d-4cf9-8ee2-0ac13724a9ac
md"""
# Visit level
"""

# ╔═╡ 3821df61-9609-4dc1-bdc9-56765cded097
md"""
Double check how visit recombination works and if uncertanties make sense.
"""

# ╔═╡ c89af488-7800-4b9e-ae2f-9faa5823e907
 apogee_allvisit = read_fits("../../all/data/allVisit-dr17-synspec_rev1.fits")

# ╔═╡ 96ae8bc4-83ac-4497-9a12-94b41a4bebab
allvisit = leftjoin(good, apogee_allvisit, on=:APOGEE_ID, makeunique=true) |> x-> rename!(x, 
	:RV => :RV_mean,
	:RV_err => :RV_mean_err,
	:VRELERR => :RV_err,
	:VHELIO => :RV)

# ╔═╡ a7674294-6841-4e71-9da7-d14537e434a9
write_fits("processed/apogee_allvisit_xmatch.fits", allvisit, overwrite=true)

# ╔═╡ dc83a043-8e92-43fd-8d1b-a5504f6649d5
starflag = 1<<0 + 1<<3 + 1<<4 #+ 1<<19 + 1<<22

# ╔═╡ 27e57506-8e50-4a10-b295-d1704dbe66de
goodspec = (allvisit.STARFLAG_1 .& starflag) .== 0

# ╔═╡ b4fd6a65-b1bc-40fd-b5cc-735d030cd702
allvisit.STARFLAG_1

# ╔═╡ 682da3a3-8e75-4f20-b99f-4e35935570be
sum( goodspec)

# ╔═╡ 42bea921-21d6-41fb-8236-9ffca2dcfd79
⊕(a, b) = sqrt(a^2 + b^2)

# ╔═╡ 98c9d2f1-eda2-4ab6-a9d0-10d39944fc83
Δv = @. abs(allvisit.RV .- allvisit.RV_mean) / (allvisit.RV_mean_err ⊕ allvisit.RV_err)

# ╔═╡ c959db9d-02c7-416d-ad22-8ae639d3dec9
Δv2 = @. abs(allvisit.RV_mean .- allvisit.RV) / allvisit.RV_sigma

# ╔═╡ dcaf2cd3-34bd-44fa-9c0a-ec29dd4aa8dc
import Distributions: cdf, Normal, Chisq

# ╔═╡ 84eb7182-f6e3-4d20-99a9-9065eb08804f
p = @. 1 - cdf(Chisq(2), Δv[goodspec])

# ╔═╡ d6c46bae-238d-4cdb-b45c-c79db84ef0a0
p2 = @. 1 - cdf(Chisq(1), Δv2[goodspec])

# ╔═╡ cb974375-e0ed-4f2f-a919-1c1616ee2391
minimum(skipmissing(Δv))

# ╔═╡ 75a77f2c-16e1-4ab3-a41b-4d8b8e8a71d6
hist(p)

# ╔═╡ 1fbd52ab-bbca-482a-80dc-35733fb534ed
hist(filter(x->!ismissing(x) & isfinite(x), allvisit.RV_err))

# ╔═╡ 095ae700-774b-4ec8-a47b-abb89a16702e
RV_mean = [lguys.mean(allvisit.RV[allvisit.source_id .== source_id],) for source_id in good.source_id]

# ╔═╡ c340afca-9bd2-44dc-898a-cf721a5af909
RV_err = [1 / sqrt(sum(allvisit.RV_err[allvisit.source_id .== source_id] .^ -2)) for source_id in good.source_id]

# ╔═╡ f8d205be-f298-40a6-ace1-858d6c5f1752
RV_sigma =  [lguys.std(allvisit.RV[allvisit.source_id .== source_id]) for source_id in good.source_id]

# ╔═╡ 768ec557-f305-48f7-86f7-fa02362b1c68
allvisit.RV[allvisit.source_id .== good.source_id[1]]

# ╔═╡ 524b33c4-f435-45fd-90e5-9277333268a1
scatter(RV_mean, good.RV .- RV_mean)

# ╔═╡ 871afe21-b267-40fd-98d5-3ca8abc2bea2
scatter(RV_err, good.RV_err ./ RV_err ./ sqrt.(good.RV_count))

# ╔═╡ 1f79d6e7-0575-45e9-93c0-9e0fbb2034ee
scatter(RV_sigma, good.RV_sigma ./ RV_sigma)

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
# ╟─811c5da0-7e70-4393-b59d-c0fdb89523ca
# ╟─72fe0dab-a127-4dd1-9bb5-d3868edccf59
# ╠═04bbc735-e0b4-4f0a-9a83-e50c8b923caf
# ╠═399bc588-515c-4a5f-88c4-cae5272d20a2
# ╠═08fc8a65-fa92-4739-879e-a32d0228d181
# ╠═0efc36d3-7958-4873-a15a-244778e94c61
# ╠═1142d331-52c9-43f2-95f4-9983dae8aee0
# ╠═ab07e08f-0d80-409a-98cd-35cedf341ba3
# ╟─4e6691d3-42f2-4929-9351-60be5d3235d5
# ╠═f26e6eb1-640b-40fa-8cbc-e87092f132f7
# ╠═5a9c76ed-eccf-42d7-96f1-814471f66d89
# ╟─0587bd4e-459e-43f6-b0b3-ef39c393d55b
# ╠═a1dd09f8-b635-48e5-8d0c-afd58dff4641
# ╠═bc8791ad-dfc3-49c2-9849-3d7bcf448978
# ╠═a620300d-1078-4c4c-9803-34f8ef10ec4d
# ╠═ff98fe21-fc5a-4bb8-bb0e-3cd2f4725fa3
# ╟─80cf2f44-1e4d-4cf9-8ee2-0ac13724a9ac
# ╟─3821df61-9609-4dc1-bdc9-56765cded097
# ╠═c89af488-7800-4b9e-ae2f-9faa5823e907
# ╠═96ae8bc4-83ac-4497-9a12-94b41a4bebab
# ╠═a7674294-6841-4e71-9da7-d14537e434a9
# ╠═dc83a043-8e92-43fd-8d1b-a5504f6649d5
# ╠═27e57506-8e50-4a10-b295-d1704dbe66de
# ╠═b4fd6a65-b1bc-40fd-b5cc-735d030cd702
# ╠═682da3a3-8e75-4f20-b99f-4e35935570be
# ╠═42bea921-21d6-41fb-8236-9ffca2dcfd79
# ╠═98c9d2f1-eda2-4ab6-a9d0-10d39944fc83
# ╠═c959db9d-02c7-416d-ad22-8ae639d3dec9
# ╠═dcaf2cd3-34bd-44fa-9c0a-ec29dd4aa8dc
# ╠═84eb7182-f6e3-4d20-99a9-9065eb08804f
# ╠═d6c46bae-238d-4cdb-b45c-c79db84ef0a0
# ╠═cb974375-e0ed-4f2f-a919-1c1616ee2391
# ╠═75a77f2c-16e1-4ab3-a41b-4d8b8e8a71d6
# ╠═1fbd52ab-bbca-482a-80dc-35733fb534ed
# ╠═095ae700-774b-4ec8-a47b-abb89a16702e
# ╠═c340afca-9bd2-44dc-898a-cf721a5af909
# ╠═f8d205be-f298-40a6-ace1-858d6c5f1752
# ╠═768ec557-f305-48f7-86f7-fa02362b1c68
# ╠═524b33c4-f435-45fd-90e5-9277333268a1
# ╠═871afe21-b267-40fd-98d5-3ca8abc2bea2
# ╠═1f79d6e7-0575-45e9-93c0-9e0fbb2034ee
# ╟─62c90fb6-1310-47d7-a0ca-b8505f9ad129
# ╠═fdb6e310-070c-41e7-b81d-b239e85c961c
# ╠═209031d1-9305-40a2-817a-8eb5a904af76
# ╠═5488c919-4a20-4918-af26-647e0f1ca716
# ╠═f7901e49-8ac6-40ea-8128-77785e72dbf7
