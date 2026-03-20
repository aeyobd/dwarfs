### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 0eb2456e-23ac-11f1-853b-af81307f0fdc
begin
	import Pkg;
	Pkg.activate()

	using PyFITS
	using Arya
	using CairoMakie
end

# ╔═╡ 1010ae2c-9d31-4388-ab0d-ecd24dc3eee2
using PythonCall

# ╔═╡ e079279a-7e34-4ce4-85de-b0c15f7ae55e


# ╔═╡ 35665422-b330-4c36-ad4b-a2ceca5f8abb
all_stars = read_fits("data/delve_dr2.fits")

# ╔═╡ dcc11561-1d6c-4647-9b71-d1e4de2287b1
md"""
PS depth from DELVE DR 2 data release paper (5sigma PSF mags) (Drlica-Wagner et al. 2022).
"""

# ╔═╡ c34b7bb3-bc7c-4e3c-ad26-83ef4351bb27
begin
	ps_depth_g = 24.3
	ps_depth_r = 23.9
	ps_depth_i = 23.5
	ps_depth_z = 22.8
end

# ╔═╡ 8606a03d-a87c-425d-a811-9e5a33089403
filt_good = (all_stars.mag_psf_g .< 99) .&
	(all_stars.magerr_psf_g .< 99) .&
	(all_stars.mag_psf_r .< 99) .&
	(all_stars.magerr_psf_r .< 99) .&
	(all_stars.extended_class_g .<= 1) .& 
	(all_stars.mag_psf_g .< ps_depth_g) .&
	(all_stars.mag_psf_r .< ps_depth_r)
	 

# ╔═╡ 178978d2-ab0b-4967-8284-735dbd5d92cc
sum(filt_good)

# ╔═╡ bcd2ddf5-e48f-4833-9fc8-d7be28f2b1c7
hist(all_stars.extended_class_g)

# ╔═╡ ebfd855c-23fe-4009-8d3e-1450b634d6e1
good_stars = all_stars[filt_good, :]

# ╔═╡ bf77a2fc-716c-4b92-bccc-ebdda4d148f9
hist(good_stars.mag_psf_g)

# ╔═╡ d7af8873-3744-4717-b9e4-2e5cbeb270ac
scatter(good_stars.mag_psf_g, good_stars.magerr_psf_g)

# ╔═╡ 01f5645a-651f-4e28-9621-59589538a0cf
scatter(good_stars.mag_psf_r, good_stars.magerr_psf_r)

# ╔═╡ b8411456-53e0-4916-9145-af5d1b9e1069
scatter(good_stars.mag_psf_g .- good_stars.mag_psf_r, good_stars.mag_psf_g)

# ╔═╡ 944893ba-7d69-48b9-8d02-2c297429dfed
Dustmaps = pyimport("dustmaps.sfd")

# ╔═╡ db00cdb7-a5ab-4548-a315-bcbef5b52278
DUSTMAP = Dustmaps.SFDQuery()

# ╔═╡ 761d8cdc-3aea-4ba6-b91a-ee74fb9032ef
SkyCoord = pyimport("astropy.coordinates").SkyCoord

# ╔═╡ 7f011860-c221-40ea-a25d-420e0524dc85
function get_extinction(ra, dec)
    icrss = SkyCoord(ra, dec, unit="degree", frame="icrs")

    ebv = pyconvert(Vector, DUSTMAP(icrss))

    A0=3.1*ebv

    Ag, Ar, Ai = let
		# from abbot et al. 2018 section 4.1 (DES DR 1)
		Rg=3.186
		Rr=2.140
		Ri=1.569
		Rz=1.196
		RY=1.048

		return Rg*ebv, Rr*ebv, Ri * ebv
    end

    return Ag, Ar, Ai
end


# ╔═╡ 3c6572f2-e135-438e-858e-1e537743c919
Ag, Ar, Ai = get_extinction(good_stars.ra, good_stars.dec)

# ╔═╡ 4906a26a-7968-46e8-a9d8-9d6c2ac3f816
df_out = let
	df = copy(good_stars)
	df[!, :gmag] = df.mag_psf_g .- Ag
	df[!, :rmag] = df.mag_psf_r .- Ar
	df[!, :imag] = df.mag_psf_i .- Ai

	df[!, :A_g] = Ag
	df[!, :A_g] = Ar
	df[!, :A_g] = Ai

	df
end

# ╔═╡ d4a91259-f9a1-4832-ad5c-04aa8855bdae
write_fits("data/delve_dr2_good.fits", df_out, overwrite=true)

# ╔═╡ 1928c2c0-6720-4060-85e4-c5002f9dc525
let
	fig = Figure()
	ax = Axis(fig[1,1])
		
	p = scatter!(good_stars.ra, good_stars.dec, color=Ag .- Ar, markersize=1)

	Colorbar(fig[1,2], p)
	fig
end

# ╔═╡ Cell order:
# ╠═0eb2456e-23ac-11f1-853b-af81307f0fdc
# ╠═e079279a-7e34-4ce4-85de-b0c15f7ae55e
# ╠═35665422-b330-4c36-ad4b-a2ceca5f8abb
# ╠═178978d2-ab0b-4967-8284-735dbd5d92cc
# ╠═8606a03d-a87c-425d-a811-9e5a33089403
# ╠═dcc11561-1d6c-4647-9b71-d1e4de2287b1
# ╠═c34b7bb3-bc7c-4e3c-ad26-83ef4351bb27
# ╠═bcd2ddf5-e48f-4833-9fc8-d7be28f2b1c7
# ╠═ebfd855c-23fe-4009-8d3e-1450b634d6e1
# ╠═bf77a2fc-716c-4b92-bccc-ebdda4d148f9
# ╠═d7af8873-3744-4717-b9e4-2e5cbeb270ac
# ╠═01f5645a-651f-4e28-9621-59589538a0cf
# ╠═b8411456-53e0-4916-9145-af5d1b9e1069
# ╠═1010ae2c-9d31-4388-ab0d-ecd24dc3eee2
# ╠═944893ba-7d69-48b9-8d02-2c297429dfed
# ╠═db00cdb7-a5ab-4548-a315-bcbef5b52278
# ╠═761d8cdc-3aea-4ba6-b91a-ee74fb9032ef
# ╠═3c6572f2-e135-438e-858e-1e537743c919
# ╠═4906a26a-7968-46e8-a9d8-9d6c2ac3f816
# ╠═d4a91259-f9a1-4832-ad5c-04aa8855bdae
# ╠═1928c2c0-6720-4060-85e4-c5002f9dc525
# ╠═7f011860-c221-40ea-a25d-420e0524dc85
