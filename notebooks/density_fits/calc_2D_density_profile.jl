### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 142a5ace-1432-4093-bee7-4a85c19b0d72
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames, CSV
	
	using CairoMakie
			
	import LilGuys as lguys
	using Arya
	end

# ╔═╡ 852717c0-aabf-4c03-9cf5-a6d91174e0f9
md"""
given a sample of points, can we centre and calculate the 2D density profile

### Inputs
- Number of bins
- Ellipticity
- Centring method
"""

# ╔═╡ 4cc4e2be-6bf6-4cbd-a2b1-121354a862bc
simulation = false

# ╔═╡ c4574ed3-f431-4c0e-a721-1c8d88dda10f
name = "exp2d_rs0.1_today"

# ╔═╡ bb644db4-7fb4-43c8-abf9-7235aa279ad6
r_centre = 3

# ╔═╡ 73f0b3a1-a4b6-422d-9f7e-be816c4a9cfc
begin 
	samplename = "/astro/dboyea/sculptor/orbit1/stars/$name.fits" # = "$(name)_sample.fits" 
	# samplename = "sculptor/fiducial_sample.fits" # = "$(name)_sample.fits" 
	#samplename = "../test_sky_recon.fits"
	if !simulation
		samplename = "sculptor/fiducial_sample.fits"
	end
end

# ╔═╡ a2465c61-ce25-42aa-8b5c-57ad7ffe16f6
outname = dirname(samplename) * "/$(name)_profile.toml"

# ╔═╡ a82d8fa5-32db-42d1-8b0a-d54ae47dc7be
begin 
	if simulation
		ell = 0
	else
		ell = 0.37
	end
	PA = 91
	centre_method="mean"
	if simulation
		mass_column = :probability
	else
		mass_column = nothing
	end
	normalize = true
end

# ╔═╡ 72d975fe-9d97-474b-ba3f-f61ba12c7c80
begin 
	f = FITS(samplename, "r")
	sample = DataFrame(f[2])
	close(f)
end

# ╔═╡ 305f79e0-a9bb-4c15-a9fd-09cb1c25db41
import StatsBase: sem, mean

# ╔═╡ 77e1461f-26c3-4b3a-8d16-a9dc095eb572
sem(sample.ra)

# ╔═╡ 3ade8bb5-8698-4d10-862c-e9685e54e570
sem(sample.dec)

# ╔═╡ 568162b2-ad89-4c72-92e4-ee3cc302bbb7
mean(sample.pmra), sem(sample.pmra)

# ╔═╡ 0272086f-ebeb-432b-b44e-7c35a17c1d58
mean(sample.pmdec), sem(sample.pmdec)

# ╔═╡ 8d276372-add5-4388-b713-b22e38d56f37
if mass_column === nothing
	weights = ones(size(sample, 1))
else
	weights = sample[:, mass_column]
end

# ╔═╡ cb38c9f9-d6ff-4bcd-a819-9b442776ccfc
ra0, dec0 = lguys.calc_centre2D(sample.ra, sample.dec, centre_method, weights )

# ╔═╡ 86dd90bc-83dd-4b1a-8e98-1bb0333c6610
xi, eta = lguys.to_tangent(sample.ra, sample.dec, ra0, dec0)

# ╔═╡ ef19dcd1-fae0-4777-a0d8-d242435f892f
r_max = sqrt(maximum(xi .^ 2 .+ eta .^ 2))

# ╔═╡ 69018984-ef00-44ef-ba6e-7cccf930aef9
let
	global r_ell
	b = sqrt(1-ell)
	a = 1/b
	
	r_ell = lguys.calc_r_ell(xi, eta, a, b, PA-90)

	r_ell_max = r_max .* b
	r_ell = r_ell[r_ell .< r_ell_max]


	r_ell .*= 60 # to arcmin
end

# ╔═╡ f476c859-ba4b-4343-8184-f6f41dc092ee
profile = lguys.calc_properties(r_ell, bins=50, weights=weights, normalization=:central, r_centre=r_centre)

# ╔═╡ b87cf54b-2c4c-46c5-9336-c13e773e29ec
begin 

	open(outname, "w") do f
		print(f, profile)
	end

	println("wrote data to ", outname)
end

# ╔═╡ 0b80353f-c698-4fae-b40f-1796f7c89792
sum(profile.counts), size(sample)

# ╔═╡ e63df7b9-1fc1-4cc0-a91f-c0f1395d7ff4
md"""
# Plots
"""

# ╔═╡ d382b455-73ba-41d5-bbc8-eca53ea2166e
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		limits = (ra0 .+ (-1, 1) ./ cosd(dec0), dec0 .+ (-1, 1)),
		aspect=1,
		xlabel="ra / degrees",
		ylabel = "dec / degrees"
	)
	
	hist2d!(sample.ra, sample.dec, bins=100, weights=weights, alpha=0.1)
	scatter!(ra0, dec0)
	

	fig
end

# ╔═╡ 619ca573-fd4b-4c65-8996-eab3869f2142
let 
	fig = Figure()
	ax = PolarAxis(fig[1,1])

	ϕ = atan.(sample.eta, sample.xi) .- deg2rad(PA)
	scatter!(ϕ, sample.r_ell, alpha=0.1)

	fig
end

# ╔═╡ 548d7186-59eb-4f08-81b6-a68127f0df6a
log_r_label = "log r / arcmin"

# ╔═╡ 91df1ddf-a197-4d43-b39c-e25409eef082
let 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=log_r_label,
		ylabel = "sigma",
	)
	
	Arya.errscatter!(profile.log_r, profile.log_Sigma, yerr=profile.log_Sigma_err)

	fig
end

# ╔═╡ dd139766-4281-4a27-af72-428dff73e4c4
profile.counts

# ╔═╡ 3589182c-ee3e-4a0d-ae1a-efc9ac98649a
let
	fig = Figure()
	ax = Axis(fig[1, 1],
		xlabel=log_r_label,
		ylabel = "counts / bin",
		#yscale = log10,
		limits = (nothing, (0., nothing)),
	)
	
	errscatter!(profile.log_r, profile.counts, yerr=sqrt.(profile.counts))

	fig
end

# ╔═╡ d6a899a6-15e8-4fea-b1b7-849f65faac55
profile.log_r

# ╔═╡ 60a223f2-ac5c-4a6a-af79-4b314d7d5509
let
	fig, ax, p = errscatter(profile.log_r, profile.M_in, profile.M_in_err)

	ax.xlabel = log_r_label
	ax.ylabel = "mass interior"


	fig
end

# ╔═╡ c4a1621e-1943-49f1-8d2f-27fa335a0a4f
let
	fig = Figure(size=(700, 300))
	ax = Axis(fig[1, 1], limits=((-0.8, 2), (-10, 10)),
		xlabel=log_r_label,
		ylabel=L"\Gamma"
	)

	
	errscatter!(profile.log_r, profile.Gamma, yerr=profile.Gamma_err)
	

	ax_lin = Axis(fig[1, 2],
		xlabel="r / arcmin",
		yticklabelsvisible=false,
	)

	
	errscatter!(Arya.value.(10 .^ profile.log_r), profile.Gamma, yerr=profile.Gamma_err)

	linkyaxes!(ax, ax_lin)

	fig
end

# ╔═╡ 59b2acc5-66b4-48c6-9507-045ea77e6914
let
	fig = Figure()
	ax = Axis(fig[1, 1], limits=((-1.1, 3), (-2.5, 2.5)),
		xlabel=log_r_label,
		ylabel=L"\Gamma_\textrm{max}"
	)

	
	errscatter!(profile.log_r, profile.Gamma_max, yerr=profile.Gamma_max_err)
	

	fig
end

# ╔═╡ 9bfe90a5-5969-4ebb-96af-72360bbced3b
md"""
# misc
"""

# ╔═╡ 8e75bc6b-2cb1-4993-92ad-037094092612
begin 
	pmra0 = lguys.mean(sample.pmra)
	pmdec0 = lguys.mean(sample.pmdec)
end

# ╔═╡ e081f84f-595f-4cbf-8613-dba2b8f69323
# let
# 	fig = Figure()
# 	ax = Axis(fig[1,1],
# 		xlabel="log r / arcmin",
# 		ylabel="mean proper motion radial"
# 	)

# 	vx = sample.pmra .- pmra0
# 	vy = sample.pmdec .- pmdec0

# 	x = sample.xi
# 	y = sample.eta

# 	r = @. sqrt(x^2 + y^2) * 60

# 	pm_rad = (vx .* x .+ vy .* y ) ./ r

# 	bins, counts = lguys.calc_histogram(log10.(r), weights=pm_rad)
# 	_, N = lguys.calc_histogram(log10.(r))

# 	errscatter!(lguys.midpoint(bins), counts, yerr=counts ./ sqrt.(N))
# 	fig
# end

# ╔═╡ Cell order:
# ╠═142a5ace-1432-4093-bee7-4a85c19b0d72
# ╟─852717c0-aabf-4c03-9cf5-a6d91174e0f9
# ╠═4cc4e2be-6bf6-4cbd-a2b1-121354a862bc
# ╠═c4574ed3-f431-4c0e-a721-1c8d88dda10f
# ╠═bb644db4-7fb4-43c8-abf9-7235aa279ad6
# ╠═73f0b3a1-a4b6-422d-9f7e-be816c4a9cfc
# ╠═a2465c61-ce25-42aa-8b5c-57ad7ffe16f6
# ╠═a82d8fa5-32db-42d1-8b0a-d54ae47dc7be
# ╠═72d975fe-9d97-474b-ba3f-f61ba12c7c80
# ╠═cb38c9f9-d6ff-4bcd-a819-9b442776ccfc
# ╠═305f79e0-a9bb-4c15-a9fd-09cb1c25db41
# ╠═77e1461f-26c3-4b3a-8d16-a9dc095eb572
# ╠═3ade8bb5-8698-4d10-862c-e9685e54e570
# ╠═568162b2-ad89-4c72-92e4-ee3cc302bbb7
# ╠═0272086f-ebeb-432b-b44e-7c35a17c1d58
# ╠═86dd90bc-83dd-4b1a-8e98-1bb0333c6610
# ╠═ef19dcd1-fae0-4777-a0d8-d242435f892f
# ╠═8d276372-add5-4388-b713-b22e38d56f37
# ╠═69018984-ef00-44ef-ba6e-7cccf930aef9
# ╠═f476c859-ba4b-4343-8184-f6f41dc092ee
# ╠═b87cf54b-2c4c-46c5-9336-c13e773e29ec
# ╠═0b80353f-c698-4fae-b40f-1796f7c89792
# ╠═e63df7b9-1fc1-4cc0-a91f-c0f1395d7ff4
# ╠═d382b455-73ba-41d5-bbc8-eca53ea2166e
# ╠═619ca573-fd4b-4c65-8996-eab3869f2142
# ╠═91df1ddf-a197-4d43-b39c-e25409eef082
# ╠═548d7186-59eb-4f08-81b6-a68127f0df6a
# ╠═dd139766-4281-4a27-af72-428dff73e4c4
# ╠═3589182c-ee3e-4a0d-ae1a-efc9ac98649a
# ╠═d6a899a6-15e8-4fea-b1b7-849f65faac55
# ╠═60a223f2-ac5c-4a6a-af79-4b314d7d5509
# ╠═c4a1621e-1943-49f1-8d2f-27fa335a0a4f
# ╠═59b2acc5-66b4-48c6-9507-045ea77e6914
# ╠═9bfe90a5-5969-4ebb-96af-72360bbced3b
# ╠═8e75bc6b-2cb1-4993-92ad-037094092612
# ╠═e081f84f-595f-4cbf-8613-dba2b8f69323
