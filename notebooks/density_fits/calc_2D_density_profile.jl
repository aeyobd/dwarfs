### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 142a5ace-1432-4093-bee7-4a85c19b0d72
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames, CSV
	
	using CairoMakie
	using Measurements
	#using KernelDensity
	
	import SciPy
	using QuadGK
	
	import LinearAlgebra: diag
	
	import LilGuys as lguys
	using Arya
	
	using JSON
end

# ╔═╡ 852717c0-aabf-4c03-9cf5-a6d91174e0f9
md"""
given a sample of points, can we centre and calculate the 2D density profile

### Inputs
- Number of bins
- Ellipticity
- Centring method
"""

# ╔═╡ 73f0b3a1-a4b6-422d-9f7e-be816c4a9cfc
begin 
	samplename = "/cosma/home/durham/dc-boye1/sculptor/orbits/orbit1/exp2d_stars_today.fits" # = "$(name)_sample.fits" 
	samplename = "sculptor/fiducial_sample.fits" # = "$(name)_sample.fits" 
end

# ╔═╡ a2465c61-ce25-42aa-8b5c-57ad7ffe16f6
outname = splitext(samplename)[1]  * "_profile.json"

# ╔═╡ a82d8fa5-32db-42d1-8b0a-d54ae47dc7be
begin 
	ecc = 0.37
	PA = 94
	centre_method="mean"
	mass_column = nothing
end

# ╔═╡ ac49a334-230a-42e4-b655-a472cd505688
function add_xi_eta!(stars, ra0, dec0)
	xi, eta = lguys.to_tangent(stars.ra, stars.dec, ra0, dec0)
	
	stars[:, "xi"] = xi
	stars[:, "eta"] = eta
	stars
end

# ╔═╡ 4b132b0d-bea3-4fe1-9a91-fdebcf2145b7
function filter_edges!(sample)
	r_max = sqrt(maximum(sample.xi .^ 2 + sample.eta .^ 2))
	r_cut = r_max * sqrt(1 - ecc)

	filt = sample.r_ell .< r_cut

	sample = sample[filt, :]

	Ncut = sum(map(!, filt))
	println("cut $Ncut stars on edge ")
end

# ╔═╡ 514027d9-6e70-41f7-b134-526c142182c6
md"""
# Centre finding
"""

# ╔═╡ 5dbfcbc7-b21b-45c0-b8d4-52f215a90d54
function mean_centre(stars)
	return (lguys.mean(stars.ra), lguys.mean(stars.dec))
end

# ╔═╡ fa6dcca1-6a54-4841-9e44-30718b7f67f9
import StatsBase: weights

# ╔═╡ ccf2ed52-439c-48fd-9ade-61cb3d6199d7
function weighted_centre(stars, w)
	w = weights(w)
	return (lguys.mean(stars.ra, w), lguys.mean(stars.dec, w))
end

# ╔═╡ f6f2d0f1-d57a-42db-99b9-02c5315ed58b
function calc_centre(sample)
	if centre_method == "mean"
		ra0, dec0 = mean_centre(sample)
	elseif center_method == "weighted"
		ra0, dec0 = weighted_centre(sample, sample.probability)
	end

	return ra0, dec0
end

# ╔═╡ 5b8e641c-29d9-42a3-8149-9e2e5b3d46cb
md"""
# Calculating Radii
"""

# ╔═╡ 9b067a18-d2eb-4112-91cf-5814fc02b388
function add_r_ell!(stars, ecc, PA)
	b = sqrt(1 - ecc)
	a = 1/b
	
	println("a, b / r_h = $a, $b")
	r_ell = lguys.calc_r_ell(stars.xi, stars.eta, a, b, PA-90)
	stars[:, "r_ell"] = r_ell;
end

# ╔═╡ 72d975fe-9d97-474b-ba3f-f61ba12c7c80
begin 
	f = FITS(samplename, "r")
	sample = DataFrame(f[2])
	close(f)

	ra0, dec0 = calc_centre(sample)
	add_xi_eta!(sample, ra0, dec0)
	add_r_ell!(sample, ecc, PA)

	if mass_column === nothing
		sample[!, "mass"] = ones(size(sample, 1))
	else
		sample[!, "mass"] = sample[mass_column]
	end

	filter_edges!(sample)
end

# ╔═╡ b063d5f4-ad59-4bf3-bd88-a93d8b27fdf4
let
	p = scatter(sample.ra, sample.dec, alpha=0.1)
	scatter!(ra0, dec0)
	p
end

# ╔═╡ d382b455-73ba-41d5-bbc8-eca53ea2166e
let
	fig = Figure()
	ax = Axis(fig[1, 1],
	limits = (ra0 .+ (-1, 1), dec0 .+ (-1, 1)))
	
	scatter!(sample.ra, sample.dec, alpha=0.1)
	scatter!(ra0, dec0)
	

	fig
end

# ╔═╡ 619ca573-fd4b-4c65-8996-eab3869f2142
let 
	fig = Figure()
	ax = PolarAxis(fig[1,1])


	ϕ = atan.(sample.eta, sample.xi) .- deg2rad(PA)
	scatter!(ϕ, sample.r_ell)

	fig
end

# ╔═╡ f3408055-ebec-41f3-8da1-e83f010da5bf
rs = sample.r_ell * 60

# ╔═╡ 76db218a-c75d-4189-9e3a-ee1bf97bab4b
md"""
# Properties
"""

# ╔═╡ a34bc111-e1ec-4f14-a759-27702d33abf2
function norm_hist(xs, bw; weights=nothing)
	bins = collect(minimum(xs):bw:(maximum(xs) + bw))
	x_h, y_h = lguys.calc_histogram(xs, bins, weights=weights)
	if weights === nothing
    	y_e = sqrt.(y_h)
		counts = y_h

	else
		_, counts = lguys.calc_histogram(xs, bins)
		y_e = y_h ./ sqrt.(counts)
	end

	area = sum(diff(bins) .* y_h)
    
	
	return x_h, y_h ./ area, y_e ./ area, counts
end

# ╔═╡ bd78829c-17be-11ef-0fac-53f9dbeb9a51
function running_hist(xs, bw, normalize=false)
	N = length(xs)
	hist = zeros(N)

	x_sort = sort(xs)
	for i in 1:N
		x0 = x_sort[i]
		count = sum(x0 + bw .> x_sort .> x0 - bw)
		hist[i] = count
	end

	if normalize
		dx = lguys.gradient(x_sort)
		area = sum(hist .* dx)
		hist ./= area
	end
	return x_sort, hist
end

# ╔═╡ eefa30dd-c368-4d4d-b043-711f231ebdb9
function calc_Σ(log_r, hist)
	r = 10 .^ log_r
	Σ = hist ./ (2π * log(10) * r .^ 2) # is equivalent because of derivative of log r
	return Σ
end

# ╔═╡ 1c30e04b-9110-4266-94a4-0c8852118cd8
function calc_Σ_mean(log_r, hist)
	r = 10 .^ lguys.midpoint(log_r)
	counts = cumsum(hist .* diff(log_r))
	Areas = @. π * r^2
	σ = counts ./ Areas
	return σ
end

# ╔═╡ 3b97ac3c-b4ae-4b4c-878b-f0aff7c8f0e9
function calc_Γ(log_rs, Σs, step=1)
	dx = lguys.gradient(log_rs)
	dρ = lguys.gradient(log10.(Σs))

	return dρ ./ dx #lguys.gradient(log10.(ρs), log_rs)
end

# ╔═╡ bdd6d321-1026-4d0f-9e45-fbc2c5a2a61c
function calc_properties(rs; weights=nothing, bw=0.1)
    log_r_bin, mass, δ_counts, counts = norm_hist(log10.(rs), weights=weights, 0.1)
    mass = mass .± δ_counts
    log_r = lguys.midpoint(log_r_bin)
    δ_log_r = diff(log_r_bin) ./ 2
    log_r = log_r .± δ_log_r
    r_bin = 10 .^ log_r_bin

    r = 10 .^ log_r

    As = π * diff((10 .^ log_r_bin) .^ 2)
    Σ = mass ./ As 
    # Σ_e = @. y_e / ys * Σ
    #Σ_m = calc_Σ_mean(xs, ys)

    M_in = cumsum(mass)
    A_in = @. π * (r_bin[2:end])^2
    Σ_m = M_in ./ A_in

	println(Σ)
    Γ = calc_Γ(log_r, Σ)
    Γ_max = @. 2*(1 - Σ / Σ_m)

	log_Σ = log10.(Σ)


    
    return Dict{String, Any}(
		"log_r"=>value.(log_r), 
		"log_r_bins"=>log_r_bin, 
		"log_r_units"=>"log arcmin",
		"counts"=>value.(counts),
		"mass"=>value.(mass),
		"mass_err"=>err.(mass),
		"surface_dens"=>value.(Σ), 
		"surface_dens_err"=>err.(Σ), 
		"log_Sigma"=>value.(log_Σ),
		"log_Sigma_err"=>err.(log_Σ),
		"surface_dens_mean"=>value.(Σ_m), 
		"Gamma"=>value.(Γ), 
		"Gamma_err"=>err.(Γ), 
		"Gamma_max"=>value.(Γ_max), 
		"Gamma_max_err"=>err.(Γ_max), 
		"M_in"=>value.(M_in),
		#"N"=>[length(rs)]
	)
end

# ╔═╡ e283c660-056c-4572-845c-7a90fdbfc79b
obs = calc_properties(rs, weights=sample.mass, bw=1)

# ╔═╡ 779b74ce-dbf9-4bea-b28c-ead21d75070c
let 
	fig = Figure()
	ax = Axis(fig[1,1])
	
	stephist!(log10.(rs), weights=sample.mass)

	x, y = lguys.calc_histogram(log10.(rs), weights=sample.mass)
	scatter!(lguys.midpoint(x), y)
	fig
end

# ╔═╡ e62f3840-11db-4123-ab14-6cebe02891af
obs["counts"]

# ╔═╡ 548d7186-59eb-4f08-81b6-a68127f0df6a
log_r_label = "log r / arcmin"

# ╔═╡ 91df1ddf-a197-4d43-b39c-e25409eef082
let 
	fig = Figure()
	ax = Axis(fig[1,1], 
		xlabel=log_r_label,
		ylabel = "sigma",
	)
	
	Arya.errscatter!(obs["log_r"], obs["log_Sigma"], yerr=obs["log_Sigma_err"])

	fig
end

# ╔═╡ 48e385fd-f532-49ab-929c-b55911095464
obs["counts"]

# ╔═╡ 3589182c-ee3e-4a0d-ae1a-efc9ac98649a
let
	fig, ax, p = errscatter(obs["log_r"], obs["counts"],)

	ax.xlabel = log_r_label
	ax.ylabel = "count / bin"

	ax.yscale=Makie.pseudolog10

	fig
end

# ╔═╡ 60a223f2-ac5c-4a6a-af79-4b314d7d5509
let
	fig, ax, p = errscatter(obs["log_r"], obs["mass"], obs["mass_err"])

	ax.xlabel = log_r_label
	ax.ylabel = "mass / bin"

	ax.yscale=Makie.pseudolog10

	fig
end

# ╔═╡ 31973f72-4ba5-4041-b501-d3a0c4fd12d0
profile = obs

# ╔═╡ c4a1621e-1943-49f1-8d2f-27fa335a0a4f
let
	fig = Figure(size=(700, 300))
	ax = Axis(fig[1, 1], limits=((-0.8, 3), nothing),
		xlabel=log_r_label,
		ylabel=L"\Gamma"
	)

	
	errscatter!(value.(profile["log_r"]), profile["Gamma"], yerr=profile["Gamma_err"])
	

	ax_lin = Axis(fig[1, 2],
		xlabel="r / arcmin",
		yticklabelsvisible=false,
	)

	
	errscatter!(value.(10 .^ profile["log_r"]), profile["Gamma"], yerr=profile["Gamma_err"])

	linkyaxes!(ax, ax_lin)

	fig
end

# ╔═╡ 59b2acc5-66b4-48c6-9507-045ea77e6914
let
	fig = Figure()
	ax = Axis(fig[1, 1], limits=((-0.8, 3), nothing),
		xlabel=log_r_label,
		ylabel=L"\Gamma_\textrm{max}"
	)

	
	errscatter!(value.(profile["log_r"]), profile["Gamma_max"], yerr=profile["Gamma_max_err"])
	

	fig
end

# ╔═╡ b87cf54b-2c4c-46c5-9336-c13e773e29ec
begin 

	open(outname, "w") do f
		JSON.print(f, obs)
	end

	println("wrote data to ", outname)
end

# ╔═╡ ab72c141-57ad-409c-914d-24fb0cd959a8
obs

# ╔═╡ 0b80353f-c698-4fae-b40f-1796f7c89792
sum(obs["counts"]) == length(rs)

# ╔═╡ 5c287839-7404-4386-bbbe-7da219e941c3
obs

# ╔═╡ Cell order:
# ╠═852717c0-aabf-4c03-9cf5-a6d91174e0f9
# ╠═142a5ace-1432-4093-bee7-4a85c19b0d72
# ╠═73f0b3a1-a4b6-422d-9f7e-be816c4a9cfc
# ╠═a2465c61-ce25-42aa-8b5c-57ad7ffe16f6
# ╠═a82d8fa5-32db-42d1-8b0a-d54ae47dc7be
# ╠═72d975fe-9d97-474b-ba3f-f61ba12c7c80
# ╠═f6f2d0f1-d57a-42db-99b9-02c5315ed58b
# ╠═ac49a334-230a-42e4-b655-a472cd505688
# ╠═4b132b0d-bea3-4fe1-9a91-fdebcf2145b7
# ╟─514027d9-6e70-41f7-b134-526c142182c6
# ╠═5dbfcbc7-b21b-45c0-b8d4-52f215a90d54
# ╠═ccf2ed52-439c-48fd-9ade-61cb3d6199d7
# ╠═fa6dcca1-6a54-4841-9e44-30718b7f67f9
# ╠═b063d5f4-ad59-4bf3-bd88-a93d8b27fdf4
# ╠═d382b455-73ba-41d5-bbc8-eca53ea2166e
# ╠═619ca573-fd4b-4c65-8996-eab3869f2142
# ╠═5b8e641c-29d9-42a3-8149-9e2e5b3d46cb
# ╠═9b067a18-d2eb-4112-91cf-5814fc02b388
# ╠═f3408055-ebec-41f3-8da1-e83f010da5bf
# ╠═76db218a-c75d-4189-9e3a-ee1bf97bab4b
# ╠═a34bc111-e1ec-4f14-a759-27702d33abf2
# ╠═bd78829c-17be-11ef-0fac-53f9dbeb9a51
# ╠═eefa30dd-c368-4d4d-b043-711f231ebdb9
# ╠═1c30e04b-9110-4266-94a4-0c8852118cd8
# ╠═3b97ac3c-b4ae-4b4c-878b-f0aff7c8f0e9
# ╠═bdd6d321-1026-4d0f-9e45-fbc2c5a2a61c
# ╠═e283c660-056c-4572-845c-7a90fdbfc79b
# ╠═779b74ce-dbf9-4bea-b28c-ead21d75070c
# ╠═e62f3840-11db-4123-ab14-6cebe02891af
# ╠═91df1ddf-a197-4d43-b39c-e25409eef082
# ╠═548d7186-59eb-4f08-81b6-a68127f0df6a
# ╠═48e385fd-f532-49ab-929c-b55911095464
# ╠═3589182c-ee3e-4a0d-ae1a-efc9ac98649a
# ╠═60a223f2-ac5c-4a6a-af79-4b314d7d5509
# ╠═31973f72-4ba5-4041-b501-d3a0c4fd12d0
# ╠═c4a1621e-1943-49f1-8d2f-27fa335a0a4f
# ╠═59b2acc5-66b4-48c6-9507-045ea77e6914
# ╠═b87cf54b-2c4c-46c5-9336-c13e773e29ec
# ╠═ab72c141-57ad-409c-914d-24fb0cd959a8
# ╠═0b80353f-c698-4fae-b40f-1796f7c89792
# ╠═5c287839-7404-4386-bbbe-7da219e941c3
