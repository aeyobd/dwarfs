### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 142a5ace-1432-4093-bee7-4a85c19b0d72
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames, CSV
	
	using GLMakie
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

# ╔═╡ 13fdf275-e2df-4385-ac9c-4bdd79e44b3f
name = "sculptor/fiducial"

# ╔═╡ 73f0b3a1-a4b6-422d-9f7e-be816c4a9cfc
samplename = "$(name)_sample.fits"

# ╔═╡ a82d8fa5-32db-42d1-8b0a-d54ae47dc7be
begin 
	ecc = 0.3
	PA = 94
	centre_method="mean"
end

# ╔═╡ 72d975fe-9d97-474b-ba3f-f61ba12c7c80
begin 
	f = FITS(samplename, "r")
	sample = DataFrame(f[2])
	close(f)
end

# ╔═╡ 514027d9-6e70-41f7-b134-526c142182c6
md"""
# Centre finding
"""

# ╔═╡ 5dbfcbc7-b21b-45c0-b8d4-52f215a90d54
function mean_centre(stars)
	return (lguys.mean(stars.ra), lguys.mean(stars.dec))
end

# ╔═╡ 7a1d304e-5be0-4b93-bded-b1158085c962
mean_centre(sample)

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

# ╔═╡ f3408055-ebec-41f3-8da1-e83f010da5bf
rs = sample.r_ell * 60

# ╔═╡ 76db218a-c75d-4189-9e3a-ee1bf97bab4b
md"""
# Properties
"""

# ╔═╡ a34bc111-e1ec-4f14-a759-27702d33abf2
function norm_hist(xs, bw)
	bins = collect(minimum(xs):bw:(maximum(xs) + bw))
	x_h, y_h = lguys.calc_histogram(xs, bins)
    y_e = sqrt.(y_h)
    
	area = length(xs)

	return x_h, y_h ./ area, y_e ./ area
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
function calc_properties(rs, bw=0.1)
    log_r_bin, counts, δ_counts = norm_hist(log10.(rs), 0.1)
    counts = counts .± δ_counts
    log_r = lguys.midpoint(log_r_bin)
    δ_log_r = diff(log_r_bin) ./ 2
    log_r = log_r .± δ_log_r
    r_bin = 10 .^ log_r_bin

    r = 10 .^ log_r

    As = π * diff((10 .^ log_r_bin) .^ 2)
    Σ = counts ./ As 
    # Σ_e = @. y_e / ys * Σ
    #Σ_m = calc_Σ_mean(xs, ys)

    M_in = cumsum(counts)
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
		"counts"=>value.(counts) * length(rs),
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
obs = calc_properties(rs)

# ╔═╡ 91df1ddf-a197-4d43-b39c-e25409eef082
Arya.errscatter(obs["log_r"], obs["log_Sigma"], yerr=obs["log_Sigma_err"])

# ╔═╡ ddded03f-e2f3-47de-aca2-398a6209c282
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel="radius",
		ylabel="gamma",
		limits=(nothing, (-4, 1))
	)
	Arya.errscatter!(obs["log_r"], obs["Gamma"], yerr=obs["Gamma_err"])

	fig
end

# ╔═╡ b87cf54b-2c4c-46c5-9336-c13e773e29ec
begin 
	out_name = name * "_profile.json"


	open(out_name, "w") do f
		JSON.print(f, obs)
	end

	println("wrote data to ", out_name)
end

# ╔═╡ ab72c141-57ad-409c-914d-24fb0cd959a8
obs

# ╔═╡ bb547672-fc6b-49cb-8053-326dc7d8d7cb
JSON.json

# ╔═╡ 0b80353f-c698-4fae-b40f-1796f7c89792
sum(obs["counts"]) == length(rs)

# ╔═╡ 5c287839-7404-4386-bbbe-7da219e941c3
obs

# ╔═╡ Cell order:
# ╠═852717c0-aabf-4c03-9cf5-a6d91174e0f9
# ╠═142a5ace-1432-4093-bee7-4a85c19b0d72
# ╠═13fdf275-e2df-4385-ac9c-4bdd79e44b3f
# ╠═73f0b3a1-a4b6-422d-9f7e-be816c4a9cfc
# ╠═a82d8fa5-32db-42d1-8b0a-d54ae47dc7be
# ╠═72d975fe-9d97-474b-ba3f-f61ba12c7c80
# ╟─514027d9-6e70-41f7-b134-526c142182c6
# ╠═5dbfcbc7-b21b-45c0-b8d4-52f215a90d54
# ╠═7a1d304e-5be0-4b93-bded-b1158085c962
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
# ╠═91df1ddf-a197-4d43-b39c-e25409eef082
# ╠═ddded03f-e2f3-47de-aca2-398a6209c282
# ╠═b87cf54b-2c4c-46c5-9336-c13e773e29ec
# ╠═ab72c141-57ad-409c-914d-24fb0cd959a8
# ╠═bb547672-fc6b-49cb-8053-326dc7d8d7cb
# ╠═0b80353f-c698-4fae-b40f-1796f7c89792
# ╠═5c287839-7404-4386-bbbe-7da219e941c3
