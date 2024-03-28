### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 59d999de-dfe8-11ee-0ba0-d7bde6f4231a
begin 
	import Pkg; Pkg.activate()
	
	using FITSIO
	using Plots; 
	using Arya
	
	using DataFrames 
	using CSV
	using LaTeXStrings
	
	import LilGuys as lguys
end

# ╔═╡ 9eaa93cd-3f49-45e8-b83c-721c10b1671b
begin 
	using Turing
	using Optim
	using StatsPlots
end

# ╔═╡ a0d149b6-8f16-4573-8c93-b670fa865e40
md"""
Given a list of densities with uncertanties (sculptor_densities.csv), which profiles best fit?
"""

# ╔═╡ 058e97e0-2284-4a9b-ab09-4d51fd01c0b8
@model function linear_reg(x, y)
	# priors
	γ ~ Normal(0, 2) # expect ~-1
	α ~ Normal(3, 2)
	σ ~ Exponential(0.4)

	# likelyhood
	for i in eachindex(y)
		y[i] ~ Normal(α + γ * x[i], σ)
	end
end

# ╔═╡ 4103cc1e-c7db-4ab2-8d45-e1922ec018f0
ch_prior = sample(linear_reg([], []), Prior(), 1000)

# ╔═╡ 73216801-4923-48fc-aedc-4d27939b2d95
plot(ch_prior)

# ╔═╡ f01dcf6f-1732-4406-a233-7b7f78d9e507
begin 
	p1 = plot()
	for i in 1:10:length(ch_prior)
		a = ch_prior[:α].data[i]
		b = ch_prior[:γ].data[i]
		y_model = @. (a + rs * b)
		plot!(log10.(rs), y_model, alpha=0.1, legend=false, color="black" )
	end
	
	p1
end

# ╔═╡ 716982e1-7371-4d43-8b8c-c90268179879
begin 
	nf = not.(isnan.(r_mid) .| isnan.(Σs))
	model = linear_reg(r_mid[nf], log10.(Σs[nf]))
end

# ╔═╡ 3bbab706-af87-47ce-a565-7d55fdd8d5ca
map_estimate = optimize(model, MAP())

# ╔═╡ b83b8d4a-e33b-40b7-a5e2-c1413de15436
chain = sample(model, NUTS(), 1000, initial_params=map_estimate.values.array)

# ╔═╡ 366373bb-02f2-406a-9edb-7ff11f587313
plot(chain)

# ╔═╡ Cell order:
# ╟─a0d149b6-8f16-4573-8c93-b670fa865e40
# ╠═59d999de-dfe8-11ee-0ba0-d7bde6f4231a
# ╠═9eaa93cd-3f49-45e8-b83c-721c10b1671b
# ╠═058e97e0-2284-4a9b-ab09-4d51fd01c0b8
# ╠═4103cc1e-c7db-4ab2-8d45-e1922ec018f0
# ╠═73216801-4923-48fc-aedc-4d27939b2d95
# ╠═f01dcf6f-1732-4406-a233-7b7f78d9e507
# ╠═716982e1-7371-4d43-8b8c-c90268179879
# ╠═3bbab706-af87-47ce-a565-7d55fdd8d5ca
# ╠═b83b8d4a-e33b-40b7-a5e2-c1413de15436
# ╠═366373bb-02f2-406a-9edb-7ff11f587313
