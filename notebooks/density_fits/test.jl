### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ e602e986-275c-11ef-3bea-edf5f1dca0de
begin
	import Pkg; Pkg.activate()

	using Arya
	using GLMakie

	import LilGuys as lguys
end

# ╔═╡ 4d2ecbf7-3e1e-4d96-98b3-eb6eb4a7492e
using FITSIO, DataFrames

# ╔═╡ fa75f03c-a35c-438f-82db-3e30b9570754
include("density_utils.jl")

# ╔═╡ c0ea6c9c-2a74-47ba-a259-ba0e245527f2
include("filter_utils.jl")

# ╔═╡ 996434e3-e5b1-44fc-8f85-ef8b49d5886e
samplename = "sculptor/fiducial_sample.fits" # = "$(name)_sample.fits" 

# ╔═╡ 4153c813-bfa5-4e8b-a07d-cc92c923b8fd
sample = load_fits(samplename)

# ╔═╡ f01e7646-3881-43d4-b24f-610106ff21ea
begin 
	r, filt = calc_radii(sample.ra, sample.dec)
	r = r[filt]
end

# ╔═╡ 60d935fd-f59c-4d97-9894-5bf0990741d6
h = Arya.histogram(log10.(r), 10, weights=ones(length(r)))

# ╔═╡ 276eb1eb-870e-4ac2-848f-f20dfb04a358
calc_Σ_mean(h.bins, h.values)

# ╔═╡ 56b33dee-8c9e-40ce-bbea-f16dc4a710ed
prof = calc_properties(r)

# ╔═╡ bf2c0946-c6ae-4984-9d75-079bd414090e
open("test.toml", "w") do f
	print(f, prof)
end

# ╔═╡ Cell order:
# ╠═e602e986-275c-11ef-3bea-edf5f1dca0de
# ╠═996434e3-e5b1-44fc-8f85-ef8b49d5886e
# ╠═4d2ecbf7-3e1e-4d96-98b3-eb6eb4a7492e
# ╠═fa75f03c-a35c-438f-82db-3e30b9570754
# ╠═c0ea6c9c-2a74-47ba-a259-ba0e245527f2
# ╠═4153c813-bfa5-4e8b-a07d-cc92c923b8fd
# ╠═f01e7646-3881-43d4-b24f-610106ff21ea
# ╠═60d935fd-f59c-4d97-9894-5bf0990741d6
# ╠═276eb1eb-870e-4ac2-848f-f20dfb04a358
# ╠═56b33dee-8c9e-40ce-bbea-f16dc4a710ed
# ╠═bf2c0946-c6ae-4984-9d75-079bd414090e
