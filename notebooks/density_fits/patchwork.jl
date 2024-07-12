### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 2a0d82b2-3fe5-11ef-224e-e352922a2e5a
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames, CSV
	
	using CairoMakie

	
	import LilGuys as lguys
	using Arya
	
end

# ╔═╡ ae6cdb42-4193-436e-ab6a-008f9c6b7b32
using KernelDensity

# ╔═╡ 218c549f-a8a9-4795-ac07-6e7aa3a5e758
filenames = ["centre", "leading", "trailing", "background"]

# ╔═╡ eb97c46b-5b9d-4ae6-82c1-2f7c409b69ab
import TOML

# ╔═╡ 13834fc5-04c5-44ed-ad34-a87408bfbe9b
scl_props = TOML.parsefile("/astro/dboyea/dwarfs/sculptor_obs_properties.toml")

# ╔═╡ 1aec8806-29bd-43e3-9b39-5719d9012c7c
function load_sample(filename)
	return lguys.load_fits("sculptor/$(filename)_sample.fits")
end

# ╔═╡ 72bd86f1-41ce-4727-979f-9262658f1873
load_sample("simple")

# ╔═╡ 5ec3651f-653d-4a43-bd09-da523c49f4cb
samples = load_sample.(filenames)

# ╔═╡ 67ed76ef-0da3-4c70-bf46-593b819f981e
ra_all = vcat([sample.ra for sample in samples]...)

# ╔═╡ 01f7fec7-ef67-45d7-a651-525f83eb0579
dec_all = vcat([sample.dec for sample in samples]...)

# ╔═╡ a8627570-50c3-48bd-a343-4ad6218ed064
allstars = unique(
	vcat(samples...), [:source_id], keep=:first)

# ╔═╡ b5593f3c-bc6d-4e8a-acc4-22eef851859f
unique(allstars.source_id)

# ╔═╡ 4683ce6c-4fce-4e05-afde-b50f2cbb3c9c
xi, eta = lguys.to_tangent(allstars.ra, allstars.dec, scl_props["ra"], scl_props["dec"])

# ╔═╡ f328d46f-8ede-4c91-97f5-724dac51064d
bandwidth = (0.05, 0.05)

# ╔═╡ 7e9fd578-2789-428f-acb9-aa700b8569b4
kd = kde([xi eta], bandwidth=bandwidth)

# ╔═╡ 146894f9-0da8-4be5-90be-677566abdfcf
let
	fig = Figure()
	ax = Axis(fig[1, 1])


	scatter!(ra_all, dec_all, alpha=0.05)

	fig
end

# ╔═╡ d7edb06b-e54a-46cc-9b74-0914931ac507
let
	fig = Figure()
	ax = Axis(fig[1, 1], aspect=DataAspect(),
		xlabel = L"$\xi$ / degrees",
		ylabel = L"$\eta$ / degrees",
		xreversed = true,
		#limits=(-r_ell_max, r_ell_max, -r_ell_max, r_ell_max) ./ 60
	)
	scale = 0.008


	h = heatmap!(kd.x, kd.y, asinh.(kd.density ./ scale))

	#lines!(orbit.xi, orbit.eta)
	Colorbar(fig[1, 2], h, label="asinh density / $scale")
	fig
end

# ╔═╡ Cell order:
# ╠═2a0d82b2-3fe5-11ef-224e-e352922a2e5a
# ╠═218c549f-a8a9-4795-ac07-6e7aa3a5e758
# ╠═eb97c46b-5b9d-4ae6-82c1-2f7c409b69ab
# ╠═13834fc5-04c5-44ed-ad34-a87408bfbe9b
# ╠═1aec8806-29bd-43e3-9b39-5719d9012c7c
# ╠═72bd86f1-41ce-4727-979f-9262658f1873
# ╠═5ec3651f-653d-4a43-bd09-da523c49f4cb
# ╠═67ed76ef-0da3-4c70-bf46-593b819f981e
# ╠═01f7fec7-ef67-45d7-a651-525f83eb0579
# ╠═a8627570-50c3-48bd-a343-4ad6218ed064
# ╠═ae6cdb42-4193-436e-ab6a-008f9c6b7b32
# ╠═b5593f3c-bc6d-4e8a-acc4-22eef851859f
# ╠═4683ce6c-4fce-4e05-afde-b50f2cbb3c9c
# ╠═7e9fd578-2789-428f-acb9-aa700b8569b4
# ╠═f328d46f-8ede-4c91-97f5-724dac51064d
# ╠═146894f9-0da8-4be5-90be-677566abdfcf
# ╠═d7edb06b-e54a-46cc-9b74-0914931ac507
