### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./style.jl")

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "test" let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y"
	)

	for i in 1:3
		scatter!(randn(10), rand(10))

	end
	
	fig
end

# ╔═╡ 82f72204-4340-41a6-a7c4-4fc8bc86d87a
@savefig "test" let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "∫x²",
		ylabel = "y"
	)

	for i in 1:3
		lines!(0.1:0.1:10, 100rand() .+ rand(100))

	end
	
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═82f72204-4340-41a6-a7c4-4fc8bc86d87a
