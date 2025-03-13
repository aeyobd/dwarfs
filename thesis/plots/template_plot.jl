### A Pluto.jl notebook ###
# v0.20.4

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

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "test" let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y"
	)

	scatter!(randn(10), rand(10))

	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
