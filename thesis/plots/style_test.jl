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

# ╔═╡ e25d463d-3fb1-4eb9-a3aa-95a64179558f
CairoMakie.activate!(pt_per_unit=5, type="svg") # display in higher resolution

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y"
	)

	scatter!(randn(10), rand(10))

	fig
end

# ╔═╡ da03103f-b36c-469a-adbb-b5d2b1b11d59
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y"
	)

	scatter!(randn(10), rand(10), label = "a")
	lines!(randn(10), rand(10), label = "b")
	errorscatter!(randn(10), rand(10), yerror=rand(10)/3, label = "c")

	axislegend()
	fig
end

# ╔═╡ 4e8012d9-020f-4fc9-bfc6-9aba9dc3e1e3
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y"
	)

	p = scatter!(randn(10), rand(10), color=randn(10))

	Colorbar(fig[1,2], p, minorticksize=2)
	fig
end

# ╔═╡ 018c7149-42ff-4777-8302-837490421308
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y",
		xminorticks = SmartMinorTicks(),
	)

	scatter!(randn(10), rand(10))

	fig
end

# ╔═╡ 2d2a0a07-3dc4-45f3-90ea-859d2389a23d
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y",
		xminorticks = SmartMinorTicks(),
		xticks = DefaultLinearTicks,
	)

	scatter!(2rand(10), rand(10))

	fig
end

# ╔═╡ 2c8355ef-7eff-48f2-90ba-dcac16d158f6
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y",
		xminorticks = SmartMinorTicks(),
		xticks = DefaultLinearTicks,
	)

	scatter!(LinRange(1, 9, 10), rand(10))

	fig
end

# ╔═╡ 95655883-722c-4982-a5b5-2e87dc8eccd8
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y",
		xminorticks = SmartMinorTicks(),
		xticks = DefaultLinearTicks,
	)

	scatter!(LinRange(1, 12, 10), rand(10))

	fig
end

# ╔═╡ 8d2ada36-fbbf-4015-81db-67fca93af7f6
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y",
		xminorticks = SmartMinorTicks(),
	)

	scatter!(LinRange(1, 14.5, 10), rand(10))

	fig
end

# ╔═╡ 8ce0cf4f-c5b8-484f-a319-c65c462aa103
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y",
		xscale = log10,
		xticks = Makie.LogTicks(Arya.DefaultLinearTicks),
		xminorticks = SmartMinorTicks(),
	)

	scatter!(10 .^ LinRange(0, 3, 100), rand(100))

	fig
end

# ╔═╡ d589b682-a116-4dc0-9829-aa4168b052be
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y",
		xscale = log10,
		xticks = Makie.LogTicks(Arya.DefaultLinearTicks),
		xminorticks = SmartMinorTicks(),
	)

	scatter!(10 .^ LinRange(0, 1, 100), rand(100))

	fig
end

# ╔═╡ eba32611-91be-4185-b409-2c06944d6821
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y",
		xscale = log10,
	)

	scatter!(10 .^ LinRange(0, 3, 100), rand(100))

	fig
end

# ╔═╡ 30a5c7f2-c24b-4b82-9290-aeb3569d2b63
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = "x",
		ylabel = "y",
		xscale = log10,
		xminorticks = SmartMinorTicks(),
	)

	scatter!(10 .^ LinRange(0, 0.5, 100), rand(100))

	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═e25d463d-3fb1-4eb9-a3aa-95a64179558f
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═da03103f-b36c-469a-adbb-b5d2b1b11d59
# ╠═4e8012d9-020f-4fc9-bfc6-9aba9dc3e1e3
# ╠═018c7149-42ff-4777-8302-837490421308
# ╠═2d2a0a07-3dc4-45f3-90ea-859d2389a23d
# ╠═2c8355ef-7eff-48f2-90ba-dcac16d158f6
# ╠═95655883-722c-4982-a5b5-2e87dc8eccd8
# ╠═8d2ada36-fbbf-4015-81db-67fca93af7f6
# ╠═8ce0cf4f-c5b8-484f-a319-c65c462aa103
# ╠═d589b682-a116-4dc0-9829-aa4168b052be
# ╠═eba32611-91be-4185-b409-2c06944d6821
# ╠═30a5c7f2-c24b-4b82-9290-aeb3569d2b63
