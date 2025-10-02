### A Pluto.jl notebook ###
# v0.20.18

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

# ╔═╡ d8fe9e91-2ed3-49f3-81d2-d33f16d74352
using Distributions

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ b0a36e69-1f67-4746-9541-f7723e92c13d
distance = 76.5

# ╔═╡ 16989e15-7cc6-43bf-86c3-10b858f094d4
function plot_pace_components(;
	R_h_inner = Normal(0.221, 0.017),
	R_h_outer = Normal(0.374, 0.050),
	f_inner = Normal(0.54, 0.10),

)
	fig = Figure()
	ax = Axis(fig[1,1])

	alpha = 0.1
	Nsamples = 100



	x = LinRange(-2, 1, 1000) # kpc
	R = 10 .^ x
	R_am = LilGuys.kpc2arcmin.(R, distance)
	x_am = log10.(R_am)

	R_trans = zeros(Nsamples)
	for i in 1:Nsamples
		R_in = rand(R_h_inner)
		R_out = rand(R_h_outer)
		f = rand(f_inner)
		prof_in = LilGuys.Plummer(r_s=R_in, M=f)
		prof_out = LilGuys.Plummer(r_s=R_out, M=1-f)

		y_in = @. log10(LilGuys.surface_density(prof_in, R))
		y_out = @. log10(LilGuys.surface_density(prof_out, R))

		lines!(x_am, y_in, alpha=alpha, color=COLORS[1])
		lines!(x_am, y_out, alpha=alpha, color=COLORS[2])

		try
			R_guess = R[argmin(abs.(y_in .- y_out))]
			R_trans[i] = LilGuys.find_zero(R->LilGuys.surface_density(prof_in, R)- LilGuys.surface_density(prof_out, R), (R_guess))
			@info R_guess
		catch Exception
			R_trans[i] = NaN
		end
	end


	ax2 = Axis(fig[0, 1])
	x = LilGuys.kpc2arcmin.(R_trans[R_trans .< 10], distance)
	hist!(log10.(x), LinRange(-2, 1, 100))

	linkxaxes!(ax, ax2)

	fig

end

# ╔═╡ 67351af0-b9a6-4829-bb01-958a7c53355a
plot_pace_components()

# ╔═╡ 7e2602e3-aa76-4785-94b9-d403d285ad09
plot_pace_components(R_h_inner=Normal(0.253, 0.018), R_h_outer=Normal(0.512, 0.145), f_inner=Normal(0.76, 0.10))

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═d8fe9e91-2ed3-49f3-81d2-d33f16d74352
# ╠═b0a36e69-1f67-4746-9541-f7723e92c13d
# ╠═16989e15-7cc6-43bf-86c3-10b858f094d4
# ╠═67351af0-b9a6-4829-bb01-958a7c53355a
# ╠═7e2602e3-aa76-4785-94b9-d403d285ad09
