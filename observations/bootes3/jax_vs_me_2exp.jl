### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ fc5c6e7c-c4aa-11f0-ac5d-2b084e6256fa
begin
	import Pkg
	Pkg.activate()

	using LilGuys
	using Arya
	using CairoMakie
	using DataFrames, CSV
end

# ╔═╡ e0c0b153-92d9-47a4-b652-a7ecef76dce1
import TOML

# ╔═╡ 060b4e91-614d-4996-95d9-ed4d4f750338
samples = CSV.read("mcmc/samples.mcmc_2exp.csv", DataFrame, ntasks=1, delim=",")

# ╔═╡ 0cea472f-f023-4a31-9493-3064be1d7c8e
obs_props = TOML.parsefile("./observed_properties.toml") |> LilGuys.collapse_errors

# ╔═╡ 875e3d04-d266-433d-9aba-8704b9a20401
α = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ e80cfad4-8893-4428-ac6e-db37cd6a85a1
r_h = @. sqrt(1 + samples.ellipticity) .* samples.R_s * α

# ╔═╡ 85a053a7-9fce-471b-8804-11f655c8ad9d
r_h_outer = @. sqrt(1 + samples.ellipticity_outer) .* samples.R_s_outer * α

# ╔═╡ 1ea2725b-5ce1-4748-9ffe-a734792f4c8d
function plot_meas!(name::String)
	meas = obs_props[name]

	plot_meas!(meas)
end

# ╔═╡ 65f1c84d-d3d8-424c-b830-0b2dd0bc2dcb
function plot_meas!(meas)
	vlines!(middle(meas), color=:black)
	vspan!(middle(meas) - lower_error(meas), middle(meas) + upper_error(meas), color=:black, alpha=0.2)
end

# ╔═╡ 27975a2a-ef57-4b00-8fa3-da25a49e4687
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "r_h / arcmin (major axis)")


	hist!(r_h, bins=100)

	plot_meas!("r_h")
	fig
end

# ╔═╡ bdcf5ec1-6942-4552-9892-2077444bb86a
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "position_angle")


	hist!(samples.position_angle .- 180, bins=100)
	stephist!(samples.position_angle_outer .- 180, bins=100, color=COLORS[2])

	plot_meas!("position_angle")
	fig
end

# ╔═╡ d86c1fed-ac83-4281-bb0e-5898a30da1e3
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "ellipticity")


	hist!(samples.ellipticity, bins=100)
	stephist!(samples.ellipticity_outer, bins=100, color=COLORS[2])

	plot_meas!("ellipticity")
	fig
end

# ╔═╡ be7de1ad-40e1-4ae1-9bbb-b02eb8b9ea12
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log10 r_h outer")


	stephist!(log10.(r_h_outer), bins=100, color=COLORS[2])

	vlines!(log10(1.901 * 60), color=:black)
	vspan!(log10.((1.901 .+ [-1.841, 0.554]) * 60)..., color=:black, alpha=0.1)
	xlims!(0, 4)
	fig
end

# ╔═╡ 84e06b07-ba0b-4764-9538-b1daeef85921
B = Measurement(0.015, 0.025, 0.06)

# ╔═╡ bb887685-0085-4902-b935-9b36d411c133
md"""
the mass of an exponential with central surface density 1 is 2π R_s^2
"""

# ╔═╡ 46734c68-ebbb-424e-81c5-92d19b9a3293
M1_jax = 2π * obs_props["r_h"].middle^2 #/ (1 + obs_props["ellipticity"].middle)

# ╔═╡ 796b612a-b117-498d-9620-6f41c14b2e99
M2_jax = 2π * (1.901 * 60 * α)^2 #/ (1 + obs_props["ellipticity"].middle)

# ╔═╡ 4b7c68f4-c91b-4d6d-8d83-602842558a32
f_outer_jax = 1 / (1 + M1_jax / (B * M2_jax))

# ╔═╡ a52319e9-5268-41e5-a734-c76dd95ccb90
1 / (1 + M1_jax / ((B.middle - B.lower) * M2_jax))

# ╔═╡ 1c466025-ff24-452f-9a1b-e60a8f3c3a38
1 / (1 + M1_jax / ((B.middle + B.upper) * M2_jax))

# ╔═╡ 4ab5b332-2792-40ea-a567-71a573b0d3fb
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "f_outer")


	hist!(samples.f_outer, bins=100)

	vlines!(f_outer_jax.middle, color=:black)
	fig
end

# ╔═╡ 9b533deb-de2e-45d5-a356-e90ec588504f
function f_to_B(f, R_s_1, R_s_2)
	return (1/2π * R_s_1^2 * f) / (1/2π * R_s_2^2 * (1-f))
end

# ╔═╡ 06b8bf03-6bf3-4f6e-bed6-0fa14f4228a5
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "B (inner / outer central surface density)")


	hist!(f_to_B.(samples.f_outer, samples.R_s, samples.R_s_outer), bins=LinRange(0, 0.3, 100))

	plot_meas!(B)
	xlims!(-0.02, 0.3)
	fig
end

# ╔═╡ 55265f89-8e1f-47b5-90cb-1397c87094f6
f_to_B.(samples.f_outer, samples.R_s, samples.R_s_outer)

# ╔═╡ eea4b2ff-9a27-4025-8199-649e4bc45229
f_to_B(f_outer_jax, obs_props["r_h"].middle, 1.901 * 60 * α)

# ╔═╡ 04e179dc-0bc6-4c09-8053-31187a4713d0
B_me = @. f_to_B(samples.f_outer, samples.R_s_outer, samples.R_s)

# ╔═╡ bdd0fd2c-4ae6-4c48-9cde-f4f44b35f5af
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "r_h / arcmin (major axis)")


	hist!(r_h)
	fig
end

# ╔═╡ f751f4c8-3211-4fe8-a614-b6c9daa41d60
LilGuys.quantile(r_h, [0.5, 0.16, 0.84])

# ╔═╡ a32d53c7-6785-4da2-b4f0-dc31e374542d
let
	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "log r_h / arcmin (major axis)", ylabel = "log r_h_outer / log r_h")
	
	scatter!(log10.(r_h), log10.(r_h_outer) .- log10.(r_h), markersize=1, alpha=0.01, color=:black)
	
	fig
end

# ╔═╡ 6f23bfe1-902d-44a7-af3b-62c7a8ccbc7a
hist(log10.(r_h_outer))

# ╔═╡ Cell order:
# ╠═fc5c6e7c-c4aa-11f0-ac5d-2b084e6256fa
# ╠═e0c0b153-92d9-47a4-b652-a7ecef76dce1
# ╠═060b4e91-614d-4996-95d9-ed4d4f750338
# ╠═0cea472f-f023-4a31-9493-3064be1d7c8e
# ╠═875e3d04-d266-433d-9aba-8704b9a20401
# ╠═e80cfad4-8893-4428-ac6e-db37cd6a85a1
# ╠═85a053a7-9fce-471b-8804-11f655c8ad9d
# ╠═1ea2725b-5ce1-4748-9ffe-a734792f4c8d
# ╠═65f1c84d-d3d8-424c-b830-0b2dd0bc2dcb
# ╠═27975a2a-ef57-4b00-8fa3-da25a49e4687
# ╠═bdcf5ec1-6942-4552-9892-2077444bb86a
# ╠═d86c1fed-ac83-4281-bb0e-5898a30da1e3
# ╠═be7de1ad-40e1-4ae1-9bbb-b02eb8b9ea12
# ╠═84e06b07-ba0b-4764-9538-b1daeef85921
# ╠═4b7c68f4-c91b-4d6d-8d83-602842558a32
# ╠═a52319e9-5268-41e5-a734-c76dd95ccb90
# ╠═1c466025-ff24-452f-9a1b-e60a8f3c3a38
# ╠═bb887685-0085-4902-b935-9b36d411c133
# ╠═46734c68-ebbb-424e-81c5-92d19b9a3293
# ╠═796b612a-b117-498d-9620-6f41c14b2e99
# ╠═4ab5b332-2792-40ea-a567-71a573b0d3fb
# ╠═06b8bf03-6bf3-4f6e-bed6-0fa14f4228a5
# ╠═55265f89-8e1f-47b5-90cb-1397c87094f6
# ╠═9b533deb-de2e-45d5-a356-e90ec588504f
# ╠═eea4b2ff-9a27-4025-8199-649e4bc45229
# ╠═04e179dc-0bc6-4c09-8053-31187a4713d0
# ╠═bdd0fd2c-4ae6-4c48-9cde-f4f44b35f5af
# ╠═f751f4c8-3211-4fe8-a614-b6c9daa41d60
# ╠═a32d53c7-6785-4da2-b4f0-dc31e374542d
# ╠═6f23bfe1-902d-44a7-af3b-62c7a8ccbc7a
