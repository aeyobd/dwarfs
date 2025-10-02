### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ b088fc88-9e6a-11f0-3095-3d9c1fe82a36
begin
	import Pkg; Pkg.activate()

	import Agama
	using LilGuys
	
	using CairoMakie
	using Arya
end

# ╔═╡ 91ab3185-a7a3-453d-ac72-9fd026d30e77
import TOML

# ╔═╡ ac195c2d-e80d-49d8-9524-9d440844f59d
module OrbitUtils
	include("../orbit_utils.jl")
end

# ╔═╡ 76bc084b-6714-47f8-9032-90baa407267f
galaxyname = "bootes3"

# ╔═╡ 309abb2f-9c6a-4562-a7da-caebbed9bc03
modelname = "vasiliev+21"

# ╔═╡ 203f2d03-15f5-4553-9eeb-2699635c212c
pot = Agama.Potential(file=joinpath(modelname, "agama_potential.ini"))

# ╔═╡ cea8a018-e315-4968-93b6-f2a59a1bcc44
pot_noreflex = Agama.Potential(file=joinpath("vasiliev+21_noreflex", "agama_potential.ini"))

# ╔═╡ 5678931c-6d87-42d8-8332-bf282c46dce8
pot_mw =  Agama.Potential(file=joinpath("vasiliev+21", "potential_mw.ini"))

# ╔═╡ ad14e47b-c5d4-4aea-aff4-5a466ece1157
pot_mw_static =  Agama.Potential(file=joinpath("vasiliev+21", "potential_mw_init.ini"))

# ╔═╡ 59946955-acdc-41d8-ad67-d7843087f309
pot_lmc =  Agama.Potential(file=joinpath("vasiliev+21", "potential_lmc.ini"))

# ╔═╡ dad0a7e3-db11-4bb3-b30f-3b68b8c1ce29
pot_reflex =  Agama.Potential(file=joinpath("vasiliev+21", "potential_reflex.ini"))

# ╔═╡ 6195c401-2dc7-425b-a082-8bbda16d0ba3
obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))

# ╔═╡ 2ccdac03-56fd-4f30-a237-b473435fe321
icrs = ICRS(obs_props)

# ╔═╡ d0cf64df-d46b-4776-9b46-d1c69a02cc61
units = Agama.VASILIEV_UNITS

# ╔═╡ 77503fde-9be2-4a7b-bed0-847297be2c76
orbit = LilGuys.agama_orbit(pot, icrs, timerange=(0, -3/T2GYR), agama_units=units)

# ╔═╡ f7c3629d-4339-4e8f-b77b-f9c48ea0bbe4
orbit_no = LilGuys.agama_orbit(pot_mw_static, icrs, timerange=(0, -3/T2GYR), agama_units=units)

# ╔═╡ 6aea4a07-3405-4e04-a4d4-b31598455043
orbit_noreflex = LilGuys.agama_orbit(pot_noreflex, icrs, timerange=(0, -3/T2GYR), agama_units=units)

# ╔═╡ 299038ee-d228-4eb9-9a97-708a95668fd2
acc_tot = Agama.acceleration(pot, orbit.positions, units, t=orbit.times)

# ╔═╡ ad95c613-aa71-4b1c-b53b-ad6a9648af32
acc_mw_total = Agama.acceleration(pot_mw, orbit.positions, units, t=orbit.times)

# ╔═╡ a8622194-93a1-4db9-8ae4-b7dde2fd1cc7
acc_mw = Agama.acceleration(pot_mw_static, orbit.positions, units, t=orbit.times)

# ╔═╡ a1b205f3-3399-4f0a-96ce-c18f724ea2bd
acc_mw_response = acc_mw_total .- acc_mw

# ╔═╡ c2dd6977-228b-49c3-bc23-922a57f22c39
acc_lmc = Agama.acceleration(pot_lmc, orbit.positions, units, t=orbit.times)

# ╔═╡ 1d61af93-6638-4c17-9f3b-1bff3d9be69b
acc_reflex = Agama.acceleration(pot_reflex, orbit.positions, units, t=orbit.times)

# ╔═╡ 7ccad48a-7895-46dd-97ff-4e2758146cb3
traj_lmc = OrbitUtils.get_lmc_orbit(modelname)

# ╔═╡ f8a15581-743f-406a-9a62-ebe14bf69ee5
LilGuys.plot_xyz(orbit.positions)

# ╔═╡ c7d41f48-4874-4524-b6c3-fd2c9e224b1b
@savefig "boo3_orbits_w_wo_reflex" LilGuys.plot_xyz(LilGuys.positions.([orbit, orbit_no, orbit_noreflex])..., labels=["+reflex", "MW only", "MW+LMC"])

# ╔═╡ 54f3841c-b5f9-4159-b11b-9695e2725650
LilGuys.plot_xyz(acc_tot, acc_mw .+ acc_lmc .+ acc_reflex)

# ╔═╡ 269cb10d-723a-4e57-ae95-0502c40e5040
times = orbit.times * T2GYR

# ╔═╡ f7e1ed18-25fe-4664-8640-16ed7304200a
let
	fig = Figure(size=(5*72, 3.5*72))
	ax = Axis(fig[1,1],
			 xlabel = "time", 
			 ylabel = "acceleration",
			  yscale = log10,
			  yticks=Makie.automatic,
			 )


	lines!(times, radii(acc_tot), label="total")
	lines!(times, radii(acc_mw), label="mw")
	lines!(times, radii(acc_lmc), label="lmc")
	lines!(times, radii(acc_reflex), label="reflex")
	lines!(times, radii(acc_mw_response), label="mw response")


	Legend(fig[1,2], ax)
	@savefig "boo3_acceleration_breakdown" 
	fig
end

# ╔═╡ afbf13a0-03ac-45e4-b529-3f30202b9107
function xyz_arrows!(ax, x0_vec, x_vec; kwargs...)
	x = @lift ($x_vec)[1:1]
	y = @lift ($x_vec)[2:2]
	z = @lift ($x_vec)[3:3]


	x0 = @lift ($x0_vec)[1:1]
	y0 = @lift ($x0_vec)[2:2]
	z0 = @lift ($x0_vec)[3:3]
	default_kwargs = (; minshaftlength=0, tiplength=2, strokemask=0.01, tipwidth=2, shaftwidth=1)

	
	arrows2d!(ax[1], x0, y0, x, y; default_kwargs..., kwargs...)
	arrows2d!(ax[2], y0, z0, y, z; default_kwargs..., kwargs...)
	arrows2d!(ax[3], x0, z0, x, z; default_kwargs..., kwargs...)
end


# ╔═╡ c27c54d5-fe4a-40a9-99d1-187b33a8b7ec


# ╔═╡ 564ca8cb-2074-43a6-8953-5af6c6eec1e7
function xyz_scatter!(ax, x_vec; kwargs...)
	x = @lift ($x_vec)[1:1]
	y = @lift ($x_vec)[2:2]
	z = @lift ($x_vec)[3:3]


	scatter!(ax[1],  x, y; kwargs...)
	scatter!(ax[2], y, z; kwargs...)
	scatter!(ax[3], x, z; kwargs...)
end


# ╔═╡ 48a9096d-64ba-46c3-800e-02ce98342740
acc_deviation = (acc_tot  .- acc_mw)

# ╔═╡ 22d61511-b22a-4938-94e2-2b98f383378a
traj_lmc_resampled = LilGuys.resample(traj_lmc, orbit.times)

# ╔═╡ f66b35ba-5c9b-4878-8625-8cb5b3b2a1f9
let
	fig = Figure()

	r_max = 150
	limits = tuple(fill((-r_max, r_max), 3)...)
	ax = LilGuys.axis_xyz(fig, limits=limits)

	d_lmc = Observable(acc_lmc[:, end] )
	d_reflex = Observable(acc_reflex[:, end])
	d_tot = Observable(acc_deviation[:, end])

	x_umi = Observable(orbit.positions[:, end])
	x_lmc = Observable(traj_lmc_resampled.positions[:, end])

	scale = 1e4
	xyz_arrows!(ax, x_umi, d_lmc, lengthscale=scale, color=COLORS[2])
	xyz_arrows!(ax, x_umi, d_reflex, lengthscale=scale, color=COLORS[3])
	xyz_arrows!(ax, x_umi, d_tot, lengthscale=scale)

	xyz_scatter!(ax, x_umi)
	xyz_scatter!(ax, x_lmc)
	xyz_scatter!(ax, Observable([0., 0, 0]), color=COLORS[3])

	record(fig, "force_animation.mp4", reverse(eachindex(times))) do i
		d_lmc[] = acc_lmc[:, i]
		d_reflex[] = acc_reflex[:, i]
		d_tot[] = acc_deviation[:, i]
		x_umi[] = orbit.positions[:, i]
		x_lmc[] = traj_lmc_resampled.positions[:, i]
	end
	
	fig
end

# ╔═╡ 9a1e95e0-2d5e-45ce-a365-62cf707c9cb2
L = LilGuys.angular_momenta(orbit)

# ╔═╡ 2189ba6c-9d20-4874-99c3-ee7b7ddfb8a1
function plot_E!(orbit)
	Φ = Agama.potential(pot_mw_static, orbit.positions, units)

	E =  Φ + 1/2 * speeds(orbit) .^2

	
	lines!(orbit.times * T2GYR, E)
end

# ╔═╡ f99f8f18-f60d-48ce-9df5-8eeb1a61dabc
let
	fig = Figure()
	ax = Axis(fig[1,1])

	plot_E!(orbit)
	plot_E!(orbit_no)
	plot_E!(orbit_noreflex)
	fig
end

# ╔═╡ 9ccaa9d4-d965-487a-8513-4cd3f43813d9
lines(times, L[3, :])

# ╔═╡ Cell order:
# ╠═b088fc88-9e6a-11f0-3095-3d9c1fe82a36
# ╠═91ab3185-a7a3-453d-ac72-9fd026d30e77
# ╠═ac195c2d-e80d-49d8-9524-9d440844f59d
# ╠═76bc084b-6714-47f8-9032-90baa407267f
# ╠═309abb2f-9c6a-4562-a7da-caebbed9bc03
# ╠═203f2d03-15f5-4553-9eeb-2699635c212c
# ╠═cea8a018-e315-4968-93b6-f2a59a1bcc44
# ╠═5678931c-6d87-42d8-8332-bf282c46dce8
# ╠═ad14e47b-c5d4-4aea-aff4-5a466ece1157
# ╠═59946955-acdc-41d8-ad67-d7843087f309
# ╠═dad0a7e3-db11-4bb3-b30f-3b68b8c1ce29
# ╠═6195c401-2dc7-425b-a082-8bbda16d0ba3
# ╠═2ccdac03-56fd-4f30-a237-b473435fe321
# ╠═d0cf64df-d46b-4776-9b46-d1c69a02cc61
# ╠═77503fde-9be2-4a7b-bed0-847297be2c76
# ╠═f7c3629d-4339-4e8f-b77b-f9c48ea0bbe4
# ╠═6aea4a07-3405-4e04-a4d4-b31598455043
# ╠═299038ee-d228-4eb9-9a97-708a95668fd2
# ╠═ad95c613-aa71-4b1c-b53b-ad6a9648af32
# ╠═a8622194-93a1-4db9-8ae4-b7dde2fd1cc7
# ╠═a1b205f3-3399-4f0a-96ce-c18f724ea2bd
# ╠═c2dd6977-228b-49c3-bc23-922a57f22c39
# ╠═1d61af93-6638-4c17-9f3b-1bff3d9be69b
# ╠═7ccad48a-7895-46dd-97ff-4e2758146cb3
# ╟─f8a15581-743f-406a-9a62-ebe14bf69ee5
# ╠═c7d41f48-4874-4524-b6c3-fd2c9e224b1b
# ╠═54f3841c-b5f9-4159-b11b-9695e2725650
# ╠═269cb10d-723a-4e57-ae95-0502c40e5040
# ╠═f7e1ed18-25fe-4664-8640-16ed7304200a
# ╠═afbf13a0-03ac-45e4-b529-3f30202b9107
# ╠═c27c54d5-fe4a-40a9-99d1-187b33a8b7ec
# ╠═564ca8cb-2074-43a6-8953-5af6c6eec1e7
# ╠═48a9096d-64ba-46c3-800e-02ce98342740
# ╠═22d61511-b22a-4938-94e2-2b98f383378a
# ╠═f66b35ba-5c9b-4878-8625-8cb5b3b2a1f9
# ╠═9a1e95e0-2d5e-45ce-a365-62cf707c9cb2
# ╠═2189ba6c-9d20-4874-99c3-ee7b7ddfb8a1
# ╠═f99f8f18-f60d-48ce-9df5-8eeb1a61dabc
# ╠═9ccaa9d4-d965-487a-8513-4cd3f43813d9
