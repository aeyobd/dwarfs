### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 156c47a4-39b1-11f1-861a-df30869de806
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using CairoMakie, Arya
end

# ╔═╡ 9f80a0a7-d83b-4a27-8038-dc375a82dc15
CairoMakie.activate!(type=:png)

# ╔═╡ fddd99d4-fd2f-4f90-bf04-8f8d20c39c48
dwarfs_dir = ENV["DWARFS_ROOT"]

# ╔═╡ cc264eb9-65df-4fc0-9be1-26142ecd2690
module Rapha
	include(joinpath(ENV["DWARFS_ROOT"], "utils/rapha_utils.jl"))
end

# ╔═╡ 18b526bf-f488-4606-85c2-8d9825bd1075
md"""
# Stellar energy
"""

# ╔═╡ 197aae18-1fe7-4d96-b7db-af2a2fd74a76
md"""
## Basic inputs (potential, stellar density, half-light radius, tidal evolution)
"""

# ╔═╡ 9ffa9087-4e74-44ad-be3f-4e663f8d69d4
halo_i = NFW(r_s=1, M_s=1)

# ╔═╡ 1ee59250-856b-4401-ac5b-f0694f9ed73c
r_max_rel = 0.1

# ╔═╡ 8e64262c-d9ff-4e6b-89c8-5291f60a5e73
α_star = 3; β_star = 3

# ╔═╡ 9c13853b-ccba-4ae6-af45-c7ea91ea823d
E_s = 0.2

# ╔═╡ a7c5d2e6-f6c9-40fe-8554-9e8c16560ead
r_min = 1e-4

# ╔═╡ 04f2c201-d973-42dc-901d-c3e5e67db3ae
r_max = 1e4

# ╔═╡ 31f63912-0b00-4227-bbd8-b9b67e561b03
md"""
# Utilities
"""

# ╔═╡ d11cd15e-f69e-405c-9b0d-2761f3d843e4
function density_from_f(f, Ψ, r)
	Ψ_0 = Ψ(r)
	# Eq 4.43 in BT2008
	
	return 4π * LilGuys.integrate(
		E -> f(E) * sqrt(2*(Ψ_0 - E)),
		0, Ψ_0
	)
end

# ╔═╡ 9a0fe901-ce83-498f-b798-c6bceaa3b82f
function density_of_states(Ψ, E)
	r_max = LilGuys.find_zero(r -> Ψ(r) + -E, (1e-5, 1e8))
	

	@assert 1 >= E > 0
	(4π)^2 * LilGuys.integrate(r-> r^2 * sqrt(2 * (-E + Ψ(r))), 
					   0, r_max)
end

# ╔═╡ 1cb35d11-6aca-4a83-8096-10e91bad45ab
struct InterpProfile <: LilGuys.SphericalProfile
	interp
end

# ╔═╡ 5c47e317-c697-47bd-a2e9-123d768e3cee
LilGuys.density(pot::InterpProfile, r) = pot.interp(r)

# ╔═╡ 173d168d-bf3a-4312-b959-e3e93be685bf
function calc_mass(ρ, r)
	4π * LilGuys.integrate(r -> ρ(r) * r^2, 0, r)
end

# ╔═╡ 90951b9c-4215-43d7-8ebb-0d10e228f82c
md"""
## Initial quantities
"""

# ╔═╡ 3ebe75bf-796e-436d-bbf3-a084bdad6d2d
logerp(x, y) = a -> 10^LilGuys.lerp(log10.(x), log10.(y))(log10(a))

# ╔═╡ 178a8242-f0e4-456f-9e37-f6d8fa6df07a
md"""
Rapha defines energy as a dimensionelsss version. 
Here I use ``\epsilon`` for the bindedness:

``
\epsilon = -\Phi - v^2 / 2 = -E
``

And I use ``x`` for rapha's dimensioness energy

``
x \equiv 1 - \frac{\epsilon}{|\Phi_0|} = 1 - \frac{E}{\Phi_0}
``

"""

# ╔═╡ 58133b1c-0e74-4569-882b-12213c5fbf55
md"""
### Dark matter potential and quantities
"""

# ╔═╡ dd96cf39-db76-4e60-86c0-02d9091e6d66
Ψ_i(r) = -LilGuys.potential(halo_i, r)

# ╔═╡ c0475038-3c18-4ef8-ac84-d9db1b4fc921
r_max_i = LilGuys.r_circ_max(halo_i)

# ╔═╡ 8417bb15-a637-4e93-8dfd-bdeb6a8f89b4
v_max_i = LilGuys.v_circ_max(halo_i)

# ╔═╡ 18716ef9-5426-4d3b-a0f9-924774f4bcb6
ϵ_max_i = Ψ_i(0)

# ╔═╡ 30106550-db04-4b1b-a96d-b0186aceb194
r_dm = logrange(r_min, r_max, 1000) |> collect

# ╔═╡ 8783310c-1f4d-4ba2-acf4-98f51aa149f9
function make_stellar_density(f_star, Ψ)
	r = r_dm
	ρ =  density_from_f.([f_star], Ψ, r)

	ρ[end-1:end] .= 0
	# ρ[1] = density_from_df.(df_stars, Φ, 0)
	return InterpProfile(LilGuys.lerp(r, ρ))
end

# ╔═╡ 899fd6ea-71ef-4584-a3c8-b42ca96b7989
ϵ_i = Ψ_i.(r_dm)

# ╔═╡ 7b5e779d-dcdb-439b-a1b7-668dfcef0472
x_i = 1 .- ϵ_i / ϵ_max_i

# ╔═╡ 2dee785a-a23c-4965-8709-bd910ddcbe81
f_dm_i_interp = LilGuys.DistributionFunction(LilGuys.density.(halo_i, r_dm), Ψ_i.(r_dm), r_dm)

# ╔═╡ 798bab2e-1a1d-45db-bd04-879282db74e2
f_dm_i = f_dm_i_interp.(ϵ_i)

# ╔═╡ 88ce15e6-189d-4bf9-8590-577d90fc2eac
g_i = density_of_states.(Ψ_i, ϵ_i)

# ╔═╡ c5cb0d96-f14e-4d48-abe8-916a020d1b4b
dNde_dm_i = f_dm_i .* g_i

# ╔═╡ 3c4cce88-6020-4bc3-ab15-d5bc8fc08d8a
A_dm_i = Rapha.normalize_dN_dE(x_i, dNde_dm_i)

# ╔═╡ 92930dbf-eb73-4066-989c-7c333bf7d438
md"""
## Final dark matter
"""

# ╔═╡ 545a5c7a-8180-4a3c-8b6e-7130e63d2a11
v_max_rel = LilGuys.v_circ_EN21(r_max_rel)

# ╔═╡ 03a5403c-91f8-4ead-8a9d-1c91a0c15d55
M_max_rel = r_max_rel * v_max_rel^2

# ╔═╡ e814f4b0-bb23-4cb0-a02d-3ba4311565db
halo_f = Rapha.ExpCusp(r_max_rel * r_max_i, v_max_rel * v_max_i)

# ╔═╡ b47aee4d-fe25-42b6-90d8-017e52f04a2c
Ψ_f(r) = -LilGuys.potential(halo_f, r)

# ╔═╡ 9d7840fd-14b4-477f-80ae-19a621aa7d68
ϵ_max_f = Ψ_f(r_min)

# ╔═╡ a5fb405e-051e-4063-a4e1-f249c6dbc6c8
md"""
## Initial stars
"""

# ╔═╡ ec59d79f-9da1-410b-8554-dbe3d883a1b7
α_exp = LilGuys.R_h(LilGuys.Exp3D())

# ╔═╡ 7c5a8005-13bd-4b85-b03d-2fd31dd98663
begin 
	dNde_star_i = Rapha.dN_dE_poly.(x_i, E_s, α_star, β_star)

	f_star_i = dNde_star_i ./ g_i


	# f_star_i = LilGuys.DistributionFunction(LilGuys.density.(prof_stars, r_dm), Ψ_i.(r_dm), r_dm).(ϵ_i)
	# dNde_star_i = f_star_i .* g_i
end

# ╔═╡ 5b52f2b2-17da-4632-a07e-2ebc202c2610
A_star_i = Rapha.normalize_dN_dE(x_i, dNde_star_i)

# ╔═╡ 27f342a3-ac37-4635-8da2-41c2dd9361f9
md"""
## Final stars
"""

# ╔═╡ d86322b1-0c08-41c8-a7fc-cd6dd30f5bcf
x_f, dNde_star_f = Rapha.final_stellar_energies(x_i, dNde_star_i, M_max_rel)

# ╔═╡ 8d605b3d-7c21-49e9-a3d3-87d62b4dbe46
ϵ_f = (1 .- x_f) .* (ϵ_max_f)

# ╔═╡ 1ca0d968-d37b-476e-b001-e8e1c81973de
f_dm_f = LilGuys.DistributionFunction(LilGuys.density.(halo_f, r_dm), Ψ_f.(r_dm), r_dm).(ϵ_f)

# ╔═╡ 78866f3e-1dbd-4e5f-bf86-ecde5eb3622c
g_f = density_of_states.(Ψ_f, ϵ_f)

# ╔═╡ 6ca92bea-a322-4591-98d0-3c617be7c5e2
dNde_dm_f = f_dm_f .* g_f

# ╔═╡ a6808503-dcda-43f7-91d9-445edbd537fb
A_dm_f = Rapha.normalize_dN_dE(x_f, dNde_dm_f)

# ╔═╡ 4954e68a-e697-40df-8224-75f6a9800ef2
f_star_f = dNde_star_f ./ g_f

# ╔═╡ ea1ed4d4-535c-46b2-a105-835bf9c90d16
A_star_f = Rapha.normalize_dN_dE(x_f, dNde_star_f)

# ╔═╡ 3e3e1293-2bb5-43f4-aae1-15959f43b9e4
md"""
# Plots
"""

# ╔═╡ 73647694-ea95-4b7e-947a-5b778076fb0d
snap = LilGuys.Snapshot(joinpath(ENV["DWARFS_ROOT"], "agama/halos/nfw_1e6_t20.0_xi3.0.hdf5"))

# ╔═╡ ad627754-8f28-4f0e-bf27-55c78c29fa74
ϵ_nbody = LilGuys.specific_energy(snap)

# ╔═╡ 973e9aa7-0409-49d7-a52c-140c53b06eb6
x_nbody = 1 .- ϵ_nbody ./ maximum(ϵ_nbody)

# ╔═╡ 317ca5e0-187b-49d4-a7f2-7ca09c41aacc
let
	f = lines(log10.(x_i), log10.(dNde_dm_i ./ A_dm_i ),
	 axis = (;
			 limits=(-1.5, 0, -4, 2))
	 )

	
	bins, val, _ = LilGuys.histogram(x_nbody, normalization=:pdf)

	x = midpoints(bins)
	A = Rapha.normalize_dN_dE(x, val)
	lines!(log10.(x), log10.(val ./ A))

	f
end

# ╔═╡ fe9d4e28-5072-4f77-a16c-4aa7b4e0b354
Rapha.log10normal(25, 25, 0.5)

# ╔═╡ a0a5bb2c-bec8-42d6-be82-6fbe70f61c0c
import Distributions: pdf, LogNormal

# ╔═╡ 2863e921-6d96-478a-844a-0c888489e2be
f_star_i_interp = LilGuys.lerp(ϵ_i, f_star_i)

# ╔═╡ b8b9d995-d623-4097-9ae7-05391bead2d0
ρ_stars(r) = density_from_f(f_star_i_interp, Ψ_i, r)

# ╔═╡ 2539649e-c016-4376-af63-7ca151fbf2b7
R_end = r_max

# ╔═╡ 4b568a40-1896-4f89-8c4a-33ce82ae3655
function calc_R_h(Σ)
	
	M(r) =  2π * LilGuys.integrate(r -> Σ(r) * r, 0, r)

	Mtot = M(R_end)
	return LilGuys.find_zero(r -> M(r) - 1/2 * Mtot, (1e-8, R_end))
end

# ╔═╡ 8dc820ba-95ab-44aa-9128-5360681c330e
function calc_r_h(ρ)
	
	M(r) =  4π * LilGuys.integrate(r -> ρ(r) * r^2, 0, r)

	Mtot = M(R_end)
	return LilGuys.find_zero(r -> M(r) - 1/2 * Mtot, (1e-8, R_end))
end

# ╔═╡ fade313d-63f5-42ae-895d-6823eb449a9b
function surface_density_from_density(profile::InterpProfile, )
	Σ(R) = surface_density_from_density(profile, R)
	R = logrange(1e-3, 1e3, 10_000)
	f = LilGuys.lerp(log10.(R), log10.(Σ.(R)))
	
	return R -> exp10(f(log10(R)))
end

# ╔═╡ 2491008c-7cd4-4bc2-a8e1-9507b64036d2
function surface_density_from_density(profile::InterpProfile, R::Real)
    integrand(r) = LilGuys.density(profile, r) * r / sqrt(r^2 - R^2)
    return 2*LilGuys.integrate(integrand, R * (1 + 1e-8), Inf)
end

# ╔═╡ 96f9f16d-323f-44d2-b707-c1e0dec84465
function calc_R_h(df_stars, Φ)
	ρ = make_stellar_density(df_stars, Φ)
	Σ = surface_density_from_density(ρ)
	calc_R_h(Σ)	
end

# ╔═╡ 45e3ee5a-be28-400a-b369-3a858c71b74a
ρ = make_stellar_density(f_star_i_interp, Ψ_i)

# ╔═╡ e2c60d16-debe-4367-a5a9-7ebdc739d351
M = 4π * LilGuys.integrate(r -> LilGuys.density(ρ, r) * r^2, 0, Inf)

# ╔═╡ 845d6da0-af32-42b9-85f0-342f6b2dd0e0
Σ = surface_density_from_density(ρ)

# ╔═╡ f43d58b1-414a-430d-a3ab-7ea5cb590632
R_h = calc_R_h(Σ)

# ╔═╡ 84385b77-a15d-47e3-8ec1-151c8b6e06ab
R_h / r_max

# ╔═╡ 37666df8-7bb6-44d3-bfd5-88b14067b1be
prof_stars = LilGuys.Exp3D(r_s= R_h / α_exp)

# ╔═╡ 840ac8b7-a4d9-4865-89e2-1e1525a8d4fe
2π * LilGuys.integrate(r -> Σ(r) * r, 0, R_end)

# ╔═╡ cf3c99f3-190c-47c6-af64-e370f27ff53f
4π * LilGuys.integrate(r -> LilGuys.density(ρ, r) * r^2, 0, 1)

# ╔═╡ 78588c60-9dd3-4235-be22-b4dde9df6d8d
let
	fig = Figure()
	ax = Axis(fig[1,1],
			  xlabel = "log r",
			  ylabel = "log ρ"
)

	r = logrange(1e-3, 1e2, 1000)

	rho = Σ.(r)

	lines!(log10.(r), log10.(rho), label="model")
	lines!(log10.(r), log10.(
		LilGuys.surface_density.(prof_stars, r) .*M), linestyle=:dot, label="analytic")

	axislegend(position=:lb)
	xlims!(-6, 2)
	ylims!(-5, 0)
	fig

end

# ╔═╡ 0d07d47c-1481-4ee6-a246-c3b2172985b7
calc_r_h(ρ.interp)

# ╔═╡ 5a83951f-cd24-415c-82ee-c3ff1328f78c
calc_R_h(Σ)

# ╔═╡ 0498aba6-92d0-4819-b4fc-aa245fd0493b
R_h

# ╔═╡ 3b5f6519-9d9a-42f1-b7b7-bcdb8a7b1367
LilGuys.density(prof_stars, 0.5)

# ╔═╡ 260c0f50-ac2b-4fc3-bfc1-9183d6408560
let
	fig = Figure()
	ax = Axis(fig[1,1],
)

	r = logrange(1e-6, 1e2, 1000)

	rho = Σ.(r)

	lines!(log10.(r), (rho) .* r)
	lines!(log10.(r), LilGuys.surface_density.(prof_stars, r) .*M .* r, linestyle=:dot)
	
	xlims!(-6, 2)
	fig

end

# ╔═╡ 04b191db-dc07-41a3-85ae-bc58ad929758
x_max_t = 0.77 * M_max_rel^0.43

# ╔═╡ 55e34425-a89d-47a3-902e-0e8cb83e538a
Rapha.x_t_map(0.0002, x_max_t)

# ╔═╡ 0c263f82-82c3-4951-9882-6d1fad61748b
0.00002 ./ x_max_t

# ╔═╡ f3a8ea5d-9fdc-4574-8b8b-72e9f86544ca
x_max_t

# ╔═╡ 00a21155-0df4-4e87-b94c-44f3f0ca8d10
dNde_star_t = dNde_star_i .* Rapha.dNde_map.(x_i, x_max_t)

# ╔═╡ 4b2eab47-d421-472e-8c05-7cd476b8294f
let
	f = lines(log10.(x_i), log10.(dNde_dm_i ./ A_dm_i ),
	 axis = (;
			 limits=(-1.5, 0, -4, 0))
	 )


	lines!(log10.(x_i), log10.(dNde_dm_i ./ A_dm_i .* Rapha.dNde_map.(x_i, x_max_t)), color=COLORS[1])

	vlines!(log10.(x_max_t))

	A_extra = 1000
	vlines!(log10.(E_s))
	lines!(log10.(x_i), log10.(dNde_star_i ./ A_star_i ./ A_extra), color=COLORS[2])
	lines!(log10.(x_i), log10.(dNde_star_t ./ A_star_i ./ A_extra), linestyle=:dot)
	f
end

# ╔═╡ 69bae1ac-c9a6-4433-aab4-dfd307c93f1d
let
	f = lines(log10.(x_i), log10.(dNde_dm_i ./ A_dm_i ),
	 axis = (;
			 limits=(-1.5, 0, -4, 0))
	 )


	lines!(log10.(x_f), log10.(dNde_dm_f))


	A_extra = 1000
	vlines!(log10.(E_s))
	lines!(log10.(x_i), log10.(dNde_star_i ./ A_star_i ./A_extra))
	lines!(log10.(x_i), log10.(dNde_star_t ./ A_star_i ./A_extra))

  	lines!(log10.(x_i), log10.(dNde_star_f ./ A_star_f ./ A_extra), linestyle=:dot)
	f
end

# ╔═╡ 8843c1ab-2e5e-4bb4-960b-3fa7bd0730f5
let
	fig = Figure()
	ax = Axis(fig[1,1])

A_extra = 1
	
	
	lines!(log10.(x_f), (dNde_star_t ./ A_star_i ./A_extra))

	
	dNde_star_t_interp = LilGuys.lerp(x_i, dNde_star_t)
	dNde_star_conv(x_f) = LilGuys.integrate(
		xx_i -> dNde_star_t_interp(xx_i) * 
			Rapha.log10normal(x_f, Rapha.x_t_map(xx_i, x_max_t), 0.02),
		1e-3, 1
	)

	
	lines!(log10.(x_f), (dNde_star_conv.(x_f) ./ A_star_f ./ A_extra ./ 2), linestyle=:dash)

	
	lines!(log10.(x_f), (dNde_star_f ./ A_star_f ./ A_extra ./ 2), linestyle=:dot)

	xlims!(-1, 0)
	fig
end

# ╔═╡ 16ae753c-24a7-4aad-bbc7-c50f868dedaa
let
	fig = Figure()
	ax = Axis(fig[1,1],
			  limits = (-2, 0.5, -2, 0),
			  aspect = 1
			 )
	
	lines!(log10.(x_i ./ x_max_t), log10.(x_f))

	fig

end

# ╔═╡ be72946c-ce17-4f0e-9aa6-4aabf26ab57f
md"""
# Tidal tracks
"""

# ╔═╡ f7f0da17-bb98-4d68-8beb-4d9130399247
r_max_i

# ╔═╡ 7f54f778-56a9-4577-80b9-0b5e95e7b526
function final_ρ_dm(r_max_rel)
	v_max_rel = LilGuys.v_circ_EN21(r_max_rel)
	M_max_rel = r_max_rel * v_max_rel^2
	halo_f = Rapha.ExpCusp(r_max_rel * r_max_i, v_max_rel * v_max_i)
end

# ╔═╡ eee6e249-93a3-44e7-8e8f-2e5dbd16a829
LilGuys.r_circ_max(final_ρ_dm(1))

# ╔═╡ 16b8c1a4-1486-4de0-a0b3-1b055dd34ed7
function final_ρ_star(r_max_rel)
	v_max_rel = LilGuys.v_circ_EN21(r_max_rel)
	M_max_rel = r_max_rel * v_max_rel^2
	halo_f = Rapha.ExpCusp(r_max_rel * r_max_i, v_max_rel * v_max_i)
	Ψ_f(r) = -LilGuys.potential(halo_f, r)
	ϵ_max_f = Ψ_f(r_min)

	
	x_f, dNde_star_f = Rapha.final_stellar_energies(x_i, dNde_star_i, M_max_rel)
	ϵ_f = (1 .- x_f) .* (ϵ_max_f)
	g_f = density_of_states.(Ψ_f, ϵ_f)

	f_star_f = dNde_star_f ./ g_f

	return make_stellar_density(LilGuys.lerp(ϵ_f, f_star_f), Ψ_f)
end


# ╔═╡ 4d3e85f8-8c76-4db5-a321-43cbfbea9590
let
	fig = Figure()
	ax = Axis(fig[1,1])

	for r_max_rel in [1, 0.3, 0.1, 0.03, 0.01]
		h = final_ρ_dm(r_max_rel)
		
		prof = final_ρ_star(r_max_rel)

		x = LinRange(-3, 1, 100)
		r = 10 .^ x
		@info prof.interp.x

		lines!(x, log10.(LilGuys.density.(h, r)) .- 4, linestyle=:dash)
		lines!(x, log10.(LilGuys.density.(prof, r)))
	end

	ylims!(-12, 0)

	fig
end

	

# ╔═╡ c0e45c43-ca27-4e8c-b547-f490fedc6476
ρ2 = final_ρ_star(1)

# ╔═╡ e927ba50-0a89-4ee4-b8b4-0957f954ca84
let
	fig = Figure()
	ax = Axis(fig[1,1],
)

	r = r_dm

	rho = ρ_stars.(r)

	lines!(log10.(r), (rho))
	# lines!(log10.(r), LilGuys.density.(ρ, r), linestyle=:dot)
	lines!(log10.(r), LilGuys.density.(prof_stars, r) .*M, linestyle=:dot)


	lines!(log10.(r), ρ2.interp.(r) .* 2.2, linestyle=:dash)
	xlims!(-6, 2)
	fig

end

# ╔═╡ 23a25536-4168-4f72-856d-a35bd1f59265
calc_r_h(ρ2.interp)

# ╔═╡ e722fb47-34d7-45ca-9b6c-6cb9fe3ee5ca
function final_properties(halo, stellar)
	r_h = calc_r_h(stellar.interp)
	σv = LilGuys.σv_star_mean(halo, stellar)
	M = calc_mass(stellar.interp, R_end)
	M, r_h, σv
	
end

# ╔═╡ 94a8c720-0051-4d7d-9db9-98a09e66df18
let
	fig = Figure()
	ax = Axis(fig[1,1])

	x = LinRange(-1.5, 1, 1000)

	lines!(x, log10.(LilGuys.v_circ.(halo_i, r_max_i .* 10 .^ x) ./ v_max_i))

	x, y = LilGuys.EN21_tidal_track(1, 1)
	lines!(log10.(x), log10.(y))

	for r_max_rel in logrange(1, 0.01, 100)
		h = final_ρ_dm(r_max_rel)
		
		prof = final_ρ_star(r_max_rel)

		M, r_h, σv = final_properties(h, prof)

		scatter!(log10.(r_h ./ r_max_i), 
				 log10.(LilGuys.v_circ(h, r_h) ./ v_max_i), color=log10.(r_max_rel), colorrange=(-2, 0))
	end

	fig
end

# ╔═╡ 066adc4c-cb74-496c-97c2-8cdca4f8a5d3


# ╔═╡ Cell order:
# ╠═156c47a4-39b1-11f1-861a-df30869de806
# ╠═9f80a0a7-d83b-4a27-8038-dc375a82dc15
# ╠═fddd99d4-fd2f-4f90-bf04-8f8d20c39c48
# ╠═cc264eb9-65df-4fc0-9be1-26142ecd2690
# ╟─18b526bf-f488-4606-85c2-8d9825bd1075
# ╠═197aae18-1fe7-4d96-b7db-af2a2fd74a76
# ╠═9ffa9087-4e74-44ad-be3f-4e663f8d69d4
# ╠═1ee59250-856b-4401-ac5b-f0694f9ed73c
# ╠═8e64262c-d9ff-4e6b-89c8-5291f60a5e73
# ╠═9c13853b-ccba-4ae6-af45-c7ea91ea823d
# ╠═a7c5d2e6-f6c9-40fe-8554-9e8c16560ead
# ╠═04f2c201-d973-42dc-901d-c3e5e67db3ae
# ╠═84385b77-a15d-47e3-8ec1-151c8b6e06ab
# ╟─31f63912-0b00-4227-bbd8-b9b67e561b03
# ╠═d11cd15e-f69e-405c-9b0d-2761f3d843e4
# ╠═8783310c-1f4d-4ba2-acf4-98f51aa149f9
# ╠═9a0fe901-ce83-498f-b798-c6bceaa3b82f
# ╠═1cb35d11-6aca-4a83-8096-10e91bad45ab
# ╠═5c47e317-c697-47bd-a2e9-123d768e3cee
# ╠═4b568a40-1896-4f89-8c4a-33ce82ae3655
# ╠═173d168d-bf3a-4312-b959-e3e93be685bf
# ╠═8dc820ba-95ab-44aa-9128-5360681c330e
# ╠═96f9f16d-323f-44d2-b707-c1e0dec84465
# ╟─90951b9c-4215-43d7-8ebb-0d10e228f82c
# ╠═3ebe75bf-796e-436d-bbf3-a084bdad6d2d
# ╟─178a8242-f0e4-456f-9e37-f6d8fa6df07a
# ╟─58133b1c-0e74-4569-882b-12213c5fbf55
# ╠═dd96cf39-db76-4e60-86c0-02d9091e6d66
# ╠═c0475038-3c18-4ef8-ac84-d9db1b4fc921
# ╠═8417bb15-a637-4e93-8dfd-bdeb6a8f89b4
# ╠═18716ef9-5426-4d3b-a0f9-924774f4bcb6
# ╠═30106550-db04-4b1b-a96d-b0186aceb194
# ╠═899fd6ea-71ef-4584-a3c8-b42ca96b7989
# ╠═7b5e779d-dcdb-439b-a1b7-668dfcef0472
# ╠═2dee785a-a23c-4965-8709-bd910ddcbe81
# ╠═798bab2e-1a1d-45db-bd04-879282db74e2
# ╠═88ce15e6-189d-4bf9-8590-577d90fc2eac
# ╠═c5cb0d96-f14e-4d48-abe8-916a020d1b4b
# ╠═3c4cce88-6020-4bc3-ab15-d5bc8fc08d8a
# ╟─92930dbf-eb73-4066-989c-7c333bf7d438
# ╠═545a5c7a-8180-4a3c-8b6e-7130e63d2a11
# ╠═03a5403c-91f8-4ead-8a9d-1c91a0c15d55
# ╠═e814f4b0-bb23-4cb0-a02d-3ba4311565db
# ╠═b47aee4d-fe25-42b6-90d8-017e52f04a2c
# ╠═9d7840fd-14b4-477f-80ae-19a621aa7d68
# ╠═1ca0d968-d37b-476e-b001-e8e1c81973de
# ╠═78866f3e-1dbd-4e5f-bf86-ecde5eb3622c
# ╠═6ca92bea-a322-4591-98d0-3c617be7c5e2
# ╠═a6808503-dcda-43f7-91d9-445edbd537fb
# ╟─a5fb405e-051e-4063-a4e1-f249c6dbc6c8
# ╠═37666df8-7bb6-44d3-bfd5-88b14067b1be
# ╠═ec59d79f-9da1-410b-8554-dbe3d883a1b7
# ╠═7c5a8005-13bd-4b85-b03d-2fd31dd98663
# ╠═5b52f2b2-17da-4632-a07e-2ebc202c2610
# ╟─27f342a3-ac37-4635-8da2-41c2dd9361f9
# ╠═d86322b1-0c08-41c8-a7fc-cd6dd30f5bcf
# ╠═8d605b3d-7c21-49e9-a3d3-87d62b4dbe46
# ╠═4954e68a-e697-40df-8224-75f6a9800ef2
# ╠═ea1ed4d4-535c-46b2-a105-835bf9c90d16
# ╟─3e3e1293-2bb5-43f4-aae1-15959f43b9e4
# ╠═73647694-ea95-4b7e-947a-5b778076fb0d
# ╠═ad627754-8f28-4f0e-bf27-55c78c29fa74
# ╠═973e9aa7-0409-49d7-a52c-140c53b06eb6
# ╠═317ca5e0-187b-49d4-a7f2-7ca09c41aacc
# ╠═4b2eab47-d421-472e-8c05-7cd476b8294f
# ╠═69bae1ac-c9a6-4433-aab4-dfd307c93f1d
# ╠═8843c1ab-2e5e-4bb4-960b-3fa7bd0730f5
# ╠═55e34425-a89d-47a3-902e-0e8cb83e538a
# ╠═0c263f82-82c3-4951-9882-6d1fad61748b
# ╠═f3a8ea5d-9fdc-4574-8b8b-72e9f86544ca
# ╠═fe9d4e28-5072-4f77-a16c-4aa7b4e0b354
# ╠═a0a5bb2c-bec8-42d6-be82-6fbe70f61c0c
# ╠═2863e921-6d96-478a-844a-0c888489e2be
# ╠═b8b9d995-d623-4097-9ae7-05391bead2d0
# ╠═f43d58b1-414a-430d-a3ab-7ea5cb590632
# ╠═00a21155-0df4-4e87-b94c-44f3f0ca8d10
# ╠═2539649e-c016-4376-af63-7ca151fbf2b7
# ╠═fade313d-63f5-42ae-895d-6823eb449a9b
# ╠═2491008c-7cd4-4bc2-a8e1-9507b64036d2
# ╠═45e3ee5a-be28-400a-b369-3a858c71b74a
# ╠═e2c60d16-debe-4367-a5a9-7ebdc739d351
# ╠═840ac8b7-a4d9-4865-89e2-1e1525a8d4fe
# ╠═845d6da0-af32-42b9-85f0-342f6b2dd0e0
# ╠═cf3c99f3-190c-47c6-af64-e370f27ff53f
# ╠═78588c60-9dd3-4235-be22-b4dde9df6d8d
# ╠═e927ba50-0a89-4ee4-b8b4-0957f954ca84
# ╠═0d07d47c-1481-4ee6-a246-c3b2172985b7
# ╠═23a25536-4168-4f72-856d-a35bd1f59265
# ╠═5a83951f-cd24-415c-82ee-c3ff1328f78c
# ╠═0498aba6-92d0-4819-b4fc-aa245fd0493b
# ╠═3b5f6519-9d9a-42f1-b7b7-bcdb8a7b1367
# ╠═260c0f50-ac2b-4fc3-bfc1-9183d6408560
# ╠═04b191db-dc07-41a3-85ae-bc58ad929758
# ╠═16ae753c-24a7-4aad-bbc7-c50f868dedaa
# ╠═be72946c-ce17-4f0e-9aa6-4aabf26ab57f
# ╠═4d3e85f8-8c76-4db5-a321-43cbfbea9590
# ╠═94a8c720-0051-4d7d-9db9-98a09e66df18
# ╠═eee6e249-93a3-44e7-8e8f-2e5dbd16a829
# ╠═f7f0da17-bb98-4d68-8beb-4d9130399247
# ╠═7f54f778-56a9-4577-80b9-0b5e95e7b526
# ╠═16b8c1a4-1486-4de0-a0b3-1b055dd34ed7
# ╠═c0e45c43-ca27-4e8c-b547-f490fedc6476
# ╠═e722fb47-34d7-45ca-9b6c-6cb9fe3ee5ca
# ╠═066adc4c-cb74-496c-97c2-8cdca4f8a5d3
