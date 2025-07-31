### A Pluto.jl notebook ###
# v0.20.13

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
include("./paper_style.jl")

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ bb0ceaa6-b714-4164-8fde-922bdafcd273
r_peri = 50

# ╔═╡ 3261234e-d81d-487f-95c4-ab4aa6c43f48
import LinearAlgebra: ×, norm

# ╔═╡ 4e1fdfd3-8144-4510-9e4b-a8febf07d0ba
struct RotatingFrame <: LilGuys.Potential
	Ω::Vector{Float64}
end

# ╔═╡ a0e7d667-5fc6-4a37-a8ff-6494b438c1f1
struct ShiftedNFW <: LilGuys.Potential
	nfw::NFW
	x0::Vector{Float64}
end

# ╔═╡ ed916ba0-ea1c-4dd8-896e-588980f4c1f7
pot_halo = ShiftedNFW(LilGuys.NFW(M_s=79.5, r_s=20.2), [0., 0., 0.])

# ╔═╡ 803fa993-0b26-4638-ac87-7f1c7479cfd9
Ω = [0, 0, v_circ(pot_halo.nfw, r_peri) /  r_peri]

# ╔═╡ 69c7334c-afbd-48fd-b96e-874b7fb77e65
pot_eff = RotatingFrame(Ω)

# ╔═╡ 71bc271c-1930-45f9-aa42-c22375213a55
pot_dwarf = ShiftedNFW(LilGuys.NFW(v_circ_max=31 / V2KMS, r_circ_max=5.9), [r_peri, 0, 0])

# ╔═╡ e5f2285a-fc86-458a-abeb-b5b82c9f12e2
LilGuys.acceleration(pot::ShiftedNFW, x::Vector{Float64}, v=nothing, t=nothing) = - LilGuys.mass(pot.nfw, radii(x, pot.x0)) / radii(x, pot.x0)^3 .* (x .- pot.x0)

# ╔═╡ ff5f6a19-6b45-4b69-b111-ed1516fb2929
LilGuys.potential(pot::RotatingFrame, x) = - 1/2 * norm(pot.Ω × x)^2

# ╔═╡ 45253000-53bc-44c6-94b9-22eabc58e423
LilGuys.potential(pot::ShiftedNFW, x::Vector{Float64}) = LilGuys.potential(pot.nfw, radii(x .- pot.x0))

# ╔═╡ 186a7cdd-e090-4d7d-bd05-ad23cde01641
LilGuys.acceleration(pot::RotatingFrame, x, v, t) = 2 * pot.Ω × v - pot.Ω × (pot.Ω × x)

# ╔═╡ 8043edc9-21de-4675-be2d-0eb08d2fbf24
LilGuys.acceleration(pot_dwarf, [r_peri .+ 0.000000001, 0., 0])

# ╔═╡ 9e61adf6-6ef9-461f-8867-f0e7a7a0626a
pot_total = LilGuys.CompositePotential([pot_halo, pot_eff, pot_dwarf])

# ╔═╡ 3bccfe4b-0836-4edb-84e7-89dc521ff0b9
LilGuys.acceleration(pot_total, [r_peri, 0., 0.], [0., 0., 0.], 0)

# ╔═╡ c0f6c623-a978-45d0-be27-bb45b9c0902e
LilGuys.potential(pot_total, [r_peri, 0., 0.])

# ╔═╡ faa44b84-e693-4d40-8c5c-72b0599ef4d5
function integrate_isopotential(pot, initial=[-8,  0.]; x_direction=1, y_direction=2, s_scale=0.00001, h_scale=0.001, h0=0.0001)

	Φ(x) = LilGuys.potential(pot, x)
	∇Φ(x, v, t) = LilGuys.acceleration(pot, x, v, t)
	
	
	x = initial[1]
	y = initial[2]

	xs = [x]
	ys = [y]
	h = h0
    s = h0

	θ = atan(y, x)
	x_vec = zeros(3)
	x_vec[x_direction] = 1

	y_vec = zeros(3)
	y_vec[y_direction] = 1

    ρ_0 = Φ(x_vec * x .+ y_vec * y)
    dlρ_max = 0

    x0 = zeros(3)
    x0[x_direction] = x
    x0[y_direction] = y
    v0 = ∇Φ(x0, [0., 0., 0.], 0.)
	dx, dy, _ = v0

    s = s_scale * (sqrt(x^2 + y^2) / sqrt(dx^2 + dy^2) )
    h = h_scale * s

	for i in 1:100000
		x0 = zeros(3)
		x0[x_direction] = x
		x0[y_direction] = y
    	v0 = ∇Φ(x0, [0., 0., 0.], 0)
		dx, dy, _ = v0

		x += s .* dy
		y += -s .* dx

		push!(xs, x)
		push!(ys, y)

        s = s_scale * (sqrt(x^2 + y^2) / sqrt(dx^2 + dy^2) )
        h = h_scale * s

		θ_new = atan(y, x)
        dlρ_max = max(dlρ_max, abs(log10(-Φ(x0)) - log10(-ρ_0)))

		if θ_new > 0 && (θ < 0)
			break
		end

		θ = θ_new
	end
    @info "max log rel error = $dlρ_max"

    return xs, ys
end


# ╔═╡ e969d427-26cf-4f67-af7f-9636728dd771
L1 = LilGuys.find_zero(x->LilGuys.acceleration(pot_total, [x, 0, 0], [0., 0., 0.], 0)[1], r_peri - 2)

# ╔═╡ ac149774-9600-4299-bd64-c5a622306f59
L2 = LilGuys.find_zero(x->LilGuys.acceleration(pot_total, [x, 0, 0], [0., 0., 0.], 0)[1], r_peri + 2)

# ╔═╡ 5a8aff7b-b4ba-4498-8570-d7462f56e950
let
	fig = Figure()
	ax = Axis(fig[1,1])

	x = LinRange(0, 100, 1000)

	y = LilGuys.v_circ.(pot_halo.nfw, x)


	lines!(x, y .* V2KMS)

	fig

end

# ╔═╡ 4ece4c56-1081-4956-95b9-32bdc48529e5
let
	fig = Figure()
	ax = Axis(fig[1,1])

	x = LinRange(0, 100, 1000)

	y = [LilGuys.potential(pot_total, [xx, 0., 0.]) for xx in x]


	lines!(x, y )

	fig

end

# ╔═╡ 88699c63-9c18-4125-8d54-3cc9faf60bd7
let
	fig = Figure()
	ax = Axis(fig[1,1])

	x = LinRange(0, 100, 1000)

	y = [LilGuys.acceleration(pot_total, [xx, 0., 0.], zeros(3), 0)[1] for xx in x]


	lines!(x, y )

	vlines!(r_peri)
	fig

end

# ╔═╡ bc32151b-1762-4232-b931-9d1a9659e9e7
let
	fig = Figure()
	ax = Axis(fig[1,1])

	x = LinRange(0, 100, 1000)

	y = [LilGuys.potential(pot_total, [xx, 0., 0.]) for xx in x]


	lines!(x, y )

	vlines!(r_peri)
	fig

end

# ╔═╡ e923b834-5d2c-47be-a65b-e16eafa4419e
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=DataAspect())


	x, y = integrate_isopotential(pot_halo, [20, 0.])
	lines!(x, y)
	x, y = integrate_isopotential(pot_total, [L2* 0.99, 0.])
	lines!(x, y, linestyle=:solid)

	x, y = integrate_isopotential(pot_total, [L1 / 0.99, 0.])
	lines!(x, y, linestyle=:solid)
	x, y = integrate_isopotential(pot_total, [L2 - 2, 0.])
	lines!(x, y, linestyle=:solid)

	x, y = integrate_isopotential(pot_total, [L2 + 1, 0.])
	lines!(x, y, linestyle=:solid)
	fig

end

# ╔═╡ 19bbb402-40e7-4ac2-8d17-b945c4c18272
Φ0 = LilGuys.potential(pot_total, [L1, 0, 0])

# ╔═╡ 77fd177e-f277-4afd-b888-77aeb00ea4fa
r_J = LilGuys.find_zero(x -> LilGuys.density(pot_dwarf.nfw, x) - 3*LilGuys.density(pot_halo.nfw, r_peri), 2)

# ╔═╡ d3faf15b-cd8c-4338-90c9-e183e042e7aa
L1 - r_peri

# ╔═╡ 143d84a6-c6e2-4c26-abf0-c5df38a4bbaa
L2 - r_peri

# ╔═╡ 696c625b-c43b-4ab0-a55d-7650a94d0bf0
@savefig "lagrange_points" let
	fig = Figure()
	R = 10
	ax = Axis(fig[1,1], 
			  xlabel = L"x",
			  ylabel = L"y",
			  aspect=DataAspect(), limits=(r_peri-R, r_peri+R, -R, R))


	xs = LinRange(r_peri - R, r_peri + R, 1000)
	ys = LinRange(-R, R, 1000)

	dΦ = 0.03
	zs = [LilGuys.potential(pot_total, [x, y, 0]) for x in xs, y in ys]
	
	p = contourf!(xs, ys, zs, levels=LinRange(Φ0-dΦ, Φ0+dΦ, 50), colormap=:bluesreds, extendlow=:black)

	scatter!([L1, L2], zeros(2), color=:black)
	text!(L1, 0, text="L₁", align=(:center, -0.2))
	text!(L2, 0, text="L₂", align=(:center, -0.2))
	Colorbar(fig[1,2], p, label="effective potential")

	#arc!([r_peri, 0], r_J, 0, 2π, color=:black)

	arrows2d!([0.1], [0.5], [-0.075], [0.0], space=:relative)
	text!([0.1], [0.5], text="MW", space=:relative, align=(:center, 1.5))
	
	fig

end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═ed916ba0-ea1c-4dd8-896e-588980f4c1f7
# ╠═71bc271c-1930-45f9-aa42-c22375213a55
# ╠═69c7334c-afbd-48fd-b96e-874b7fb77e65
# ╠═bb0ceaa6-b714-4164-8fde-922bdafcd273
# ╠═803fa993-0b26-4638-ac87-7f1c7479cfd9
# ╠═e5f2285a-fc86-458a-abeb-b5b82c9f12e2
# ╠═45253000-53bc-44c6-94b9-22eabc58e423
# ╠═8043edc9-21de-4675-be2d-0eb08d2fbf24
# ╠═3261234e-d81d-487f-95c4-ab4aa6c43f48
# ╠═4e1fdfd3-8144-4510-9e4b-a8febf07d0ba
# ╠═a0e7d667-5fc6-4a37-a8ff-6494b438c1f1
# ╠═ff5f6a19-6b45-4b69-b111-ed1516fb2929
# ╠═186a7cdd-e090-4d7d-bd05-ad23cde01641
# ╠═3bccfe4b-0836-4edb-84e7-89dc521ff0b9
# ╠═c0f6c623-a978-45d0-be27-bb45b9c0902e
# ╠═9e61adf6-6ef9-461f-8867-f0e7a7a0626a
# ╠═faa44b84-e693-4d40-8c5c-72b0599ef4d5
# ╠═e969d427-26cf-4f67-af7f-9636728dd771
# ╠═ac149774-9600-4299-bd64-c5a622306f59
# ╠═5a8aff7b-b4ba-4498-8570-d7462f56e950
# ╠═4ece4c56-1081-4956-95b9-32bdc48529e5
# ╠═88699c63-9c18-4125-8d54-3cc9faf60bd7
# ╠═bc32151b-1762-4232-b931-9d1a9659e9e7
# ╠═e923b834-5d2c-47be-a65b-e16eafa4419e
# ╠═19bbb402-40e7-4ac2-8d17-b945c4c18272
# ╠═77fd177e-f277-4afd-b888-77aeb00ea4fa
# ╠═d3faf15b-cd8c-4338-90c9-e183e042e7aa
# ╠═143d84a6-c6e2-4c26-abf0-c5df38a4bbaa
# ╠═696c625b-c43b-4ab0-a55d-7650a94d0bf0
