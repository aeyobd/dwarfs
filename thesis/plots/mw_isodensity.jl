### A Pluto.jl notebook ###
# v0.20.19

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

# ╔═╡ 48655feb-9ccc-4c61-a0b7-1c8c816be4fd
using Agama

# ╔═╡ 7ea4d01f-31f3-4822-bed3-2428789e2914
using CSV, DataFrames

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ bf905cc5-50e7-4789-abe1-95e74ad93ecd
potname = "EP2020"

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ a0dfb6b2-31a9-4ac5-b935-4375385bc16e
potfile = Dict(
	"EP2020" => "EP2020.ini",
	"L3M11" => "vasiliev24/L3M11/potential_mw_init.ini"
)[potname]

# ╔═╡ 6b9f554c-dfb4-42f5-9cfc-d51d4a696b2b
pot = Agama.Potential(file=joinpath(ENV["DWARFS_ROOT"], "agama/potentials", potfile))

# ╔═╡ 7717e4e9-0231-4132-a5c2-6da577090033
X_SUN = [-8.122, 0., 0.]

# ╔═╡ 5161ee09-2651-47fd-b703-05bfb9ef13a0
function integrate_isodensity_2d(pot, initial=[-X_SUN[1], 0.]; s_scale=0.0003, h_scale=0.0001, h0=0.0001, kwargs...)
	x = initial[1]
	y = initial[2]

	xs = [x]
	ys = [y]
	h = h0
    s = h0

	θ = atan(y, x)
	x_vec = [1, 0]

	y_vec = [0, 1]
	ρ(x) = pot._py.projectedDensity(x; kwargs...) |> py2f
	
    ρ_0 = ρ(x_vec * x .+ y_vec * y)
    dlρ_max = 0

    x0 = [x, y]
    dx = (ρ(x0 .+ x_vec * h) - ρ(x0))/h
    dy = (ρ(x0 .+ y_vec * h) - ρ(x0))/h

    s = s_scale * (sqrt(x^2 + y^2) / sqrt(dx^2 + dy^2) )
    h = h_scale * s

	for i in 1:100000
		x0 = [x, y]
		dx = (ρ(x0 .+ x_vec * h) - ρ(x0))/h
		dy = (ρ(x0 .+ y_vec * h) - ρ(x0))/h
		x += s .* dy
		y += -s .* dx

		push!(xs, x)
		push!(ys, y)

        s = s_scale * (sqrt(x^2 + y^2) / sqrt(dx^2 + dy^2) )
        h = h_scale * s

		θ_new = atan(y, x)
        dlρ_max = max(dlρ_max, abs(log10(ρ(x0)) - log10(ρ_0)))

		if θ_new > 0 && (θ < 0)
			break
		end

		θ = θ_new
	end
    @info "max log rel error = $dlρ_max"

    return xs, ys
end


# ╔═╡ e6636b6d-57fc-4463-ad2a-1fe0cf7abe52
isodensity_xy_2d = integrate_isodensity_2d(pot)

# ╔═╡ 2858a1a4-964f-46bc-b323-5a1fa9b20994
isodensity_xz_2d = integrate_isodensity_2d(pot, beta = π/2)

# ╔═╡ 65ab78a2-75d3-4ac8-b69e-358e8d368756
function integrate_isodensity(pot, initial=[-X_SUN[1], 0.]; x_direction=2, y_direction=3, s_scale=0.0001, h_scale=0.0001, h0=0.0001)
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
	ρ(x) = Agama.density(pot, x)
    ρ_0 = ρ(x_vec * x .+ y_vec * y)
    dlρ_max = 0

    x0 = zeros(3)
    x0[x_direction] = x
    x0[y_direction] = y
    dx = (ρ(x0 .+ x_vec * h) - ρ(x0))/h
    dy = (ρ(x0 .+ y_vec * h) - ρ(x0))/h

    s = s_scale * (sqrt(x^2 + y^2) / sqrt(dx^2 + dy^2) )
    h = h_scale * s

	for i in 1:100000
		x0 = zeros(3)
		x0[x_direction] = x
		x0[y_direction] = y
		dx = (ρ(x0 .+ x_vec * h) - ρ(x0))/h
		dy = (ρ(x0 .+ y_vec * h) - ρ(x0))/h
		x += s .* dy
		y += -s .* dx

		push!(xs, x)
		push!(ys, y)

        s = s_scale * (sqrt(x^2 + y^2) / sqrt(dx^2 + dy^2) )
        h = h_scale * s

		θ_new = atan(y, x)
        dlρ_max = max(dlρ_max, abs(log10(ρ(x0)) - log10(ρ_0)))

		if θ_new > 0 && (θ < 0)
			break
		end

		θ = θ_new
	end
    @info "max log rel error = $dlρ_max"

    return xs, ys
end


# ╔═╡ 1f4e145a-0a68-4f62-b9f0-ecb4a78a46c3
isodensity_yz = integrate_isodensity(pot)

# ╔═╡ 00845d5f-0f2e-4638-b39e-e5e7d8ca4ea4
isodensity_xy = integrate_isodensity(pot, x_direction=1, y_direction=2)

# ╔═╡ d875656e-5c8b-46d6-8bdb-3d6612b5c9e0
isodensity_xz = integrate_isodensity(pot, x_direction=1, y_direction=3)

# ╔═╡ 222b1789-7f92-4e42-a3dc-6e38e77ef19e
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=DataAspect())

	lines!(isodensity_xy...)

	fig

end

# ╔═╡ 5f8d41bd-64d5-4e51-bd7e-bf30a33ea8fe
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=DataAspect())

	lines!(isodensity_xy_2d...)

	fig

end

# ╔═╡ bcf10081-03a5-476f-a4cc-3e4e3ed34665
let
	fig = Figure()
	ax = Axis(fig[1,1], aspect=DataAspect())

	poly!(isodensity_xz...)

	lines!(isodensity_xz_2d..., color=COLORS[2])
	fig

end

# ╔═╡ a6b4ace3-275a-4443-ac3e-d8191cb29f80
let
	fig = Figure()
	ax = Axis(fig[1,1], autolimitaspect=1)

	poly!(isodensity_yz...)
	lines!(isodensity_xz_2d..., color=COLORS[2])

	fig

end

# ╔═╡ cf0ea249-4d48-482b-b850-baa33cf9644f
lines(Agama.density(pot, hcat(isodensity_xy..., zeros(length(isodensity_xy[1])))'))

# ╔═╡ 148a5569-c4bc-46d9-9f36-b740adc8d988
lines(Agama.density(pot, hcat(isodensity_xz[1], zeros(length(isodensity_xz[1])), isodensity_xz[2])'))

# ╔═╡ ffb38d9e-c4de-4787-b0f1-03f0daa038a1
lines(Agama.density(pot, hcat(zeros(length(isodensity_yz[1])), isodensity_yz..., )'))

# ╔═╡ 63be0abf-0fd0-4f5b-b779-ddb2a52f0b77
CSV.write("resources/$(potname)_iso_xy.csv", DataFrame("x"=>isodensity_xy_2d[1], "y"=>isodensity_xy_2d[2]))

# ╔═╡ e1d3fa27-77f0-4ad4-8504-990848886a3d
CSV.write("resources/$(potname)_iso_yz.csv", DataFrame("y"=>isodensity_xz_2d[1], "z"=>isodensity_xz_2d[2]))

# ╔═╡ 76ba4985-4bef-4965-b572-7593427e469e
CSV.write("resources/$(potname)_iso_xz.csv", DataFrame("x"=>isodensity_xz_2d[1], "z"=>isodensity_xz_2d[2]))

# ╔═╡ Cell order:
# ╠═bf905cc5-50e7-4789-abe1-95e74ad93ecd
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═48655feb-9ccc-4c61-a0b7-1c8c816be4fd
# ╠═a0dfb6b2-31a9-4ac5-b935-4375385bc16e
# ╠═6b9f554c-dfb4-42f5-9cfc-d51d4a696b2b
# ╠═7717e4e9-0231-4132-a5c2-6da577090033
# ╠═e6636b6d-57fc-4463-ad2a-1fe0cf7abe52
# ╠═2858a1a4-964f-46bc-b323-5a1fa9b20994
# ╠═5161ee09-2651-47fd-b703-05bfb9ef13a0
# ╠═65ab78a2-75d3-4ac8-b69e-358e8d368756
# ╠═1f4e145a-0a68-4f62-b9f0-ecb4a78a46c3
# ╠═00845d5f-0f2e-4638-b39e-e5e7d8ca4ea4
# ╠═d875656e-5c8b-46d6-8bdb-3d6612b5c9e0
# ╠═222b1789-7f92-4e42-a3dc-6e38e77ef19e
# ╠═5f8d41bd-64d5-4e51-bd7e-bf30a33ea8fe
# ╠═bcf10081-03a5-476f-a4cc-3e4e3ed34665
# ╠═a6b4ace3-275a-4443-ac3e-d8191cb29f80
# ╠═cf0ea249-4d48-482b-b850-baa33cf9644f
# ╠═148a5569-c4bc-46d9-9f36-b740adc8d988
# ╠═ffb38d9e-c4de-4787-b0f1-03f0daa038a1
# ╠═7ea4d01f-31f3-4822-bed3-2428789e2914
# ╠═63be0abf-0fd0-4f5b-b779-ddb2a52f0b77
# ╠═e1d3fa27-77f0-4ad4-8504-990848886a3d
# ╠═76ba4985-4bef-4965-b572-7593427e469e
