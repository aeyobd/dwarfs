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

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ ca9e7089-6617-4bff-996b-b16de0587c5b
Mstar = 3e-4

# ╔═╡ 1b202d9d-92ab-4798-8851-f3d8998b5c5f
stars = LilGuys.read_hdf5_table(joinpath(ENV["DWARFS_ROOT"], "analysis/sculptor/1e6_v48_r7_oblate_0.5/stars/exp2d_rs0.20/probabilities_stars.hdf5"))

# ╔═╡ 4cfdc774-f065-4d71-bdd9-ffd26dbe64a7
out = Output(joinpath(ENV["DWARFS_ROOT"], "analysis/sculptor/1e6_v48_r7_oblate_0.5/orbit_smallperi"), weights=stars.probability * Mstar)

# ╔═╡ e9dcf228-7c2d-4deb-bbe5-b9215512bc70
function moment_inertia_tensor(snap)
	snap_bound = snap[LilGuys.bound_particles(snap)]

	I = zeros(3, 3)

	x = snap_bound.positions .- snap_bound.x_cen
	r = radii(x)

	for i in 1:3
		for j in i:3
			I[i, j] = sum(x[i, :] .* x[j, :] ./ r .^ 3)
		end
	end

	I[2, 1] = I[1,2]
	I[3, 1] = I[1, 3]
	I[3, 2] = I[2, 3]
	return I
end

# ╔═╡ 3d2f16ce-3d97-429a-9acb-004dc5162575
import LinearAlgebra: eigen

# ╔═╡ c554e3d1-7ccc-4075-86c1-cf1df69afaf6
function principal_axes(snap)
	vecs =  eigen(moment_inertia_tensor(snap)).vectors

	return vecs
end

# ╔═╡ 6f24e979-b2ce-4590-b8d1-d442924ab2fd
function project_principal(snap)
	eig = principal_axes(snap)
	
	x = snap.positions .- snap.x_cen

	return eig' * x
end

# ╔═╡ 9a281d8f-748d-4632-b7a6-458ee8eb9a53
snap_i = out[1] 

# ╔═╡ ece3e463-bc42-40c6-bcfd-be39222f433a
snap_f = out[end]

# ╔═╡ f19eeb21-3067-48c6-ba8d-5c7757bfe2ae
moment_inertia_tensor(snap_i)

# ╔═╡ 549fb8ce-3e4a-4fa8-b6ce-bde7a06a135b
eig = principal_axes(snap_f)

# ╔═╡ 46f5aa3b-8d8b-452f-8ed9-7ab6b60e1513
x_i = project_principal(snap_i)

# ╔═╡ 2d0de394-c79a-4870-b7a6-9228a47fb769
x_f = project_principal(snap_f)

# ╔═╡ b801dbdc-9e7e-49c4-ad04-4462aa22a30d
LilGuys.plot_xyz(x_i, plot=:scatter, markersize=1, alpha=0.01)

# ╔═╡ 93019bc6-cf42-46de-90c0-33ad4ce3e4de
LilGuys.plot_xyz(x_f[:, LilGuys.bound_particles(snap_f)], plot=:scatter, markersize=1, alpha=0.01, limits=((-10, 10), (-10, 10), (-10, 10)))

# ╔═╡ 0ca0c4df-6f59-4191-b14c-461d5ad906c8
import DensityEstimators: histogram2d

# ╔═╡ ff983edb-7461-487b-8ac9-0ad197f3d661
import KernelDensity: kde

# ╔═╡ a177f182-07d9-4520-89b9-a7c37ba35f24
import DensityEstimators: density2d

# ╔═╡ 332f418d-6b41-4b1a-9396-1a994660e5d2
function plot_isocontours(gs, positions; bandwidth=0.4, filter_bound=false, )
	
	r_max = 10
	x = positions[3, :]
	y = positions[1, :]

	h = kde((x, y), boundary=((-r_max, r_max), (-r_max, r_max)), bandwidth=(bandwidth, bandwidth))
	h_min = 1e-10

	h_max = log10(maximum(h.density) )
	
	ax = Axis(gs, xlabel=L"$x'$ / kpc", ylabel=L"$z'$ / kpc", )
	contour!(h.x, h.y, h.density, colorscale=log10, levels=10 .^ LinRange(h_max-2, h_max, 20), linewidth=theme(:linewidth)[]/2, linestyle=:solid)


	ax
end

# ╔═╡ 8e60f605-a2fb-4223-9cb9-7f11d8e2f24c
let
	fig = Figure()

	ax = plot_isocontours(fig[1,1], x_i)
	ax.title = "initial"

	ax = plot_isocontours(fig[1, 2], x_f)
	ax.title = "final"
	hideydecorations!(ticks=false, minorticks=false)

	rowsize!(fig.layout, 1, Aspect(1, 1))

	resize_to_layout!()
	@savefig "oblate_projected_2d"
	fig

end

# ╔═╡ 7c81c3b6-e51a-40f2-9b4a-006366e7bc34
function plot_isocontours_2(gs, positions; filter_bound=false, )
	
	r_max = 10
	x = positions[3, :]
	y = positions[1, :]

	# h = kde((x, y), boundary=((-r_max, r_max), (-r_max, r_max)))
	h = histogram2d(x, y, limits=((-r_max, r_max), (-r_max, r_max)))
	h_min = 1e-10

	h_max = log10(maximum(h.values) )
	
	ax = Axis(gs, xlabel=L"$x'$ / kpc", ylabel=L"$z'$ / kpc", )
	contour!(midpoints(h.xbins), midpoints(h.ybins), h.values, colorscale=log10, levels=10 .^ LinRange(h_max-2, h_max, 20), linewidth=theme(:linewidth)[]/2, linestyle=:solid)


	h
end

# ╔═╡ 102f5701-dafa-4c32-a463-644f41228e24
let
	fig = Figure()

	plot_isocontours_2(fig[1,1], x_i)

	plot_isocontours_2(fig[1, 2], x_f)
	hideydecorations!(ticks=false, minorticks=false)

	rowsize!(fig.layout, 1, Aspect(1, 1))
	fig

end

# ╔═╡ 1610d536-b75b-454a-83b5-dbbe09f0c4ad
import StatsBase: median

# ╔═╡ c3f8fd10-eae6-4f91-a078-008b214beb0e
[median(abs.(row)) for row in eachrow(x_f[:, LilGuys.bound_particles(snap_f)])]

# ╔═╡ 67b8052a-e0d0-41a0-8f49-0f20623365f4
[median(abs.(row)) for row in eachrow(x_i[:, LilGuys.bound_particles(snap_i)])]

# ╔═╡ bf26847d-e61d-488d-bbb2-e9a5023959f9
part_mass = snap_i.masses[1]

# ╔═╡ 5853e3b1-2534-4691-90ed-af99b3064355
function density_1d(positions, axis=1; weights=snap_i.masses, θ=π/12, num_per_bin=100)
	r = positions[axis, :]
	ax_2 = (axis + 1) % 3 + 1
	ax_3 = (axis) % 3 + 1
	@info "axes: $axis, $ax_2, $ax_3"

	r_tan = radii(positions[[ax_2, ax_3], :])
	filt = abs.(atan.(r_tan ./ r)) .< θ

	r_filt = abs.(r[filt])
	@info LilGuys.mean(filt)

	bins = LilGuys.bins_equal_number(r_filt, nothing, num_per_bin=num_per_bin)


	Nbins = length(bins) - 1

	ρ = zeros(Nbins)
	for i in 1:Nbins
		M = sum(weights[filt][bins[i] .< r_filt .< bins[i+1]])
		Ω = 4π * sin(θ/2)^2 * 2
		volume = Ω/3 * (bins[i+1]^3 - bins[i]^3)

		ρ[i] = M ./ volume 
	end

	return midpoints(bins), ρ
end

# ╔═╡ 33b20c20-69c1-491a-8b82-9e48048650c4
filt_bound = LilGuys.bound_particles(snap_f)

# ╔═╡ 56f556f6-b865-47b8-ace8-780e947828c9
halo = LilGuys.TruncNFW(r_circ_max=7, v_circ_max=48/V2KMS, trunc=20, xi=3)

# ╔═╡ a17fa4df-3bde-4cbc-8ff1-c29d47d84361
prof = LilGuys.DensityProfile(snap_i)

# ╔═╡ 2bbf3e1e-bbee-4b41-b06c-7863251ad334
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log radius / kpc",
		ylabel = "log density",
		limits=(-1, 2, -8, -1.5)
			 )

	x = LinRange(-1, 2, 1000)
	y = LilGuys.density.(halo, 10 .^ x)

	lw = theme(:linewidth)[]/2
	lines!(x, log10.(y), color=:black,  linewidth=lw)
	lines!(x .- log10.(2), log10.(y), color=:black, linewidth=lw)
	
	r, ρ = density_1d(x_i, 3)
	lines!(log10.(r), log10.(ρ), linestyle=:dot, color=COLORS[1], label="major axis, initial")


	r, ρ = density_1d(x_i)
	lines!(log10.(r), log10.(ρ), linestyle=:dot, color=COLORS[2], label="minor axis, initial")


	r, ρ = density_1d(x_f[:, filt_bound], weights=snap_f.masses[filt_bound], 3, num_per_bin=100)
	lines!(log10.(r), log10.(ρ), color=COLORS[1], label="major axis, final")
	r, ρ = density_1d(x_f[:, filt_bound], num_per_bin=50, weights=snap_f.masses[filt_bound])
	lines!(log10.(r), log10.(ρ), color=COLORS[2], label="minor axis, final")

	axislegend(position=:lb)

	@savefig "oblate_density_i_f"
	
	fig
end

# ╔═╡ 1caf2f60-224f-449f-9f85-8a3b07c191f6
halo_stars = LilGuys.Exp2D(R_s=0.10, M=Mstar)

# ╔═╡ e1dd8759-db92-4298-b89c-d5259ceb1ed2
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log radius / kpc",
		ylabel = "log density",
		limits=(-1, 2, -12, -1.5)
			 )

	x = LinRange(-1, 2, 1000)
	y = LilGuys.density.(halo_stars, 10 .^ x)

	lw = theme(:linewidth)[]/2
	lines!(x, log10.(y), color=:black,  linewidth=lw)
	lines!(x .- log10(2), log10.(y), color=:black, linewidth=lw)
	
	r, ρ = density_1d(x_i, 3, weights=snap_i.weights,θ=π/12)
	lines!(log10.(r), log10.(ρ), linestyle=:dot, color=COLORS[1], label="major axis, initial")


	r, ρ = density_1d(x_i, 1, weights=snap_i.weights, θ=π/6)
	lines!(log10.(r), log10.(ρ), linestyle=:dot, color=COLORS[2], label="minor axis, initial")


	r, ρ = density_1d(x_f[:, filt_bound], weights=snap_f.weights[filt_bound], 3, num_per_bin=100)
	lines!(log10.(r), log10.(ρ), color=COLORS[1], label="major axis, final")
	r, ρ = density_1d(x_f[:, filt_bound], num_per_bin=50, weights=snap_f.weights[filt_bound])
	lines!(log10.(r), log10.(ρ), color=COLORS[2], label="minor axis, final")

	axislegend(position=:lb)

	
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═ca9e7089-6617-4bff-996b-b16de0587c5b
# ╠═4cfdc774-f065-4d71-bdd9-ffd26dbe64a7
# ╠═1b202d9d-92ab-4798-8851-f3d8998b5c5f
# ╠═e9dcf228-7c2d-4deb-bbe5-b9215512bc70
# ╠═3d2f16ce-3d97-429a-9acb-004dc5162575
# ╠═c554e3d1-7ccc-4075-86c1-cf1df69afaf6
# ╠═6f24e979-b2ce-4590-b8d1-d442924ab2fd
# ╠═9a281d8f-748d-4632-b7a6-458ee8eb9a53
# ╠═ece3e463-bc42-40c6-bcfd-be39222f433a
# ╠═f19eeb21-3067-48c6-ba8d-5c7757bfe2ae
# ╠═549fb8ce-3e4a-4fa8-b6ce-bde7a06a135b
# ╠═46f5aa3b-8d8b-452f-8ed9-7ab6b60e1513
# ╠═2d0de394-c79a-4870-b7a6-9228a47fb769
# ╠═b801dbdc-9e7e-49c4-ad04-4462aa22a30d
# ╠═93019bc6-cf42-46de-90c0-33ad4ce3e4de
# ╠═0ca0c4df-6f59-4191-b14c-461d5ad906c8
# ╠═ff983edb-7461-487b-8ac9-0ad197f3d661
# ╠═8e60f605-a2fb-4223-9cb9-7f11d8e2f24c
# ╠═102f5701-dafa-4c32-a463-644f41228e24
# ╠═a177f182-07d9-4520-89b9-a7c37ba35f24
# ╠═332f418d-6b41-4b1a-9396-1a994660e5d2
# ╠═7c81c3b6-e51a-40f2-9b4a-006366e7bc34
# ╠═1610d536-b75b-454a-83b5-dbbe09f0c4ad
# ╠═c3f8fd10-eae6-4f91-a078-008b214beb0e
# ╠═67b8052a-e0d0-41a0-8f49-0f20623365f4
# ╠═bf26847d-e61d-488d-bbb2-e9a5023959f9
# ╠═5853e3b1-2534-4691-90ed-af99b3064355
# ╠═33b20c20-69c1-491a-8b82-9e48048650c4
# ╠═56f556f6-b865-47b8-ace8-780e947828c9
# ╠═a17fa4df-3bde-4cbc-8ff1-c29d47d84361
# ╠═2bbf3e1e-bbee-4b41-b06c-7863251ad334
# ╠═1caf2f60-224f-449f-9f85-8a3b07c191f6
# ╠═e1dd8759-db92-4298-b89c-d5259ceb1ed2
