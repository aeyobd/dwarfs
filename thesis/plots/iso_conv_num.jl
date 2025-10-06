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

# ╔═╡ 2abed1bf-9861-49b0-9bea-aa9e9f67829e
include("./paper_style.jl")

# ╔═╡ 7ae386b0-8c0b-43c1-b414-ac0518b73c70
models_dir = joinpath(ENV["DWARFS_ROOT"], "analysis/isolation")

# ╔═╡ 0f08434b-6f55-4d73-90eb-750e5d979a47
function load_profile(name) 
    path = joinpath(models_dir, "$name/profiles.hdf5")
    profiles = LilGuys.read_structs_from_hdf5(path, LilGuys.MassProfile)
    idx = parse.(Int, first.(profiles))
    profiles = last.(profiles)
    return profiles[sortperm(idx)]
    
end

# ╔═╡ 2a0d9629-0f16-4c51-bcca-93d9d97ea915
halo = LilGuys.load_profile("$models_dir/1e6_new/fiducial/halo.toml")

# ╔═╡ f095ddc2-4d80-49bb-8a55-ce0f2e6524b7
nfw = NFW(M_s=halo.M_s, r_s=halo.r_s, c=halo.c')

# ╔═╡ 2c77e958-c57c-4bb0-94f3-f29414e275cb


# ╔═╡ ba744355-8da6-4caa-a252-ad83ad4b8b13
log_r_label = "log radius / kpc"

# ╔═╡ bd023583-28b3-4c4d-8ced-7161350a5d7a
md"""
## Ploting utilities
"""

# ╔═╡ 57818793-f612-42ed-a003-9e4a4810a163
RES_SIZE = 0.3

# ╔═╡ c50be256-ab72-4d87-90a8-0fac98ef35b7
function vcirc_axes(; xlims, ylims)
    
    fig = Figure()
    ax = Axis(fig[1, 1],
        xlabel=log_r_label,
        ylabel=L"$\textrm{v}_\textrm{circ}$ / km\,s$^{-1}$",
        limits=(xlims[1], xlims[2], ylims[1], ylims[2]),
		yscale=log10,
		yticks=[1:10; 10:10:100]
        )



    ax_res = Axis(fig[2, 1],
        xlabel=log_r_label,
        ylabel=L"\Delta\,\textrm{v}\,/\textrm{v}_\textrm{NFW}",
        limits=(xlims[1], xlims[2], -0.2, 0.2),
    )



	linkxaxes!(ax, ax_res, )
	rowsize!(fig.layout, 2, Auto(RES_SIZE))
	hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
	rowgap!(fig.layout, 0)

	return fig, ax, ax_res
end

# ╔═╡ 95af5b9d-e613-4d7f-849e-f6ddb221bd29
function plot_vcirc_profiles!(profiles, v_circ_0=x->0; errskip=1)
    for i in eachindex(profiles)
        label, profs = profiles[i]
        profile = profs[1]

		y0 = v_circ_0(profile.radii)
		y = (LilGuys.circular_velocity(profile) * V2KMS)
        lines!(log10.(profile.radii), y .- y0,
            linestyle=:dot,
            color=COLORS[i]
        )
        
        profile = profs[end]
		v0 = v_circ_0(profile.radii)

		y0 = v_circ_0(profile.radii)
		y = (LilGuys.circular_velocity(profile) * V2KMS)
		x = log10.(profile.radii)
        lines!(x, y .- y0 ,
            color=COLORS[i],
            label=label
        )

        idx = 1:errskip:length(x)
        errorbars!(x[idx], middle.(y .- y0)[idx], LilGuys.error_interval.(y .- y0)[idx], color=COLORS[i])

	end
end

# ╔═╡ 533eb462-eb9a-4452-b861-253e78b87490
function plot_marks!(profiles, marks::Real, ylims; y0=nothing, residual=false, kwargs...)
	x = log10(marks)
	
	if isnothing(y0)
		y = (LilGuys.v_circ(nfw, marks) * V2KMS)
	elseif y0 == :low
		y = ylims[1]
	else
		y = y0
	end

	dy = 0.05 * (ylims[2] - ylims[1])
	if residual
		dy /= RES_SIZE
		@info "marking ($x, $y, $dy)"
		arrows2d!([x], [y] .+ 2*dy, [0], [-dy]; kwargs...)
	else
		dy = (ylims[2] / ylims[1])^0.05
		@info "marking ($x, $y, $dy)"

		arrows2d!([x], [y * dy^2], [0], [y*dy-y*dy^2]; kwargs...)

	end

end

# ╔═╡ ffa1eff7-57d9-4df4-868c-7cb48d5f1909
function plot_marks!(profiles, marks, ylims; kwargs...)
    for i in eachindex(profiles)
		label, _ = profiles[i]

		if label ∈ keys(marks)
			@info "marking $label"
			
			plot_marks!(profiles, marks[label], ylims; color=COLORS[i], kwargs...)
		end
    end
end

# ╔═╡ 27ee63c5-7923-4bbc-b6e6-690968df190d
function plot_analytic!(nfw, xlims)
    x = LinRange(xlims[1], xlims[2], 1000)
    y = LilGuys.v_circ.(nfw, 10 .^ x) 
    lines!(x, (y * V2KMS), linestyle=:dash, color=:black, label="NFW")
end

# ╔═╡ 9f1116ee-1fe2-4088-a242-f65a0f22f5f5
function get_r_res(profile)
	x = log10.(profile.radii)
	
	y_exp = LilGuys.v_circ.(nfw, 10 .^ x)
	dy = LilGuys.circular_velocity(profile) .- y_exp
	res = dy ./ y_exp
	res_err = LilGuys.sym_error.(res)
	res = LilGuys.middle.(res)
	return x, res, res_err
end

# ╔═╡ 64efb625-9fed-4ec9-9546-79842e7e6b90
function plot_residual_profiles!(profiles; errskip=1)  
    for i in eachindex(profiles)
        label, profs = profiles[i]

        profile = profs[1]
		x, res, res_err = get_r_res(profile)
        lines!(x, res, color=COLORS[i], linestyle=:dot)

        profile = profs[end]
		x, res, res_err = get_r_res(profile)
        scatterlines!(x, res, color=COLORS[i], markersize=3)

        idx = 1:errskip:length(res)
        errorbars!(x[idx], res[idx], res_err[idx], color=COLORS[i])

    end
end

# ╔═╡ 69516e98-4e15-44a8-b0d4-dc8d8bacd7e8
CairoMakie.activate!(type=:png)

# ╔═╡ 00dc1954-0298-4f20-8ad7-adf2d886c74c
smallfontsize=0.8*theme(:fontsize)[]

# ╔═╡ 685c9c53-4c4e-4d1f-82da-a64488562662
smalllinewidth = theme(:linewidth)[] /  2

# ╔═╡ f7cee4fe-e8df-429b-9cb0-f508d60c9788
load_profile("1e5/fiducial")

# ╔═╡ a4b6b203-a42a-48dd-8944-ea150eee376b
profiles = [
    "10⁴" => load_profile("1e4_new/fiducial"),
    "10⁵" => load_profile("1e5_new/fiducial"),
	"10⁶" => load_profile("1e6_new/fiducial"),
     "10⁷" => load_profile("1e7_new/fiducial")
    ];

# ╔═╡ 6456950b-d95a-4fc9-942a-0d90e140ca39
LilGuys.circular_velocity(profiles[1].second[1])

# ╔═╡ 9c3c8684-de16-4601-a6e0-7bcb9c7a1f84
softening = Dict(
	"10⁴" => 0.44,
	"10⁵" => 0.14,
	"10⁶" => 0.044,
	"10⁷" => 0.014,
)

# ╔═╡ cc50d80b-78b8-4f39-8a41-18ae79bf7287
function compare_vcirc(profiles; 
					   errskip=1, xlims = (-2, 1.6), vlims=(3, 36), marks=Dict())
	fig, ax, ax_res = vcirc_axes(xlims=xlims, ylims=vlims)
	
	Makie.current_axis!(ax)

	plot_vcirc_profiles!(profiles)
	plot_marks!(profiles, marks, vlims)
	plot_marks!(profiles, softening, vlims, tiplength=0, y0=:low, shaftwidth=smalllinewidth)
	plot_analytic!(nfw, xlims)
	text!(log10(0.014), 3.5*1.1266415859202983^2, text="softening", offset=(0, 6), align=(:left, :center), fontsize=smallfontsize, rotation=π/2, color=COLORS[4])
	text!(-0.9301245845065439, 9.53748778279495*1.1266415859202983^2, text="converged", offset=(0, 6), align=(:left, :center), fontsize=smallfontsize, rotation=π/2, color=COLORS[4])

    axislegend(position=:rb, patchsize=[24, 6])
	li = lines!([NaN], [NaN], linestyle=:dot, label="initial", color=:grey)
	lf = lines!([NaN], [NaN], linestyle=:solid, label="final", color=:grey)
	axislegend(ax, [li, lf], ["initial", "final"], position=:lt, patchsize=[24, 6])

    # residual
	Makie.current_axis!(ax_res)
	plot_residual_profiles!(profiles)

	plot_marks!(profiles, marks, (-0.2, 0.2), y0=0, residual=true)

    hlines!(0, color=:black, linestyle=:dash)

    
    fig
end

# ╔═╡ 3aacc4eb-5e7b-4433-adac-dd83d41b62e8
part_num = Dict(
	"10⁴" => 1e4,
	"10⁵" => 1e5,
	"10⁶" => 1e6,
	"10⁷" => 1e7,

)

# ╔═╡ 3c0bb9b7-d1b5-4a67-a9ad-8a89c615110f
md"""
# Convergence tests
"""

# ╔═╡ 7f0a1887-bf40-4fe6-930e-cee87bf4d32a
Ntot = 1e7

# ╔═╡ 7341e90b-355d-449b-a913-c2bae8ec9249
N(r, halo, Ntot=1e7) = Ntot * LilGuys.mass(halo, r) / LilGuys.mass(halo)

# ╔═╡ 6d118749-a1ec-4ca3-bd7d-2bc36ddf864c
ϵ = 0.014 * sqrt(Ntot / 1e7)^-1

# ╔═╡ 89259f41-1b8a-42ee-842b-0cde1749ed6d
function t_circ(r)
	return 2π*r / v_circ(halo, r)
end

# ╔═╡ cb0478f4-8399-42eb-b3b9-538c1de69385
ρ_rel(r, halo) = LilGuys.mean_density(halo, r) / LilGuys.ρ_crit

# ╔═╡ a4f8a47a-e708-4674-b5dd-d3616e0e5d1d
t_relax(r, halo, Ntot=1e7) = t_circ(LilGuys.R200(halo)) * sqrt(200) / 8 * N(r, halo, Ntot) / log(N(r, halo, Ntot)) * (ρ_rel(r, halo))^(-1/2)

# ╔═╡ 1ce1a5e2-7269-41f8-a633-d29b7f07855f
r_conv(halo, N=1e7) = LilGuys.Interface.find_zero(r -> t_relax.(r, halo, N) .-  1060, [ϵ, 100])

# ╔═╡ 9ea7fe71-67a0-41ab-bb5b-3d191eb4275c
@savefig "iso_converg_num" compare_vcirc(profiles, marks=Dict(k => r_conv(halo, v) for (k, v) in part_num), vlims=(3.5, 38))

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═2abed1bf-9861-49b0-9bea-aa9e9f67829e
# ╠═0f08434b-6f55-4d73-90eb-750e5d979a47
# ╠═7ae386b0-8c0b-43c1-b414-ac0518b73c70
# ╠═2a0d9629-0f16-4c51-bcca-93d9d97ea915
# ╠═f095ddc2-4d80-49bb-8a55-ce0f2e6524b7
# ╠═2c77e958-c57c-4bb0-94f3-f29414e275cb
# ╠═ba744355-8da6-4caa-a252-ad83ad4b8b13
# ╠═bd023583-28b3-4c4d-8ced-7161350a5d7a
# ╠═57818793-f612-42ed-a003-9e4a4810a163
# ╠═c50be256-ab72-4d87-90a8-0fac98ef35b7
# ╠═6456950b-d95a-4fc9-942a-0d90e140ca39
# ╠═95af5b9d-e613-4d7f-849e-f6ddb221bd29
# ╠═533eb462-eb9a-4452-b861-253e78b87490
# ╠═ffa1eff7-57d9-4df4-868c-7cb48d5f1909
# ╠═27ee63c5-7923-4bbc-b6e6-690968df190d
# ╠═9f1116ee-1fe2-4088-a242-f65a0f22f5f5
# ╠═64efb625-9fed-4ec9-9546-79842e7e6b90
# ╠═69516e98-4e15-44a8-b0d4-dc8d8bacd7e8
# ╠═cc50d80b-78b8-4f39-8a41-18ae79bf7287
# ╠═00dc1954-0298-4f20-8ad7-adf2d886c74c
# ╠═685c9c53-4c4e-4d1f-82da-a64488562662
# ╠═f7cee4fe-e8df-429b-9cb0-f508d60c9788
# ╠═a4b6b203-a42a-48dd-8944-ea150eee376b
# ╠═9c3c8684-de16-4601-a6e0-7bcb9c7a1f84
# ╠═3aacc4eb-5e7b-4433-adac-dd83d41b62e8
# ╠═9ea7fe71-67a0-41ab-bb5b-3d191eb4275c
# ╟─3c0bb9b7-d1b5-4a67-a9ad-8a89c615110f
# ╠═7f0a1887-bf40-4fe6-930e-cee87bf4d32a
# ╠═7341e90b-355d-449b-a913-c2bae8ec9249
# ╠═6d118749-a1ec-4ca3-bd7d-2bc36ddf864c
# ╠═89259f41-1b8a-42ee-842b-0cde1749ed6d
# ╠═cb0478f4-8399-42eb-b3b9-538c1de69385
# ╠═a4f8a47a-e708-4674-b5dd-d3616e0e5d1d
# ╠═1ce1a5e2-7269-41f8-a633-d29b7f07855f
