### A Pluto.jl notebook ###
# v0.20.8

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
log_r_label = "log r / kpc"

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
        ylabel=L"$v_\textrm{circ}$ / km\,s$^{-1}$",
        limits=(xlims[1], xlims[2], ylims[1], ylims[2]),
		yscale=log10,
		yticks=[1:10; 10:10:100]
        )



    ax_res = Axis(fig[2, 1],
        xlabel=log_r_label,
        ylabel=L"\Delta\,v\,/v_\textrm{exp}",
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
function plot_marks!(profiles, marks::Real, ylims; residual=false, kwargs...)
	x = log10(marks)
	
	if !residual
		y = (LilGuys.v_circ(nfw, marks) * V2KMS)
	else
		y = 0
	end

	dy = 0.05 * (ylims[2] - ylims[1])
	if residual
		dy /= RES_SIZE
		@info "marking ($x, $y)"
		arrows!([x], [y] .+ 2*dy, [0], [-dy]; kwargs...)
	else
		dy = (ylims[2] / ylims[1])^0.05
		arrows!([x], [y * dy^2], [0], [y*dy-y*dy^2]; kwargs...)

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

# ╔═╡ cc50d80b-78b8-4f39-8a41-18ae79bf7287
function compare_vcirc(profiles; 
					   errskip=1, xlims = (-2, 1.6), vlims=(3, 36), marks=Dict())
	fig, ax, ax_res = vcirc_axes(xlims=xlims, ylims=vlims)
	
	Makie.current_axis!(ax)

	plot_vcirc_profiles!(profiles)
	plot_marks!(profiles, marks, vlims)
	plot_analytic!(nfw, xlims)
	

    axislegend(position=:rb, patchsize=[24, 6])
	li = lines!([NaN], [NaN], linestyle=:dot, label="initial", color=:grey)
	lf = lines!([NaN], [NaN], linestyle=:solid, label="final", color=:grey)
	axislegend(ax, [li, lf], ["initial", "final"], position=:lt, patchsize=[24, 6])

    # residual
	Makie.current_axis!(ax_res)
	plot_residual_profiles!(profiles)

	plot_marks!(profiles, marks, (-0.2, 0.2), residual=true)

    hlines!(0, color=:black, linestyle=:dash)

    
    fig
end

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

# ╔═╡ 265f7694-c4e9-4e19-8c48-7d9afc40c1ad
profiles[1].second

# ╔═╡ 9c3c8684-de16-4601-a6e0-7bcb9c7a1f84
softening = Dict(
	"10⁴" => 0.44,
	"10⁵" => 0.14,
	"10⁶" => 0.044,
	"10⁷" => 0.014,
)

# ╔═╡ 3aacc4eb-5e7b-4433-adac-dd83d41b62e8
part_num = Dict(
	"10⁴" => 1e4,
	"10⁵" => 1e5,
	"10⁶" => 1e6,
	"10⁷" => 1e7,

)

# ╔═╡ 69516e98-4e15-44a8-b0d4-dc8d8bacd7e8
CairoMakie.activate!(type=:png)

# ╔═╡ 1453c857-dfdd-4afc-98ec-2ec794c129d4
10^1.4

# ╔═╡ a1d729da-b041-467e-8ec7-1d88cda2d5a2
log10(LilGuys.R200(halo))

# ╔═╡ 0cba9678-4469-4a14-93cb-417c01beb18d
profiles_softening = [

	L"$h_\textrm{grav} = 0.14\,$kpc" => load_profile("1e6_new/s0.14"),
    L"$h_\textrm{grav} = 0.044\,$kpc" => load_profile("1e6_new/fiducial"),
    L"$h_\textrm{grav} = 0.014\,$kpc" => load_profile("1e6_new/s0.014"),
    ];

# ╔═╡ 49c4798b-59e3-4908-90f8-e0a6801a64fd
[prof[end].time for (k, prof) in profiles_softening]

# ╔═╡ 1e109d8c-d7e4-4f46-9e39-3adfc881f563
profiles_methods = [
    "fiducial" => load_profile("1e5_new/fiducial"),
 	"gadget2" => load_profile("1e5_new/gadget2"),
    "dt0.1" => load_profile("1e5_new/dt_0.1"),
    "high accuracy" => load_profile("1e5_new/f_0.001"),
   "geometric" => load_profile("1e5_new/geometric"),
];

# ╔═╡ 886b0f81-dca7-47e2-847b-bcd3450d4d07
[prof[end].time for (k, prof) in profiles_methods]

# ╔═╡ 42696671-2828-4d38-a3dd-2a140a271d58
1 / sqrt(1330) 

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

# ╔═╡ 755823d8-4e56-4c30-9199-dcd79df79e86
t_circ(0.01)

# ╔═╡ 05eb20ad-ce15-4bf8-b819-7a766f3ebddb
N(0.01, halo)

# ╔═╡ 85a345df-3647-4003-ba4f-15c775a94658
t_relax(10ϵ, halo)

# ╔═╡ 1ce1a5e2-7269-41f8-a633-d29b7f07855f
r_conv(halo, N=1e7) = LilGuys.Interface.find_zero(r -> t_relax.(r, halo, N) .-  1176, [ϵ, 100])

# ╔═╡ 9ea7fe71-67a0-41ab-bb5b-3d191eb4275c
@savefig "iso_converg_num" compare_vcirc(profiles, marks=Dict(k => r_conv(halo, v) for (k, v) in part_num), vlims=(3.5, 38))

# ╔═╡ c25d0e36-7d7c-438e-ad45-61abd3278a75
@savefig "iso_converg_softening" compare_vcirc(profiles_softening, marks = r_conv(halo, 1e6), xlims=(-1.2, 1.6), vlims=(5.5, 38))

# ╔═╡ 46d7290c-a546-4910-90b3-a6184760a29d
@savefig "iso_converg_methods" compare_vcirc(profiles_methods, errskip=1, marks=r_conv(halo, 1e5), xlims=(-1, 1.6), vlims=(7.5, 36))

# ╔═╡ bc661866-5f31-4193-8ea7-d3b0a145fcdd
r_conv_0 = r_conv(halo)

# ╔═╡ fa43c524-40c4-490c-a92b-1a88bbd7827f
log10(ρ_rel(r_conv_0, halo))

# ╔═╡ 62563803-0853-4c9c-b7a1-904af9392339
t_relax(r_conv_0, halo)

# ╔═╡ 927bc186-67af-483c-bc8b-5cb855706f66
t_circ(r_conv_0) * N(r_conv_0, halo) / (8log(N(r_conv_0, halo)))

# ╔═╡ 6df842d4-436d-4473-b8ce-c3a82a7346cf
r_conv_0 / ϵ

# ╔═╡ 7324abb7-9496-44d7-932e-a6d273477a4c
0.044 * 10

# ╔═╡ 8207080f-5b93-4d23-ac98-3580cefcb9fe
N(r_conv_0, halo)

# ╔═╡ 8ed9880b-25d7-437a-9144-b61dc399e535
v_circ(halo, 0.01)

# ╔═╡ 05e1b82c-c931-4e62-af27-2d9501e59877
t_circ(0.01)

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
# ╠═cc50d80b-78b8-4f39-8a41-18ae79bf7287
# ╠═f7cee4fe-e8df-429b-9cb0-f508d60c9788
# ╠═a4b6b203-a42a-48dd-8944-ea150eee376b
# ╠═265f7694-c4e9-4e19-8c48-7d9afc40c1ad
# ╠═9c3c8684-de16-4601-a6e0-7bcb9c7a1f84
# ╠═3aacc4eb-5e7b-4433-adac-dd83d41b62e8
# ╠═69516e98-4e15-44a8-b0d4-dc8d8bacd7e8
# ╠═1453c857-dfdd-4afc-98ec-2ec794c129d4
# ╠═9ea7fe71-67a0-41ab-bb5b-3d191eb4275c
# ╠═a1d729da-b041-467e-8ec7-1d88cda2d5a2
# ╠═0cba9678-4469-4a14-93cb-417c01beb18d
# ╠═49c4798b-59e3-4908-90f8-e0a6801a64fd
# ╠═c25d0e36-7d7c-438e-ad45-61abd3278a75
# ╠═1e109d8c-d7e4-4f46-9e39-3adfc881f563
# ╠═886b0f81-dca7-47e2-847b-bcd3450d4d07
# ╠═46d7290c-a546-4910-90b3-a6184760a29d
# ╠═42696671-2828-4d38-a3dd-2a140a271d58
# ╟─3c0bb9b7-d1b5-4a67-a9ad-8a89c615110f
# ╠═7f0a1887-bf40-4fe6-930e-cee87bf4d32a
# ╠═7341e90b-355d-449b-a913-c2bae8ec9249
# ╠═6d118749-a1ec-4ca3-bd7d-2bc36ddf864c
# ╠═89259f41-1b8a-42ee-842b-0cde1749ed6d
# ╠═cb0478f4-8399-42eb-b3b9-538c1de69385
# ╠═fa43c524-40c4-490c-a92b-1a88bbd7827f
# ╠═a4f8a47a-e708-4674-b5dd-d3616e0e5d1d
# ╠═62563803-0853-4c9c-b7a1-904af9392339
# ╠═927bc186-67af-483c-bc8b-5cb855706f66
# ╠═bc661866-5f31-4193-8ea7-d3b0a145fcdd
# ╠═755823d8-4e56-4c30-9199-dcd79df79e86
# ╠═05eb20ad-ce15-4bf8-b819-7a766f3ebddb
# ╠═85a345df-3647-4003-ba4f-15c775a94658
# ╠═1ce1a5e2-7269-41f8-a633-d29b7f07855f
# ╠═6df842d4-436d-4473-b8ce-c3a82a7346cf
# ╠═7324abb7-9496-44d7-932e-a6d273477a4c
# ╠═8207080f-5b93-4d23-ac98-3580cefcb9fe
# ╠═8ed9880b-25d7-437a-9144-b61dc399e535
# ╠═05e1b82c-c931-4e62-af27-2d9501e59877
