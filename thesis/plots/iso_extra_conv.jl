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
r_label = "radius / kpc"

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
        xlabel=r_label,
		xscale=log10,
		xticks=[0.01, 0.1, 1, 10, 100],
        ylabel=L"$\textrm{v}_\textrm{circ}$ / km\,s$^{-1}$",
        limits=(xlims[1], xlims[2], ylims[1], ylims[2]),
		yscale=log10,
		yticks=[1:10; 10:10:100]
        )



    ax_res = Axis(fig[2, 1],
        xlabel=r_label,
		xscale=log10,
		xticks=[0.01, 0.1, 1, 10, 100],
        ylabel=L"\Delta\,\textrm{v}\,/\textrm{v}_\textrm{NFW}",
        limits=(xlims[1], xlims[2], -0.2, 0.2),
    )



	linkxaxes!(ax, ax_res, )
	rowsize!(fig.layout, 2, Auto(RES_SIZE))
	hidexdecorations!(ax, grid=false, ticks=false, minorticks=false)
	rowgap!(fig.layout, 0)

	return fig, ax, ax_res
end

# ╔═╡ c625ea14-662e-4dc5-851a-10dfba4df25f
function vcirc_ax(gs, ; xlims, ylims)    
    ax = Axis(gs,
        xlabel=r_label,
		xscale=log10,
		xticks=[0.01, 0.1, 1, 10, 100],
        ylabel=L"circular velocity / km\,s$^{-1}$",
        limits=(xlims[1], xlims[2], ylims[1], ylims[2]),
		yscale=log10,
		yticks=[1:10; 10:10:100]
        )

	return ax
end

# ╔═╡ 95af5b9d-e613-4d7f-849e-f6ddb221bd29
function plot_vcirc_profiles!(profiles, v_circ_0=x->0; errskip=1)
    for i in eachindex(profiles)
        label, profs = profiles[i]
        profile = profs[1]

		y0 = v_circ_0(profile.radii)
		y = (LilGuys.circular_velocity(profile) * V2KMS)
        lines!(profile.radii, y .- y0,
            linestyle=:dot,
            color=COLORS[i]
        )
        
        profile = profs[end]
		v0 = v_circ_0(profile.radii)

		y0 = v_circ_0(profile.radii)
		y = (LilGuys.circular_velocity(profile) * V2KMS)
		x = profile.radii
        lines!(x, y .- y0 ,
            color=COLORS[i],
            label=label
        )

        idx = 1:errskip:length(x)
        errorbars!(x[idx], middle.(y .- y0)[idx], LilGuys.error_interval.(y .- y0)[idx], color=COLORS[i])

	end
end

# ╔═╡ e019b3e3-a6d5-4430-a655-7bd785440a27
function get_lws(profiles)
	LinRange(3, 1, length(profiles))
end

# ╔═╡ 531169b0-55a7-428d-bf52-1a6a371336c9
function plot_vcirc_profiles_end!(profiles, v_circ_0=x->0; errskip=1)
	lws = get_lws(profiles)
    for i in eachindex(profiles)
        label, profs = profiles[i]

		profile = profs[end]
		y0 = v_circ_0(profile.radii)
		y = (LilGuys.circular_velocity(profile) * V2KMS)
		x = profile.radii
        lines!(x, y .- y0 ,
            color=COLORS[i],
            label=label,
			linewidth=lws[i]
        )

        idx = 1:errskip:length(x)
        errorbars!(x[idx], middle.(y .- y0)[idx], LilGuys.error_interval.(y .- y0)[idx], color=COLORS[i], linewidth=lws[i]
)

	end
end

# ╔═╡ 533eb462-eb9a-4452-b861-253e78b87490
function plot_marks!(profiles, marks::Real, ylims; residual=false, y0=nothing, kwargs...)
	x = marks
	
	if isnothing(y0)
		y = (LilGuys.v_circ(nfw, marks) * V2KMS)
	elseif y0 == :low
		y = ylims[1]
	else
		y = y0
	end
	
	dy = 0.05 * (ylims[2] - ylims[1])
	@info "marking ($x, $y)"

	if residual
		dy /= RES_SIZE
		arrows2d!([x], [y] .+ 2*dy, [0], [-dy]; kwargs...)
	else
		dy = (ylims[2] / ylims[1])^0.05
		arrows2d!([x], [y * dy^2], [0], [y*dy-y*dy^2]; kwargs...)

	end

end

# ╔═╡ ffa1eff7-57d9-4df4-868c-7cb48d5f1909
function plot_marks!(profiles, marks, ylims; kwargs...)
    for i in eachindex(profiles)
		label, _ = profiles[i]
		lw = get_lws(profiles[i])

		if label ∈ keys(marks)
			@info "marking $label"
			
			plot_marks!(profiles, marks[label], ylims; color=COLORS[i], shaftwidth=lw, kwargs...)
		end
    end
end

# ╔═╡ 27ee63c5-7923-4bbc-b6e6-690968df190d
function plot_analytic!(nfw, xlims)
    x = 10 .^ LinRange(xlims[1], xlims[2], 1000)
    y = LilGuys.v_circ.(nfw, x) 
    lines!(x, (y * V2KMS), linestyle=:dash, color=:black, label="NFW")
end

# ╔═╡ 9f1116ee-1fe2-4088-a242-f65a0f22f5f5
function get_r_res(profile)
	x = profile.radii
	
	y_exp = LilGuys.v_circ.(nfw, x)
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

# ╔═╡ 96d0bc0a-42e8-4127-9c5a-55a66925b14f
function compare_vcirc_orbit(profiles, nfw; 
					   errskip=1, xlims = (0.01, 50), vlims=(3, 36), marks=Dict())
	fig = Figure()
		
	ax = vcirc_ax(fig[1,1], xlims=xlims, ylims=vlims)
	
	plot_vcirc_profiles_end!(profiles)
	plot_marks!(profiles, marks, vlims)
		

    axislegend(position=:rb, patchsize=[24, 6])
	li = lines!([NaN], [NaN], linestyle=:dot, label="initial", color=:grey)
	lf = lines!([NaN], [NaN], linestyle=:solid, label="final", color=:grey)
	axislegend(ax, [li, lf], ["initial", "final"], position=:lt, patchsize=[24, 6])

    
    fig
end

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

# ╔═╡ 265f7694-c4e9-4e19-8c48-7d9afc40c1ad
profiles[1].second

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
	

    axislegend(position=:rb, patchsize=[24, 6])
	li = lines!([NaN], [NaN], linestyle=:dot, label="initial", color=:grey)
	lf = lines!([NaN], [NaN], linestyle=:solid, label="final", color=:grey)
	axislegend(ax, [li, lf], ["initial", "final"], position=:lt, patchsize=[24, 6])

    # residual
	Makie.current_axis!(ax_res)
	plot_residual_profiles!(profiles)

	plot_marks!(profiles, marks, (-0.2, 0.2), residual=true)
	plot_marks!(profiles, softening, (-0.2, 0.2), residual=true,  y0=:low, tiplength=0, shaftwidth=smalllinewidth)

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

# ╔═╡ 46d7290c-a546-4910-90b3-a6184760a29d
# @savefig "iso_converg_methods" compare_vcirc(profiles_methods, errskip=1, marks=r_conv(halo, 1e5), xlims=(-1, 1.6), vlims=(7.5, 36))

# ╔═╡ b0131440-cea1-4791-a5fa-883f5746269a
profiles_soft_orbit = [
    "benchmark" => load_profile("../sculptor/1e7_new_v31_r3.2/orbit_smallperi"),
	"larger softening" => load_profile("../sculptor/1e5_s0.44_v31_r3.2/orbit_smallperi"),
    "fiducial" => load_profile("../sculptor/1e5_new_v31_r3.2/orbit_smallperi"),
    "smaller softening" => load_profile("../sculptor/1e5_s0.044_v31_r3.2/orbit_smallperi"),
];

# ╔═╡ 6c8633e7-3460-400e-b7bd-a544672904dd


# ╔═╡ bfa9eeab-0e2a-411f-85a1-8c0008084d41
profiles_methods_orbit = [
    "benchmark" => load_profile("../sculptor/1e7_new_v31_r3.2/orbit_smallperi"),
    "fiducial" => load_profile("../sculptor/1e5_new_v31_r3.2/orbit_smallperi"),
	"small timestep" => load_profile("../sculptor/1e5_acc0.003_v31_r3.2//orbit_smallperi"),
    "high acc. force" => load_profile("../sculptor/1e5_f0.001_v31_r3.2//orbit_smallperi"),
];

# ╔═╡ 07d70f01-b3b6-4183-8d74-5b7ae56a2496
profiles_orbit = [
    "benchmark" => load_profile("../sculptor/1e7_new_v31_r3.2/orbit_smallperi"),
    "1e6" => load_profile("../sculptor/1e6_new_v31_r3.2/orbit_smallperi"),
    "1e5" => load_profile("../sculptor/1e5_new_v31_r3.2/orbit_smallperi"),
];

# ╔═╡ 5dc1fe46-8c1d-4c83-b1db-1f2d82a1f948
halo_smallperi = LilGuys.TruncNFW(r_circ_max=3.2, v_circ_max=31/V2KMS, trunc=20, xi=3)

# ╔═╡ 487b087b-1cd5-40c7-a481-3594dc0ffca6
compare_vcirc_orbit(profiles_orbit, halo_smallperi)

# ╔═╡ 10e4958b-95ae-4245-a046-be18273b0f41
0.14*halo_smallperi.r_s/ halo.r_s

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

# ╔═╡ c25d0e36-7d7c-438e-ad45-61abd3278a75
@savefig "iso_converg_softening" compare_vcirc(profiles_softening, marks = r_conv(halo, 1e6), xlims=(0.1, 50), vlims=(5.5, 38))

# ╔═╡ 27a3628b-35dd-450f-8267-16349b341a64
@savefig "iso_converg_methods" compare_vcirc(profiles_methods, marks = r_conv(halo, 1e5), xlims=(0.05, 50), vlims=(5.5, 38))

# ╔═╡ 0bf50703-33fc-493f-8c93-e7a3e6d8d35e
@savefig "orbit_converg_methods" let
	fig = Figure()
	h0 = 0.14*halo_smallperi.r_s/ halo.r_s
	h_soft_orbit = Dict(
	    "benchmark" => h0/10,
		"larger softening" => h0*sqrt(10),
	    "fiducial" => h0,
	    "smaller softening" => h0/sqrt(10),
	)
	
	xlims=(0.022, 10)
	vlims=(6, 25)
	ax = vcirc_ax(fig[1,1], xlims=xlims, ylims=vlims)
	plot_vcirc_profiles_end!(profiles_soft_orbit)
	axislegend(position=:rb, patchsize=(30, 5))
	
	plot_marks!(profiles_soft_orbit, r_conv(halo_smallperi, 1e5), vlims)
	plot_marks!(profiles_soft_orbit, h_soft_orbit, vlims, tiplength=0, shaftwidth=smalllinewidth, y0=5.6)
	text!(0.07505790373441329, 5.6, text="softening", fontsize=0.8*theme(:fontsize)[], offset=(0, 20), align=(:left, :center), rotation=π/2, color=COLORS[3])


	
	ax2 = vcirc_ax(fig[1, 2], xlims=xlims, ylims=vlims)
	plot_vcirc_profiles_end!(profiles_methods_orbit)
	plot_marks!(profiles_methods_orbit, r_conv(halo_smallperi, 1e5), vlims)
	text!(0.6132391588137316, 19.580283390724112, text="converged", fontsize=0.8*theme(:fontsize)[], offset=(0, 20), align=(:center, :bottom))

	plot_marks!(profiles_methods_orbit, h0, vlims, tiplength=0, shaftwidth=smalllinewidth, y0=5.6)


	hideydecorations!(ticks=false, minorticks=false)


	axislegend(position=:rb, patchsize=(30, 5))

	rowsize!(fig.layout, 1, Aspect(1, 1))

	resize_to_layout!()
	fig
end

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
# ╠═c625ea14-662e-4dc5-851a-10dfba4df25f
# ╠═6456950b-d95a-4fc9-942a-0d90e140ca39
# ╠═95af5b9d-e613-4d7f-849e-f6ddb221bd29
# ╠═e019b3e3-a6d5-4430-a655-7bd785440a27
# ╠═531169b0-55a7-428d-bf52-1a6a371336c9
# ╠═533eb462-eb9a-4452-b861-253e78b87490
# ╠═ffa1eff7-57d9-4df4-868c-7cb48d5f1909
# ╠═27ee63c5-7923-4bbc-b6e6-690968df190d
# ╠═9f1116ee-1fe2-4088-a242-f65a0f22f5f5
# ╠═64efb625-9fed-4ec9-9546-79842e7e6b90
# ╠═cc50d80b-78b8-4f39-8a41-18ae79bf7287
# ╠═96d0bc0a-42e8-4127-9c5a-55a66925b14f
# ╠═685c9c53-4c4e-4d1f-82da-a64488562662
# ╠═f7cee4fe-e8df-429b-9cb0-f508d60c9788
# ╠═a4b6b203-a42a-48dd-8944-ea150eee376b
# ╠═265f7694-c4e9-4e19-8c48-7d9afc40c1ad
# ╠═9c3c8684-de16-4601-a6e0-7bcb9c7a1f84
# ╠═3aacc4eb-5e7b-4433-adac-dd83d41b62e8
# ╠═69516e98-4e15-44a8-b0d4-dc8d8bacd7e8
# ╠═1453c857-dfdd-4afc-98ec-2ec794c129d4
# ╠═a1d729da-b041-467e-8ec7-1d88cda2d5a2
# ╠═0cba9678-4469-4a14-93cb-417c01beb18d
# ╠═49c4798b-59e3-4908-90f8-e0a6801a64fd
# ╠═c25d0e36-7d7c-438e-ad45-61abd3278a75
# ╠═1e109d8c-d7e4-4f46-9e39-3adfc881f563
# ╠═27a3628b-35dd-450f-8267-16349b341a64
# ╠═886b0f81-dca7-47e2-847b-bcd3450d4d07
# ╠═46d7290c-a546-4910-90b3-a6184760a29d
# ╠═b0131440-cea1-4791-a5fa-883f5746269a
# ╠═6c8633e7-3460-400e-b7bd-a544672904dd
# ╠═bfa9eeab-0e2a-411f-85a1-8c0008084d41
# ╠═07d70f01-b3b6-4183-8d74-5b7ae56a2496
# ╠═487b087b-1cd5-40c7-a481-3594dc0ffca6
# ╠═5dc1fe46-8c1d-4c83-b1db-1f2d82a1f948
# ╠═10e4958b-95ae-4245-a046-be18273b0f41
# ╠═0bf50703-33fc-493f-8c93-e7a3e6d8d35e
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
