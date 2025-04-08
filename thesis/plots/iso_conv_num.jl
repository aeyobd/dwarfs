### A Pluto.jl notebook ###
# v0.20.5

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
    profiles = LilGuys.read_structs_from_hdf5(path, LilGuys.MassProfile3D)
    idx = parse.(Int, first.(profiles))
    profiles = last.(profiles)
    return profiles[sortperm(idx)]
    
end

# ╔═╡ 2c77e958-c57c-4bb0-94f3-f29414e275cb
halo = LilGuys.load_profile("$models_dir/1e6/fiducial/halo.toml")

# ╔═╡ ba744355-8da6-4caa-a252-ad83ad4b8b13
log_r_label = "log r / kpc"

# ╔═╡ cc50d80b-78b8-4f39-8a41-18ae79bf7287
function compare_vcirc(profiles; errskip=1, xlims = (-2, 3), vlims=(0.3, 1.6))
    
    fig = Figure()
    ax = Axis(fig[1, 1],
        xlabel=log_r_label,
        ylabel=L"$\log\,v_\textrm{circ}$ / km\,s$^{-1}$",
        limits=(xlims[1], xlims[2], vlims[1], vlims[2]),
        )
    pi = 1

    
    for i in eachindex(profiles)
        label, profs = profiles[i]
        profile = profs[1]
        lines!(log10.(profile.r_circ), log10.(profile.v_circ * V2KMS),
            linestyle=:dot,
            color=COLORS[i]
        )
        
        profile = profs[end]
        lines!(log10.(profile.r_circ), log10.(profile.v_circ* V2KMS), 
            color=COLORS[i],
            label=label
        )

        println(label, " number per bin: ", LilGuys.mean(diff(profile.n_circ)))

    
    end

    x = LinRange(xlims[1], xlims[2], 1000)
    y = LilGuys.v_circ.(halo, 10 .^ x)
    lines!(x, log10.(y * V2KMS), linestyle=:dash, color=:black, label="NFW")

    axislegend(position=:rb)

    # residual

    ax_res = Axis(fig[2, 1],
        xlabel=log_r_label,
        ylabel=L"\Delta\,v\,/v_\textrm{exp}",
        limits=(xlims[1], xlims[2], -0.2, 0.2),
    )
    
    for i in eachindex(profiles)
        label, profs = profiles[i]
        profile = profs[end]

        x = log10.(profile.r_circ)
        
        y_exp = LilGuys.v_circ.(halo, 10 .^ x)
        dy = profile.v_circ .- y_exp
        res = dy ./ y_exp
        res_err = profile.v_circ_err ./ y_exp
        
        scatterlines!(x, res, color=COLORS[i], markersize=3)

        idx = 1:errskip:length(res)
        errorbars!(x[idx], res[idx], res_err[idx], color=COLORS[i])
    end
    hlines!(0, color=:black)


	linkxaxes!(ax, ax_res, )
	rowsize!(fig.layout, 2, Auto(0.3))
	hidexdecorations!(ax, grid=false)
    
    fig
end

# ╔═╡ f7cee4fe-e8df-429b-9cb0-f508d60c9788
load_profile("1e5/fiducial")

# ╔═╡ a4b6b203-a42a-48dd-8944-ea150eee376b
profiles = [
    "1e4" => load_profile("1e4/fiducial"),
    "1e5" => load_profile("1e5/fiducial"),
    "1e6" => load_profile("1e6/fiducial"),
    "1e7" => load_profile("1e7/fiducial")
    ];

# ╔═╡ 69516e98-4e15-44a8-b0d4-dc8d8bacd7e8
CairoMakie.activate!(type=:png)

# ╔═╡ 9ea7fe71-67a0-41ab-bb5b-3d191eb4275c
@savefig "iso_converg_num" compare_vcirc(profiles)

# ╔═╡ 0cba9678-4469-4a14-93cb-417c01beb18d
profiles_softening = [

	"s0.14" => load_profile("1e6/s0.14"),
    "s0.044" => load_profile("1e6/fiducial")[1:6],
    "s0.014" => load_profile("1e6/s0.014"),
 #    "1e7" => load_profile("1e7/fiducial")
    ];

# ╔═╡ 49c4798b-59e3-4908-90f8-e0a6801a64fd
[prof[end].time for (k, prof) in profiles_softening]

# ╔═╡ c25d0e36-7d7c-438e-ad45-61abd3278a75
@savefig "iso_converg_softening" compare_vcirc(profiles_softening)

# ╔═╡ 641ee083-0640-42e6-9822-a1ff1d5aecac
# profiles_accuracy = [

# 	"s0.14" => load_profile("1e6/s0.14"),
#     "s0.044" => load_profile("1e6/fiducial")[1:6],
#     "s0.014" => load_profile("1e6/s0.014"),
#  #    "1e7" => load_profile("1e7/fiducial")
#     ];

# ╔═╡ 1e109d8c-d7e4-4f46-9e39-3adfc881f563
profiles_methods = [
    "fiducial" => load_profile("1e6/fiducial")[1:13],
 	"gadget2" => load_profile("1e6/gadget2"),
#     "fixed timestep" => load_profile("1e6/s0.014"),
#  #    "high accuracy" => load_profile("1e7/fiducial")
];

# ╔═╡ 886b0f81-dca7-47e2-847b-bcd3450d4d07
[prof[end].time for (k, prof) in profiles_methods]

# ╔═╡ 46d7290c-a546-4910-90b3-a6184760a29d
@savefig "iso_converg_methods" compare_vcirc(profiles_methods)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═2abed1bf-9861-49b0-9bea-aa9e9f67829e
# ╠═0f08434b-6f55-4d73-90eb-750e5d979a47
# ╠═7ae386b0-8c0b-43c1-b414-ac0518b73c70
# ╠═2c77e958-c57c-4bb0-94f3-f29414e275cb
# ╠═ba744355-8da6-4caa-a252-ad83ad4b8b13
# ╠═cc50d80b-78b8-4f39-8a41-18ae79bf7287
# ╠═f7cee4fe-e8df-429b-9cb0-f508d60c9788
# ╠═a4b6b203-a42a-48dd-8944-ea150eee376b
# ╠═69516e98-4e15-44a8-b0d4-dc8d8bacd7e8
# ╠═9ea7fe71-67a0-41ab-bb5b-3d191eb4275c
# ╠═0cba9678-4469-4a14-93cb-417c01beb18d
# ╠═49c4798b-59e3-4908-90f8-e0a6801a64fd
# ╠═c25d0e36-7d7c-438e-ad45-61abd3278a75
# ╠═641ee083-0640-42e6-9822-a1ff1d5aecac
# ╠═1e109d8c-d7e4-4f46-9e39-3adfc881f563
# ╠═886b0f81-dca7-47e2-847b-bcd3450d4d07
# ╠═46d7290c-a546-4910-90b3-a6184760a29d
