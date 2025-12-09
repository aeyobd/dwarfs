### A Pluto.jl notebook ###
# v0.20.21

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

# ╔═╡ 8aab441c-3268-449d-b401-94f8e15eb5f3
using CSV, DataFrames

# ╔═╡ cf655228-1528-4e77-9629-07e99984951f
using OrderedCollections

# ╔═╡ 6d2bff6d-b347-49c4-87df-ad58d8a27ff3
include("./style.jl")

# ╔═╡ b47ea4ad-9f71-44a5-abc0-df6a748d5247
CairoMakie.activate!(type=:png)

# ╔═╡ 67dd3488-79e1-4ba2-b9ac-f417765d55de
import TOML

# ╔═╡ 4397ebbd-9928-4a79-9d9f-ad8deb6f6fb4
theme(:Lines)[:cycle] = Cycle([:color]=>:color)

# ╔═╡ 68fdf1ff-2463-4d2a-bd86-29158af57262
theme(:Lines)

# ╔═╡ 257236d1-37f8-4e45-81c4-0407216010f4
function get_R_h(galaxyname)                                                      
	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))
	R_h = obs_props["R_h"]                                                        
	return R_h
end          

# ╔═╡ 40ef5263-bf52-4ad7-8c39-ed1e11c45fc4
function load_profile(galaxyname; algname="fiducial")
	@info galaxyname
	
	filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, 
		"density_profiles/$(algname)_profile.toml")
	
    prof = LilGuys.SurfaceDensityProfile(filename) |> LilGuys.filter_empty_bins

	R_h = get_R_h(galaxyname)
	
	@info "counts = $(sum(prof.counts))"

	Σ_h = 10 .^ LilGuys.lerp(prof.log_R, middle.(prof.log_Sigma))(log10(R_h))
	M_s = Σ_h * R_h .^ 2

	
    prof = LilGuys.scale(prof, 1/R_h, 1/M_s)

	filt = maximum.(LilGuys.error_interval.(prof.log_Sigma)).< 1

	prof.log_R = prof.log_R[filt]
	prof.log_Sigma = prof.log_Sigma[filt]

	prof
end

# ╔═╡ 4ee54995-a31f-4ae1-8204-e55884d786bb
α = LilGuys.R_h(LilGuys.Exp2D(R_s=1))

# ╔═╡ 8adc28b1-b0d0-4735-a1f1-ee3bb82e8ef2
prof_scl = load_profile("sculptor")

# ╔═╡ c3b52d5a-a8b4-4207-a8c7-9d32914aca93
prof_umi = load_profile("ursa_minor")

# ╔═╡ 740d04a4-4327-47cd-b5aa-b71ed94a0610
prof_fornax = load_profile("fornax")

# ╔═╡ dea5ec9f-7f31-42eb-b8f7-6ca5d22e9f0b
begin 
	plummer = LilGuys.Plummer()
	Σ_0 = surface_density(plummer, 1)
	plummer = LilGuys.Plummer(M = 1/surface_density(plummer, 1))
end

# ╔═╡ 8821d235-5622-48c6-ab49-587d0382e583
galaxies = OrderedDict(
	"fornax" => (
		color = COLORS[1],
		marker = :dot,
		label = "Fornax",
		linestyle=:solid,
		#linewidth=2
	),
	"sculptor" => (
		color = COLORS[2], 
		marker=:rect,
		label = "Sculptor",
		linestyle=:solid,
		#linewidth=1.75
	),
	"ursa_minor" => (
		color = COLORS[3],
		marker = :star5,
		label = "Ursa Minor",
		linestyle=:solid,
		#markersize=8,
		#linewidth=1.5
	),
)

# ╔═╡ 552e9438-862f-4710-a7d4-c8d798b5f1aa
galaxynames = [
	"fornax",
	"leo1",
	"leo2",
	"carina",
	"sextans1",
	"draco",	
	"canes_venatici1",
	"crater2",
]

# ╔═╡ 4db3ed1d-ef02-4f28-8c9b-b6c6d491199f
[print(galaxy, ", ") for galaxy in galaxynames]

# ╔═╡ f4e8b66c-5f18-45fd-8859-32479d7227bc
begin 
	profiles = OrderedDict(name => load_profile(name) for name in galaxynames)
	profiles["sculptor"] = prof_scl
	profiles["ursa_minor"] = prof_umi
end

# ╔═╡ d34dfbe4-c80a-46e8-9082-30148cbcbd0b
profiles["sculptor"]

# ╔═╡ ae186cf2-6fd5-4958-acc0-8a87eab45358
p08_notides = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/all/data", "penarrubia08_notides.csv"), DataFrame)

# ╔═╡ 180d22b5-21f0-4559-9696-03887675d9f9
p08_r_scale = 0.7

# ╔═╡ fe81daf4-d2a9-484f-8ba4-9d4f7348d3fd
p08_m_scale = 5 * p08_r_scale

# ╔═╡ eb85be07-d72e-4d22-8cef-a8521d465a4c
p08_tides = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/all/data", "penarrubia08_apo.csv"), DataFrame)

# ╔═╡ 5eb61a4e-063f-4c11-8dea-07be62cee089
begin 
	prof_i = LilGuys.SurfaceDensityProfile(
		log_R=p08_notides.x,
		log_Sigma=p08_notides." y",
		R_units="", 
		log_R_bins=[]
	) 
	
	prof_f = LilGuys.SurfaceDensityProfile(
		log_R=p08_tides.x,
		log_Sigma=p08_tides." y",
		R_units="", 
		log_R_bins=[]
	)

	prof_i = LilGuys.scale(prof_i, p08_r_scale, p08_m_scale)

	prof_f = LilGuys.scale(prof_f, p08_r_scale, p08_m_scale)

end

# ╔═╡ ec3d92b1-76b1-41cd-b52c-0268c4a4a584
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		ylabel = L"log surface density / $\Sigma_h$",
				xticks = -2:1:1,
		limits=(-1.5, nothing, -4.5, 1.1)

	)

    x = LinRange(-1.5, 1.2, 1000)
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y, color=:black, label="exponential")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
   # lines!(x, y, color=:black, label="Plummer", linestyle=:dash)

	for galaxy in galaxynames
		prof = profiles[galaxy]
		lines!(prof.log_R, prof.log_Sigma, color=COLORS[1], alpha=0.5, label="classical dwarfs", linestyle=:solid)
	end
	axislegend(position=:lb, merge=true, unique=true, patchsize=(48, 6))


	# residuals
	ax_res = Axis(fig[2,1], 
		ylabel = "residual",
		xlabel = L"log Radius / $R_h$",
		xticks = -2:1:1,
		limits=(-1.5, 1.1, -1.5, 1.5),
		yticks = -2:1:2
	)

	f(x) = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
	
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Exp2D")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
  #  lines!(x, y .- f(x), color=:black, label="Plummer", linestyle=:dash)


	for galaxy in galaxynames
		prof = profiles[galaxy]
		lines!(prof.log_R, prof.log_Sigma .- f(prof.log_R), color=COLORS[1], alpha=0.5, linestyle=:solid)
	end


	linkxaxes!(ax, ax_res)#, ax_res_plummer)
	rowgap!(fig.layout, 0)
	rowsize!(fig.layout, 2, Relative(1/5))
	#rowsize!(fig.layout, 3, Relative(1/5))
	hidexdecorations!(ax, ticks=false, minorticks=false)
	#hidexdecorations!(ax_res, ticks=false, minorticks=false)

	@savefig "classical_dwarfs"

	fig
end


# ╔═╡ 99426526-4e5b-42ae-8e08-b1e93d5f7fe6
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		ylabel = L"log surface density / $\Sigma_h$",
				xticks = -2:1:1,
		limits=(-1.5, nothing, -4.5, 1.1)

	)

    x = LinRange(-1.5, 1.2, 1000)
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y, color=:black, label="exponential")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
   #lines!(x, y, color=:black, label="Plummer", linestyle=:dash)

	for galaxy in galaxynames
		prof = profiles[galaxy]
		lines!(prof.log_R, prof.log_Sigma, color=COLORS[1], alpha=0.5, label="classical dwarfs", linestyle=:solid)
	end
	for galaxy in ["sculptor", "ursa_minor"]
		prof = profiles[galaxy]
		scatterlines!(prof.log_R, prof.log_Sigma; galaxies[galaxy]...)
	end

	axislegend(position=:lb, merge=true, unique=true, patchsize=(48, 6))

	# residuals
	ax_res = Axis(fig[2,1], 
		ylabel = "residual",
		xlabel = L"log Radius / $R_h$",
		xticks = -2:1:1,
		limits=(-1.5, 1.1, -1.5, 1.5),
		yticks = -2:1:2
	)

	f(x) = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
	
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Exp2D")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    #lines!(x, y .- f(x), color=:black, label="Plummer", linestyle=:dash)


	for galaxy in galaxynames
		prof = profiles[galaxy]
		lines!(prof.log_R, prof.log_Sigma .- f(prof.log_R), color=COLORS[1], alpha=0.5, linestyle=:solid)
	end
	
	for galaxy in ["sculptor", "ursa_minor"]
		prof = profiles[galaxy]
		scatterlines!(prof.log_R, prof.log_Sigma .- f(prof.log_R); galaxies[galaxy]...)
	end

	linkxaxes!(ax, ax_res)#, ax_res_plummer)
	rowgap!(fig.layout, 0)
	rowsize!(fig.layout, 2, Relative(1/5))
	#rowsize!(fig.layout, 3, Relative(1/5))
	hidexdecorations!(ax, ticks=false, minorticks=false)
	#hidexdecorations!(ax_res, ticks=false, minorticks=false)

	@savefig "classical_dwarfs_vs_scl_umi"

	fig
end


# ╔═╡ 86ecd227-8524-4a4b-8f7d-bb0a0fc97257
# let
# 	fig = Figure()

# 	ax = Axis(fig[1,1], 
# 		ylabel = L"log surface density / $\Sigma_h$",
# 				xticks = -2:1:1,
# 		limits=(-1.5, nothing, -4.5, 1.1)

# 	)

#     x = LinRange(-1.5, 1.2, 1000)
#     y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
#    # lines!(x, y, color=:black, label="exponential")
    

#     y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
#    #lines!(x, y, color=:black, label="Plummer", linestyle=:dash)

# 	for galaxy in galaxynames[1:1]
# 		prof = profiles[galaxy]
# 		lines!(prof.log_R, prof.log_Sigma, color=COLORS[1], alpha=0.5, label="Fornax")
# 	end
# 	for galaxy in ["sculptor", "ursa_minor"]
# 		prof = profiles[galaxy]
# 		scatterlines!(prof.log_R, prof.log_Sigma; galaxies[galaxy]...)
# 	end

# 	lines!(prof_i.log_R, prof_i.log_Sigma, label="P+08 initial", linestyle=:dot, color=:black)

	
# 	axislegend(position=:lb, merge=true, unique=true, patchsize=(48, 6))

# 	# residuals
# 	ax_res = Axis(fig[2,1], 
# 		ylabel = "residual",
# 		xlabel = L"log Radius / $R_h$",
# 		xticks = -2:1:1,
# 		limits=(-1.5, 1.1, -1.5, 1.5),
# 		yticks = -2:1:2
# 	)

# 	f(x) = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
	
#     y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
#     lines!(x, y .- f(x), color=:black, label="Exp2D")
    

#     y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
#     #lines!(x, y .- f(x), color=:black, label="Plummer", linestyle=:dash)


# 	for galaxy in galaxynames[1:1]
# 		prof = profiles[galaxy]
# 		lines!(prof.log_R, prof.log_Sigma .- f(prof.log_R), color=COLORS[1], alpha=0.5)
# 	end
	
# 	for galaxy in ["sculptor", "ursa_minor"]
# 		prof = profiles[galaxy]
# 		scatterlines!(prof.log_R, prof.log_Sigma .- f(prof.log_R); galaxies[galaxy]...)
# 	end

# 	lines!(prof_i.log_R, prof_i.log_Sigma .- f(prof_i.log_R), label="no tides", linestyle=:dot, color=:black)

# 	linkxaxes!(ax, ax_res)#, ax_res_plummer)
# 	rowgap!(fig.layout, 0)
# 	rowsize!(fig.layout, 2, Relative(1/5))
# 	#rowsize!(fig.layout, 3, Relative(1/5))
# 	hidexdecorations!(ax, ticks=false, minorticks=false)
# 	#hidexdecorations!(ax_res, ticks=false, minorticks=false)

# 	@savefig "classical_dwarfs_vs_scl_umi_p08_i"

# 	fig
# end


# ╔═╡ ed718819-40a0-4568-92b3-0e86808729e2
# let
# 	fig = Figure()

# 	ax = Axis(fig[1,1], 
# 		ylabel = L"log surface density / $\Sigma_h$",
# 				xticks = -2:1:1,
# 		limits=(-1.5, nothing, -4.5, 1.1)

# 	)

#     x = LinRange(-1.5, 1.2, 1000)
#     y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
#    # lines!(x, y, color=:black, label="exponential")
    

#     y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
#    #lines!(x, y, color=:black, label="Plummer", linestyle=:dash)

# 	for galaxy in galaxynames[1:1]
# 		prof = profiles[galaxy]
# 		lines!(prof.log_R, prof.log_Sigma, color=COLORS[1], alpha=0.5, label="Fornax")
# 	end
# 	for galaxy in ["sculptor", "ursa_minor"]
# 		prof = profiles[galaxy]
# 		scatterlines!(prof.log_R, prof.log_Sigma; galaxies[galaxy]...)
# 	end

# 	lines!(prof_i.log_R, prof_i.log_Sigma, label="P+08 initial", linestyle=:dot, color=:black)
# 	lines!(prof_f.log_R, prof_f.log_Sigma, label="P+08 with tides", linewidth=6, color=:black)

	
# 	axislegend(position=:lb, merge=true, unique=true, patchsize=(48, 6))

# 	# residuals
# 	ax_res = Axis(fig[2,1], 
# 		ylabel = "residual",
# 		xlabel = L"log Radius / $R_h$",
# 		xticks = -2:1:1,
# 		limits=(-1.5, 1.1, -1.5, 1.5),
# 		yticks = -2:1:2
# 	)

# 	f(x) = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
	
#     y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
#     lines!(x, y .- f(x), color=:black, label="Exp2D")
    

#     y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
#     #lines!(x, y .- f(x), color=:black, label="Plummer", linestyle=:dash)


# 	for galaxy in galaxynames[1:1]
# 		prof = profiles[galaxy]
# 		lines!(prof.log_R, prof.log_Sigma .- f(prof.log_R), color=COLORS[1], alpha=0.5)
# 	end
	
# 	for galaxy in ["sculptor", "ursa_minor"]
# 		prof = profiles[galaxy]
# 		scatterlines!(prof.log_R, prof.log_Sigma .- f(prof.log_R); galaxies[galaxy]...)
# 	end

# 	lines!(prof_i.log_R, prof_i.log_Sigma .- f(prof_i.log_R), label="no tides", linestyle=:dot, color=:black)
# 	lines!(prof_f.log_R, prof_f.log_Sigma .- f(prof_f.log_R), label="tides", linewidth=6, color=:black)

# 	linkxaxes!(ax, ax_res)#, ax_res_plummer)
# 	rowgap!(fig.layout, 0)
# 	rowsize!(fig.layout, 2, Relative(1/5))
# 	#rowsize!(fig.layout, 3, Relative(1/5))
# 	hidexdecorations!(ax, ticks=false, minorticks=false)
# 	#hidexdecorations!(ax_res, ticks=false, minorticks=false)

# 	@savefig "classical_dwarfs_vs_scl_umi_p08_i_f"

# 	fig
# end


# ╔═╡ 0599cad0-64b4-4d6d-9669-b660dc659313
md"""
# Extra
"""

# ╔═╡ 45bb2bf2-66de-4834-adb2-725277fa6180
function plot_fit(prof, title="")
	fig = Figure()

	ax = Axis(fig[1,1], 
		ylabel = L"\log\,\Sigma\ / \ \Sigma_h",
			  title=title
				# xticks = -2:1:1,
		#limits=(-1, nothing, -4.2, 1.1)

	)

    x = LinRange(extrema(prof.log_R_bins)..., 1000)
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y, color=:black, label="Exp2D")

	vlines!(log10(3/α), color=COLORS[2], label = "fit R max")
	LilGuys.plot_log_Σ!(ax, prof, label = "observations")

	axislegend(position=:lb)

	
	# residuals
	ax_res = Axis(fig[2,1], 
		ylabel = L"\delta \log\,\Sigma",
		xlabel = L"\log\,R\ / \ R_h",
		# xticks = -2:1:1,
		limits=(-1, 1.0, -1.2, 1.2)
	)
	f(x) = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))

	y = prof.log_Sigma .- f(prof.log_R)

	errorscatter!(prof.log_R, y, yerror=error_interval.(y))
	hlines!(0, color=:black)
	vlines!(log10(3/α), color=COLORS[2])

	linkxaxes!(ax, ax_res,)

	rowsize!(fig.layout, 2, Relative(1/4))
	hidexdecorations!(ax, ticks=false, minorticks=false)
	fig

end

# ╔═╡ e16b1441-f184-4fd2-8b68-52801d0a0960
plot_fit(profiles["fornax"], "fornax")

# ╔═╡ 356c06df-79ad-4313-a7a2-a4506a651700
plot_fit(prof_scl, "sculptor")

# ╔═╡ 89a375e4-b810-4899-aff1-5da81e5ca067
plot_fit(profiles["leo1"], "Leo I")

# ╔═╡ cf78a1a7-1f03-436c-9c76-45047e499e84
plot_fit(profiles["carina"], "Carina")

# ╔═╡ eb5b28a9-7873-4039-822f-7f40be719098
plot_fit(profiles["leo2"], "Leo II",  )

# ╔═╡ 8286e2a4-2bd0-47a0-badc-b55697e2a8fa
plot_fit(profiles["sextans1"], "Sextans1")

# ╔═╡ 4514e8fe-abd6-4b0b-9e64-a9cb7f2c782d
plot_fit(prof_umi, "Ursa Minor")

# ╔═╡ c20eb47f-e042-49f2-90a5-39056f5f5284
plot_fit(profiles["draco"], "Draco")

# ╔═╡ 1e81485d-283d-4a83-944f-55c9e3e9883d
plot_fit(profiles["crater2"], "Crater II")

# ╔═╡ 4eb83e39-1988-4a20-90a8-5c0a822dd162
plot_fit(profiles["canes_venatici1"], "CVe II")

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═8aab441c-3268-449d-b401-94f8e15eb5f3
# ╠═b47ea4ad-9f71-44a5-abc0-df6a748d5247
# ╠═67dd3488-79e1-4ba2-b9ac-f417765d55de
# ╠═4397ebbd-9928-4a79-9d9f-ad8deb6f6fb4
# ╠═68fdf1ff-2463-4d2a-bd86-29158af57262
# ╠═257236d1-37f8-4e45-81c4-0407216010f4
# ╠═40ef5263-bf52-4ad7-8c39-ed1e11c45fc4
# ╠═4ee54995-a31f-4ae1-8204-e55884d786bb
# ╠═8adc28b1-b0d0-4735-a1f1-ee3bb82e8ef2
# ╠═c3b52d5a-a8b4-4207-a8c7-9d32914aca93
# ╠═740d04a4-4327-47cd-b5aa-b71ed94a0610
# ╠═6d2bff6d-b347-49c4-87df-ad58d8a27ff3
# ╠═dea5ec9f-7f31-42eb-b8f7-6ca5d22e9f0b
# ╠═d34dfbe4-c80a-46e8-9082-30148cbcbd0b
# ╠═8821d235-5622-48c6-ab49-587d0382e583
# ╠═552e9438-862f-4710-a7d4-c8d798b5f1aa
# ╠═4db3ed1d-ef02-4f28-8c9b-b6c6d491199f
# ╠═cf655228-1528-4e77-9629-07e99984951f
# ╠═f4e8b66c-5f18-45fd-8859-32479d7227bc
# ╠═ae186cf2-6fd5-4958-acc0-8a87eab45358
# ╠═5eb61a4e-063f-4c11-8dea-07be62cee089
# ╠═180d22b5-21f0-4559-9696-03887675d9f9
# ╠═fe81daf4-d2a9-484f-8ba4-9d4f7348d3fd
# ╠═eb85be07-d72e-4d22-8cef-a8521d465a4c
# ╠═ec3d92b1-76b1-41cd-b52c-0268c4a4a584
# ╠═99426526-4e5b-42ae-8e08-b1e93d5f7fe6
# ╠═86ecd227-8524-4a4b-8f7d-bb0a0fc97257
# ╠═ed718819-40a0-4568-92b3-0e86808729e2
# ╟─0599cad0-64b4-4d6d-9669-b660dc659313
# ╠═45bb2bf2-66de-4834-adb2-725277fa6180
# ╠═e16b1441-f184-4fd2-8b68-52801d0a0960
# ╠═356c06df-79ad-4313-a7a2-a4506a651700
# ╠═89a375e4-b810-4899-aff1-5da81e5ca067
# ╠═cf78a1a7-1f03-436c-9c76-45047e499e84
# ╠═eb5b28a9-7873-4039-822f-7f40be719098
# ╠═8286e2a4-2bd0-47a0-badc-b55697e2a8fa
# ╠═4514e8fe-abd6-4b0b-9e64-a9cb7f2c782d
# ╠═c20eb47f-e042-49f2-90a5-39056f5f5284
# ╠═1e81485d-283d-4a83-944f-55c9e3e9883d
# ╠═4eb83e39-1988-4a20-90a8-5c0a822dd162
