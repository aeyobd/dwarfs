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
	),
	"sculptor" => (
		color = COLORS[2], 
		marker=:rect,
		label = "Sculptor",
		linestyle=:solid,
		linewidth = 1.5*theme(:linewidth)[],
		markersize=1.5*theme(:markersize)[]

		#linewidth=1.75
	),
	"ursa_minor" => (
		color = COLORS[3],
		marker = :star5,
		label = "Ursa Minor",
		linestyle=:solid,
		linewidth = 1.5*theme(:linewidth)[],
		markersize=1.5*theme(:markersize)[]

		#markersize=8,
		#linewidth=1.5
	),
)

# ╔═╡ b05a2d84-2b5f-4d37-8c45-48e55152c43a
theme(:linewidth)

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
		limits=(-1.5, nothing, -3.5, 0.9)

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
	axislegend(position=:lb, merge=true, unique=true, patchsize=(48, 6), margin=theme(:Legend)[:margin][])


	# residuals
	ax_res = Axis(fig[2,1], 
		ylabel = "residual",
		xlabel = L"log Radius / $R_h$",
		xticks = -2:1:1,
		limits=(-1.5, 1.0, -1.5, 1.5),
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
		limits=(-1.5, nothing, -3.5, 0.9)

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

	axislegend(position=:lb, merge=true, unique=true, patchsize=(48, 6), margin=theme(:Legend)[:margin][])

	# residuals
	ax_res = Axis(fig[2,1], 
		ylabel = "residual",
		xlabel = L"log Radius / $R_h$",
		xticks = -2:1:1,
		limits=(-1.5, 1.0, -1.5, 1.5),
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
# ╠═b05a2d84-2b5f-4d37-8c45-48e55152c43a
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
