### A Pluto.jl notebook ###
# v0.20.18

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

# ╔═╡ cf655228-1528-4e77-9629-07e99984951f
using OrderedCollections

# ╔═╡ 6d2bff6d-b347-49c4-87df-ad58d8a27ff3
include("./paper_style.jl")

# ╔═╡ 67dd3488-79e1-4ba2-b9ac-f417765d55de
import TOML

# ╔═╡ 9f06a5c4-a8ae-49d2-991e-1c1576361700
function get_R_h(galaxyname)
	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml")) |> LilGuys.collapse_errors
	R_h = obs_props["R_h"]
	return R_h
end

# ╔═╡ c869bbaa-bed3-48aa-8e8f-5dcfc833e2ac
get_R_h("leo1")

# ╔═╡ be13884c-fe1b-4b30-95a0-983ff73eaa0c
get_R_h("fornax")

# ╔═╡ 0d5d8a0f-97df-4a46-876e-3f11f52ceecf
function get_R_h_inner(galaxyname)
	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml")) |> LilGuys.collapse_errors
	R_h = obs_props["R_h_inner"]
	return R_h
end

# ╔═╡ 40ef5263-bf52-4ad7-8c39-ed1e11c45fc4
function load_profile(galaxyname; algname="jax", inner=false)
	@info galaxyname
	
	filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, 
		"density_profiles/$(algname)_eqw_profile.toml")
	
    prof = LilGuys.SurfaceDensityProfile(filename) |> LilGuys.filter_empty_bins

	if inner
		R_h = get_R_h_inner(galaxyname)
	else
		R_h = get_R_h(galaxyname)
	end
	
	@info "R_h = $R_h arcmin"
	@info "counts = $(sum(prof.counts))"

	Σ_h = 10 .^ LilGuys.lerp(prof.log_R, middle.(prof.log_Sigma))(log10(middle(R_h)))
	M_s = Σ_h * middle(R_h) .^ 2

	
    prof = LilGuys.scale(prof, 1/middle(R_h), 1/M_s)

	filt = maximum.(LilGuys.error_interval.(prof.log_Sigma)).< 1

	prof.log_Sigma[.!filt] .= Measurement(NaN, NaN)

	prof |> LilGuys.filter_empty_bins
end

# ╔═╡ 4ee54995-a31f-4ae1-8204-e55884d786bb
α = LilGuys.R_h(LilGuys.Exp2D(R_s=1))

# ╔═╡ 8adc28b1-b0d0-4735-a1f1-ee3bb82e8ef2
prof_scl = load_profile("sculptor", algname="jax_2c")

# ╔═╡ c3b52d5a-a8b4-4207-a8c7-9d32914aca93
prof_umi = load_profile("ursa_minor", algname="jax_2c")

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
		marker = :circle,
		label = "Fornax",
		linewidth=2
	),
	"sculptor" => (
		color = COLORS[2], 
		marker=:rect,
		label = "Sculptor",
		linewidth=1.75
	),
	"ursa_minor" => (
		color = COLORS[3],
		marker = :star5,
		label = "Ursa Minor",
		markersize=8,
		linewidth=1.5
	),
)

# ╔═╡ 07e1390f-2789-4e20-80d6-b2cb4631be5c
galaxies_scatter = OrderedDict(
	"fornax" => (
		color = COLORS[1],
		marker = :circle,
		label = "Fornax",
	),
	"sculptor" => (
		color = COLORS[2], 
		marker=:rect,
		label = "Sculptor",
	),
	"ursa_minor" => (
		color = COLORS[3],
		marker = :utriangle,
		label = "Ursa Minor",
		markersize=8,
	),
)

# ╔═╡ 552e9438-862f-4710-a7d4-c8d798b5f1aa
galaxynames = [
	"fornax",
	"leo1",
	# "antlia2", this one is unreliable...
	 "leo2",
	 "carina",
	 "draco",	
	"canes_venatici1",
	 "sextans1",
	"crater2",
]

# ╔═╡ f4e8b66c-5f18-45fd-8859-32479d7227bc
begin 
	profiles = OrderedDict(name => load_profile(name) for name in galaxynames)
	profiles["sculptor"] = prof_scl
	profiles["ursa_minor"] = prof_umi
end

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		ylabel = L"\log\,\Sigma\ / \ \Sigma_h",
				xticks = -2:1:1,
		limits=(-1.5, nothing, -4.2, 1.1)

	)

    x = LinRange(-1.4, 1, 1000)
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y, color=:black, label="Exp2D")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    lines!(x, y, color=:black, label="Plummer", linestyle=:dash)


	for (galaxy, prof) in profiles
		lines!(prof.log_R, prof.log_Sigma, color=COLORS[1], alpha=0.5, label="classical dwarfs")
	end


	scatterlines!(prof_scl.log_R, prof_scl.log_Sigma, color=COLORS[2], label="Sculptor", )
	scatterlines!(prof_umi.log_R, prof_umi.log_Sigma, color=COLORS[3], label="Ursa Minor", )

	axislegend(position=:lb, merge=true, unique=true)



	# residuals
	ax_res = Axis(fig[2,1], 
		ylabel = L"\delta \log\,\Sigma",
		xlabel = L"\log\,R\ / \ R_h",
		xticks = -2:1:1,
		limits=(-1.5, 1.1, -1.2, 3)
	)

	f(x) = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
	
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Exp2D")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Plummer", linestyle=:dash)

	
	for (galaxy, prof) in profiles
		lines!(prof.log_R, prof.log_Sigma .- f(prof.log_R), color=COLORS[1], alpha=0.5)
	end

	scatterlines!(prof_scl.log_R, prof_scl.log_Sigma .- f(prof_scl.log_R), color=COLORS[2])
	scatterlines!(prof_umi.log_R, prof_umi.log_Sigma .- f(prof_umi.log_R), color=COLORS[3])


	linkxaxes!(ax, ax_res)
	rowgap!(fig.layout, 0)
	rowsize!(fig.layout, 2, Relative(1/4))
	hidexdecorations!(ax, ticks=false, minorticks=false)

	@savefig "classical_dwarf_profiles"
	fig
end

# ╔═╡ d34dfbe4-c80a-46e8-9082-30148cbcbd0b
profiles["sculptor"]

# ╔═╡ ceb62c28-85b7-424e-a6e4-243905333a86
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		ylabel = L"\log\,\Sigma\ / \ \Sigma_h",
				xticks = -2:1:1,
		limits=(-1.5, nothing, -4.2, 1.1)

	)

    x = LinRange(-1.5, 1.2, 1000)
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y, color=:black, label="Exp2D")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    lines!(x, y, color=:black, label="Plummer", linestyle=:dash)

	for (galaxy, kwargs) in galaxies
		prof = profiles[galaxy]
		scatterlines!(prof.log_R, prof.log_Sigma; kwargs...)
	end

	axislegend(position=:lb, merge=true, unique=true, patchsize=(24, 6))



	# residuals
	ax_res = Axis(fig[2,1], 
		ylabel = L"\delta\log\,\Sigma_\textrm{Exp}",
		xlabel = L"\log\,R\ / \ R_h",
		xticks = -2:1:1,
		limits=(-1.5, 1.6, -1.2, 1.2),
				  ylabelsize=10
	)

	f(x) = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
	
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Exp2D")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Plummer", linestyle=:dash)


	for (galaxy, kwargs) in galaxies
		prof = profiles[galaxy]
		scatterlines!(prof.log_R, prof.log_Sigma .- f(prof.log_R); kwargs...)
	end


	# residuals plummer
	ax_res_plummer = Axis(fig[3,1], 
		ylabel = L"\delta\log\,\Sigma_\textrm{plummer}",
		xlabel = L"\log\,R\ / \ R_h",
		xticks = -2:1:1,
		limits=(-1.5, 1.2, -1.2, 1.2),
						  ylabelsize=10
	)

	f(x) = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
	
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Exp2D")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Plummer", linestyle=:dash)

	for (galaxy, kwargs) in galaxies
		prof = profiles[galaxy]
		scatterlines!(prof.log_R, prof.log_Sigma .- f(prof.log_R); kwargs...)
	end

	
	linkxaxes!(ax, ax_res, ax_res_plummer)
	rowgap!(fig.layout, 0)
	rowsize!(fig.layout, 2, Relative(1/5))
	rowsize!(fig.layout, 3, Relative(1/5))
	hidexdecorations!(ax, ticks=false, minorticks=false)
	hidexdecorations!(ax_res, ticks=false, minorticks=false)

	@savefig "scl_umi_fornax_exp_fit"
	fig
end

# ╔═╡ 5b5c426e-ce27-48bd-a3ac-7f23e7663ed7
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		ylabel = L"\log\,\Sigma\ / \ \Sigma_h",
		xlabel = L"\log\,R\ / \ R_h",

				xticks = -2:1:1,
		limits=(-1.5, 1.2, -5, 1.1)

	)

    x = LinRange(-1.5, 1.2, 1000)
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y, color=:black, label="Exp2D")
    

	for (galaxy, kwargs) in galaxies_scatter
		prof = profiles[galaxy]
			errorscatter!(ax, prof.log_R, prof.log_Sigma, yerror=LilGuys.error_interval.(prof.log_Sigma), ; kwargs...)
	end

	axislegend(position=:lb, merge=true, unique=true, patchsize=(24, 6))



	@savefig "scl_umi_fornax_exp"
	fig
end

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
		# limits=(-1, 1.0, -1.2, 1.2)
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

# ╔═╡ dece1fde-c036-4b37-814b-e63b55bb9b4f
plot_fit(profiles["crater2"], "Crater II")

# ╔═╡ 13de9235-2abb-4ba5-8e5b-87111191722c
plot_fit(profiles["canes_venatici1"], "CVn I")

# ╔═╡ f706fc01-fade-4d4f-97d3-ce08264ec680
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		ylabel = L"\log\,\Sigma\ / \ \Sigma_h",
				xticks = -2:1:1,
		limits=(nothing, nothing, -4.2, 1.2)

	)


	for (galaxy, prof) in profiles
		lines!(prof.log_R, prof.log_Sigma, color=COLORS[1], alpha=0.2, label="classical dwarfs")
	end


    x = LinRange(-1.4, 1, 1000)
    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y, color=:black, label="Sérsic (n=1)")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    lines!(x, y, color=:black, label="Plummer", linestyle=:dash)

	axislegend(position=:lb, merge=true, unique=true)



	# residuals
	ax_res = Axis(fig[2,1], 
		ylabel = L"\delta \log\,\Sigma",
		xlabel = L"\log\,R\ / \ R_h",
		xticks = -2:1:1,
		limits=(-1.4, 1, -1, 2)
	)

	f(x) = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
	
	for (galaxy, prof) in profiles
		lines!(prof.log_R, prof.log_Sigma .- f(prof.log_R), color=COLORS[1], alpha=0.2)
	end


    y = log10.(LilGuys.surface_density.(LilGuys.Sersic(n=1), exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Sérsic (n=1)")
    

    y = log10.(LilGuys.surface_density.(plummer, exp10.(x)))
    lines!(x, y .- f(x), color=:black, label="Plummer", linestyle=:dash)


	linkxaxes!(ax, ax_res)
	rowgap!(fig.layout, 0)
	rowsize!(fig.layout, 2, Relative(1/4))
	hidexdecorations!(ax, ticks=false, minorticks=false)

	# @savefig "classical_dwarf_profiles"
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═67dd3488-79e1-4ba2-b9ac-f417765d55de
# ╠═9f06a5c4-a8ae-49d2-991e-1c1576361700
# ╠═c869bbaa-bed3-48aa-8e8f-5dcfc833e2ac
# ╠═be13884c-fe1b-4b30-95a0-983ff73eaa0c
# ╠═0d5d8a0f-97df-4a46-876e-3f11f52ceecf
# ╠═40ef5263-bf52-4ad7-8c39-ed1e11c45fc4
# ╠═4ee54995-a31f-4ae1-8204-e55884d786bb
# ╠═8adc28b1-b0d0-4735-a1f1-ee3bb82e8ef2
# ╠═c3b52d5a-a8b4-4207-a8c7-9d32914aca93
# ╠═740d04a4-4327-47cd-b5aa-b71ed94a0610
# ╠═6d2bff6d-b347-49c4-87df-ad58d8a27ff3
# ╠═dea5ec9f-7f31-42eb-b8f7-6ca5d22e9f0b
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═d34dfbe4-c80a-46e8-9082-30148cbcbd0b
# ╠═8821d235-5622-48c6-ab49-587d0382e583
# ╠═07e1390f-2789-4e20-80d6-b2cb4631be5c
# ╠═552e9438-862f-4710-a7d4-c8d798b5f1aa
# ╠═cf655228-1528-4e77-9629-07e99984951f
# ╠═f4e8b66c-5f18-45fd-8859-32479d7227bc
# ╠═ceb62c28-85b7-424e-a6e4-243905333a86
# ╠═5b5c426e-ce27-48bd-a3ac-7f23e7663ed7
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
# ╠═dece1fde-c036-4b37-814b-e63b55bb9b4f
# ╠═13de9235-2abb-4ba5-8e5b-87111191722c
# ╠═f706fc01-fade-4d4f-97d3-ce08264ec680
