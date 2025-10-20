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

	using DataFrames, CSV

end

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
function load_Sigmas(galaxyname)
	outdir = joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "mcmc")
	
	df_chains = CSV.read(joinpath(outdir, "samples.mcmc_hist_fast.csv"), DataFrame)
	df_summary = CSV.read(joinpath(outdir, "summary.mcmc_hist_fast.csv"), DataFrame)
	
	bins = 10 .^ [df_summary.log_R_low; df_summary.log_R_high[end]]
	
	Nbins = length(bins) - 1

	Nc = size(df_chains, 1)
	
	log_Sigmas_fast = Matrix{Float64}(undef, Nbins, Nc)
	areas = diff(π * bins .^2)
	counts = df_summary.N_stars

	for i in 1:Nc
		params = [df_chains[i, "params[$j]"] for j in 1:Nbins]
		f_sat = exp10.(params) ./(1 .+ exp10.(params))
		log_Sigmas_fast[:, i] .= log10.( f_sat .* counts ./ areas)
	end

	return bins, log_Sigmas_fast

end

# ╔═╡ e37864c0-3af8-45df-af73-f5a4852a76ef
function load_jax_prof(galaxyname)	
	local profile = nothing
	
	for filename in ["jax_2c_eqw_profile", "jax_eqw_profile", "jax_profile"]
		path = joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname,  "density_profiles/", filename * ".toml")
		if isfile(path)
			profile = LilGuys.SurfaceDensityProfile(path)
			break
		end
	end
	
	profile |> LilGuys.filter_empty_bins
end

# ╔═╡ 8d39a682-acd5-41d9-af5e-36731987d4a7
galaxy_labels = Dict(
	"fornax" => "Fornax",
	"leo1" => "Leo I",
	"sculptor" => "Sculptor",
	"antlia2" => "Antlia II",
	"leo2" => "Leo II",
	"carina" => "Carina",
	"draco" => "Draco", 
	"ursa_minor" => "Ursa Minor",
	"canes_venatici1" => "Canes Venatici I",
	"sextans1" => "Sextans I",
	"crater2" => "Crater II",
)

# ╔═╡ 4358f380-a61a-4bb5-8b7e-073b386d5080
function plot_galaxy(gs, galaxyname; title="")
	ax = Axis(gs,
		# xlabel = L"$\log\,R_\textrm{ell}$ / arcmin",
		# ylabel = "log surface density",
		limits = (nothing, nothing, -5, nothing),
		title = title
	)

	prof_jax_cut = load_jax_prof(galaxyname)
	bins, log_Sigmas = load_Sigmas(galaxyname)
	skip = 48
	jitter = 0.005
	
	for i in 1:skip:size(log_Sigmas, 2)
		y = log_Sigmas[:, i]
		x = midpoints(log10.(bins))
		scatter!(x .+ jitter * randn(length(x)), y, color=:black, alpha=0.03, markersize=1, label="MCMC samples" => (; alpha=0.5), rasterize=2)
	end

	errorscatter!(prof_jax_cut.log_R, prof_jax_cut.log_Sigma, markersize=3, color=COLORS[2], label="J+24")

	text!(0.0, 0.0, offset=(6, 6), space=:relative, text=galaxy_labels[galaxyname])
	return ax
end

# ╔═╡ 27e1ef30-31bc-4688-b314-36570d9e540e
load_jax_prof("antlia2")

# ╔═╡ 5c437e54-8f9b-4276-99dd-e85a8bfd09de
hidexlabels!() =  hidexdecorations!(ticks=false, minorticks=false)

# ╔═╡ b6bb981b-ea26-4df9-b85c-13e2b7505ada
hideylabels!() =  hideydecorations!(ticks=false, minorticks=false)

# ╔═╡ 8deb9d05-ca0e-4cf3-901e-22f0d768630d
hidelabels!() =  hidedecorations!(ticks=false, minorticks=false)

# ╔═╡ b77c2b3a-0961-4ea8-b909-4ed29bd994c0
theme(:size)[] ./ 72

# ╔═╡ ac2be006-902d-4221-88e5-2f6257dc3593
let
	fig = Figure(size=(6*72, 6*72))
	plot_galaxy(fig[1,1], "fornax")
	hidexlabels!()
	
	plot_galaxy(fig[1,2], "leo1")
	hidelabels!()

	plot_galaxy(fig[2,1], "sculptor")
	hidexlabels!()
	
	ax = plot_galaxy(fig[2,2], "antlia2")
	hidelabels!()

	plot_galaxy(fig[3,1], "leo2")
	
	plot_galaxy(fig[3,2], "carina")
	hideylabels!()


	linkaxes!(fig.content...)


	Label(fig[end+1, :], L"log $R_\textrm{ell}$ / arcmin", )
	Label(fig[:, 0], "log surface density", rotation=π/2, )

	axislegend(ax, position=:rt, unique=true, backgroundcolor=(:white, 0.8))
	
	@savefig "mcmc_histograms"
end

# ╔═╡ 5a2f99c6-416c-4101-a17d-1bbd4e40562b
let 
	fig = Figure(size=(6*72, 6*72))


	ax = plot_galaxy(fig[0,1], "draco")
	hidexlabels!()

	
	plot_galaxy(fig[1,1], "ursa_minor")
	hidexlabels!()
	
	plot_galaxy(fig[1,2], "canes_venatici1")
	hidelabels!()
	
	plot_galaxy(fig[2,1], "sextans1")
	# hidexlabels!()
	
	plot_galaxy(fig[2, 2], "crater2")
	hideylabels!()


	
	linkaxes!(fig.content...)

	Label(fig[end+1, :], L"log $R_\textrm{ell}$ / arcmin", )
	Label(fig[:, 0], "log surface density", rotation=π/2, )
	Legend(fig[0,2], ax, tellwidth=false, unique=true)


	@savefig "mcmc_histograms2"

end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═e37864c0-3af8-45df-af73-f5a4852a76ef
# ╠═4358f380-a61a-4bb5-8b7e-073b386d5080
# ╠═8d39a682-acd5-41d9-af5e-36731987d4a7
# ╠═27e1ef30-31bc-4688-b314-36570d9e540e
# ╠═5c437e54-8f9b-4276-99dd-e85a8bfd09de
# ╠═b6bb981b-ea26-4df9-b85c-13e2b7505ada
# ╠═8deb9d05-ca0e-4cf3-901e-22f0d768630d
# ╠═b77c2b3a-0961-4ea8-b909-4ed29bd994c0
# ╠═ac2be006-902d-4221-88e5-2f6257dc3593
# ╠═5a2f99c6-416c-4101-a17d-1bbd4e40562b
