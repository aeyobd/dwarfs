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

# ╔═╡ 0c60047d-469d-49bb-b6cc-3ecc5a5ad8c2
include("./paper_style.jl")

# ╔═╡ f926569c-5d14-4236-be29-c11869bcb7d0
CairoMakie.activate!(type="svg", pt_per_unit=2)

# ╔═╡ 619d2e6b-f3e0-41bd-8235-5db53ed32517
import TOML

# ╔═╡ c8981c9d-dbd2-411a-bf1e-a769592f0929


# ╔═╡ 507bf30f-9b76-4255-8df6-ca96ccbd541b
obs_dir = joinpath(ENV["DWARFS_ROOT"], "observations")

# ╔═╡ 24505198-6110-4398-abe3-158499d7894d
galaxies = [
	#"antlia2",
 #   "bootes1",
  #  "bootes3",
   # "crater2",
    "carina",
    #"canes_venatici1",
    "draco",
    "fornax",
    "leo2",
    "leo1",
    #"reticulum2",
    "sextans1",
    "sculptor",
    "ursa_minor",
    ]

# ╔═╡ c5a1bf8b-76c1-4fa4-bc46-16d5be93f74b
function load_density_fit(galaxyname; algname=nothing)
	if galaxyname ∈ ["crater2", "antlia2"]
		algname = "mcmc_hist_fast"
	else
		algname = "mcmc_hist"
	end
	filename = joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, 
		"density_profiles/jax_profile.toml")


	density_fit = TOML.parsefile(filename * "_density_fits.toml")
end

# ╔═╡ 4cca1a03-d542-4cf5-a7fa-d734e40e7ec6
properties = [TOML.parsefile(obs_dir * "/$galaxy/observed_properties.toml") for galaxy in galaxies]

# ╔═╡ a448099c-cada-4a28-9403-1b35148aeda8
density_fits = load_density_fit.(galaxies)

# ╔═╡ bd76e24f-e722-4d99-9bbd-74b99a74c3fd
properties[1]

# ╔═╡ ebc625bd-f0e7-4b64-979f-11182d7713ae
Mv = [prop["Mv"] for prop in properties]

# ╔═╡ e13268f7-b7b0-4429-95e2-3ee3c1ae9452
ellipticity = [prop["ellipticity"] for prop in properties]

# ╔═╡ 9439fcd0-8a1a-4ac9-9e80-8a6b5a8e9a0c
α = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 19af5904-5d5a-4a7a-81b0-edf9fef78441
R_h = α * [fit["log_R_s_exp2d_inner"] for fit in density_fits]

# ╔═╡ a8afccf2-a300-4e0c-b807-33153b1dbe81
n_sersic = [fit["n_sersic"] for fit in density_fits]

# ╔═╡ 346bcb0f-0f84-4152-9f1b-43f6d5eae576
function get_err(properties, key)
	errs = []

	for i in eachindex(properties)
		props = properties[i]

		if key*"_em" ∈ keys(props)
			push!(errs, (props[key*"_em"], props[key*"_ep"]) )
		else
			push!(errs, (props[key*"_err"], props[key*"_err"]) )
		end

	end

	return errs
end

# ╔═╡ 50765eea-f657-4c84-86ee-50617897a22b
n_sersic_err = get_err(density_fits, "n_sersic")

# ╔═╡ 733ba2e9-6e0b-450d-b76d-3a29b2dcf004
R_h_err =  [α .* e for e in  get_err(density_fits, "log_R_s_exp2d_inner")]

# ╔═╡ 8683840e-a7b2-47e9-ac15-a0baf45dff48
Mv_err = get_err(properties, "Mv")

# ╔═╡ 1cd44536-044f-4e62-a2d1-ada28e046444
let
	fig = Figure()
	ax = Axis(fig[1,1],
	    xlabel = "absolute magnitude (Mv)",
	    ylabel = L"Sérsic $n$",
		limits = (-15, nothing, 0.5, 1.5)
	    )
	hlines!(1, color=(:black, 0.2), linestyle=:dash)

	errorscatter!(Mv, n_sersic, xerror=Mv_err, yerror=n_sersic_err, color=(:black, 0.2))
	
	p = scatter!(Mv, n_sersic, color=ellipticity, markersize=5)
	
	text!(Mv, n_sersic, text=galaxies, fontsize=10, color=p.color, rotation=0)
	LilGuys.hide_grid!(ax)


	Colorbar(fig[1,2], p, label="ellipticity")

	ax.xreversed = true

	@savefig "classical_dwarf_scatter"
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═0c60047d-469d-49bb-b6cc-3ecc5a5ad8c2
# ╠═f926569c-5d14-4236-be29-c11869bcb7d0
# ╠═619d2e6b-f3e0-41bd-8235-5db53ed32517
# ╠═c8981c9d-dbd2-411a-bf1e-a769592f0929
# ╠═507bf30f-9b76-4255-8df6-ca96ccbd541b
# ╠═24505198-6110-4398-abe3-158499d7894d
# ╠═c5a1bf8b-76c1-4fa4-bc46-16d5be93f74b
# ╠═4cca1a03-d542-4cf5-a7fa-d734e40e7ec6
# ╠═a448099c-cada-4a28-9403-1b35148aeda8
# ╠═bd76e24f-e722-4d99-9bbd-74b99a74c3fd
# ╠═ebc625bd-f0e7-4b64-979f-11182d7713ae
# ╠═e13268f7-b7b0-4429-95e2-3ee3c1ae9452
# ╠═9439fcd0-8a1a-4ac9-9e80-8a6b5a8e9a0c
# ╠═19af5904-5d5a-4a7a-81b0-edf9fef78441
# ╠═a8afccf2-a300-4e0c-b807-33153b1dbe81
# ╠═50765eea-f657-4c84-86ee-50617897a22b
# ╠═733ba2e9-6e0b-450d-b76d-3a29b2dcf004
# ╠═346bcb0f-0f84-4152-9f1b-43f6d5eae576
# ╠═8683840e-a7b2-47e9-ac15-a0baf45dff48
# ╠═1cd44536-044f-4e62-a2d1-ada28e046444
