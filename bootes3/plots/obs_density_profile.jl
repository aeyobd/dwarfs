### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 02612c68-042b-11f0-2f84-e9431cc0ee83
begin
	using Pkg; Pkg.activate()

	using PlutoUI
	using CairoMakie
	using Arya

end

# ╔═╡ cb95f753-b139-4150-bed9-75cd791ceffc
using DataFrames: DataFrame

# ╔═╡ abead192-161f-4375-899f-2c769d40cbfb
using OrderedCollections

# ╔═╡ 898d2c8b-b761-4a69-b561-658a644f44df
begin
	using LilGuys
	FIGDIR = "./figures"
	FIGSUFFIX = ".compare_profiles"
end

# ╔═╡ fbf40c95-37e6-4ae5-8f78-532fc3da4d77
import CSV

# ╔═╡ e1c0f486-2917-4c3f-a909-f0091ee27c58
CairoMakie.activate!(type=:png)

# ╔═╡ cec06b83-84de-4bb3-b41c-5dffcd6fe0f3
galaxyname = "bootes3"

# ╔═╡ ad174507-779f-4117-9d71-10843d42981d
function add_profile!(profiles, filename; label=filename)
	filepath = joinpath("../../observations", galaxyname, "density_profiles", filename)
	if isfile(filepath)
		@info "reading $filepath"
		prof = LilGuys.SurfaceDensityProfile(filepath)
		# idxs = LilGuys.find_longest_consecutive_true(LilGuys.sym_error.(prof.log_Sigma) .< 1.0)
		
		# if length(idxs) < length(prof.log_R)
		# 	prof =  LilGuys.filter_by_bin(prof, idxs)
		# end

		profiles[label] = prof
		return prof
	else
		@info "$filepath not found, skipping"
	end
end

# ╔═╡ cc9ac531-ef25-4573-9d35-4e4a25418adc
[prof.log_m_scale for (key, prof) in profiles]

# ╔═╡ 33a4d4a2-8f2b-4339-9e89-cf3919c56918
log_r_label = L"$\log\,R_\text{ell}$ /  arcmin"

# ╔═╡ e04789ef-2f38-4233-b2bb-426deebf451f
log_Sigma_label = L"$\log\,\Sigma_\star$ / arcmin$^{-2}$"

# ╔═╡ 161937ea-277c-49a8-8648-0793e52e7414
update_theme!(ErrorScatter=(; cycle=Cycle([:color=>:color, :marker=>:marker], covary=true)))

# ╔═╡ ea2beacd-d28a-4f7f-89a5-d462e4405b7b
md"""
# Plots
"""

# ╔═╡ b6baf3ed-9929-4252-b88d-9819d78be35f


# ╔═╡ 494cae7d-c78d-4680-b140-7441c0b8507b
markersize=5

# ╔═╡ 2300f4f3-bb27-4f16-9001-a9c626759a6b
function get_ylims(profiles)
	lower = Inf
	upper = -Inf
	for (key, prof) in profiles
		ys = filter(isfinite, prof.log_Sigma)
		y_l = middle.(ys) .- 1.1*abs.(lower_error.(ys))
		lower = min(minimum(y_l), lower)
		y_u = middle.(ys) .+ 1.1*upper_error.(ys)
		upper = max(maximum(y_u), upper)

	end

	return lower, upper

end

# ╔═╡ 21ee9c3c-6884-4240-8bd6-96a2bac74498
function plot_profiles(profiles_eq; jitter=0.005, markersize=markersize, legend=true, ylims=(nothing, nothing))
	fig = Figure(size=(5*72, 3.5*72))
	
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
		yticks = Makie.automatic,
		limits=(nothing, nothing, get_ylims(profiles_eq)...)
	)

	for (i, (key, prof)) in enumerate(profiles_eq)
		filt = isfinite.(prof.log_Sigma) 

		y =  LilGuys.middle.(prof.log_Sigma)
		ye = LilGuys.error_interval.(prof.log_Sigma)
		if i == 0
			offset = 0
		elseif i % 2 == 0
			offset = i/2 * jitter
		else
			offset = -(i-1)/2 * jitter
		end
		errorscatter!(prof.log_R[filt] .+ offset, y[filt], 
					  yerror=ye[filt], label=key, markersize=markersize
		 )	
	end
	if legend
		axislegend(position=:lb)
	end

	ylims!(ylims...)

	fig
end

# ╔═╡ 6ccf3382-3df8-4411-8350-4979ce49525a
mcmcdir = joinpath(ENV["DWARFS_ROOT"], "observations", "bootes3", "final_derivations/mcmc")

# ╔═╡ 50e23b20-bf23-4475-8f69-cfafc42b4024
samples_plummer = CSV.read("$mcmcdir/samples.j24_1c.mcmc_ell.csv", DataFrame)

# ╔═╡ 61dc2e5c-c02e-4838-b796-e968f81ce1d5
samples_sersic = CSV.read("$mcmcdir/samples.j24_1c.mcmc_sersic_ell.csv", DataFrame)

# ╔═╡ 2ad9c124-0343-4f14-98fb-1c5e753e09a1
function plot_plummer!(ax, df_out; color=COLORS[2])

	x = LinRange(0.5, 2.0, 1000)

	R = 10 .^ x

	for i in 1:100
		prof = LilGuys.Plummer(r_s=df_out.R_h[i], M=1)
		M =  df_out.N_sat[i]
		
		y = @. log10(LilGuys.surface_density(prof, R) * M)

		lines!(ax, x, y, color=color, alpha=0.1)
	end

end

# ╔═╡ abef5257-45f4-432e-9359-53a717cf6f38
function plot_sersic!(ax, df_out)


	x = LinRange(0.5, 2.0, 1000)

	R = 10 .^ x

	for i in 1:100
		prof = LilGuys.Sersic(R_h=df_out.R_h[i], n=df_out.n[i])
		M = df_out.N_sat[i] / LilGuys.mass_2D(prof, 600)
		
		y = @. log10(LilGuys.surface_density(prof, R) * M)

		lines!(ax, x, y, color=COLORS[3], alpha=0.1)
	end



end

# ╔═╡ e6ed390a-9da9-4d7b-b344-5f506411141d
function plot_exp!(ax, df_out)


	x = LinRange(0.5, 2.0, 1000)

	R = 10 .^ x

	for i in 1:100
		prof = LilGuys.Sersic(R_h=df_out.R_h[i], n=df_out.n[i])
		M = df_out.N_sat[i] / LilGuys.mass_2D(prof, 600)
		
		y = @. log10(LilGuys.surface_density(prof, R) * M)

		lines!(ax, x, y, color=COLORS[3], alpha=0.1)
	end



end

# ╔═╡ 866bbcf1-75e7-4556-b782-32ba92cb8035
let 
	profs = OrderedDict{String, LilGuys.SurfaceDensityProfile}()
	add_profile!(profs, "jax_LLR_1_eqw_profile.toml", label="selected stars")
	# add_profile!(profs, "jax_LLR_1_eqw_sub_profile.toml", label="LL>1")

	
	fig = plot_profiles(profs, legend=false)
	plot_plummer!(fig.content[1], samples_plummer)
	plot_sersic!(fig.content[1], samples_sersic )
	ylims!(-3, -1)
	lines!([NaN], [NaN], color=COLORS[2], label="Plummer samples")
	lines!([NaN], [NaN], color=COLORS[3], label="Sérsic samples")

	axislegend(position=:lb)
	@savefig "obs_density_profile"
	fig
end

# ╔═╡ 9b5e2dae-cc92-4778-aa77-d79fc1c91008
let
	profiles_mcmc = OrderedDict{String, LilGuys.SurfaceDensityProfile}()
	add_profile!(profiles_mcmc, "jax_LLR_1_eqw_profile.toml", label="fiducial")
	add_profile!(profiles_mcmc, "../mcmc/hist_fast_profile.toml", label="mcmc")

	plot_profiles(profiles_mcmc, ylims=(-5, 0))

end

# ╔═╡ Cell order:
# ╠═02612c68-042b-11f0-2f84-e9431cc0ee83
# ╠═fbf40c95-37e6-4ae5-8f78-532fc3da4d77
# ╠═cb95f753-b139-4150-bed9-75cd791ceffc
# ╠═e1c0f486-2917-4c3f-a909-f0091ee27c58
# ╠═abead192-161f-4375-899f-2c769d40cbfb
# ╠═cec06b83-84de-4bb3-b41c-5dffcd6fe0f3
# ╠═898d2c8b-b761-4a69-b561-658a644f44df
# ╠═ad174507-779f-4117-9d71-10843d42981d
# ╠═cc9ac531-ef25-4573-9d35-4e4a25418adc
# ╠═33a4d4a2-8f2b-4339-9e89-cf3919c56918
# ╠═e04789ef-2f38-4233-b2bb-426deebf451f
# ╠═161937ea-277c-49a8-8648-0793e52e7414
# ╟─ea2beacd-d28a-4f7f-89a5-d462e4405b7b
# ╠═b6baf3ed-9929-4252-b88d-9819d78be35f
# ╠═494cae7d-c78d-4680-b140-7441c0b8507b
# ╠═2300f4f3-bb27-4f16-9001-a9c626759a6b
# ╠═21ee9c3c-6884-4240-8bd6-96a2bac74498
# ╠═6ccf3382-3df8-4411-8350-4979ce49525a
# ╠═50e23b20-bf23-4475-8f69-cfafc42b4024
# ╠═61dc2e5c-c02e-4838-b796-e968f81ce1d5
# ╠═2ad9c124-0343-4f14-98fb-1c5e753e09a1
# ╠═abef5257-45f4-432e-9359-53a717cf6f38
# ╠═e6ed390a-9da9-4d7b-b344-5f506411141d
# ╠═866bbcf1-75e7-4556-b782-32ba92cb8035
# ╠═9b5e2dae-cc92-4778-aa77-d79fc1c91008
