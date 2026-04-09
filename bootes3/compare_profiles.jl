### A Pluto.jl notebook ###
# v0.20.23

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
	filepath = joinpath("../observations", galaxyname, "density_profiles", filename)
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

# ╔═╡ 33a4d4a2-8f2b-4339-9e89-cf3919c56918
log_r_label = L"$\log\,R$ /  arcmin"

# ╔═╡ e04789ef-2f38-4233-b2bb-426deebf451f
log_Sigma_label = L"$\log\,\Sigma$"

# ╔═╡ 161937ea-277c-49a8-8648-0793e52e7414
update_theme!(ErrorScatter=(; cycle=Cycle([:color=>:color, :marker=>:marker], covary=true)))

# ╔═╡ ea2beacd-d28a-4f7f-89a5-d462e4405b7b
md"""
# Plots
"""

# ╔═╡ b6baf3ed-9929-4252-b88d-9819d78be35f


# ╔═╡ 494cae7d-c78d-4680-b140-7441c0b8507b
markersize=5

# ╔═╡ 6ccf3382-3df8-4411-8350-4979ce49525a
mcmcdir = joinpath(ENV["DWARFS_ROOT"], "observations", "bootes3", "mcmc")

# ╔═╡ 50e23b20-bf23-4475-8f69-cfafc42b4024
samples_plummer = CSV.read("$mcmcdir/samples.best_sample.mcmc_ell.csv", DataFrame)

# ╔═╡ 61dc2e5c-c02e-4838-b796-e968f81ce1d5
samples_sersic = CSV.read("$mcmcdir/samples.best_sample.mcmc_sersic_ell.csv", DataFrame)

# ╔═╡ ad5475f8-c80a-4bc8-8aee-636b584a4904
samples_sersic_delve = CSV.read("$mcmcdir/samples.delve_matched_filter_6deg.mcmc_sersic.csv", DataFrame)

# ╔═╡ 9eeb6dd0-18ef-48d2-af40-3b60e87d874f
samples_plummer_delve = CSV.read("$mcmcdir/samples.delve_matched_filter_6deg.mcmc_plummer.csv", DataFrame)

# ╔═╡ 04168b26-6e36-41da-9162-45440554a9c9
samples_plummer_delve.N_sat ./ samples_plummer_delve.f_sat

# ╔═╡ 2ad9c124-0343-4f14-98fb-1c5e753e09a1
function plot_plummer!(ax, df_out, N_stars; color=COLORS[2])

	x = LinRange(0.5, 2.0, 1000)

	R = 10 .^ x

	for i in 1:100
		prof = LilGuys.Plummer(r_s=df_out.R_h[i], M=1)
		M = N_stars * df_out.f_sat[i]
		
		y = @. log10(LilGuys.surface_density(prof, R) * M)

		lines!(ax, x, y, color=color, alpha=0.1)
	end

end

# ╔═╡ abef5257-45f4-432e-9359-53a717cf6f38
function plot_sersic!(ax, df_out, N_stars)


	x = LinRange(0.5, 2.0, 1000)

	R = 10 .^ x

	for i in 1:100
		prof = LilGuys.Sersic(R_h=df_out.R_h[i], n=df_out.n[i])
		M = N_stars * df_out.f_sat[i] / LilGuys.mass_2D(prof, 600)
		
		y = @. log10(LilGuys.surface_density(prof, R) * M)

		lines!(ax, x, y, color=COLORS[3], alpha=0.1)
	end



end

# ╔═╡ 978487c1-ca31-4840-b7e2-a8253fd29d51
N_gaia = 19440.0

# ╔═╡ bae678d8-a093-4dd7-b061-cb2f671a09c8
N_delve = 352678.0 ./ 3

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
function plot_profiles(profiles_eq; jitter=0.005, markersize=markersize, legend=true)
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
	key_ref = collect(keys(profiles_eq))[begin]

	# ax2 = Axis(fig[2,1],
	# 	xlabel = log_r_label,
	# 	ylabel = L"\log\Sigma / \Sigma_\textrm{%$key_ref}",
	# 	limits=(nothing, nothing, -1, 1)
	# )
	# hlines!(0, color=:black)

	# prof_ref = profiles_eq[key_ref]
	# for (key, prof) in profiles_eq
	# 	ym = LilGuys.lerp(prof_ref.log_R, LilGuys.middle.(prof_ref.log_Sigma)).(prof.log_R)

	# 	filt = isfinite.(prof.log_Sigma) .& isfinite.(ym)

	# 	y =  LilGuys.middle.(prof.log_Sigma) .- ym
	# 	ye = LilGuys.error_interval.(prof.log_Sigma)
	# 	errorscatter!(prof.log_R[filt].+ jitter * LilGuys.randu(-1, 1), y[filt], 
	# 				  yerror=ye[filt], markersize=markersize
	# 				 )
	# end

	# hidexdecorations!(ax)

	# linkxaxes!(ax, ax2)
	# rowsize!(fig.layout, 2, Relative(1/4))
	# rowgap!(fig.layout, 0.)

	fig
end

# ╔═╡ 866bbcf1-75e7-4556-b782-32ba92cb8035
let 
	profs = OrderedDict{String, LilGuys.SurfaceDensityProfile}()
	add_profile!(profs, "jax_LLR_1_eqw_profile.toml", label="LL>1 w bg")
	# add_profile!(profs, "jax_LLR_0_eqw_sub_profile.toml", label="LL>0")
	add_profile!(profs, "jax_LLR_1_eqw_sub_profile.toml", label="LL>1")
	# add_profile!(profs, "jax_LLR_2_eqw_sub_profile.toml", label="LL>2")
	# add_profile!(profs, "delve_mf_1_eqw_sub_profile.toml", label="delve")
	fig = plot_profiles(profs)


	plot_plummer!(fig.content[1], samples_plummer, N_gaia)
	plot_sersic!(fig.content[1], samples_sersic, N_gaia)
	ylims!(-3.5, -1)
	@savefig "LLR_density_comparison"
	fig
end

# ╔═╡ 53f920f6-5bbd-4f70-9730-57dc28aa2c40
let 
	profs = OrderedDict{String, LilGuys.SurfaceDensityProfile}()
	add_profile!(profs, "mf_0_sub_profile.toml", label="LL>0")
	add_profile!(profs, "mf_0.5_sub_profile.toml", label="LL>0.5")

	add_profile!(profs, "mf_1_sub_profile.toml", label="LL>1")

	f = plot_profiles(profs)

	plot_plummer!(f.content[1], samples_sersic, N_gaia .* 6, color=COLORS[4])
	plot_plummer!(f.content[1], samples_plummer_delve, N_delve .* 4 )
	plot_sersic!(f.content[1], samples_sersic_delve, N_delve .* 4)
	f

end

# ╔═╡ fb9d2db8-012c-4951-b58e-a9ba9fa80078
let 
	profs = OrderedDict{String, LilGuys.SurfaceDensityProfile}()
	add_profile!(profs, "mf_1_profile.toml", label="all")
	add_profile!(profs, "mf_bhb_-2_profile.toml", label="BHB")
	add_profile!(profs, "mf_rgb_1_profile.toml", label="RGB")
	add_profile!(profs, "mf_ms_1_profile.toml", label="MS")
	plot_profiles(profs, jitter=0.003)
end

# ╔═╡ 2c2b136d-f22e-460a-91c3-ce8ec858dbc7
md"""
This plot uses my fiducial binning method and compares a variety of different methods.
"""

# ╔═╡ 5295c0da-f77f-45b9-84bf-e174e2d5a521
begin
	profiles = OrderedDict{String, LilGuys.SurfaceDensityProfile}()
	add_profile!(profiles, "fiducial_profile.toml", label="J+24")
	add_profile!(profiles, "jax_LLR_1_eqw_sub_profile.toml", label="LL cut")
	add_profile!(profiles, "delve_mf_1_eqw_sub_profile.toml", label="delve")
	add_profile!(profiles, "simple_sub_profile.toml", label="simple")
	add_profile!(profiles, "../mcmc/hist_fast_profile.toml", label="MCMC")


end

# ╔═╡ cc9ac531-ef25-4573-9d35-4e4a25418adc
[prof.log_m_scale for (key, prof) in profiles]

# ╔═╡ 035abb97-067a-4b3b-abf3-d2713a2da8b6
[(key => prof.Gamma) for (key, prof) in profiles]

# ╔═╡ 1fcec5cc-fa17-40f4-bd5d-3dbff5db50dc
function tot_num_stars(prof)
	areas = π * diff(LilGuys.radii_bins(prof) .^ 2)
	counts = middle.(LilGuys.surface_density(prof)) .* areas
	return sum(counts[isfinite.(counts)])
end

# ╔═╡ 4a6e21c8-8deb-497c-ba3a-fa90869d1b6d
for (label, prof) in profiles
	@info label, tot_num_stars(prof)
end

# ╔═╡ 27d95410-65dd-48b1-a9eb-e67f3e91cf18
@savefig "log_Sigma" let
	fig = plot_profiles(profiles, legend=false)
	x = LinRange(0.5, 2.5, 1000)
	y = -2x
	lines!(x, y .+ 0.7, label=L"\Sigma\sim R^{-2}", color=:black)
	ylims!(-5, -1)
	xlims!(0.6, 2.5)
	axislegend(position=:rt)
	fig
end

# ╔═╡ f8924d9d-6a88-4180-8222-3714ee4e0e37
begin
	profiles_methods = OrderedDict{String, LilGuys.SurfaceDensityProfile}()
	add_profile!(profiles_methods, "jax_eqw_profile.toml", label="fiducial")
	add_profile!(profiles_methods, "jax_2c_eqw_profile.toml", label="2c")
	add_profile!(profiles_methods, "jax_circ_eqw_profile.toml", label="circ")
	# add_profile!(profiles_methods, "jax_LLR_0_eqw_sub_profile.toml", label="LL")
	add_profile!(profiles_methods, "best_eqw_sub_profile.toml", label="BG sub")
	#add_profile!(profiles_methods, "../mcmc/hist_struct_profile.toml", label="hist struct")
	add_profile!(profiles_methods, "../mcmc/hist_fast_profile.toml", label="mcmc")
	#add_profile!(profiles_methods, "../mcmc/hist_fast_approx_profile.toml", label="hist approx")

	# add_profile!(profiles, "jax_LL_0_profile.toml", label="LLR cut")
	#add_profile!(profiles, "../processed/profile.mcmc_hist.toml", label="piecewise w/ ell&PA")

end

# ╔═╡ e835372a-3ef6-4ede-82b0-0f2d24c798c5
plot_profiles(profiles_methods)

# ╔═╡ 9b5e2dae-cc92-4778-aa77-d79fc1c91008
let
	profiles_mcmc = OrderedDict{String, LilGuys.SurfaceDensityProfile}()
	add_profile!(profiles_mcmc, "fiducial_profile.toml", label="fiducial")
	add_profile!(profiles_mcmc, "../mcmc/hist_fast_profile.toml", label="hist")

	plot_profiles(profiles_mcmc)

end

# ╔═╡ 5e72529d-4ac2-49e6-b881-67f5184f7594
md"""
# Density profile slope
"""

# ╔═╡ 2f1a42a3-e100-4413-a039-a20244ada256
let
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = L"\Gamma",
		yticks = Makie.automatic,
		limits=(nothing, nothing, -12, 4)
	)

	for (key, prof) in profiles
		if length(prof.Gamma) > 0
			lines!(radii(prof), middle.(prof.Gamma), label=key)
		end
	end

	axislegend(position=:lb)

	@savefig "Gamma_by_assumptions"
	fig
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
# ╠═21ee9c3c-6884-4240-8bd6-96a2bac74498
# ╠═6ccf3382-3df8-4411-8350-4979ce49525a
# ╠═50e23b20-bf23-4475-8f69-cfafc42b4024
# ╠═61dc2e5c-c02e-4838-b796-e968f81ce1d5
# ╠═ad5475f8-c80a-4bc8-8aee-636b584a4904
# ╠═9eeb6dd0-18ef-48d2-af40-3b60e87d874f
# ╠═04168b26-6e36-41da-9162-45440554a9c9
# ╠═2ad9c124-0343-4f14-98fb-1c5e753e09a1
# ╠═abef5257-45f4-432e-9359-53a717cf6f38
# ╠═978487c1-ca31-4840-b7e2-a8253fd29d51
# ╠═bae678d8-a093-4dd7-b061-cb2f671a09c8
# ╠═866bbcf1-75e7-4556-b782-32ba92cb8035
# ╠═53f920f6-5bbd-4f70-9730-57dc28aa2c40
# ╠═fb9d2db8-012c-4951-b58e-a9ba9fa80078
# ╠═035abb97-067a-4b3b-abf3-d2713a2da8b6
# ╠═2300f4f3-bb27-4f16-9001-a9c626759a6b
# ╠═2c2b136d-f22e-460a-91c3-ce8ec858dbc7
# ╠═5295c0da-f77f-45b9-84bf-e174e2d5a521
# ╠═1fcec5cc-fa17-40f4-bd5d-3dbff5db50dc
# ╠═4a6e21c8-8deb-497c-ba3a-fa90869d1b6d
# ╠═27d95410-65dd-48b1-a9eb-e67f3e91cf18
# ╠═f8924d9d-6a88-4180-8222-3714ee4e0e37
# ╠═e835372a-3ef6-4ede-82b0-0f2d24c798c5
# ╠═9b5e2dae-cc92-4778-aa77-d79fc1c91008
# ╠═5e72529d-4ac2-49e6-b881-67f5184f7594
# ╠═2f1a42a3-e100-4413-a039-a20244ada256
