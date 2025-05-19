### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 02612c68-042b-11f0-2f84-e9431cc0ee83
begin
	using Pkg; Pkg.activate()

	using PlutoUI
	using CairoMakie
	using Arya

end

# ╔═╡ abead192-161f-4375-899f-2c769d40cbfb
using OrderedCollections

# ╔═╡ e1c0f486-2917-4c3f-a909-f0091ee27c58
CairoMakie.activate!(type=:png)

# ╔═╡ cec06b83-84de-4bb3-b41c-5dffcd6fe0f3
@bind galaxyname confirm(TextField(default="draco"))

# ╔═╡ 898d2c8b-b761-4a69-b561-658a644f44df
begin
	using LilGuys
	FIGDIR = joinpath(galaxyname, "figures")
	FIGSUFFIX = ".compare_profiles"
end

# ╔═╡ ad174507-779f-4117-9d71-10843d42981d
function add_profile!(profiles, filename; label=filename)
	filepath = joinpath(galaxyname, "density_profiles", filename)
	if isfile(filepath)
		prof = LilGuys.StellarDensityProfile(filepath)
		if prof.log_m_scale != 0
			prof.log_Sigma .-= prof.log_m_scale
		end
		filt = maximum.(error_interval.(prof.log_Sigma)) .< 1
		prof.log_R = prof.log_R[filt]
		prof.log_Sigma = prof.log_Sigma[filt]
		if length(prof.Gamma) > 0
			prof.Gamma = prof.Gamma[filt]
		end
		profiles[label] = prof
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

# ╔═╡ 2c2b136d-f22e-460a-91c3-ce8ec858dbc7
md"""
This plot uses my fiducial binning method and compares a variety of different methods.
"""

# ╔═╡ b6baf3ed-9929-4252-b88d-9819d78be35f
jitter=0.005

# ╔═╡ 494cae7d-c78d-4680-b140-7441c0b8507b
markersize=3

# ╔═╡ 21ee9c3c-6884-4240-8bd6-96a2bac74498
function plot_profiles(profiles_eq; jitter=0.005, markersize=markersize)
	fig = Figure(size=(4*72, 3*72))
	
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
		yticks = Makie.automatic,
		limits=(nothing, nothing, get_ylims(profiles_eq)...)
	)

	for (key, prof) in profiles_eq
		filt = isfinite.(prof.log_Sigma) 

		y =  LilGuys.middle.(prof.log_Sigma)
		ye = LilGuys.error_interval.(prof.log_Sigma)
		errorscatter!(prof.log_R[filt] .+ jitter * LilGuys.randu(-1, 1), y[filt], 
					  yerror=ye[filt], label=key, size=markersize
		 )	
	end

	axislegend(position=:lb)
	key_ref = collect(keys(profiles_eq))[begin]

	ax2 = Axis(fig[2,1],
		xlabel = log_r_label,
		ylabel = L"\log\Sigma / \Sigma_\textrm{%$key_ref}",
		limits=(nothing, nothing, -1, 1)
	)
	hlines!(0, color=:black)

	prof_ref = profiles_eq[key_ref]
	for (key, prof) in profiles_eq
		ym = LilGuys.lerp(prof_ref.log_R, LilGuys.middle.(prof_ref.log_Sigma)).(prof.log_R)

		filt = isfinite.(prof.log_Sigma) .& isfinite.(ym)

		y =  LilGuys.middle.(prof.log_Sigma) .- ym
		ye = LilGuys.error_interval.(prof.log_Sigma)
		errorscatter!(prof.log_R[filt].+ jitter * LilGuys.randu(-1, 1), y[filt], 
					  yerror=ye[filt], size=markersize
					 )
	end

	hidexdecorations!(ax)

	linkxaxes!(ax, ax2)
	rowsize!(fig.layout, 2, Relative(1/4))
	rowgap!(fig.layout, 0.)

	fig
end

# ╔═╡ 5295c0da-f77f-45b9-84bf-e174e2d5a521
begin
	profiles = OrderedDict{String, LilGuys.StellarDensityProfile}()
	add_profile!(profiles, "jax_profile.toml", label="fiducial")
	add_profile!(profiles, "jax_2c_profile.toml", label="2exp")
	add_profile!(profiles, "jax_LLR_0_sub_profile.toml", label="LL cut")
	add_profile!(profiles, "simple_sub_profile.toml", label="simple")
	add_profile!(profiles, "../mcmc/hist_struct_profile.toml", label="MCMC struct")
	add_profile!(profiles, "../mcmc/hist_fast_profile.toml", label="MCMC")

	# add_profile!(profiles, "jax_LL_0_profile.toml", label="LLR cut")
	add_profile!(profiles, "jax_circ_eqw_profile.toml", label="circ.")
	#add_profile!(profiles, "../processed/profile.mcmc_hist.toml", label="piecewise w/ ell&PA")

end

# ╔═╡ cc9ac531-ef25-4573-9d35-4e4a25418adc
[prof.log_m_scale for (key, prof) in profiles]

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
			LilGuys.gammaplot!(ax, prof, label=key, size=3, cycle=[:color])
		end
	end

	axislegend(position=:lb)

	@savefig "Gamma_by_assumptions"
	fig
end

# ╔═╡ 035abb97-067a-4b3b-abf3-d2713a2da8b6
[(key => prof.Gamma) for (key, prof) in profiles]

# ╔═╡ d2e1374d-978d-41b9-99fb-e1d025dc9bbc
get_ylims(profiles)

# ╔═╡ 27d95410-65dd-48b1-a9eb-e67f3e91cf18
@savefig "log_Sigma" plot_profiles(profiles, jitter=jitter)

# ╔═╡ 27164128-5562-48e0-ab4b-e13965287946
begin
	profiles_eq = OrderedDict{String, LilGuys.StellarDensityProfile}()
	add_profile!(profiles_eq, "jax_2c_eqw_profile.toml", label="2exp")
	add_profile!(profiles_eq, "jax_eqw_profile.toml", label="1c")
	add_profile!(profiles_eq, "jax_LLR_0_eqw_profile.toml", label="LL")
	add_profile!(profiles_eq, "simple_eqw_sub_profile.toml", label="simple")
	add_profile!(profiles_eq, "best_eqw_profile.toml", label="best")
	add_profile!(profiles_eq, "best_eqw_sub_profile.toml", label="best - bg")


	# add_profile!(profiles, "jax_LL_0_profile.toml", label="LLR cut")
	# add_profile!(profiles, "jax_circ_profile.toml", label="2exp circ.")
	#add_profile!(profiles, "../processed/profile.mcmc_hist.toml", label="piecewise w/ ell&PA")

end

# ╔═╡ 99732187-a038-424e-b31d-15aff4a2afaa
profiles_eq["LL"]

# ╔═╡ e6373852-3932-4d6f-9d1b-cde32a1f73f5
@savefig "log_Sigma_eqw" plot_profiles(profiles_eq)

# ╔═╡ f8924d9d-6a88-4180-8222-3714ee4e0e37
begin
	profiles_methods = OrderedDict{String, LilGuys.StellarDensityProfile}()
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
@savefig "log_Sigma_methods" plot_profiles(profiles_methods)

# ╔═╡ 9b5e2dae-cc92-4778-aa77-d79fc1c91008
begin
	profiles_mcmc = OrderedDict{String, LilGuys.StellarDensityProfile}()
	add_profile!(profiles_mcmc, "jax_eqw_profile.toml", label="fiducial")
	# add_profile!(profiles_methods, "jax_LLR_0_eqw_sub_profile.toml", label="LL")
	add_profile!(profiles_mcmc, "../mcmc/hist_struct_profile.toml", label="hist struct")
	add_profile!(profiles_mcmc, "../mcmc/hist_fast_profile.toml", label="hist")
	add_profile!(profiles_mcmc, "../mcmc/hist_fast_approx_profile.toml", label="hist approx")
	add_profile!(profiles_mcmc, "../mcmc/hist_fast_nostruct_profile.toml", label="hist nostruct")

	#add_profile!(profiles, "../processed/profile.mcmc_hist.toml", label="piecewise w/ ell&PA")

end

# ╔═╡ e7c839df-97b2-45f6-a577-4a9ba8d2b39d
@savefig "log_Sigma_mcmc" plot_profiles(profiles_mcmc)

# ╔═╡ Cell order:
# ╠═02612c68-042b-11f0-2f84-e9431cc0ee83
# ╠═e1c0f486-2917-4c3f-a909-f0091ee27c58
# ╠═abead192-161f-4375-899f-2c769d40cbfb
# ╠═cec06b83-84de-4bb3-b41c-5dffcd6fe0f3
# ╠═898d2c8b-b761-4a69-b561-658a644f44df
# ╠═ad174507-779f-4117-9d71-10843d42981d
# ╠═cc9ac531-ef25-4573-9d35-4e4a25418adc
# ╠═33a4d4a2-8f2b-4339-9e89-cf3919c56918
# ╠═e04789ef-2f38-4233-b2bb-426deebf451f
# ╠═161937ea-277c-49a8-8648-0793e52e7414
# ╠═2f1a42a3-e100-4413-a039-a20244ada256
# ╠═035abb97-067a-4b3b-abf3-d2713a2da8b6
# ╠═99732187-a038-424e-b31d-15aff4a2afaa
# ╠═2300f4f3-bb27-4f16-9001-a9c626759a6b
# ╠═d2e1374d-978d-41b9-99fb-e1d025dc9bbc
# ╠═21ee9c3c-6884-4240-8bd6-96a2bac74498
# ╠═2c2b136d-f22e-460a-91c3-ce8ec858dbc7
# ╠═b6baf3ed-9929-4252-b88d-9819d78be35f
# ╠═494cae7d-c78d-4680-b140-7441c0b8507b
# ╠═5295c0da-f77f-45b9-84bf-e174e2d5a521
# ╠═27d95410-65dd-48b1-a9eb-e67f3e91cf18
# ╠═27164128-5562-48e0-ab4b-e13965287946
# ╠═e6373852-3932-4d6f-9d1b-cde32a1f73f5
# ╠═f8924d9d-6a88-4180-8222-3714ee4e0e37
# ╠═e835372a-3ef6-4ede-82b0-0f2d24c798c5
# ╠═9b5e2dae-cc92-4778-aa77-d79fc1c91008
# ╠═e7c839df-97b2-45f6-a577-4a9ba8d2b39d
