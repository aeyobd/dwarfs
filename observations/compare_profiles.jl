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
@bind galaxyname confirm(TextField(default="sculptor"))

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
		
		profiles[label] = prof
	else
		@info "$filepath not found, skipping"
	end
end

# ╔═╡ 5295c0da-f77f-45b9-84bf-e174e2d5a521
begin
	profiles = OrderedDict{String, LilGuys.StellarDensityProfile}()
	add_profile!(profiles, "jax_2c_profile.toml", label="2exp")
	add_profile!(profiles, "jax_profile.toml", label="1c")
	add_profile!(profiles, "jax_LLR_0_sub_profile.toml", label="LL")
	add_profile!(profiles, "simple_sub_profile.toml", label="simple")
	# add_profile!(profiles, "jax_LL_0_profile.toml", label="LLR cut")
	# add_profile!(profiles, "jax_circ_profile.toml", label="2exp circ.")
	# add_profile!(profiles, "../processed/profile.mcmc_hist_nostruct.toml", label="piecewise")
	#add_profile!(profiles, "../processed/profile.mcmc_hist.toml", label="piecewise w/ ell&PA")

end

# ╔═╡ cc9ac531-ef25-4573-9d35-4e4a25418adc
[prof.log_m_scale for (key, prof) in profiles]

# ╔═╡ 33a4d4a2-8f2b-4339-9e89-cf3919c56918
log_r_label = L"$\log\,R$ /  arcmin"

# ╔═╡ e04789ef-2f38-4233-b2bb-426deebf451f
log_Sigma_label = L"$\log\,\Sigma$"

# ╔═╡ 161937ea-277c-49a8-8648-0793e52e7414
update_theme!(ErrorScatter=(; cycle=Cycle([:color=>:color, :marker=>:marker], covary=true)))

# ╔═╡ 21814a70-ce10-43f9-bb02-91bf15da860f
let
	fig = Figure(size=(4*72, 3*72))
	
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
		yticks = Makie.automatic,
	)

	for (key, prof) in profiles
		LilGuys.plot_log_Σ!(ax, prof, label=key)
	end

	axislegend(position=:lb)

	ax2 = Axis(fig[2,1],
		xlabel = log_r_label,
		ylabel = L"\log\Sigma / \Sigma_\textrm{%$(collect(keys(profiles))[begin])}",
		limits=(nothing, nothing, -1, 1)
	)

	prof_ref = profiles[collect(keys(profiles))[begin]]
	for (key, prof) in profiles
		ym = LilGuys.lerp(prof_ref.log_R, LilGuys.middle.(prof_ref.log_Sigma)).(prof.log_R)
		
		errorscatter!(prof.log_R, LilGuys.middle.(prof.log_Sigma) .- ym, 
					  yerror=LilGuys.error_interval.(prof.log_Sigma)
					 )
	end

	hidexdecorations!(ax)

	linkxaxes!(ax, ax2)
	rowsize!(fig.layout, 2, Relative(1/4))
	rowgap!(fig.layout, 0.)

	@savefig "log_Sigma_by_assumptions"
	fig
end

# ╔═╡ 2f1a42a3-e100-4413-a039-a20244ada256
let
	fig = Figure()
	
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
		yticks = Makie.automatic,
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

# ╔═╡ 21ee9c3c-6884-4240-8bd6-96a2bac74498
LilGuys.projecteddensity

# ╔═╡ Cell order:
# ╠═02612c68-042b-11f0-2f84-e9431cc0ee83
# ╠═e1c0f486-2917-4c3f-a909-f0091ee27c58
# ╠═abead192-161f-4375-899f-2c769d40cbfb
# ╠═cec06b83-84de-4bb3-b41c-5dffcd6fe0f3
# ╠═898d2c8b-b761-4a69-b561-658a644f44df
# ╠═ad174507-779f-4117-9d71-10843d42981d
# ╠═5295c0da-f77f-45b9-84bf-e174e2d5a521
# ╠═cc9ac531-ef25-4573-9d35-4e4a25418adc
# ╠═33a4d4a2-8f2b-4339-9e89-cf3919c56918
# ╠═e04789ef-2f38-4233-b2bb-426deebf451f
# ╠═161937ea-277c-49a8-8648-0793e52e7414
# ╠═21814a70-ce10-43f9-bb02-91bf15da860f
# ╠═2f1a42a3-e100-4413-a039-a20244ada256
# ╠═035abb97-067a-4b3b-abf3-d2713a2da8b6
# ╠═21ee9c3c-6884-4240-8bd6-96a2bac74498
