### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ dc9f8e5a-289d-11f1-b227-375ff5c871ba
let
	import Pkg;
	Pkg.activate()

	using LilGuys
	using Arya, CairoMakie
	import TOML
end

# ╔═╡ db1ea197-fa80-4049-9d8e-966afb448c20
using OrderedCollections

# ╔═╡ aea3559d-ab36-4499-81ec-7533994d0b39
CairoMakie.activate!(type=:png)

# ╔═╡ df11aa76-0bc5-4728-97dd-da074c518e04
obs_dir = joinpath(ENV["DWARFS_ROOT"], "observations/bootes3")

# ╔═╡ 13dbbaff-ae35-4817-9f53-3a50a4c44205
function add_profile!(profiles, filename; label=filename)
	filepath = joinpath(obs_dir, "density_profiles", filename)
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

# ╔═╡ 27084f67-b307-4eda-944a-63dc8a7df03c
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

# ╔═╡ 9e290938-5d7a-4f31-9398-5e4b771b274e
begin
	profiles_jax = OrderedDict{String, LilGuys.SurfaceDensityProfile}()
	add_profile!(profiles_jax, "jax_LLR_1_eqw_profile.toml", label="fiducial")
	# add_profile!(profiles_jax, "jax_2c_eqw_profile.toml", label="J+24")
	add_profile!(profiles_jax, "jax_eqw_profile.toml", label="J+24 1c")

	add_profile!(profiles_jax, "../mcmc/hist_fast_profile.toml", label="MCMC")

	# add_profile!(profiles_jax, "jax_circ_eqw_profile.toml", label="circ.")
	#add_profile!(profiles, "../processed/profile.mcmc_hist.toml", label="piecewise w/ ell&PA")
	add_profile!(profiles_jax, "delve_mf_sub_profile.toml", label="DELVE MF")

	profiles_jax["DELVE MF"].log_Sigma .-= 1.8

	profiles_jax

end

# ╔═╡ 00542773-26bf-4145-83fe-be5f094a06dc
markersize=3

# ╔═╡ e855a4de-c484-435c-9ed1-f628887f1e7b
log_r_label = L"$\log\,R$ /  arcmin"

# ╔═╡ 87a78e46-8c9f-4518-a003-757a2fe0c838
log_Sigma_label = L"$\log\,\Sigma$"

# ╔═╡ e02f4269-f6b2-47bc-aac9-8c0d9d60012c
function plot_profiles(profiles_eq; jitter=0.005, markersize=markersize, ylims=nothing)
	fig = Figure(figsize=(6*72, 4*72))

	if isnothing(ylims)
		ylims = get_ylims(profiles_eq)
	end
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
		yticks = Makie.automatic,
		limits=(nothing, ylims)
	)

	for (key, prof) in profiles_eq
		filt = isfinite.(prof.log_Sigma) 

		y =  LilGuys.middle.(prof.log_Sigma)
		ye = LilGuys.error_interval.(prof.log_Sigma)
		errorscatter!(prof.log_R[filt] .+ jitter * LilGuys.randu(-1, 1), y[filt], 
					  yerror=ye[filt], label=key, markersize=markersize
		 )	
	end

	axislegend(position=:lb)
	# key_ref = collect(keys(profiles_eq))[begin]

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

# ╔═╡ 6af4557f-5368-4c54-994e-83eba56b0b0b
update_theme!(ErrorScatter=(; cycle=Cycle([:color=>:color, :marker=>:marker], covary=true)))

# ╔═╡ b353f864-cc88-4f22-af1e-84f903ff444f
begin
	profiles = OrderedDict{String, LilGuys.SurfaceDensityProfile}()
	# add_profile!(profiles, "fiducial_profile.toml", label="J+24")
	# add_profile!(profiles, "jax_LLR_0_eqw_profile.toml", label="LLR0")
	add_profile!(profiles, "jax_LLR_0_eqw_sub_profile.toml", label="LLR0-bkg")
	add_profile!(profiles, "delve_mf_sub_profile.toml", label="DELVE MF")
	profiles["DELVE MF"].log_Sigma .-= 1.8

	profiles
end

# ╔═╡ 1c140656-0fa9-435c-ba22-8237e9eecb11
plot_profiles(profiles)

# ╔═╡ 4235c9be-87d3-4ffa-a19d-5db97bdefda3
let 
	f = plot_profiles(profiles_jax, ylims=(-6, -1))
	xlims!(0.7, 2.4)

	f
end

# ╔═╡ Cell order:
# ╠═dc9f8e5a-289d-11f1-b227-375ff5c871ba
# ╠═aea3559d-ab36-4499-81ec-7533994d0b39
# ╠═db1ea197-fa80-4049-9d8e-966afb448c20
# ╠═df11aa76-0bc5-4728-97dd-da074c518e04
# ╠═13dbbaff-ae35-4817-9f53-3a50a4c44205
# ╠═27084f67-b307-4eda-944a-63dc8a7df03c
# ╠═9e290938-5d7a-4f31-9398-5e4b771b274e
# ╠═e02f4269-f6b2-47bc-aac9-8c0d9d60012c
# ╠═00542773-26bf-4145-83fe-be5f094a06dc
# ╠═e855a4de-c484-435c-9ed1-f628887f1e7b
# ╠═87a78e46-8c9f-4518-a003-757a2fe0c838
# ╠═6af4557f-5368-4c54-994e-83eba56b0b0b
# ╠═b353f864-cc88-4f22-af1e-84f903ff444f
# ╠═1c140656-0fa9-435c-ba22-8237e9eecb11
# ╠═4235c9be-87d3-4ffa-a19d-5db97bdefda3
