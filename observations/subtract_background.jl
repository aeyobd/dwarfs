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

# ╔═╡ 4a6522fd-39ec-4cba-8494-b4a50a19d9fb
md"""
The goal of this notebook is to determine the background density of a given density profile and subtract this value. 
This is likely better done with a bayesian analysis
"""

# ╔═╡ 4a985c25-0077-4eae-b2f9-9e7d4caa3c28
import NaNMath as nm

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

# ╔═╡ 7557d9c9-891a-4bce-88d6-f04e06e81738
@bind profilename confirm(TextField(default="jax_LL_0_profile.toml"))

# ╔═╡ fe484bae-bd19-41c0-9d9d-7a6446b99039
@bind n_bins_bg confirm(NumberField(1:1:100))

# ╔═╡ cb844ce4-334c-45ce-8343-d38722ec185b
@bind n_bins_cut confirm(NumberField(0:1:100, default=1))

# ╔═╡ ed0c091b-d638-428e-800e-d0aacda2b01c
profpath = joinpath(galaxyname, "density_profiles", profilename)


# ╔═╡ 7b07b18f-8048-4947-b72c-ffaf69306a44
prof = LilGuys.StellarDensityProfile(profpath)

# ╔═╡ 858c7c17-8e6e-47f6-aa29-30573c4a2b8a
Nb = length(prof.log_R)

# ╔═╡ c2206a3a-f41e-4b23-93ad-5ec23874d943
filt_sat = 1:(Nb-n_bins_bg-n_bins_cut)

# ╔═╡ 22b1035d-43f5-48e8-ac9e-5bb5a954887e
filt_bg = (Nb-n_bins_bg-n_bins_cut+1):(Nb-n_bins_cut)

# ╔═╡ 5816a3be-20b9-447d-8efa-889e5603bed3
Sigmas_bg = prof.log_Sigma[end-n_bins_bg-n_bins_cut+1:end-n_bins_cut]

# ╔═╡ 3a35f811-11cf-4154-86c0-c016a6253631
Sigmas_bg_err = maximum.(LilGuys.ci_of.(Sigmas_bg))

# ╔═╡ 3cc7bd35-0f9f-48c9-aeb1-d29cc8381869
Sigmas_bg_m = LilGuys.value.(Sigmas_bg)

# ╔═╡ f56378dd-e0dd-4856-a192-d7c2b7fa0222
Sigma_bg = LilGuys.mean(Sigmas_bg_m, 1 ./ Sigmas_bg_err .^ 2)

# ╔═╡ 63ffa29a-9157-4af5-ba08-dc945e70240d
Sigma_bg_err = sqrt(1/LilGuys.sum(Sigmas_bg_err .^ -2))

# ╔═╡ 7f4bad59-ec23-4add-aba8-58cb6fcbf9e6
Sigma_bg_u = LilGuys.Measurement(Sigma_bg, Sigma_bg_err)

# ╔═╡ a0b6e938-be54-45e1-9433-3f68f2cbb8da
begin
	Sigma_sub = 10 .^ prof.log_Sigma .- 10 .^ Sigma_bg_u
	log_Sigma_sub = log10.(Sigma_sub[1:end-n_bins_bg-n_bins_cut])
	prof_sub = LilGuys.StellarDensityProfile(
		R_units=prof.R_units,
		log_R=prof.log_R[1:end-n_bins_bg-n_bins_cut],
		log_R_bins=prof.log_R_bins[1:end-n_bins_bg-n_bins_cut],
		log_Sigma=log_Sigma_sub,
		log_R_scale=prof.log_R_scale,
		log_m_scale=prof.log_m_scale,
		annotations=prof.annotations
	)
end

# ╔═╡ 33a4d4a2-8f2b-4339-9e89-cf3919c56918
log_r_label = L"$\log\,R$ /  arcmin"

# ╔═╡ e04789ef-2f38-4233-b2bb-426deebf451f
log_Sigma_label = L"$\log\,\Sigma$"

# ╔═╡ 8a1b3f9e-4a88-4527-b9ee-ebb97b3cfc8d
let
	fig = Figure(size=(4*72, 3*72))
	
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
		yticks = Makie.automatic,
	)

	errorscatter!(prof.log_R, LilGuys.value.(prof.log_Sigma), 
		yerror=LilGuys.ci_of.(prof.log_Sigma),
	)

	scatter!(prof.log_R[filt_bg], LilGuys.value.(prof.log_Sigma)[filt_bg], color=COLORS[3])
	fig
end

# ╔═╡ 21814a70-ce10-43f9-bb02-91bf15da860f
let
	fig = Figure(size=(4*72, 3*72))
	
	ax = Axis(fig[1,1],
		xlabel = log_r_label,
		ylabel = log_Sigma_label,
		yticks = Makie.automatic,
	)

	errorscatter!(prof.log_R, LilGuys.value.(prof.log_Sigma), 
		yerror=LilGuys.ci_of.(prof.log_Sigma),
	)

	scatter!(prof.log_R[filt_bg], LilGuys.value.(prof.log_Sigma)[filt_bg], color=COLORS[3])


	errorscatter!(prof_sub.log_R, LilGuys.value.(prof_sub.log_Sigma), 
		yerror=LilGuys.ci_of.(prof_sub.log_Sigma),
	)

	hlines!(Sigma_bg)
	hspan!(Sigma_bg - Sigma_bg_err, Sigma_bg + Sigma_bg_err, alpha=0.2)

	fig
end

# ╔═╡ c35fcdeb-1165-4fa6-882f-d2a24eff41fb
profilename_out = split(profilename, "_profile")[1] * "_sub_profile.toml"

# ╔═╡ b74aee7a-6f7b-43ca-8c30-baa68d679c6f
profilename_bg = split(profilename, "_profile")[1] * "_background.toml"

# ╔═╡ 3a866a51-d79b-4722-8243-0043192d9ce5
open(joinpath(galaxyname, "density_profiles", profilename_out), "w") do f
	print(f, prof_sub)
end

# ╔═╡ 3ecf11e6-1dab-476d-a012-d492cce17795
df_bg = OrderedDict(
	"log_Sigma" => Sigma_bg,
	"log_Sigma_err" => Sigma_bg_err,
	"log_R_bg" => prof.log_R[filt_bg]
)

# ╔═╡ a760ebfc-204c-449d-b83b-c4d2282d2f6e
open(joinpath(galaxyname, "density_profiles", profilename_bg), "w") do f
	print(f, df_bg)
end

# ╔═╡ Cell order:
# ╠═4a6522fd-39ec-4cba-8494-b4a50a19d9fb
# ╠═02612c68-042b-11f0-2f84-e9431cc0ee83
# ╠═4a985c25-0077-4eae-b2f9-9e7d4caa3c28
# ╠═e1c0f486-2917-4c3f-a909-f0091ee27c58
# ╠═abead192-161f-4375-899f-2c769d40cbfb
# ╠═cec06b83-84de-4bb3-b41c-5dffcd6fe0f3
# ╠═7557d9c9-891a-4bce-88d6-f04e06e81738
# ╠═fe484bae-bd19-41c0-9d9d-7a6446b99039
# ╠═cb844ce4-334c-45ce-8343-d38722ec185b
# ╠═898d2c8b-b761-4a69-b561-658a644f44df
# ╠═ed0c091b-d638-428e-800e-d0aacda2b01c
# ╠═7b07b18f-8048-4947-b72c-ffaf69306a44
# ╠═858c7c17-8e6e-47f6-aa29-30573c4a2b8a
# ╠═c2206a3a-f41e-4b23-93ad-5ec23874d943
# ╠═22b1035d-43f5-48e8-ac9e-5bb5a954887e
# ╠═5816a3be-20b9-447d-8efa-889e5603bed3
# ╠═3a35f811-11cf-4154-86c0-c016a6253631
# ╠═3cc7bd35-0f9f-48c9-aeb1-d29cc8381869
# ╠═f56378dd-e0dd-4856-a192-d7c2b7fa0222
# ╠═63ffa29a-9157-4af5-ba08-dc945e70240d
# ╠═7f4bad59-ec23-4add-aba8-58cb6fcbf9e6
# ╠═a0b6e938-be54-45e1-9433-3f68f2cbb8da
# ╠═33a4d4a2-8f2b-4339-9e89-cf3919c56918
# ╠═e04789ef-2f38-4233-b2bb-426deebf451f
# ╠═8a1b3f9e-4a88-4527-b9ee-ebb97b3cfc8d
# ╠═21814a70-ce10-43f9-bb02-91bf15da860f
# ╠═c35fcdeb-1165-4fa6-882f-d2a24eff41fb
# ╠═b74aee7a-6f7b-43ca-8c30-baa68d679c6f
# ╠═3a866a51-d79b-4722-8243-0043192d9ce5
# ╠═3ecf11e6-1dab-476d-a012-d492cce17795
# ╠═a760ebfc-204c-449d-b83b-c4d2282d2f6e
