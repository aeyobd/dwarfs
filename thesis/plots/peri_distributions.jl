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

# ╔═╡ f0e19873-1835-4b47-9389-f4034ec4b032
using OrderedCollections

# ╔═╡ 665d33fd-1d6f-40b8-b03b-b56a604d9862
using DataFrames

# ╔═╡ e11d4bde-923f-4e50-b7b0-ed4eb9ad2ae4
using Printf

# ╔═╡ 7b8afb33-83c0-444d-8333-ecd1f96c3e1b
using PyFITS

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 1368d513-2631-4adb-bf3a-cbb33cd93b34
import TOML

# ╔═╡ ad86ba28-77a0-42ed-ac09-7f0d295f4ebf
import Agama

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png)

# ╔═╡ 2dfea445-50ff-4983-8c33-44e2d05b0cdc
function calc_r_J(halo, rho_peri)
    return 10 ^ LilGuys.find_zero(lr -> LilGuys.mean_density(halo, 10^lr) - 3*rho_peri, log10(halo.r_s))
end

# ╔═╡ 4fa535d8-0098-439b-b0c1-6ec31c110d13
pot_ep = Agama.Potential(file = joinpath(ENV["DWARFS_ROOT"], "agama/potentials/EP2020.ini"))

# ╔═╡ a25c4187-8409-4efb-9884-10ba9b4af3e9
calc_r_J(NFW(r_circ_max=4.0, v_circ_max=38/V2KMS), 0.00030350705117854145)

# ╔═╡ 0744d9a2-256f-45ed-8ba5-36f549a08c2b
function print_quantity(label, r)
    @printf "%s kpc   \t%0.2f ± %0.2f\n" label LilGuys.mean(r) LilGuys.std(r)
end

# ╔═╡ 2b03efa2-f91a-4e50-a90d-8c31a929ce1e
function print_radius(label, r, df)
    r_arcmin = LilGuys.kpc2arcmin.(r, df.distance)
    @printf "%s kpc   \t%0.2f ± %0.2f\n" label LilGuys.mean(r) LilGuys.std(r)
    @printf "%s arcmin\t%0.2f ± %0.2f\n" label LilGuys.mean(r_arcmin) LilGuys.std(r_arcmin)
end

# ╔═╡ d03e54bb-6311-40d3-a669-b72d017b0901
function mean_density(pot, r, units)
    M =  Agama.enclosed_mass(pot, r, units) 
    V = (4π/3 * r .^ 3)

    return M ./ V
end

# ╔═╡ d537bae8-8eee-4fa9-b22b-7a933d0a40a6
log10(Measurement(27,6, 7))

# ╔═╡ 6fdf86dc-0130-4c6e-a8e9-923be6315afb
function sample_scl_halo()
	σv = 0
	local h
	for i in 1:100
	  	vmax = 31 + 3*randn()
		rmax = 6 * 10 ^ (0.15 * randn())
		h = LilGuys.NFW(v_circ_max=vmax / V2KMS, r_circ_max=rmax)
		σv = LilGuys.σv_star_mean(h, LilGuys.Exp2D(R_s=0.13))
		if σv * V2KMS > 9.7 - 0.2*3# - 1
			return h
		end

	end
	@warn "iterma"

	return h

end

# ╔═╡ fdd9469f-140d-4118-af2c-2c389c7c8f34
function sample_umi_halo()
	σv = 0
	local h
	for i in 1:100
		vmax = 27 * 10 .^ (0.1 * randn())
		rmax = 5 * 10 .^ (0.12 * randn())	
		
		h = LilGuys.NFW(v_circ_max=vmax / V2KMS, r_circ_max=rmax)
		σv = LilGuys.σv_star_mean(h, LilGuys.Exp2D(R_s=0.13))
		if σv * V2KMS > 8.6 #- 0.3*3# - 1
			return h
		end

	end
	@warn "itermax"

	return h

end

# ╔═╡ 8d1daca0-5d79-47c2-9cca-1974b0c3ec65
function sample_scl_halos(N)
	return [sample_scl_halo() for _ in 1:N]
end


# ╔═╡ d8533d10-63f6-4f9e-a7e2-49809d19f6f2
function sample_umi_halos(N)
	return [sample_umi_halo() for _ in 1:N]

end


# ╔═╡ a9a983f5-4570-40d3-9d6e-e22d085b332a
replace_missings(x) = ifelse.(ismissing.(x), NaN, x)

# ╔═╡ 72736bff-eccc-491e-8c31-63cc52ea95ff
function analyze_lmc(galaxyname, modelname, halo; units=Agama.VASILIEV_UNITS)

    obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))


	modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, modelname)

    df = read_fits(joinpath(modeldir, "orbital_properties.fits"))

    pot = Agama.Potential(file = joinpath(modeldir, "potential_mw_init.ini"))
    pot_lmc = Agama.Potential(file = joinpath(modeldir, "potential_lmc_init.ini"))

    rho_peri = mean_density(pot, df.pericentre, units)
    rho_peri_lmc = mean_density(pot_lmc, df.pericentre_lmc, units)

    r_J = calc_r_J.(halo, rho_peri)
    r_J_lmc = calc_r_J.(halo, rho_peri_lmc)
    print_radius("r_J", r_J, df)
    print_radius("r_J_lmc", r_J_lmc, df)


    σv = obs_props["sigma_v"] .+ randn(size(df, 1)) * obs_props["sigma_v_err"]

	
    r_break = LilGuys.break_radius.(σv / V2KMS, df.time_last_peri)
    r_break_lmc = LilGuys.break_radius.(σv / V2KMS, replace_missings(df.time_last_peri_lmc))
	
    print_radius("r_break", r_break, df)
    print_radius("r_break_lmc", r_break_lmc, df)
    print_radius("r_break_lmc", r_break_lmc, df)

	print_quantity("pericentre", df.pericentre)
	print_quantity("pericentre lmc", df.pericentre_lmc)
	print_quantity("apocentre", df.apocentre)
	print_quantity("apocentre_lmc", df.apocentre_lmc)
	print_quantity("time last", df.time_last_peri * T2GYR)
	print_quantity("time last lmc", replace_missings(df.time_last_peri_lmc) * T2GYR)
	# print_quantity("periods", floor.(10 ./ (df.period * T2GYR)))


	return DataFrame(
		:r_break => r_break,
		:r_break_lmc => r_break_lmc,
		:r_J => r_J,
		:r_J_lmc => r_J_lmc
	)
end

# ╔═╡ 1b474924-ab4f-4e8f-b18d-b28dd116ef77
function analyze_mw(galaxyname, modelname, halo; units=Agama.AgamaUnits())
    obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))


	modeldir = joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, modelname)

    df = read_fits(joinpath(modeldir, "orbital_properties.fits"))

    pot = Agama.Potential(file = joinpath(modeldir, "agama_potential.ini"))

    rho_peri = mean_density(pot, df.pericentre, units)

    r_J = calc_r_J.(halo, rho_peri)
    print_radius("r_J", r_J, df)


    σv = obs_props["sigma_v"] .+ randn(size(df, 1)) * obs_props["sigma_v_err"]

	t_last_peri =  df.time_last_peri
	t_last_peri[DataFrames.ismissing.(t_last_peri)] .== -1
	
    r_break = LilGuys.break_radius.(σv / V2KMS, t_last_peri)
    print_radius("r_break", r_break, df)
	print_quantity("pericentre", df.pericentre)
	print_quantity("apocentre", df.apocentre)
	print_quantity("time last", df.time_last_peri * T2GYR)
	print_quantity("periods", floor.(10 ./ (df.period * T2GYR)))

	return DataFrame(
		:r_break => r_break,
		:r_J => r_J,
		:pericentre => df.pericentre,
	)
end

# ╔═╡ a19f8018-30aa-4cb9-a433-d7ff0c938af9
r_kink_scl = LilGuys.arcmin2kpc(Measurement(25, 5), Measurement(83.2, 2))

# ╔═╡ 442c6efc-38cf-4f69-b6c2-74dd3859574c
r_kink_umi = LilGuys.arcmin2kpc(Measurement(30, 5), Measurement(70.1, 3.6))

# ╔═╡ 9ab135a5-9bff-4d0b-93f3-979fd41f52d6
umi_halo = NFW(r_circ_max=4.0, v_circ_max=38/V2KMS)

# ╔═╡ e8e0236b-dac4-4196-bed6-d1122a1c4377
df_lmc_umi = analyze_lmc("ursa_minor", "L3M11_5Gyr", umi_halo)

# ╔═╡ 0a1b2a5e-6cca-44a4-bb8f-759bff2de019
df_umi = analyze_mw("ursa_minor", "EP2020", umi_halo)

# ╔═╡ e1f5f252-7d60-4320-b638-2b0a0fa34af5
scl_halo =  NFW(r_circ_max=3.2, v_circ_max=31/V2KMS)

# ╔═╡ 2aabaa14-256c-4500-9262-f9c24ce4b5bc
df_lmc = analyze_lmc("sculptor", "vasiliev24_L3M11", scl_halo)

# ╔═╡ c8f21530-9e47-471b-a8d9-80ec395c2bcb
df = analyze_mw("sculptor", "EP2020", scl_halo)

# ╔═╡ 9b0252d2-6253-4215-8570-dc7da3db369e
samples = OrderedDict(
	"MW only" => df.r_break,
	"MW+LMC (mw)" => df_lmc.r_break,
	"MW+LMC (lmc)" => df_lmc.r_break_lmc,
)

# ╔═╡ a2057cdb-febc-4fd3-a4bc-fda4c5d91cfa
samples_umi = OrderedDict(
	"MW only" => df_umi.r_break,
	"MW+LMC (mw)" => df_lmc_umi.r_break,
	#"MW+LMC (lmc)" => df_lmc_umi.r_break_lmc,
)

# ╔═╡ 41fd69fd-9836-4635-b009-8848dd6ae8ad
samples_rJ = OrderedDict(
	"MW only" => df.r_J,
	"MW+LMC (mw)" => df_lmc.r_J,
	"MW+LMC (lmc)" => df_lmc.r_J_lmc,
)

# ╔═╡ d98b72d8-ae5d-44e2-93a2-1e1e5348f677
samples_rJ_umi = OrderedDict(
	"MW only" => df_umi.r_J,
	"MW+LMC (mw)" => df_lmc_umi.r_J,
	#"MW+LMC (lmc)" => df_lmc_umi.r_J_lmc,
)

# ╔═╡ 17c66ef5-8a27-4da1-8d1d-4fe92f709771
theme(:Annotation)[:style]

# ╔═╡ 946f8dbf-287f-454f-be57-975f1c3d1209
function annotate_radii!(ax, r_jacobi, r_break)
	annotation!(ax, 0, 36, r_jacobi, 0, color=COLORS[1] )

	annotation!(ax, 0, 36, r_break, 0, color=COLORS[2])
end

# ╔═╡ af7bf2ce-5972-4917-b634-adcbd5f0a702
calc_r_J(scl_halo, mean_density(pot_ep, 42., Agama.AgamaUnits()))

# ╔═╡ e60f1827-e1de-4c9f-a60b-0f658352f19c
function plot_samples(samples, samples_rJ, r_kink; limits=(-0.5, 1))
	fig = Figure()
	bins = limits[1]:0.01:limits[end]

	axs = []
	for (i, (label, sample)) in enumerate(samples)
		ax = Axis(fig[i,1], xlabel = "Radius / kpc")
		
		hist!(log10.(abs.(samples_rJ[label])), bins=bins, normalization=:pdf, label="Jacobi radii")
		stephist!(log10.(abs.(sample)), bins=bins, normalization=:pdf, color=COLORS[2], linewidth=theme(:linewidth)[]/2, label="break radii")

		text!(0, 1, text=label, space=:relative, align=(:left, :top), offset=(6, -6))
		ylims!(0, nothing)
		xlims!(limits...)
		hideydecorations!()

		push!(axs, ax)

		if i < length(samples)
			hidexdecorations!(ax, ticks=false, minorticks=false)
		end

		vspan!(log10(r_kink.middle - r_kink.lower), log10(r_kink.middle + r_kink.upper), alpha=0.2, color=COLORS[3], label="observed break")


		if i == 2
			axislegend(position=:lb)
		end
	end

	linkxaxes!(axs...)
	rowgap!(fig.layout, 0)

	
	fig
end

# ╔═╡ eaaf7b31-8e49-4b5d-b13c-9cbc84f25f29
let 
	fig = plot_samples(samples, samples_rJ, r_kink_scl, limits=(-0.5, 0.9))
	fig.content[1].title = "Sculptor"

	annotate_radii!(fig.content[1], log10(3.5), log10(2.3))


	annotate_radii!(fig.content[2], log10(2.8), log10(1.6))
	annotate_radii!(fig.content[4], log10(2.6), log10(0.5))

	@savefig "scl_break_radii_hists"
	fig

end

# ╔═╡ 3eb3e487-4f06-466e-88a2-0aab8e4ce647
let
	fig = plot_samples(samples_umi, samples_rJ_umi, r_kink_scl, limits=(-1, 1.2))

	annotate_radii!(fig.content[1], log10(2.9), log10(3.7))

	fig.content[1].title = "Ursa Minor"
	@savefig "umi_break_radii_hists"

	fig

end

# ╔═╡ d4fae6f9-7355-42b5-b74f-22026749fd67
# df_lmc_halos = analyze_lmc("sculptor", "vasiliev24_L3M11", sample_scl_halos(100_000))

# ╔═╡ ad31388a-dc6f-4a05-853e-4e37c812446e
# df_umi_halos = analyze_mw("ursa_minor", "EP2020", sample_umi_halos(100_000))

# ╔═╡ 3b85c1df-dc0a-4ca0-8c7a-707b46cabaf9
plot_samples(OrderedDict(
	"MW only" => df_lmc_halos.r_break,
	"MW+LMC (mw)" => df_lmc.r_break,
) , OrderedDict(
	"MW only" => df_lmc_halos.r_J,
	"MW+LMC (mw)" => df_lmc.r_J,
), r_kink_scl,limits=(0, 7))

# ╔═╡ 63123436-8f01-4783-a7c8-41da3fe8e556
plot_samples(OrderedDict(
	"MW only" => df_umi.r_break,
	"MW+LMC (mw)" => df_umi_halos.r_break,
) , OrderedDict(
	"MW only" => df_umi.r_J,
	"MW+LMC (mw)" => df_umi_halos.r_J,
), r_kink_scl,limits=(0, 7))

# ╔═╡ b3957136-6f04-4528-93aa-6c355153901c
LilGuys.error_interval(r_kink_scl)

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═1368d513-2631-4adb-bf3a-cbb33cd93b34
# ╠═f0e19873-1835-4b47-9389-f4034ec4b032
# ╠═665d33fd-1d6f-40b8-b03b-b56a604d9862
# ╠═ad86ba28-77a0-42ed-ac09-7f0d295f4ebf
# ╠═e11d4bde-923f-4e50-b7b0-ed4eb9ad2ae4
# ╠═7b8afb33-83c0-444d-8333-ecd1f96c3e1b
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═2dfea445-50ff-4983-8c33-44e2d05b0cdc
# ╠═4fa535d8-0098-439b-b0c1-6ec31c110d13
# ╠═a25c4187-8409-4efb-9884-10ba9b4af3e9
# ╠═0744d9a2-256f-45ed-8ba5-36f549a08c2b
# ╠═2b03efa2-f91a-4e50-a90d-8c31a929ce1e
# ╠═d03e54bb-6311-40d3-a669-b72d017b0901
# ╠═d537bae8-8eee-4fa9-b22b-7a933d0a40a6
# ╠═6fdf86dc-0130-4c6e-a8e9-923be6315afb
# ╠═fdd9469f-140d-4118-af2c-2c389c7c8f34
# ╠═8d1daca0-5d79-47c2-9cca-1974b0c3ec65
# ╠═d8533d10-63f6-4f9e-a7e2-49809d19f6f2
# ╠═72736bff-eccc-491e-8c31-63cc52ea95ff
# ╠═a9a983f5-4570-40d3-9d6e-e22d085b332a
# ╠═1b474924-ab4f-4e8f-b18d-b28dd116ef77
# ╠═2aabaa14-256c-4500-9262-f9c24ce4b5bc
# ╠═e8e0236b-dac4-4196-bed6-d1122a1c4377
# ╠═c8f21530-9e47-471b-a8d9-80ec395c2bcb
# ╠═0a1b2a5e-6cca-44a4-bb8f-759bff2de019
# ╠═a19f8018-30aa-4cb9-a433-d7ff0c938af9
# ╠═442c6efc-38cf-4f69-b6c2-74dd3859574c
# ╠═9ab135a5-9bff-4d0b-93f3-979fd41f52d6
# ╠═e1f5f252-7d60-4320-b638-2b0a0fa34af5
# ╠═9b0252d2-6253-4215-8570-dc7da3db369e
# ╠═a2057cdb-febc-4fd3-a4bc-fda4c5d91cfa
# ╠═41fd69fd-9836-4635-b009-8848dd6ae8ad
# ╠═d98b72d8-ae5d-44e2-93a2-1e1e5348f677
# ╠═17c66ef5-8a27-4da1-8d1d-4fe92f709771
# ╠═946f8dbf-287f-454f-be57-975f1c3d1209
# ╠═af7bf2ce-5972-4917-b634-adcbd5f0a702
# ╠═eaaf7b31-8e49-4b5d-b13c-9cbc84f25f29
# ╠═3eb3e487-4f06-466e-88a2-0aab8e4ce647
# ╠═e60f1827-e1de-4c9f-a60b-0f658352f19c
# ╠═d4fae6f9-7355-42b5-b74f-22026749fd67
# ╠═ad31388a-dc6f-4a05-853e-4e37c812446e
# ╠═3b85c1df-dc0a-4ca0-8c7a-707b46cabaf9
# ╠═63123436-8f01-4783-a7c8-41da3fe8e556
# ╠═b3957136-6f04-4528-93aa-6c355153901c
