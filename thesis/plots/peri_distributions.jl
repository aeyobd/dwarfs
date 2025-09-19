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
    r_break_lmc = LilGuys.break_radius.(σv / V2KMS, df.time_last_peri_lmc)
    print_radius("r_break", r_break, df)
    print_radius("r_break_lmc", r_break_lmc, df)

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

    r_break = LilGuys.break_radius.(σv / V2KMS, df.time_last_peri)
    print_radius("r_break", r_break, df)

	return DataFrame(
		:r_break => r_break,
		:r_J => r_J,
	)
end

# ╔═╡ 2aabaa14-256c-4500-9262-f9c24ce4b5bc
df_lmc = analyze_lmc("sculptor", "vasiliev24_L3M11", NFW(r_circ_max=3.2, v_circ_max=31/V2KMS))

# ╔═╡ e8e0236b-dac4-4196-bed6-d1122a1c4377
df_lmc_umi = analyze_lmc("ursa_minor", "vasiliev24_L3M11", NFW(r_circ_max=3.2, v_circ_max=31/V2KMS))

# ╔═╡ c8f21530-9e47-471b-a8d9-80ec395c2bcb
df = analyze_mw("sculptor", "EP2020", NFW(r_circ_max=3.2, v_circ_max=31/V2KMS))

# ╔═╡ a19f8018-30aa-4cb9-a433-d7ff0c938af9
r_kink_scl = LilGuys.arcmin2kpc(Measurement(25, 5), Measurement(83.2, 2))

# ╔═╡ 442c6efc-38cf-4f69-b6c2-74dd3859574c
r_kink_umi = LilGuys.arcmin2kpc(Measurement(30, 5), Measurement(70.1, 3.6))

# ╔═╡ 0a1b2a5e-6cca-44a4-bb8f-759bff2de019
df_umi = analyze_mw("ursa_minor", "EP2020", NFW(r_circ_max=3.2, v_circ_max=31/V2KMS))

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
	"MW+LMC (lmc)" => df_lmc_umi.r_break_lmc,
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
	"MW+LMC (lmc)" => df_lmc_umi.r_J_lmc,
)

# ╔═╡ 88cd177e-2066-4833-b7a4-ef8b3ede3077
LilGuys.mean(samples_rJ_umi["MW+LMC (lmc)"] .< 0.7)

# ╔═╡ e60f1827-e1de-4c9f-a60b-0f658352f19c
function plot_samples(samples, samples_rJ, r_kink; limits=(0, 10))
	fig = Figure()
	bins = limits[1]:0.05:limits[end]

	axs = []
	for (i, (label, sample)) in enumerate(samples)
		ax = Axis(fig[i,1], xlabel = "Radius / kpc")
		
		hist!(abs.(samples_rJ[label]), bins=bins, normalization=:pdf, label="Jacobi radii")
		stephist!(abs.(sample), bins=bins, normalization=:pdf, color=COLORS[2], linewidth=theme(:linewidth)[]/2, label="break radii")

		text!(0, 1, text=label, space=:relative, align=(:left, :top), offset=(6, -6))
		ylims!(0, nothing)
		xlims!(limits...)
		hideydecorations!()

		push!(axs, ax)

		if i < length(samples)
			hidexdecorations!(ax, ticks=false, minorticks=false)
		end

		vspan!(r_kink.middle - r_kink.lower, r_kink.middle + r_kink.upper, alpha=0.2, color=COLORS[3], label="observed break")


		if i == 1
			axislegend(position=:rt)
		end
	end

	linkxaxes!(axs...)
	rowgap!(fig.layout, 0)

	
	fig
end

# ╔═╡ eaaf7b31-8e49-4b5d-b13c-9cbc84f25f29
let 
	fig = plot_samples(samples, samples_rJ, r_kink_scl, limits=(0, 7))
	fig.content[1].title = "Sculptor"

	@savefig "scl_break_radii_hists"
	fig

end

# ╔═╡ 3eb3e487-4f06-466e-88a2-0aab8e4ce647
let
	fig = plot_samples(samples_umi, samples_rJ_umi, r_kink_scl, limits=(0, 30))

	fig.content[1].title = "Ursa Minor"
	@savefig "umi_break_radii_hists"

	fig

end

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
# ╠═2b03efa2-f91a-4e50-a90d-8c31a929ce1e
# ╠═d03e54bb-6311-40d3-a669-b72d017b0901
# ╠═72736bff-eccc-491e-8c31-63cc52ea95ff
# ╠═1b474924-ab4f-4e8f-b18d-b28dd116ef77
# ╠═2aabaa14-256c-4500-9262-f9c24ce4b5bc
# ╠═e8e0236b-dac4-4196-bed6-d1122a1c4377
# ╠═c8f21530-9e47-471b-a8d9-80ec395c2bcb
# ╠═a19f8018-30aa-4cb9-a433-d7ff0c938af9
# ╠═442c6efc-38cf-4f69-b6c2-74dd3859574c
# ╠═0a1b2a5e-6cca-44a4-bb8f-759bff2de019
# ╠═9b0252d2-6253-4215-8570-dc7da3db369e
# ╠═a2057cdb-febc-4fd3-a4bc-fda4c5d91cfa
# ╠═41fd69fd-9836-4635-b009-8848dd6ae8ad
# ╠═d98b72d8-ae5d-44e2-93a2-1e1e5348f677
# ╠═eaaf7b31-8e49-4b5d-b13c-9cbc84f25f29
# ╠═3eb3e487-4f06-466e-88a2-0aab8e4ce647
# ╠═88cd177e-2066-4833-b7a4-ef8b3ede3077
# ╠═e60f1827-e1de-4c9f-a60b-0f658352f19c
# ╠═b3957136-6f04-4528-93aa-6c355153901c
