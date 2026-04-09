using DataFrames, CSV
import TOML
using LilGuys
using StatsBase: median, quantile

module RVUtils
include(joinpath(ENV["DWARFS_ROOT"], "observations/rv_utils.jl"))
end

obs_dir = joinpath(ENV["DWARFS_ROOT"], "observations", "bootes3")
obs_props = TOML.parsefile(joinpath(obs_dir, "observed_properties.toml"))

Δv_gsr = RVUtils.rv_gsr_shift(obs_props["ra"], obs_props["dec"])


function print_properties_kinematic(samplename="mcmc_samples_vz.csv")
    filename = joinpath(obs_dir, "final_derivations/mcmc", samplename)

    samples = CSV.read(filename, DataFrame)

    println(samplename)
    print_quantity("μ (gsr)", samples.μ)
    print_quantity("μ (icrs)", samples.μ .+ Δv_gsr)
    print_quantity("σ", samples.σ )
end


function print_properties(samplename="samples.G_21.mcmc_plummer.csv")
    filename = joinpath(obs_dir, "mcmc", samplename)

    ra0, dec0 = obs_props["ra_original"], obs_props["dec_original"]

    samples = CSV.read(filename, DataFrame)
    N = size(samples, 1)
    ras = Vector{Float64}(undef, N)
    decs = Vector{Float64}(undef, N)

    for i in 1:N
        ras[i], decs[i] = LilGuys.from_tangent(samples.d_xi[i]/60, samples.d_eta[i]/60, ra0, dec0)
    end

    println(samplename)
    print_quantity("ra", ras)
    print_quantity("dec", decs)
    print_quantity("d_xi", samples.d_xi)
    print_quantity("d_eta", samples.d_eta)
    print_quantity("R_h", samples.R_h)
    print_quantity("N_sat", samples.N_sat)
    print_quantity("f_sat", samples.f_sat)

    for col in ["ellipticity", "position_angle", "n"]
        if col ∈ names(samples)
            print_quantity(col, samples[!, col])
        end
    end

    println()
end



function print_quantity(key, x::AbstractVector)
    x = disallowmissing(replace(x, missing => NaN)) # skip missings
    n_missing = LilGuys.mean(isnan.(x))
    if n_missing > 0
        @info "missing rate", n_missing
    end

    x = x[.!isnan.(x)]
    if length(x) == 0
        println(key, "\tNAN")
        return
    end

    l, m, h = LilGuys.quantile(x, [0.16, 0.5, 0.84])
    print_quantity(key, Measurement(m, m-l, h-m))
end


function print_quantity(key, m::Measurement)
    println(key, "\t", latex_string(m))
end

function print_quantity(key, x::Real)
    println(key, "\t", x)
end


function latex_string(m::Measurement)
    lo = m.lower
    hi = m.upper
    mid = m.middle

    # find order of magnitude of uncertainties
    sig_lo = lo == 0 ? -3 : floor(Int, log10(lo))
    sig_hi = hi == 0 ? -3 : floor(Int, log10(hi))
    sig = min(sig_lo, sig_hi)

    # we want 1 or 2 sig figs for uncertainty
    round_digits = -sig + 1

    lo_r = round(lo; sigdigits=2)
    hi_r = round(hi; sigdigits=2)
    mid_r = round(mid; digits=round_digits)

    return "\$$(mid_r)^{+$(hi_r)}_{-$(lo_r)}\$"
end


function (@main)(ARGS)
    print_properties_kinematic()

    println()

    print_properties("../final_derivations/mcmc/samples.j24_1c.mcmc_exp_ell.csv")
    print_properties("../final_derivations/mcmc/samples.j24_1c.mcmc_ell.csv")
    print_properties("../final_derivations/mcmc/samples.j24_1c.mcmc_sersic_ell.csv")

    println()
    println("delve")
    print_properties("samples.delve_matched_filter_6deg.mcmc_ell.csv")
    print_properties("samples.delve_matched_filter_6deg.g22.mcmc_ell.csv")
end

