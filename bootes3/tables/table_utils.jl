import TOML
using LilGuys
import DataFrames: disallowmissing

modelnames = TOML.parsefile("model_key.toml")


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

function get_obs_props(galaxyname)
    return TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))
end

