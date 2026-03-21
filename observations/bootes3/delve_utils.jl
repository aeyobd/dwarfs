
import LilGuys
using CSV, DataFrames

ISO_HEADER = "Zini     MH   logAge Mini        int_IMF         Mass   logL    logTe  logg  label   McoreTP C_O  period0  period1  period2  period3  period4  pmode  Mloss  tau1m   X   Y   Xc  Xn  Xo  Cexcess  Z 	 mbolmag  umag    gmag    rmag    imag    zmag    Ymag"


"""
    get_isochrone(M_H::Real, age::Real=12)

Loads and filters a Padova isochrone of a given metallicity and age.
Returns a DataFrame containing stellar evolution quantities and magnitudes.
"""
function get_isochrone(M_H::Real, age::Real=12)
    iso_columns = string.(split(ISO_HEADER, r"\s+"))

    all_isochrones = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/padova/isochrone.decam.$(age)Gyrs.dat"),DataFrame,
					  comment="#", ignorerepeated=true, delim = ' ', header=iso_columns)

	
	M_Hs = unique(all_isochrones.MH)
	M_H_adopted = M_Hs[argmin(abs.(M_H .- M_Hs))]
	@info "using isochrone with metallicity $M_H_adopted, age $age"

	filt = isapprox.(all_isochrones.MH, M_H_adopted)
	filt .&= all_isochrones.label .< 5 # only keep through HB

	all_isochrones[filt, :]
end



"""
    kroupa_imf(m) -> Float64

Kroupa (2001) broken power-law IMF. Returns dN/dm ∝ m^{-α}.
"""
function kroupa_imf(m::Real)
    m <= 0 && return 0.0
    if m < 0.08
        return m^(-0.3)
    elseif m < 0.50
        return 0.08^(1.3 - 0.3) * m^(-1.3)   # continuous at 0.08
    else
        return 0.08^(1.3 - 0.3) * 0.5^(2.3 - 1.3) * m^(-2.3)
    end
end

"""
    chabrier_imf(m) -> Float64

Chabrier (2003) IMF: log-normal below 1 M☉, Salpeter power-law above.
"""
function chabrier_imf(m::Real)
    m <= 0 && return 0.0
    if m <= 1.0
        return exp(-((log10(m) - log10(0.22))^2) / (2 * 0.57^2)) / m
    else
        return chabrier_imf(1.0) * m^(-2.3)
    end
end

# ╔═╡ bb0600b7-cf8a-484b-9d4a-7dab2c963da6
err_model(x, params) = @. (params[1] + x*params[2] + x^2 * params[3])

function fit_color_err(stars)

    gr_err = @. sqrt(stars.magerr_psf_g^2 + stars.magerr_psf_r^2)

    popt, covt = LilGuys.curve_fit(err_model, stars.gmag, log10.(gr_err), [0., 1., 1.])

    color_err(gmag) = 10 .^ err_model(gmag, popt)

    return color_err
end

function fit_color_err(gmag, gr_err)
    popt, covt = LilGuys.curve_fit(err_model, gmag, log10.(gr_err), [0., 1., 1.])

    color_err(gmag) = 10 .^ err_model(gmag, popt)

    return color_err
end

